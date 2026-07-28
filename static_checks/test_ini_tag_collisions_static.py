"""Guards the INI tag-matching contract.

The parser used to match a wanted tag against a file key by SUBSTRING, taking
the first hit in file order. That silently bound 'err_label' to
'fluxnet_err_label', so the missing-value token written into every output file
came from the wrong key. It also made indexed keys such as 'gas_1_col' unsafe
to introduce, because 'col_co2' is a substring of 'col_co2_something'.

SearchLocalTags now matches by exact equality. These checks keep it that way
and keep the tag tables free of the ambiguities that motivated the change.
"""

import re
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def read(path):
    return (ROOT / path).read_text(encoding="utf-8", errors="replace")


# e.g.:  data EPPrjNTags(3)%Label / 'col_ts' / &
#             EPPrjCTags(36)%Label / 'err_label'        / &
TAG_RE = re.compile(
    r"^\s*(?:data\s+)?([A-Za-z]\w*Tags)\((\d+)\)%Label\s*/\s*'([^']*)'\s*/"
)

# module file -> tag tables declared in it, and the INI scope each is parsed with.
# Tables sharing a scope are searched against the same set of file lines.
TABLES = {
    "src/src_common/m_common_global_var.f90": {
        "EPPrjNTags": "project",
        "EPPrjCTags": "project",
        "ANTags": "metadata",
        "ACTags": "metadata",
    },
    "src/src_rp/m_rp_global_var.f90": {
        "SNTags": "rawprocess",
        "SCTags": "rawprocess",
    },
    "src/src_fcc/m_fx_global_var_mod.f90": {
        "SNTags": "fluxcorrection",
        "SCTags": "fluxcorrection",
    },
}


def tag_tables():
    """(module, table) -> {slot index: label}, skipping commented-out slots."""
    out = defaultdict(dict)
    for rel, names in TABLES.items():
        for line in read(rel).splitlines():
            if line.lstrip().startswith("!"):
                continue
            m = TAG_RE.match(line)
            if m and m.group(1) in names:
                out[(rel, m.group(1))][int(m.group(2))] = m.group(3)
    return out


def test_tag_tables_are_discoverable():
    """If the data-block format changes, the other checks here go silently vacuous."""
    tables = tag_tables()
    assert len(tables) == 8, sorted(tables)
    for key, slots in tables.items():
        assert slots, f"no labels parsed for {key}"


def test_search_local_tags_matches_tag_names_exactly():
    source = read("src/src_common/parse_ini_file.f90")
    body = source[source.index("subroutine SearchLocalTags"):]

    assert "trim(adjustl(Tags(j)%Label)) == trim(adjustl(NumTags(i)%label))" in body
    assert "trim(adjustl(Tags(j)%Label)) == trim(adjustl(CharTags(i)%label))" in body
    # the substring form must not come back
    assert "index(Tags(j)%Label" not in body


def test_search_local_tags_does_not_declare_wanted_tags_intent_out():
    """%label is read inside the routine, so intent(out) would be undefined on entry."""
    body = read("src/src_common/parse_ini_file.f90")
    body = body[body.index("subroutine SearchLocalTags"):]
    assert "type(Numerical), intent(inout) :: NumTags(nnum)" in body
    assert "type(Text), intent(inout) :: CharTags(nchar)" in body


def test_no_duplicate_labels_within_a_tag_table():
    """Two slots sharing a label both bind to the same key -- always a mistake."""
    for (rel, name), slots in tag_tables().items():
        seen = defaultdict(list)
        for index, label in slots.items():
            if label:
                seen[label].append(index)
        dupes = {k: v for k, v in seen.items() if len(v) > 1}
        assert not dupes, f"{name} in {rel} has duplicate labels: {dupes}"


def test_no_tag_label_is_a_substring_of_another_in_the_same_scope():
    """Exact matching makes these harmless, but they signal a naming mistake.

    Known and accepted: the legacy 'prof_<gas>' profile tags are prefixes of
    'prof_<gas>_z1'..'_z7'. The GUI writes no 'prof_' keys at all, so they are
    dead weight rather than a live hazard; they are listed explicitly so that a
    NEW collision cannot slip in unnoticed.
    """
    accepted = {
        (f"prof_{gas}", f"prof_{gas}_z{z}")
        for gas in ("co2", "h2o", "ch4", "gas4")
        for z in range(1, 8)
    }
    # 'err_label' / 'fluxnet_err_label' is the bug this module exists for. It is
    # deliberately NOT accepted here: exact matching is what makes it safe, and
    # test_search_local_tags_matches_tag_names_exactly is what enforces that.
    accepted.add(("err_label", "fluxnet_err_label"))
    # Same shape in the metadata file: a standalone 'head_corr' alongside the
    # per-instrument 'instr_<K>_head_corr'. Under the old substring matching the
    # standalone tag would bind to whichever instrument block came first in the
    # file; exact matching is what makes it harmless.
    n_instr = int(
        re.search(
            r"integer, parameter :: MaxNumInstruments = (\d+)",
            read("src/src_common/m_typedef.f90"),
        ).group(1)
    )
    accepted |= {("head_corr", f"instr_{k}_head_corr") for k in range(1, n_instr + 2)}

    by_scope = defaultdict(set)
    for rel, names in TABLES.items():
        for name, scope in names.items():
            by_scope[scope] |= {
                label for label in tag_tables()[(rel, name)].values() if label
            }

    unexpected = []
    for scope, labels in by_scope.items():
        for short in labels:
            for long in labels:
                if short != long and short in long and (short, long) not in accepted:
                    unexpected.append((scope, short, long))
    assert not unexpected, f"new tag-name collisions: {unexpected}"


def test_metadata_tag_tables_match_the_generator():
    """The .metadata tables are generated; a hand edit or a changed
    MaxNumInstruments without re-running gen_metadata_tags.py silently
    desynchronises the tag indices from the stride arithmetic that reads them."""
    import subprocess
    import sys

    gen = ROOT / "prj" / "gen_metadata_tags.py"
    assert gen.exists(), gen

    n = int(
        re.search(
            r"integer, parameter :: MaxNumInstruments = (\d+)",
            read("src/src_common/m_typedef.f90"),
        ).group(1)
    )
    r = subprocess.run(
        [sys.executable, str(gen), "--instruments", str(n), "--check"],
        capture_output=True, text=True, cwd=str(ROOT / "prj"),
    )
    assert r.returncode == 0, (
        f"metadata tag tables are stale for MaxNumInstruments={n}; "
        f"re-run prj/gen_metadata_tags.py --instruments {n}\n{r.stdout}{r.stderr}"
    )


def test_project_tag_tables_match_the_generator():
    """The [Project]/RP/FCC tables are generated too; a hand edit desynchronises
    the tag indices from the positional code that consumes them."""
    import subprocess
    import sys

    gen = ROOT / "prj" / "gen_project_tags.py"
    assert gen.exists(), gen
    r = subprocess.run(
        [sys.executable, str(gen), "--check"],
        capture_output=True, text=True, cwd=str(ROOT / "prj"),
    )
    assert r.returncode == 0, (
        f"project tag tables are stale; re-run prj/gen_project_tags.py"
        f"\n{r.stdout}{r.stderr}"
    )


def test_metadata_group_origins_match_the_generated_layout():
    """read_metadata_file.f90 addresses the groups by arithmetic, so its
    origins must match where the generator actually put them."""
    src = read("src/src_common/read_metadata_file.f90")
    tables = tag_tables()
    an = tables[("src/src_common/m_common_global_var.f90", "ANTags")]
    ac = tables[("src/src_common/m_common_global_var.f90", "ACTags")]

    an_col = min(i for i, l in an.items() if l.startswith("col_1_"))
    ac_col = min(i for i, l in ac.items() if l.startswith("col_1_"))
    assert f"init_an_col = {an_col} - leap_an_col" in src
    assert f"init_ac_col = {ac_col} - leap_ac_col" in src

    data_label = next(i for i, l in ac.items() if l == "data_label")
    assert f"ACTags({data_label})%value" in src


def test_full_spectra_gas4_tag_matches_the_key_the_gui_writes():
    """The GUI writes 'out_full_sp_gas4'; the engine used to look for
    'out_full_sp_n2o', so the setting was silently dropped and the flag always
    read false. Every other 4th-gas tag already uses the '_gas4' spelling."""
    source = read("src/src_rp/m_rp_global_var.f90")
    assert "'out_full_sp_gas4'" in source
    assert "'out_full_sp_n2o'" not in source
    # still consumed from the same slot
    assert "RPsetup%out_full_sp(gas4) = SCTags(34)%value(1:1) == '1'" in read(
        "src/src_rp/read_ini_rp.f90"
    )
