#!/usr/bin/env python3
"""Regenerate the project-file tag tables (EPPrj*, RP S*Tags, FCC S*Tags).

Companion to gen_metadata_tags.py, same contract: these tables are positional,
so a tag's identity is its array index and the blocks must stay contiguous and
consistent with the code that reads them by index.

This step is deliberately a NO-OP on the tag set: it takes over ownership of the
blocks and re-emits exactly what is already there, one assignment per line and
chunked under the Fortran 255-continuation limit. Widening comes separately, so
that "the generator reproduces reality" is proved before anything changes.

Gaps are preserved: several slots are deliberately commented out (retired keys),
and only assigned indices are emitted, so those stay unassigned.

Usage:  python gen_project_tags.py [--check]

--check exits non-zero if a file is not what the generator would produce.
"""

import argparse
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

# (path, table name, marker name, size parameter). Marker names are unique per
# block because two different modules both declare tables called SNTags/SCTags;
# the size parameter is rewritten in the same file, so Nsn/Nsc are unambiguous.
BLOCKS = [
    ("src/src_common/m_common_global_var.f90", "EPPrjNTags", "EPPrjNTags", "Npn"),
    ("src/src_common/m_common_global_var.f90", "EPPrjCTags", "EPPrjCTags", "Npc"),
    ("src/src_rp/m_rp_global_var.f90", "SNTags", "RP.SNTags", "Nsn"),
    ("src/src_rp/m_rp_global_var.f90", "SCTags", "RP.SCTags", "Nsc"),
    ("src/src_fcc/m_fx_global_var_mod.f90", "SNTags", "FCC.SNTags", "Nsn"),
    ("src/src_fcc/m_fx_global_var_mod.f90", "SCTags", "FCC.SCTags", "Nsc"),
]

MAX_CONT = 200

# --------------------------------------------------------------------------
# Record schema. One indexed record per measurement, mirroring the way the
# metadata file already stores instr_<K>_* and col_<N>_*, so the engine can
# walk them with the same stride arithmetic.
#
# New slots are APPENDED after each table's current maximum index. Existing
# indices must never move: consumers address these tables positionally
# (SCTags(34)%value and friends), so an insertion silently re-wires settings
# with no compile error.
# --------------------------------------------------------------------------

# [Project] scalars introduced alongside the records
PROJECT_COUNTS = ["gas_num", "cell_num", "diag_num"]

#: `fluxnet_default` marks the record whose FLUXNET columns carry the bare
#: species name. CO2, H2O and CH4 are required FLUXNET variables, so the
#: standard spelling has to exist: the designated record is FC/CO2 and further
#: records of that species are CO2_2, CO2_3. Unflagged, the lowest record index
#: is designated, which is what every project written before this key expects.
GAS_NUMERIC = ["col", "moist", "cell", "mw", "diff", "fluxnet_default"]
GAS_TEXT = ["var", "instr"]
#: One Conditional Eddy Covariance pairing. `co2` and `h2o` are gas RECORD
#: indices, not raw column numbers, so a project that is re-ordered keeps its
#: pairings. `extra` is a comma-separated list of further record indices whose
#: species are partitioned in the same pair's octants - carbonyl sulfide, for
#: instance - which is why it is a text key rather than a fixed run of numbers.
#:
#: Appended AFTER the diag block, so the gas/cell/diag origins do not move.
CEC_NUMERIC = ["meth", "co2", "h2o"]
CEC_TEXT = ["extra"]

CELL_NUMERIC = ["col"]
CELL_TEXT = ["var", "instr"]
DIAG_NUMERIC = ["col"]
DIAG_TEXT = ["var", "instr"]

# Per-gas processing settings, replacing the single *_gas4 key of each family.
#> Append only. Every name here is a slot at a fixed offset from the record
#> origin, and read_ini_rp.f90 addresses them by that offset, so inserting in
#> the middle silently repoints every setting after it.
#>
#> drift_dir_0..6 and drift_inv_0..6 are the direct and inverse
#> absorptance<->density calibration polynomials. They exist per gas because
#> they are per-channel instrument calibrations: the legacy tag table carries
#> a set for each of the four historical slots and nothing beyond, so a fifth
#> gas had no way to be drift-corrected. A gas that supplies none keeps
#> DriftCorr's `error` and is skipped - drift stays opt-in per gas, since
#> there is no general polynomial for an arbitrary species on an arbitrary
#> analyser.
RP_GAS_NUMERIC = [
    "sr_lim", "al_min", "al_max", "ds_hf", "ds_sf", "tl_def",
    "to_min_flux", "to_min_lag", "to_max_lag", "pwb_min_lag", "pwb_max_lag",
] + [f"drift_dir_{k}" for k in range(7)] \
  + [f"drift_inv_{k}" for k in range(7)]
RP_GAS_TEXT = ["out_full_sp", "out_full_cosp_w", "out_raw"]
FCC_GAS_NUMERIC = [
    "sa_fmin", "sa_fmax", "sa_hfn_fmin", "sa_min_st", "sa_min_un", "sa_max",
]
#> The months this gas pools before a transfer function is fitted, written as
#> a group list: `1-12` is one group over the calendar, `1-6,7-12` is two. A
#> group's ordinal in the list is its class index.
#>
#> Text rather than 24 numeric slots per gas. The flat form this replaces was
#> sa_<slot>_g<k>_start / _stop, twelve pairs for each of three slots - 72
#> tags to express a setting every real project leaves at "all months". Per
#> gas that shape would be 1536 tags, almost all of them empty. The engine
#> parses the list in ParseMonthGrouping.
FCC_GAS_TEXT = ["sa_months"]


TYPEDEF = (ROOT / "src" / "src_common" / "m_typedef.f90").read_text(
    encoding="utf-8", errors="surrogateescape")


def const(name):
    """A capacity constant, read from m_typedef.f90 so it cannot drift."""
    m = re.search(rf"integer, parameter :: {name} = (\d+)", TYPEDEF)
    if not m:
        raise SystemExit(f"could not read {name} from m_typedef.f90")
    return int(m.group(1))


#: The most groups a gas can pool its months into - twelve, because there are
#: twelve months. Read rather than written here so the retired-tag set below
#: cannot fall out of step with the engine's own bound.
MAX_GAS_CLASSES = const("MaxGasClasses")

#: How many CO2/water pairings a project may declare.
CEC_PAIRS = const("MaxNumCecPairs")


# --------------------------------------------------------------------------
# Hand-placed tags.
#
# A new setting is normally appended, which is safe for the record blocks
# because they are appended too. It is NOT safe for anything else: the record
# block starts at max(kept) + 1, so a tag written past the last hand-written
# one lifts that origin and re-emits all 1600 per-gas slots at new indices -
# rpGasOriginN went 425 -> 2033 and Nsn 2024 -> 3632 the first time these eight
# were added that way.
#
# So they go into slots that were never assigned, BELOW the origin, and
# process() refuses to write them if that origin moves anyway.
#
# 357..369 is the longer of the two free runs in the RP numerical table
# (the other is 285..289 and 291..299). Eight of the thirteen are used here.
# --------------------------------------------------------------------------
INSTR_LACK_ORIGIN = 357

#: How much of its OWN expected data an instrument may be missing, in percent.
#: Absent means "use the project-wide max_lack", which is every project written
#: before this key existed. Keyed by the same 1-based index the .metadata uses
#: for instr_<K>_*, so instr_3_max_lack is the allowance of instr_3_model.
#: Tags that occupy a stated index rather than being appended. A new [Project]
#: scalar goes into one of the blanks left by a retired key, never onto the end
#: of the table: the record blocks begin at gasNumTag and every per-gas setting
#: after it is addressed by its offset from that origin, so appending would
#: silently repoint all of them.
FIXED_TAGS = {
    "RP.SNTags": {
        INSTR_LACK_ORIGIN + k - 1: f"instr_{k}_max_lack"
        for k in range(1, const("MaxNumInstruments") + 1)
    },
    "EPPrjNTags": {
        4: "cec_singular_band",
        5: "cec_stationarity_mode",
    },
}


def limits():
    def derived(name, factor_of):
        """MaxNumCellCols / MaxNumDiagCols are declared as multiples."""
        m = re.search(rf"integer, parameter :: {name} = {factor_of} \* (\d+)",
                      TYPEDEF)
        if not m:
            raise SystemExit(f"could not read {name} from m_typedef.f90")
        return const(factor_of) * int(m.group(1))

    return {
        "gases": const("MaxNumGases"),
        "cells": derived("MaxNumCellCols", "MaxNumInstruments"),
        "diags": derived("MaxNumDiagCols", "MaxNumInstruments"),
    }


def appended(marker, lim):
    """Labels to append for a block, in order."""
    g, c, d = lim["gases"], lim["cells"], lim["diags"]
    out = []
    if marker == "EPPrjNTags":
        out += PROJECT_COUNTS
        out += [f"gas_{i}_{s}" for i in range(1, g + 1) for s in GAS_NUMERIC]
        out += [f"cell_{i}_{s}" for i in range(1, c + 1) for s in CELL_NUMERIC]
        out += [f"diag_{i}_{s}" for i in range(1, d + 1) for s in DIAG_NUMERIC]
        out += ["cec_num"]
        out += [f"cec_{i}_{s}" for i in range(1, CEC_PAIRS + 1) for s in CEC_NUMERIC]
    elif marker == "EPPrjCTags":
        out += [f"gas_{i}_{s}" for i in range(1, g + 1) for s in GAS_TEXT]
        out += [f"cell_{i}_{s}" for i in range(1, c + 1) for s in CELL_TEXT]
        out += [f"diag_{i}_{s}" for i in range(1, d + 1) for s in DIAG_TEXT]
        out += [f"cec_{i}_{s}" for i in range(1, CEC_PAIRS + 1) for s in CEC_TEXT]
    elif marker == "RP.SNTags":
        out += [f"gas_{i}_{s}" for i in range(1, g + 1) for s in RP_GAS_NUMERIC]
    elif marker == "RP.SCTags":
        out += [f"gas_{i}_{s}" for i in range(1, g + 1) for s in RP_GAS_TEXT]
    elif marker == "FCC.SNTags":
        out += [f"gas_{i}_{s}" for i in range(1, g + 1) for s in FCC_GAS_NUMERIC]
    elif marker == "FCC.SCTags":
        out += [f"gas_{i}_{s}" for i in range(1, g + 1) for s in FCC_GAS_TEXT]
    return out


def begin(name):
    return f"    !> BEGIN GENERATED {name} - edit gen_project_tags.py, not this block"


def end(name):
    return f"    !> END GENERATED {name}"


def region(text, name):
    b, e = begin(name), end(name)
    if b not in text or e not in text:
        raise SystemExit(f"marker pair missing for {name}")
    head, rest = text.split(b, 1)
    body, tail = rest.split(e, 1)
    return head, body, tail


#: Tags retired with the 5.0.0 record format. Gases, cell measurements and
#: diagnostics are described by records now, which name the analyser as well
#: as the column and can carry the same species more than once.
#:
#: Blanked rather than removed - see process(). col_ts, col_air_t and
#: col_air_p are absent deliberately: they are one per project, not one per
#: instrument, and are still live.
RETIRED_LABELS = {
    #: The fixed full-output format, which named four gas blocks whatever the
    #: project held. It described a column set no reader in this fork has been
    #: checked against, its two header branches had already drifted apart, and
    #: a fifth gas simply fell out of the file. The full output covers every
    #: configured gas now, so there is nothing left for the flag to select.
    "fix_out_format",
    "col_co2", "col_h2o", "col_ch4", "col_gas4",
    "col_cell_t", "col_int_t_1", "col_int_t_2", "col_int_p",
    "col_diag_72", "col_diag_75", "col_diag_77", "col_diag_anem",
    "gas_mw", "gas_diff",
}

#: Retired in one table only. The ru_* keys are genuine EPPrjNTags entries that
#: the engine reads; the copies in the RP table were duplicates nothing ever
#: read, and their presence there made the keys look like RawProcess settings -
#: which is where the interface wrote them, and so where the engine never
#: looked. Scoped per table, because blanking them everywhere would delete the
#: live ones.
#: The flat per-gas settings, retired with the 5.0.0 record format.
#:
#: Each of these named one of the four legacy slots - co2, h2o, ch4, gas4 -
#: and the engine read it into that slot before letting a record override it.
#: Nothing reads them now: a gas states its own thresholds through
#: gas_<i>_<setting>, which reaches every gas rather than the first four.
#:
#: Built rather than typed, because the list is a hundred labels and a typo
#: would blank a live one. Note what is deliberately absent: sa_max_h,
#: sa_max_le and sa_max_ustar are whole-run thresholds that merely share the
#: sa_max_ prefix, and out_full_cosp_w_u/v/ts are anemometric.
_SLOTS = ("co2", "h2o", "ch4", "gas4")


def _flat_per_gas(templates, slots=_SLOTS):
    return {t.format(s=s) for t in templates for s in slots}


RETIRED_RP_NUMERIC = _flat_per_gas([
    "sr_lim_{s}",
    "al_{s}_min", "al_{s}_max",
    "ds_hf_{s}", "ds_sf_{s}",
    "tl_def_{s}",
    "to_{s}_min_lag", "to_{s}_max_lag",
    "pwb_{s}_min_lag", "pwb_{s}_max_lag",
] + ["drift_dir_{s}_%d" % k for k in range(7)]
  + ["drift_inv_{s}_%d" % k for k in range(7)])
#: Water is judged by LE, so these three never had an h2o member.
RETIRED_RP_NUMERIC |= _flat_per_gas(["to_{s}_min_flux"],
                                    ("co2", "ch4", "gas4"))

RETIRED_RP_TEXT = _flat_per_gas([
    "out_full_sp_{s}", "out_full_cosp_w_{s}", "out_raw_{s}",
])

RETIRED_FCC_NUMERIC = _flat_per_gas([
    "sa_fmin_{s}", "sa_fmax_{s}", "sa_hfn_{s}_fmin",
])
RETIRED_FCC_NUMERIC |= _flat_per_gas(
    ["sa_min_st_{s}", "sa_min_un_{s}", "sa_max_{s}"], ("co2", "ch4", "gas4"))
#: The month grouping, twelve start/stop pairs for each of three slots. A gas
#: states its own now, as the single string gas_<i>_sa_months, so these 72
#: tags describe nothing: water never had a table, and every gas past the
#: fourth could only inherit CO2's.
#:
#: Built rather than typed, for the same reason as the sets above.
RETIRED_FCC_NUMERIC |= {
    "sa_%s_g%d_%s" % (s, k, e)
    for s in ("co2", "ch4", "gas4")
    for k in range(1, MAX_GAS_CLASSES + 1)
    for e in ("start", "stop")
}

#: The biomet gas profile, retired with the storage block that was its only
#: reader. prof_t_z1..z7 and one set per gas fed bSetup%zT/zCO2/... and the dz
#: built from them; that block has been commented out throughout this fork's
#: history, so the settings were parsed and consumed by nothing.
#:
#: biom_ta..biom_rg are NOT here - bSetup%sel reads those and they are live.
RETIRED_PROFILE = (
    {"prof_ts", "prof_ta"}
    | {"prof_%s" % s for s in _SLOTS}
    | {"biom_%s" % s for s in _SLOTS}
    | {"prof_t_z%d" % k for k in range(1, 8)}
    | {"prof_%s_z%d" % (s, k) for s in _SLOTS for k in range(1, 8)}
)

RETIRED_LABELS_BY_TABLE = {
    "RP.SNTags": {"ru_meth", "ru_its_meth", "ru_its_sec_factor",
                  "ru_tlag_max"} | RETIRED_RP_NUMERIC | RETIRED_PROFILE,
    "RP.SCTags": RETIRED_RP_TEXT,
    "FCC.SNTags": RETIRED_FCC_NUMERIC,
}


def parse_block(body, table):
    """index -> label for one block, skipping commented-out slots."""
    pat = re.compile(rf"{re.escape(table)}\((\d+)\)%Label\s*/\s*'([^']*)'\s*/")
    out = {}
    for line in body.splitlines():
        if line.lstrip().startswith("!"):
            continue
        for m in pat.finditer(line):
            out[int(m.group(1))] = m.group(2)
    return out


def emit(table, entries, width=26):
    lines = []
    items = sorted(entries.items())
    for start in range(0, len(items), MAX_CONT):
        chunk = items[start:start + MAX_CONT]
        for i, (idx, label) in enumerate(chunk):
            lead = (f"    data {table}({idx})%Label" if i == 0
                    else f"         {table}({idx})%Label")
            cont = " /" if i == len(chunk) - 1 else " / &"
            lines.append(f"{lead:<{width}} / '{label}'{cont}")
    return lines


def process(path, table, marker, size_param, lim, check):
    p = ROOT / path
    text = p.read_text(encoding="utf-8", errors="surrogateescape")
    head, body, tail = region(text, marker)
    entries = parse_block(body, table)
    if not entries:
        raise SystemExit(f"no {table} entries parsed in {path} - marker misplaced?")

    # Drop any previously appended record slots, then re-append. Without this
    # the generator would not be idempotent: a second run would stack another
    # full set of records on top of the first.
    record = re.compile(r"^(gas|cell|diag|cec)_\d+_|^(gas|cell|diag|cec)_num$")
    kept = {i: l for i, l in entries.items() if not record.match(l)}

    # Retired tags keep their slot and lose their label.
    #
    # These tables are positional: the reader addresses them by index, so
    # DELETING an entry would renumber every tag after it and silently
    # rebind hundreds of settings. Blanking leaves the slot in place and
    # unmatchable, which is the convention the tables already use for the
    # other retired slots.
    retired = RETIRED_LABELS | RETIRED_LABELS_BY_TABLE.get(marker, set())
    for i, label in list(kept.items()):
        if label in retired:
            kept[i] = ""

    # Where the record block starts before the hand-placed tags go in. Taking
    # it first is what lets the check below tell "filled a gap" from "appended
    # past the end", which look identical in the finished table.
    origin = max(kept) + 1
    for idx, label in sorted(FIXED_TAGS.get(marker, {}).items()):
        # Re-placing the tag the previous run wrote is the idempotent case and
        # has to be allowed; anything else at that index is a live tag being
        # overwritten.
        if kept.get(idx) not in (None, "", label):
            raise SystemExit(
                f"{marker}: slot {idx} already holds {kept[idx]!r}, so "
                f"{label} would overwrite a live tag - pick a free slot")
        kept[idx] = label

    base = max(kept) + 1
    if base != origin:
        raise SystemExit(
            f"{marker}: the hand-placed tags moved the record origin from "
            f"{origin} to {base}, which re-indexes every appended record slot. "
            f"They must go in slots below {origin}.")
    nxt = base
    for label in appended(marker, lim):
        kept[nxt] = label
        nxt += 1

    new_body = "\n" + "\n".join(emit(table, kept)) + "\n"
    new = head + begin(marker) + new_body + end(marker) + tail

    # Never shrink a table. Some tables carry unlabelled spare slots past the
    # last assigned index; those are unread today, but reducing the declared
    # size to the last label would break any code that reaches one.
    cur = re.search(rf"integer, parameter :: {size_param} = (\d+)", new)
    size = max(max(kept), int(cur.group(1))) if cur else max(kept)
    new = re.sub(rf"integer, parameter :: {size_param} = \d+",
                 f"integer, parameter :: {size_param} = {size}", new)

    changed = new != text
    if not check and changed:
        p.write_text(new, encoding="utf-8", errors="surrogateescape")
    return kept, len(appended(marker, lim)), changed, base


ORIGIN_MARKER = "ProjectRecordOrigins"


def origins_block(origins, lim):
    """Fortran parameters naming where each record group starts.

    The reader walks these groups by stride arithmetic, exactly as
    read_metadata_file.f90 does for instr_<K>_* and col_<N>_*. Emitting the
    origins here keeps the reader free of magic numbers that would silently
    rot the next time a slot is appended.
    """
    L = [
        "    !> Slot origins for the appended gas/cell/diag records. The value",
        "    !> is the index of the FIRST field of record 1, so record i field f",
        "    !> is <origin> + (i-1)*<leap> + f, with f zero-based.",
    ]
    for name, value in origins:
        L.append(f"    integer, parameter :: {name} = {value}")
    L.append(f"    integer, parameter :: gasRecLeapN   = {len(GAS_NUMERIC)}")
    L.append(f"    integer, parameter :: gasRecLeapC   = {len(GAS_TEXT)}")
    L.append(f"    integer, parameter :: cellRecLeapN  = {len(CELL_NUMERIC)}")
    L.append(f"    integer, parameter :: cellRecLeapC  = {len(CELL_TEXT)}")
    L.append(f"    integer, parameter :: diagRecLeapN  = {len(DIAG_NUMERIC)}")
    L.append(f"    integer, parameter :: diagRecLeapC  = {len(DIAG_TEXT)}")
    L.append(f"    integer, parameter :: rpGasLeapN    = {len(RP_GAS_NUMERIC)}")
    L.append(f"    integer, parameter :: rpGasLeapC    = {len(RP_GAS_TEXT)}")
    L.append(f"    integer, parameter :: cecRecLeapN   = {len(CEC_NUMERIC)}")
    L.append(f"    integer, parameter :: cecRecLeapC   = {len(CEC_TEXT)}")
    L.append(f"    integer, parameter :: fccGasLeapN   = {len(FCC_GAS_NUMERIC)}")
    L.append(f"    integer, parameter :: fccGasLeapC   = {len(FCC_GAS_TEXT)}")
    return L


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    args = ap.parse_args()

    lim = limits()
    print(f"capacity: {lim['gases']} gases, {lim['cells']} cell columns, "
          f"{lim['diags']} diagnostic columns\n")

    stale = []
    origins = []
    for path, table, marker, size_param in BLOCKS:
        entries, n_new, changed, base = process(
            path, table, marker, size_param, lim, args.check)
        gaps = sorted(set(range(1, max(entries) + 1)) - set(entries))
        print(f"{marker:<12} {len(entries):>5} tags  {size_param}={max(entries):<5}"
              f" (+{n_new} record slots)"
              f"{'  gaps: ' + str(len(gaps)) if gaps else ''}")
        if changed:
            stale.append(marker)

        g, c = lim["gases"], lim["cells"]
        if marker == "EPPrjNTags":
            d = lim["diags"]
            diag_n = base + 3 + g * len(GAS_NUMERIC) + c * len(CELL_NUMERIC)
            cec_num_n = diag_n + d * len(DIAG_NUMERIC)
            origins += [("gasNumTag", base), ("cellNumTag", base + 1),
                        ("diagNumTag", base + 2),
                        ("gasRecOriginN", base + 3),
                        ("cellRecOriginN", base + 3 + g * len(GAS_NUMERIC)),
                        ("diagRecOriginN", diag_n),
                        ("cecNumTag", cec_num_n),
                        ("cecRecOriginN", cec_num_n + 1)]
        elif marker == "EPPrjCTags":
            d = lim["diags"]
            diag_c = base + g * len(GAS_TEXT) + c * len(CELL_TEXT)
            origins += [("gasRecOriginC", base),
                        ("cellRecOriginC", base + g * len(GAS_TEXT)),
                        ("diagRecOriginC", diag_c),
                        ("cecRecOriginC", diag_c + d * len(DIAG_TEXT))]
        elif marker == "RP.SNTags":
            origins.append(("rpGasOriginN", base))
            origins.append(("rpInstrMaxLackN", INSTR_LACK_ORIGIN))
        elif marker == "RP.SCTags":
            origins.append(("rpGasOriginC", base))
        elif marker == "FCC.SNTags":
            origins.append(("fccGasOriginN", base))
        elif marker == "FCC.SCTags":
            origins.append(("fccGasOriginC", base))

    # The origins live beside the [Project] tables, where every reader can see
    # them via m_common_global_var.
    p = ROOT / "src" / "src_common" / "m_common_global_var.f90"
    text = p.read_text(encoding="utf-8", errors="surrogateescape")
    head, _, tail = region(text, ORIGIN_MARKER)
    new = (head + begin(ORIGIN_MARKER) + "\n"
           + "\n".join(origins_block(origins, lim)) + "\n"
           + end(ORIGIN_MARKER) + tail)
    if new != text:
        if args.check:
            stale.append(ORIGIN_MARKER)
        else:
            p.write_text(new, encoding="utf-8", errors="surrogateescape")

    if args.check:
        if stale:
            print(f"\nstale: {', '.join(stale)} - re-run gen_project_tags.py")
            return 1
        print("\nproject tag tables up to date")
    return 0


if __name__ == "__main__":
    sys.exit(main())
