#!/usr/bin/env python3
"""Regenerate the .metadata tag tables (ANTags / ACTags) in m_common_global_var.f90.

These tables are positional: a tag's identity IS its array index, so the blocks
must stay contiguous and the group origins must match the stride arithmetic in
read_metadata_file.f90. Widening the instrument block therefore shifts every
column entry after it - about 1,500 hand-numbered lines. Hand-maintaining that
is how the err_label collision and the missing out_full_sp_gas4 key got in, so
it is generated instead.

Layout produced (1-based, contiguous):

    ANTags : header(9) | instruments(N x 15) | columns(100 x 8)
    ACTags : header(24) | instruments(N x 8) | data_label(1) | columns(100 x 7)
             | instr_<N+1>_manufacturer   <- overflow sentinel

The sentinel lets the engine notice that the file describes more instruments
than it can hold, instead of silently dropping them.

Usage:  python gen_metadata_tags.py [--instruments 8] [--check]

--check exits non-zero if the file is not what the generator would produce,
which is what makes this safe to re-run in CI.
"""

import argparse
import re
import sys
from pathlib import Path

SRC = Path(__file__).resolve().parents[1] / "src" / "src_common" / "m_common_global_var.f90"

NUM_COLS = 100
TAG_RE = re.compile(r"(A[NC]Tags)\((\d+)\)%Label\s*/\s*'([^']*)'\s*/")

# Markers delimiting the two generated regions. The generator rewrites strictly
# between them, so hand-written code above and below is preserved.
BEGIN_AN = "    !> BEGIN GENERATED ANTags - edit gen_metadata_tags.py, not this block"
END_AN = "    !> END GENERATED ANTags"
BEGIN_AC = "    !> BEGIN GENERATED ACTags - edit gen_metadata_tags.py, not this block"
END_AC = "    !> END GENERATED ACTags"


def parse_existing(text):
    """index -> label, per table, ignoring commented-out lines."""
    tables = {"ANTags": {}, "ACTags": {}}
    for line in text.splitlines():
        if line.lstrip().startswith("!"):
            continue
        for m in TAG_RE.finditer(line):
            tables[m.group(1)][int(m.group(2))] = m.group(3)
    return tables


def suffixes(labels, prefix):
    """Ordered suffixes of e.g. 'instr_1_*' / 'col_1_*', in index order."""
    out = []
    for idx in sorted(labels):
        lab = labels[idx]
        if lab.startswith(prefix):
            out.append(lab[len(prefix):])
    return out


# Fortran allows at most 255 continuation lines in one statement. gfortran only
# warns, but other compilers reject it, so the block is split into several
# `data` statements rather than one enormous continued one.
MAX_CONT = 200


def emit(table, entries, width):
    """Render `data` lines, one assignment per line, in <=MAX_CONT chunks."""
    lines = []
    for start in range(0, len(entries), MAX_CONT):
        chunk = entries[start:start + MAX_CONT]
        for i, (idx, label) in enumerate(chunk):
            lead = (f"    data {table}({idx})%Label" if i == 0
                    else f"         {table}({idx})%Label")
            cont = " /" if i == len(chunk) - 1 else " / &"
            lines.append(f"{lead:<{width}} / '{label}'{cont}")
    return lines


def build(n_instr, tables):
    an, ac = tables["ANTags"], tables["ACTags"]

    # Group boundaries are derived from the CURRENT file, so the generator does
    # not hard-code how many header tags there are.
    an_instr = suffixes(an, "instr_1_")
    an_col = suffixes(an, "col_1_")
    ac_instr = suffixes(ac, "instr_1_")
    ac_col = suffixes(ac, "col_1_")

    first_an_instr = min(i for i, l in an.items() if l.startswith("instr_1_"))
    first_ac_instr = min(i for i, l in ac.items() if l.startswith("instr_1_"))
    first_ac_col = min(i for i, l in ac.items() if l.startswith("col_1_"))

    # How many instruments the file currently describes, EXCLUDING the overflow
    # sentinel (which is deliberately one past the limit and has only the
    # manufacturer key). Deriving this rather than assuming the previous count
    # is what makes the generator idempotent: re-running it must not mistake
    # already-widened instrument blocks for content that belongs after them.
    instr_nums = sorted({
        int(re.match(r"instr_(\d+)_", l).group(1))
        for l in ac.values() if re.match(r"instr_(\d+)_", l)
    })
    full = [k for k in instr_nums
            if sum(1 for l in ac.values() if l.startswith(f"instr_{k}_")) == len(ac_instr)]
    cur_instr = max(full) if full else 0

    an_header = [an[i] for i in range(1, first_an_instr)]
    ac_header = [ac[i] for i in range(1, first_ac_instr)]
    # Whatever sits between the instrument block and the column block (today
    # just 'data_label'); kept so its position stays correct after widening.
    old_ac_instr_end = first_ac_instr + cur_instr * len(ac_instr)
    ac_middle = [ac[i] for i in range(old_ac_instr_end, first_ac_col)]

    an_entries, idx = [], 1
    for lab in an_header:
        an_entries.append((idx, lab)); idx += 1
    for k in range(1, n_instr + 1):
        for s in an_instr:
            an_entries.append((idx, f"instr_{k}_{s}")); idx += 1
    an_col_origin = idx
    for c in range(1, NUM_COLS + 1):
        for s in an_col:
            an_entries.append((idx, f"col_{c}_{s}")); idx += 1
    nan = idx - 1

    ac_entries, idx = [], 1
    for lab in ac_header:
        ac_entries.append((idx, lab)); idx += 1
    for k in range(1, n_instr + 1):
        for s in ac_instr:
            ac_entries.append((idx, f"instr_{k}_{s}")); idx += 1
    for lab in ac_middle:
        ac_entries.append((idx, lab)); idx += 1
    ac_col_origin = idx
    for c in range(1, NUM_COLS + 1):
        for s in ac_col:
            ac_entries.append((idx, f"col_{c}_{s}")); idx += 1
    # Overflow sentinel: present in the file only when there are more
    # instruments than we can store.
    sentinel_idx = idx
    ac_entries.append((idx, f"instr_{n_instr + 1}_manufacturer")); idx += 1
    nac = idx - 1

    return {
        "an_entries": an_entries, "ac_entries": ac_entries,
        "nan": nan, "nac": nac,
        "an_instr_leap": len(an_instr), "ac_instr_leap": len(ac_instr),
        "an_col_leap": len(an_col), "ac_col_leap": len(ac_col),
        "an_instr_origin": first_an_instr, "ac_instr_origin": first_ac_instr,
        "an_col_origin": an_col_origin, "ac_col_origin": ac_col_origin,
        "ac_middle_origin": first_ac_instr + n_instr * len(ac_instr),
        "sentinel_idx": sentinel_idx,
    }


def splice(text, begin, end, body):
    """Replace between markers, inserting them if absent is not supported."""
    if begin not in text or end not in text:
        raise SystemExit(
            f"marker missing: {begin!r}. Add the BEGIN/END marker pair around the "
            "existing data block once, by hand, then re-run."
        )
    head, rest = text.split(begin, 1)
    _, tail = rest.split(end, 1)
    return head + begin + "\n" + "\n".join(body) + "\n" + end + tail


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--instruments", type=int, default=8)
    ap.add_argument("--check", action="store_true")
    args = ap.parse_args()

    text = SRC.read_text(encoding="utf-8", errors="surrogateescape")
    info = build(args.instruments, parse_existing(text))

    new = splice(text, BEGIN_AN, END_AN, emit("ANTags", info["an_entries"], 24))
    new = splice(new, BEGIN_AC, END_AC, emit("ACTags", info["ac_entries"], 24))

    new = re.sub(r"integer, parameter :: Nan = \d+",
                 f"integer, parameter :: Nan = {info['nan']}", new)
    new = re.sub(r"integer, parameter :: Nac = \d+",
                 f"integer, parameter :: Nac = {info['nac']}", new)

    if args.check:
        if new != text:
            print("m_common_global_var.f90 is out of date; re-run gen_metadata_tags.py")
            return 1
        print("metadata tag tables up to date")
        return 0

    SRC.write_text(new, encoding="utf-8", errors="surrogateescape")
    print(f"instruments      : {args.instruments}")
    print(f"Nan              : {info['nan']}   (columns start at {info['an_col_origin']})")
    print(f"Nac              : {info['nac']}   (columns start at {info['ac_col_origin']})")
    print(f"overflow sentinel: ACTags({info['sentinel_idx']})")
    print()
    print("read_metadata_file.f90 must use:")
    print(f"    leap_an_instr = {info['an_instr_leap']} ; init_an_instr = "
          f"{info['an_instr_origin']} - leap_an_instr")
    print(f"    leap_ac_instr = {info['ac_instr_leap']} ; init_ac_instr = "
          f"{info['ac_instr_origin']} - leap_ac_instr")
    print(f"    leap_an_col   = {info['an_col_leap']} ; init_an_col   = "
          f"{info['an_col_origin']} - leap_an_col")
    print(f"    leap_ac_col   = {info['ac_col_leap']} ; init_ac_col   = "
          f"{info['ac_col_origin']} - leap_ac_col")
    print(f"    data_label is ACTags({info['ac_middle_origin']})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
