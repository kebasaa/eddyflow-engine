#!/usr/bin/env python3
"""Build the slow-instrument fixture: an analyser sampling under the row rate.

The motivating case for the per-instrument allowance. Every other fixture has
one acquisition frequency, so `expected == nrow` everywhere and the per-column
completeness test reduces exactly to the project-wide one - which is what makes
those fixtures the "nothing moved" gate and this one the only thing that
exercises the feature at all.

What it builds, from the same three hours of CH-LAE the other fixtures use:

  * a copy of the six raw files in which the MIRO's six columns carry -9999 on
    nine rows out of ten. CleanUpE2Set maps anything below -300 to the error
    code, so those rows are error rows between the analyser's real samples -
    exactly the shape a 1 Hz instrument writes into a 10 Hz file. The value has
    to be numeric: ImportAscii reads a whole row list-directed, so a text token
    in one field would drop the entire row, sonic and all.

  * three projects over that data, each isolating one half of the mechanism:

      base_slow_naive   nothing declared. The MIRO's columns are 90 % error
                        against the row grid, over the 40 % global allowance,
                        so they are dropped and their fluxes are -9999. This
                        is what the engine did before the feature, reproduced
                        with the current binary.
      base_slow_lack    instr_3_max_lack=95, no rate. Proves the project key
                        is read and reaches the drop test: 90 % missing is now
                        under the instrument's own allowance, so the columns
                        survive without anything knowing the MIRO is slow.
      base_slow         instr_3_ac_freq=1.0 with a tight instr_3_max_lack=10.
                        Proves the rate is read: measured against what a 1 Hz
                        instrument should have produced, the column is
                        complete, so a 10 % allowance is ample.

Usage:  python gen_slow.py [--data-out DIR]

Writes the projects and metadata beside this script, and the decimated raw
files to a sibling of the source data directory.
"""

import argparse
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent

#: The instrument block the MIRO occupies in base_site.metadata, and the rate
#: it is declared to run at in the `base_slow` variant. 1 Hz against the
#: station's 10 is the ratio the decimation below writes.
MIRO_SLOT = 3
MIRO_AC_FREQ = 1.0
#: 1-based raw-file columns belonging to that instrument: its four gases and
#: its cell temperature and pressure. Its flow rate and status word are left
#: alone deliberately - they are on the same instrument and stay complete, so
#: the run also shows that a column with more data than expected is not
#: mistaken for a corrupt one.
MIRO_COLS = [7, 8, 9, 10, 11, 12]
#: Keep one row in ten. Anything below -300 is error to CleanUpE2Set.
KEEP_EVERY = 10
FILL = "-9999"

HEADER_ROWS = 4
SOURCE_PROJECT = "base_n_gas.eddyflow"
SOURCE_METADATA = "base_site.metadata"

#: The three-hour subset base_n_gas processes, one file per half hour.
FILES = [
    f"CH-LAE_ec-cos_20250601-{h:02d}{m:02d}.csv"
    for h in range(3)
    for m in (0, 30)
]


def ini_value(text, key):
    m = re.search(rf"^{re.escape(key)}=(.*)$", text, re.M)
    if not m:
        raise SystemExit(f"{key} not found")
    return m.group(1).strip()


def set_key(text, key, value, after=None):
    """Replace a key if present, otherwise insert it.

    `after` names the key to insert below, which is how a new key lands in the
    right INI section: the engine parses each table against one section prefix,
    so a per-instrument allowance appended at the end of the file would be read
    only if the last section happened to be a RawProcess one.
    """
    pat = re.compile(rf"^{re.escape(key)}=.*$", re.M)
    if pat.search(text):
        return pat.sub(f"{key}={value}", text)
    if after is None:
        return text.rstrip("\n") + f"\n{key}={value}\n"
    anchor = re.search(rf"^{re.escape(after)}=.*$", text, re.M)
    if not anchor:
        raise SystemExit(f"cannot place {key}: no {after} to insert after")
    return text[:anchor.end()] + f"\n{key}={value}" + text[anchor.end():]


def decimate(src, dst):
    """Copy one raw file, blanking the MIRO's columns on the rows it skips."""
    idx = [c - 1 for c in MIRO_COLS]
    kept = 0
    with open(src, encoding="utf-8", errors="replace") as fi, \
         open(dst, "w", encoding="utf-8", newline="") as fo:
        for n, line in enumerate(fi):
            if n < HEADER_ROWS:
                fo.write(line)
                continue
            if (n - HEADER_ROWS) % KEEP_EVERY == 0:
                fo.write(line)
                kept += 1
                continue
            f = line.rstrip("\r\n").split(",")
            for i in idx:
                if i < len(f):
                    f[i] = FILL
            fo.write(",".join(f) + "\n")
    return kept


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-out", default=None,
                    help="where to write the decimated raw files "
                         "(default: a '_slow' sibling of the source data dir)")
    args = ap.parse_args()

    project = (HERE / SOURCE_PROJECT).read_text(encoding="utf-8")
    metadata = (HERE / SOURCE_METADATA).read_text(encoding="utf-8")

    src_dir = Path(ini_value(project, "data_path"))
    out_dir = Path(args.data_out) if args.data_out else \
        src_dir.with_name(src_dir.name + "_slow")
    out_dir.mkdir(parents=True, exist_ok=True)

    for name in FILES:
        src = src_dir / name
        if not src.exists():
            raise SystemExit(f"missing source file {src}")
        kept = decimate(src, out_dir / name)
        print(f"{name}  {kept} rows kept of every {KEEP_EVERY}th")

    #> Two metadata files: one that says nothing about the MIRO's rate, which
    #> is every metadata file written before the feature, and one that states
    #> it. The engine reads instr_<K>_ac_freq at a slot that was reserved and
    #> unread, so the first is not merely an older file - it is the same file.
    slow_meta = set_key(metadata, f"instr_{MIRO_SLOT}_ac_freq",
                        f"{MIRO_AC_FREQ:.3f}",
                        after=f"instr_{MIRO_SLOT}_ko")
    (HERE / "base_slow.metadata").write_text(slow_meta, encoding="utf-8")
    (HERE / "base_slow_naive.metadata").write_text(metadata, encoding="utf-8")

    #> The same rate, declared as integrating over the interval rather than
    #> reporting an instant. Only the w pairing differs, so only the COSPECTRA
    #> may move - the gas's own spectrum is built from the same samples either
    #> way, and if it moves the flag is reaching something it should not.
    integr_meta = set_key(slow_meta, f"instr_{MIRO_SLOT}_integrates", "1",
                          after=f"instr_{MIRO_SLOT}_ac_freq")
    (HERE / "base_slow_integr.metadata").write_text(integr_meta, encoding="utf-8")

    variants = {
        "base_slow_naive": ("base_slow_naive.metadata", {}),
        "base_slow_lack": ("base_slow_naive.metadata",
                           {f"instr_{MIRO_SLOT}_max_lack": "95"}),
        "base_slow": ("base_slow.metadata",
                      {f"instr_{MIRO_SLOT}_max_lack": "10"}),
        "base_slow_integr": ("base_slow_integr.metadata",
                             {f"instr_{MIRO_SLOT}_max_lack": "10"}),
    }
    for name, (meta, keys) in variants.items():
        p = set_key(project, "data_path", out_dir.as_posix())
        p = set_key(p, "proj_file", (HERE / meta).as_posix())
        for k, v in keys.items():
            p = set_key(p, k, v, after="max_lack")
        (HERE / f"{name}.eddyflow").write_text(p, encoding="utf-8")
        print(f"wrote {name}.eddyflow  ({meta}, {keys or 'nothing declared'})")


if __name__ == "__main__":
    main()
