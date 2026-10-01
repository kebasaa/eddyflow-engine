#!/usr/bin/env python3
"""Build the gappy fixtures: a pre-pass range with holes in its raw data.

Every other fixture reads an unbroken run of half-hourly files, so a parallel
pre-pass over it cuts into pieces that all carry the same work. Real records
do not look like that - the Yatir season thins out in June and stops for two
weeks in July - and the pre-pass now cuts its range by work, giving a period
with no raw file almost no weight. These fixtures are what exercises that.

What it builds, from 36 hours of CH-LAE starting 2025-06-01 00:00:

  * a copy of the raw files with holes punched in it (the files are simply not
    copied):

      00:00-03:00 on the 1st   a leading gap. The parent runs the first piece
                               and must read a raw file in it, so the first
                               piece has to stretch past this.
      05:30 on the 1st         one missing file.
      10:00-20:00 on the 1st   a long gap, a third of the range.
      02:00 on the 2nd         one missing file.

  * two projects over that data:

      base_gappy_tlag   time-lag optimisation and the planar fit - both
                        pre-passes, so both splits run over the holes - over
                        the whole 36 hours, processing only the last three.
      base_gappy_pwb    PWB cache generation, processing the whole 36 hours:
                        a PWB cache covers the processing range, not the
                        pre-pass one, so with three hours processed it would
                        pre-pass six periods and never split.

The gate is check_parallel.sh on each, at several -j: the split runs must be
byte-identical to -j 1.

Usage:  python gen_gappy.py [--data-out DIR]

Writes the projects beside this script and the raw files to data_gappy/
beside it, which is gitignored. Never beside the source: the source is the
Lagern repository's own data folder, read here and never written into.
"""

import argparse
import re
import shutil
from datetime import datetime, timedelta
from pathlib import Path

HERE = Path(__file__).resolve().parent

START = datetime(2025, 6, 1, 0, 0)
SLOTS = 72  # 36 hours of half-hours
PROTOTYPE = "CH-LAE_ec_preproc_10hz_{:%Y%m%d-%H%M}.csv"

#: Half-open [from, to) stretches with no raw file.
HOLES = [
    (datetime(2025, 6, 1, 0, 0), datetime(2025, 6, 1, 3, 0)),
    (datetime(2025, 6, 1, 5, 30), datetime(2025, 6, 1, 6, 0)),
    (datetime(2025, 6, 1, 10, 0), datetime(2025, 6, 1, 20, 0)),
    (datetime(2025, 6, 2, 2, 0), datetime(2025, 6, 2, 2, 30)),
]

#: The pre-pass covers the whole range; the time-lag fixture processes only
#: its last 3 hours.
PREPASS = (datetime(2025, 6, 1, 0, 0), datetime(2025, 6, 2, 12, 0))
PROCESS = (datetime(2025, 6, 2, 9, 0), datetime(2025, 6, 2, 12, 0))

VARIANTS = {
    "base_gappy_tlag": ("base_tlag_par.eddyflow", {"rot_meth": "3"}, PROCESS),
    "base_gappy_pwb": ("base_pwb_par.eddyflow", {}, PREPASS),
}


def ini_value(text, key):
    m = re.search(rf"^{re.escape(key)}=(.*)$", text, re.M)
    if not m:
        raise SystemExit(f"{key} not found")
    return m.group(1).strip()


def set_key(text, key, value):
    pat = re.compile(rf"^{re.escape(key)}=.*$", re.M)
    if not pat.search(text):
        raise SystemExit(f"{key} not found")
    return pat.sub(f"{key}={value}", text)


def set_window(text, prefix, span):
    a, b = span
    text = set_key(text, f"{prefix}_start_date", f"{a:%Y-%m-%d}")
    text = set_key(text, f"{prefix}_start_time", f"{a:%H:%M}")
    text = set_key(text, f"{prefix}_end_date", f"{b:%Y-%m-%d}")
    return set_key(text, f"{prefix}_end_time", f"{b:%H:%M}")


def in_hole(t):
    return any(a <= t < b for a, b in HOLES)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-out", default=None,
                    help="where to write the raw files "
                         "(default: data_gappy/ beside this script)")
    args = ap.parse_args()

    template = (HERE / VARIANTS["base_gappy_tlag"][0]).read_text(encoding="utf-8")
    src_dir = Path(ini_value(template, "data_path"))
    out_dir = Path(args.data_out) if args.data_out else HERE / "data_gappy"
    out_dir.mkdir(parents=True, exist_ok=True)
    for old in out_dir.glob("CH-LAE_ec_preproc_10hz_*.csv"):
        old.unlink()

    kept = 0
    for i in range(SLOTS):
        t = START + timedelta(minutes=30 * i)
        if in_hole(t):
            continue
        src = src_dir / PROTOTYPE.format(t)
        if not src.exists():
            raise SystemExit(f"missing source file {src}")
        shutil.copyfile(src, out_dir / src.name)
        kept += 1
    print(f"{kept} of {SLOTS} half-hours copied to {out_dir}")

    for name, (source, keys, process) in VARIANTS.items():
        p = (HERE / source).read_text(encoding="utf-8")
        p = set_key(p, "data_path", out_dir.as_posix())
        p = set_window(p, "pr", process)
        p = set_window(p, "to", PREPASS)
        p = set_window(p, "pf", PREPASS)
        for k, v in keys.items():
            p = set_key(p, k, v)
        (HERE / f"{name}.eddyflow").write_text(p, encoding="utf-8")
        print(f"wrote {name}.eddyflow  (from {source}{', ' + str(keys) if keys else ''})")


if __name__ == "__main__":
    main()
