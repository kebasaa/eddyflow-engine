#!/usr/bin/env python3
"""Build base_lastgas: a run whose last half-hour has no data for one gas.

What the end of an RP run writes partly comes from the run's last period:
which gases were present then decides, among other things, which gases the
PWB time-lag summary lists. A production pass split across workers once took
that from the parent's own last period instead, and a Yatir run split four ways
left a gas out of its summary that the serial run listed. Every other fixture
has every gas in every period, so nothing caught it. Here COS has no data in the
last half-hour; check_parallel.sh PAR_KIND=pr on this fixture must still be
byte-identical to -j 1.

What it builds, from CH-LAE 2025-06-01 00:00-03:00 (base_rec's range):
  * a copy of the raw files in data_lastgas/, with the COS_DRY column emptied
    in the file of the last half-hour (02:30);
  * base_lastgas.eddyflow: base_rec reading that copy.

Usage:  python gen_lastgas.py [--data-out DIR]
The raw files go to data_lastgas/ beside this script, which is gitignored.
The source is the Lagern repository's own data folder, read and never written.
"""

import argparse
import re
import shutil
from datetime import datetime, timedelta
from pathlib import Path

HERE = Path(__file__).resolve().parent
START = datetime(2025, 6, 1, 0, 0)
SLOTS = 8  # 00:00 to 03:30: base_rec's range and the file after it
EMPTIED = datetime(2025, 6, 1, 2, 30)
COLUMN = "COS_DRY"
PROTOTYPE = "CH-LAE_ec_preproc_10hz_{:%Y%m%d-%H%M}.csv"


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


def empty_column(src, dst, name):
    with open(src, newline="") as f:
        lines = f.read().splitlines(True)
    header = lines[0].rstrip("\r\n").split(",")
    col = header.index(name)
    out = [lines[0]]
    for line in lines[1:]:
        end = line[len(line.rstrip("\r\n")):]
        fields = line.rstrip("\r\n").split(",")
        if col < len(fields):
            fields[col] = ""
        out.append(",".join(fields) + end)
    with open(dst, "w", newline="") as f:
        f.writelines(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-out", default=None,
                    help="where to write the raw files "
                         "(default: data_lastgas/ beside this script)")
    args = ap.parse_args()

    template = (HERE / "base_rec.eddyflow").read_text(encoding="utf-8")
    src_dir = Path(ini_value(template, "data_path"))
    out_dir = Path(args.data_out) if args.data_out else HERE / "data_lastgas"
    out_dir.mkdir(parents=True, exist_ok=True)
    for old in out_dir.glob("CH-LAE_ec_preproc_10hz_*.csv"):
        old.unlink()

    for i in range(SLOTS):
        t = START + timedelta(minutes=30 * i)
        src = src_dir / PROTOTYPE.format(t)
        if not src.exists():
            raise SystemExit(f"missing source file {src}")
        if t == EMPTIED:
            empty_column(src, out_dir / src.name, COLUMN)
        else:
            shutil.copyfile(src, out_dir / src.name)
    print(f"{SLOTS} half-hours copied to {out_dir}, {COLUMN} emptied at {EMPTIED:%H:%M}")

    p = set_key(template, "data_path", out_dir.as_posix())
    (HERE / "base_lastgas.eddyflow").write_text(p, encoding="utf-8")
    print("wrote base_lastgas.eddyflow  (from base_rec.eddyflow)")


if __name__ == "__main__":
    main()
