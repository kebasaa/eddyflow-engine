#!/usr/bin/env python3
"""Build the GHG fixture: the same three hours, as LI-COR archives.

`file_type=0` - the compressed LI-COR GHG archive - had no fixture at all
before this. That is the one raw format whose per-file cost is not parsing:
UnZipArchive spends five shell invocations and a full decompress before a
single record is read, and ReadLicorGhgArchive a sixth to clean up. Measured
here, a GHG file costs about 350 ms more than the same data as CSV.

It is also the only path that rewrites `Metadata%ac_freq`, `NumCol` and
`FileInterpreter` per file, because each archive carries its own metadata -
so it is the one place where those globals move under the period loop. None of
that was exercised by anything.

What this builds, from the same raw files the other fixtures read:

  * one `.ghg` per half-hour, each a zip of `<name>.data` (the raw record file,
    unchanged) and `<name>.metadata` (a copy of base_site.metadata). That is
    exactly what UnZipArchive looks for - it searches the extracted files by
    extension, not by name.

  * `base_ghg.eddyflow`, the same project as `base_tlag_opt` but pointed at
    them with `file_type=0` and a `.ghg` prototype.

The gate is that it must produce **the same fluxes as the CSV fixture over the
same window**, because it is the same data. That is checked directly: run
base_ghg and base_tlag_opt over one window and diff the FLUXNET files with the
filename column normalised. They match today.

Needs 7-Zip on PATH, which is what the engine shells out to. sweep.sh skips
the fixture when it is missing rather than failing, because a machine without
it can still run every other fixture.

Usage:
    python gen_ghg.py [--data DIR] [--out DIR] [--hours N]
"""

import argparse
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
DEFAULT_DATA = Path(r"C:\Users\jonmuell\Documents\_data\CH-LAE COS\Data_202506")
DEFAULT_OUT = Path(r"C:\Users\jonmuell\Documents\_data\CH-LAE COS\Data_202506_ghg")

#: Where 7-Zip is looked for when it is not already on PATH. The portable
#: distribution ships one, and that is the same binary the engine would use.
FALLBACKS = [
    Path(r"C:\Users\jonmuell\Documents\GitHub\eddyflow-portable\bin\7z.exe"),
    Path(r"C:\Program Files\7-Zip\7z.exe"),
]


def find_7z():
    found = shutil.which("7z") or shutil.which("7za")
    if found:
        return found
    for candidate in FALLBACKS:
        if candidate.exists():
            return str(candidate)
    sys.exit("7-Zip not found. Put 7z.exe on PATH - the engine shells out to it.")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", type=Path, default=DEFAULT_DATA)
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--hours", type=float, default=3.0,
                    help="how much of the dataset to convert")
    args = ap.parse_args()

    sevenzip = find_7z()
    metadata = HERE / "base_site.metadata"
    if not metadata.exists():
        sys.exit("missing %s" % metadata)
    if not args.data.is_dir():
        sys.exit("raw data not found at %s - it lives outside the repo" % args.data)

    files = sorted(args.data.glob("*.csv"))[: int(args.hours * 2)]
    if not files:
        sys.exit("no .csv files under %s" % args.data)

    args.out.mkdir(parents=True, exist_ok=True)
    for stale in args.out.glob("*.ghg"):
        stale.unlink()

    for src in files:
        stem = src.stem
        with tempfile.TemporaryDirectory() as stage:
            stage = Path(stage)
            #> The names inside matter only by extension, but keeping the
            #> archive's own stem is what LI-COR does and what a reader
            #> debugging this would expect to see.
            shutil.copy2(src, stage / (stem + ".data"))
            shutil.copy2(metadata, stage / (stem + ".metadata"))
            archive = args.out / (stem + ".ghg")
            #> -mx=1 because these are regenerated often and read once each;
            #> the engine does not care how they were compressed.
            subprocess.run(
                [sevenzip, "a", "-tzip", "-mx=1", str(archive),
                 str(stage / (stem + ".data")), str(stage / (stem + ".metadata"))],
                check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        print("  %s" % archive.name)

    #> The project, derived from the one whose numbers it has to match.
    base = (HERE / "base_tlag_opt.eddyflow").read_text(encoding="utf-8",
                                                       errors="replace")
    out_lines = []
    for line in base.splitlines(keepends=True):
        key = line.split("=", 1)[0]
        if key == "file_type":
            line = "file_type=0" + line[len(line.rstrip()):]
        elif key == "file_prototype":
            line = ("file_prototype=CH-LAE_ec-cos_yyyymmdd-HHMM.ghg"
                    + line[len(line.rstrip()):])
        elif key == "data_path":
            line = ("data_path=" + args.out.as_posix()
                    + line[len(line.rstrip()):])
        out_lines.append(line)
    (HERE / "base_ghg.eddyflow").write_text("".join(out_lines), encoding="utf-8")

    print("\n%d archives in %s" % (len(files), args.out))
    print("base_ghg.eddyflow written")


if __name__ == "__main__":
    main()
