#!/usr/bin/env python3
"""Build the mixed-acquisition-frequency GHG fixtures from the LI-COR archives.

Each .ghg carries its own .metadata, and a site can change rate part way through
a project: an analyser swapped for a faster one, a logger reconfigured. RP used
to size its period buffer once, from the first file. After a change to a higher
rate, only the first half of each period fitted and every period was dropped for
"not enough samples"; the pre-passes silently used the half that fitted.

The two committed archives (data_ghg/, 10 Hz, 01:30 and 02:00) are stretched to
six half-hours, 01:00 to 03:30, alternating between the two and renamed, with
every name inside them renamed to match:

data_ghg_mixed/ - 01:00 to 02:00 are DECIMATED to 5 Hz (every second record,
    acquisition_frequency=5.0); 02:30 to 03:30 stay at 10 Hz. That is a change
    UP in rate, which is the case that broke.
      base_ghg_mixed     30-min periods. All six are processed, each at its own
                         rate: 9000 records, then 18000.
      base_ghg_mixed_60  60-min periods. 02:00-03:00 holds a 5 Hz and a 10 Hz
                         file, so it is skipped with Warning(116); 01:00-02:00
                         (5 Hz) and 03:00-04:00 (10 Hz) are processed.

data_ghg_mixed_instr/ - all six at 10 Hz, but from 02:30 on they state that the
    LI-7700 runs at 1 Hz (instr_3_ac_freq=1.0). The row rate never changes, so
    only the per-instrument comparison can see it.
      base_ghg_mixed_instr  60-min periods; 02:00-03:00 is skipped with
                            Warning(116), the other two are processed.

The 10 Hz periods of base_ghg_mixed are the committed archives unchanged but for
their names, so their fluxes must equal base_ghg_licor's to the last digit - the
direct check that a period after a change of rate is read whole.

The renamed copies keep their original timestamps INSIDE the data (the Date and
Time columns, the "Timestamp:" header line and the embedded biomet). The engine
places a file by its name, so that is harmless for fluxes; the embedded biomet
of a period simply finds fewer records in its window.

Generated rather than committed, like data_ghg_ext/.

Usage:  python gen_ghg_mixed.py     (needs 7-Zip on PATH, same as the fixtures)
"""

import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
SRC = HERE / "data_ghg"
BASE = HERE / "base_ghg_licor.eddyflow"

#: (source stem time, target stem time, decimate?, slow LI-7700?)
MIXED = [
    ("013000", "010000", True, False),
    ("020000", "013000", True, False),
    ("013000", "020000", True, False),
    ("020000", "023000", False, False),
    ("013000", "030000", False, False),
    ("020000", "033000", False, False),
]
MIXED_INSTR = [
    ("013000", "010000", False, False),
    ("020000", "013000", False, False),
    ("013000", "020000", False, False),
    ("020000", "023000", False, True),
    ("013000", "030000", False, True),
    ("020000", "033000", False, True),
]

#: fixture: (data directory, averaging interval, extra project keys)
#:
#: base_ghg_mixed_sa lowers sa_min_smpl to 2 so that the six periods - three
#: per rate - are enough for a spectral assessment at each rate. The fits are
#: statistically meaningless; what the fixture exercises is the machinery: one
#: assessment pass per rate, one assessment file with a column set per rate,
#: and every period corrected with its own rate's result.
FIXTURES = {
    "base_ghg_mixed": ("data_ghg_mixed", 30, {}),
    "base_ghg_mixed_60": ("data_ghg_mixed", 60, {}),
    "base_ghg_mixed_instr": ("data_ghg_mixed_instr", 60, {}),
    "base_ghg_mixed_sa": ("data_ghg_mixed", 30, {"sa_min_smpl": "2"}),
}

HEADER_ROWS = 8


def sevenzip():
    for exe in ("7z", "7za"):
        if shutil.which(exe):
            return exe
    sys.exit("7-Zip not on PATH; the GHG fixtures need it too")


def decimate(raw):
    """Keep the header and every second record. Bytes, so the line endings
    and the encoding come out exactly as they went in."""
    lines = raw.splitlines(keepends=True)
    return b"".join(lines[:HEADER_ROWS] + lines[HEADER_ROWS::2])


def rewrite_metadata(text, decimated, slow_7700):
    if decimated:
        text, n = re.subn(r"(?m)^acquisition_frequency=.*$",
                          "acquisition_frequency=5.0", text)
        assert n == 1, "no acquisition_frequency line"
    if slow_7700:
        text, n = re.subn(r"(?m)^(instr_3_model=.*)$",
                          r"\1\ninstr_3_ac_freq=1.0", text)
        assert n == 1, "no instr_3_model line"
    return text


def build(name, plan, z):
    dst = HERE / name
    if dst.exists():
        shutil.rmtree(dst)
    dst.mkdir()
    for src_time, dst_time, decimated, slow_7700 in plan:
        src = next(SRC.glob(f"*T{src_time}_*.ghg"))
        src_stem = src.stem
        dst_stem = src_stem.replace(f"T{src_time}_", f"T{dst_time}_")
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            #: Unpack everything, rename, repack: a renamed archive whose
            #: members kept the old stem would still be found by extension,
            #: but would not look like anything LI-COR writes.
            subprocess.run([z, "x", str(src), "-o" + str(tmp), "-y"],
                           check=True, stdout=subprocess.DEVNULL)
            for f in list(tmp.rglob("*")):
                if f.is_file() and src_stem in f.name:
                    f.rename(f.with_name(f.name.replace(src_stem, dst_stem)))
            data = tmp / (dst_stem + ".data")
            meta = tmp / (dst_stem + ".metadata")
            if decimated:
                data.write_bytes(decimate(data.read_bytes()))
            #: cp1252 and bytes, as gen_ghg_ext.py does, for the same reason.
            meta.write_bytes(rewrite_metadata(
                meta.read_bytes().decode("cp1252"), decimated, slow_7700
            ).encode("cp1252"))
            target = dst / (dst_stem + ".ghg")
            members = [p.name for p in tmp.iterdir()]
            subprocess.run([z, "a", "-tzip", "-mx=1", str(target)] + members,
                           check=True, cwd=tmp, stdout=subprocess.DEVNULL)
        tags = []
        if decimated:
            tags.append("5 Hz")
        if slow_7700:
            tags.append("LI-7700 at 1 Hz")
        print(f"  {name}/{target.name}" + (f"  ({', '.join(tags)})" if tags else ""))


def write_fixture(fixture, data_dir, avrg_len, extra):
    text = BASE.read_text(encoding="utf-8", errors="replace")
    data_path = (HERE / data_dir).as_posix()
    text, n = re.subn(r"(?m)^data_path=.*$", "data_path=" + data_path, text)
    assert n == 1
    text, n = re.subn(r"(?m)^avrg_len=.*$", f"avrg_len={avrg_len}", text)
    assert n == 1
    for key, value in extra.items():
        text, n = re.subn(r"(?m)^%s=.*$" % re.escape(key), f"{key}={value}", text)
        assert n == 1, key
    (HERE / (fixture + ".eddyflow")).write_text(text, encoding="utf-8")
    print(f"  {fixture}.eddyflow  ({data_dir}, {avrg_len}-min periods)")


def main():
    z = sevenzip()
    if not BASE.exists():
        sys.exit(f"missing {BASE}")
    build("data_ghg_mixed", MIXED, z)
    build("data_ghg_mixed_instr", MIXED_INSTR, z)
    for fixture, (data_dir, avrg_len, extra) in FIXTURES.items():
        write_fixture(fixture, data_dir, avrg_len, extra)


if __name__ == "__main__":
    main()
