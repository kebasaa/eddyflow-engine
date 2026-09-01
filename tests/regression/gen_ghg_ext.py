#!/usr/bin/env python3
"""Build the extended-.ghg fixture from the committed LI-COR archives.

Writes data_ghg_ext/ beside data_ghg/, holding the SAME two archives with their
embedded .metadata rewritten into the extended format documented as "Extended
GHG Files":

  * ghg_format_version=2.0
  * the LI-7500A block declares itself as a generic_open_path stand-in with real
    geometry - what an EddyPro that did not know the instrument would need -
    and states its true identity in instr_<k>_ef_model
  * every col_N_instrument that named the LI-7500A is repointed at the stand-in,
    because that is the only name EddyPro can resolve

The point of the fixture is an EQUALITY: base_ghg_ext must produce byte-identical
output to base_ghg_licor. EddyFlow reads ef_model, so it is processing the same
LI-7500A either way, and any difference means the override did not take - which
is the failure this cannot otherwise be caught doing, since a stand-in that is
silently believed still produces plausible fluxes.

Generated rather than committed: the archives are 3.5 MB of binary and the
interesting part is the ~10 lines of metadata that differ, which belong in a
diffable script rather than inside a second zip.

Usage:  python gen_ghg_ext.py     (needs 7-Zip on PATH, same as the fixtures)
"""

import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
SRC = HERE / "data_ghg"
DST = HERE / "data_ghg_ext"

#: The instrument to stand in for, and what to declare it as. Only the analyser
#: is swapped: the LI-7700 and the Metek are models EddyPro knows, and swapping
#: instruments EddyPro can already name would test nothing.
REAL = "li7500a_1"
STANDIN = "generic_open_path_1"
#: Required for any generic model - EddyPro refuses one whose
#: hpath * vpath * tau is zero.
#:
#: These are the LI-7500A's REAL numbers, taken from the table
#: RetrieveSensorParams assigns by model (0.127 m vertical, 0.0095 m
#: horizontal, 0.1 s), converted to the centimetres a .metadata states. Stating
#: the truth is the point rather than a nicety: a stand-in carrying plausible
#: but wrong geometry - 12.5 cm horizontal instead of 0.95 - moves LE by 0.2 %
#: in both engines, which is the whole of the measured "cost of standing in".
#: With these, EddyFlow reading the stand-in and EddyFlow reading the LI-7500A
#: agree to the digit.
#:
#: EddyFlow does not read them here anyway - ef_model puts the model back to
#: li7500a_1 before the block that would - so they exist for EddyPro's sake.
GEOMETRY = [("hpath_length", "0.95"), ("vpath_length", "12.7"), ("tau", "0.1")]


def sevenzip():
    for exe in ("7z", "7za"):
        if shutil.which(exe):
            return exe
    sys.exit("7-Zip not on PATH; the GHG fixtures need it too")


def rewrite(text):
    """Classic .metadata -> extended .metadata. Byte-preserving elsewhere."""
    out, block, did_geom = [], None, False
    for line in text.split("\n"):
        key = line.split("=", 1)[0].strip()

        #: The version key goes at the end of [FileDescription], which is the
        #: section that describes the file rather than the site.
        if line.strip() == "[FileDescription]":
            out.append(line)
            out.append("ghg_format_version=2.0")
            continue

        if key.endswith("_model") and line.strip().endswith("=" + REAL):
            block = key[: -len("_model")]          # e.g. 'instr_2'
            out.append(f"{block}_model={STANDIN}")
            continue

        #: Geometry and the true identity are appended to the block being
        #: stood in for, after its last original key. sw_version is the key
        #: that follows model in every file LI-COR writes.
        if block and key == f"{block}_sw_version" and not did_geom:
            out.append(line)
            out += [f"{block}_{k}={v}" for k, v in GEOMETRY]
            out.append(f"{block}_ef_model={REAL}")
            did_geom = True
            continue

        if key.endswith("_instrument") and line.strip().endswith("=" + REAL):
            out.append(f"{key}={STANDIN}")
            continue

        out.append(line)

    if not did_geom:
        sys.exit(f"no instrument block declares model={REAL}; nothing to stand in for")
    return "\n".join(out)


def main():
    z = sevenzip()
    archives = sorted(SRC.glob("*.ghg"))
    if not archives:
        sys.exit(f"no archives in {SRC}")

    DST.mkdir(exist_ok=True)
    for ghg in archives:
        target = DST / ghg.name
        shutil.copy2(ghg, target)
        with tempfile.TemporaryDirectory() as tmp:
            #: The archive holds TWO .metadata files - the site's, and the
            #: biomet one beside it. Extracted by exact name rather than by
            #: glob, because <stem>-biomet.metadata sorts first and describes
            #: no instruments at all: taking it produces a rewrite that
            #: changes nothing and a fixture that proves nothing.
            want = ghg.stem + ".metadata"
            subprocess.run([z, "e", str(target), "-o" + tmp, want, "-y"],
                           check=True, stdout=subprocess.DEVNULL)
            md = Path(tmp) / want
            if not md.is_file():
                sys.exit(f"{ghg.name} holds no {want}")
            #: cp1252 and LF, read and written as BYTES. Python's text mode
            #: would translate the line endings on Windows and every line in
            #: the file would differ from the archive it came out of.
            md.write_bytes(
                rewrite(md.read_bytes().decode("cp1252")).encode("cp1252"))
            subprocess.run([z, "u", str(target), md.name, "-y"],
                           check=True, cwd=tmp, stdout=subprocess.DEVNULL)
        print(f"  {target.name}")
    print(f"{len(archives)} extended archive(s) in {DST}")


if __name__ == "__main__":
    main()
