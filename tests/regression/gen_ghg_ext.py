#!/usr/bin/env python3
"""Build the extended-.ghg fixtures from the committed LI-COR archives.

Writes two directories beside data_ghg/, each holding the SAME two archives with
their embedded .metadata rewritten into the extended format documented as
"Extended GHG Files". Both add ghg_format_version=2.0, declare an instrument the
way EddyPro needs to see it, and state what it really is in
instr_<k>_ef_model - and both repoint every col_N_instrument that named that
instrument, because the name in the columns is the only one EddyPro can resolve.

data_ghg_ext/ - the STAND-IN case. The LI-7500A declares itself a
    generic_open_path carrying the analyser's real geometry: what an archive
    describing an instrument EddyPro cannot name has to look like.

data_ghg_campbell/ - the RENAME case, which the LI-COR archives cannot show on
    their own because every instrument in them is spelt the same on both sides.
    The Metek is swapped for a Campbell CSAT3B, whose EddyPro spelling is
    `csat3b` and whose EddyFlow spelling is `csi_csat3b`.

    This is the case that was broken. Nothing canonicalises col_N_instrument, so
    an instrument written as csat3b_1 and canonicalised to csi_csat3b_1 left
    every column pointing at a name the instrument list no longer had: no column
    bound to the sonic at all, and the run died with "exactly one selected u, v,
    w and one selected ts or sos are required", which does not mention
    instruments.

Generated rather than committed: the archives are 3.5 MB of binary and the
interesting part is the handful of metadata lines that differ, which belong in a
diffable script rather than inside another zip.

Usage:  python gen_ghg_ext.py     (needs 7-Zip on PATH, same as the fixtures)
"""

import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
SRC = HERE / "data_ghg"

#: Each case renames ONE instrument. `real` is what the archive says today,
#: `standin` what EddyPro must be shown, `ef` what EddyFlow should read, `firm`
#: the manufacturer that goes with the stand-in - EddyPro validates firm and
#: model against separate lists, and its IRGA firms are only `licor` and
#: `other_irga`, its sonic firms only gill/metek/young/csi/other_sonic.
#:
#: The geometry is REQUIRED for a generic model - EddyPro refuses one whose
#: hpath * vpath * tau is zero - and must be the instrument's real geometry.
#: Stating the truth is the point rather than a nicety: a stand-in carrying
#: plausible but wrong numbers - 12.5 cm of horizontal path where the LI-7500A
#: has 0.95 - moves LE by 0.24 % in EddyPro and 0.21 % in EddyFlow. With the
#: real numbers both engines are bit-identical to reading the instrument
#: directly. A named model needs none, since both engines look its geometry up.
CASES = {
    "data_ghg_ext": dict(
        real="li7500a_1",
        standin="generic_open_path_1",
        firm="other_irga",
        ef="li7500a_1",
        #: RetrieveSensorParams' own numbers for the LI-7500A - 0.127 m
        #: vertical, 0.0095 m horizontal, 0.1 s - in the centimetres a
        #: .metadata states.
        geom=[("hpath_length", "0.95"), ("vpath_length", "12.7"),
              ("tau", "0.1")],
    ),
    "data_ghg_campbell": dict(
        real="usoni3_classa_mp_1",
        #: EddyPro's own spelling, and the only one it takes: upstream's
        #: metadata_file_validation.f90 has `csat3` and `csat3b` and no
        #: prefixed form of either.
        standin="csat3b_1",
        firm="csi",
        ef="csi_csat3b_1",
        #: None: csat3b is a model EddyPro knows, so it is a rename rather than
        #: a stand-in and both engines supply the geometry themselves.
        geom=[],
    ),
}


def sevenzip():
    for exe in ("7z", "7za"):
        if shutil.which(exe):
            return exe
    sys.exit("7-Zip not on PATH; the GHG fixtures need it too")


def rewrite(text, case):
    """Classic .metadata -> extended .metadata. Byte-preserving elsewhere."""
    lines = text.split("\n")

    #: Which instr_<k> block to rewrite, found in a FIRST pass. The block
    #: cannot be learnt from the model line as it goes past, because
    #: instr_<k>_manufacturer comes BEFORE instr_<k>_model in every file
    #: LI-COR writes - so a single pass reaches the manufacturer without yet
    #: knowing whether this is the block being renamed, and silently leaves it
    #: naming the old maker.
    block = None
    for line in lines:
        key = line.split("=", 1)[0].strip()
        if key.endswith("_model") and line.strip().endswith("=" + case["real"]):
            block = key[: -len("_model")]
            break
    if block is None:
        sys.exit("no instrument block declares model=%s" % case["real"])

    out, did = [], False
    for line in lines:
        key = line.split("=", 1)[0].strip()

        #: The version key goes at the end of [FileDescription], the section
        #: that describes the file rather than the site.
        if line.strip() == "[FileDescription]":
            out += [line, "ghg_format_version=2.0"]
            continue

        if key == f"{block}_model":
            out.append(f"{block}_model=" + case["standin"])
            continue

        #: The manufacturer has to move with the model. EddyPro checks the two
        #: against separate lists, and a Campbell sonic still calling itself
        #: `metek`, or a generic analyser still calling itself `licor`, fails
        #: the firm test rather than the model one - and the message names
        #: neither.
        if key == f"{block}_manufacturer":
            out.append(f"{block}_manufacturer=" + case["firm"])
            continue

        #: Geometry and the true identity are appended to the block being stood
        #: in for. sw_version is the key that follows model in every file
        #: LI-COR writes.
        if key == f"{block}_sw_version" and not did:
            out.append(line)
            out += [f"{block}_{k}={v}" for k, v in case["geom"]]
            out.append(f"{block}_ef_model=" + case["ef"])
            did = True
            continue

        if key.endswith("_instrument") and line.strip().endswith("=" + case["real"]):
            out.append(f"{key}=" + case["standin"])
            continue

        out.append(line)

    if not did:
        sys.exit("no instrument block declares model=%s" % case["real"])
    return "\n".join(out)


def build(name, case, z, archives):
    dst = HERE / name
    dst.mkdir(exist_ok=True)
    for ghg in archives:
        target = dst / ghg.name
        shutil.copy2(ghg, target)
        with tempfile.TemporaryDirectory() as tmp:
            #: The archive holds TWO .metadata files - the site's, and the
            #: biomet one beside it. Extracted by exact name rather than by
            #: glob, because <stem>-biomet.metadata sorts first and describes
            #: no instruments at all: taking it produces a rewrite that changes
            #: nothing and a fixture that proves nothing.
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
                rewrite(md.read_bytes().decode("cp1252"), case).encode("cp1252"))
            #: str(target) is ABSOLUTE. With cwd set to the temp directory a
            #: relative archive path makes 7-Zip silently CREATE a new archive
            #: there and exit 0, leaving the real one untouched.
            subprocess.run([z, "u", str(target), want, "-y"],
                           check=True, cwd=tmp, stdout=subprocess.DEVNULL)
    print(f"  {name}: {case['real']} -> {case['standin']} "
          f"(ef_model={case['ef']}), {len(archives)} archive(s)")


def main():
    z = sevenzip()
    archives = sorted(SRC.glob("*.ghg"))
    if not archives:
        sys.exit(f"no archives in {SRC}")
    for name, case in CASES.items():
        build(name, case, z, archives)


if __name__ == "__main__":
    main()
