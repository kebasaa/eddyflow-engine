"""RP and FCC write one metadata file, in one layout and one set of units.

Which program writes `*_metadata_*.csv` depends only on whether FCC follows
RP (an in-situ spectral correction or ensemble spectra), so the user cannot
tell which one they got - and the two disagreed:

  - RP's header had a DOY column after time, FCC's did not, so every column
    after time sat under the wrong name in one of the two.
  - RP wrote lengths in cm, tube diameter in mm and flow rate in l min-1 -
    the units of the .metadata file - and FCC wrote the SI values the ex
    reader had converted them to: 0.25 from FCC where RP wrote 25.0.
  - RP scaled the error code along with the value: a missing tube length
    came out as -999900 and a missing flow rate as -599940000.
  - FCC wrote the row list-directed, which may lead with a blank.

EddyPro 7 has the same split; it is fixed here, not inherited. Latitude and
longitude still differ in precision: FCC reads them back from the ex file,
which carries six significant digits (G20.6).
"""

from pathlib import Path
import re
import unittest


ROOT = Path(__file__).resolve().parents[1]

RP_WRITER = "src/src_rp/write_out_metadata.f90"
FCC_WRITER = "src/src_fcc/write_out_metadata_fcc.f90"
RP_INIT = "src/src_rp/init_outfiles_rp.f90"
FCC_INIT = "src/src_fcc/init_out_files.f90"

#> The same quantity under each program's name for it.
SAME = {"d": "disp_height", "z0": "rough_length"}


def code(path):
    """Source with comment-only lines dropped."""
    text = (ROOT / path).read_text(encoding="utf-8")
    return "\n".join(ln for ln in text.splitlines()
                     if not ln.lstrip().startswith("!"))


def metadata_header(path):
    """The fixed header literal and the per-gas suffixes after it."""
    body = code(path)
    block = body[body.index("out_md) then"):]
    block = block[:block.index("end if\n", block.index("end do"))]
    lit = re.search(r"AddDatum\(\w+,\s*'(filename,.*?)', separator\)",
                    block, re.S)
    assert lit, path
    fixed = re.sub(r"&\s*\n\s*&", "", lit.group(1))
    per_gas = re.findall(r"'(irga_\w+|measure_type),?'", block)
    return fixed, per_gas


def written_fields(path):
    """[(field, gain)] in the order the row writes them, gain None if
    unscaled. Field is the last component (lat, hpath_length, ...)."""
    out = []
    for m in re.finditer(r"write\(field_val, \*\) (.+)", code(path)):
        expr = m.group(1).strip()
        unit = re.fullmatch(r"OutUnits\((.+), (\w+)\)", expr)
        gain = None
        if unit:
            expr, gain = unit.group(1), unit.group(2)
        expr = re.sub(r"\(1:len_trim\(.*\)\)$", "", expr)
        field = re.split(r"%", expr)[-1]
        field = re.sub(r"\(.*\)$", "", field)
        out.append((SAME.get(field, field), gain))
    return out


class OneHeader(unittest.TestCase):

    def test_fixed_part_is_the_same(self):
        rp, _ = metadata_header(RP_INIT)
        fcc, _ = metadata_header(FCC_INIT)
        self.assertEqual(rp, fcc)
        self.assertTrue(rp.startswith("filename,date,time,DOY,"), rp)

    def test_gas_blocks_are_the_same(self):
        _, rp = metadata_header(RP_INIT)
        _, fcc = metadata_header(FCC_INIT)
        self.assertEqual(len(rp), 14)
        self.assertEqual(rp, fcc)


class OneRow(unittest.TestCase):

    def test_fcc_writes_doy_after_time(self):
        body = code(FCC_WRITER)
        t = body.index("lEx%end_time(1:5)")
        doy = body.index("float_doy", t)
        self.assertLess(doy, body.index("lEx%lat", t))
        self.assertIn("index(field_val, '.') + 4", body[doy:])

    def test_same_fields_same_units(self):
        rp = written_fields(RP_WRITER)
        fcc = written_fields(FCC_WRITER)
        #> FCC writes its own timestamp fields; RP gets them as one string.
        start = [f for f, _ in fcc].index("lat")
        self.assertEqual(rp, fcc[start:])

    def test_lengths_are_scaled(self):
        """The SI fields are the ones scaled, so a field dropped from the
        scaling cannot pass by being dropped from both."""
        scaled = {f: g for f, g in written_fields(RP_WRITER) if g}
        for f in ("hpath_length", "vpath_length", "nsep", "esep", "vsep",
                  "tube_l"):
            self.assertEqual(scaled.get(f), "1d2", f)
        self.assertEqual(scaled.get("tube_d"), "1d3")
        self.assertEqual(scaled.get("tube_f"), "6d4")


class ErrorCodeIsNotScaled(unittest.TestCase):

    def test_no_bare_scaling(self):
        for path in (RP_WRITER, FCC_WRITER):
            self.assertNotRegex(code(path), r"\* *(1d2|1d3|6d4)\s*\n", path)

    def test_outunits_passes_error_through(self):
        for path in (RP_WRITER, FCC_WRITER):
            self.assertRegex(
                code(path),
                r"OutUnits = merge\(si \* gain, si, si /= error\)", path)

    def test_rows_are_written_verbatim(self):
        for path in (RP_WRITER, FCC_WRITER):
            self.assertIn("write(umd, '(a)')", code(path), path)


if __name__ == "__main__":
    unittest.main()
