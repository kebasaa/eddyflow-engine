"""The smoothed covariance peak search after innFLUX, opt-in.

innFLUX (Striednig et al. 2020, Atmos. Meas. Tech. 13, 1447-1465;
find_covariance_peak.m, v1.2.0) smooths the cross-covariance with a Hamming
window before taking the peak (after Taipale et al. 2010) and can decide the
peak's sign from the function's derivative. Two keys bring that to CovMax:

  covmax_smooth_len  RP.SNTags(425)  Hamming window in rows, 0 = off
  covmax_peak_sign   RP.SNTags(426)  0 |max| (as before), 1 max, 2 innFLUX's rule

Both default to off, and off must leave CovMax's own search untouched: the
innFLUX branch takes its own, wider series and returns before anything the
default path reads. The arithmetic is checked against a port of the MATLAB
code in tests/innflux/compare.py and test_innflux_compare_static.py.

Part of the EddyFlow engine's static checks.
"""

from pathlib import Path
import re
import unittest

ROOT = Path(__file__).resolve().parents[1]


def read(*parts):
    return ROOT.joinpath(*parts).read_text(encoding="utf-8", errors="replace")


def code(*parts):
    return "\n".join(ln for ln in read(*parts).splitlines()
                     if not ln.lstrip().startswith("!"))


TAGS = read("src", "src_rp", "m_rp_global_var.f90")
GEN = read("prj", "gen_project_tags.py")
READER = code("src", "src_rp", "read_ini_rp.f90")
TYPEDEF = code("src", "src_common", "m_typedef.f90")
TL = code("src", "src_rp", "timelag_handle.f90")
IMPORT = read("src", "src_common", "m_eddypro_import.f90")


def covmax_body():
    m = re.search(r"subroutine CovMax\(.*?end subroutine CovMax", TL, re.S)
    assert m, "CovMax not found"
    return m.group(0)


class TheKeys(unittest.TestCase):
    def test_slots(self):
        self.assertRegex(TAGS, r"SNTags\(425\)%Label\s*/\s*'covmax_smooth_len'")
        self.assertRegex(TAGS, r"SNTags\(426\)%Label\s*/\s*'covmax_peak_sign'")
        self.assertIn('425: "covmax_smooth_len"', GEN)
        self.assertIn('426: "covmax_peak_sign"', GEN)

    def test_setup_fields(self):
        self.assertRegex(TYPEDEF, r"integer :: covmax_smooth_len")
        self.assertRegex(TYPEDEF, r"integer :: covmax_peak_sign")

    def test_read_with_literal_defaults_and_found_guards(self):
        self.assertIn("RPSetup%covmax_smooth_len = 0", READER)
        self.assertIn("RPSetup%covmax_peak_sign = 0", READER)
        self.assertIn("if (SNTagFound(425)) then", READER)
        self.assertIn("if (SNTagFound(426)) then", READER)

    def test_out_of_range_values_fall_back_to_off(self):
        self.assertIn("if (RPSetup%covmax_smooth_len < 2) RPSetup%covmax_smooth_len = 0",
                      READER)
        self.assertRegex(READER, r"covmax_peak_sign < 0 \.or\. RPSetup%covmax_peak_sign > 2")

    def test_not_combined_with_stochastic_detrending(self):
        self.assertRegex(READER, r"if \(RPSetup%covmax_stocdet \.and\. "
                                 r"\(RPSetup%covmax_smooth_len > 0")

    def test_the_eddypro_import_does_not_state_them(self):
        """Absent means off, which is the reader's default: nothing to pin."""
        self.assertNotIn("covmax_smooth_len", IMPORT)
        self.assertNotIn("covmax_peak_sign", IMPORT)


class CovMaxBranch(unittest.TestCase):
    def setUp(self):
        self.body = covmax_body()

    def test_the_branch_is_taken_only_when_asked(self):
        self.assertIn("if (RPSetup%covmax_smooth_len > 1 .or. "
                      "RPSetup%covmax_peak_sign /= 0) then", self.body)

    def test_the_branch_returns_before_the_default_path(self):
        branch = self.body.index("call FindCovPeak(")
        ret = self.body.index("return", branch)
        default = self.body.index("allocate(CovSeries(nlag))")
        self.assertLess(ret, default)

    def test_it_reads_a_padded_series(self):
        self.assertIn("PeakPadding(RPSetup%covmax_smooth_len", self.body)
        self.assertIn("call CrossCovarianceSeries(Col1, Col2, nrow, elo, ehi, error, "
                      "ExtSeries)", self.body)

    def test_the_default_path_is_as_it_was(self):
        self.assertIn("call CrossCovarianceSeries(Col1, Col2, nrow, lagmin, lagmax, "
                      "error, CovSeries)", self.body)
        self.assertIn("score = dabs(CovSeries(k))", self.body)


if __name__ == "__main__":
    unittest.main()
