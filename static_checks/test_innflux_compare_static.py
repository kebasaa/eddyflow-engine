"""The innFLUX-derived peak search and its comparison harness.

innFLUX: Striednig, M., Graus, M., Maerk, T. D., and Karl, T. G., Atmos. Meas.
Tech., 13, 1447-1465, https://doi.org/10.5194/amt-13-1447-2020, 2020; MATLAB
code v1.2.0, University of Innsbruck. Smoothing the cross-covariance before
taking its extremum follows Taipale et al. (2010).

``m_lagpeak_core`` reimplements innFLUX's find_covariance_peak.m: a Hamming
window smoothed through MATLAB filtfilt, and the sign of the peak decided from
a 100-point moving average of the function's derivative. filtfilt of an FIR
filter is, away from the ends, a convolution with the zero-phase kernel
sum_{j-k=m} b(j) b(k); the engine computes the covariance PeakPadding lags
past the search window so only interior values are searched.

tests/innflux/compare.py compares the module with a line-by-line Python port
of the MATLAB code (needs scipy - the conda env ``dp``). This file pins what
does not need scipy, and runs the driver on hand-computable cases when it has
been built (``cmake --build --preset release --target innfluxcmp``).

Part of the EddyFlow engine's static checks.
"""

from pathlib import Path
import re
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
DRIVER = ROOT / "obj" / "win" / "innflux_compare.exe"


def read(*parts):
    return ROOT.joinpath(*parts).read_text(encoding="utf-8", errors="replace")


def code(*parts):
    return "\n".join(ln for ln in read(*parts).splitlines()
                     if not ln.lstrip().startswith("!"))


CORE = code("src", "src_common", "m_lagpeak_core.f90")
CMAKE = read("CMakeLists.txt")


class TestModule(unittest.TestCase):
    def test_public_interface(self):
        for name in ("HammingKernel", "ZeroPhaseKernel", "ZeroPhaseSmooth",
                     "PeakPadding", "FindCovPeak"):
            self.assertRegex(CORE, rf"public :: {name}\b")

    def test_it_is_pure(self):
        """A function of its arguments only, so the driver can check it."""
        for name in ("HammingKernel", "ZeroPhaseKernel", "ZeroPhaseSmooth",
                     "FindCovPeak"):
            self.assertRegex(CORE, rf"pure subroutine {name}\b")
        self.assertNotIn("use m_rp_global_var", CORE)

    def test_matlab_hamming_is_the_symmetric_window(self):
        """hamming(L): 0.54 - 0.46 cos(2 pi k / (L - 1)), k = 0..L-1."""
        self.assertIn("0.54d0 - 0.46d0 * cos(2d0 * pi * dble(k) / dble(l - 1))", CORE)
        self.assertIn("b = b / sum(b)", CORE)

    def test_derivative_filter_is_innflux_s(self):
        self.assertRegex(CORE, r"DerivFilterLen = 100\b")

    def test_padding_covers_both_filters(self):
        self.assertIn("pad = smooth_len - 1", CORE)
        self.assertIn("pad = max(pad, DerivFilterLen)", CORE)

    def test_the_source_is_credited(self):
        header = read("src", "src_common", "m_lagpeak_core.f90")[:3000]
        self.assertIn("Striednig et al. 2020", header)
        self.assertIn("Taipale et al. (2010)", header)
        self.assertIn("find_covariance_peak.m", header)


class TestHarness(unittest.TestCase):
    def test_fixture_target(self):
        self.assertIn("eddyflow_fixture(innfluxcmp innflux_compare "
                      "innflux_compare_main.f90", CMAKE)
        self.assertRegex(CMAKE, r"innflux_core_objects m_numeric_kinds "
                                r"m_covmax_core m_lagpeak_core")

    def test_harness_files(self):
        for p in ("tests/innflux/compare.py", "tests/innflux/innflux_port.py",
                  "tests/innflux/README.md",
                  "src/src_tools/innflux_compare_main.f90"):
            self.assertTrue((ROOT / p).is_file(), p)

    def test_port_pads_like_matlab(self):
        """scipy's default pad is one filter length longer than MATLAB's."""
        port = read("tests", "innflux", "innflux_port.py")
        self.assertIn("padlen=3 * (len(b) - 1)", port)
        self.assertIn("hamming(filter_length, sym=True)", port)


@unittest.skipUnless(DRIVER.exists(), "innflux_compare not built")
class TestDriver(unittest.TestCase):
    def run_cases(self, text):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "c.txt"
            p.write_text(text)
            r = subprocess.run([str(DRIVER), str(p)], capture_output=True,
                               text=True, timeout=60)
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        return [l.split() for l in r.stdout.splitlines() if l.strip()]

    def triangle(self, centre, sign, lo=-300, hi=300):
        return [sign * max(0.0, 50.0 - abs(k - centre)) for k in range(lo, hi + 1)]

    def peak_case(self, values, smooth, mode, lo=-300, hi=300, slo=-100, shi=100):
        body = "\n".join(f"{v:.17e}" for v in values)
        return f"peak\n{lo} {hi} {slo} {shi} {smooth} {mode} 0\n{body}\n"

    def test_positive_and_negative_peaks(self):
        cases = (self.peak_case(self.triangle(37, 1.0), 15, 2)
                 + self.peak_case(self.triangle(-23, -1.0), 15, 2)
                 + self.peak_case(self.triangle(-23, -1.0), 15, 1)
                 + self.peak_case(self.triangle(-23, -1.0), 0, 0))
        out = self.run_cases(cases)
        # Derivative rule: the maximum of a positive peak, the minimum of a
        # negative one. Positive-only finds no maximum inside a negative
        # triangle's flat zero, so it returns the latest of the tied zeros.
        self.assertEqual(out[0][1], "37")
        self.assertEqual(out[1][1], "-23")
        self.assertEqual(out[2][1], "100")
        self.assertEqual(out[3][1], "-23")

    def test_smoothing_is_zero_phase(self):
        """A symmetric peak stays where it was however hard it is smoothed."""
        out = self.run_cases(self.peak_case(self.triangle(12, 1.0), 41, 1))
        self.assertEqual(out[0][1], "12")


if __name__ == "__main__":
    unittest.main()
