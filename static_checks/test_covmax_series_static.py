"""Covariance maximisation computes the one covariance it reads, and moves no bit.

``CovMax`` is what time-lag optimisation runs for every gas of every period,
and what PWB runs for its terminal fallback. For **every lag** it used to
allocate two arrays, copy the shifted pair into one and that into the other,
and hand the pair to ``CovarianceMatrixNoError`` - which computes the full 2x2
matrix, four passes over the data, of which only ``Cov(1,2)`` was read.

``CrossCovarianceSeries`` (m_covmax_core) computes ``Cov(1,2)`` alone, reading
the shifted elements in place. For that entry it performs exactly what the
matrix routine performs: the same pairwise test against the error code, the
same products, sums and count accumulated in the same order, and the same
``C/N - (si/N)(sj/N)`` at the end. So every value is bitwise the same.

Measured at the Yatir shape (36 000 rows at 20 Hz, a 0-25 s window, a gas on
one row in twenty): 0.27 s -> 0.036 s per gas per period, 7.3x. With ten gases
that is about 2.7 s of every period of a time-lag optimisation pre-pass.

The proof is ``covmax_bitexact`` (``mingw32-make covmaxcheck``): the old
per-lag path frozen verbatim, compared on raw bits over 10 812 cases at
error-code densities from none to all. It ran clean.

The stochastic-detrending path is left exactly as it was: there each shifted
window is detrended on its own, so the copies are real work.

Part of the EddyFlow engine's static checks.
"""

from pathlib import Path
import subprocess
import unittest


ROOT = Path(__file__).resolve().parents[1]

CORE = ROOT / "src" / "src_common" / "m_covmax_core.f90"
HANDLE = ROOT / "src" / "src_rp" / "timelag_handle.f90"
STATS = ROOT / "src" / "src_common" / "stats_operator_no_error.f90"
HARNESS = ROOT / "src" / "src_tools" / "covmax_bitexact_main.f90"
MAKEFILE = ROOT / "prj" / "Makefile"
DRIVER = ROOT / "obj" / "win" / "covmax_bitexact.exe"


def read(path):
    return path.read_text(encoding="utf-8", errors="replace")


def code(path):
    return (chr(10)).join(ln for ln in read(path).splitlines()
                          if not ln.lstrip().startswith("!"))


SERIES = code(CORE)
COVMAX = code(HANDLE)
COVMAX = COVMAX[COVMAX.index("subroutine CovMax("):COVMAX.index("end subroutine CovMax")]


class OnlyTheEntryThatIsReadIsComputed(unittest.TestCase):

    def test_covmax_calls_the_series_routine_when_not_detrending(self):
        self.assertIn("if (.not. RPSetup%covmax_stocdet) then", COVMAX)
        self.assertIn("call CrossCovarianceSeries(Col1, Col2, nrow, lagmin, lagmax, error, CovSeries)",
                      COVMAX)

    def test_the_detrending_path_still_builds_and_detrends_each_window(self):
        tail = COVMAX[COVMAX.index("call CrossCovarianceSeries("):]
        self.assertIn("call VariableStochasticDetrending(", tail)
        self.assertIn("call CovarianceMatrixNoError(", tail)

    def test_the_routine_is_pure_and_copies_nothing(self):
        self.assertIn("pure subroutine CrossCovarianceSeries", SERIES)
        self.assertNotIn("allocate", SERIES)


class ItDoesWhatTheMatrixRoutineDidForThatEntry(unittest.TestCase):
    """The bitwise claim rests on these lines matching the matrix routine's."""

    def test_the_same_pairwise_test(self):
        self.assertIn("if (Set(k, i) /= err_float .and. Set(k, j) /= err_float) then", code(STATS))
        self.assertIn("if (a /= err .and. b /= err) then", SERIES)

    def test_the_same_accumulations(self):
        for line in ("nact = nact + 1", "c = c + a * b", "si = si + a", "sj = sj + b"):
            self.assertIn(line, SERIES)

    def test_the_same_normalisation_in_the_same_order(self):
        steps = ["si = si / dble(nact)", "sj = sj / dble(nact)",
                 "c = c / dble(nact)", "c = c - si * sj"]
        pos = [SERIES.index(s) for s in steps]
        self.assertEqual(pos, sorted(pos))

    def test_the_same_shift_convention(self):
        self.assertIn("a = col1(k - lag)", SERIES)
        self.assertIn("b = col2(k + lag)", SERIES)

    def test_no_valid_pair_gives_the_error_code(self):
        self.assertIn("c = err", SERIES)


class TheProofIsKeptAndRuns(unittest.TestCase):

    def test_the_harness_freezes_the_old_path(self):
        h = read(HARNESS)
        self.assertIn("subroutine FrozenSeries(", h)
        self.assertIn("subroutine FrozenCovMatrix(", h)
        self.assertIn("out(i) = CovMat(1, 2)", h)

    def test_it_covers_sparse_columns(self):
        self.assertIn("0.95d0", read(HARNESS))

    def test_it_compares_raw_bits(self):
        self.assertIn("transfer(a(k), 0_int64) /= transfer(b(k), 0_int64)", read(HARNESS))

    def test_it_is_linked_against_the_shipped_object(self):
        self.assertIn("covmaxcheck : m_numeric_kinds.o m_covmax_core.o", read(MAKEFILE))

    @unittest.skipUnless(DRIVER.exists(), "covmax_bitexact not built (mingw32-make covmaxcheck)")
    def test_the_shipped_routine_is_bit_identical(self):
        r = subprocess.run([str(DRIVER)], capture_output=True, text=True, timeout=600)
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertIn("BIT-IDENTICAL", r.stdout)
        self.assertIn("mismatched elements: 0", r.stdout)


if __name__ == "__main__":
    unittest.main()
