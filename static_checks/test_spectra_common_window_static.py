"""Every (co)spectrum is computed on one window, the same rows for each column.

FixDatasetForSpectra interpolated each column over all its rows, and
ReplaceGapWithLinearInterpolation drops a column's LEADING error rows by moving
the rest of it up. The full-rate pass pairs w and a gas row by row, so any
column that opened on a gap was out of step with w by that many rows - in its
cospectrum and in the covariance that cospectrum was divided by. A negative
lag opens a column with exactly such a gap, so for it the move undid the lag
compensation; data missing or filtered at the start of a period did the same.
On CH-LAE the MIRO's normalising covariance matched its flux in a fifth to a
half of the half-hours. A gas identical to w and opening with a 40-row gap
came out with a correlation of 0.06 against w; on one window, 0.9994.

Pinned: each column is interpolated over its own span and nothing moves; the
caller takes one window from FixDatasetForSpectra rather than trimming by the
spread of RowLags after the columns had moved; and a window under two thirds
of the period skips the spectra rather than transforming a fragment.
"""

import unittest
from pathlib import Path

SRC = Path(__file__).resolve().parent.parent / "src" / "src_rp"
FIX = (SRC / "fix_dataset_for_spectra.f90").read_text(encoding="utf-8", errors="replace")
MAIN = (SRC / "eddyflow-rp_main.f90").read_text(encoding="utf-8", errors="replace")


def spectral_block():
    a = MAIN.index("if (RPsetup%do_spectral_analysis) then")
    return MAIN[a:MAIN.index("call SpectralAnalysis(", a) + 200]


class EachColumnStaysWhereItIs(unittest.TestCase):

    def test_interpolated_over_its_own_span(self):
        self.assertIn("call ReplaceGapWithLinearInterpolation(Set(first:last, j),", FIX)
        self.assertNotIn("ReplaceGapWithLinearInterpolation(Set(1:nrow, j)", FIX)

    def test_the_window_is_where_every_column_has_data(self):
        self.assertIn("lo = max(lo, first)", FIX)
        self.assertIn("hi = min(hi, last)", FIX)
        self.assertIn("nrow2 = max(0, hi - lo + 1)", FIX)

    def test_nothing_is_shifted(self):
        self.assertIn("SpecLead = 0", FIX)
        self.assertNotIn("SpecLead(j) = ", FIX)

    def test_the_window_is_handed_back(self):
        self.assertIn("subroutine FixDatasetForSpectra(Set, nrow, ncol, lo, nrow2)", FIX)
        self.assertIn("integer, intent(out) :: lo", FIX)


class TheCallerUsesThatWindow(unittest.TestCase):

    def test_the_set_is_cut_from_the_window(self):
        block = spectral_block()
        self.assertIn("size(E2Primes, 1), size(E2Primes, 2), SpecLo, N2)", block)
        self.assertIn("SpecSet(1: SpecRow, :) = E2Primes(SpecLo: SpecLo + SpecRow - 1, :)", block)
        self.assertIn("SpecRowOffset = SpecLo - 1", block)

    def test_no_second_trim_by_the_lags(self):
        block = spectral_block()
        self.assertNotIn("maxval(RowLags)", block)
        self.assertNotIn("minval(RowLags)", block)
        self.assertNotIn("Nmin = - Nmin", block)

    def test_a_short_window_skips_the_spectra(self):
        block = spectral_block()
        guard = block.index("if (N2 < (2 * size(E2Primes, 1)) / 3) then")
        self.assertLess(guard, block.index("call SpectralAnalysis("))
        self.assertIn("Essentials%degH(:) = error", block[guard:])


if __name__ == "__main__":
    unittest.main()
