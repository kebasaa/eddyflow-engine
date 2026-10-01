"""A gas sampled slower than the file is detected at its own rate, refined at the file's.

The Yatir inputs carry a 1 Hz laser on the 20 Hz sonic grid, nineteen rows in
twenty empty. PwbDetectGas counted that gas's completeness against the row
grid - 5 % - so every laser gas, COS included, failed its validity test on
every period and ended in the terminal fallback.

PwbDetectSlowGas now takes any gas whose instrument is slower than the file:

* completeness is counted against what that instrument owes;
* stage 1 runs the unchanged pre-whitening chain on the gas's REAL samples
  (SlowColumnSampleRows: its true rows - the phase drifts at Yatir - with a
  slot for each missed sample), with w and Ts taken at the same rows, point or
  interval-averaged as the instrument samples, and windows, block length and
  smoothing converted at the measured spacing;
* stage 2 evaluates every row lag within one gas interval of the coarse peak,
  pre-whitened with the winner's filter and bootstrapped on the winner's own
  stream, so the lag and HDI come out at the file's resolution and the lag in
  file rows - what the streaming pass shifts by;
* where stage 1's replicates spread over more than an interval and a half, the
  reported HDI covers both stages, so stage 2's narrow range cannot make an
  uncertain period look certain.

Every full-rate gas takes the old path, bit for bit: RunPwbCombination gained
the rate and smoothing width as arguments and the full-rate caller passes the
file's rate and the configured width. Gate: check_slow_lag.sh - COS made from
w delayed by exactly 3.35 s on a drifting 0.987 Hz grid is recovered as 3.3 or
3.4 s in every period (the previous build fell back and said 8.1 or 10.0).

Part of the EddyFlow engine's static checks.
"""

from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[1]


def code(rel):
    text = (ROOT / rel).read_text(encoding="utf-8", errors="replace")
    return "\n".join(ln for ln in text.splitlines() if not ln.lstrip().startswith("!"))


PWB = code("src/src_rp/pwb_timelag_handle.f90")
SAMPLING = code("src/src_rp/column_sampling.f90")


def body(text, name, kind="subroutine"):
    return text[text.index("%s %s(" % (kind, name)):text.index("end %s %s" % (kind, name))]


DETECT = body(PWB, "PwbDetectGas")
SLOW = body(PWB, "PwbDetectSlowGas")
RUN = body(PWB, "RunPwbCombination")


class OnlySlowGasesTakeTheNewPath(unittest.TestCase):

    def test_the_branch_is_on_the_column_rate(self):
        self.assertIn("if (ColumnAcFreq(gas) < Metadata%ac_freq) then", DETECT)
        i = DETECT.index("if (ColumnAcFreq(gas) < Metadata%ac_freq) then")
        self.assertIn("call PwbDetectSlowGas(", DETECT[i:i + 200])
        self.assertIn("return", DETECT[i:i + 250])

    def test_full_rate_combinations_pass_todays_rate_and_width(self):
        self.assertEqual(DETECT.count("Metadata%ac_freq, PWBSetup%smoothing_width, candidate("), 4)

    def test_the_combination_runner_uses_only_its_arguments(self):
        self.assertNotIn("Metadata%ac_freq", RUN)
        self.assertNotIn("PWBSetup%smoothing_width", RUN)
        self.assertIn("requested_block_len = nint(PWBSetup%block_length_s * rate)", RUN)


class StageOneRunsOnTheRealSamples(unittest.TestCase):

    def test_samples_are_found_not_assumed(self):
        self.assertIn("call SlowColumnSampleRows(Set(:, gas), nrow, stride, error, rows, isreal, ns, nreal)", SLOW)
        helper = body(SAMPLING, "SlowColumnSampleRows")
        self.assertIn("if (2 * gap > 3 * stride) then", helper)
        self.assertIn("nmiss = nint(dble(gap) / dble(stride)) - 1", helper)

    def test_completeness_counts_against_the_instruments_own_samples(self):
        self.assertIn("expected = max(1, nint(dble(nrow) * ColumnAcFreq(gas) / Metadata%ac_freq))", SLOW)
        self.assertIn("dble(nreal) < min_valid * dble(expected)", SLOW)

    def test_conversions_use_the_measured_spacing(self):
        self.assertIn("dbar = dble(rows(ns) - rows(1)) / dble(ns - 1)", SLOW)
        self.assertIn("rate1 = Metadata%ac_freq / dbar", SLOW)
        self.assertIn("min_s = floor(dble(min_rl) / dbar)", SLOW)
        self.assertIn("max_s = ceiling(dble(max_rl) / dbar)", SLOW)

    def test_the_unchanged_chain_runs_at_the_gas_rate(self):
        self.assertIn("call PwbPreWhiten(ss, ww, tt, ns, min_s, max_s,", SLOW)
        self.assertEqual(SLOW.count("rate1, swidth1, candidate("), 4)

    def test_the_driver_is_taken_the_way_the_gas_was(self):
        driver = SLOW[SLOW.index("function DriverAt("):SLOW.index("end function DriverAt")]
        self.assertIn("if (E2Col(gas)%instr%integrates) then", driver)
        self.assertIn("lo = max(1, r - stride + 1)", driver)


class StageTwoRefinesAtTheRows(unittest.TestCase):

    def test_within_one_interval_of_the_coarse_peak(self):
        self.assertIn("half = ceiling(dbar)", SLOW)
        self.assertIn("lo2 = max(min_rl, LocResult%row_lag - half)", SLOW)
        self.assertIn("hi2 = min(max_rl, LocResult%row_lag + half)", SLOW)

    def test_gas_sample_k_pairs_with_the_driver_l_rows_earlier(self):
        """The CCF pairs x(i) with y(i + lag): the gas lags the driver."""
        self.assertIn("xk(k) = DriverAt(wfull, rows(k1 + k - 1) - l)", SLOW)

    def test_prewhitened_with_the_winners_filter_on_the_winners_stream(self):
        self.assertIn("call ApplyArFilter(Differenced(xk, n2, pw%differenced), ne2, phi, p, xf)", SLOW)
        self.assertIn("call ApplyArFilter(Differenced(yk, n2, pw%differenced), ne2, phi, p, yf)", SLOW)
        self.assertIn("state = PwbStreamSeed(gas, combo(best))", SLOW)

    def test_the_lag_leaves_in_file_rows(self):
        self.assertIn("LocResult%row_lag = lag", SLOW)
        self.assertIn("LocResult%selected_lag = dble(lag) / Metadata%ac_freq", SLOW)

    def test_the_refinement_edge_counts_as_pinned(self):
        self.assertIn("LocResult%edge_pinned = lag == lo2 .or. lag == hi2", SLOW)

    def test_stage_one_uncertainty_reaches_the_hdi(self):
        self.assertIn("if (h1hi - h1lo > 1.5d0 * dbar / Metadata%ac_freq) then", SLOW)
        self.assertIn("LocResult%hdi_low = min(LocResult%hdi_low, h1lo)", SLOW)
        self.assertIn("LocResult%hdi_high = max(LocResult%hdi_high, h1hi)", SLOW)


if __name__ == "__main__":
    unittest.main()
