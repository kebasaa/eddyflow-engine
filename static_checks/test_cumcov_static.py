"""Time lag from cumulated cross-covariances (tlag_meth 6), and clock drift.

After innFLUX (Striednig et al. 2020, Atmos. Meas. Tech. 13, 1447-1465,
Sect. 2.3.2; innFLUX_step2.m), who credit Park et al. (2013) with summing
covariance functions; the clock-drift correction after innFLUX v1.2.0's
driftfix tool (2023, unpublished; M. Graus, M. Striednig). Opt-in: a new
time-lag method, and drift only with tlag_drift_meth=1.

What this pins, because a mistake in any of it gives a plausible lag rather
than an error:

- The pre-pass writes one record per period and gas to a scratch file per
  process and sums nothing; the parent orders them by period before summing,
  so a split pre-pass gives the same sums as a serial one.
- newunit= returns NEGATIVE unit numbers: "is it open" is a flag, never a
  sign test. The first build tested `unit < 0` and gathered nothing.
- Production applies the lag by the period's own stamp (Stats%date/time), the
  end of the period, and segments place a period by its start, as innFLUX
  compares interval start times.
- Drift is taken relative to its median, so a drift-shifted function keeps its
  peak inside the gas's search window.

tests/innflux/compare.py stage 2 checks the drift estimate against a port of
the MATLAB tool; base_cumcov and base_cumcov_drift run the method end to end.

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
CUMCOV = code("src", "src_rp", "m_cumcov.f90")
CUMCOV_RAW = read("src", "src_rp", "m_cumcov.f90")
DRIFT = read("src", "src_common", "m_driftfix_core.f90")
MAIN = code("src", "src_rp", "eddyflow-rp_main.f90")
TL = code("src", "src_rp", "timelag_handle.f90")
FLX = code("src", "src_rp", "write_out_fluxnet.f90")


class TheKeys(unittest.TestCase):
    def test_slots(self):
        for idx, key in ((429, "tlag_cum_filter_len"), (430, "tlag_cum_min_ustar"),
                         (431, "tlag_cum_min_compl_sonic"), (432, "tlag_cum_min_compl_gas"),
                         (433, "tlag_cum_max_sst"), (434, "tlag_drift_meth"),
                         (435, "tlag_drift_span_d"), (436, "tlag_drift_jump_thr"),
                         (437, "tlag_drift_max_s"), (438, "out_ccf_max_lag_s"),
                         (439, "out_ccf_decim")):
            self.assertRegex(TAGS, rf"SNTags\({idx}\)%Label\s*/\s*'{key}'")
        self.assertRegex(TAGS, r"SCTags\(102\)%Label\s*/\s*'tlag_cum_segments'")

    def test_per_gas_keys_are_appended(self):
        self.assertIn('+ ["tlag_drift_ref", "out_ccf"]', GEN)
        self.assertIn("rpGasLeapN + 26", READER)
        self.assertIn("rpGasLeapN + 27", READER)

    def test_method_six(self):
        self.assertRegex(READER, r"case \('6'\)\s*\n\s*Meth%tlag = 'cumcov'")
        self.assertIn("CumCovGenerate = Meth%tlag == 'cumcov'", READER)

    def test_innflux_defaults(self):
        for literal in ("RPSetup%tlag_cum_filter_len = 15",
                        "RPSetup%tlag_cum_min_ustar = 0.3d0",
                        "RPSetup%tlag_cum_min_compl_sonic = 0.95d0",
                        "RPSetup%tlag_cum_min_compl_gas = 0.5d0",
                        "RPSetup%tlag_cum_max_sst = 0.3d0",
                        "RPSetup%tlag_drift_meth = 0",
                        "RPSetup%tlag_drift_span_d = 5d0",
                        "RPSetup%tlag_drift_jump_thr = -0.8d0"):
            self.assertIn(literal, READER)

    def test_fluxnet_method_code(self):
        self.assertRegex(FLX, r"case\('cumcov'\)\s*\n\s*call AddIntDatumToDataline\(6,")


class ThePrepass(unittest.TestCase):
    def test_gathered_before_the_lag_is_applied(self):
        g = MAIN.index("call CumCovGatherPeriod(E2Set,")
        t = MAIN.index("call TimeLagHandle('maxcov', E2Set,", g)
        self.assertLess(t - g, 400)

    def test_flushed_after_ustar(self):
        a = MAIN.index("Stats7 = Stats")
        b = MAIN.index("call CumCovFlushPeriod(", a)
        self.assertLess(b - a, 600)

    def test_slice_closed_before_a_worker_stops(self):
        a = MAIN.index("end do to_periods_loop")
        b = MAIN.index("if (CumCovGenerate) call CumCovEndSlice()")
        stop = MAIN.index("call FinishBatchWorker()")
        self.assertLess(a, b)
        self.assertLess(b, stop)

    def test_resolved_by_the_parent_over_every_piece(self):
        self.assertIn("call CumCovResolve(PrepassChunkCount())", MAIN)

    def test_openness_is_a_flag_not_a_sign(self):
        self.assertNotRegex(CUMCOV, r"uScratch\s*[<>]=?\s*0")
        self.assertNotRegex(CUMCOV, r"usum\s*[<>]=?\s*0")
        self.assertIn("scratchOpen = ios == 0", CUMCOV)

    def test_records_sorted_by_period_then_gas(self):
        self.assertIn("if (recs(order(k))%ord < recs(j)%ord) exit", CUMCOV)

    def test_innflux_filters(self):
        for literal in ("r%ustar <= RPSetup%tlag_cum_min_ustar",
                        "r%csonic <= RPSetup%tlag_cum_min_compl_sonic",
                        "r%cgas <= RPSetup%tlag_cum_min_compl_gas",
                        "r%sst >= RPSetup%tlag_cum_max_sst"):
            self.assertIn(literal, CUMCOV)

    def test_summed_absolute_and_positive_peak(self):
        self.assertIn("abs(dble(recs(j)%ccf(L + d)))", CUMCOV)
        self.assertIn("PeakPositive, .false., error, lag2, found", CUMCOV)


class Production(unittest.TestCase):
    def test_applied_by_the_periods_own_stamp(self):
        self.assertIn("if (CumCovGenerate) call CumCovPeriodLags(Stats%date, Stats%time)",
                      MAIN)

    def test_segments_by_start_time(self):
        self.assertIn("- dble(RPsetup%avrg_len) / 1440d0", CUMCOV)

    def test_lag_handler_case(self):
        m = re.search(r"case \('cumcov'\)(.*?)case \('maxcov', 'maxcov&default'\)", TL, re.S)
        self.assertIsNotNone(m)
        self.assertIn("RowLags(j) = nint(CumCovLagNow(j) * Metadata%ac_freq)", m.group(1))
        self.assertIn("DefTlagUsed(j) = .false.", m.group(1))

    def test_scratch_and_table_in_the_parents_folder(self):
        """A worker's TmpDir is its own; the parent's is BatchTmpDir."""
        self.assertIn("dir = trim(BatchTmpDir) // slash", CUMCOV)
        self.assertIn("path = trim(SharedDir()) // 'cumcov_b'", CUMCOV)
        self.assertIn("path = trim(SharedDir()) // 'cumcov_lags.bin'", CUMCOV)

    def test_workers_read_the_table(self):
        self.assertIn("if (.not. tableReady) call ReadTable()", CUMCOV)


class Drift(unittest.TestCase):
    def test_relative_to_its_median(self):
        self.assertIn("perDrift(i, g2) = drift(i) - dmed", CUMCOV)

    def test_rlowess_like_matlab(self):
        self.assertIn("6d0 * mad", DRIFT)
        self.assertIn("(1d0 - rw(i) ** 2) ** 2", DRIFT)
        self.assertIn("(1d0 - wj ** 3) ** 3", DRIFT)

    def test_credited(self):
        self.assertIn("driftfix", DRIFT)
        self.assertIn("Cleveland 1979", DRIFT)
        self.assertIn("Striednig et al. 2020", CUMCOV_RAW)
        self.assertIn("Park et", CUMCOV_RAW)


class CcfOutput(unittest.TestCase):
    """Per-period cross-covariance functions, after innFLUX's step 1 files.
    Off unless a gas asks (gas_<i>_out_ccf); a worker's files go to the
    parent's folder like the cospectra."""

    def test_only_when_a_gas_asks(self):
        self.assertIn("if (any(RPSetup%out_ccf)) call WriteCcfOutput(E2Set,", MAIN)
        self.assertIn("if (.not. any(RPSetup%out_ccf)) return", CUMCOV)
        self.assertIn("RPSetup%out_ccf = .false.", READER)

    def test_before_the_lag_is_applied(self):
        a = MAIN.index("call WriteCcfOutput(E2Set,")
        b = MAIN.index("call TimeLagHandle(Meth%tlag(1:len_trim(Meth%tlag)), E2Set,")
        self.assertLess(a, b)

    def test_workers_write_to_the_parent_tree(self):
        par = code("src", "src_rp", "production_parallel.f90")
        self.assertIn("if (any(RPSetup%out_ccf)) call ToParentTree(CcfDir)", par)
        init = code("src", "src_rp", "init_outfiles_rp.f90")
        self.assertIn("CcfDir = Dir%main_out(1:len_trim(Dir%main_out)) // SubDirCcf // slash",
                      init)


if __name__ == "__main__":
    unittest.main()
