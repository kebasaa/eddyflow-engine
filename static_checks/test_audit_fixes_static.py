"""Defects found by the 8.1.8 audit of the RP and FCC chains, one class each.

Each of these gave a wrong number or flag without saying so, in a
configuration people use. The classes pin the corrected code shapes, so a
later edit that brings one back fails here rather than in someone's data.

Part of the EddyFlow engine's static checks.
"""

from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"


def code(*parts):
    text = SRC.joinpath(*parts).read_text(encoding="utf-8", errors="replace")
    return "\n".join(ln for ln in text.splitlines()
                     if not ln.lstrip().startswith("!"))


def body(src, start, end):
    return src[src.index(start):src.index(end)]


class MauderDespikingDrawsFromGoodSamples(unittest.TestCase):
    """Mauder (2013) drew its interpolation from the row before the run,
    which could be the error code: a spike beside one filtered sample became
    about (x - 9999) / 2."""

    MAUDER = code("src_rp", "test_spike_detection_mauder_13.f90")

    def test_no_line_from_the_row_before_the_run(self):
        self.assertNotIn("Set(i - (cnt + 1), j)", self.MAUDER)
        self.assertNotIn("IsSpike(i-k, j) = .true.", self.MAUDER)

    def test_the_line_starts_at_the_last_good_sample(self):
        self.assertIn("m = (Set(i, j) - Set(prev_good, j)) / dble(i - prev_good)", self.MAUDER)
        self.assertIn("if (prev_good /= 0) then", self.MAUDER)
        self.assertIn("Set(k, j) = m * dble(k - prev_good) + q", self.MAUDER)
        # Both arms track the good samples and flag samples only.
        self.assertEqual(self.MAUDER.count("prev_good = i"), 2)
        self.assertEqual(self.MAUDER.count("if (Set(k, j) /= error) IsSpike(k, j) = .true."), 2)
        self.assertEqual(self.MAUDER.count("new_spike = .not. IsSpike(run_last, j)"), 2)


class TheIterationCompoundsNothing(unittest.TestCase):
    """The iterative correction calls Level 1 once per pass. RP's u* and
    both programs' oxygen correction scaled their inputs in place, so they
    were applied once per pass."""

    RP1 = code("src_rp", "fluxes1_rp.f90")
    FCC1 = code("src_fcc", "fluxes1.f90")
    RPMAIN = code("src_rp", "eddyflow-rp_main.f90")
    FCCMAIN = code("src_fcc", "eddyflow-fcc_main.f90")

    def test_ustar_is_built_from_level_zero(self):
        self.assertIn("Ambient%us = Flux0%ustar * dsqrt(BPCF%of(w_u))", self.RP1)
        self.assertNotIn("Ambient%us = Ambient%us * dsqrt", self.RP1)

    def test_the_oxygen_correction_is_its_own_routine(self):
        for src, name, level1 in ((self.RP1, "OxygenCorrectionKrypton_rp", "Fluxes1_rp"),
                                  (self.FCC1, "OxygenCorrectionKrypton", "Fluxes1")):
            l1 = body(src, "subroutine %s(" % level1, "end subroutine %s\n" % level1)
            self.assertNotIn("Cox", l1)
            self.assertIn("subroutine %s(" % name, src)

    def test_it_runs_once_before_the_passes(self):
        fcc = self.FCCMAIN
        self.assertLess(fcc.index("call OxygenCorrectionKrypton(lEx)"),
                        fcc.index("do corr_pass = 1, corr_passes"))
        rp = self.RPMAIN
        self.assertLess(rp.index("call OxygenCorrectionKrypton_rp("),
                        rp.index("do corr_pass = 1, corr_passes"))
        # And ahead of the FCC-only branch's single Level 1.
        self.assertEqual(rp.count("call OxygenCorrectionKrypton_rp("), 2)


class TheOxygenCorrectionHasItsHeatCapacityAndReachesTheFlux(unittest.TestCase):
    """Cox = 1 + 0.23 (ko/kw) beta lambda / (cp T). The 1/cp was missing, and
    the correction scaled the covariance and variance columns only, never
    FH2O, LE or ET."""

    RP1 = code("src_rp", "fluxes1_rp.f90")
    FCC1 = code("src_fcc", "fluxes1.f90")
    RPMAIN = code("src_rp", "eddyflow-rp_main.f90")

    def test_the_formula_divides_by_cp(self):
        self.assertIn("/ (Ambient%RhoCp / RHO%a * Ambient%Ta)", self.RP1)
        self.assertIn("* lEx%Bowen * lEx%lambda / (lEx%Cp * lEx%Ta)", self.FCC1)
        self.assertNotIn("* Ambient%Bowen * Ambient%lambda / Ambient%Ta", self.RP1)
        self.assertNotIn("* lEx%Bowen * lEx%lambda / lEx%Ta", self.FCC1)

    def test_the_water_flux_is_corrected(self):
        fcc = body(self.FCC1, "subroutine OxygenCorrectionKrypton(", "end subroutine OxygenCorrectionKrypton")
        self.assertIn("lEx%Flux0%gas(gas) = Cox * lEx%Flux0%gas(gas)", fcc)
        self.assertIn("lEx%Flux0%LE = Cox * lEx%Flux0%LE", fcc)
        rp = body(self.RP1, "subroutine OxygenCorrectionKrypton_rp(", "end subroutine OxygenCorrectionKrypton_rp")
        self.assertIn("if (scale_fluxes) then", rp)
        self.assertIn("Flux0%gas(msl) = Cox * Flux0%gas(msl)", rp)
        self.assertIn("if (Flux0%LE /= error) Flux0%LE = Cox * Flux0%LE", rp)

    def test_once_only_when_fcc_follows(self):
        # RP corrects the flux when it is the last word on it; when FCC
        # follows, the ex record carries Level 0 uncorrected for FCC.
        self.assertIn("call OxygenCorrectionKrypton_rp(.true.)\n                do corr_pass", self.RPMAIN)
        self.assertIn("call OxygenCorrectionKrypton_rp(.false.)\n                    call Fluxes1_rp()", self.RPMAIN)


class DetrendingRunsOnRowsNotSamples(unittest.TestCase):
    """A column with missing rows - slower than the file, or filtered -
    was fitted on sample counts and evaluated on rows (linear), updated at
    the wrong weight (running mean), or decayed one step per sample
    (exponential). The latter two also started from uninitialised memory
    after a leading gap."""

    FLUCT = code("src_common", "fluctuations.f90")

    def test_linear_regresses_on_the_row_offset(self):
        lin = body(self.FLUCT, "subroutine LinDetrend", "end subroutine LinDetrend")
        self.assertIn("sumx1(j) = sumx1(j) + (Set(i, j)*(dble(i - imin)))", lin)
        self.assertIn("(dble(i - imin) - sumtime(j) / dble(nn(j)))", lin)
        self.assertNotIn("dble(nn(j) - 1)", lin)
        self.assertNotIn("mm(j)", lin)

    def test_running_mean_averages_the_valid_samples_of_a_gappy_column(self):
        rm = body(self.FLUCT, "subroutine RunningMean", "end subroutine RunningMean")
        self.assertIn("if (any(Set(1:N, k) == error)) then", rm)
        self.assertIn("if (win_cnt > 0) Trend(i, k) = win_sum / dble(win_cnt)", rm)
        self.assertIn("Trend = error", rm)

    def test_exponential_decays_by_the_gap_and_has_a_start(self):
        ew = body(self.FLUCT, "subroutine ExpWeightAvrg", "end subroutine ExpWeightAvrg")
        self.assertIn("if (last_valid > 0) bk = b ** (i - last_valid)", ew)
        self.assertIn("Trend(i, k) = bk * Trend(i - back, k) &", ew)
        self.assertIn("if (Trend(i, k) == error) Trend(i, k) = Set(i, k)", ew)
        self.assertIn("Trend = error", ew)


class TheHaarTestComparesLikeWithLike(unittest.TestCase):
    """The half windows were cut before the molar-density gases were
    converted, and only the whole window's standard deviation was squared."""

    DS = code("src_rp", "test_discontinuities.f90")

    def test_conversion_precedes_the_half_windows(self):
        self.assertLess(self.DS.index("XX(1:nn, gas) = XX(1:nn, gas) * StdVair * 1d3"),
                        self.DS.index("XX_dw(i, :) = XX(i, :)"))
        self.assertIn("where (XX(1:nn, gas) /= error) XX(1:nn, gas) = XX(1:nn, gas) * StdVair", self.DS)

    def test_all_three_are_variances(self):
        for v in ("Var", "Var_dw", "Var_up"):
            self.assertIn("where (%s /= error) %s = %s**2" % (v, v, v), self.DS)


class PwbCachedLagsFollowThePeriodsRate(unittest.TestCase):
    """The cached row lags were converted once at the last pre-pass rate."""

    STREAM = code("src_rp", "pwb_stream.f90")

    def test_rows_come_from_seconds_at_this_rate(self):
        self.assertIn("RowLags(j) = nint(ev%cache_used(j) * ev%ac_freq)", self.STREAM)
        self.assertIn("pwb_raw_OwnRowLags(j) = RowLags(j)", self.STREAM)
        self.assertNotIn("pwb_raw_OwnRowLags(j) = ev%cache_row(j)", self.STREAM)


if __name__ == "__main__":
    unittest.main()
