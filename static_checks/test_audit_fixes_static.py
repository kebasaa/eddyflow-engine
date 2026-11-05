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


class NoHygrometerMeansNoWaterColumns(unittest.TestCase):
    """PrimaryWaterOutSlot falls back to slot six, which holds another gas
    in a project without water: LE/ET were gated on it and ET was filled
    with that gas's flux times 0.0648."""

    WRITERS = [("src_fcc", "write_out_full_fcc.f90"), ("src_fcc", "init_out_files.f90"),
               ("src_fcc", "write_out_fluxnet_fcc.f90"), ("src_rp", "write_out_full.f90"),
               ("src_rp", "write_out_fluxnet.f90"), ("src_rp", "init_outfiles_rp.f90")]

    def test_writers_ask_the_real_water(self):
        for parts in self.WRITERS:
            src = code(*parts)
            self.assertNotIn("PrimaryWaterOutSlot()", src, parts[1])
            self.assertIn("PrimaryWater", src, parts[1])

    def test_the_helpers_decline_without_water(self):
        slots = code("src_common", "gas_slot_resolution.f90")
        for name in ("PrimaryWaterIn", "PrimaryWaterOf", "PrimaryWaterQC"):
            fn = body(slots, "function %s(" % name, "end function %s" % name)
            self.assertIn("if (wsl < firstGas .or. wsl > n) return", fn)

    def test_nothing_else_takes_slot_six_for_water(self):
        self.assertIn("HumiditySlot = gas == PrimaryWaterSlot()", code("src_rp", "pwb_timelag_handle.f90"))
        self.assertIn("if (wsl >= firstGas) E2Col(wsl)%present = .true.", code("src_rp", "optimize_timelags.f90"))
        self.assertIn("wsl = PrimaryWaterSlot()", code("src_rp", "add_to_timelag_opt_dataset.f90"))
        self.assertIn("'H2O RH classes: no hygrometer configured... NOT APPLIED'",
                      code("src_fcc", "spectral_assessment_diagnostics.f90"))


class PostFluxDespikingSeesThePublishedFluxes(unittest.TestCase):
    """It banked Level 2, so its series did not match the outputs."""

    MAIN = code("src_fcc", "eddyflow-fcc_main.f90")

    def test_level_three_is_banked(self):
        self.assertIn("pfd_nee, Flux3%H, Flux3%LE)", self.MAIN)
        self.assertNotIn("Flux2%gas(PrimaryCarbonSlot()), Flux2%H, Flux2%LE", self.MAIN)
        self.assertIn("if (pfd_csl >= firstGas) pfd_nee = Flux3%gas(pfd_csl)", self.MAIN)


class NoFootprintWithoutStability(unittest.TestCase):
    """Hsieh, Kormann-Meixner and Kljun ran on z/L or L = -9999."""

    FOOT = code("src_common", "footprint_handle.f90")

    def test_it_refuses_before_any_model(self):
        guard = "if (ustar == error .or. ustar <= 0d0 .or. zL == error &"
        self.assertIn(guard, self.FOOT)
        self.assertLess(self.FOOT.index(guard), self.FOOT.index("if (foot_model_used == 'kljun_04' .and. &"))


class AShortPathAveragesToOne(unittest.TestCase):
    """A tiny path length made the path-averaging function NaN, then the
    error code, which was multiplied into the band-pass function."""

    TF = code("src_common", "bpcf_analytic_transfer_functions.f90")
    AUX = code("src_common", "bpcf_aux_subs.f90")

    def test_the_limit_is_taken_below_the_cancellation(self):
        self.assertIn("where (2d0 * p * fp_irga(1:N) < 1d-7)", self.TF)
        self.assertIn("BPTF(1:N)%LP(var)%wirga = dsqrt(1d0 - 2d0 * p * fp_irga(1:N) / 6d0)", self.TF)

    def test_the_error_code_never_enters_the_product(self):
        bp = body(self.AUX, "subroutine BandPassTransferFunction", "end subroutine BandPassTransferFunction")
        self.assertIn("where (BPTF%LP(var1)%wirga == error .or. BPTF%LP(var2)%wirga == error)", bp)


class AbsoluteLimitsCountWhatTheyRemove(unittest.TestCase):

    def test_zero_only_when_nothing_is_removed(self):
        al = code("src_rp", "test_absolute_limits.f90")
        guard = al.index("if (.not. RPsetup%filter_al) then")
        self.assertIn("Essentials%al_s(u:lastGas) = 0", al[guard:guard + 200])
        self.assertNotIn("    if (RPsetup%filter_al) then\n        where (E2Col(u:lastGas)%present)", al)


class VpdIsConvertedOnItsOwnValue(unittest.TestCase):

    def test_the_guard_reads_vpd(self):
        rex = code("src_common", "read_ex_record.f90")
        self.assertIn("if (lEx%VPD /= error) lEx%VPD = lEx%VPD * 1d2", rex)
        self.assertNotIn("if (lEx%Pa /= error) lEx%VPD", rex)


class GhgRatesAreSetByTheFirstFileRead(unittest.TestCase):

    def test_the_reference_follows_the_first_read_file(self):
        imp = code("src_rp", "import_current_period.f90")
        self.assertIn("FirstRead = .false.\n    rawfile_loop: do", imp)
        self.assertIn("if (.not. FirstRead) then\n                        rate_changed = .true.", imp)
        self.assertIn("if (.not. FirstRead) then\n                        PeriodInstrFreq = LocCol%instr%ac_freq", imp)
        self.assertIn("                FirstRead = .true.\n", imp)
        self.assertNotIn("if (CurrentFile == FirstFile) then\n                        PeriodInstrFreq", imp)


if __name__ == "__main__":
    unittest.main()
