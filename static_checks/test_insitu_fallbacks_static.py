"""What the in-situ spectral corrections do when part of them is missing.

Five ways a period's in-situ correction can lack an ingredient, and what each
must lead to. Each of these once led somewhere else, silently:

1. A stalled Ibrom et al. (2007) fit. The solver is judged not to have moved
   when it still sits at its starting point; the stable fit was compared with
   (1, 0.02), a point it never starts from, so a stalled one passed its own
   starting guess off as a result. Start and check now read one constant.

2. A failed fit. Its parameters are the error code, and the Ibrom formula
   turned -9999 into a factor of about WS + 1. Now it yields no factor, and the
   hub corrects that gas with the analytic method, as it does any gas the
   in-situ method cannot take.

3. A Fratini et al. (2012) period with no full w/T cospectrum. Every factor
   was set to 1, Tau's and H's included, which the hub had already computed.
   Now the gases take the method's own Ibrom 2007 fallback and Tau and H are
   left alone.

4. A second hygrometer whose fit is its own while the primary's RH regression
   failed. The gate and the readiness report asked the primary's coefficients
   of every hygrometer, while the correction evaluated each one's own. All
   three now ask one rule.

5. A period with no cut-off frequency. A hygrometer with no humidity had its
   RH regression evaluated at the error code, and a gas with no cut-off kept
   the factor of 1 its caller starts from. Now it has no factor, the hub
   corrects it analytically, and the FCC run ends by listing every such
   period, gas by gas, with the reason.

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


MODEL = code("src_fcc", "correction_factor_model.f90")
AUX = code("src_common", "bpcf_aux_subs.f90")
IBROM = code("src_common", "bpcf_Ibrom_07.f90")
FRATINI = code("src_common", "bpcf_fratini_12.f90")
HUB = code("src_common", "bpcf_bandpass_spectral_corrections.f90")
SLOTS = code("src_common", "gas_slot_resolution.f90")
REPORT = code("src_fcc", "spectral_assessment_diagnostics.f90")
HORST = code("src_common", "bpcf_Horst_97.f90")
PARALLEL = code("src_fcc", "fcc_parallel.f90")
MAIN = code("src_fcc", "eddyflow-fcc_main.f90")


class AStalledFitIsCaught(unittest.TestCase):

    def test_start_and_check_read_one_constant(self):
        self.assertIn("real(kind = dbl), parameter :: UnStart(2) = (/ 2.5d-3, 6.0d-4 /)", MODEL)
        self.assertIn("real(kind = dbl), parameter :: StStart(2) = (/ 6.5d-3, 3.0d-3 /)", MODEL)
        self.assertIn("UnPar(1:2) = UnStart", MODEL)
        self.assertIn("StPar(1:2) = StStart", MODEL)
        self.assertIn("if (all(UnPar(1:2) == UnStart) .or. info < 1 .or. info > 3) then", MODEL)
        self.assertIn("if (all(StPar(1:2) == StStart) .or. info < 1 .or. info > 3) then", MODEL)

    def test_no_check_against_a_point_the_fit_never_starts_from(self):
        self.assertNotIn("StPar(1) == 1d0", MODEL)
        self.assertNotIn("0.02d0", MODEL)


class AFailedFitGivesNoFactor(unittest.TestCase):

    def test_the_formula_refuses_error_parameters(self):
        fn = body(AUX, "subroutine CorrectionFactorsIbrom07", "end subroutine CorrectionFactorsIbrom07")
        for par in ("StPar", "UnPar"):
            self.assertIn("if (any(%s(1:2) == error)) then" % par, fn)
        # One for each stability, and one for a missing cut-off (case 5).
        self.assertEqual(fn.count("lBPCF%of(gas) = error"), 3)
        # The guard comes before the formula it protects.
        self.assertLess(fn.index("if (any(StPar(1:2) == error)) then"),
                        fn.index("lBPCF%of(gas) = StPar(1) * lEx%WS"))
        self.assertLess(fn.index("if (any(UnPar(1:2) == error)) then"),
                        fn.index("lBPCF%of(gas) = UnPar(1) * lEx%WS"))

    def test_ibrom_carries_the_missing_factor(self):
        # Not the low-frequency part alone, which would read as corrected.
        self.assertIn("where (BPCF_Ibrom%of(firstGas:lastGas) == error)\n"
                      "        BPCF%of(firstGas:lastGas) = error", IBROM)

    def test_the_separation_term_leaves_a_missing_factor_missing(self):
        self.assertEqual(HUB.count(".and. BPCF%of(firstGas:lastGas) /= error)"), 3)

    def test_the_hub_corrects_such_a_gas_analytically(self):
        route = HUB.index("if (BPCF%of(gas) /= error) cycle")
        moncrieff = HUB.index("if (n_analytic > 0) then")
        self.assertLess(route, moncrieff)
        block = HUB[route - 400:route + 200]
        self.assertIn("if (.not. insitu_ok(gas)) cycle", block)
        self.assertIn("insitu_ok(gas) = .false.", block)
        self.assertIn("n_analytic = n_analytic + 1", block)


class AFratiniGapKeepsTauAndH(unittest.TestCase):

    def test_no_factor_is_reset_to_one(self):
        self.assertNotIn("BPCF%of(:) = 1d0", FRATINI)

    def test_the_gases_take_the_ibrom_fallback(self):
        gap = FRATINI[FRATINI.rindex("    else\n", 0, FRATINI.index("contains")):
                      FRATINI.index("contains")]
        self.assertIn("call RetrieveLPTFpars(lEx, 'iir', LocSetup)", gap)
        self.assertIn("if (loc_var_present(gas)) call CorrectionFactorsIbrom07(gas, BPCF, lEx)", gap)
        self.assertNotIn("w_u", gap)
        self.assertNotIn("w_ts", gap)


class EveryHygrometerAsksOneRule(unittest.TestCase):

    def test_the_rule(self):
        rule = body(SLOTS, "subroutine HygrometerRhCoefficients", "end subroutine HygrometerRhCoefficients")
        self.assertLess(rule.index("A = RegPar(gas_slot, dum)%e1"),
                        rule.index("A = RegPar(dum, dum)%e1"))
        self.assertIn("logical function HygrometerHasRhFit(gas_slot)", SLOTS)

    def test_the_correction_the_gate_and_the_report_all_ask_it(self):
        self.assertIn("call HygrometerRhCoefficients(gas, A, B, C)", AUX)
        # Nested rather than one .and., which gfortran warns need not evaluate
        # the function.
        self.assertIn("if (insitu_ok(gas)) insitu_ok(gas) = HygrometerHasRhFit(gas)", HUB)
        self.assertEqual(REPORT.count("if (gas_ok) gas_ok = HygrometerHasRhFit(gas)"), 2)

    def test_nothing_asks_the_primary_for_every_hygrometer(self):
        self.assertNotIn("rh_ok", HUB)
        self.assertNotIn("h2o_ok", REPORT)


class AMissingCutoffGivesNoFactor(unittest.TestCase):
    """5. A period with no cut-off frequency, such as a hygrometer with no
    humidity to evaluate its RH regression at. The regression was evaluated
    at the error code, and the formulas left the 1 their callers start from;
    either way the period read as corrected. Now there is no factor, the hub
    corrects the gas analytically, and the run ends by naming the periods."""

    def test_no_humidity_no_cutoff(self):
        arm = body(AUX, "case('iir')", "case('sigma')")
        self.assertLess(arm.index("if (lRH == error) lRH = lEx%RH"),
                        arm.index("if (lRH == error) cycle"))
        self.assertLess(arm.index("if (lRH == error) cycle"),
                        arm.index("f_c(gas) = dexp(A * lRH**2 + B * lRH + C)"))
        self.assertIn("if (f_c(gas) <= 0d0 .or. f_c(gas) > huge(1d0)) f_c(gas) = error", arm)

    def test_both_formulas_give_no_factor(self):
        ibrom = body(AUX, "subroutine CorrectionFactorsIbrom07", "end subroutine CorrectionFactorsIbrom07")
        horst = body(AUX, "subroutine CorrectionFactorsHorst97", "end subroutine CorrectionFactorsHorst97")
        for fn in (ibrom, horst):
            self.assertIn("if (f_c(gas) == error .or. f_c(gas) <= 0d0) then", fn)
        self.assertNotIn("if (f_c(gas) == error) return", ibrom)
        self.assertEqual(ibrom.count("lBPCF%of(gas) = error"), 3)
        self.assertIn("lBPCF%of(gas) = error\n            cycle", horst)

    def test_horst_carries_the_missing_factor(self):
        self.assertIn("where (BPCF_Horst%of(firstGas:lastGas) == error)\n"
                      "        BPCF%of(firstGas:lastGas) = error", HORST)

    def test_fratini_keeps_no_direct_factor_without_a_cutoff(self):
        loop = FRATINI[FRATINI.index("if (f_c(gas) == error) then"):]
        self.assertTrue(loop.startswith("if (f_c(gas) == error) then\n"
                                        "                BPCF%of(gas) = error\n"
                                        "                cycle"))
        self.assertLess(FRATINI.index("if (f_c(gas) == error) then"),
                        FRATINI.index("BPCF%of(gas) <= min_bpcf_f12(gas)"))

    def test_the_hub_notes_each_period_with_its_reason(self):
        block = HUB[HUB.index("if (BPCF%of(gas) /= error) cycle"):HUB.index("if (n_analytic > 0) then")]
        self.assertIn("call StoreInsituFallback(gas, lEx%end_date, lEx%end_time, InsituNoCutoff)", block)
        self.assertIn("call StoreInsituFallback(gas, lEx%end_date, lEx%end_time, InsituFitFailed)", block)

    def test_workers_hand_their_periods_back(self):
        finish = body(PARALLEL, "subroutine FinishFccWorker", "end subroutine FinishFccWorker")
        merge = body(PARALLEL, "subroutine MergeFccPieces", "end subroutine MergeFccPieces")
        begins = body(PARALLEL, "subroutine FccPieceBegins", "end subroutine FccPieceBegins")
        self.assertIn("InsituFallbackN = 0", begins)
        # Written and read in the same place: after the despiking rows.
        self.assertLess(finish.index("write(ud) PfdCacheN"), finish.index("write(ud) InsituFallbackN"))
        self.assertIn("write(ud) InsituFallbacks(1:InsituFallbackN)", finish)
        self.assertLess(merge.index("call StorePfdCache"), merge.index("read(ud) fallbacks"))
        self.assertIn("call StoreInsituFallback(fallbacks(i)%gas, fallbacks(i)%date, &", merge)

    def test_the_run_ends_by_naming_them(self):
        self.assertLess(MAIN.index("call MergeFccPieces()"), MAIN.index("call ReportInsituFallbacks()"))
        report = body(REPORT, "subroutine ReportInsituFallbacks", "end subroutine ReportInsituFallbacks")
        self.assertIn("if (InsituFallbackN < 1) return", report)
        self.assertIn("call Put(InsituFallbacks(i)%date // ' ' // InsituFallbacks(i)%time)", report)


if __name__ == "__main__":
    unittest.main()
