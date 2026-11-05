"""What the in-situ spectral corrections do when part of them is missing.

Four ways a period's in-situ correction can lack an ingredient, and what each
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
        self.assertEqual(fn.count("lBPCF%of(gas) = error"), 2)
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


if __name__ == "__main__":
    unittest.main()
