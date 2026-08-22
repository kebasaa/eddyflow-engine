"""Lenschow et al. (2000) instrument noise, and what it must not be confused with.

White instrument noise is uncorrelated between samples, so it lands entirely in
the autocovariance at zero lag. Fit a line through the first few lags,
extrapolate to zero, and the gap is the noise variance
(``EC_Software_Preproc/EddyUH_unc_Preproc.m:157-185``).

The whole risk in this method is that it looks like the three arms beside it
and answers a different question. It is the **analyser's own noise**, not a
sampling error (Finkelstein & Sims, Mann & Lenschow) and not a resolvability
floor (Billesbach), and it is systematically the smallest of the four. Three
things guard against them being blurred:

1. ``ru_meth = 5``, because 4 is Billesbach - the row and the method number
   stopped running in step when 3 became Mahrt.
2. The method name is ``lenschow_00``, distinct from the pre-existing
   ``mann_lenschow_94``. They share a surname and nothing else.
3. It **declines** a period rather than reporting a small number when its
   assumption fails.
"""

import re
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
HANDLE = ROOT / "src" / "src_rp" / "random_error_handle.f90"
DECODER = ROOT / "src" / "src_common" / "write_processing_project_variables.f90"
EDDYUH = (ROOT.parent / "EddyUH_testing" / "EddyUH" / "EddyUH_1.7b_COS"
          / "EC_Software_Preproc" / "EddyUH_unc_Preproc.m")


def read(path):
    return path.read_text(encoding="utf-8", errors="replace")


SRC = read(HANDLE)
#: The loop that fills Essentials, and the estimator it calls. Separated so
#: the conditional time-lag borrowing can reach the same number long before
#: the random-uncertainty stage runs - see the function's own note.
LOOP = SRC[SRC.index("subroutine RU_Lenschow_00("):
           SRC.index("end subroutine RU_Lenschow_00")]
ROUTINE = SRC[SRC.index("function LenschowFluxNoise(Set, N, M, gas)"):
              SRC.index("end function LenschowFluxNoise")]


def code(text):
    """Comments stripped - the house style names in prose what it forbids."""
    return re.sub(r"^\s*!.*$", "", text, flags=re.M)


class TheMethodNumberDoesNotMoveTheOthers(unittest.TestCase):

    def test_five_decodes_to_lenschow_00(self):
        src = read(DECODER)
        block = src[src.index("RUsetup%meth = 'none'") - 2000:]
        block = block[:block.index("end select")]
        self.assertRegex(block, r"case\(5\)[\s\S]{0,400}?'lenschow_00'")

    def test_the_four_that_were_there_keep_their_numbers(self):
        src = read(DECODER)
        for n, name in ((1, "finkelstein_sims_01"), (2, "mann_lenschow_94"),
                        (3, "mahrt_98"), (4, "billesbach_11")):
            self.assertRegex(
                src, r"case\(%d\)[\s\S]{0,400}?RUsetup%%meth = '%s'" % (n, name),
                "ru_meth %d no longer decodes to %s" % (n, name))

    def test_the_dispatch_has_an_arm_for_it(self):
        block = SRC[SRC.index("select case (RUsetup%meth)"):]
        block = block[:block.index("case default")]
        self.assertRegex(
            block, r"case\('lenschow_00'\)[\s\S]{0,400}?call RU_Lenschow_00")

    def test_it_is_not_the_1994_paper(self):
        #> mann_lenschow_94 is a sampling error from a different paper. The
        #> shared surname is the only thing they have in common, and a reader
        #> skimming a select case is exactly who would merge them.
        self.assertIn("RU_Mann_Lenschow_04", SRC)
        self.assertIn("RU_Lenschow_00", SRC)
        self.assertNotIn("mann_lenschow_00", SRC)
        self.assertNotIn("lenschow_94", code(ROUTINE))

    def test_it_does_not_ask_for_an_integral_turbulence_scale(self):
        #> The two arms that integrate over lags need it; this fits a fixed
        #> five and does not. Computing it anyway would be wasted work that
        #> looks load-bearing.
        block = SRC[SRC.index("case('lenschow_00')"):]
        block = block[:block.index("case('billesbach_11')")]
        self.assertNotIn("IntegralTurbulenceScale", block)


class TheEstimateItself(unittest.TestCase):

    def test_the_window_is_five_lags_in_samples(self):
        self.assertIn("integer, parameter :: first_lag = 1", ROUTINE)
        self.assertIn("integer, parameter :: last_lag = 5", ROUTINE)
        #> In SAMPLES, not seconds. Deriving it from the acquisition rate
        #> would be a different method wearing the same citation, so the
        #> reason is written where the constants are.
        self.assertIn("SAMPLES", ROUTINE)

    def test_the_noise_is_the_gap_at_lag_zero(self):
        self.assertIn("w_noise = acov_w(0) - InterceptAtZero(acov_w)", ROUTINE)
        self.assertIn("g_noise = acov_g(0) - InterceptAtZero(acov_g)", ROUTINE)

    def test_the_flux_noise_uses_the_total_variance_of_w(self):
        #> Not w's noise variance. EddyUH's form: the scalar's own noise
        #> beaten against the full vertical wind signal. Using w_noise here
        #> would be a plausible-looking and different quantity.
        self.assertIn(
            "LenschowFluxNoise = dsqrt(g_noise * acov_w(0) / dble(N))", ROUTINE)
        tail = ROUTINE[ROUTINE.index("LenschowFluxNoise = dsqrt"):]
        self.assertNotIn("w_noise", tail[:120])

    def test_the_intercept_is_a_closed_form_not_a_fit_call(self):
        helper = ROUTINE[ROUTINE.index("function InterceptAtZero"):]
        self.assertIn("slope = sxy / sxx", helper)
        self.assertIn("InterceptAtZero = ybar - slope * xbar", helper)

    def test_the_estimator_is_the_only_implementation(self):
        #> The loop that fills Essentials and the borrowing that judges a
        #> covariance must divide by the same number. Two copies of this
        #> arithmetic would drift, and the drift would show up as a lag
        #> population rather than as an uncertainty column.
        self.assertIn("LenschowFluxNoise(Set, N, M, var)", LOOP)
        self.assertNotIn("InterceptAtZero", LOOP)
        self.assertIn("LaggedCovarianceNoError(Set(:, w), Set(:, w)", ROUTINE)

    def test_it_refuses_a_column_that_is_not_there(self):
        #> Called per gas from two places now, one of which does not pre-check.
        self.assertIn("if (.not. E2Col(gas)%present) return", ROUTINE)


class ItDeclinesRatherThanGuesses(unittest.TestCase):

    def test_a_non_positive_intercept_rejects_the_period(self):
        self.assertIn("if (w_noise <= 0d0 .or. g_noise <= 0d0) return", ROUTINE)

    def test_the_vertical_wind_gates_every_gas(self):
        #> EddyUH rejects on either intercept. A w autocovariance with no
        #> noise-like step means the assumption fails for the period, whatever
        #> the scalar did - so w_noise is in the same condition as g_noise and
        #> not merely reported.
        cond = ROUTINE[ROUTINE.index("if (w_noise"):]
        self.assertIn("w_noise <= 0d0", cond[:60])
        #> And an error-coded autocovariance on either side refuses too.
        self.assertIn(
            "if (any(acov_w == error) .or. any(acov_g == error)) return",
            ROUTINE)

    def test_the_declined_value_is_the_error_code_not_a_zero(self):
        #> Set before every early return, so a refusal reads as "not
        #> available" rather than as a very small noise estimate.
        head = ROUTINE[:ROUTINE.index("do lag = 0, last_lag")]
        self.assertIn("LenschowFluxNoise = error", head)
        #> And the loop clears the slot before asking, so a gas that is not
        #> there does not keep the previous one's answer.
        self.assertIn("Essentials%rand_uncer(var) = error", LOOP)


class ItMatchesEddyUH(unittest.TestCase):

    @unittest.skipUnless(EDDYUH.is_file(), "EddyUH tree not beside this one")
    def test_eddyuh_still_fits_lags_one_to_five(self):
        src = read(EDDYUH)
        self.assertIn("polyfit([1:5]',autocov_w(lagsapu<=5 & lagsapu>=1),1)", src)

    @unittest.skipUnless(EDDYUH.is_file(), "EddyUH tree not beside this one")
    def test_eddyuh_still_gates_on_both_intercepts(self):
        src = read(EDDYUH)
        self.assertIn("if wwerr==0 | ggerr==0", src)

    @unittest.skipUnless(EDDYUH.is_file(), "EddyUH tree not beside this one")
    def test_eddyuh_still_uses_the_total_w_variance(self):
        src = read(EDDYUH)
        self.assertIn("(1/N)^0.5.*(ggerr*autocov_w(lagsapu==0))^0.5", src)


if __name__ == "__main__":
    unittest.main()
