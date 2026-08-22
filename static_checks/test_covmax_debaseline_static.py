"""Static checks for baseline-subtracted time-lag selection.

Instead of the largest absolute cross-covariance in the search window, take
the largest departure from the straight line joining the window's two ends. A
weak flux often sits on a sloping cross-covariance - from a trend, or from a
neighbouring stronger correlation - and the plain maximum then lands on
whichever end the slope is highest at rather than on the peak.

EddyUH's ``EddyUH_SC_Flux.m:281``. Worth knowing that in the shipped EddyUH
this is what the *"standard"* menu entry runs: the plain ``max|cov|`` line is
commented out immediately above it.

Three things are pinned.

**The default path must be untouched.** Selecting the lag by the largest
absolute covariance is what every existing project does, and rewriting
``CovMax`` into two passes had to leave it computing the same numbers in the
same order.

**Two latent faults were fixed in passing and must stay fixed.** ``RLag`` was
``intent(out)`` and assigned only inside the comparison, so a window in which
nothing beat the initial zero returned an undefined row lag that the caller
then shifted the series by. And an error-coded covariance was compared like
any other - the code being -9999, its magnitude beat every real covariance in
the window, so a lag at which the two series shared no valid sample won
outright.

**It is a modifier, not a method.** A new ``tlag_meth`` value would force a
decision in the EddyPro importer, which must keep mapping EddyPro's
``tlag_meth=2`` onto plain ``maxcov``.
"""

import re
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / "src"


def read(rel):
    return (SRC / rel).read_text(encoding="utf-8", errors="replace")


TLAG = read("src_rp/timelag_handle.f90")
READER = read("src_rp/read_ini_rp.f90")
TAGS = read("src_rp/m_rp_global_var.f90")
TYPEDEF = read("src_common/m_typedef.f90")
IMPORT = read("src_common/m_eddypro_import.f90")

#: The slot it was given, out of the blanks at 79-82.
SLOT = 79


def covmax_body():
    start = TLAG.index("subroutine CovMax(")
    return TLAG[start:TLAG.index("end subroutine CovMax", start)]


class TheSettingReachesTheEngine(unittest.TestCase):

    def test_the_key_holds_its_slot(self):
        self.assertRegex(
            TAGS, r"SCTags\(%d\)%%Label\s*/\s*'covmax_debaseline'\s*/" % SLOT)

    def test_it_is_read_under_its_found_guard(self):
        #> Unlike covmax_var and covmax_stocdet beside it, which every project
        #> states. A new character key must not be read blind: the tag array
        #> is saved between parses.
        self.assertRegex(
            READER,
            r"RPSetup%%covmax_debaseline = SCTagFound\(%d\)" % SLOT,
            "covmax_debaseline is read without its found guard")

    def test_it_defaults_to_off(self):
        #> The guard is the default: absent tag, false.
        self.assertNotRegex(
            READER,
            r"RPSetup%covmax_debaseline\s*=\s*\.true\.",
            "the modifier is being switched on somewhere")

    def test_the_setup_type_carries_it(self):
        self.assertRegex(TYPEDEF, r"logical :: covmax_debaseline")


class ItStaysAModifierNotAMethod(unittest.TestCase):

    def test_no_new_timelag_method_was_added(self):
        #> EddyPro has four methods and the importer maps its tlag_meth=2 onto
        #> maxcov. A fifth EddyFlow method is pwb; a sixth would need an
        #> importer decision that nothing has made.
        methods = re.findall(r"Meth%tlag = '([a-z_&]+)'", READER)
        self.assertEqual(
            sorted(set(methods)),
            sorted({"none", "constant", "maxcov&default", "maxcov",
                    "tlag_opt", "pwb"}),
            "the set of time-lag methods changed: %s" % sorted(set(methods)))

    def test_the_importer_does_not_write_it(self):
        #> An imported EddyPro project must behave as EddyPro did, and EddyPro
        #> has no such criterion.
        self.assertNotIn(
            "covmax_debaseline", IMPORT,
            "the EddyPro importer now states covmax_debaseline, which would "
            "change what a converted project computes")


class TheDefaultPathIsUnchanged(unittest.TestCase):

    def setUp(self):
        self.body = covmax_body()

    def test_the_plain_criterion_is_still_the_absolute_covariance(self):
        self.assertIn("score = dabs(CovSeries(k))", self.body)

    def test_ties_still_go_to_the_earliest_lag(self):
        #> Strictly greater, as before. Changing this to >= would silently
        #> move every lag that has a flat maximum.
        self.assertIn("if (score > MaxCov) then", self.body)

    def test_the_covariance_is_computed_once_per_lag_as_before(self):
        #> Pass one still allocates, aligns, optionally detrends and calls the
        #> covariance exactly as the single loop did, so the arithmetic and
        #> its floating-point association are untouched.
        self.assertIn("call CovarianceMatrixNoError(ShPrimes", self.body)
        self.assertIn("if (RPSetup%covmax_stocdet) then", self.body)
        self.assertEqual(
            self.body.count("call CovarianceMatrixNoError("), 1,
            "the covariance is computed more than once per lag")


class TheTwoLatentFaultsStayFixed(unittest.TestCase):

    def setUp(self):
        self.body = covmax_body()

    def test_the_row_lag_is_initialised(self):
        self.assertIn(
            "RLag = lagmin", self.body,
            "RLag is intent(out) and would be returned undefined when no lag "
            "beats the initial zero")

    def test_an_error_coded_covariance_is_skipped(self):
        self.assertIn(
            "if (CovSeries(k) == error) cycle", self.body,
            "an error-coded covariance is being compared again; its magnitude "
            "is 9999 and beats every real covariance in the window")

    def test_the_skip_precedes_the_comparison(self):
        skip = self.body.index("if (CovSeries(k) == error) cycle")
        cmp_ = self.body.index("if (score > MaxCov) then")
        self.assertLess(skip, cmp_)


class TheChordIsDrawnBetweenTheWindowEnds(unittest.TestCase):

    def setUp(self):
        self.body = covmax_body()

    def test_it_interpolates_between_the_first_and_last_lag(self):
        self.assertIn("baseline = CovSeries(1)", self.body)
        self.assertIn("(CovSeries(nlag) - CovSeries(1))", self.body)
        self.assertIn("dble(k - 1) / dble(nlag - 1)", self.body)

    def test_the_score_is_the_departure_from_it(self):
        self.assertIn("score = dabs(CovSeries(k) - baseline)", self.body)

    def test_it_declines_a_window_too_short_to_have_an_interior(self):
        #> Two points make a line through themselves and nothing else.
        self.assertIn("nlag >= 3", self.body)

    def test_it_declines_when_either_end_is_missing(self):
        self.assertIn(
            "CovSeries(1) /= error .and. CovSeries(nlag) /= error", self.body,
            "the chord is drawn from an error-coded endpoint")

    def test_the_whole_window_is_kept_before_anything_is_chosen(self):
        #> The reason for two passes at all: by the time the far end has been
        #> computed, a running scalar has lost the near one.
        fill = self.body.index("CovSeries(i - lagmin + 1) = Cov")
        choose = self.body.index("do k = 1, nlag")
        self.assertLess(fill, choose)


if __name__ == "__main__":
    unittest.main()
