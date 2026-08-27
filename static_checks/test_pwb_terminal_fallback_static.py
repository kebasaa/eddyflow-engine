"""The terminal fallback is the period's own covariance maximum.

When nothing at all reaches a period - no detection of its own, nothing to
interpolate between, nothing to carry, no same-analyser neighbour, not even a
per-gas median - `PostProcessPwbTimelagCache` ends with an arm that labels
itself `maxcov_default`. This file is about making that label true.

It used to hand back `fallback_lag`, captured from whatever the STREAMING pass
in timelag_handle had settled on. For a gas the rule rejected everywhere that
is the same number: the streaming pass never had a previous lag to carry
either, so it fell back to covariance maximisation too. The two part company
in one case - the HDI pre-filter runs in the post-pass and NOT in the streaming
pass, so a gas whose every detection is discarded there had settled and carried
during the streaming walk. What came back was a carried lag wearing the
`maxcov_default` label.

It shows up plainly on `base_pwb_prefilt`, which tightens the pre-filter to
0.10 s so every row takes this arm: before the fix, h2o came back with 18.8 s
for three consecutive periods, which is not something a per-period covariance
maximum can do. 10 of 21 rows changed.

The second reason is pass order. `pwb_last_optimal_lag` depends on where the
walk began, so a worker process starting cold carries a different lag and this
arm hands back a different answer. That is one of the values that has to stop
carrying pass order before the PWB pre-pass can be split across processes.

So the covariance maximum is now taken per period at detection time, whether
or not anything needs it, and stored on the row. That costs one CovMax per gas
per period on every PWB run - the price of the answer being a property of the
period rather than of the walk.

Part of the EddyFlow engine's static checks.
"""

from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]

MODULE = "src/src_rp/pwb_timelag_handle.f90"
HANDLER = "src/src_rp/timelag_handle.f90"
TYPES = "src/src_common/m_typedef.f90"


def read(path):
    return (ROOT / path).read_text(encoding="utf-8")


def code(path):
    """Source with full-line comments removed."""
    return (chr(10)).join(ln for ln in read(path).splitlines()
                          if not ln.lstrip().startswith("!"))


def body_of(source, opener, closer):
    return source[source.index(opener):source.index(closer)]


MOD = code(MODULE)
HAND = code(HANDLER)
TYP = code(TYPES)

POST = body_of(MOD, "subroutine PostProcessPwbTimelagCache",
               "end subroutine PostProcessPwbTimelagCache")


class ThePeriodCarriesItsOwnCovarianceMaximum(unittest.TestCase):

    def test_the_row_has_somewhere_to_put_it(self):
        self.assertIn("real(kind = dbl) :: maxcov_lag", TYP)

    def test_it_is_taken_at_detection_time(self):
        self.assertIn("lPwbResult%maxcov_lag = mc_used", HAND)

    def test_it_is_taken_for_every_period_not_only_the_ones_that_need_it(self):
        """Whether a period needs it cannot be known while the walk is still
        going - that is decided by the settled table afterwards. So it is taken
        straight after detection, before any classification has happened, and
        before the branch that used to be the only caller."""
        detect = HAND.index("call PwbDetectGas(Set, nrow, ncol, j, lPwbResult, pwb_success)")
        taken = HAND.index("lPwbResult%maxcov_lag = mc_used")
        classified = HAND.index("lPwbResult%reliability_class = 'S1_optimal'")
        self.assertLess(detect, taken)
        self.assertLess(taken, classified,
                        "the maximum must be taken before anything is classified")

    def test_it_does_not_depend_on_a_previous_lag(self):
        """pwb_has_previous is the streaming carry flag. If the new call sat
        under it, the value would carry pass order and nothing would be
        fixed."""
        seg = HAND[HAND.index("call PwbDetectGas(Set, nrow, ncol, j, lPwbResult, pwb_success)"):
                   HAND.index("lPwbResult%maxcov_lag = mc_used")]
        self.assertNotIn("pwb_has_previous", seg)

    def test_it_starts_unset(self):
        self.assertIn("res%maxcov_lag = error", MOD)


class TheTerminalArmUsesIt(unittest.TestCase):

    def test_step_eight_reads_the_row(self):
        self.assertIn("PwbTimelagCache(i)%used_lag = PwbTimelagCache(i)%result%maxcov_lag",
                      POST)

    def test_the_arm_still_says_what_it_is(self):
        #> The label was always 'maxcov_default'. It is the value that changed.
        self.assertIn("'maxcov_default'", POST)

    def test_the_old_source_survives_only_as_a_guard(self):
        """fallback_lag is kept for a row carrying no maximum of its own, which
        a table written by this build cannot produce - but a guard that silently
        did the old thing for every row would be the bug reinstated, so it must
        sit under the error test rather than replace it."""
        arm = POST[POST.index("PwbTimelagCache(i)%result%maxcov_lag /= error"):]
        arm = arm[:arm.index("end if")]
        self.assertIn("PwbTimelagCache(i)%result%maxcov_lag", arm)
        self.assertIn("fallback_lag(i)", arm)


class TheCacheCarriesItToTheNextRun(unittest.TestCase):
    """The table is written out and read back, and a lag that survived only in
    memory would be lost the moment a run reused a cache."""

    def test_the_version_moved(self):
        #> The reader accepts exactly one version and refuses the rest, so a
        #> column added without a bump would be read off the end of an old row.
        self.assertEqual(MOD.count("PWB_TIMELAG_CACHE_VERSION=5"), 2)
        self.assertNotIn("PWB_TIMELAG_CACHE_VERSION=4", MOD)

    def test_the_column_is_written(self):
        self.assertIn("maxcov_lag_s", MOD)
        self.assertIn("PwbTimelagCache(i)%result%maxcov_lag", MOD)

    def test_the_column_is_read_back(self):
        self.assertIn("res%maxcov_lag = maxcov_lag", MOD)

    def test_it_is_the_last_column(self):
        """The record is parsed positionally by a list-directed read, so a new
        field anywhere but the end shifts every one after it."""
        header = [ln for ln in read(MODULE).splitlines() if "maxcov_lag_s" in ln][0]
        self.assertTrue(header.rstrip().endswith("maxcov_lag_s'"),
                        "maxcov_lag_s must be the last column: " + header.strip())


if __name__ == "__main__":
    unittest.main()
