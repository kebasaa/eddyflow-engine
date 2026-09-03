"""The aggregate time-lag dataset is built from the settled table.

Before it computes any flux, a PWB cache-generation run walks every averaging
period and builds a second thing alongside the half-hourly table: the aggregate
dataset, one row per period, which `OptimizeTimelags` reduces to a per-gas
window and a set of H2O relative-humidity classes.

`AddPwbTimelagSummaryDataset` built it as the walk went, and decided each row's
membership from the STREAMING classification - a guess made having read only
the periods before this one. `PostProcessPwbTimelagCache` then overrules that
guess for every row, having read the whole run. So the dataset was assembled
from classifications the table had already discarded, and the two disagree in
practice: on `base_pwb_cache` the streaming walk settles 6 cos periods that the
finished table records as `S4_instrument_filled` instead.

The dataset is now filled afterwards, from the finished table. Only the
humidity cannot come from there - it is a property of the period's own
statistics and is nowhere in the cache - so that is recorded during the walk,
ungated, and the rebuild applies the water gate once the table says whether
water settled.

**What this does not do.** On both PWB fixtures the run output is unchanged,
and that is expected rather than lucky: production overwrites the aggregate
file with its own resolve, and reads its lags from the cache rather than from
the windows this dataset produces. The pre-pass aggregate is a by-product, as
the comment in `optimize_timelags.f90` says. Two things still make this worth
doing - the by-product should agree with the table it describes, and the
streaming version could not have been transported across worker processes,
because a worker starting cold classifies its first periods differently. The
rebuild needs nothing from a worker but the humidity.

The live path is untouched: with no post-pass to defer to, its streaming
classification is the cache's own, since a cache hit fills `PWBResult` from
the table.

Part of the EddyFlow engine's static checks.
"""

from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]

MODULE = "src/src_rp/pwb_timelag_handle.f90"
RP_MAIN = "src/src_rp/eddyflow-rp_main.f90"


def read(path):
    return (ROOT / path).read_text(encoding="utf-8")


def code(path):
    """Source with full-line comments removed."""
    return (chr(10)).join(ln for ln in read(path).splitlines()
                          if not ln.lstrip().startswith("!"))


def body_of(source, opener, closer):
    return source[source.index(opener):source.index(closer)]


MOD = code(MODULE)
MAIN = code(RP_MAIN)

RECORD = body_of(MOD, "subroutine RecordPwbTimelagOptPeriod",
                 "end subroutine RecordPwbTimelagOptPeriod")
REBUILD = body_of(MOD, "subroutine RebuildPwbTimelagOptFromCache",
                  "end subroutine RebuildPwbTimelagOptFromCache")
ADD = body_of(MOD, "subroutine AddPwbTimelagSummaryDataset",
              "end subroutine AddPwbTimelagSummaryDataset")


class TheWalkRecordsOnlyWhatTheTableCannotSay(unittest.TestCase):

    def test_it_does_not_classify(self):
        """The whole point. A reliability class read during the walk is the
        guess the post-pass exists to replace."""
        self.assertNotIn("reliability_class", RECORD)

    def test_it_records_the_humidity(self):
        self.assertIn("Stats%RH", RECORD)

    def test_it_records_which_period_this_is(self):
        """The dataset has no time axis of its own, so without this the
        rebuild cannot find a period's rows in the table."""
        self.assertIn("PwbOptDate(n) = PwbPeriodDate", RECORD)
        self.assertIn("PwbOptTime(n) = PwbPeriodTime", RECORD)

    def test_it_leaves_the_lags_empty(self):
        self.assertIn("TimelagOpt(n)%tlag = error", RECORD)


class TheRebuildReadsTheSettledTable(unittest.TestCase):

    def test_membership_comes_from_the_settled_class(self):
        self.assertIn("PwbTimelagCache(i)%result%reliability_class", REBUILD)
        for arm in ("'S1_optimal'", "'S2_optimal'"):
            self.assertIn(arm, REBUILD)

    def test_the_lag_comes_from_the_settled_row(self):
        self.assertIn("TimelagOpt(k)%tlag(gas) = PwbTimelagCache(i)%used_lag",
                      REBUILD)

    def test_it_clears_before_it_fills(self):
        """Otherwise a row the table no longer settles would keep the lag the
        walk had put there."""
        self.assertLess(REBUILD.index("TimelagOpt(k)%tlag = error"),
                        REBUILD.index("TimelagOpt(k)%tlag(gas) ="))

    def test_the_water_gate_is_applied_here(self):
        """RH travels with the water record's own lag, and only the table
        knows whether water settled - which is precisely what the streaming
        version could not wait for."""
        self.assertIn("PrimaryWaterOutSlot", REBUILD)
        self.assertIn("TimelagOpt(k)%tlag(wsl) = error", REBUILD)
        self.assertIn("TimelagOpt(k)%RH = error", REBUILD)

    def test_the_search_is_a_cursor_not_a_scan(self):
        """Rows and periods are both in period order. A season is some twenty
        thousand rows against four thousand periods, and the nested form of
        that is quadratic for nothing."""
        self.assertIn("PeriodSlot(i, cursor, n)", REBUILD)
        self.assertIn("cursor = k", REBUILD)


class TheOrderInMainIsWhatMakesItCorrect(unittest.TestCase):

    def test_the_walk_records_and_does_not_accumulate(self):
        loop = MAIN.index("call RecordPwbTimelagOptPeriod")
        self.assertGreater(loop, 0)

    def test_the_rebuild_follows_the_post_pass(self):
        self.assertLess(MAIN.index("call PostProcessPwbTimelagCache()"),
                        MAIN.index("call RebuildPwbTimelagOptFromCache"))

    def test_the_rebuild_precedes_what_consumes_the_dataset(self):
        """FixTimelagOptDataset reduces it and OptimizeTimelags fits the
        windows; a rebuild after either would be fitting the empty set."""
        self.assertLess(MAIN.index("call RebuildPwbTimelagOptFromCache"),
                        MAIN.index("call FixTimelagOptDataset(PwbTimelagOpt"))


class TheLivePathIsLeftAlone(unittest.TestCase):
    """It has no post-pass to defer to, and needs none: a cache hit fills
    PWBResult from the table, so its streaming class IS the settled one."""

    def test_it_still_accumulates_as_it_goes(self):
        self.assertIn("call AddPwbTimelagSummaryDataset", MAIN)

    def test_that_routine_still_gates_on_the_class_it_has(self):
        self.assertIn("PWBResult(gas)%reliability_class", ADD)


if __name__ == "__main__":
    unittest.main()
