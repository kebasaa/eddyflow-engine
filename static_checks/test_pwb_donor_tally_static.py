"""The aggregate summary counts donors off the settled table, not the guesses.

`ResolvePwbAggregateSummary` picks each gas's lender by count: whichever
same-analyser, non-water candidate lent to it most often wins. So that count
decides a real number - the time-lag window a gas inherits, which
`SetTimelags` writes into `E2Col(gas)%min_tl/max_tl`.

It used to be accumulated by `AddPwbTimelagSummaryDataset` as the pre-pass
walked, off the STREAMING classification in timelag_handle. Two things follow
from that, and both are why the count now comes from the finished table.

**It was simply wrong.** The streaming classifier decides a period having read
only the periods before it; `PostProcessPwbTimelagCache` then re-decides every
row having read the whole run, and overwrites `reliability_class` on all of
them. The tally was never revisited, so it described guesses the table had
since overruled. Measured on `base_pwb_cache`, three hours of CH-LAE: the
streaming pass counted **2**, the settled table **7**, for the same gas pair.
That is the same reason the run-log tallies are recounted from the table three
lines above - this array was left behind when they were fixed.

**It leaked pass order into the result.** The streaming chain depends on where
the pass began, so a worker process starting cold classifies its first periods
differently and tallies a different lender. That is what stops the PWB
pre-pass being split, and this is one of the values that has to stop carrying
it. The others have their own checks.

A note on what this does NOT cover: production accumulates into the same array
during the run proper, and resolves again at the end. That path is already
table-derived - a cache hit fills `PWBResult` from the cache - so it is left
alone.

Part of the EddyFlow engine's static checks.
"""

from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]

MODULE = "src/src_rp/pwb_timelag_handle.f90"


def read(path):
    return (ROOT / path).read_text(encoding="utf-8")


def code(path):
    """Source with full-line comments removed."""
    return (chr(10)).join(ln for ln in read(path).splitlines()
                          if not ln.lstrip().startswith("!"))


def body_of(source, opener, closer):
    return source[source.index(opener):source.index(closer)]


SRC = read(MODULE)
CODE = code(MODULE)

POST = body_of(CODE, "subroutine PostProcessPwbTimelagCache",
               "end subroutine PostProcessPwbTimelagCache")
ADD = body_of(CODE, "subroutine AddPwbTimelagSummaryDataset",
              "end subroutine AddPwbTimelagSummaryDataset")


class TheTallyIsRebuiltFromTheSettledTable(unittest.TestCase):

    def test_it_is_cleared_before_being_counted(self):
        """Or the streaming accumulation would be added to, not replaced."""
        self.assertIn("PwbSummaryDonorCount = 0", POST)

    def test_it_is_counted_from_the_cache_rows(self):
        self.assertIn("origin = PwbTimelagCache(i)%result%origin_gas", POST)
        self.assertIn("PwbSummaryDonorCount(gas, origin)", POST)

    def test_the_clear_comes_before_the_count(self):
        self.assertLess(POST.index("PwbSummaryDonorCount = 0"),
                        POST.index("origin = PwbTimelagCache(i)%result%origin_gas"))

    def test_it_happens_inside_the_post_pass(self):
        """Not in the caller. ResolvePwbAggregateSummary runs immediately
        after PostProcessPwbTimelagCache in the cache-generation branch, and
        anything between the two would have to be kept in step by hand."""
        self.assertIn("PwbSummaryDonorCount = 0", POST)


class ItCountsTheSameArmTheStreamingPassDid(unittest.TestCase):
    """A different rule would be a behaviour change dressed as a refactor."""

    def test_a_settled_row_is_not_a_borrowing(self):
        for arm in ("'S1_optimal'", "'S2_optimal'"):
            self.assertIn(arm, POST)

    def test_a_gas_does_not_lend_to_itself(self):
        #> origin_gas is the gas itself on every native row, so without this
        #> every settled period would read as a loan from nobody.
        self.assertIn("origin /= gas", POST)
        self.assertIn("origin /= gas", ADD)

    def test_the_donor_slot_is_range_checked(self):
        self.assertIn("origin >= firstGas .and. origin <= lastGas", POST)
        self.assertIn("origin >= firstGas .and. origin <= lastGas", ADD)

    def test_an_absent_gas_is_skipped(self):
        self.assertIn("if (.not. E2Col(gas)%present) cycle", POST)


class TheReasonIsWrittenDown(unittest.TestCase):

    def test_the_comment_names_the_routine_that_reads_the_tally(self):
        """A bare `PwbSummaryDonorCount = 0` in a 400-line routine reads like
        a reset that could be moved or dropped. It cannot: the resolve runs
        straight after and consumes it."""
        self.assertIn("ResolvePwbAggregateSummary", SRC)

    def test_it_says_the_streaming_tally_was_what_it_replaced(self):
        self.assertIn("AddPwbTimelagSummaryDataset", SRC)


if __name__ == "__main__":
    unittest.main()
