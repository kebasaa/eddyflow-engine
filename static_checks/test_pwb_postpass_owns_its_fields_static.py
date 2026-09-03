"""Every field saying HOW a period was settled is written by the post-pass.

`PostProcessPwbTimelagCache` re-decides every row of the table once the whole
run has been read, and overwrites `reliability_class` on all of them. But the
row carries more than a class: `donor_gas`, `carry_hours`, `fallback_source`
and `fallback_used` all describe how the lag was arrived at, and they are what
the half-hourly file reports.

Not one arm sets all of them. Interpolate, carry-forward, back-fill and median
leave `donor_gas` alone; median leaves `carry_hours` too. A field an arm does
not set keeps what the STREAMING pass left in the row - a value the table has
just finished overruling, describing a decision that no longer stands.

**Why a serial run cannot show this.** The leftover is at least the same
leftover every time, so the file is wrong but stable, and a diff against a
previous run of the same build says nothing. It took running the pre-pass
across worker processes to make it visible: a worker starting cold leaves
different leftovers, and `base_pwb_par` disagreed with its own serial run on
three `S3_interpolated` rows, whose `donor_gas` read `h2o` and `co2` one way
and `co2` and `none` the other - with the settled lag, the class and
`origin_gas` identical in both. Neither answer was right. They were merely
differently wrong.

So the fields are blanked for every row before the tiers run, and each arm
overrides what it decides. Blanking rather than patching the four arms is the
point: an arm added later inherits a defined value instead of a stale one,
which is the mistake this file exists to stop repeating.

Safe because nothing reads them in between - every use of all four inside the
routine is a write.

Part of the EddyFlow engine's static checks.
"""

from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]

MODULE = "src/src_rp/pwb_timelag_handle.f90"

FIELDS = ("donor_gas", "carry_hours", "fallback_source", "fallback_used")


def read(path):
    return (ROOT / path).read_text(encoding="utf-8")


def code(path):
    return (chr(10)).join(ln for ln in read(path).splitlines()
                          if not ln.lstrip().startswith("!"))


def body_of(source, opener, closer):
    return source[source.index(opener):source.index(closer)]


MOD = code(MODULE)
POST = body_of(MOD, "subroutine PostProcessPwbTimelagCache",
               "end subroutine PostProcessPwbTimelagCache")

BLANK = "PwbTimelagCache(i)%result%donor_gas = 'none'"


class TheFieldsStartBlank(unittest.TestCase):

    def test_all_four_are_cleared(self):
        for field in FIELDS:
            self.assertIn("PwbTimelagCache(i)%result%" + field + " =",
                          POST[:POST.index("do gas = firstGas, lastGas")],
                          field + " is not blanked before the tiers run")

    def test_it_happens_before_any_arm_decides(self):
        """After the first tier it would erase what that tier had just set."""
        self.assertLess(POST.index(BLANK),
                        POST.index("do gas = firstGas, lastGas"))

    def test_it_happens_after_the_prefilter(self):
        """The pre-filter is step 1 and sets hdi_prefiltered, which is not one
        of these - the ordering only has to put the blanking ahead of the arms."""
        self.assertLess(POST.index("hdi_prefiltered = .false."),
                        POST.index(BLANK))


class TheArmsStillSayWhatTheyDecide(unittest.TestCase):
    """Blanking is a floor, not a replacement. An arm that borrows a lag has to
    name the lender, or the file says 'none' for a borrowed period."""

    def test_a_settled_period_has_no_donor(self):
        self.assertIn("PwbTimelagCache(i)%result%donor_gas = 'none'", POST)

    def test_a_borrowed_period_names_its_lender(self):
        self.assertIn(
            "PwbTimelagCache(i)%result%donor_gas = GasLabel(PwbTimelagCache(shared)%gas)",
            POST)

    def test_the_terminal_arm_still_raises_the_fallback_flag(self):
        """It is the one arm where the blanked default is the wrong answer."""
        self.assertIn("PwbTimelagCache(i)%result%fallback_used = .true.", POST)


class NothingReadsThemInBetween(unittest.TestCase):
    """Which is what makes blanking safe rather than destructive."""

    def test_every_use_inside_the_routine_is_a_write(self):
        for field in FIELDS:
            token = "%result%" + field
            for line in POST.split(chr(10)):
                if token not in line:
                    continue
                after = line.split(token, 1)[1].lstrip()
                self.assertTrue(
                    after.startswith("="),
                    "%s is read, not written, at: %s" % (field, line.strip()))


if __name__ == "__main__":
    unittest.main()
