"""The RH-class counts reach the report instead of being thrown away.

``WriteOutTimelagOptimization`` prints, per relative-humidity class, the median
lag and the number of determinations behind it. The count is computed by
``OptimizeTimelags``, which fills the caller's array; the caller then hands the
same array to the writer.

The writer used to declare it ``intent(out)``.

**An ``intent(out)`` dummy is undefined on entry.** So the declaration threw
away the counts the caller had just computed, and the ``class_num`` column
printed whatever happened to be on the stack - values like ``1818717765``,
which is ``0x6C696D45``: text, read as an integer. Every other column of the
table was correct, which is what made it survive since ``c74a11b``, the
initial EddyPro 6.2.2 fork.

Nothing about the fluxes changed - this is a diagnostic report - but a count
column that is silently garbage is worse than no column, because a reader
compares it against the threshold printed directly above it and concludes that
the classes were rejected for a reason.

Two properties are pinned here, and a third that the fix exposed:

1. **The writer only reads the counts**, so they arrive intact.
2. **The optimiser writes them**, so there is something to arrive.
3. **The threshold in the sentence is the threshold in the gate.** The report
   claimed "numerosity < 30" while the optimiser used 15, and had done since
   the fork. Nobody could see it while the column beside it was noise. Both now
   come from one parameter.
"""

import re
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
WRITER = ROOT / "src" / "src_rp" / "writeout_timelag_optimization.f90"
OPTIMISER = ROOT / "src" / "src_rp" / "optimize_timelags.f90"
TYPEDEF = ROOT / "src" / "src_common" / "m_typedef.f90"
RP_MAIN = ROOT / "src" / "src_rp" / "eddyflow-rp_main.f90"


def read(p):
    return p.read_text(encoding="utf-8", errors="replace")


W = read(WRITER)
O = read(OPTIMISER)
T = read(TYPEDEF)
M = read(RP_MAIN)


class TheWriterOnlyReadsTheCounts(unittest.TestCase):
    """The bug itself, in one declaration."""

    def test_the_dummy_is_intent_in(self):
        self.assertRegex(W, r"integer, intent\(in\) :: h2o_n\(ncls\)")

    def test_it_is_not_intent_out(self):
        #> The whole defect. intent(out) discards what the caller passed.
        self.assertNotRegex(W, r"intent\(out\) *:: *h2o_n")

    def test_the_writer_never_assigns_it(self):
        #> If it ever did, intent(in) would fail to compile - but a future
        #> edit might "fix" that by widening the intent instead, which is
        #> how this started.
        for m in re.finditer(r"^\s*h2o_n\s*\(?[^=!]*\)?\s*=[^=]", W, re.M):
            self.fail("the writer assigns h2o_n at %r" % m.group(0).strip())

    def test_it_is_still_printed(self):
        #> A guard that dropped the column would also pass the tests above.
        self.assertIn("h2o_n(cls)", W)


class TheOptimiserProducesThem(unittest.TestCase):

    def test_the_optimiser_declares_them_out(self):
        self.assertRegex(O, r"integer, intent\(out\) :: h2o_n\(MM\)")

    def test_the_optimiser_counts_into_them(self):
        self.assertIn("h2o_n(cls) = 0", O)
        self.assertIn("h2o_n(cls) = h2o_n(cls) + 1", O)

    def test_the_caller_fills_before_it_writes(self):
        #> Every WriteOutTimelagOptimization call must be preceded by an
        #> OptimizeTimelags call on the same array - otherwise the counts are
        #> undefined again, by a different route.
        fills = [m.start() for m in
                 re.finditer(r"call OptimizeTimelags\(", M)]
        writes = [m.start() for m in
                  re.finditer(r"call WriteOutTimelagOptimization\(", M)]
        self.assertTrue(writes, "no call sites found")
        for w in writes:
            self.assertTrue(any(f < w for f in fills),
                            "a write at %d has no fill before it" % w)

    def test_the_array_the_caller_passes_is_the_one_it_filled(self):
        self.assertIn("call OptimizeTimelags(toSet, size(toSet), tlagn, "
                      "E2NumVar, toH2On,", M)
        self.assertIn("toH2On", M)


class TheStatedThresholdIsTheRealOne(unittest.TestCase):
    """Exposed by the fix: the sentence said 30, the gate used 15."""

    def test_the_threshold_is_one_shared_parameter(self):
        self.assertRegex(T, r"integer, parameter :: toMinH2OClassN = \d+")

    def test_the_optimiser_uses_it(self):
        self.assertIn("min_numerosity = toMinH2OClassN", O)

    def test_the_report_prints_it_rather_than_a_literal(self):
        #> Spelling the number out in the sentence is exactly how the two
        #> drifted apart in the first place.
        self.assertIn("toMinH2OClassN", W)
        self.assertNotIn("numerosity < 30", W)
        self.assertNotRegex(W, r"numerosity < \d")

    def test_the_gate_is_the_one_the_sentence_describes(self):
        #> "inferred" means the class got no direct fit, which is the
        #> N < min_numerosity branch - not the > used to find the first and
        #> last good class for the edges.
        self.assertIn("if (N < min_numerosity) cycle", O)


if __name__ == "__main__":
    unittest.main()
