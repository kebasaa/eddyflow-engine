"""The summary's columns partition the periods they claim to.

``ReportPwbDiagnostics`` prints one line per gas - attempts, then how many
periods each arm of the classification settled. The columns are meant to add
up to attempts, so every reliability class has to be counted once and by name.

``CountPwbDiagnostic`` decided the last of them with a bare ``else``, and that
is the defect. Any class the ladder above did not claim landed in
``pwb_successes`` and printed under a column headed ``S1/S2`` - a heading that
asserts a lag was measured in that period. ``S4_instrument_filled``, which the
post-pass assigns when a gas takes a settled lag from a neighbour on the same
analyser, reaches exactly there. On ``base_pwb_cache`` that printed
``cos: S1/S2=6`` while the settled table holds no S1 or S2 row for cos at all;
the six were borrowed.

Two things follow, and both are pinned here.

**The successes branch names its classes.** A column headed S1/S2 means S1 or
S2, not "whatever was left over".

**Nothing falls off the end silently.** A class added later would otherwise go
missing from a line whose columns are supposed to sum to attempts. It is
counted and reported instead, so the next one is loud.

``S4_instrument_filled`` is counted with ``S4_instrument_shared``. To a reader
of a summary line they are the same answer - this gas did not detect a lag and
took one measured down the same tube - and they differ only in when the settled
table reached for the neighbour. The S3 arms are already lumped for that
reason, and the column is named ``S4_borrowed`` rather than ``S4_shared`` so
the heading says what it now holds.

Part of the EddyFlow engine's static checks.
"""

from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]

MODULE = "src/src_rp/pwb_timelag_handle.f90"


def read(path):
    return (ROOT / path).read_text(encoding="utf-8")


def code(path):
    return (chr(10)).join(ln for ln in read(path).splitlines()
                          if not ln.lstrip().startswith("!"))


def body_of(source, opener, closer):
    return source[source.index(opener):source.index(closer)]


MOD = code(MODULE)
SRC = read(MODULE)

COUNT = body_of(MOD, "subroutine CountPwbDiagnostic",
                "end subroutine CountPwbDiagnostic")
REPORT = body_of(MOD, "subroutine ReportPwbDiagnostics",
                 "end subroutine ReportPwbDiagnostics")


class ADetectionIsCountedByNameNotByElimination(unittest.TestCase):

    def test_the_successes_branch_names_its_classes(self):
        arm = COUNT[:COUNT.index("pwb_successes(gas) = pwb_successes(gas) + 1")]
        self.assertIn("'S1_optimal'", arm)
        self.assertIn("'S2_optimal'", arm)

    def test_it_is_not_the_catch_all(self):
        """The bare `else` is the bug. Whatever ends up there must not be a
        detection."""
        tail = COUNT[COUNT.rindex("    else"):]
        self.assertNotIn("pwb_successes", tail)


class NothingFallsOffTheEndQuietly(unittest.TestCase):

    def test_there_is_a_counter_for_the_unnamed(self):
        self.assertIn("pwb_unclassified(gas) = pwb_unclassified(gas) + 1", COUNT)

    def test_it_is_reset_with_the_others(self):
        reset = body_of(MOD, "subroutine ResetPwbDiagnostics",
                        "end subroutine ResetPwbDiagnostics")
        self.assertIn("pwb_unclassified = 0", reset)

    def test_it_is_reported(self):
        """A counter nobody prints is the silence this replaces."""
        self.assertIn("pwb_unclassified(gas) > 0", REPORT)


class BorrowedIsBorrowedWhicheverArmReachedForIt(unittest.TestCase):

    def test_both_instrument_arms_count_together(self):
        arm = COUNT[:COUNT.index(
            "pwb_instrument_shared(gas) = pwb_instrument_shared(gas) + 1")]
        self.assertIn("'S4_instrument_shared'", arm)
        self.assertIn("'S4_instrument_filled'", arm)

    def test_the_heading_says_what_it_holds(self):
        self.assertIn("S4_borrowed=", SRC)
        self.assertNotIn("S4_shared=", SRC)

    def test_the_heading_is_written_to_both_streams(self):
        """The console line and the run log are written separately, and the
        log is the copy the regression harness compares."""
        self.assertEqual(SRC.count("', S4_borrowed=', pwb_instrument_shared(gas), &"), 2)


if __name__ == "__main__":
    unittest.main()
