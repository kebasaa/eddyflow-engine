"""Static checks for detection-limit-conditional time-lag borrowing.

Nemitz et al. (2018): a gas whose cross-covariance peak cannot be told from
noise takes the lag of a gas that shares its tube, because gases drawn down
one tube share a transport delay and the weak one has nothing of its own to
detect.

**The call site and the signature must agree.** This is an external
subroutine with no interface block, so a mismatched argument list compiles
without complaint. During development the call still passed ``ActTLag`` after
it had been dropped from the signature, which shifted every following
argument by one: ``DefTlagUsed(gas) = .true.`` then wrote a logical into the
real ``TLag`` array, and one gas's published time lag came out as the
denormal 0.2122E-313 - in periods where nothing had been borrowed for it.
Nothing raised. The argument list is pinned here because the compiler will
not do it.

**It must not run without a detection limit.** There is nothing to compare a
covariance against, so the routine returns early and the interface greys the
control.

**A borrowed lag must not be borrowed again.** Otherwise one detection walks
down a whole tube. The trusted set is snapshotted before anything is taken.

**The donor is the best-resolved tube-mate, not the first one.** Taking the
first put whichever gas sat lowest in the metadata in charge, which on a tube
where nothing clears the threshold comfortably produced pairings that read
backwards.
"""

import re
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / "src"


def read(rel):
    return (SRC / rel).read_text(encoding="utf-8", errors="replace")


BORROW = read("src_rp/borrow_timelag.f90")
TLAG = read("src_rp/timelag_handle.f90")
READER = read("src_rp/read_ini_rp.f90")
TAGS = read("src_rp/m_rp_global_var.f90")
TYPEDEF = read("src_common/m_typedef.f90")

ROUTINE = "BorrowTimelagBelowDetectionLimit"


class TheCallSiteMatchesTheSignature(unittest.TestCase):
    """No interface block, so nothing else checks this."""

    @staticmethod
    def arglist(text, anchor):
        i = text.index(anchor)
        chunk = text[i:i + 400]
        chunk = chunk[chunk.index("(") + 1: chunk.index(")")]
        chunk = chunk.replace("&", " ").replace("\n", " ")
        return [a.strip() for a in chunk.split(",") if a.strip()]

    def test_the_names_and_order_agree(self):
        call = self.arglist(TLAG, "call " + ROUTINE + "(")
        sig = self.arglist(BORROW, "subroutine " + ROUTINE + "(")
        self.assertEqual(
            call, sig,
            "the call passes %s and the subroutine declares %s; with no "
            "interface block this compiles and silently writes each argument "
            "into the next one's storage" % (call, sig))

    def test_the_actual_lag_is_not_among_them(self):
        #> ActTLag reports what each gas's own maximisation found. Borrowing
        #> must leave it alone, so the two differing is the record that a lag
        #> was taken from elsewhere.
        sig = self.arglist(BORROW, "subroutine " + ROUTINE + "(")
        self.assertNotIn("ActTLag", sig)


class TheSettingsReachTheEngine(unittest.TestCase):

    def test_the_keys_hold_their_slots(self):
        self.assertRegex(TAGS, r"SCTags\(80\)%Label\s*/\s*'tlag_borrow_meth'\s*/")
        self.assertRegex(TAGS, r"SNTags\(60\)%Label\s*/\s*'tlag_borrow_snr'\s*/")

    def test_both_are_read_under_found_guards(self):
        self.assertIn("SCTagFound(80)", READER)
        self.assertIn("SNTagFound(60)", READER)

    def test_the_threshold_defaults_to_the_papers_three(self):
        self.assertRegex(READER, r"RPSetup%tlag_borrow_snr = 3d0")

    def test_a_non_positive_threshold_is_refused(self):
        #> Zero would make every gas borrow, which nobody means.
        self.assertRegex(
            READER, r"if \(RPSetup%tlag_borrow_snr <= 0d0\)")

    def test_the_setup_type_carries_both(self):
        self.assertRegex(TYPEDEF, r"logical :: tlag_borrow_meth")
        self.assertRegex(TYPEDEF, r"real\(kind = dbl\) :: tlag_borrow_snr")


class ItIsInertWithoutADetectionLimit(unittest.TestCase):

    def test_it_returns_when_not_asked_for(self):
        self.assertIn("if (.not. RPSetup%tlag_borrow_meth) return", BORROW)

    def test_it_returns_when_there_is_nothing_to_compare_against(self):
        self.assertIn(
            "if (RPSetup%detlim_meth /= 'wienhold_94') return", BORROW,
            "the routine would divide by a detection limit that was never "
            "computed")

    def test_both_guards_precede_any_work(self):
        body = BORROW[BORROW.index("subroutine " + ROUTINE):]
        first_work = body.index("ColW(1:nrow) = Set")
        self.assertLess(body.index("if (.not. RPSetup%tlag_borrow_meth) return"),
                        first_work)
        self.assertLess(
            body.index("if (RPSetup%detlim_meth /= 'wienhold_94') return"),
            first_work)


class TheBorrowingRuleHolds(unittest.TestCase):

    def test_the_trusted_set_is_decided_before_anything_is_taken(self):
        #> Otherwise a borrowed lag becomes a donor and one detection walks
        #> down the whole tube.
        settle = BORROW.index("trusted(gas) = snr(gas) >=")
        take = BORROW.index("RowLags(gas) = RowLags(chosen)")
        self.assertLess(settle, take)

    def test_a_donor_must_share_the_analyser(self):
        self.assertIn("if (.not. SameAnalyser(gas, donor)) cycle", BORROW)

    def test_the_donor_is_the_best_resolved_not_the_first(self):
        self.assertIn("if (snr(donor) > best) then", BORROW)
        self.assertNotRegex(
            BORROW,
            r"chosen = donor\s*\n\s*exit",
            "the donor is being taken in slot order again, which puts "
            "whichever gas sits lowest in the metadata in charge")

    def test_water_is_excluded_on_both_sides(self):
        #> Its lag is the one every other gas's water covariance is taken at.
        self.assertIn("if (GasSlotIsWater(g)) return", BORROW)

    def test_open_paths_are_excluded(self):
        self.assertIn("path_type /= 'closed') return", BORROW)

    def test_a_borrowed_lag_is_flagged_as_not_its_own(self):
        self.assertIn("DefTlagUsed(gas) = .true.", BORROW)

    def test_the_donor_is_named_in_the_log(self):
        self.assertIn("' taken from '", BORROW)


class TheCallSiteSeesWhatItNeeds(unittest.TestCase):

    def test_it_runs_after_the_detection_limit(self):
        detlim = TLAG.index("call FluxDetectionLimit(")
        borrow = TLAG.index("call " + ROUTINE + "(")
        self.assertLess(
            detlim, borrow,
            "borrowing reads a detection limit that has not been computed yet")

    def test_it_runs_before_the_series_are_shifted(self):
        borrow = TLAG.index("call " + ROUTINE + "(")
        shift = TLAG.index("Align data according to relevant time-lags")
        self.assertLess(
            borrow, shift,
            "a lag changed after the shift would not take effect")


class TheBigArraysAreOnTheHeap(unittest.TestCase):
    """Two 288 kB automatic arrays, from deep inside the period loop."""

    def test_the_columns_are_allocatable(self):
        self.assertRegex(BORROW, r"real\(kind = dbl\), allocatable :: ColW\(:\)")
        self.assertRegex(BORROW, r"real\(kind = dbl\), allocatable :: ColGas\(:\)")

    def test_they_are_released(self):
        self.assertIn("deallocate(ColW)", BORROW)
        self.assertIn("deallocate(ColGas)", BORROW)


if __name__ == "__main__":
    unittest.main()
