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
        #> Only the detection limit needs a second switch. It is computed
        #> elsewhere and read from Essentials, so asking for it without
        #> switching it on would divide by a number nothing produced. The
        #> Lenschow noise is measured here from the series in hand, so the
        #> guard is conditional on which floor was chosen - which is what
        #> lets EddyUH's combination be selected on its own.
        self.assertIn(
            "if (RPSetup%tlag_borrow_noise == 'detlim' .and. &", BORROW,
            "the routine would divide by a detection limit that was never "
            "computed")
        self.assertIn("RPSetup%detlim_meth /= 'wienhold_94') return", BORROW)

    def test_both_guards_precede_any_work(self):
        body = BORROW[BORROW.index("subroutine " + ROUTINE):]
        first_work = body.index("ColW(1:nrow) = Set")
        self.assertLess(body.index("if (.not. RPSetup%tlag_borrow_meth) return"),
                        first_work)
        self.assertLess(
            body.index("if (RPSetup%tlag_borrow_noise == 'detlim' .and. &"),
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


class TheTwoChoicesDefaultToOurOwn(unittest.TestCase):
    """Which noise floor, and which donor.

    Both were ported wrong the first time. EddyUH tests a covariance against
    the Lenschow instrument noise and borrows specifically from carbon
    dioxide; the first pass tested against the Wienhold detection limit and
    borrowed from the best-resolved tube-mate. EddyUH's own inline comment at
    EddyUH_SC_Flux2.m:325 calls unc3 a detection limit, which is how that
    happened.

    Both rules are now selectable, and both DEFAULT to this engine's own,
    so a project that already had borrowing switched on keeps what it had.
    """

    def reader(self):
        #> READER is already the file's text, not a path.
        i = READER.index("RPSetup%tlag_borrow_noise = 'detlim'")
        return READER[i:i + 700]

    def test_the_defaults_are_ours_not_eddyuhs(self):
        block = self.reader()
        self.assertIn("RPSetup%tlag_borrow_noise = 'detlim'", block)
        self.assertIn("RPSetup%tlag_borrow_donor = 'best_resolved'", block)

    def test_only_an_explicit_one_selects_eddyuhs(self):
        #> Guarded reads, and an equality against '1' rather than a negation,
        #> so an absent tag and an unknown value both leave the default.
        block = self.reader()
        self.assertIn("if (SCTagFound(81) .and. SCTags(81)%value(1:1) == '1')",
                      block)
        self.assertIn("if (SCTagFound(82) .and. SCTags(82)%value(1:1) == '1')",
                      block)
        self.assertIn("RPSetup%tlag_borrow_noise = 'lenschow_00'", block)
        self.assertIn("RPSetup%tlag_borrow_donor = 'carbon_dioxide'", block)

    def test_the_noise_floor_is_chosen_once_where_it_is_divided_by(self):
        #> Inside SignalToNoise, which is the only place either floor is
        #> used. Choosing at the call sites instead would need the same
        #> branch in two places.
        snr = BORROW[BORROW.index("function SignalToNoise"):]
        snr = snr[:snr.index("end function SignalToNoise")]
        self.assertIn("if (RPSetup%tlag_borrow_noise == 'lenschow_00') then",
                      snr)
        self.assertIn("floor_ = LenschowFluxNoise(Set, nrow, ncol, g)", snr)
        self.assertIn("floor_ = Essentials%detlim(g)", snr)
        self.assertIn("SignalToNoise = dabs(cov) / floor_", snr)

    def test_either_floor_is_refused_when_absent_or_non_positive(self):
        snr = BORROW[BORROW.index("function SignalToNoise"):]
        snr = snr[:snr.index("end function SignalToNoise")]
        self.assertIn("if (floor_ == error) return", snr)
        self.assertIn("if (floor_ <= 0d0) return", snr)

    def test_the_carbon_donor_refuses_rather_than_falling_back(self):
        #> CarbonOnAnalyserOf returns 0 when that analyser measures none, and
        #> nothing here substitutes another instrument's carbon dioxide. The
        #> point of the rule is a shared tube; another analyser's gas shares
        #> only a clock.
        block = BORROW[BORROW.index("if (RPSetup%tlag_borrow_donor =="):]
        block = block[:block.index("else")]
        self.assertIn("donor = CarbonOnAnalyserOf(gas)", block)
        self.assertIn("if (donor /= gas .and. trusted(donor)) chosen = donor",
                      block)
        self.assertNotIn("PrimaryCarbon", block)

    def test_the_carbon_donor_still_has_to_be_trusted(self):
        #> EddyUH takes the carbon dioxide lag unconditionally. Requiring it
        #> to have cleared the threshold itself is a deliberate departure: a
        #> donor that could not resolve its own peak is not evidence, and
        #> copying it would launder one bad detection onto every gas on the
        #> tube.
        block = BORROW[BORROW.index("if (RPSetup%tlag_borrow_donor =="):]
        block = block[:block.index("else")]
        self.assertIn("trusted(donor)", block)

    def test_the_best_resolved_rule_is_untouched(self):
        #> The default arm. It still ranks by signal-to-noise over the
        #> analyser's own gases rather than taking the first in slot order.
        block = BORROW[BORROW.index("best = 0d0"):]
        block = block[:block.index("end do") + 6]
        self.assertIn("if (.not. SameAnalyser(gas, donor)) cycle", block)
        self.assertIn("if (snr(donor) > best) then", block)


if __name__ == "__main__":
    unittest.main()
