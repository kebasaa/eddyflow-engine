"""Consecutive-difference despiking, and the two things it must not break.

EddyUH's ``spi_method = 1`` (``EC_Software_Common/EddyUH_despike.m:66-76``): if
a sample steps further from the one before it than a stated limit, replace it
with its predecessor and count a spike. A rate-of-change limit, not a
statistical outlier test - nothing is scaled by a standard deviation and
nothing iterates, which is what separates it from the two methods beside it.

Two properties carry the design:

1. **The historical values of ``despike_vm`` decode exactly as before.** ``'0'``
   was Vickers & Mahrt and *anything else* was Mauder; the new arm is ``'2'``
   and the default arm stays Mauder rather than becoming the nominal default.
   That is what makes an existing project byte-identical.
2. **The three methods publish their outcome the same way.** A spike count that
   reaches the output through a fourth spelling is how one of them would
   quietly stop reporting.
"""

import re
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
NEW = ROOT / "src" / "src_rp" / "despike_consecutive_diff.f90"
VICKERS = ROOT / "src" / "src_rp" / "test_spike_detection_vickers_97.f90"
MAUDER = ROOT / "src" / "src_rp" / "test_spike_detection_mauder_13.f90"
SCREEN = ROOT / "src" / "src_rp" / "statistical_screening.f90"
READER = ROOT / "src" / "src_rp" / "read_ini_rp.f90"
TAGS = ROOT / "src" / "src_rp" / "m_rp_global_var.f90"
GEN = ROOT / "prj" / "gen_project_tags.py"
TYPEDEF = ROOT / "src" / "src_common" / "m_typedef.f90"


def read(path):
    return path.read_text(encoding="utf-8", errors="replace")


NEW_SRC = read(NEW)


class TheHistoricalValuesStillDecodeAsBefore(unittest.TestCase):

    def block(self):
        src = read(READER)
        i = src.index("select case (SCTags(90)%value(1:1))")
        return src[i:src.index("end select", i)]

    def test_zero_is_still_vickers(self):
        self.assertRegex(
            self.block(),
            r"case \('0'\)\s*\r?\n\s*RPSetup%despike_meth = 'vickers_97'")

    def test_the_default_arm_is_mauder_not_the_nominal_default(self):
        #> The old form was `despike_vickers97 = value == '0'`, so every value
        #> that was not '0' meant Mauder - including a malformed one. Making
        #> the default arm Vickers would look tidier and would change what a
        #> malformed project computes.
        self.assertRegex(
            self.block(),
            r"case default\s*\r?\n\s*RPSetup%despike_meth = 'mauder_13'")

    def test_two_is_the_new_arm(self):
        self.assertRegex(
            self.block(),
            r"case \('2'\)(?:\s*\r?\n\s*!>[^\r\n]*)*"
            r"\s*\r?\n\s*RPSetup%despike_meth = 'consecutive_diff'")

    def test_the_dispatch_covers_all_three_and_defaults_to_mauder(self):
        src = read(SCREEN)
        block = src[src.index("select case (trim(adjustl(RPSetup%despike_meth)))"):]
        block = block[:block.index("end select")]
        self.assertIn("call TestSpikeDetectionVickers97", block)
        self.assertIn("call DespikeConsecutiveDiff", block)
        self.assertRegex(
            block, r"case default\s*\r?\n\s*call TestSpikeDetectionMauder13")


class TheThreeMethodsPublishTheSameWay(unittest.TestCase):

    #: Each routine's tail, from the flag string to the end of the where.
    def tail(self, path):
        src = read(path)
        i = src.index("call PackFlagString")
        j = src.index("endwhere", i)
        body = src[i:j]
        #> Local variable names differ by design; compare the shape, not the
        #> spelling of the counters.
        body = re.sub(r"\b(tot_spikes_sng|nreplaced)\b", "REPLACED", body)
        body = re.sub(r"\b(tot_spikes|nspikes)\b", "COUNTED", body)
        #> Only the new routine carries a comment here.
        body = re.sub(r"^\s*!>.*$", "", body, flags=re.M)
        return re.sub(r"\s+", " ", body).strip()

    def test_the_new_routine_matches_vickers(self):
        self.assertEqual(self.tail(NEW), self.tail(VICKERS))

    def test_vickers_and_mauder_still_match_each_other(self):
        #> If this fails the two existing routines have drifted and the
        #> comparison above is anchored to whichever one moved.
        self.assertEqual(self.tail(VICKERS), self.tail(MAUDER))

    def test_the_new_routine_takes_the_hard_flag_per_column(self):
        #> Dividing by N dilutes a slower column's spike count by its stride,
        #> so a 1 Hz analyser on a 10 Hz file needs ten times the spikes to
        #> reach sr%hf_lim. Vickers was fixed for this; the new routine
        #> follows it.
        for path in (NEW, VICKERS):
            self.assertIn("ColumnAcFreq(j) / Metadata%ac_freq", read(path))

    def test_mauders_hard_flag_is_still_per_file_row(self):
        #> NOT a parity assertion - a record of a defect this work found and
        #> deliberately did not fix here. Mauder still divides by dble(N),
        #> which is the stride dilution Vickers had corrected with a written
        #> rationale. Fixing it changes what every project using Mauder
        #> computes, so it belongs in its own commit with its own regression,
        #> not as a rider on a new method.
        #>
        #> When that happens, this test fails and should be replaced by the
        #> parity assertion above extended to all three.
        self.assertIn("dble(tot_spikes(j)) / dble(N)", read(MAUDER))
        self.assertNotIn("ColumnAcFreq(j)", read(MAUDER))


class TheMethodItself(unittest.TestCase):

    def test_it_compares_against_the_last_valid_sample(self):
        #> Not simply i-1. A gap between two good samples would otherwise
        #> manufacture a step out of the values either side of it.
        self.assertIn("prev = 0", NEW_SRC)
        self.assertIn("if (Set(i, j) == error) cycle", NEW_SRC)
        self.assertIn("dabs(Set(i, j) - Set(prev, j)) > step_lim(j)", NEW_SRC)

    def test_a_column_without_a_limit_is_left_alone_and_named(self):
        #> Zero means "not stated", the same shape as EddyUH's dlim where NaN
        #> means the variable is not despiked. A method that silently does
        #> nothing when nobody filled in its numbers is worse than one that
        #> says so.
        self.assertIn("if (step_lim(j) <= 0d0) then", NEW_SRC)
        self.assertIn("No step limit stated for:", NEW_SRC)

    def test_counting_and_replacing_are_separated_by_filter_sr(self):
        #> The count is taken either way; only the replacement is gated, and
        #> the tail then zeroes the replacement count exactly as the other two
        #> routines do. Asserted here because "detected" and "removed" are
        #> different numbers and the output reports both.
        i = NEW_SRC.index("nspikes(j) = nspikes(j) + 1")
        j = NEW_SRC.index("prev = i", i)
        between = NEW_SRC[i:j]
        self.assertIn("if (RPsetup%filter_sr) Set(i, j) = Set(prev, j)",
                      between)
        self.assertIn("if (.not. RPsetup%filter_sr) nreplaced(u:pe) = 0",
                      NEW_SRC)

    def test_nothing_here_is_scaled_by_a_deviation(self):
        #> The moment a sigma multiplier appears this has become a third
        #> variant of the two tests beside it rather than EddyUH's method.
        #> Comments stripped: the header explains what it is NOT by naming
        #> the very things this forbids.
        code = re.sub(r"^\s*!.*$", "", NEW_SRC, flags=re.M)
        for token in ("sr%lim_u", "sr%lim_w", "sr%lim_gas", "StDev", "MAD"):
            self.assertNotIn(token, code)

    def test_the_sonic_components_have_limits_of_their_own(self):
        #> These are absolute limits in each variable's own units. Sharing one
        #> key across u and Ts, as the Vickers sigma multipliers do, would put
        #> a kelvin and a metre per second behind the same number.
        for f in ("sr%step_u", "sr%step_v", "sr%step_w", "sr%step_ts"):
            self.assertIn(f, NEW_SRC)


class TheKeysAreWhereTheyWereSaidToBe(unittest.TestCase):

    def test_the_four_sonic_limits_hold_their_slots(self):
        tags = read(TAGS)
        for slot, name in ((61, "sr_step_u"), (62, "sr_step_v"),
                           (63, "sr_step_w"), (64, "sr_step_ts")):
            self.assertRegex(
                tags, r"SNTags\(%d\)%%Label\s*/\s*'%s'\s*/" % (slot, name))

    def test_the_per_gas_limit_is_appended_not_inserted(self):
        #> RP_GAS_NUMERIC is positional: inserting anywhere but the end
        #> silently re-points every later per-gas setting of every gas.
        gen = read(GEN)
        block = gen[gen.index("RP_GAS_NUMERIC = ["):gen.index("RP_GAS_TEXT")]
        self.assertRegex(block, r'\+ \["step_lim"\]\s*$',
                         "step_lim is no longer the last entry")
        self.assertLess(block.index("sr_lim"), block.index("step_lim"))

    def test_the_reader_takes_the_offset_from_the_generated_stride(self):
        #> Not a hand-counted index: rpGasLeapN grows with the list, and a
        #> literal here would go stale the next time a per-gas key is added.
        self.assertIn("rpGasOriginN + (i - 1) * rpGasLeapN + 25", read(READER))

    def test_every_step_limit_read_is_guarded(self):
        src = read(READER)
        block = src[src.index("sr%step_u = 0d0"):src.index("!> Dropout test")]
        for slot in (61, 62, 63, 64):
            self.assertIn("if (SNTagFound(%d))" % slot, block)
        self.assertIn("if (SNTagFound(rpGasOriginN", block)

    def test_the_selector_is_no_longer_a_logical(self):
        #> A boolean named for one method cannot describe a third.
        self.assertNotIn("despike_vickers97", read(TYPEDEF))
        self.assertIn("character(32) :: despike_meth", read(TYPEDEF))


if __name__ == "__main__":
    unittest.main()
