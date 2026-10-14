"""The preliminary time-lag assessment walks a random sample of the periods
and writes every lag it finds.

The interface's "Timelag preliminary assessment" dialog runs a throw-away copy
of the project as an assessment-only time-lag pre-pass over a wide window. Four
keys, which only that copy ever carries, ask for it:

    to_sample_frac   share of the periods with raw data to walk (0 = all)
    to_sample_min_n  at least this many of them
    to_sample_seed   the draw's seed
    to_dump_lags     write <id>_tlag_preassess_<ts>.csv

The sample has to be the same in every process of a split pre-pass, or the
workers walk periods the parent never meant to; and the periods outside it
must cost nothing, or the share buys no time.
"""

import re
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def read(rel):
    return (ROOT / rel).read_text(encoding="utf-8", errors="replace").replace("\r\n", "\n")


GLOBALS = read("src/src_rp/m_rp_global_var.f90")
READER = read("src/src_rp/read_ini_rp.f90")
MAIN = read("src/src_rp/eddyflow-rp_main.f90")
POOL = read("src/src_common/batch_pool.f90")
MOD = read("src/src_rp/tlag_preassess.f90")
ADD = read("src/src_rp/add_to_timelag_opt_dataset.f90")
TYPEDEF = read("src/src_common/m_typedef.f90")


def body(text, kind, name):
    i = text.index("%s %s(" % (kind, name))
    return text[i:text.index("end %s %s" % (kind, name), i)]


class TheKeysSitInFreeSlots(unittest.TestCase):

    def test_each_key_has_its_own_slot_and_is_read_from_it(self):
        for idx, label in ((365, "to_sample_frac"), (366, "to_sample_min_n"),
                           (367, "to_sample_seed"), (368, "to_dump_lags")):
            self.assertRegex(GLOBALS, r"SNTags\(%d\)%%Label\s*/\s*'%s'" % (idx, label))
            self.assertIn("SNTagFound(%d)" % idx, READER)

    def test_they_mean_nothing_outside_an_on_the_fly_pre_pass(self):
        i = READER.index("TOSetup%sample_frac = 0d0")
        block = READER[i:READER.index("end if", READER.index("SNTagFound(368)", i))]
        self.assertIn("if (RPsetup%to_onthefly) then", block)
        self.assertIn("TOSetup%dump_lags = .false.", block)

    def test_absent_is_off(self):
        self.assertIn("real(kind = dbl) :: sample_frac = 0d0", TYPEDEF)
        self.assertIn("logical :: dump_lags = .false.", TYPEDEF)


class EveryProcessDrawsTheSameSample(unittest.TestCase):

    def test_drawn_over_the_whole_range_before_a_worker_narrows_it(self):
        draw = MAIN.index("call DrawTlagSample(toStartTimestampIndx, toEndTimestampIndx")
        narrow = MAIN.index("toStartTimestampIndx = BatchSliceStart")
        self.assertLess(draw, narrow)
        self.assertLess(draw, MAIN.index("call PlanPrepassBatches(toWalk"))

    def test_the_draw_is_a_fixed_sequence_not_the_compilers(self):
        draw = body(MOD, "subroutine", "DrawTlagSample")
        self.assertNotIn("random_number", draw)
        self.assertNotIn("random_seed", draw)
        self.assertIn("state = mod(PmMultiplier * state, PmModulus)", draw)

    def test_only_periods_with_a_file_are_drawn_and_at_least_the_minimum(self):
        draw = body(MOD, "subroutine", "DrawTlagSample")
        self.assertIn("call PeriodsWithFile(", draw)
        self.assertIn("n = max(TOSetup%sample_min_n, nint(TOSetup%sample_frac * dble(nAvail)))", draw)
        self.assertIn("n = min(n, nAvail)", draw)

    def test_nothing_is_drawn_when_nothing_was_asked(self):
        draw = body(MOD, "subroutine", "DrawTlagSample")
        self.assertLess(draw.index("if (TOSetup%sample_frac <= 0d0"),
                        draw.index("call SetPeriodSample("))


class AnUndrawnPeriodCostsNothing(unittest.TestCase):

    def test_it_is_passed_before_its_file_is_looked_for(self):
        i = MAIN.index("to_periods_loop: do")
        loop = MAIN[i:MAIN.index("end do to_periods_loop", i)]
        self.assertLess(loop.index("if (.not. PeriodSampled(pcount)) cycle to_periods_loop"),
                        loop.index("call FirstFileOfCurrentPeriod("))
        # After the exit test, so the loop still ends where it did.
        self.assertLess(loop.index("exit to_periods_loop"),
                        loop.index("PeriodSampled(pcount)"))

    def test_the_sample_ends_with_the_pre_pass(self):
        i = MAIN.index("to_periods_loop: do")
        end = MAIN.index("call ProgressEnd()", i)
        self.assertLess(end, MAIN.index("call ClearPeriodSample()", i))
        self.assertLess(MAIN.index("call ClearPeriodSample()", i),
                        MAIN.index("FinishBatchWorker()", i))

    def test_progress_counts_walked_periods(self):
        say = body(POOL, "subroutine", "ProgressSay")
        self.assertIn("SpanWork(ProgFirst, ProgEnd)", say)
        span = body(POOL, "function", "SpanWork")
        self.assertIn("if (.not. SampleOn) then", span)
        self.assertIn("SpanWork = max(0, b - a)", span)


class EveryLagFoundIsTabled(unittest.TestCase):

    def test_the_record_carries_the_raw_lag_and_its_verdicts(self):
        i = TYPEDEF.index("type :: TimeLagOptType")
        t = TYPEDEF[i:TYPEDEF.index("end type TimeLagOptType", i)]
        for field in ("RawLag(E2NumVar)", "WinMin(E2NumVar)", "WinMax(E2NumVar)",
                      "Flux(E2NumVar)", "AtEdge(E2NumVar)", "FluxOk(E2NumVar)",
                      "character(10) :: date", "character(5) :: time"):
            self.assertIn(field, t)

    def test_both_arms_note_the_lag_after_their_own_gate(self):
        self.assertEqual(ADD.count("call NoteRawLag("), 2)
        note = body(ADD, "subroutine", "NoteRawLag")
        self.assertIn("lag = Essentials%actual_timelag(gas)", note)

    def test_the_gates_themselves_are_untouched(self):
        self.assertEqual(ADD.count(".and. Essentials%used_timelag(gas) /= E2Col(gas)%max_tl"), 2)

    def test_written_before_the_dataset_is_compacted(self):
        i = MAIN.index("call WriteTlagPreassessTable(TimelagOpt, TimelagOptSize, ton)")
        self.assertLess(i, MAIN.index("call FixTimelagOptDataset(TimelagOpt, TimelagOptSize", i - 2000))
        self.assertLess(MAIN.index("call MergeTlagBatchDumps('to'"), i)

    def test_pwb_tables_its_own_detection_after_the_table_is_settled(self):
        i = MAIN.index("call WritePwbPreassessTable()")
        self.assertLess(MAIN.index("call PostProcessPwbTimelagCache()"), i)
        self.assertIn("PwbTimelagCache(i)%result%selected_lag",
                      body(MOD, "subroutine", "WritePwbPreassessTable"))

    def test_the_table_says_which_record_each_gas_is(self):
        self.assertIn("'date,time,gas,gas_record,method,detected_lag_s,'", MOD)
        self.assertIn("gas - firstGas + 1", body(MOD, "subroutine", "WriteRow"))

    def test_nothing_is_written_unless_asked(self):
        for name in ("WriteTlagPreassessTable", "WritePwbPreassessTable"):
            self.assertIn("if (.not. TOSetup%dump_lags) return",
                          body(MOD, "subroutine", name))


if __name__ == "__main__":
    unittest.main()
