"""The interface is told how far a run has got, split or not, so its estimate
of the time left holds.

The interface estimated the time left from the parent's own period lines. A
split pass prints only the parent's first piece - one of every few of the
range's periods - and then waits, saying no more than how many pieces are
done; so the estimate was the serial time while that piece ran, and then
stood still. Nor could the interface tell which passes lay ahead.

For the interface only (-c gui), the parent now says which passes the run
will walk:
    "  Run plan: to pf rp (17520 periods)"     (RP)
    "  Run plan: sa fx"                         (FCC)
and, every five seconds at most, how far the pass under way has got over
every process:
    "  Progress [rp]: 1234 of 17520 periods done."
Each worker notes its own progress beside its return code; the parent adds
those, the pieces done in full, and its own. None of it goes into the run log,
and a console run prints none of it.
"""

import re
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def read(rel):
    return (ROOT / rel).read_text(encoding="utf-8", errors="replace").replace("\r\n", "\n")


POOL = read("src/src_common/batch_pool.f90")
LOG = read("src/src_common/m_log.f90")
RP = read("src/src_rp/eddyflow-rp_main.f90")
FCC = read("src/src_fcc/eddyflow-fcc_main.f90")
PREPASS = read("src/src_rp/prepass_parallel.f90")


def body(text, kind, name):
    i = text.index("%s %s(" % (kind, name))
    return text[i:text.index("end %s %s" % (kind, name), i)]


class TheLinesAreForTheInterfaceOnly(unittest.TestCase):

    def test_console_say_writes_the_console_and_not_the_log(self):
        say = body(LOG, "subroutine", "ConsoleSay")
        self.assertIn("write(*, '(a)') text", say)
        self.assertNotIn("ulog", say)
        # A progress line may be open; this ends it first, as LogSay does.
        self.assertIn("call EndOpenLine()", say)

    def test_only_the_parent_of_a_gui_run_reports(self):
        begin = body(POOL, "subroutine", "ProgressBegin")
        self.assertIn("ProgOn = BatchIndex == 0 .and. EddyFlowProj%caller == 'gui'", begin)
        self.assertIn("if (BatchIndex == 0 .and. EddyFlowProj%caller == 'gui') then", RP)
        self.assertIn("if (BatchIndex == 0 .and. EddyFlowProj%caller == 'gui') then", FCC)

    def test_every_line_goes_through_console_say(self):
        self.assertIn("call ConsoleSay('  Progress [' // trim(tag)",
                      body(POOL, "subroutine", "ProgressSay"))
        self.assertNotIn("LogSay('  Progress", POOL)
        self.assertIn("call ConsoleSay('  Run plan:' // trim(RunPlan)", RP)
        self.assertIn("call ConsoleSay('  Run plan: sa fx')", FCC)

    def test_at_most_every_five_seconds(self):
        self.assertIn("now - ProgLastSay < 5 * rate", body(POOL, "subroutine", "ProgressSay"))


class EveryProcessIsCounted(unittest.TestCase):

    def test_a_worker_notes_its_progress_beside_its_return_code(self):
        note = body(POOL, "subroutine", "NotePieceProgress")
        self.assertIn("trim(BatchTmpDir) // slash // 'batch_'", note)
        self.assertIn("'.prog'", note)
        self.assertIn("now - PieceLastNote < 2 * rate", note)
        self.assertIn("ProgPath = trim(TmpDir) // 'batch_' // trim(tag) // '.prog'",
                      body(POOL, "function", "ProgPath"))

    def test_the_parent_adds_the_pieces_done_running_and_its_own(self):
        done = body(POOL, "function", "ProgressDone")
        #> Counted in walked periods: all of them, or the preliminary
        #> time-lag assessment's sample (SpanWork).
        self.assertIn("ProgressDone = SpanWork(ProgFirst, min(ProgOwnNext, ChunkEnd(1)))", done)
        self.assertIn("ProgressDone + SpanWork(ChunkStart(k), ChunkEnd(k))", done)
        self.assertIn("SpanWork(ChunkStart(k), ChunkStart(k) + PieceProg(k))", done)
        self.assertIn("write(TotalString, '(i0)') SpanWork(ProgFirst, ProgEnd)",
                      body(POOL, "subroutine", "ProgressSay"))
        self.assertIn("call ReadPieceProgress(k)", done)
        # A file caught mid-write keeps the last value.
        self.assertIn("PieceProg(k) = max(PieceProg(k),", body(POOL, "subroutine", "ReadPieceProgress"))

    def test_a_stale_progress_file_is_cleared_with_the_return_code(self):
        start = body(POOL, "subroutine", "StartPrepassBatches")
        self.assertIn("rcPath = ProgPath(kind, k)", start)

    def test_the_wait_reports_and_completes_the_pass(self):
        wait = body(POOL, "subroutine", "WaitPrepassBatches")
        self.assertIn("call ProgressSay(.false.)", wait)
        self.assertIn("if (ProgressUsesPool() .and. kind /= 'pd') then", wait)

    def test_the_detection_round_is_tagged_apart(self):
        say = body(POOL, "subroutine", "ProgressSay")
        self.assertIn("if (PoolKind == 'pd') tag = 'pd'", say)

    def test_top_up_serves_workers_and_parent_alike(self):
        top = body(POOL, "subroutine", "TopUpPrepassBatches")
        self.assertIn("call NotePieceProgress(next)", top)
        self.assertIn("call ProgressAt(next)", top)


class EveryPassReports(unittest.TestCase):

    def test_rp_passes(self):
        self.assertIn("call ProgressBegin('to', toStartTimestampIndx, toEndTimestampIndx, 'periods')", RP)
        self.assertIn("call ProgressBegin('pf', pfStartTimestampIndx, pfEndTimestampIndx, 'periods')", RP)
        self.assertIn("call ProgressBegin('rp', rpStartTimestampIndx, rpEndTimestampIndx, 'periods')", RP)
        self.assertEqual(RP.count("call ProgressEnd()"), 3)
        self.assertIn("public :: ProgressBegin, ProgressEnd", PREPASS)

    def test_the_pre_pass_begins_over_the_whole_range_before_it_splits(self):
        self.assertLess(RP.index("call ProgressBegin('to',"), RP.index("call StartPrepassBatches('to',"))
        self.assertLess(RP.index("call ProgressBegin('pf',"), RP.index("call StartPrepassBatches('pf',"))

    def test_every_loop_reports_where_it_stands(self):
        self.assertEqual(RP.count("call TopUpPrepassBatches(pcount + 1)"), 2)
        # Production: unconditional now, workers and serial parent too.
        self.assertIn("        call TopUpPrepassBatches(pcount)\n", RP)
        self.assertNotIn("if (ProdSplit) call TopUpPrepassBatches", RP)

    def test_fcc_passes(self):
        self.assertIn("call ProgressBegin('sa', saStartTimestampIndx, saEndTimestampIndx + 1, 'files')", FCC)
        self.assertIn("call ProgressAt(fcount)", FCC)
        self.assertIn("call ProgressBegin('fx', 1, NumExRecords + 1, 'records')", FCC)
        i = FCC.index("ex_loop: do i = 1, NumExRecords")
        loop = FCC[i:FCC.index("call ReadExRecord('', uex, -1,", i)]
        self.assertIn("call TopUpPrepassBatches(i)", loop)
        self.assertEqual(FCC.count("call ProgressEnd()"), 2)

    def test_the_rp_plan_lists_the_passes_the_run_will_walk(self):
        i = RP.index("RunPlan = ''")
        plan = RP[i:RP.index("call ConsoleSay('  Run plan:'", i)]
        for stage in ("' to'", "' pf'", "' rp'"):
            self.assertIn(stage, plan)
        self.assertIn("RPsetup%to_onthefly", plan)
        self.assertIn("RPsetup%pf_onthefly", plan)
        self.assertIn("if (.not. AssessmentOnly) RunPlan = trim(RunPlan) // ' rp'", plan)
        # Before the first pass begins.
        self.assertLess(i, RP.index("call ProgressBegin('to',"))


if __name__ == "__main__":
    unittest.main()
