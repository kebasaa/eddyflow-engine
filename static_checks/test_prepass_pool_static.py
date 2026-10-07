"""A parallel pre-pass hands out pieces of similar work to whichever worker is free.

The range used to be cut into one slice per worker, by period count. On the
Yatir run (6664 half-hours, 3387 raw files) that gave one slice 801 files and
another none: the data thin out from June and stop for two weeks in July. The
worker with no data finished in 79 s and idled; the run waited on the
heaviest. And on a machine with fast and slow cores, a slice that lands on a
slow one finishes last whatever its size.

So the range is now cut into several pieces per worker, at equal shares of
estimated work - a period weighs 1 if a raw file covers it, almost nothing if
none does - and a fixed number run at once. The parent runs piece 1 (it has
to, for the state the loop establishes), then dispatches: each worker that
finishes makes room for the next piece. A fast core simply finishes more.

What does not change: the pieces tile the range, are merged in piece order,
and a failure stops the run - now at once rather than after the rest. Results
are identical to a serial run (``check_parallel.sh``). Without parallel
processing (``-j 1``, the interface's tickbox off) none of this runs.

Part of the EddyFlow engine's static checks.
"""

import re
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / "src"


def read(rel):
    return (SRC / rel).read_text(encoding="utf-8", errors="replace")


def code(text):
    return "\n".join(ln for ln in text.splitlines() if not ln.lstrip().startswith("!"))


# The worker pool itself moved to src_common/batch_pool.f90, shared with FCC;
# what is checked here is read from both.
PARALLEL = code(read("src_rp/prepass_parallel.f90") + read("src_common/batch_pool.f90"))
MAIN = code(read("src_rp/eddyflow-rp_main.f90"))


def body(name, kind="subroutine"):
    return PARALLEL[PARALLEL.index("%s %s(" % (kind, name)):
                    PARALLEL.index("end %s %s" % (kind, name))]


PLAN = body("PlanPrepassChunks")
START = body("StartPrepassBatches")
WAIT = body("WaitPrepassBatches")
POLL = body("PollPrepassPool")
TOPUP = body("TopUpPrepassBatches")
LAUNCH = body("LaunchChunks")


class PiecesAreCutByWorkNotByCount(unittest.TestCase):

    def test_a_period_with_a_raw_file_weighs_one_and_an_empty_one_little(self):
        self.assertRegex(PARALLEL, r"EmptyPeriodWeight = 0\.\d+d0")
        self.assertIn("weight(p) = EmptyPeriodWeight", PLAN)
        self.assertIn("weight(p) = 1d0", PLAN)

    def test_coverage_comes_from_the_file_list_in_one_sweep(self):
        """Files are in time order; a file covers a period if it starts
        before the period ends and ends after it starts."""
        self.assertIn("Files(j)%timestamp + DatafileDateStep > Series(p)", PLAN)
        self.assertIn("Files(j)%timestamp < Series(p + 1)", PLAN)
        self.assertEqual(PLAN.count("do p = iStart, iEnd - 1"), 1,
                         "one pass over the periods, not one per file")

    def test_cuts_fall_at_equal_shares_of_the_total(self):
        self.assertIn("total = sum(weight)", PLAN)
        self.assertIn("cum >= total * dble(nCuts + 1) / dble(want)", PLAN)

    def test_several_pieces_per_worker_within_the_file_name_width(self):
        self.assertRegex(PARALLEL, r"integer, parameter :: ChunksPerWorker = \d+")
        self.assertIn("integer, parameter :: MaxChunks = 99", PARALLEL)
        self.assertIn("min(ChunksPerWorker * nEff, MaxChunks,", PLAN)
        #> Two digits in every per-piece file name; 99 is the most that fits.
        self.assertIn("'_b', k", PARALLEL)
        self.assertIn("i2.2)') trim(kind), '_b', k", PARALLEL)

    def test_piece_one_reaches_the_first_period_with_a_file(self):
        """The parent runs piece 1 and must read a raw file in it, or the
        finalisation runs on state the loop never set."""
        self.assertIn("if (firstData == 0) firstData = p", PLAN)
        self.assertIn("if (cuts(1) > firstData) exit", PLAN)

    def test_the_count_passed_to_each_worker_is_the_piece_count(self):
        """--batch k:count:... is validated as k <= count."""
        self.assertIn("call WriteChildScript(kind, k, NumChunks, ChunkStart(k), ChunkEnd(k),",
                      START)


class AFreeWorkerTakesTheNextPiece(unittest.TestCase):

    def test_while_the_parent_runs_piece_one_the_others_start(self):
        self.assertIn("call LaunchChunks(kind, min(nEff - 1, NumChunks - 1))", START)

    def test_the_wait_keeps_neff_workers_running_while_pieces_are_left(self):
        loop = WAIT[WAIT.index("        do\n"):WAIT.index("deallocate(Running)")]
        self.assertIn("call PollPrepassPool(nEff, .false., progressed)", loop)
        self.assertIn("nRunning = count(Running)", POLL)
        self.assertIn("if (nRunning < slots .and. NextChunk <= NumChunks) then", POLL)
        self.assertIn("call LaunchChunks(PoolKind, min(slots - nRunning, NumChunks - NextChunk + 1))",
                      POLL)

    def test_a_piece_is_launched_once(self):
        self.assertIn("NextChunk = NextChunk + n", LAUNCH)
        self.assertIn("do k = NextChunk, NextChunk + n - 1", LAUNCH)

    def test_a_half_written_return_code_reads_as_not_yet(self):
        rc = body("ChunkReturnCode")
        i = rc.index("read(u, *, iostat = io_status) rc")
        self.assertIn("if (io_status /= 0) then", rc[i:])
        self.assertLess(rc[i:].index("return"), rc[i:].index("finished = .true."))

    def test_the_wait_is_bounded_by_progress_not_by_total_time(self):
        """A season legitimately takes hours; what must not happen is a day
        with no piece finishing."""
        self.assertIn("progressed = .true.", POLL[POLL.index("Running(k) = .false."):])
        self.assertIn("if (progressed) ticks = 0", WAIT)
        self.assertIn("if (ticks > MaxWaitTicks) then", WAIT)


class AFailureStillStopsTheRun(unittest.TestCase):

    def test_at_once_inside_the_dispatch_loop(self):
        self.assertIn("error stop 'A parallel pre-pass worker failed.'", POLL)
        self.assertIn("error stop 'A parallel pre-pass worker produced no output.'", POLL)

    def test_a_failure_noticed_mid_period_ends_the_open_line_first(self):
        i = POLL.index("if (.not. delivered) then")
        self.assertLess(POLL.index("if (quiet) call LogSay('')", i),
                        POLL.index("call AppendWorkerLog(PoolKind, k)", i))

    def test_worker_logs_join_the_run_log_in_piece_order(self):
        tail = WAIT[WAIT.index("deallocate(Running)"):]
        self.assertIn("do k = 2, NumChunks", tail)
        self.assertIn("call AppendWorkerLog(kind, k)", tail)


class ProgressIsReported(unittest.TestCase):

    def test_each_finished_piece_is_counted_once(self):
        self.assertIn("' pieces done.'", body("SayPiecesDone"))
        i = POLL.index("PoolDone = PoolDone + 1")
        self.assertLess(i, POLL.index("if (.not. quiet) call SayPiecesDone()"))

    def test_the_parent_s_own_piece_counts_once_it_waits(self):
        i = WAIT.index("call LogSay('  Waiting for the workers:')")
        self.assertLess(i, WAIT.index("PoolDone = PoolDone + 1"))
        self.assertIn("if (PoolDone > 1) call SayPiecesDone()", WAIT)


class NoCoreIdlesWhileTheParentWorks(unittest.TestCase):
    """A worker that finishes while the parent is still on piece 1 is given
    the next piece then, not when the parent reaches the wait."""

    def test_the_parent_tops_up_quietly_leaving_its_own_core(self):
        self.assertIn("if (.not. PoolActive) return", TOPUP)
        self.assertIn("call PollPrepassPool(PoolEff - 1, .true., progressed)", TOPUP)

    def test_every_loop_the_parent_shares_tops_up_once_a_period(self):
        for loop in ("to_periods_loop: do", "pf_periods_loop: do"):
            i = MAIN.index(loop)
            self.assertIn("call TopUpPrepassBatches()", MAIN[i:i + 600], loop)
        i = MAIN.index("if (ProdSplit .and. pcount >= ProdParentEnd) exit periods_loop")
        self.assertIn("if (ProdSplit) call TopUpPrepassBatches()", MAIN[i:i + 200])

    def test_the_pool_is_set_up_when_the_workers_start(self):
        self.assertLess(START.index("allocate(Running(NumChunks))"),
                        START.index("call LaunchChunks(kind, min(nEff - 1, NumChunks - 1))"))
        self.assertIn("Running(2:NextChunk - 1) = .true.", START)
        self.assertIn("PoolActive = .true.", START)
        self.assertIn("PoolActive = .false.", WAIT)


class NothingOfThisRunsWithoutParallelProcessing(unittest.TestCase):

    def test_one_job_means_no_pieces_and_no_workers(self):
        plan = body("PlanPrepassBatches")
        self.assertIn("if (requested <= 1) return", plan)
        for kind, flag in (("to", "toParallel = toWorkers > 1"),
                           ("pf", "pfParallel = pfWorkers > 1")):
            self.assertIn(flag, MAIN)
            i = MAIN.index(flag)
            self.assertIn("call StartPrepassBatches('%s'" % kind, MAIN[i:i + 300])

    def test_the_old_equal_count_slicer_is_gone(self):
        self.assertNotIn("PrepassSlice", PARALLEL)
        self.assertNotIn("PrepassSlice", MAIN)


if __name__ == "__main__":
    unittest.main()
