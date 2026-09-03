"""Static checks for the parallel assessment pre-passes.

The planar fit and the time-lag optimiser each walk every period in their
range and append one independent record to a flat array, so the range can be
cut into slices and each slice run in a copy of ``eddyflow_rp``. Four things
about that arrangement are load-bearing and none of them fails loudly.

**A worker must never spawn workers.** The growth is exponential: the first
attempt at this shipped a command line whose switches were silently dropped,
so every child believed it was a parent, and the machine was carrying 1909
processes within a minute. ``BatchIndex > 0`` is the ordinary guard; the
interlock behind it asks the *raw command line* instead, because the parsed
result is exactly what fails when switch handling is wrong.

**Switches must only consume a value when they take one.** The argument loop
used to read a value after every token, so the project path swallowed
whatever followed it - which is why a ``-e`` written after the path was
silently ignored, and why ``-j`` written there did nothing at all.

**The slices are half-open.** Both period loops exit on ``pcount >= endIndex``,
so an end index is one past the last period processed. Treating it as
inclusive lost one period per worker, which showed up only as a slightly
different planar fit.

**The parent runs a slice itself.** The code after the period loop reads
global state the loop established - ``SortWindBySector`` takes the north
offset out of ``E2Col``, which is filled when a raw file's metadata is read.
A parent that had skipped the loop reached that code with the offset unset
and binned every period into a sector rotated by the site's true offset, 209
degrees at the test site. Nothing enumerates what the finalisation depends
on, so the parent must not skip the loop.

**The PWB cache pre-pass is not split at all.** Its classifier decides a
period partly from the last settled detection before it, and that chain has
no time limit. A lead-in was tried and does not save it: over two days of
CH-LAE, COS reached no settled detection in 103 periods, so no lead-in of any
length rebuilds its state. Splitting it left every applied lag identical but
moved donor attribution at the seams, and the aggregate summary resolved from
those donor votes moved with it.
"""

import re
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / "src"


def read(rel):
    return (SRC / rel).read_text(encoding="utf-8", errors="replace")


PARALLEL = read("src_rp/prepass_parallel.f90")
MAIN = read("src_rp/eddyflow-rp_main.f90")
ENV = read("src_common/init_env.f90")
RUNLOG = read("src_common/init_run_log.f90")
GLOBALS = read("src_common/m_common_global_var.f90")
OSENV = read("src_common/set_os_environment.f90")
PWB = read("src_rp/pwb_timelag_handle.f90")


class AWorkerNeverSpawnsWorkers(unittest.TestCase):
    """The one failure mode that takes the machine down with it."""

    def test_the_ordinary_guard_is_the_first_thing_the_planner_does(self):
        body = PARALLEL[PARALLEL.index("subroutine PlanPrepassBatches"):
                        PARALLEL.index("end subroutine PlanPrepassBatches")]
        guard = body.index("if (BatchIndex > 0) return")
        detect = body.index("DetectCoreCount()")
        self.assertLess(
            guard, detect,
            "PlanPrepassBatches must refuse before it works out a worker count")

    def test_the_interlock_reads_the_raw_command_line(self):
        body = PARALLEL[PARALLEL.index("subroutine PlanPrepassBatches"):
                        PARALLEL.index("end subroutine PlanPrepassBatches")]
        self.assertIn("call get_command(cmdline)", body)
        self.assertIn("index(cmdline, '--batch')", body)
        self.assertIn("error stop", body[body.index("index(cmdline, '--batch')"):])

    def test_the_interlock_precedes_any_launch(self):
        body = PARALLEL[PARALLEL.index("subroutine PlanPrepassBatches"):
                        PARALLEL.index("end subroutine PlanPrepassBatches")]
        self.assertLess(body.index("index(cmdline, '--batch')"),
                        body.index("DetectCoreCount()"))

    def test_children_are_told_to_stay_serial(self):
        self.assertIn("' -j 1'", PARALLEL,
                      "a worker must carry -j 1 as well as --batch")

    def test_the_worker_count_is_capped(self):
        self.assertRegex(PARALLEL, r"integer, parameter :: MaxWorkers = \d+")
        self.assertIn("min(requested, MaxWorkers)", PARALLEL)


class SwitchesOnlyConsumeAValueWhenTheyTakeOne(unittest.TestCase):
    """Otherwise the project path swallows whatever follows it."""

    #: Every switch InitEnv reads a value for.
    VALUED = ("-s", "--system", "-e", "--environment", "-m", "--mode",
              "-c", "--caller", "-j", "--jobs",
              "--batch", "--batch-out")

    def test_the_loop_asks_before_reading_a_value(self):
        loop = ENV[ENV.index("arg_loop: do"):ENV.index("end do arg_loop")]
        self.assertIn("if (SwitchTakesValue(switch)) then", loop)
        self.assertEqual(
            loop.count("call get_command_argument(i, value=arg"), 1,
            "the value is read in exactly one place, behind the test")

    def test_every_valued_switch_is_declared_as_such(self):
        body = ENV[ENV.index("logical function SwitchTakesValue"):
                   ENV.index("end function SwitchTakesValue")]
        for switch in self.VALUED:
            self.assertIn("'%s'" % switch, body,
                          "%s takes a value but is not listed" % switch)

    def test_help_and_version_take_no_value(self):
        body = ENV[ENV.index("logical function SwitchTakesValue"):
                   ENV.index("end function SwitchTakesValue")]
        true_arm = body[:body.index("case default")]
        for switch in ("'-h'", "'--help'", "'-v'", "'--version'"):
            self.assertNotIn(switch, true_arm)

    def test_the_lowercase_loop_has_its_own_counter(self):
        """It used to run on `i` and leave the argument loop past the end."""
        loop = ENV[ENV.index("arg_loop: do"):ENV.index("end do arg_loop")]
        self.assertNotIn("do i = 1, len_trim(lowerPath)", loop)
        self.assertIn("do ch = 1, len_trim(lowerPath)", loop)

    def test_jobs_is_documented(self):
        self.assertRegex(ENV, r"\[-j \| --jobs")

    def test_the_batch_switches_are_not_documented(self):
        """They are written by the parent, not typed by anyone."""
        help_body = ENV[ENV.index("subroutine CommandLineHelp"):]
        self.assertNotIn("--batch", help_body)


class SlicesTileTheRangeExactly(unittest.TestCase):
    """Half-open, because that is what the period loops are."""

    def test_the_slice_arithmetic_is_half_open(self):
        body = PARALLEL[PARALLEL.index("subroutine PrepassSlice"):
                        PARALLEL.index("end subroutine PrepassSlice")]
        self.assertIn("total = iEnd - iStart", body)
        self.assertNotIn("total = iEnd - iStart + 1", body)
        self.assertIn("sliceEnd = sliceStart + len", body)
        self.assertNotIn("sliceEnd = sliceStart + len - 1", body)

    def test_the_tiling_is_asserted_before_anything_is_launched(self):
        body = PARALLEL[PARALLEL.index("subroutine StartPrepassBatches"):
                        PARALLEL.index("end subroutine StartPrepassBatches")]
        self.assertIn("if (covered /= iEnd - iStart) &", body)
        self.assertLess(body.index("covered /= iEnd - iStart"),
                        body.index("call system(trim(cmd))"))

    def test_no_slice_reads_outside_its_own_range(self):
        """There is no lead-in: nothing that is split carries state."""
        body = PARALLEL[PARALLEL.index("subroutine PrepassSlice"):
                        PARALLEL.index("end subroutine PrepassSlice")]
        self.assertNotIn("warmup", body)
        self.assertIn("covered", PARALLEL)


class TheParentRunsASliceItself(unittest.TestCase):
    """So it reaches the finalisation with the state the loop establishes."""

    def test_neither_period_loop_is_skipped(self):
        """The loops used to sit in an else arm the parent took."""
        for loop in ("to_periods_loop: do", "pf_periods_loop: do"):
            before = MAIN[:MAIN.index(loop)]
            tail = before[-400:]
            self.assertNotRegex(
                tail, r"\n\s+else\s*\n\s*$",
                "%s must not be reachable only through an else arm" % loop)

    def test_the_parent_takes_slice_one(self):
        for kind, first in (("'to'", "toWorkers"), ("'pf'", "pfWorkers")):
            i = MAIN.index("call StartPrepassBatches(%s" % kind)
            chunk = MAIN[i:i + 500]
            self.assertIn("call PrepassSlice(", chunk)
            self.assertIn("%s, 1, " % first, chunk,
                          "the parent must ask for slice 1 of %s" % kind)

    def test_workers_are_launched_before_the_parent_starts_its_slice(self):
        for kind in ("'to'", "'pf'"):
            self.assertLess(MAIN.index("call StartPrepassBatches(%s" % kind),
                            MAIN.index("call WaitPrepassBatches(%s" % kind))

    def test_the_merge_starts_at_the_second_slice(self):
        for name in ("MergeTlagBatchDumps", "MergePfBatchDumps"):
            body = PARALLEL[PARALLEL.index("subroutine %s" % name):
                            PARALLEL.index("end subroutine %s" % name)]
            self.assertIn("do k = 2, nEff", body,
                          "%s must not re-read the parent's own slice" % name)
            self.assertNotIn("do k = 1, nEff", body)

    def test_the_merge_appends_rather_than_resets(self):
        for name in ("MergeTlagBatchDumps", "MergePfBatchDumps"):
            body = PARALLEL[PARALLEL.index("subroutine %s" % name):
                            PARALLEL.index("end subroutine %s" % name)]
            self.assertIn("intent(inout) :: n", body,
                          "%s must add to the parent's own count" % name)

    def test_the_launcher_only_starts_the_other_slices(self):
        body = PARALLEL[PARALLEL.index("subroutine StartPrepassBatches"):
                        PARALLEL.index("end subroutine StartPrepassBatches")]
        self.assertIn("call WriteChildScript(", body)
        i = body.index("call WriteChildScript(")
        self.assertIn("do k = 2, nEff", body[:i])


class AWorkerStopsBeforeItCanWriteAnything(unittest.TestCase):
    """It hands back records; the run's output belongs to the parent."""

    def test_each_prepass_dumps_then_stops(self):
        for dump in ("call WriteTlagBatchDump(", "call WritePfBatchDump("):
            i = MAIN.index(dump)
            chunk = MAIN[i:i + 700]
            self.assertIn("stop ''", chunk,
                          "%s must be followed by a stop" % dump)

    def test_a_worker_stops_before_the_output_files_are_opened(self):
        last_stop = max(MAIN.index("call WriteTlagBatchDump("),
                        MAIN.index("call WritePfBatchDump("))
        self.assertLess(last_stop, MAIN.index("call InitOutFiles_rp"))

    def test_a_worker_gets_the_resolved_project_path(self):
        """An EddyPro project is imported once, by the parent, not raced for."""
        body = PARALLEL[PARALLEL.index("subroutine WriteChildScript"):
                        PARALLEL.index("end subroutine WriteChildScript")]
        self.assertIn("trim(PrjPath)", body)

    def test_switches_precede_the_project_path_on_the_child_command_line(self):
        body = PARALLEL[PARALLEL.index("subroutine WriteChildScript"):
                        PARALLEL.index("end subroutine WriteChildScript")]
        cmd = body[body.index("cmd = "):body.index("open(newunit = u")]
        self.assertLess(cmd.index("--batch-out"), cmd.index("trim(PrjPath)"),
                        "the project path must come last")

    def test_a_worker_is_told_to_report_as_a_console_caller(self):
        body = PARALLEL[PARALLEL.index("subroutine WriteChildScript"):
                        PARALLEL.index("end subroutine WriteChildScript")]
        self.assertIn("' -c console'", body)


class WorkersDoNotCollideWithTheParentOrEachOther(unittest.TestCase):
    """They share their parent's start second, so the stamp is not enough."""

    def test_the_temporary_directory_carries_the_batch_index(self):
        self.assertIn("tmpDirPadding = trim(tmpDirPadding) // trim(batchPadding)",
                      ENV)

    def test_the_output_stamp_does_not(self):
        """It is fixed-width and concatenated whole into every file name."""
        self.assertNotIn(
            "Timestamp_FilePadding = trim(Timestamp_FilePadding) // trim(batchPadding)",
            ENV)
        self.assertRegex(GLOBALS, r"character\(22\)\s+:: Timestamp_FilePadding")

    def test_a_worker_logs_beside_its_records(self):
        self.assertIn("LogPath = BatchOutPath(1:len_trim(BatchOutPath)) // LogExt",
                      RUNLOG)
        self.assertIn("if (BatchIndex > 0) then", RUNLOG)

    def test_the_parent_folds_worker_logs_into_the_run_log(self):
        self.assertIn("call AppendWorkerLog(", PARALLEL)


class TheWaitIsPortableAndBounded(unittest.TestCase):
    """gfortran's SLEEP is a GNU extension and this builds -std=f2008."""

    def test_the_pause_is_a_shell_command(self):
        self.assertIn("comm_sleep", GLOBALS)
        self.assertEqual(OSENV.count("comm_sleep        ="), 3,
                         "one per operating system arm")
        self.assertIn("call system(comm_sleep)", PARALLEL)

    def test_the_wait_gives_up_eventually(self):
        body = PARALLEL[PARALLEL.index("subroutine WaitPrepassBatches"):
                        PARALLEL.index("end subroutine WaitPrepassBatches")]
        self.assertIn("if (ticks > MaxWaitTicks) then", body)
        self.assertIn("error stop", body)

    def test_stale_return_codes_are_cleared_before_launching(self):
        body = PARALLEL[PARALLEL.index("subroutine StartPrepassBatches"):
                        PARALLEL.index("end subroutine StartPrepassBatches")]
        self.assertIn("comm_del", body)
        self.assertLess(body.index("comm_del"), body.index("call system(trim(cmd))"))


class AFailedWorkerStopsTheRun(unittest.TestCase):
    """A partial fit is a wrong answer that looks like a right one."""

    def test_a_non_zero_return_code_is_fatal(self):
        body = PARALLEL[PARALLEL.index("subroutine WaitPrepassBatches"):
                        PARALLEL.index("end subroutine WaitPrepassBatches")]
        self.assertIn("if (rc /= 0) then", body)
        self.assertIn("error stop 'A parallel pre-pass worker failed.'", body)

    def test_a_missing_dump_is_fatal(self):
        body = PARALLEL[PARALLEL.index("subroutine WaitPrepassBatches"):
                        PARALLEL.index("end subroutine WaitPrepassBatches")]
        self.assertIn("exited cleanly but wrote no records", body)

    def test_the_diagnostic_points_at_what_is_actually_there(self):
        """A worker that died never closed its log, so the console capture is"""
        """the only record of what it managed to say - and it follows, not"""
        """precedes, the message."""
        body = PARALLEL[PARALLEL.index("subroutine WaitPrepassBatches"):
                        PARALLEL.index("end subroutine WaitPrepassBatches")]
        self.assertNotIn("is above", body)
        i = body.index("before it stopped:")
        self.assertLess(i, body.index("call DumpWorkerStdout(", i - 400))

    def test_the_failing_worker_s_output_is_shown(self):
        body = PARALLEL[PARALLEL.index("subroutine WaitPrepassBatches"):
                        PARALLEL.index("end subroutine WaitPrepassBatches")]
        self.assertEqual(body.count("call DumpWorkerStdout("), 2)

    def test_a_foreign_dump_is_refused(self):
        self.assertIn("if (magic /= BatchMagic) &", PARALLEL)
        #> One per merge, and there are three of them now: time-lag records,
        #> planar-fit means, and the PWB cache rows. A merge that skipped the
        #> check would read a dump from a crashed earlier run as data.
        self.assertEqual(PARALLEL.count("if (magic /= BatchMagic) &"), 3,
                         "one per merge")


class ThePwbCachePrepassIsSplitToo(unittest.TestCase):
    """This class used to assert the opposite, and the reasoning it gave was
    sound about the wrong thing.

    The streaming classifier does decide a period partly from the last settled
    detection before it, and that chain has no time limit, and the weak species
    never settle - over two days of CH-LAE, COS reached no settled detection at
    all. All true, and none of it decides the output of a cache-generation run.
    `PostProcessPwbTimelagCache` does, once the whole run has been read, and it
    overrules every row.

    What actually blocked the split was four values carrying the classifier's
    verdict into the settled table anyway. Each now comes from the table, and
    each has its own file:

      * the terminal fallback's lag - `test_pwb_terminal_fallback_static`
      * the aggregate dataset's membership - `test_pwb_aggregate_dataset_static`
      * the donor tally - `test_pwb_donor_tally_static`
      * every field saying how a period was settled -
        `test_pwb_postpass_owns_its_fields_static`

    The transport itself is pinned in `test_pwb_prepass_split_static`, and the
    end-to-end claim is `check_parallel.sh base_pwb_par.eddyflow`: 24 periods
    across 6 workers, 85 files byte-identical to the same run under -j 1.
    """

    def test_the_planner_still_honours_a_caller_that_refuses(self):
        """The switch remains, and remains ahead of the core count - a worker
        must never spawn workers of its own, whatever the caller allows."""
        body = PARALLEL[PARALLEL.index("subroutine PlanPrepassBatches"):
                        PARALLEL.index("end subroutine PlanPrepassBatches")]
        self.assertIn("if (.not. allowed) return", body)
        self.assertLess(body.index("if (.not. allowed) return"),
                        body.index("DetectCoreCount()"))

    def test_the_time_lag_prepass_is_no_longer_gated_on_the_cache_mode(self):
        i = MAIN.index("call PlanPrepassBatches(toEndTimestampIndx")
        self.assertNotIn(".not. PwbCacheGenerate", MAIN[i:i + 200])
        self.assertIn(".true.", MAIN[i:i + 200])

    def test_the_planar_fit_is_always_allowed(self):
        """A planar-fit period contributes three means and carries nothing."""
        i = MAIN.index("call PlanPrepassBatches(pfEndTimestampIndx")
        self.assertIn(".true.", MAIN[i:i + 200])

    def test_the_cache_now_does_cross_the_process_boundary(self):
        """It is the evidence the parent settles from, so it has to."""
        self.assertIn("PwbTimelagCache", PARALLEL)
        self.assertIn("PwbOptDate", PARALLEL)

    def test_but_the_diagnostics_still_do_not(self):
        """They are recounted from the settled table, so transporting them
        would be carrying a number that is about to be thrown away."""
        for name in ("PwbSummaryDonorCount", "GetPwbDiagnostics",
                     "AddPwbDiagnostics"):
            self.assertNotIn(name, PARALLEL,
                             "%s should not be reachable here" % name)

    def test_the_accessors_added_for_the_first_attempt_are_still_gone(self):
        for name in ("StorePwbTimelagCacheEntry", "GetPwbDiagnostics",
                     "AddPwbDiagnostics", "nPwbDiagnostics"):
            self.assertNotIn(name, PWB, "%s is dead code" % name)

    def test_the_dump_format_was_versioned_when_it_changed(self):
        """04 carries the PWB sections; an 03 dump left by an older build has
        none of them and must be refused rather than read short."""
        self.assertIn("BatchMagic = 'EDDYFLOW_PREPASS_04 '", PARALLEL)


if __name__ == "__main__":
    unittest.main()
