"""Static checks for the output-disk guard and the progress-aware watchdog.

On CH-LAE the network share a run wrote to filled up while 24 workers wrote
raw datasets at all seven levels, 77 MB a period. From then on every file came
out empty and nothing reported an error - the share accepted the writes and
lost them - while each period took sixteen minutes instead of twenty seconds.
The parent, which only reset its one-day timeout when a whole piece finished,
then stopped the run a day later, with every worker still going.

Three things are pinned.

**The disk is watched, not the writes.** Every process writing a pass looks
at the free space of the folders it writes to at the top of every period and
stops below MinFreeOutputBytes - a production worker at its parent's folder
too, since its per-period files go there.

**A run that cannot fit is refused before any worker starts.** The parent
measures the raw datasets of the first period it writes, projects them over
the periods still to come, and stops before the production pass is split.

**The watchdog counts progress, not only finished pieces.** A worker that gets
one period further resets the wait, so a slow pass is told from a stuck one,
and a stuck one is reported with where each worker had got to.
"""

import re
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / "src"


def read(rel):
    return (SRC / rel).read_text(encoding="utf-8", errors="replace")


POOL = read("src_common/batch_pool.f90")
WIN = read("src_os/win/m_process_os.f90")
POSIX = read("src_os/posix/m_process_os.f90")
RP = read("src_rp/eddyflow-rp_main.f90")
RAW = read("src_rp/out_raw_data.f90")
FCC = read("src_fcc/eddyflow-fcc_main.f90")
PROD = read("src_rp/production_parallel.f90")


def body(src, kind, name):
    m = re.search(r"\n\s*(?:[\w()]+\s+)*%s\s+%s\b(.*?)\n\s*end %s %s"
                  % (kind, name, kind, name), src, re.S | re.I)
    assert m, "%s %s not found" % (kind, name)
    return m.group(1)


class TheDiskIsAskedOnBothPlatforms(unittest.TestCase):

    def test_both_os_files_offer_disk_free_bytes(self):
        for src in (WIN, POSIX):
            self.assertRegex(src, r"public :: ProcessAlive, DiskFreeBytes")
            self.assertIn("integer(8) function DiskFreeBytes(dir)", src)

    def test_windows_asks_for_what_the_user_may_write(self):
        #> The first of the three counts, which honours a quota; the total
        #> free on the volume does not.
        self.assertIn("bind(C, name = 'GetDiskFreeSpaceExA')", WIN)
        fn = body(WIN, "function", "DiskFreeBytes")
        self.assertIn("DiskFreeBytes = int(avail, 8)", fn)

    def test_posix_reads_df_from_the_capacity_column_back(self):
        fn = body(POSIX, "function", "DiskFreeBytes")
        self.assertIn("df -Pk", fn)
        self.assertIn("line(i:i) == '%'", fn)

    def test_an_unanswerable_disk_stops_nothing(self):
        for src in (WIN, POSIX):
            self.assertIn("DiskFreeBytes = -1", body(src, "function", "DiskFreeBytes"))
        self.assertIn("if (free < 0 .or. free >= MinFreeOutputBytes) return",
                      body(POOL, "subroutine", "StopIfDirLow"))


class EveryWriterWatchesTheDisk(unittest.TestCase):

    def test_the_floor_is_a_gibibyte(self):
        #> tests/test_output_space_static.py in the GUI holds its
        #> MinFreeBytes to this number.
        self.assertIn("integer(8), parameter :: MinFreeOutputBytes = 1073741824_8", POOL)

    def test_it_stops_with_a_stop_code_the_interface_recognises(self):
        fn = body(POOL, "subroutine", "StopIfDirLow")
        self.assertIn("The disk this run writes to has only", fn)
        self.assertIn("error stop 'The output disk is full.'", fn)

    def test_rp_checks_every_period_before_it_writes(self):
        #> After the worker skip-ahead, so a worker checks only once it is
        #> writing; with its parent's folder, where its per-period files go.
        call = "call StopIfOutputDiskFull(Dir%main_out, ProdParentOutDir())"
        self.assertEqual(RP.count(call), 1)
        at = RP.index(call)
        before = RP[:at]
        #> In the production loop: after its batch block, before its logging.
        self.assertGreater(at, RP.index("\n    periods_loop: do"))
        self.assertTrue(before.rstrip().endswith("!> production worker's per-period files go to its parent's folder."))
        self.assertIn("call ResetPwbOutsideWindow()", before[-1200:])
        self.assertIn("processing new", RP[at:at + 800])

    def test_a_worker_knows_its_parents_folder(self):
        fn = body(PROD, "function", "ProdParentOutDir")
        self.assertIn("ProdParentOutDir = ParentMainOut", fn)

    def test_fcc_checks_every_record(self):
        loop = FCC[FCC.index("ex_loop: do"):]
        self.assertLess(loop.index("call StopIfOutputDiskFull(Dir%main_out)"),
                        loop.index("call ReadExRecord("))

    def test_the_check_is_rate_limited(self):
        fn = body(POOL, "subroutine", "StopIfOutputDiskFull")
        self.assertIn("now - LastDiskCheck < 2 * rate", fn)


class ARunThatCannotFitIsRefused(unittest.TestCase):

    def test_raw_datasets_are_measured_as_written(self):
        tail = RAW[RAW.index("close(udf)"):]
        self.assertIn("call NoteRawOutputBytes(RawPath)", tail)
        fn = body(POOL, "subroutine", "NoteRawOutputBytes")
        self.assertIn("if (BatchIndex > 0 .or. OutputFitChecked) return", fn)

    def test_the_projection_comes_before_the_split(self):
        self.assertLess(RP.index("if (BatchIndex == 0) call CheckRawOutputWillFit()"),
                        RP.index("call TryProdSplit()"))

    def test_it_projects_over_the_periods_that_have_data(self):
        fn = body(RP, "subroutine", "CheckRawOutputWillFit")
        self.assertIn("call PeriodsWithFile(pcount + 1, rpEndTimestampIndx", fn)
        self.assertIn("nLeft = count(has)", fn)

    def test_it_stops_only_when_the_projection_will_not_fit(self):
        fn = body(POOL, "subroutine", "CheckOutputWillFit")
        self.assertIn("need = perPeriod * max(0, nLeft)", fn)
        self.assertIn("if (need <= free - MinFreeOutputBytes) return", fn)
        self.assertIn("will not fit on the", fn)
        self.assertIn("error stop 'Not enough room on the output disk for the raw datasets.'", fn)


class TheWatchdogCountsProgress(unittest.TestCase):

    def test_every_wait_resets_on_worker_progress(self):
        for name in ("WaitPrepassBatches", "WaitForBatchPiece", "CancelBatchPool"):
            fn = body(POOL, "subroutine", name)
            self.assertIn("if (WorkersAdvanced()) ticks = 0", fn, name)
            self.assertIn("if (ticks > MaxWaitTicks) then", fn, name)
            self.assertIn("call StopStalledPool(", fn, name)

    def test_progress_is_read_from_the_pieces_own_notes(self):
        fn = body(POOL, "function", "WorkersAdvanced")
        self.assertIn("call ReadPieceProgress(k)", fn)
        self.assertIn("WorkersAdvanced = s > LastProgSum", fn)

    def test_a_new_pool_starts_from_no_progress(self):
        self.assertIn("LastProgSum = -1", body(POOL, "subroutine", "StartPrepassBatches"))

    def test_a_stall_names_where_each_worker_got_to(self):
        fn = body(POOL, "subroutine", "StopStalledPool")
        self.assertIn("call DumpWorkerStdout(kind, k, 3)", fn)
        self.assertNotIn("pre-pass worker", fn)


if __name__ == "__main__":
    unittest.main()
