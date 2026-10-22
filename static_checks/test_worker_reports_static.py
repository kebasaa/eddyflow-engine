"""A run never waits for a worker that can no longer report, and two sessions'
runs keep out of each other's temporary folders.

A Yatir run waited a day and gave up: every closed file in its temporary
folder had been removed while its workers ran - by another session's Stop,
which wiped the whole tmp folder of the version. cmd reads a worker's script
as it runs it, so with the script gone the worker finished and its return code
was never written.

Each worker now writes its process id beside its return code; the pool looks
every ten seconds, and a worker that has exited without a code - or has
neither a process id nor its script a minute after it was started - stops the
run with a message saying what happened. Each run's temporary folder names its
owner, so the interface removes only those of runs that have ended.
"""

import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def read(rel):
    return (ROOT / rel).read_text(encoding="utf-8", errors="replace").replace("\r\n", "\n")


POOL = read("src/src_common/batch_pool.f90")
ENV = read("src/src_common/init_env.f90")


def body(text, kind, name):
    i = text.index("%s %s(" % (kind, name))
    return text[i:text.index("end %s %s" % (kind, name), i)]


class EachRunNamesItsOwner(unittest.TestCase):

    def test_the_temporary_folder_carries_the_process_id(self):
        self.assertIn("write(pidPadding, '(a,i0)') '_p', ProcessSelfId()", ENV)
        self.assertIn("tmpDirPadding = trim(tmpDirPadding) // trim(pidPadding)", ENV)
        self.assertIn("character(48) :: tmpDirPadding", ENV)

    def test_a_desktop_run_writes_owner_pid_into_it(self):
        i = ENV.index("if (EddyFlowProj%run_env == 'desktop') then\n        open(newunit = upid, file = trim(TmpDir) // 'owner.pid'")
        self.assertGreater(i, ENV.index("make_dir = CreateDir('\"' // trim(TmpDir) // '\"')"))

    def test_a_worker_writes_its_process_id_beside_its_return_code(self):
        self.assertIn("if (BatchIndex > 0 .and. len_trim(BatchTmpDir) > 0) then", ENV)
        self.assertIn("// trim(BatchKind) // trim(batchPadding) // '.pid'", ENV)
        self.assertIn("PidPath = trim(TmpDir) // 'batch_' // trim(tag) // '.pid'",
                      body(POOL, "function", "PidPath"))


class AWorkerThatCannotReportStopsTheRun(unittest.TestCase):

    C = body(POOL, "subroutine", "CheckWorkersAlive")

    def test_every_poll_looks_every_ten_seconds(self):
        self.assertIn("call CheckWorkersAlive(quiet)", body(POOL, "subroutine", "PollPrepassPool"))
        self.assertIn("if (now - LastLivenessCheck < 10 * rate) return", self.C)

    def test_an_exited_worker_gets_ten_seconds_to_show_its_code(self):
        self.assertIn("if (.not. ProcessAlive(pid)) SuspectAt(k) = now", self.C)
        self.assertIn("if (now - SuspectAt(k) >= 10 * rate) &", self.C)

    def test_no_process_id_and_no_script_after_a_minute(self):
        self.assertIn("if (now - LaunchedAt(k) < 60 * rate) cycle", self.C)
        self.assertIn("if (.not. ex) call WorkerFilesGone(PoolKind, k, quiet)", self.C)

    def test_a_missing_script_is_not_launched(self):
        launch = body(POOL, "subroutine", "LaunchChunks")
        i = launch.index("if (.not. ex) call WorkerFilesGone(kind, k, .false.)")
        self.assertLess(i, launch.index("open(newunit = u, file = trim(masterPath)"))

    def test_stale_process_ids_are_cleared_with_stale_return_codes(self):
        start = body(POOL, "subroutine", "StartPrepassBatches")
        self.assertIn("rcPath = PidPath(kind, k)", start)

    def test_the_message_says_what_happened_and_shows_the_end_of_its_output(self):
        gone = body(POOL, "subroutine", "WorkerFilesGone")
        self.assertIn("were removed while the run was using them.", gone)
        self.assertIn("call DumpWorkerStdout(kind, k, 15)", gone)
        self.assertIn("error stop 'A parallel worker ended without reporting back.'", gone)


if __name__ == "__main__":
    unittest.main()
