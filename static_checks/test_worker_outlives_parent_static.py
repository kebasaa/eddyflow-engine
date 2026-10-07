"""A pre-pass worker never outlives the process that started it.

A worker is a separate process, launched through a shell script that exits as
soon as it has started it, so nothing ties the worker's life to its parent's.
Killing the parent - the interface's Stop button, Task Manager, an error stop
while it waits for the workers - used to leave every worker running to the end
of its slice, writing records into a directory nobody would read again.

Found on the Yatir run: stopped from the interface at about 13:10, and at 13:12
the parent was gone while six workers and their six cmd.exe launchers were
still computing. One of them finished at 13:19 and wrote a 154 KB dump that no
process would ever merge.

So the parent hands each worker its process ID (``--batch-parent``), the worker
opens a handle to it at start-up, and checks it at the top of every period. If
the parent is gone, the worker tidies up as it would on finishing, writes
nothing, and exits with 3. Measured: parent killed, all five workers gone two
seconds later, every ``.rc`` reading 3, no dump; a worker started with a parent
ID that does not exist stops at once, before creating anything.

The interface side - a job object, so that Stop and Pause reach every process
of the run - is in the GUI repository. This is the engine side, which also
covers every way a parent can end that the interface never sees.

Part of the EddyFlow engine's static checks.
"""

from pathlib import Path
import re
import unittest


ROOT = Path(__file__).resolve().parents[1]

WIN = ROOT / "src" / "src_os" / "win" / "m_process_os.f90"
POSIX = ROOT / "src" / "src_os" / "posix" / "m_process_os.f90"
INIT = ROOT / "src" / "src_common" / "init_env.f90"
PARALLEL = ROOT / "src" / "src_rp" / "prepass_parallel.f90"
# The worker pool itself moved to src_common/batch_pool.f90, shared with FCC;
# what is checked here is read from both.
POOL = ROOT / "src" / "src_common" / "batch_pool.f90"
MAIN = ROOT / "src" / "src_rp" / "eddyflow-rp_main.f90"
MAKEFILE = ROOT / "prj" / "Makefile"


def read(path):
    return path.read_text(encoding="utf-8", errors="replace")


def code(path):
    return (chr(10)).join(ln for ln in read(path).splitlines()
                          if not ln.lstrip().startswith("!"))


def publics(path):
    names = set()
    for m in re.finditer(r"(?im)^\s*public\s*::\s*(.+)$", code(path)):
        names |= {n.strip().lower() for n in m.group(1).split(",")}
    return names


class TheTwoPlatformFilesAgree(unittest.TestCase):
    """One is compiled per platform and gen_makefile_deps.py reads only the
    Windows one, so they must declare the same module and the same names."""

    def test_same_module(self):
        for path in (WIN, POSIX):
            self.assertRegex(code(path), r"(?im)^\s*module\s+m_process_os\s*$")

    def test_same_public_names(self):
        self.assertEqual(publics(WIN), publics(POSIX))
        self.assertTrue({"processselfid", "watchparent", "parentgone"} <= publics(WIN))

    def test_neither_uses_another_engine_module(self):
        """Their Makefile rule is the source alone; a `use` of an engine module
        in either would need a dependency the generator never sees for posix."""
        for path in (WIN, POSIX):
            uses = re.findall(r"(?im)^\s*use\b[\s,]*(?:intrinsic\s*::\s*)?(\w+)", code(path))
            self.assertEqual({u.lower() for u in uses}, {"iso_c_binding"}, path.name)

    def test_the_makefile_picks_one_by_platform(self):
        mk = read(MAKEFILE)
        self.assertIn("OS_SRC_DIR_f90d1 = ../src/src_os/$(OS_FAMILY)/", mk)
        self.assertIn("OS_FAMILY = win", mk)
        self.assertIn("OS_FAMILY = posix", mk)


class TheParentIsNamedAndWatched(unittest.TestCase):

    def test_the_parent_passes_its_id(self):
        self.assertIn("// ' --batch-parent ' // trim(parentId)", (code(PARALLEL) + code(POOL)))
        self.assertIn("write(parentId, '(i0)') ProcessSelfId()", (code(PARALLEL) + code(POOL)))

    def test_the_switch_is_read_and_takes_a_value(self):
        init = code(INIT)
        self.assertIn("case('--batch-parent')", init)
        #> The argument loop drops every switch after one it does not know
        #> takes a value - which is how --batch was once silently lost and a
        #> worker spawned workers of its own.
        self.assertIn("'--batch-parent')", init[init.index("function SwitchTakesValue"):])

    def test_it_is_watched_before_anything_is_created(self):
        init = code(INIT)
        watch = init.index("WatchParent(BatchParentPid)")
        self.assertLess(init.index("end do arg_loop"), watch)
        self.assertLess(watch, init.index("call SetOSEnvironment()"))

    def test_both_pre_pass_loops_check_every_period(self):
        main = code(MAIN)
        for label in ("to_periods_loop: do", "pf_periods_loop: do"):
            i = main.index(label)
            self.assertIn("call StopIfParentGone()", main[i:i + 120], label)


class AnOrphanLeavesNothingBehind(unittest.TestCase):

    def stop(self):
        par = (code(PARALLEL) + code(POOL))
        return par[par.index("subroutine StopIfParentGone"):
                   par.index("end subroutine StopIfParentGone")]

    def test_it_writes_no_records(self):
        self.assertNotIn("Dump", self.stop())

    def test_it_tidies_as_a_finishing_worker_does(self):
        self.assertIn("call FinishBatchWorker()", self.stop())

    def test_it_exits_distinctly(self):
        """3, so an .rc file read later tells an orphan from a failure."""
        self.assertIn("stop 3", self.stop())

    def test_only_a_worker_checks(self):
        self.assertIn("if (BatchIndex <= 0) return", self.stop())


if __name__ == "__main__":
    unittest.main()
