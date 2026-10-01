"""The engine asks Windows not to run it as background work.

Windows classes a windowless program as background work, and on a hybrid
processor keeps background work on the efficiency cores. The engine is always
windowless - a child of the interface with its output captured, or a pre-pass
worker started through a script - and on the Yatir run every engine process
sat at 100 % on the four low-power cores of a Core Ultra 7 268V while the four
performance cores idled at 3-8 %.

``RequestFullSpeed`` turns execution-speed throttling off for the process
(ControlMask = EXECUTION_SPEED, StateMask = 0), the documented opt-out. It
changes where the program runs, never what it computes: every regression
fixture was byte-identical before and after. Read back from outside on a
running engine, the state is ControlMask=1 StateMask=0; the build before it
reads 0/0.

Part of the EddyFlow engine's static checks.
"""

from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[1]
WIN = ROOT / "src" / "src_os" / "win" / "m_process_os.f90"
POSIX = ROOT / "src" / "src_os" / "posix" / "m_process_os.f90"
INIT = ROOT / "src" / "src_common" / "init_env.f90"


def code(path):
    return "\n".join(ln for ln in path.read_text(encoding="utf-8", errors="replace").splitlines()
                     if not ln.lstrip().startswith("!"))


class TheProcessAsksForFullSpeed(unittest.TestCase):

    def test_execution_speed_throttling_is_turned_off(self):
        w = code(WIN)
        body = w[w.index("subroutine RequestFullSpeed"):w.index("end subroutine RequestFullSpeed")]
        self.assertIn("state%control_mask = THROTTLING_EXECUTION_SPEED", body)
        self.assertIn("state%state_mask = 0_c_int32_t", body)
        self.assertIn("w_SetProcessInformation(w_GetCurrentProcess(), PROCESS_POWER_THROTTLING", body)

    def test_the_constants_are_windows_own(self):
        w = code(WIN)
        #> ProcessPowerThrottling in PROCESS_INFORMATION_CLASS, and
        #> PROCESS_POWER_THROTTLING_CURRENT_VERSION / _EXECUTION_SPEED.
        self.assertIn("PROCESS_POWER_THROTTLING = 4_c_int", w)
        self.assertIn("THROTTLING_VERSION = 1_c_int32_t", w)
        self.assertIn("THROTTLING_EXECUTION_SPEED = 1_c_int32_t", w)
        self.assertIn("bind(C, name = 'SetProcessInformation')", w)

    def test_a_failure_is_ignored(self):
        """Scheduling, never results: an old Windows just runs as before."""
        w = code(WIN)
        body = w[w.index("subroutine RequestFullSpeed"):w.index("end subroutine RequestFullSpeed")]
        self.assertNotIn("error stop", body)
        self.assertNotIn("stop ", body)

    def test_every_platform_has_it(self):
        self.assertIn("subroutine RequestFullSpeed()", code(POSIX))
        self.assertIn("RequestFullSpeed", code(POSIX).split("public ::")[1].splitlines()[0])

    def test_every_process_asks_once_at_start_up(self):
        """Workers run the same InitEnv, so each asks for itself."""
        init = code(INIT)
        self.assertEqual(init.count("call RequestFullSpeed()"), 1)
        self.assertLess(init.index("end do arg_loop"), init.index("call RequestFullSpeed()"))


if __name__ == "__main__":
    unittest.main()
