"""The dyco comparison drives the engine's own PWB chain, on given resamples.

tests/pwb_dyco/compare.py compares the engine's PWB detector with dyco's,
number for number, by handing both the same bootstrap block starts. That only
means something if the driver it runs, src/src_tools/pwb_compare_main.f90,
calls the routines PwbDetectGas calls - so the per-combination bootstrap, its
summary and the choice of combination live in m_pwb_core as functions of their
arguments, the block starts among them, and both the engine and the driver call
them. This pins that wiring; the comparison itself needs scipy and dyco and is
run by hand (tests/pwb_dyco/README.md).

Part of the EddyFlow engine's static checks.
"""

from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]


def read(path):
    return (ROOT / path).read_text(encoding="utf-8")


def code(path):
    return (chr(10)).join(ln for ln in read(path).splitlines()
                          if not ln.lstrip().startswith("!"))


CORE = code("src/src_common/m_pwb_core.f90")
PWB = code("src/src_rp/pwb_timelag_handle.f90")
DRIVER = code("src/src_tools/pwb_compare_main.f90")
MAKE = read("prj/Makefile")
COMPARE = read("tests/pwb_dyco/compare.py")


class TheBootstrapTakesItsResamples(unittest.TestCase):

    def test_the_core_routines_are_public(self):
        for name in ("PwbDrawBlockStarts", "PwbBootstrapCombination",
                     "PwbSummariseBootstrap", "PwbBestCombination",
                     "FillMissingLinear"):
            self.assertIn(name, CORE[:CORE.index("contains")], name)

    def test_the_starts_are_an_argument(self):
        self.assertIn("integer, intent(in) :: starts(nblocks, nboot)", CORE)

    def test_the_engine_draws_them_and_hands_them_over(self):
        run = PWB[PWB.index("subroutine RunPwbCombination("):
                  PWB.index("end subroutine RunPwbCombination")]
        self.assertIn("call PwbDrawBlockStarts(state, n, block_len, nblocks, nboot, starts)", run)
        self.assertIn("call PwbBootstrapCombination(x, y, n, min_rl, max_rl, eval_lo, eval_hi,", run)
        self.assertIn("call PwbSummariseBootstrap(", run)
        self.assertNotIn("RandBelow", run)

    def test_the_choice_of_combination_is_the_cores(self):
        self.assertIn("best = PwbBestCombination(candidate(:)%ccf_at_mode, ok, 4)", PWB)


class TheDriverRunsTheSameChain(unittest.TestCase):

    def test_it_calls_what_pwb_detect_gas_calls(self):
        for call in ("call FillMissingLinear(", "call PwbPreWhiten(",
                     "call PwbBootstrapCombination(", "call PwbSummariseBootstrap(",
                     "PwbBestCombination("):
            self.assertIn(call, DRIVER, call)

    def test_it_uses_the_engines_evaluated_range(self):
        self.assertIn("margin = max(trail, nint(2d0 * hz))", DRIVER)
        self.assertIn("margin = max(trail, nint(2d0 * Metadata%ac_freq))", PWB)

    def test_it_links_against_the_core_alone(self):
        self.assertIn("pwbcmp : m_numeric_kinds.o m_pwb_core.o", MAKE)
        self.assertNotIn("use m_rp_global_var", DRIVER)


class TheComparisonFailsOnWhatShouldMatch(unittest.TestCase):

    def test_a_mode_mismatch_fails(self):
        """The mode is dyco's estimator now; only edge-pinned choices may differ."""
        self.assertIn('fail("%s mode %d vs dyco %d"', COMPARE)
        self.assertIn("return np.zeros(size)", COMPARE)

    def test_replicate_lags_must_be_identical(self):
        self.assertIn("if elags != dlags:", COMPARE)


class TheModeIsDycos(unittest.TestCase):

    def test_scott_bandwidth_on_a_512_point_grid(self):
        mode = CORE[CORE.index("integer function MapLagEstimate"):
                    CORE.index("end function MapLagEstimate")]
        self.assertIn("integer, parameter :: ngrid = 512", mode)
        self.assertIn("bw2 = var_s * (dble(n)**(-0.2d0))**2", mode)
        self.assertNotIn("1.06d0", mode)

    def test_rounded_half_to_even(self):
        mode = CORE[CORE.index("integer function MapLagEstimate"):
                    CORE.index("end function MapLagEstimate")]
        self.assertIn("elseif (mod(int(f), 2) == 0) then", mode)


if __name__ == "__main__":
    unittest.main()
