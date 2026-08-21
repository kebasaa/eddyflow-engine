"""The driver, and the one part of the chain no frozen R value reached.

This file used to carry the numerical pin as well: the engine's own
pre-whitening chain run over the two fixtures dyco ships, compared against
RFlux v3.2.0's frozen output for the same input. That is gone, because dyco is
no longer checked out beside this repository and the fixtures and the expected
values both lived in it, so the checks could only ever skip.

What was lost is worth stating plainly, because nothing here replaces it. It
was the only check in this directory that would notice the port computing the
wrong thing correctly - everything else reads source text - and it caught three
such defects the first time it ran:

  * the covariance divided by each lag's overlap count rather than by N,
    inflating it by 2.9% at lag 169 of 6000;
  * differencing that kept N samples by fabricating a zero, where R and numpy
    both return N-1, which moved the first AR coefficient in the sixth
    significant digit;
  * the AR-initialisation samples kept as zeros in the full-data CCF, where R
    drops them with na.action = na.omit, which at AR order 67 moved the
    correlation in the fourth.

All three fixes are still in the source, and test_pwb_borrowing_static.py says
where in it each one lives. What no longer exists is the check that would
notice any of them being undone. To restore it, put dyco back beside this
repository and recover this file from git history: it read the expected values
out of dyco's own test module rather than copying them here, so it needs both
the fixtures and that module.

What remains needs the driver but not dyco: the smoother against zoo's rule
written out independently, and the layout checks that keep m_pwb_core free of
engine state. Skipped when the driver has not been built
(`mingw32-make pwbref`).
"""

import subprocess
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
#> In obj/, not bin/. bin/ holds what a user runs; this is a test fixture.
DRIVER = ROOT / "obj" / "win" / "pwb_reference.exe"
if not DRIVER.exists():
    DRIVER = ROOT / "obj" / "linux" / "pwb_reference"

def run_smooth(width, n):
    """The engine's smoother over a fixed synthetic series."""
    proc = subprocess.run([str(DRIVER), "smooth", str(width), str(n)],
                          capture_output=True, text=True, check=True)
    out = {}
    for line in proc.stdout.splitlines():
        if line.startswith("y["):
            key, val = line.split("=", 1)
            out[int(key[2:-1])] = float(val)
    return out


def zoo_rolling_mean(v, width):
    """zoo's rollapply(align="center") plus the two-pass na.locf edge fill.

    An odd window is symmetric; an even one puts the extra sample after the
    centre. Written out rather than imported so the expectation does not come
    from the same place as the thing being checked.
    """
    n = len(v)
    lead, trail = (width - 1) // 2, width // 2
    out = [None] * n
    for i in range(lead, n - trail):
        out[i] = sum(v[i - lead:i + trail + 1]) / width
    first = next((i for i, x in enumerate(out) if x is not None), None)
    if first is None:
        return list(v)
    last = max(i for i, x in enumerate(out) if x is not None)
    for i in range(first):
        out[i] = out[first]
    for i in range(last + 1, n):
        out[i] = out[last]
    return out


def synthetic_series(n):
    """The ramp-with-a-spike the driver builds, mirrored here."""
    v = [float(i) for i in range(1, n + 1)]
    v[max(1, n // 2) - 1] += 100.0
    return v


@unittest.skipUnless(DRIVER.exists(),
                     "pwb_reference driver not built (mingw32-make pwbref)")
class TheSmootherFollowsZoo(unittest.TestCase):
    """The one part of the chain no frozen R value reaches.

    pww and cor_pww come from the UNSMOOTHED cross-correlation and cov_mcw
    from the raw cross-covariance, so the reference fixtures cannot see
    SmoothAndFill at all -- which is how a window that summed width+1 terms
    for an even width survived the pin. dyco verified the convention checked
    here against R for widths 4, 5 and 6; those are the three widths below.
    """

    def test_the_window_matches_zoo_for_both_parities(self):
        n = 9
        for width in (4, 5, 6):
            got = run_smooth(width, n)
            want = zoo_rolling_mean(synthetic_series(n), width)
            self.assertEqual(len(got), n, "width %d" % width)
            for i in range(1, n + 1):
                self.assertAlmostEqual(
                    got[i], want[i - 1], places=9,
                    msg="width %d, position %d: %r against zoo's %r"
                        % (width, i, got[i], want[i - 1]))

    def test_an_even_width_is_not_a_symmetric_window(self):
        #> Guard the guard. If SmoothAndFill went back to summing 2*(w/2)+1
        #> terms, the test above would catch it - but only because the two
        #> differ, so assert that they DO differ and this is not a tautology.
        n = 9
        width = 6
        symmetric = [sum(v) / len(v) for v in
                     [synthetic_series(n)[max(0, i - 3):i + 4] for i in range(n)]]
        zoo = zoo_rolling_mean(synthetic_series(n), width)
        self.assertNotEqual([round(x, 6) for x in symmetric],
                            [round(x, 6) for x in zoo])


class TheDriverProvesTheCoreIsIndependent(unittest.TestCase):
    """m_pwb_core must not reach for engine state.

    The driver links it against m_numeric_kinds and nothing else, so this is
    already enforced by the build - but only for whoever runs `make pwbref`.
    Stated here it is enforced for everyone.
    """

    def test_the_core_uses_only_the_kinds_module(self):
        raw = (ROOT / "src/src_common/m_pwb_core.f90").read_text(
            encoding="utf-8")
        #> Comments name the very things they promise are absent, so the scan
        #> has to see code only.
        source = (chr(10)).join(ln for ln in raw.splitlines()
                                if not ln.lstrip().startswith("!"))
        uses = [ln.strip() for ln in source.splitlines()
                if ln.strip().lower().startswith("use ")]
        self.assertEqual(uses, ["use m_numeric_kinds"])
        for forbidden in ("PWBSetup", "E2Col", "Metadata%", "ulog", "LogSay"):
            self.assertNotIn(forbidden, source)

    def test_the_driver_is_not_part_of_either_executable(self):
        #> src_common is a wildcard in the makefile, so a program unit there
        #> would be linked into both binaries and neither would have one main.
        self.assertTrue((ROOT / "src/src_tools/pwb_reference_main.f90").exists())
        self.assertFalse((ROOT / "src/src_common/pwb_reference_main.f90").exists())
        makefile = (ROOT / "prj/makefile").read_text(encoding="utf-8")
        self.assertIn("pwbref :", makefile)
        self.assertNotIn("pwbref", makefile.split("all : ")[1].splitlines()[0])

    def test_the_driver_does_not_land_in_bin(self):
        #> bin/ holds what a user runs. One unshipped test binary sitting
        #> beside the two real ones invites somebody to ship it.
        makefile = (ROOT / "prj/makefile").read_text(encoding="utf-8")
        recipe = makefile.split("pwbref :")[1].split(chr(10) + chr(10))[0]
        self.assertIn("$(OBJS_DIR)$(PWBREF_EXE)", recipe)
        self.assertNotIn("$(EXE_DIR)$(PWBREF_EXE)", recipe)


if __name__ == "__main__":
    unittest.main()
