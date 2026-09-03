"""PWB's arithmetic against the original R implementation.

Everything else in this directory reads source text. This runs the engine's
own pre-whitening chain over the two fixtures dyco ships and compares the
numbers with RFlux v3.2.0's output for the same input -- the AR order chosen
by AIC, the AR coefficients, the peak of the pre-whitened cross-correlation
and the raw cross-covariance, on both branches of the unit-root test.

It is the only check here that would notice the port computing the wrong
thing correctly. It caught three such defects when it was first run:

  * the covariance divided by each lag's overlap count rather than by N,
    inflating it by 2.9% at lag 169 of 6000;
  * differencing that kept N samples by fabricating a zero, where R and numpy
    both return N-1, which moved the first AR coefficient in the sixth
    significant digit;
  * the AR-initialisation samples kept as zeros in the full-data CCF, where R
    drops them with na.action = na.omit, which at AR order 67 moved the
    correlation in the fourth.

The expected values are READ OUT of dyco's own test module rather than copied
into this one, so the two cannot drift apart: if dyco regenerates its fixtures
it must update that table, and this check follows.

The bootstrap is deliberately not compared. It has no RNG stream in common
with R, which is exactly how dyco's own reference test treats it.

Skipped when dyco is not checked out beside this repository, or when the
driver has not been built (`mingw32-make pwbref`).
"""

import ast
import gzip
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DYCO = ROOT / "dyco-main"
FIXTURES = DYCO / "tests" / "data"
DYCO_TEST = DYCO / "tests" / "test_pwb_reference.py"
#> In obj/, not bin/. bin/ holds what a user runs; this is a test fixture.
DRIVER = ROOT / "obj" / "win" / "pwb_reference.exe"
if not DRIVER.exists():
    DRIVER = ROOT / "obj" / "linux" / "pwb_reference"

#: The fixtures were generated at 20 Hz with a 10 s search half-width; the R
#: run that produced the frozen numbers used the same two values.
HZ = 20
LAG_MAX_S = 10.0


def rflux_expectations():
    """The `_R` table out of dyco's test module, without importing it.

    Importing would pull in numpy, pandas and dyco itself. The table is a
    plain literal, so it can be read with ast alone.
    """
    tree = ast.parse(DYCO_TEST.read_text(encoding="utf-8"))
    for node in tree.body:
        if not isinstance(node, ast.Assign):
            continue
        if not any(getattr(t, "id", None) == "_R" for t in node.targets):
            continue
        out = {}
        for key, val in zip(node.value.keys, node.value.values):
            # dict(...) call, one per case
            out[ast.literal_eval(key)] = {
                kw.arg: ast.literal_eval(kw.value) for kw in val.keywords
            }
        return out
    raise AssertionError("no _R table found in %s" % DYCO_TEST)


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


def run_driver(fixture_gz):
    with tempfile.TemporaryDirectory() as tmp:
        csv = Path(tmp) / "fixture.csv"
        with gzip.open(fixture_gz, "rb") as fin, open(csv, "wb") as fout:
            shutil.copyfileobj(fin, fout)
        proc = subprocess.run(
            [str(DRIVER), str(csv), str(HZ), str(LAG_MAX_S)],
            capture_output=True, text=True, check=True)
    out = {}
    for line in proc.stdout.splitlines():
        if "=" not in line:
            continue
        k, v = line.split("=", 1)
        out[k.strip()] = v.strip()
    return out


@unittest.skipUnless(FIXTURES.exists(),
                     "dyco not checked out beside this repository")
@unittest.skipUnless(DRIVER.exists(),
                     "pwb_reference driver not built (mingw32-make pwbref)")
class PwbMatchesRflux(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.expected = rflux_expectations()
        cls.got = {
            case: run_driver(FIXTURES / ("pwb_reference_%s.csv.gz" % case))
            for case in ("stationary", "differencing")
        }

    def assertClose(self, got, want, sigfigs, what):
        """Agreement to a number of significant digits, not to an absolute."""
        tol = abs(want) * 10.0 ** (-sigfigs)
        self.assertLessEqual(
            abs(float(got) - want), tol,
            "%s: %r against R's %r, wanted %d significant digits"
            % (what, got, want, sigfigs))

    def test_the_table_was_found_at_all(self):
        #> Guard the guard: a parse that quietly returned nothing would make
        #> every assertion below vacuously true.
        self.assertEqual(set(self.expected), {"stationary", "differencing"})
        for case in self.got.values():
            self.assertIn("cov_mcw", case)

    def test_the_unit_root_decision_matches(self):
        #> Which branch runs is the first thing that has to agree; everything
        #> downstream is computed on different series otherwise.
        for case, want in self.expected.items():
            self.assertEqual(self.got[case]["differenced"] == "T",
                             want["differenced"], case)

    def test_the_aic_order_search_matches_exactly(self):
        #> An integer, so it agrees exactly or the search is wrong. The
        #> differencing fixture selects 67/56/56, which is deep enough into
        #> the Levinson-Durbin recursion to be a real test of it.
        for case, want in self.expected.items():
            got = self.got[case]
            self.assertEqual(
                (int(got["ar_order_scalar"]), int(got["ar_order_w"]),
                 int(got["ar_order_t"])),
                tuple(want["ar_orders"]), case)

    def test_the_ar_coefficients_match_to_ten_significant_digits(self):
        for case, want in self.expected.items():
            got = self.got[case]
            for key, expect in zip(("phi1_scalar", "phi1_w", "phi1_t"),
                                   want["phi1"]):
                self.assertClose(got[key], expect, 10, "%s %s" % (case, key))

    def test_the_prewhitened_ccf_peak_matches(self):
        #> tl_pww is a lag in records and agrees exactly; cor_pww is the CCF
        #> there. This pair is what the AR-initialisation trim decides.
        for case, want in self.expected.items():
            got = self.got[case]
            self.assertEqual(int(got["pww"]), want["pww"], case)
            self.assertClose(got["cor_pww"], want["cor_pww"], 10,
                             "%s cor_pww" % case)

    def test_the_raw_cross_covariance_matches(self):
        #> The quantity that has to be read off the UNDIFFERENCED series, with
        #> R's biased divisor. Both mistakes are silent in the detected lag
        #> and loud here.
        for case, want in self.expected.items():
            got = self.got[case]
            self.assertEqual(int(got["mcw"]), want["mcw"], case)
            self.assertClose(got["cov_mcw"], want["cov_mcw"], 10,
                             "%s cov_mcw" % case)


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
