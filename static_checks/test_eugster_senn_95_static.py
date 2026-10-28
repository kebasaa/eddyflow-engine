"""The Eugster and Senn (1995) high-frequency correction, hf_meth = 7.

Eugster, W. and Senn, W. (1995): A cospectral correction model for measurement
of turbulent NO2 flux. Boundary-Layer Meteorol. 74, 321-340.

Four things hold the method together, and each fails silently if it drifts:

1. The constants are the paper's. They are read out of the Fortran here and
   checked against the paper's own Table III, so a transposed digit fails
   this test instead of moving every flux by a percent.

2. The stable branch is Eq. 33 with a common root cancelled. Checked against
   the equation as printed at points where the printed form is defined.

3. It is analytic. RP runs it, and FCC is not forced, so it has to be in
   FccFollowsProject's list beside Moncrieff and Massman.

4. Its discard reaches the flux. The 0.4 threshold sets a mask that both
   Level-3 writers read, so a gas the paper would discard is not reported.
"""

import math
import re
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
COMMON = ROOT / "src" / "src_common"
ROUTINE = COMMON / "bpcf_eugster_senn_95.f90"
HUB = COMMON / "bpcf_bandpass_spectral_corrections.f90"
DECODER = COMMON / "write_processing_project_variables.f90"
FOLLOWS = COMMON / "reported_random_errors.f90"
TYPEDEF = COMMON / "m_typedef.f90"
GLOBALS = COMMON / "m_common_global_var.f90"
TAGS = ROOT / "prj" / "gen_project_tags.py"
FLUXNET = ROOT / "src" / "src_rp" / "write_out_fluxnet.f90"
LEVEL3 = (ROOT / "src" / "src_rp" / "fluxes23_rp.f90",
          ROOT / "src" / "src_fcc" / "fluxes23.f90")

#: Eugster and Senn (1995), Table III: correction factors 1/xi for unstable
#: and neutral stratification. Rows are z/u [s], columns are L [s].
TABLE_L = (0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.75, 1.0)
TABLE_III = {
    0.1: (2.462, 2.675, 2.869, 3.045, 3.206, 3.354, 3.939, 4.350),
    0.2: (1.830, 1.969, 2.101, 2.227, 2.348, 2.462, 2.959, 3.354),
    0.5: (1.373, 1.438, 1.501, 1.563, 1.624, 1.684, 1.969, 2.227),
    1.0: (1.204, 1.240, 1.274, 1.307, 1.341, 1.373, 1.532, 1.684),
    2.0: (1.110, 1.130, 1.150, 1.168, 1.187, 1.204, 1.291, 1.373),
    5.0: (1.041, 1.052, 1.062, 1.071, 1.080, 1.089, 1.130, 1.168),
    10.0: (1.010, 1.017, 1.023, 1.030, 1.035, 1.041, 1.066, 1.089),
}


def read(path):
    return path.read_text(encoding="utf-8", errors="replace")


SRC = read(ROUTINE)
FRACTION = SRC[SRC.index("function EugsterSenn95Fraction("):]

NUM = r"([0-9.]+)d0"
UNSTABLE = re.compile(
    r"xi = " + NUM + r" \* datan\(" + NUM
    + r" \* dlog\(z_over_u / tau_L \+ 1d0\) \+ " + NUM + r"\)")
N0 = re.compile(
    r"n0 = " + NUM + r" \* \(1d0 \+ " + NUM + r" \* zeta\)\*\*" + NUM)


def constants():
    a, b, c = (float(v) for v in UNSTABLE.search(FRACTION).groups())
    n0a, n0b, n0e = (float(v) for v in N0.search(FRACTION).groups())
    return a, b, c, n0a, n0b, n0e


def unstable(z_over_u, tau_l):
    a, b, c = constants()[:3]
    return min(a * math.atan(b * math.log(z_over_u / tau_l + 1.0) + c), 1.0)


def stable(z_over_u, zl, tau_l):
    n0a, n0b, n0e = constants()[3:]
    n0 = n0a * (1.0 + n0b * min(zl, 2.0)) ** n0e
    return z_over_u / (z_over_u + 4.0 * math.pi / math.sqrt(6.0) * tau_l * n0)


def stable_as_printed(z_over_u, zl, tau_l):
    """Eq. 33 exactly as the paper prints it."""
    n0 = 0.23 * (1.0 + 6.4 * zl) ** 0.75
    x, s6 = z_over_u, math.sqrt(6.0)
    return (3.0 / s6 * x * (s6 * x - 4.0 * math.pi * tau_l * n0)
            / (3.0 * x * x - 8.0 * math.pi ** 2 * tau_l ** 2 * n0 ** 2))


class TheConstants(unittest.TestCase):

    def test_they_are_the_papers(self):
        self.assertEqual(constants(), (0.725, 1.24, 0.21, 0.23, 6.4, 0.75))

    def test_they_reproduce_table_iii(self):
        for z_over_u, row in TABLE_III.items():
            for tau_l, factor in zip(TABLE_L, row):
                with self.subTest(z_over_u=z_over_u, L=tau_l):
                    self.assertAlmostEqual(1.0 / unstable(z_over_u, tau_l),
                                           factor, delta=0.0015)

    def test_the_unstable_fit_is_capped_at_one(self):
        #> The fit tends to 0.725 pi/2 = 1.14 as L goes to zero.
        self.assertIn("xi = min(xi, 1d0)", FRACTION)
        self.assertEqual(unstable(10.0, 1e-6), 1.0)

    def test_stability_is_held_at_two(self):
        self.assertIn("zeta = min(zL, 2d0)", FRACTION)

    def test_an_undamped_sensor_measures_everything(self):
        self.assertRegex(FRACTION, r"if \(tau_L <= 0d0\) then\s+xi = 1d0")


class TheStableBranch(unittest.TestCase):

    def test_it_is_eq_33_with_the_common_root_cancelled(self):
        for zl in (0.1, 0.5, 1.0, 2.0):
            for tau_l in (0.1, 0.35, 1.0):
                for z_over_u in (0.1, 0.3, 1.0, 3.0, 10.0):
                    root = 4.0 * math.pi * tau_l * 0.23 * (1 + 6.4 * zl) ** 0.75 / math.sqrt(6)
                    if abs(z_over_u - root) < 1e-3:
                        continue
                    with self.subTest(zl=zl, L=tau_l, z_over_u=z_over_u):
                        self.assertAlmostEqual(
                            stable(z_over_u, zl, tau_l),
                            stable_as_printed(z_over_u, zl, tau_l), places=12)

    def test_it_is_written_in_the_cancelled_form(self):
        self.assertIn(
            "xi = z_over_u / (z_over_u + 4d0 * p / dsqrt(6d0) * tau_L * n0)",
            FRACTION)


class TheWiring(unittest.TestCase):

    def test_seven_decodes_to_the_method(self):
        src = read(DECODER)
        self.assertRegex(src, r"case \('7'\)[\s\S]{0,200}?"
                              r"hf_meth = 'eugster_senn_95'[\s\S]{0,80}?"
                              r"hf_meth_in_situ = \.false\.")

    def test_the_hub_dispatches_it(self):
        src = read(HUB)
        self.assertRegex(src, r"case\('eugster_senn_95'\)\s+[\s\S]{0,120}?"
                              r"call BPCF_EugsterSenn95\(")

    def test_it_does_not_force_fcc(self):
        self.assertIn("EddyFlowProj%hf_meth /= 'eugster_senn_95'", read(FOLLOWS))

    def test_the_fluxnet_column_has_a_code_for_it(self):
        src = read(FLUXNET)
        block = src[src.index("!> Spectral correction method"):]
        block = block[:block.index("end select")]
        self.assertRegex(block, r"case\('eugster_senn_95'\)\s+"
                                r"call AddIntDatumToDataline\(7,")
        self.assertIn("case default", block)


class TheTimeConstant(unittest.TestCase):

    def test_es_tau_lives_in_the_tail_block(self):
        #> Not in GAS_NUMERIC: growing that moves every later record origin.
        gen = read(TAGS)
        head = gen[gen.index("GAS_NUMERIC = ["):gen.index("GAS_TEXT")]
        self.assertNotIn("es_tau", head)
        self.assertIn('GAS_TAIL_NUMERIC = ["es_tau"]', gen)

    def test_the_tail_block_comes_after_every_other_record_block(self):
        tags = read(GLOBALS)
        origin = lambda name: int(re.search(
            r"integer, parameter :: %s\s*=\s*(\d+)" % name, tags).group(1))
        self.assertGreater(origin("gasTailOriginN"), origin("agcRecOriginN"))
        self.assertIn("'gas_1_es_tau'", tags)

    def test_it_is_read_into_the_gas_record(self):
        self.assertIn("real(kind = dbl) :: es_tau", read(TYPEDEF))
        src = read(DECODER)
        self.assertIn("b = gasTailOriginN + (i - 1) * gasTailLeapN", src)
        self.assertIn("EddyFlowProj%gas(i)%es_tau = &", src)

    def test_the_analysers_tau_is_the_fallback(self):
        self.assertRegex(SRC, r"tau_L = EddyFlowProj%gas\(gas - firstGas \+ 1\)%es_tau\s+"
                              r"if \(tau_L < 0d0\) tau_L = LocInstr\(gas\)%tau")


class TheDiscard(unittest.TestCase):

    def test_below_forty_percent_the_gas_is_marked(self):
        self.assertIn("if (xi < 0.4d0) BPCF_discard(gas) = .true.", SRC)

    def test_the_mask_is_reset_every_period(self):
        src = read(HUB)
        reset = src.index("BPCF_discard = .false.")
        self.assertLess(reset, src.index("select case(trim(adjustl(actual_hf_method)))"))
        self.assertIn("logical :: BPCF_discard(GHGNumVar)", read(GLOBALS))

    def test_both_level_3_writers_honour_it(self):
        for path in LEVEL3:
            with self.subTest(path=path.name):
                self.assertRegex(read(path),
                                 r"if \(\.not\. BPCF_discard\(msl\)\) cycle\s+"
                                 r"Flux3%gas\(msl\) = error")


if __name__ == "__main__":
    unittest.main()
