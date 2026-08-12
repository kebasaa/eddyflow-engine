"""One humidity in the air, or one per hygrometer - but not one of each.

When a biomet RH value is available, `FluxParams` takes the site scalars from
it and then recomputes the reported mole fraction, mixing ratio and molar
density of the hygrometer. It used to do that for the *primary* hygrometer
alone, so a site with two reported biomet for one of them and the instrument's
own measurement for the other - and which was which followed the primary
designation, a naming choice with no business deciding whose numbers are real.

Seen on CH-LAE across two runs that differed only in gas-record order: a mixing
ratio of 19.9081 followed the primary slot and matched neither instrument. The
LI-7200 measured 17.1089 and the MIRO 16.354; 19.9081 was the biomet value,
attached to whichever record happened to be first.

Now every water record takes it. The consequence reaches further than the
reported columns, and deliberately so: `RHO%w_at` is built from `Stats%chi`
later in this same routine, so `sigma_at`, `Q_at`, `rho_a_at`, `RhoCp_at` and
`RH_at` all follow, and a gas is WPL-corrected with the humidity of the air
rather than with whichever hygrometer it was paired to. That propagation is the
reason the change is one loop rather than a rewrite, and it is why these checks
pin the derivation as well as the loop - break the link and the loop silently
stops mattering.

What did *not* change: with no biomet RH the site scalars still come from the
primary hygrometer's own channel. That arm is untouched, and a project without
a biomet RH column sees nothing of this.
"""

from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]

PARAMS = "src/src_rp/flux_params.f90"


def read(path):
    return (ROOT / path).read_text(encoding="utf-8")


def code(path):
    """Source with comment-only lines dropped.

    The routine explains the retired single-slot behaviour in prose right where
    it stood, so a naive search matches the explanation rather than live code.
    """
    return "\n".join(ln for ln in read(path).splitlines()
                     if not ln.lstrip().startswith("!"))


def biomet_branch():
    """The `biomet RH is available` arm, up to the arm that follows it."""
    src = code(PARAMS)
    start = src.index("if (biomet%val(bRH) > 0d0")
    return src[start: src.index("elseif (wsl >= firstGas) then", start)]


class EveryHygrometerTakesIt(unittest.TestCase):
    def setUp(self):
        self.block = biomet_branch()

    def test_the_concentrations_are_written_in_a_loop(self):
        self.assertIn("do msl = firstGas, lastGas", self.block)
        self.assertIn("if (.not. GasSlotIsWater(msl)) cycle", self.block)

    def test_no_concentration_is_written_to_the_primary_slot_alone(self):
        """`wsl` is the primary. Writing chi/r/d through it is the defect."""
        for field in ("Stats%chi(wsl)", "Stats%r(wsl)", "Stats%d(wsl)"):
            self.assertNotIn(field, self.block,
                             "%s makes the primary designation decide whose "
                             "humidity is measured and whose is biomet" % field)

    def test_all_three_quantities_are_written_per_slot(self):
        for field in ("Stats%chi(msl)", "Stats%r(msl)", "Stats%d(msl)"):
            self.assertIn(field, self.block)

    def test_the_cell_terms_are_the_slot_own(self):
        """A closed-path hygrometer's molar density goes through its own cell,
        and two analysers do not share one."""
        self.assertIn("E2Col(msl)%instr%path_type == 'closed'", self.block)
        self.assertIn("E2Col(msl)%Va", self.block)
        self.assertNotIn("E2Col(wsl)%Va", self.block)

    def test_an_absent_record_is_skipped(self):
        self.assertIn("if (.not. E2Col(msl)%present) cycle", self.block)


class TheDerivationIsNotBypassed(unittest.TestCase):
    """The loop only matters because everything downstream is built from
    `Stats%chi`. Compute RHO%w_at from anything else and the per-gas moisture
    terms go back to describing whichever instrument they came from, while this
    loop goes on looking correct."""

    def test_the_per_water_density_still_comes_from_chi(self):
        self.assertIn(
            "RHO%w_at(msl) = (Stats%chi(msl) / Ambient%Va) * MW_H2O * 1d-3",
            code(PARAMS))

    def test_it_is_computed_after_the_biomet_branch(self):
        src = code(PARAMS)
        self.assertLess(src.index("if (biomet%val(bRH) > 0d0"),
                        src.index("RHO%w_at(msl) ="),
                        "the per-water densities must be built after the "
                        "override, or they carry the pre-override values")

    def test_the_rest_of_the_regime_is_built_from_those_densities(self):
        src = code(PARAMS)
        for field in ("Ambient%e_at(msl) = RHO%w_at(msl)",
                      "Ambient%rho_a_at(msl) = Ambient%rho_d_at(msl) + RHO%w_at(msl)"):
            self.assertIn(field, src)


class WithoutBiometNothingChanges(unittest.TestCase):
    """A project with no biomet RH column must be untouched by any of this."""

    def test_the_fallback_arm_still_reads_the_primary(self):
        src = code(PARAMS)
        start = src.index("elseif (wsl >= firstGas) then")
        block = src[start: start + 900]
        self.assertIn("RHO%w = (Stats%chi(wsl) / Ambient%Va) * MW_H2O * 1d-3",
                      block,
                      "with no biomet humidity the site scalars come from the "
                      "designated hygrometer's own channel, as they always have")

    def test_the_override_is_still_gated_on_a_usable_value(self):
        self.assertIn("if (biomet%val(bRH) > 0d0 .and. biomet%val(bRH) < RHmax) then",
                      code(PARAMS))


if __name__ == "__main__":
    unittest.main()
