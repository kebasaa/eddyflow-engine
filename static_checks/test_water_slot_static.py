"""Nothing may read a slot number and mean "water".

Slots are assigned by record order, so the h2o constant identifies record two
and nothing else - and with one hygrometer per analyser there is more than one
water record, so even "record two" does not name a unique measurement.

The fixture base_h2o_late is base_n_gas with records 2 and 5 swapped. Every
column it produces must equal base_n_gas's, because the same physical columns
feed them. Getting there took six separate fixes, each pinned below; on the
build before them all 25 water-derived scalars were wrong while looking
entirely plausible - RH 70.3% -> 3.6%, dewpoint 11.5C -> -25.6C, Bowen -12 ->
8976, latent heat sign-flipped.
"""

from pathlib import Path
import re
import unittest


ROOT = Path(__file__).resolve().parents[1]


def read(path):
    return (ROOT / path).read_text(encoding="utf-8")


def code(path):
    """Source with comment lines dropped, so a note *about* a removed
    construct does not read as the construct itself."""
    return "\n".join(ln for ln in read(path).splitlines()
                     if not ln.lstrip().startswith("!"))


#> Every file that computes a water-derived quantity.
FLUX_FILES = (
    "src/src_rp/flux_params.f90",
    "src/src_rp/fluxes0_rp.f90",
    "src/src_rp/fluxes1_rp.f90",
    "src/src_rp/fluxes23_rp.f90",
    "src/src_rp/storage.f90",
    "src/src_rp/molefractions_and_mixingratios.f90",
    "src/src_common/point_by_point_to_mixing_ratio.f90",
    "src/src_fcc/fluxes23.f90",
)


class WaterIsResolvedNotAssumed(unittest.TestCase):
    def test_no_flux_file_indexes_the_h2o_slot(self):
        pattern = re.compile(r"\(h2o\)|\bw_h2o\b|== h2o\b|/= h2o\b")
        for path in FLUX_FILES:
            hit = pattern.search(code(path))
            self.assertIsNone(
                hit,
                "%s still reads the h2o slot directly; it must resolve the "
                "primary water (PrimaryWaterSlot) or the gas's own moisture "
                "reference (E2Col(gas)%%moist_ref)" % path)

    def test_the_resolver_is_record_derived(self):
        source = read("src/src_common/gas4_output_units.f90")
        body = source[source.index("integer function PrimaryWaterSlot"):]
        body = body[:body.index("end function PrimaryWaterSlot")]
        self.assertIn("GasSlotIsWater(gas)", body)
        self.assertIn("do gas = firstGas, lastGas", body)
        self.assertIn("PrimaryWaterSlot = 0", body,
                      "a project with no water must resolve to 0, so callers "
                      "can report the quantity as not performed rather than "
                      "computing it from whatever occupies the slot")

    def test_the_water_constant_is_not_a_slot_lookup(self):
        """MW(slot) is superseded per record - including slots one to four.

        MW(h2o) was used at 27 sites meaning the molecular weight of water. A
        record naming another species at slot 6 made it 44.01e-3 and scaled
        every LE, E and ET by 2.44.
        """
        self.assertIn("real(kind = sgl), parameter :: MW_H2O = 18.02e-3",
                      read("src/src_common/m_common_global_var.f90"),
                      "MW_H2O must stay `sgl` with that literal, or the "
                      "promotion to double changes and every flux moves")
        for path in FLUX_FILES + ("src/src_common/define_all_var_set.f90",
                                  "src/src_rp/set_timelags.f90"):
            self.assertNotIn("MW(h2o)", code(path),
                             "%s reads MW by slot where it means water" % path)

    def test_the_dilution_correction_uses_each_gas_own_water(self):
        """Two passes: water first, then each gas by the water it names.

        Both this and the mole-fraction conversion were one hard-coded water
        block followed by a loop diluting every gas by Stats%r(h2o).
        """
        for path in ("src/src_common/point_by_point_to_mixing_ratio.f90",
                     "src/src_rp/molefractions_and_mixingratios.f90"):
            source = code(path)
            self.assertIn("E2Col(gas)%moist_ref", source,
                          "%s must dilute each gas by the water it names"
                          % path)
            self.assertIn("GasSlotIsWater", source)

    def test_the_same_instrument_guard_survives(self):
        """moist_ref falls back to "the first H2O anywhere" when a gas's own
        analyser has none. Diluting analyser B's gas with analyser A's water
        would be a new defect, not a fix."""
        source = code("src/src_common/point_by_point_to_mixing_ratio.f90")
        self.assertIn("E2Col(gas)%instr%model /= E2Col(msl)%instr%model",
                      source)


class SpeciesPropertiesAreKeyedOnSpecies(unittest.TestCase):
    def test_tube_adsorption_is_chosen_by_species(self):
        """Lambda is a property of the molecule.

        Keyed on the slot, slots past the fourth matched no case at all and
        left it undefined; and a project with water elsewhere gave the gas at
        slot 6 water's adsorption curve.
        """
        source = code("src/src_common/bpcf_massman_00.f90")
        self.assertNotIn("case (co2, gas4)", source)
        self.assertIn("case ('H2O')", source)
        self.assertIn("case default", source,
                      "an unrecognised species must fall back, not leave "
                      "Lambda undefined")

    def test_the_measured_flow_rate_reaches_every_gas(self):
        """It drives tube velocity, Reynolds number and the transfer function.

        Bounded at the fourth slot, gases on one analyser were corrected with
        two different flow rates depending on which slot they occupied: the
        MIRO's N2O got 1.04877 where its CO2, H2O and COS got 1.04798.
        """
        source = read("src/src_rp/eddyflow-rp_main.f90")
        start = source.index("replace instrument")
        block = source[start:source.index("end do", source.index(
            "E2Col(i)%instr%tube_f = UserStats%Mean(j)", start))]
        self.assertIn("do i = firstGas, lastGas", block)


if __name__ == "__main__":
    unittest.main()
