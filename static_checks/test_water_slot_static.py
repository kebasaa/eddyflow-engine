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
#>
#> fluxes1.f90 was missing from this list, and that omission is the whole
#> reason it kept the defect: it corrected LE, E and ET with BPCF%of(w_h2o)
#> and invalidated them on Flux0%gas(h2o) while its own comment claimed to
#> use "the primary H2O slot". The list is the check - a file absent from it
#> is not being checked at all, so adding a water-derived computation
#> anywhere means adding its file here.
FLUX_FILES = (
    "src/src_rp/flux_params.f90",
    "src/src_rp/fluxes0_rp.f90",
    "src/src_rp/fluxes1_rp.f90",
    "src/src_rp/fluxes23_rp.f90",
    "src/src_rp/storage.f90",
    "src/src_rp/molefractions_and_mixingratios.f90",
    "src/src_common/point_by_point_to_mixing_ratio.f90",
    "src/src_fcc/fluxes1.f90",
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
        source = read("src/src_common/gas_slot_resolution.f90")
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

    def test_molecular_weight_and_diffusivity_default_by_species(self):
        """A record that names a species but does not quantify it.

        The fallback was gated on `slot > gas4`, so the first four slots kept
        whatever the compile-time table held: water declared at slot 9 was
        given N2O's molecular weight and diffusivity, and a gas declared at
        slot 6 was given water's.
        """
        source = read("src/src_common/write_processing_project_variables.f90")
        block = source[source.index("g mol-1 -> kg mol-1"):]
        block = block[:block.index("end do")]
        self.assertNotIn("else if (slot > gas4)", block,
                         "the MW/Dc fallback still keys on slot position")
        self.assertIn("DefaultMolecularWeight(EddyFlowProj%gas(i)%var)", block)
        self.assertIn("DefaultDiffusivity(EddyFlowProj%gas(i)%var)", block)
        for name in ("DefaultMolecularWeight", "DefaultDiffusivity"):
            body = source[source.index("function %s(var)" % name):]
            body = body[:body.index("end function %s" % name)]
            self.assertIn("case ('H2O')", body)
            self.assertIn("case default", body,
                          "%s must give an unrecognised species a usable "
                          "number, never zero" % name)

    def test_a_gas_converts_with_its_own_records_weight(self):
        """Not with the weight of the first record naming that species.

        HistoricGasSlot maps a column's species name to a fixed slot, so a
        site measuring CO2 on two analysers converted both with record one's
        molecular weight. base_mw proves it: its second CO2 scales by
        44.01/30.0, its own record's weight, not by 44.01/90.0.
        """
        source = read("src/src_common/define_all_var_set.f90")
        block = source[source.index("case('co2', 'ch4', 'n2o')"):]
        block = block[:block.index("case('h2o')")]
        self.assertIn("if (gasSlot > 0) then", block)
        collapsed = " ".join(block.split())
        self.assertIn("N, j, & gasSlot)", collapsed,
                      "the conversion must be handed the column's own slot")

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


class EveryHygrometerIsTreatedAsWater(unittest.TestCase):
    """Not just the one in record two.

    None of these is exercised by a fixture, and each says why in its own
    docstring. They are pinned here because a static check is the only thing
    standing between them and a silent revert.
    """

    def test_rh_class_sorting_covers_every_hygrometer(self):
        """Water is binned by humidity, everything else by month, because
        water's tube attenuation drifts with humidity - a property of the
        species.

        Unprovable until the binned (co)spectra file carries more than four
        gases: a second hygrometer sits past slot 8 and so has no binned
        spectra to sort at all.
        """
        source = code("src/src_fcc/spectra_sorting_and_averaging.f90")
        self.assertNotIn("if (gas /= h2o) then", source)
        self.assertIn("if (.not. GasSlotIsWater(gas)) then", source)

    def test_the_rh_cutoff_fit_resolves_its_slot(self):
        """The fit belongs to the primary water record.

        Its *result* - the exponential RegPar(dum, dum) - is one set of
        coefficients for the whole project, so a second hygrometer reuses the
        primary's RH dependence. That is a data-model limit, not a loop bound,
        and widening it is a separate change.
        """
        source = code("src/src_fcc/fit_rh_to_cutoff.f90")
        self.assertNotIn("RegPar(h2o,", source)
        self.assertIn("wsl = PrimaryWaterSlot()", source)
        self.assertIn("if (wsl < firstGas) return", source,
                      "with no water there is nothing to fit")

    def test_the_active_gas_search_window_follows_the_species(self):
        """Water's time-lag window is ten times the transit time, not two.

        It adsorbs on the tube wall, so its lag runs long - a property of the
        molecule. Written as mult(h2o) it widened whatever record two held,
        and a hygrometer declared elsewhere got the passive window instead.
        The lag it was looking for could then sit outside the range searched,
        and an unfound lag is reported as the default one, so the flux moves
        without anything reporting a failure.

        No fixture carries a hygrometer outside record two *and* runs the
        time-lag optimiser, so this is pinned here and nowhere else.
        """
        source = code("src/src_rp/adjust_timelag_opt_settings.f90")
        self.assertNotIn("mult(h2o)", source)
        self.assertIn("GasSlotIsWater(gas)", source,
                      "the active-gas window must be chosen by species")

    def test_the_oxygen_correction_covers_every_hygrometer(self):
        """Krypton and Lyman-alpha instruments, each with its own ko/kw.

        Genuinely a water-vapour correction - oxygen absorbs in the band the
        instrument uses for water - so the H2O assumption stays and only the
        slot is resolved. No fixture carries a krypton, so this is pinned
        here and nowhere else.
        """
        source = code("src/src_rp/fluxes1_rp.f90")
        block = source[source.index("open_path_krypton"):]
        block = source[:source.index("open_path_krypton")].rsplit("do msl", 1)
        self.assertEqual(len(block), 2,
                         "the correction must loop over slots, not apply to "
                         "a single resolved one")
        body = source[source.index("do msl = firstGas, lastGas"):]
        body = body[:body.index("end do")]
        self.assertIn("GasSlotIsWater(msl)", body)
        self.assertIn("E2Col(msl)%Instr%ko", body,
                      "each hygrometer must use its own extinction "
                      "coefficients, never another instrument's")


if __name__ == "__main__":
    unittest.main()
