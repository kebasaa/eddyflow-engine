"""The full output's unit basis is per gas, not per slot.

These checks used to pin the opposite: a single set of scalars resolved for
the fourth slot and applied only there. That was the reason the full output
carried no columns at all for gases past the fourth - there was nothing to
label them with. Each assertion below is the turned-around form of the one it
replaces, so re-introducing a scalar gas4 scale fails here.
"""

from pathlib import Path
import re
import unittest


ROOT = Path(__file__).resolve().parents[1]


def read(path):
    return (ROOT / path).read_text(encoding="utf-8")


class GasFullOutputUnitStaticTests(unittest.TestCase):
    def test_shared_helper_defines_pmol_and_nmol_scales(self):
        source = read("src/src_common/gas_slot_resolution.f90")
        self.assertIn("subroutine GasFullOutputUnits", source)
        self.assertIn("case ('ppb', 'nmol_mol', 'nmol/mol')", source)
        self.assertIn("flux_scale = 1d3", source)
        self.assertIn("dens_scale = 1d6", source)
        self.assertIn("case ('pmol_mol', 'pmol/mol')", source)
        self.assertIn("flux_scale = 1d6", source)
        self.assertIn("dens_scale = 1d9", source)
        self.assertIn("flux_scale = 1d0", source)
        self.assertIn("dens_scale = 1d0", source)

    def test_helper_resolves_every_configured_gas(self):
        source = read("src/src_common/gas_slot_resolution.f90")
        self.assertIn("subroutine GasFullOutputUnitsAll", source)
        self.assertIn("do gas = firstGas, lastGas", source)
        # Water is on the mmol basis internally, so it must not fall through to
        # the helper's umol default arm.
        self.assertIn("if (GasSlotIsWater(gas)) then", source)
        self.assertIn("flux_label(gas) = '[mmol+1s-1m-2]'", source)
        self.assertIn("call GasFullOutputUnits(GasUnitIn(gas)", source)

    def test_species_predicate_and_tags_come_from_the_record(self):
        source = read("src/src_common/gas_slot_resolution.f90")
        # Answered from the gas record, not from the slot number: a second
        # water record sits well past the historical h2o slot.
        self.assertIn("logical function GasSlotIsWater(gas_slot)", source)
        self.assertIn("species = EddyFlowProj%gas(rec)%var", source)
        # Repeated species are disambiguated, or a project measuring CO2 on two
        # analysers would emit two identical column families.
        self.assertIn("subroutine FullOutputGasTags(tags)", source)
        self.assertIn("repeat = repeat + 1", source)

    def test_label_and_unit_helpers_take_a_slot(self):
        source = read("src/src_common/define_all_var_set.f90")
        self.assertIn("function GasOutputLabel(gas_slot) result(label)", source)
        self.assertIn("function GasUnitIn(gas_slot) result(unit_in)", source)
        self.assertNotIn("function FourthGasLabel()", source)
        self.assertNotIn("function FourthGasUnitIn()", source)

    def test_rp_full_output_resolves_units_per_gas(self):
        header = read("src/src_rp/init_outfiles_rp.f90")
        writer = read("src/src_rp/write_out_full.f90")
        # Header and row writer share one call, so a column's label can never
        # describe a scale that was not applied to its value.
        for source in (header, writer):
            self.assertIn("call GasFullOutputUnitsAll(gas_flux_sc, gas_dens_sc", source)
            self.assertNotIn("Gas4FullOutputUnits(", source)
            self.assertNotIn("case ('ppb', 'nmol_mol')", source)
            self.assertNotIn("case ('pmol_mol')", source)

    def test_fcc_initializes_units_for_every_configured_gas(self):
        source = read("src/src_fcc/read_ini_fcc.f90")
        self.assertIn("call InitializeGas4FullOutputUnitsFcc()", source)
        self.assertIn("call ReadMetadataFile(MetadataCol, AuxFile%metadata", source)
        # Per record, not the retired singular col_gas4 slot.
        self.assertIn("gas_col = EddyFlowProj%gas(rec)%col", source)
        self.assertNotIn("EddyFlowProj%col(gas4)", source)
        self.assertIn("gas_unit = MetadataCol(gas_col)%unit_out", source)
        self.assertIn("gas_unit = MetadataCol(gas_col)%unit_in", source)
        self.assertIn("call GasFullOutputUnits(gas_unit, gas_full_flux_sc(gas)", source)

    def test_fcc_globals_are_per_slot(self):
        source = read("src/src_fcc/m_fx_global_var_mod.f90")
        for token in (
            "gas_full_flux_sc(GHGNumVar)",
            "gas_full_dens_sc(GHGNumVar)",
            "gas_full_flux_label(GHGNumVar)",
            "gas_full_conc_label(GHGNumVar)",
            "gas_full_mixr_label(GHGNumVar)",
            "gas_full_dens_label(GHGNumVar)",
        ):
            self.assertIn(token, source)

    def test_no_scale_is_conditional_on_a_gas_slot(self):
        """The defect class this whole change removes.

        A scale selected by `gas == gas4` names a position, not a species, so
        the same gas was reported on one basis from slot 8 and another from
        slot 9. Any surviving instance is that bug returning.
        """
        pattern = re.compile(r"gas\s*==\s*(gas4|ch4|n2o)\b")
        for path in (
            "src/src_rp/write_out_full.f90",
            "src/src_rp/init_outfiles_rp.f90",
            "src/src_fcc/write_out_full_fcc.f90",
            "src/src_fcc/init_out_files.f90",
        ):
            source = read(path)
            self.assertIsNone(
                pattern.search(source),
                f"{path} still selects a unit scale by gas slot",
            )
            self.assertNotIn("gas4_full_flux_sc", source)
            self.assertNotIn("gas4_full_dens_sc", source)

    def test_fcc_presence_covers_every_configured_gas(self):
        """Every FCC output loop gates on fcc_var_present.

        Left bounded at the fourth slot it made gases 5+ absent from the full
        output however wide the loops were; run to lastGas instead of to the
        record count it marked all 64 slots present and emitted an empty column
        family for each.
        """
        source = read("src/src_fcc/init_ex_vars.f90")
        self.assertIn("do k = 1, min(EddyFlowProj%gas_num, MaxNumGases)", source)
        self.assertNotIn("do gas = co2, gas4", source)

    def test_fcc_presence_is_configured_not_measured(self):
        """A gas is present because a record names a column for it.

        Derived from the essentials record instead, presence collapsed to
        whatever survived filtering: a gas whose data was discarded reported
        the error code as its measure type, so it lost its columns rather than
        carrying -9999 in them. A filtered gas is written as the error label;
        absence is not how this file says "no data".
        """
        source = read("src/src_fcc/init_ex_vars.f90")
        self.assertIn("EddyFlowProj%gas(k)%col > 0", source)
        self.assertNotIn("lEx%measure_type_int(gas) /= ierror", source)

    def test_every_unit_conversion_leaves_the_error_code_alone(self):
        """The error code is a number, so scaling it makes it data.

        `error` is -9999, and three arms of ConvertTraceGasUnits used to scale
        the whole column: a missing sample in a ppb column came out as -9.999,
        which is past CleanUpE2Set's -300 test, past every `/= error` check
        downstream, and into the flux as a plausible mixing ratio. Measured on
        base_slow_naive, N2O reported a mole fraction of -8954 - the mean of a
        column nine tenths full of the logger's fill - and a flux computed from
        it.

        Checked structurally rather than arm by arm: every assignment to fRaw
        inside the routine must sit under a `where`, so a new unit added later
        cannot reintroduce the hole.
        """
        source = read("src/src_common/define_all_var_set.f90")
        body = source[source.index("subroutine ConvertTraceGasUnits"):]
        body = body[:body.index("end subroutine ConvertTraceGasUnits")]

        unguarded = []
        guard_depth = 0
        for line in body.splitlines():
            stripped = line.strip()
            if stripped.startswith("!"):
                continue
            if stripped.startswith("where("):
                guard_depth += 1
            elif stripped.startswith("end where"):
                guard_depth = max(0, guard_depth - 1)
            elif "fRaw(1:N, j) =" in stripped and guard_depth == 0:
                unguarded.append(stripped)

        #> The one legitimate bare assignment is the round trip through
        #> LinearConversion, which carries the error code itself.
        unguarded = [u for u in unguarded if "DumVec" not in u]
        self.assertEqual(unguarded, [], f"unguarded unit scalings: {unguarded}")


if __name__ == "__main__":
    unittest.main()
