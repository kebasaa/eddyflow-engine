from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]


def read(path):
    return (ROOT / path).read_text(encoding="utf-8")


class Gas4FullOutputUnitStaticTests(unittest.TestCase):
    def test_shared_helper_defines_pmol_and_nmol_scales(self):
        source = read("src/src_common/gas4_output_units.f90")
        self.assertIn("subroutine Gas4FullOutputUnits", source)
        self.assertIn("case ('ppb', 'nmol_mol', 'nmol/mol')", source)
        self.assertIn("flux_scale = 1d3", source)
        self.assertIn("dens_scale = 1d6", source)
        self.assertIn("case ('pmol_mol', 'pmol/mol')", source)
        self.assertIn("flux_scale = 1d6", source)
        self.assertIn("dens_scale = 1d9", source)
        self.assertIn("flux_scale = 1d0", source)
        self.assertIn("dens_scale = 1d0", source)

    def test_rp_full_output_uses_shared_helper(self):
        header = read("src/src_rp/init_outfiles_rp.f90")
        writer = read("src/src_rp/write_out_full.f90")
        # The argument is FourthGasUnitIn(), not E2Col(gas4)%unit_in. Both sites
        # run before the gas records have filled E2Col, so reading E2Col here
        # silently fell through to the default umol branch; the helper consults
        # the metadata column first and keeps E2Col only as its fallback.
        self.assertIn("call Gas4FullOutputUnits(FourthGasUnitIn()", header)
        self.assertIn("call Gas4FullOutputUnits(FourthGasUnitIn()", writer)
        self.assertNotIn("case ('ppb', 'nmol_mol')", header)
        self.assertNotIn("case ('ppb', 'nmol_mol')", writer)
        self.assertNotIn("case ('pmol_mol')", header)
        self.assertNotIn("case ('pmol_mol')", writer)

    def test_fcc_initializes_gas4_units_from_metadata(self):
        source = read("src/src_fcc/read_ini_fcc.f90")
        self.assertIn("call InitializeGas4FullOutputUnitsFcc()", source)
        self.assertIn("subroutine InitializeGas4FullOutputUnitsFcc", source)
        self.assertIn("call ReadMetadataFile(MetadataCol, AuxFile%metadata", source)
        self.assertIn("gas4_col = EddyFlowProj%col(gas4)", source)
        self.assertIn("gas4_unit = MetadataCol(gas4_col)%unit_out", source)
        self.assertIn("gas4_unit = MetadataCol(gas4_col)%unit_in", source)
        self.assertIn("call Gas4FullOutputUnits(gas4_unit, gas4_full_flux_sc, gas4_full_dens_sc", source)

    def test_fcc_full_output_scales_gas4_fields_only(self):
        source = read("src/src_fcc/write_out_full_fcc.f90")
        # The per-gas fluxes are an array now - Flux3%gas(gas4), not Flux3%gas4 -
        # because a project can configure more gases than the four that once had
        # named members.
        for token in (
            "Flux3%gas(gas4) * gas4_full_flux_sc",
            "lEx%rand_uncer(gas4) * gas4_full_flux_sc",
            "lEx%Stor%of(gas) * 1d-3 * gas4_full_flux_sc",
            "lEx%rot_w * lEx%d(gas) * 1d3 * gas4_full_flux_sc",
            "lEx%d(gas) * gas4_full_dens_sc",
            "lEx%chi(gas) * gas4_full_flux_sc",
            "lEx%r(gas) * gas4_full_flux_sc",
            "lEx%Flux0%gas(gas4) * gas4_full_flux_sc",
        ):
            self.assertIn(token, source)

        # The leak these guards exist to catch - another gas picking up the
        # fourth gas's scale - would now be spelled %gas(ch4). Written against
        # the retired %ch4 member they matched text that can no longer exist in
        # any form, so they could not have caught anything.
        self.assertNotIn("Flux3%gas(ch4) * gas4_full_flux_sc", source)
        self.assertNotIn("lEx%Flux0%gas(ch4) * gas4_full_flux_sc", source)


if __name__ == "__main__":
    unittest.main()
