"""Which columns an analyser's diagnostic word may invalidate.

The list used to be spelled `case (co2:gas4, pi:pe)`, and both halves stopped
meaning what they said once the slots widened:

  - `co2:gas4` is the first four gas slots, so a gas past the fourth kept every
    record its analyser's diagnostic rejected.
  - `pi:pe` was instrument 1's cell pressure through to air pressure - three
    slots. With one cell block per instrument it spans instruments 2..8
    entirely, so their cell *temperatures* started being filtered on a rule
    instrument 1's `tc` has never been subject to. On a two-analyser site that
    wiped the second analyser's cell temperature outright, and the physics then
    silently fell back to instrument 1's conditions for every gas.

The predicate names quantities instead of a slot span, which is what the
original list described back when there was only one instrument's cell block.
"""

from pathlib import Path
import re
import unittest


ROOT = Path(__file__).resolve().parents[1]
SRC = "src/src_rp/filter_dataset_for_diagnostics.f90"


def read(path):
    return (ROOT / path).read_text(encoding="utf-8")


class DiagnosticFilterScope(unittest.TestCase):
    def test_the_slot_span_is_gone(self):
        source = read(SRC)
        self.assertNotIn(
            "case (co2:gas4, pi:pe)", source,
            "the diagnostic filter is back to a slot span: it misses gases 5+ "
            "and swallows instruments 2..8's cell temperatures")

    def test_the_loop_covers_every_gas(self):
        source = read(SRC)
        self.assertIn("do var = firstGas, pe", source)

    def test_the_predicate_names_quantities(self):
        source = read(SRC)
        self.assertIn("logical function DiagFilterable(slot)", source)
        # Gas slots.
        self.assertIn("slot >= firstGas .and. slot <= lastGas", source)
        # One cell pressure per instrument - offset 3, where `pi` sits.
        self.assertIn("mod(slot - firstCell, NumCellPerInstr) == 3", source)
        # Ambient T/P.
        self.assertIn("slot == te .or. slot == pe", source)

    def test_cell_temperatures_are_not_filterable(self):
        """The asymmetry that caused the bug.

        Offsets 0, 1 and 2 of an instrument's block are cell_t, int_t_1 and
        int_t_2. Instrument 1's are tc/ti1/ti2, which sit below `pi` and were
        never in the historical list; making instrument 2's filterable while
        instrument 1's are not is what wiped the second analyser.
        """
        source = read(SRC)
        for offset in (0, 1, 2):
            self.assertNotIn(
                f"mod(slot - firstCell, NumCellPerInstr) == {offset}", source,
                f"cell block offset {offset} is a temperature; filtering it "
                f"would treat a second analyser differently from the first")

    def test_the_missing_gas_report_covers_every_gas(self):
        source = read(SRC)
        self.assertIn("E2Col(firstGas:lastGas)%present", source)
        self.assertNotIn("E2Col(co2:gas4)%present", source)

    def test_the_loop_index_is_declared(self):
        """`i` was only ever used inside the unrolled arms.

        Widening the loop moved the arms under a new guard; without the
        declaration the file relies on implicit typing, which `implicit none`
        would reject and which no test would otherwise notice.
        """
        source = read(SRC)
        decls = re.findall(r"^\s*integer :: (.+)$", source, re.M)
        names = {n.strip() for line in decls for n in line.split(",")}
        self.assertIn("i", names)
        self.assertIn("var", names)



class CellConditionsCrossIntoFcc(unittest.TestCase):
    """FCC computes the published fluxes; the cell terms must reach it per gas.

    `lEx%Tcell`/`lEx%Pcell` are instrument 1's *and* carry the writer's degC and
    kPa gains, which the reader never inverts - so the two closed-path WPL cell
    terms were dividing by about 27 instead of 300, and by 0.07 instead of 70.
    The per-gas columns are SI and read back unchanged, which is the rule the
    NUM_GAS_INSTR block already follows for the same reason.
    """

    HEADER = "src/src_rp/init_fluxnet_file_rp.f90"
    WRITER = "src/src_rp/write_out_fluxnet.f90"
    READER = "src/src_common/read_ex_record.f90"
    FCC_EMIT = "src/src_fcc/write_out_fluxnet_fcc.f90"
    FCC_PHYS = "src/src_fcc/fluxes23.f90"

    def test_the_header_declares_the_three_groups(self):
        source = read(self.HEADER)
        for tag in ("'T_CELL_'", "'PA_CELL_'", "'W_PA_CELL_'"):
            self.assertIn(tag, source)

    def test_the_writer_emits_si(self):
        """No gain/offset: a converted value would be double-converted on read."""
        source = read(self.WRITER)
        for expr in (
            "AddFloatDatumToDataline(Ambient%Tcell_at(gas), csv_row, EddyFlowProj%err_label)",
            "AddFloatDatumToDataline(Ambient%Pcell_at(gas), csv_row, EddyFlowProj%err_label)",
        ):
            self.assertIn(expr, source)
        self.assertIn("Stats%cov(w, cellPressureSlot(gas))", source)

    def test_the_reader_accounts_for_them(self):
        source = read(self.READER)
        self.assertIn("+ 3 * nExGas", source)
        self.assertIn("lEx%Tcell_at(firstGas:lastCfg)", source)
        self.assertIn("lEx%Pcell_at(firstGas:lastCfg)", source)
        self.assertIn("lEx%cov_w_pcell(firstGas:lastCfg)", source)

    def test_fcc_echoes_them(self):
        source = read(self.FCC_EMIT)
        for member in ("lEx%Tcell_at(gas)", "lEx%Pcell_at(gas)", "lEx%cov_w_pcell(gas)"):
            self.assertIn(member, source)

    def test_fcc_physics_uses_the_per_gas_values(self):
        source = read(self.FCC_PHYS)
        self.assertIn("lEx%RhoCp * lEx%Tcell_at(gas)", source)
        self.assertIn("lEx%cov_w_pcell(gas)", source)
        self.assertIn("lEx%Pcell_at(gas)", source)
        # The scalars are instrument 1's and in the wrong units for physics.
        self.assertNotIn("lEx%RhoCp * lEx%Tcell)", source)
        self.assertNotIn("lEx%cov_w(pi)", source)

    def test_the_slot_helper_is_shared(self):
        """One definition, because the writer and the physics must agree.

        The writer puts this covariance in the file and the flux code consumes
        it; a second copy of the slot arithmetic is a silent mismatch.
        """
        self.assertIn(
            "integer function cellPressureSlot(gas) result(slot)",
            read("src/src_common/define_e2_set.f90"))
        for path in ("src/src_rp/fluxes23_rp.f90", "src/src_rp/write_out_fluxnet.f90"):
            self.assertIn("integer, external :: cellPressureSlot", read(path))


if __name__ == "__main__":
    unittest.main()


class SpectralCorrectionsReachEveryGas(unittest.TestCase):
    """Four independent gates kept gases 5+ from getting a correction factor.

    Each one alone was enough to leave `BPCF%of` at the error sentinel, so the
    symptom - `SCF = -9999` - looked the same however many were fixed. They are
    pinned together because that is how they have to be removed.
    """

    def test_transfer_functions_are_initialised_for_every_variable(self):
        """The one that made it fail under *every* method.

        BPTF is intent(out); a slot SetTransferFunctionsToValue skips is left
        undefined, and SpectralCorrectionFactors then finds no usable band-pass
        value. Nothing to do with the cospectra file.
        """
        source = read("src/src_common/bpcf_aux_subs.f90")
        self.assertIn("do var = u, lastGas", source)
        self.assertNotIn("do var = u, gas4", source)

    def test_the_analytic_transfer_function_covers_every_gas(self):
        source = read("src/src_common/bpcf_analytic_transfer_functions.f90")
        self.assertIn("case (firstGas:lastGas)", source)
        self.assertNotIn("case (co2, h2o, ch4, gas4)", source)

    def test_the_cospectral_model_covers_every_gas(self):
        """The analytic cospectrum was copied to three named slots."""
        source = read("src/src_common/bpcf_cospectral_models.f90")
        self.assertIn("do gas = firstGas, lastGas", source)
        self.assertNotIn("Cospectrum(:)%of(w_gas4) = Cospectrum(:)%of(w_co2)", source)

    def test_sensor_parameters_are_retrieved_for_every_gas(self):
        """Path lengths and response time drive the transfer functions.

        A slot left at the error sentinel gives a response time of -9999, the
        dynamic-response term collapses, and the factor comes out around 5500
        instead of about 1.05 - plausible-looking only if nobody looks.
        """
        source = read("src/src_rp/retrieve_sensor_params.f90")
        self.assertIn("do gas = firstGas, lastGas", source)
        self.assertNotIn("do gas = co2, gas4", source)

    def test_the_cospectra_file_names_every_gas(self):
        """Writer and reader must spell the column the same way.

        A name they disagree on is simply not imported, and the gas silently
        gets no correction factor.
        """
        helper = read("src/src_common/gas4_output_units.f90")
        self.assertIn("subroutine SpectralVarTags(tags)", helper)
        # The historical eight are the shipped file format and must not move.
        for name in ("'u'", "'v'", "'w'", "'ts'", "'co2'", "'h2o'", "'ch4'", "'gas4'"):
            self.assertIn("= " + name, helper)
        for path in ("src/src_rp/spectral_analysis.f90",
                     "src/src_common/bpcf_read_full_cos_wt.f90"):
            self.assertIn("call SpectralVarTags(", read(path))
        # The reader's compile-time table is gone.
        reader = read("src/src_common/bpcf_read_full_cos_wt.f90")
        self.assertNotIn("data covlabs(1:8)", reader)
        self.assertIn("covlabs(j) = 'cov(w_' // trim(vartags(j)) // ')'", reader)

    def test_level3_does_not_multiply_by_the_sentinel(self):
        for path in ("src/src_rp/fluxes23_rp.f90", "src/src_fcc/fluxes23.f90"):
            source = read(path)
            self.assertIn("if (BPCF%of(msl) /= error) then", source)
