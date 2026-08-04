"""The spectral assessment covers every configured gas, not the first four.

The chain that fits transfer functions - bin, sort, ensemble, fit, write,
read back, apply - was bounded at the fourth gas end to end. Gases past it
were never assessed, so they fell through to an analytic transfer function
while the output reported a correction factor either way. Nothing in the file
or the log said which gases had actually been fitted.

Widening it uncovered three defects that were inert only because the loops
stopped early, and each is pinned here:

  * the month/class table was built for CO2, CH4 and the fourth gas only, so
    a fifth gas indexed RegPar(gas, 0);
  * var_present was derived over every slot from a Flux0 array that is not
    reset between records, so a slot the project never declared tested as
    present;
  * the Horst & Lenschow separation correction indexed the by-role
    instrument array as `gas - 3`, which runs off its 8-entry end at the
    fifth gas and addresses an unrelated analyser before that.

The regression fixtures `base_n_gas_sa` / `base_n_gas_sa_short` exercise the
same ground at runtime: the first proves gases 5+ take a fitted transfer
function (SCF 2.47 against the analytic 1.05), the second proves a file with
too few blocks falls back rather than inventing one.
"""

from pathlib import Path
import re
import unittest


ROOT = Path(__file__).resolve().parents[1]

# The chain from binned spectra through to the applied correction factor.
CHAIN = (
    "src/src_fcc/normalize_mean_spectra_cospectra.f90",
    "src/src_fcc/cospectra_qaqc.f90",
    "src/src_fcc/fit_tf_models.f90",
    "src/src_fcc/fit_cospectral_models.f90",
    "src/src_fcc/available_mean_spectra_cospectra.f90",
    "src/src_fcc/spectra_sorting_and_averaging.f90",
    "src/src_fcc/add_to_cospectra_fit_dataset.f90",
    "src/src_fcc/ensemble_cospectra_by_stability.f90",
    "src/src_fcc/subtract_high_freq_noise.f90",
    "src/src_fcc/read_binned_file.f90",
    "src/src_fcc/spectral_assessment_diagnostics.f90",
    "src/src_fcc/output_spectral_assessment_results.f90",
    "src/src_fcc/read_spectral_assessment_file.f90",
)

WRITER = "src/src_fcc/output_spectral_assessment_results.f90"
READER = "src/src_fcc/read_spectral_assessment_file.f90"


def read(path):
    return (ROOT / path).read_text(encoding="utf-8")


def strip_comments(source):
    """Fortran comment lines, so a note *about* a removed construct does not
    read as the construct itself."""
    return "\n".join(ln for ln in source.splitlines()
                     if not ln.lstrip().startswith("!"))


class TheChainIsNotBoundedAtTheFourthGas(unittest.TestCase):
    def test_no_gas_loop_stops_at_gas4(self):
        # `co2, gas4` and `co2:gas4` as a loop bound or an array section.
        pattern = re.compile(r"\b(histCO2|u)\s*[,:]\s*histGas4\b")
        for path in CHAIN:
            hit = pattern.search(read(path))
            self.assertIsNone(
                hit,
                f"{path} still stops at the fourth gas, so gases 5+ are never "
                f"assessed and silently fall back to an analytic transfer "
                f"function",
            )


class TheWriterAndReaderAgreeOnTheBlockCount(unittest.TestCase):
    """The assessment file is a writer/reader pair over the same range.

    A count mismatch does not fail loudly: the next thing in the file is the
    exponential-fit section and a blind `read` of its title line succeeds, so
    the reader would consume it as transfer-function parameters.
    """

    def test_both_skip_water_by_species(self):
        """Both sides must leave out the same gases.

        The skip used to be `gas == h2o`, the historical slot - so a second
        hygrometer was written and read as a month-classed trace gas, while
        spectra_sorting_and_averaging RH-sorts every hygrometer. Asked of the
        record now, on both sides, so they still agree.

        The reader no longer walks the gas list at all - it is driven by the
        block headers in the file - so only the writer is checked for the
        range; see test_spectral_assessment_blocks_static.py for the reader.
        """
        for path in (WRITER, READER):
            source = read(path)
            self.assertNotIn("if (gas == h2o) cycle", source,
                             f"{path} must skip water by species, not by slot")
            self.assertIn("GasSlotIsWater(gas)", source,
                          f"{path} must skip water, whose cut-offs come from "
                          f"the RH class table, not from a gas block")
        writer = read(WRITER)
        self.assertIn("do gas = firstGas, lastGas", writer, WRITER)
        self.assertIn(
            "min(EddyFlowProj%gas_num, MaxNumGases)", writer,
            "the writer must stop at the declared gas count")

    def test_the_reader_checks_the_block_header(self):
        source = read(READER)
        self.assertIn("index(dataline, 'TFP')", source,
                      "the reader must verify it is looking at a block "
                      "header; skipping the line blindly consumes the next "
                      "section as parameters")
        self.assertIn("backspace(udf)", source,
                      "on a mismatch the peeked lines must be put back, or "
                      "the sections below no longer parse")

    def test_an_absent_gas_reads_as_unfitted_not_as_zero(self):
        """The was-this-configured guard for this widening.

        RegPar is zeroed before the read. A cut-off of zero is not a missing
        value - it is an infinitely aggressive correction, and it produced a
        correction factor of 2.6 for gases a short file did not carry.
        """
        source = read(READER)
        self.assertIn("RegPar(gas, JAN:DEC)%fc = error", source)
        self.assertIn("RegPar(gas, JAN:DEC)%Fn = error", source)

    def test_the_headers_name_the_species(self):
        """Every slot, named from its own record - no slot is a species.

        Slots are assigned by record order (slot = firstGas + i - 1), so the
        constants co2/h2o/ch4 name records one to three and say nothing about
        what those records declare. This used to pin the first three and name
        everything past the fourth "Gas 4".
        """
        self.assertIn("call SpectralGasNames", read(WRITER))
        self.assertIn("call SpectralGasNames",
                      read("src/src_fcc/spectral_assessment_diagnostics.f90"))

        source = read("src/src_common/gas4_output_units.f90")
        body = source[source.index("subroutine SpectralGasNames"):]
        body = body[:body.index("end subroutine SpectralGasNames")]
        self.assertIn("do gas = firstGas, lastGas", body,
                      "the naming loop must cover every configured slot")
        for literal in ("'co2'", "'h2o'", "'ch4'", "'gas4'"):
            self.assertNotIn(
                literal, body,
                f"SpectralGasNames names slot {literal} as a fixed species; "
                f"a project that orders its records differently puts a "
                f"different gas there")

    def test_the_diagnostics_report_names_no_fixed_species(self):
        source = read("src/src_fcc/spectral_assessment_diagnostics.f90")
        body = source[source.index("function GasName(gas)"):]
        body = body[:body.index("end function GasName")]
        for literal in ("'CO2'", "'H2O'", "'CH4'", "'Gas 4'"):
            self.assertNotIn(literal, body,
                             f"GasName still hard-codes {literal}")

    def test_the_fourth_gas_label_is_not_parsed_out_of_the_header(self):
        """g4lab searched the FLUXNET header for ',FCH4,' and took the next
        column, which is a garbage substring on a project without methane."""
        for path in ("src/src_fcc/init_ex_vars.f90",
                     "src/src_fcc/m_fx_global_var_mod.f90",
                     "src/src_fcc/fit_cospectral_models.f90",
                     WRITER):
            code = strip_comments(read(path))
            self.assertNotIn("g4lab", code, f"{path} still uses g4lab")
        self.assertNotIn(
            "FCH4", strip_comments(read("src/src_fcc/init_ex_vars.f90")),
            "the fourth gas's label is being recovered by searching the "
            "FLUXNET header for a methane column again")


class EveryConfiguredGasCanBeClassified(unittest.TestCase):
    """Gases past the fourth have no month-grouping tags of their own.

    Left at class 0 they are not merely unfitted: 0 is not a valid RegPar
    index, and the assessment could never fit them either, because every
    month would be written as `error`.
    """

    def test_gases_past_the_fourth_inherit_a_grouping(self):
        source = read("src/src_fcc/read_ini_fcc.f90")
        self.assertIn(
            "FCCsetup%SA%class(gas, JAN:DEC) = FCCsetup%SA%class(histCO2, JAN:DEC)",
            source,
            "gases past the fourth must inherit CO2's month grouping; the "
            "grouping bins the calendar, not the species",
        )

    def test_the_lookup_rejects_an_invalid_class(self):
        source = read("src/src_common/bpcf_bandpass_spectral_corrections.f90")
        self.assertIn("LocSetup%SA%class(gas, month) < 1", source,
                      "an unclassified gas must be treated as a missing fit; "
                      "class 0 indexes RegPar out of bounds")


class PresenceIsBoundedByTheDeclaredGasCount(unittest.TestCase):
    """A slot the project never declared must not test as present.

    Flux0 is not reset between records, so an undeclared slot holds 0 rather
    than the error sentinel. The phantom gas reached every var_present-gated
    loop; it surfaced as the whole spectral correction falling back to
    Moncrieff, because a phantom has no class and so no cut-off to look up.
    """

    def test_var_present_stops_at_the_declared_count(self):
        source = read("src/src_common/read_ex_record.f90")
        block = source[source.index("lEx%var_present = .false."):]
        block = block[:block.index("!> Units adjustments")]
        self.assertIn("min(EddyFlowProj%gas_num, MaxNumGases)", block)


class SeparationCorrectionUsesThePerGasInstrument(unittest.TestCase):
    """lEx%instr is the fixed by-role array; lEx%gas_instr is per slot."""

    def test_horst_lenschow_indexes_by_slot(self):
        path = "src/src_common/bpcf_additional_horst_lenschow_09.f90"
        source = read(path)
        self.assertNotIn("lEx%instr(igas)", source,
                         f"{path} indexes the 8-entry by-role instrument "
                         f"array by gas slot; it runs off the end at the "
                         f"fifth gas and mis-addresses before that")
        self.assertIn("lEx%gas_instr(gas)%vsep", source)
        self.assertIn("lEx%gas_instr(gas)%nsep", source)
        # The sonic is genuinely a role, and stays one.
        self.assertIn("lEx%instr(sonic)%height", source)


if __name__ == "__main__":
    unittest.main()
