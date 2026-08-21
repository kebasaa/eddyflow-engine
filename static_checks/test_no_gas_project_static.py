"""An anemometer measuring no gas is a site, not an error.

`base_no_gas` sets `gas_num=0`: wind and sonic temperature, nothing else. RP
handled it and wrote a column-aligned FLUXNET file. FCC then rejected every
record and died with "No valid data records found in the essentials file".

The cause was one guard. `SelectFluxnetGasSlots` built its layout only

    if (EddyFlowProj%gas_num > 0) then
        call FluxnetLayoutGasSlots(FluxnetLayoutSlots, nFluxnetLayoutSlots)

but `FluxnetLayoutGasSlots` returns three synthetic CO2/H2O/CH4 slots whatever
`gas_num` is - FLUXNET requires those columns to exist however little the site
measures - and `ReadExRecord` calls it with no guard at all. So RP emitted no
gas column families and FCC expected three of everything. Every field after the
`_NR` block was out of step, and the read failed on the first integer item,
`measure_type`, landing on the real `W_SKW`: "Bad integer for item 78 in list
input".

Header and row were short by the same block, so they agreed with each other and
the column-count check passed. The fault only became visible one executable
later, as a message about data records that named nothing to do with columns.
That is the shape worth guarding against: not a mismatch between two files, but
two writers agreeing on something the third party reads differently.

The site is told, too. A project with no analyser gets Warning(108) once,
naming what is computed and what is not, so the error-code columns are
explained before they are found.
"""

from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]

FLX_HDR = "src/src_rp/init_fluxnet_file_rp.f90"
RESOLUTION = "src/src_common/gas_slot_resolution.f90"
RP_MAIN = "src/src_rp/eddyflow-rp_main.f90"
HANDLER = "src/src_common/exception_handler.f90"
LEGEND = "src/src_rp/write_column_legend.f90"
OUTFILES = "src/src_rp/init_outfiles_rp.f90"


def read(path):
    return (ROOT / path).read_text(encoding="utf-8")


def code(path):
    return "\n".join(ln for ln in read(path).splitlines()
                     if not ln.lstrip().startswith("!"))


def routine(path, name):
    src = code(path)
    start = src.index("subroutine %s" % name)
    return src[start: src.index("end subroutine %s" % name)]


class TheLayoutIsBuiltWhateverTheGasCount(unittest.TestCase):
    def setUp(self):
        self.block = routine(FLX_HDR, "SelectFluxnetGasSlots")

    def test_the_layout_call_is_not_guarded_on_gas_num(self):
        """The reader calls FluxnetLayoutGasSlots unguarded. A writer that
        guards it disagrees with the reader about how wide the row is."""
        idx = self.block.index("call FluxnetLayoutGasSlots(")
        before = self.block[:idx]
        #> Nothing may open a gas_num test that is still open at the call.
        opens = before.count("if (EddyFlowProj%gas_num > 0) then")
        closes = before.count("end if")
        self.assertEqual(
            0, opens,
            "the layout call sits inside a gas_num guard again; a no-gas "
            "project then writes no gas columns while FCC expects three")

    def test_the_present_gas_list_stays_guarded(self):
        """That one is the gases actually measured, and is legitimately empty
        - unlike the layout, which FLUXNET requires either way."""
        self.assertIn("if (EddyFlowProj%gas_num > 0) then", self.block)
        self.assertLess(self.block.index("call FluxnetLayoutGasSlots("),
                        self.block.index("if (EddyFlowProj%gas_num > 0) then"))

    def test_a_record_less_slot_still_gets_a_species(self):
        """A required variable no record names is carried on a slot past the
        records, so the species cannot be read off a record."""
        self.assertIn("FluxnetRequiredOrder(n)", self.block)


class TheResolverAnswersForZeroGases(unittest.TestCase):
    """The writer can only stop guarding because the resolver copes."""

    def test_the_required_species_do_not_depend_on_a_record(self):
        block = routine(RESOLUTION, "FluxnetLayoutGasSlots")
        self.assertIn("required = (/ 'CO2', 'H2O', 'CH4' /)", block)
        #> The synthetic branch: nothing names it, so it takes a slot past the
        #> configured records rather than being dropped.
        self.assertIn("slots(nslots) = firstGas + nrec + nsynth", block)


class TheSiteIsTold(unittest.TestCase):
    def test_the_warning_is_raised_for_a_no_gas_project(self):
        self.assertIn("if (EddyFlowProj%gas_num <= 0) call ExceptionHandler(108)",
                      code(RP_MAIN))

    def test_it_is_raised_once_and_not_per_period(self):
        """Outside the averaging-period loop, so no save-guard is needed and
        none is relied on. Ninety-seven copies of a warning is a silence."""
        src = code(RP_MAIN)
        at = src.index("call ExceptionHandler(108)")
        loop = src.index("do while (.not. EndOfPeriod)") \
            if "do while (.not. EndOfPeriod)" in src else len(src)
        self.assertLess(at, loop,
                        "the warning must be raised before the period loop")

    def test_the_message_separates_what_runs_from_what_does_not(self):
        text = read(HANDLER)
        block = text[text.index("case(108)"):]
        block = block[: block.index("end select")]
        self.assertIn("Computed:", block)
        self.assertIn("Not computed:", block)
        self.assertIn("u*", block)
        self.assertIn("error label", block,
                      "the message has to explain the columns of -9999, which "
                      "is what the user will actually go looking for")


class TheLegendNamesTheColumns(unittest.TestCase):
    """Once a suffix depends on record order, nothing in `h2o_2_flux` says
    which analyser it came from."""

    def setUp(self):
        self.block = routine(LEGEND, "WriteColumnLegend")

    def test_every_tag_comes_from_the_naming_helpers(self):
        """Rebuilding the names here would put a header/row divergence in the
        one file whose job is to explain the headers."""
        self.assertIn("call FullOutputGasTags(gas_tag)", self.block)
        self.assertIn("FluxnetLayoutTags(j)", self.block)
        self.assertIn("call WaterOutSlots(waterSlots, waterTags, nWater)",
                      self.block)

    def test_it_does_not_spell_a_suffix_itself(self):
        for invented in ("'_1'", "'_2'", "// '_' //"):
            self.assertNotIn(invented, self.block,
                             "the legend must report the tags, not derive "
                             "them a second time")

    def test_it_is_written_unconditionally(self):
        """A legend that might not be there is one a reader cannot rely on."""
        src = code(OUTFILES)
        at = src.index("call WriteColumnLegend()")
        line_start = src.rindex("\n", 0, at)
        self.assertNotIn("if ", src[line_start:at],
                         "the legend call must not be gated on an output flag")

    def test_it_says_which_record_is_primary(self):
        self.assertIn("designatedWater = PrimaryWaterSlot()", self.block)
        self.assertIn("designatedCarbon = PrimaryCarbonSlot()", self.block)

    def test_a_required_species_nobody_measures_is_still_explained(self):
        """Those columns exist because FLUXNET requires them and carry the
        error label throughout. Omitting the row leaves a reader wondering."""
        self.assertIn("'not measured'", self.block)


if __name__ == "__main__":
    unittest.main()
