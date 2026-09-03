"""An analyser is found by its model, never by the slot it happens to occupy.

Three gates name a *particular instrument* rather than a particular gas: the
LI-7700 multiplier block, and the AGC/RSSI columns whose own heading says
LI-7200 or LI-7500. All three asked a fixed slot - E2Col(ch4) for the 7700,
E2Col(co2) for the other two - which is that analyser only when the project
happens to order its records that way.

Three ways that went wrong, none of which reports anything:

  two analysers   A site running a 7200 and a 7500 of different firmware had
                  both columns decided by whichever one held slot five, so one
                  of the two was headed AGC while carrying RSSI, or the sign
                  of its value was inverted.

  neither at five A project whose first record sits on a third instrument -
                  a MIRO, say - had both decided by something that is neither.

  header vs row   init_outfiles_rp writes the heading and write_out_fluxnet
                  writes the value. Both read slot five, so they agreed by
                  coincidence; resolving only one of them would have made the
                  column name disagree with the column.

The resolution is GasSlotByInstrModel and its companion InstrSwVerFor, which
hold the absent-instrument rule in one place: no such analyser gives version
zero, which compares older than any threshold - the arm an unpopulated
E2Col(co2)%instr%sw_ver already took.

No fixture carries two LI-COR analysers of different firmware, so nothing here
is provable by regression and this file is the only guard.
"""

from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]

HELPERS = "src/src_common/gas_slot_resolution.f90"

#> Every file that decides something from a named analyser's identity.
CALLERS = (
    "src/src_rp/write_out_fluxnet.f90",
    "src/src_rp/init_outfiles_rp.f90",
    "src/src_rp/interpret_diagnostics.f90",
)


def read(path):
    return (ROOT / path).read_text(encoding="utf-8")


def code(path):
    """Source with comment lines dropped, so a note *about* a removed
    construct does not read as the construct itself."""
    return "\n".join(ln for ln in read(path).splitlines()
                     if not ln.lstrip().startswith("!"))


class TheResolversExist(unittest.TestCase):
    def test_the_slot_search_is_record_derived(self):
        source = read(HELPERS)
        body = source[source.index("integer function GasSlotByInstrModel"):]
        body = body[:body.index("end function GasSlotByInstrModel")]
        self.assertIn("do gas = firstGas, lastGas", body)
        self.assertIn("index(E2Col(gas)%instr%model, fragment)", body,
                      "matching must be by substring - the model string "
                      "carries a trailing revision")
        self.assertIn("GasSlotByInstrModel = 0", body,
                      "a site with no such analyser must resolve to 0, so the "
                      "caller can take the absent arm rather than reading "
                      "whatever occupies the slot")

    def test_the_version_lookup_defaults_to_zero(self):
        source = read(HELPERS)
        body = source[source.index("function InstrSwVerFor(fragment)"):]
        body = body[:body.index("end function InstrSwVerFor")]
        self.assertIn("ver = SwVerType(0, 0, 0)", body,
                      "version zero compares older than any threshold, which "
                      "is what an unpopulated sw_ver did")
        self.assertIn("GasSlotByInstrModel(fragment)", body,
                      "one search, not two")

    def test_both_are_declared_for_callers(self):
        interfaces = read("src/src_common/interfaces_1.inc")
        self.assertIn("integer function GasSlotByInstrModel(fragment)", interfaces)
        self.assertIn("function InstrSwVerFor(fragment) result(ver)", interfaces)


class NoCallerReadsAFixedSlot(unittest.TestCase):
    def test_no_caller_asks_a_slot_for_a_firmware_version(self):
        for path in CALLERS:
            source = code(path)
            for slot in ("co2", "h2o", "ch4", "gas4"):
                self.assertNotIn(
                    "E2Col(%s)%%instr%%sw_ver" % slot, source,
                    "%s reads firmware from a fixed slot; it must ask the "
                    "analyser the column names, through InstrSwVerFor" % path)

    def test_the_methane_multipliers_are_written_per_gas(self):
        """No gate at all now: each gas writes its own multipliers.

        This began as a slot problem - keyed on E2Col(ch4), the writer and the
        scan that computes the multipliers could disagree, and a 7700 on any
        other record produced multipliers the writer suppressed. Asking
        GasSlotByInstrModel instead fixed that, but it still asked ONE
        question of the whole site, and the answer was one set of A/B/C for
        however many LI-7700s there were.

        Per gas there is nothing left to disagree: the writer emits
        Mul7700(gas) for each gas in the layout, and a gas with no LI-7700
        holds the error value, which AddFloatDatumToDataline writes as the
        missing-value token. A site gate would now be the bug.
        """
        source = code("src/src_rp/write_out_fluxnet.f90")
        self.assertNotIn("E2Col(ch4)%Instr%model", source)
        self.assertNotIn("GasSlotByInstrModel('li7700')", source,
                         "a single site-wide question cannot answer for two "
                         "LI-7700s; the writer emits Mul7700(gas) per gas")
        self.assertIn("Mul7700(FluxnetLayoutSlots(gas))%A", source)
        self.assertIn("Mul7700(FluxnetLayoutSlots(gas))%B", source)
        self.assertIn("Mul7700(FluxnetLayoutSlots(gas))%C", source)

    def test_each_analyser_column_is_labelled_from_its_own_firmware(self):
        """Each analyser named from its own firmware.

        The retired fixed-format header built both column names from one
        test, so the two could not disagree even when the instruments did.
        It stays two columns, but each is named from the analyser it
        describes.
        """
        source = code("src/src_rp/init_outfiles_rp.f90")
        for model in ("li7200", "li7500"):
            self.assertIn("InstrSwVerFor('%s')" % model, source)
        self.assertNotIn(
            "'RSSI_LI-7200,RSSI_LI-7500,variances", source,
            "the header names both analysers from one test")

    def test_the_row_and_the_header_ask_the_same_question(self):
        """They are in different files and were agreeing by coincidence."""
        row = code("src/src_rp/write_out_fluxnet.f90")
        header = code("src/src_rp/init_outfiles_rp.f90")
        for model in ("li7200", "li7500"):
            self.assertIn("InstrSwVerFor('%s')" % model, row)
            self.assertIn("InstrSwVerFor('%s')" % model, header)


if __name__ == "__main__":
    unittest.main()
