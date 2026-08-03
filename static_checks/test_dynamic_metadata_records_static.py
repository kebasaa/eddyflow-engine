"""Every configured gas can be named in a dynamic metadata file.

RetrieveDynamicMetadata read its per-analyser overrides from a fixed index
table - co2_irga_manufacturer = 20 through gas4_irga_tau = 75 - via four
near-identical fifteen-statement blocks. The table stops there, so a project
with more than four gases had no name it could give the fifth analyser:
nothing in the file could reach it.

The four blocks had already drifted apart, which is what four copies of the
same fifteen statements do: the fourth read gas4_measure_type into
%instr(gas4)%nsep, a real, with a list-directed read of a character token.

The header is now matched per gas record, under the record's own label, so a
COS record answers to `cos_irga_model`. The four historical spellings keep
working as aliases through GasSlotFromDynMDTag - the same resolver the drift
subsystem uses for `<gas>_ref`, so the two cannot disagree about what
`n2o_irga_model` names.

The field list lives in one place, DynMDGasFieldNames, because the matcher
and the reader both walk it and a field added to one and not the other is
silently ignored.

NOT covered by this change, and worth knowing: ExtractUsableMetadataFromDynamic
propagates a whole instrument record between gases that report the same model
("consideration of gases from same analyser"). So a per-gas override on a
shared analyser is homogenised across every gas on it - setting
cos_irga_vertical_separation on a project whose gases all sit on one analyser
moves all of them. That predates this work; it only became reachable now that
those gases can be named at all.
"""

from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]

READER = "src/src_rp/retrieve_dynamic_metadata.f90"
HEADER = "src/src_rp/init_dynamic_medata.f90"
HELPER = "src/src_common/gas4_output_units.f90"

FIXTURES = ("tests/regression/base_dynmd.eddyflow",
            "tests/regression/base_dynmd_n_gas.eddyflow")


def code(path):
    return "\n".join(ln for ln in (ROOT / path).read_text(encoding="utf-8",
                                                          errors="replace").splitlines()
                     if not ln.lstrip().startswith("!"))


class TheFieldListLivesInOnePlace(unittest.TestCase):
    def test_the_helper_names_every_field(self):
        src = code(HELPER)
        self.assertIn("subroutine DynMDGasFieldNames", src)
        body = src.split("subroutine DynMDGasFieldNames")[1].split("end subroutine")[0]
        for suffix in ("_irga_manufacturer", "_irga_model", "_measure_type",
                       "_irga_northward_separation", "_irga_eastward_separation",
                       "_irga_vertical_separation", "_irga_tube_length",
                       "_irga_tube_diameter", "_irga_tube_flowrate",
                       "_irga_kw", "_irga_ko", "_irga_hpath_length",
                       "_irga_vpath_length", "_irga_tau"):
            self.assertIn("'%s'" % suffix, body,
                          "%s is not in the shared field list" % suffix)

    def test_both_sides_walk_it(self):
        for path in (READER, HEADER):
            self.assertIn("call DynMDGasFieldNames", code(path),
                          "%s must take the field list from the helper, or a "
                          "field added to one side is silently ignored" % path)


class TheGasBlocksAreOneLoop(unittest.TestCase):
    def test_no_unrolled_per_slot_blocks_remain(self):
        src = code(READER)
        for slot in ("co2", "h2o", "ch4", "gas4"):
            self.assertNotIn("DynamicMetadataOrder(%s_irga" % slot, src,
                             "%s_irga_* is a fixed index into a table that "
                             "stops at four gases" % slot)
            self.assertNotIn("DynamicMetadata%%instr(%s)" % slot, src,
                             "%s names a slot, not a record" % slot)

    def test_the_reader_is_record_driven(self):
        src = code(READER)
        self.assertIn("DynMDGasOrder(gas, fld)", src)
        self.assertIn("do gas = firstGas, lastGas", src)

    def test_measure_type_lands_in_measure_type(self):
        """The fourth block read it into %nsep, a real, from a character
        token - the kind of divergence four copies produce."""
        src = code(READER)
        self.assertIn("DynamicMetadata%measure_type(gas)", src)

    def test_the_header_matcher_resolves_by_record(self):
        src = code(HEADER)
        self.assertIn("GasSlotFromDynMDTag", src,
                      "the same resolver the drift ref columns use, so the "
                      "two agree about what n2o_irga_model names")
        self.assertIn("DynMDGasOrder", src)


class FixturesCoverBothDirections(unittest.TestCase):
    def test_the_fixtures_exist_and_enable_the_file(self):
        """There was no dynamic-metadata fixture at all, which is how the
        fourth block's measure_type defect survived."""
        for rel in FIXTURES:
            path = ROOT / rel
            self.assertTrue(path.exists(), "%s must exist" % rel)
            text = path.read_text(encoding="utf-8", errors="replace")
            self.assertRegex(text, r"(?m)^use_dyn_md_file=1$")

    def test_the_backward_fixture_uses_the_legacy_spellings(self):
        md = (ROOT / "tests/regression/base_dynmd.metadata").read_text(
            encoding="utf-8", errors="replace")
        self.assertIn("co2_irga_model", md)
        self.assertIn("gas4_measure_type", md,
                      "the fourth slot's legacy spelling must keep working")

    def test_the_forward_fixture_names_a_gas_past_the_fourth(self):
        md = (ROOT / "tests/regression/base_dynmd_n_gas.metadata").read_text(
            encoding="utf-8", errors="replace")
        self.assertIn("co2_2_irga", md,
                      "a second record of the same species is addressed by "
                      "the _2 suffix, and no legacy spelling could reach it")


if __name__ == "__main__":
    unittest.main()
