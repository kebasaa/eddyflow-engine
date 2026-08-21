"""The AGC/RSSI columns the conditional eddy covariance screen needs.

Three faults met on one project, and this file pins all three.

A column named AGC or RSSI never reached UserCol: IsCustomOutputColumn listed
`agc` and `rssi` beside `ignore` and `flag_1`, so the engine dropped exactly the
columns two of its own loops then went looking for. CecSignalColumnFor and
SetLicorDiagnostics were dead code - the screen ran on nothing and RSSI77 was
always the error value, in both cases silently.

Which analyser a signal column belonged to was inferred, and how it was spelled
decided whether it counted at all: the comparison was case-sensitive, so `agc`
from another tool was not a signal strength and nothing said so. The project
file now carries `agc_<i>_*` records naming the column and its analyser, on the
same shape as the cell and diagnostic records.

And a record whose column has since been re-declared was still honoured. A
diagnostic column re-declared as AGC kept its diag_72 record, so a real
diagnostic elsewhere on that analyser made two records competing for the one
diagnostic slot and MetadataFileValidation aborted the run - over a record the
interface does not show and the user cannot remove. RecordNamesColumn makes such
a record inert instead.
"""

from pathlib import Path
import re
import unittest


ROOT = Path(__file__).resolve().parents[1]
GUI = ROOT.parent / "eddyflow-gui"

USER_SET = "src/src_common/define_user_set.f90"
USED = "src/src_common/define_used_variables.f90"
VARS = "src/src_rp/define_vars.f90"
NAMES = "src/src_common/record_names_column.f90"
E2SET = "src/src_common/define_e2_set.f90"
VALIDATION = "src/src_common/metadata_file_validation.f90"
SLOTS = "src/src_common/gas_slot_resolution.f90"
LICOR = "src/src_rp/set_licor_diagnostics.f90"
PROJVARS = "src/src_common/write_processing_project_variables.f90"
GLOBALS = "src/src_common/m_common_global_var.f90"
TYPEDEF = "src/src_common/m_typedef.f90"

#: Every copy of the predicate that decides what becomes a user column.
CUSTOM_COLUMN_FILES = (USER_SET, USED, VARS)


def read(path, base=ROOT):
    return (base / path).read_text(encoding="utf-8", errors="replace")


def unit(source, opener, closer):
    """The text of one program unit, by its opening and closing line."""
    start = source.index(opener)
    end = source.index(closer, start)
    return source[start:end]


class SignalStrengthColumnsReachTheUserSet(unittest.TestCase):
    """A signal-strength column is ordinary custom data.

    It is not a diagnostic, whatever its name suggests. Excluding it here is
    what made CecSignalColumnFor and SetLicorDiagnostics unreachable.
    """

    def test_no_copy_of_the_predicate_excludes_them(self):
        for path in CUSTOM_COLUMN_FILES:
            source = read(path)
            self.assertIn("logical function IsCustomOutputColumn", source,
                          path + " no longer holds the predicate")
            self.assertNotIn("'agc', 'rssi'", source,
                             path + " drops signal-strength columns again")

    def test_the_other_exclusions_are_untouched(self):
        """Only the two names moved. `ignore`, `not_numeric`, `none` and the
        two flags are still not custom output."""
        for path in CUSTOM_COLUMN_FILES:
            self.assertIn(
                "case ('ignore', 'not_numeric', 'none', 'flag_1', 'flag_2')",
                read(path), path)

    def test_the_reason_is_written_down(self):
        """A comment naming the consumers, so the exclusion is not reinstated
        the next time this list is read as a list of diagnostics."""
        for path in CUSTOM_COLUMN_FILES:
            self.assertIn("CecSignalColumnFor", read(path), path)


class TheProjectNamesItsSignalColumns(unittest.TestCase):
    """Records, not a variable name matched case-sensitively."""

    def test_the_records_are_sized_and_held(self):
        typedef = read(TYPEDEF)
        self.assertIn(
            "integer, parameter :: MaxNumAgcCols = MaxNumInstruments * 2",
            typedef)
        self.assertIn("type(MeasRecordType) :: agc(MaxNumAgcCols)", typedef)
        self.assertIn("integer :: agc_num", typedef)

    def test_the_tags_exist_at_their_own_origins(self):
        globals_ = read(GLOBALS)
        for tag in ("agc_num", "agc_1_col", "agc_1_var", "agc_1_instr",
                    "agc_16_col", "agc_16_instr"):
            self.assertIn("'%s'" % tag, globals_, tag)
        for name in ("agcNumTag", "agcRecOriginN", "agcRecOriginC",
                     "agcRecLeapN", "agcRecLeapC"):
            #> Whitespace-tolerant: the generator pads the leap
            #> constants into a column, so a literal " = " misses them.
            self.assertRegex(
                globals_, r"integer, parameter :: %s +=" % name, name)

    def test_appending_did_not_move_the_older_origins(self):
        """The tables are positional: a tag's identity is its index.

        agc_* is APPENDED, after the cec block, so nothing before it moves.
        Had it been inserted, every gas, cell and diagnostic record after the
        insertion point would read a different key with no compile error.
        """
        globals_ = read(GLOBALS)
        for name, value in (("gasNumTag", 33), ("cellNumTag", 34),
                            ("diagNumTag", 35), ("diagRecOriginC", 243),
                            ("cecRecOriginC", 275)):
            self.assertIn(
                "integer, parameter :: %s = %d" % (name, value), globals_,
                name + " moved - every record after it is now misaddressed")

    def test_the_count_is_read_and_clamped(self):
        source = read(PROJVARS)
        self.assertIn(
            "EddyFlowProj%agc_num = nint(EPPrjNTags(agcNumTag)%value)", source)
        self.assertIn(
            "EddyFlowProj%agc_num  = min(max(EddyFlowProj%agc_num,  0), "
            "MaxNumAgcCols)", source)


class TheScreenFindsTheColumnByRecord(unittest.TestCase):
    """Per gas, on that gas's own analyser."""

    def signal_column_for(self):
        return unit(read(SLOTS), "integer function CecSignalColumnFor",
                    "end function CecSignalColumnFor")

    def test_the_records_are_consulted_first(self):
        fn = self.signal_column_for()
        self.assertLess(fn.index("EddyFlowProj%agc(i)%col"),
                        fn.index("UserCol(j)%var /= 'AGC'"),
                        "the name scan runs before the records")

    def test_the_name_scan_survives_as_the_fallback(self):
        """A project written before the records existed states nothing but the
        variable name, so that is what is left to match on. Removing the
        fallback would silently unscreen every such project."""
        self.assertIn("UserCol(j)%var /= 'AGC' .and. UserCol(j)%var /= 'RSSI'",
                      read(SLOTS))

    def test_the_analyser_comes_from_the_metadata_not_the_record(self):
        """The record's `instr` and the column's metadata say the same thing
        when the interface wrote the file, and the metadata is the authority.

        Comparing id spellings is how this went wrong before: the label the
        table shows is not the id the records store, the two never match, and
        the failure is silent.
        """
        fn = self.signal_column_for()
        self.assertIn("UserCol(j)%instr%slot /= E2Col(gas_slot)%instr%slot", fn)
        self.assertNotIn("%agc(i)%instr", fn)

    def test_which_of_the_two_it_is_ignores_case(self):
        """RSSI is high-is-clean and AGC is high-is-dirty, compared against the
        same number in opposite directions, so getting this wrong keeps exactly
        the samples that should go."""
        fn = unit(read(SLOTS), "logical function CecSignalIsRssi",
                  "end function CecSignalIsRssi")
        self.assertIn("call lowercase(name)", fn)
        self.assertIn("== 'rssi'", fn)
        self.assertIn("/= 'agc'", fn)

    def test_the_essentials_read_the_same_records(self):
        source = read(LICOR)
        self.assertIn("EddyFlowProj%agc(k)%col /= UserCol(i)%orig_col", source)
        self.assertIn("call lowercase(name)", source)

    def test_each_analyser_gets_its_own_column(self):
        """Every branch used to `exit`, not only the one that matched its own
        analyser, so the first signal column of any kind ended the search and a
        site with two analysers filled one slot and left the other at error."""
        loop = unit(read(LICOR), "do i = 1, M", "end do")
        self.assertNotIn("exit", loop,
                         "the first signal column of any kind ends the search "
                         "again, so a second analyser gets none")
        for slot in ("AGC72", "AGC75", "RSSI77"):
            self.assertIn("Essentials%%%s == error" % slot, loop, slot)


class ARecordMustStillNameItsColumn(unittest.TestCase):
    """The metadata is the authority on what a column measures."""

    def test_the_predicate_exists_and_is_case_insensitive(self):
        source = read(NAMES)
        self.assertIn(
            "logical function RecordNamesColumn(col_var, record_var)", source)
        self.assertEqual(source.count("call lowercase("), 2,
                         "one side of the comparison is not normalised")

    def test_it_still_rejects_the_ignored_column(self):
        """The case this generalises. It was already exempt, and for the same
        reason: the metadata says the column holds nothing usable."""
        self.assertIn(
            "if (trim(col) == 'ignore' .or. trim(col) == 'not_numeric') return",
            read(NAMES))

    def test_both_spellings_of_the_anemometer_record_agree(self):
        """The one measurement whose record slug is not the metadata's own
        spelling of it. Both are accepted rather than one declared correct,
        because project files carrying either are already in the field."""
        source = read(NAMES)
        self.assertIn(
            "if (trim(col) == 'diag_anem') col = 'anemometer_diagnostic'",
            source)
        self.assertIn(
            "if (trim(rec) == 'diag_anem') rec = 'anemometer_diagnostic'",
            source)

    def test_every_consumer_of_a_record_asks(self):
        """Validation alone is not enough. A stale record that collides with
        nothing wins its slot instead, and the engine then decodes whatever the
        column now holds as a diagnostic bitfield."""
        for path, call in ((E2SET, "EddyFlowProj%diag(i)%var)) cycle"),
                           (USED, "RecordNamesColumn(LocCol(rec%col)%var, rec%var)"),
                           (VALIDATION, "RecordNamesColumn(Cols(rec%col)%var, rec%var)")):
            self.assertIn("RecordNamesColumn", read(path), path)
            self.assertIn(call, read(path), path + " declares it but never asks")

    def test_the_slot_mapping_asks_before_it_maps(self):
        """DefineUsedVariables is where a record takes a slot, and it is
        last-writer-wins."""
        source = read(USED)
        block = unit(source, "if (EddyFlowProj%diag_num > 0) then", "end if")
        self.assertIn(
            "if (.not. RecordStillNamesIt(EddyFlowProj%diag(i))) cycle", block)


class TheAnemometerDiagnosticSurvivesItsRecords(unittest.TestCase):
    """ApplyCellDiagRecords took only `anemometer_diagnostic` while the
    interface writes `diag_anem`, so on every project carrying diagnostic
    records it cleared the anemometer slot and then cycled past the record that
    would have refilled it. The flags were read, the presence flag was set, and
    the data was gone.
    """

    def test_the_slug_the_interface_writes_is_accepted(self):
        self.assertIn(
            "case ('diag_anem', 'anemometer_diagnostic'); slot = diagAnem",
            read(E2SET))

    @unittest.skipUnless((GUI / "src/ecproject.cpp").exists(),
                         "eddyflow-gui not checked out beside this repository")
    def test_that_is_the_slug_the_interface_writes(self):
        self.assertIn('addPlain(g.diagColumns, QStringLiteral("diag_anem")',
                      read("src/ecproject.cpp", GUI))
        self.assertIn('rec.slug = QStringLiteral("diag_anem");',
                      read("src/basicsettingspage.cpp", GUI))


@unittest.skipUnless((GUI / "src/ecproject.cpp").exists(),
                     "eddyflow-gui not checked out beside this repository")
class TheInterfaceWritesWhatTheEngineReads(unittest.TestCase):
    """The two halves of the same file format."""

    def test_the_key_prefix_matches(self):
        source = read("src/ecproject.cpp", GUI)
        self.assertIn('writePlain(QStringLiteral("agc"), g.agcColumns);',
                      source)
        self.assertIn('readPlain(QStringLiteral("agc"), g.agcColumns);', source)

    def test_stale_agc_keys_are_cleared_before_writing(self):
        """QSettings preserves what it is not asked to overwrite, so a
        shrinking list would leave orphaned keys for the reader to pick up."""
        self.assertIn('QStringLiteral("^(gas|cell|diag|agc)_")',
                      read("src/ecproject.cpp", GUI))

    def test_the_slug_is_lower_case(self):
        """Lower case, unlike the display name, because the record is what the
        engine compares now and it compares it lower-cased. Deriving the slug
        from the display name keeps the two in step."""
        page = read("src/basicsettingspage.cpp", GUI)
        fn = unit(page, "void BasicSettingsPage::syncSignalStrengthRecords()",
                  "\n}\n")
        self.assertIn("rec.slug = name.toLower();", fn)
        for n in (35, 36):
            self.assertIn("VariableDesc::getVARIABLE_VAR_STRING_%d()" % n, fn)


if __name__ == "__main__":
    unittest.main()
