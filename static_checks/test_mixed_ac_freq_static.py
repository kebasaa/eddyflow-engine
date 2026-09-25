"""Guards GHG projects whose files are not all at one acquisition frequency.

Each GHG archive carries its own .metadata, and a site can switch from 10 to
20 Hz part way through a project, or an analyser in it from 1 to 0.33 Hz. The
buffer a period is read into, and every record count checked against it, used
to be sized once from the first file. After a switch, half of each 20 Hz period
did not fit, and every period was dropped as "not enough samples". The
assessment pre-passes had no such check and silently used the half that fit.

None of this shows in a single-rate regression run, where every new path
reduces to the old one. What has to stay true:

  1. a file's rate is compared with the period's BEFORE its data are read,
     because the record window and the buffer were both sized for the period;
  2. every import goes through one wrapper, which resizes and reads the period
     again when its first file is at a new rate - the pre-passes included;
  3. a period whose files disagree, on the file rate or on any used column's
     instrument rate, is skipped with a warning rather than joined;
  4. the binned-spectra grid is built for the highest rate the survey finds,
     and FCC's Nyquist is the highest of any record, not the first one's.
"""

import re
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def read(path):
    return (ROOT / path).read_text(encoding="utf-8", errors="replace")


MAIN = read("src/src_rp/eddyflow-rp_main.f90")
IMPORT = read("src/src_rp/import_current_period.f90")
GHG = read("src/src_rp/read_licor_ghg_archive.f90")
SURVEY = read("src/src_rp/survey_ghg_ac_freq.f90")
FCC_INIT = read("src/src_fcc/init_ex_vars.f90")
DYNMD = read("src/src_rp/retrieve_dynamic_metadata.f90")


class TheRateIsCheckedBeforeTheData(unittest.TestCase):
    def test_the_reader_compares_before_importing(self):
        check = GHG.index("if (ExpectedAcFreq > 0d0 .and. DataIsNeeded) then")
        meta = GHG.index("call ReadMetadataFile(LocCol, MetaFile")
        data = GHG.index("call ImportNativeData(DataFile")
        self.assertLess(meta, check)
        self.assertLess(check, data)

    def test_the_record_window_uses_the_period_rate(self):
        """A file skipped for its rate leaves that rate in Metadata."""
        fix = IMPORT.index("if (ExpectedAcFreq > 0d0) Metadata%ac_freq = ExpectedAcFreq")
        portion = IMPORT.index("call PortionOfFileInCurrentPeriod(")
        self.assertLess(fix, portion)

    def test_the_metadata_retriever_sees_every_file(self):
        self.assertIn("ExpectedAcFreq = -1d0", IMPORT)


class OneWrapperForEveryImport(unittest.TestCase):
    def test_only_the_wrapper_calls_the_importer(self):
        self.assertEqual(MAIN.count("call ImportCurrentPeriod("), 1)
        wrapper = MAIN[MAIN.index("subroutine ImportPeriod("):]
        self.assertIn("call ImportCurrentPeriod(", wrapper)
        self.assertEqual(MAIN.count("call ImportPeriod("), 4,
                         "time lag, planar fit, drift and the main loop")

    def test_a_new_rate_resizes_and_reads_again(self):
        wrapper = MAIN[MAIN.index("subroutine ImportPeriod("):]
        self.assertIn("if (.not. RateChanged) exit", wrapper)
        self.assertIn("call SetPeriodRate(Metadata%ac_freq)", wrapper)

    def test_the_buffer_follows_the_rate(self):
        body = MAIN[MAIN.index("subroutine SetPeriodRate("):
                    MAIN.index("end subroutine SetPeriodRate")]
        self.assertIn("if (size(Raw, 1) /= MaxPeriodNumRecords) then", body)
        self.assertIn("allocate(Raw(MaxPeriodNumRecords, ncolRaw))", body)

    def test_nothing_else_recomputes_the_record_limits(self):
        self.assertEqual(MAIN.count("MaxPeriodNumRecords = "), 1)
        self.assertEqual(MAIN.count("MaxNumFileRecords   = "), 1)


class AMixedPeriodIsSkipped(unittest.TestCase):
    def test_a_later_file_at_another_rate(self):
        block = IMPORT[IMPORT.index("if (rate_mismatch) then"):]
        block = block[:block.index("return")]
        self.assertIn("rate_changed = .true.", block)
        self.assertIn("call ExceptionHandler(116)", block)
        self.assertIn("NextFile = CurrentFile", block)

    def test_an_instrument_changing_rate(self):
        """The file rate can stay while an analyser in it changes."""
        self.assertIn("PeriodInstrFreq = LocCol%instr%ac_freq", IMPORT)
        self.assertIn("if (.not. LocCol(j)%useit) cycle", IMPORT)

    def test_dynamic_metadata_does_not_override_a_ghg_file(self):
        body = DYNMD[DYNMD.index("subroutine ExtractUsableMetadataFromDynamic"):]
        guard = body.index("if (EddyFlowProj%ftype == 'licor_ghg' .and. .not. EddyFlowProj%use_extmd_file) then")
        override = body.index("Metadata%ac_freq = DynamicMetadata%ac_freq")
        self.assertLess(guard, override)


class TheSpectraGridCoversTheHighestRate(unittest.TestCase):
    def test_the_grid_is_built_from_the_survey(self):
        self.assertIn("BinGridAcFreq = max(SurveyAcFreq, Metadata%ac_freq)", MAIN)
        self.assertRegex(MAIN, r"call BinnedFrequencyVector\(bf, Meth%spec%nbins, &\s*"
                               r"RPsetup%avrg_len, BinGridAcFreq\)")

    def test_the_survey_leaves_the_globals_alone(self):
        peek = SURVEY[SURVEY.index("subroutine PeekGhgAcFreq("):]
        self.assertNotIn("ReadMetadataFile", peek)
        self.assertNotIn("Metadata%", peek)

    def test_the_survey_is_silent_for_one_rate(self):
        """A single-rate project's log must not change."""
        self.assertIn("if (nchanges > 0) then", SURVEY)

    def test_each_change_is_narrowed_to_the_file(self):
        """The list names when each rate starts, which a sparse sample alone
        cannot say: two files read that disagree are halved until they are
        neighbours."""
        self.assertIn("if (i > a + 1 .and. Differ(rate(a), rate(i))) then", SURVEY)
        self.assertIn("m = UntriedNear((a + i) / 2, a, i)", SURVEY)

    def test_the_change_is_dated_by_its_data_start(self):
        self.assertIn("if (EddyFlowLog%tstamp_end) tsFrom = tsFrom - DatafileDateStep",
                      SURVEY)

    def test_fcc_takes_the_highest_record_rate(self):
        self.assertIn("if (lEx%ac_freq > FCCMetadata%ac_freq) &", FCC_INIT)


if __name__ == "__main__":
    unittest.main()
