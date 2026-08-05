"""The automatic spectral configuration writes keys the engine reads.

That feature is not a report. It reads the flux distribution, decides better
per-gas thresholds, and writes them with EditIniFile into a copy of the
project that the next run reads. So the key it composes has to be a key
ReadIniFCC looks for, and nothing checked that.

It composed `sa_min_un_co2`, `sa_min_st_ch4`, `sa_max_gas4` for the first
four slots and `sa_min_un_gas_5_record` beyond them. None of those exist: the
flat `sa_min_*_<gas>` tags were retired with the record format, and the
spelling ReadIniFCC reads is `gas_<i>_sa_min_un`. The whole feature was inert
- for every gas, including the four it had names for - and it reported changes
that never took effect.

This is checked against the generated tag table rather than against a fixture.
The suggestion path needs more data than the regression window holds: on the
three hours available the feature correctly reports "no qualifying
recommendations", so a fixture can show the run completes but not that the
key it would write is readable. Comparing writer and table is the stronger
claim anyway - it is the contract that was broken.
"""

from pathlib import Path
import re
import unittest


ROOT = Path(__file__).resolve().parents[1]

WRITER = "src/src_fcc/spectral_assessment_diagnostics.f90"
READER = "src/src_fcc/read_ini_fcc.f90"
TABLE = "src/src_fcc/m_fx_global_var_mod.f90"

#: The three limits the feature can rewrite.
SUFFIXES = ("min_un", "min_st", "max")


def code(rel):
    return "\n".join(ln for ln in (ROOT / rel).read_text(
        encoding="utf-8", errors="replace").splitlines()
        if not ln.lstrip().startswith("!"))


class TheKeyIsARecordKey(unittest.TestCase):
    def test_the_writer_composes_the_record_spelling(self):
        src = code(WRITER)
        self.assertIn("'gas_', gas_slot - firstGas + 1, '_sa_'", src,
                      "the key must be gas_<i>_sa_<suffix>, which is what "
                      "ReadIniFCC reads")

    def test_the_retired_spellings_are_not_composed(self):
        src = code(WRITER)
        for gone in ("'sa_min_un_' //", "'sa_min_st_' //", "'sa_max_' //"):
            self.assertNotIn(
                gone, src,
                "%s builds a flat key that was retired with the record "
                "format, so the setting it writes is never read" % gone)
        self.assertNotIn(
            "_record'", src,
            "the `gas_<i>_record` suffix names no tag at all")

    def test_every_key_it_can_write_exists_in_the_tag_table(self):
        """The contract that was broken, stated directly."""
        table = code(TABLE)
        for suffix in SUFFIXES:
            self.assertIn(
                "'gas_1_sa_%s'" % suffix, table,
                "the writer can emit gas_<i>_sa_%s but the FCC tag table has "
                "no such label, so the engine will never read it" % suffix)

    def test_the_reader_reads_them_from_the_record_origin(self):
        """And that the reader has not drifted the other way."""
        src = code(READER)
        self.assertIn("fccGasOriginN + (gas - 1) * fccGasLeapN", src)


class TheReportNamesTheSameKey(unittest.TestCase):
    """A suggestion the user cannot act on is worse than none.

    The report prints the key beside the value it suggests. If that key and
    the key written to the project ever diverge, the text tells the user to
    edit a setting that does nothing.
    """

    def test_the_report_and_the_writer_share_one_label(self):
        src = code(WRITER)
        self.assertIn("trim(min_label)", src)
        self.assertIn("trim(max_label)", src)
        self.assertIn("WriteAutomaticSpectralSetting(output_project, trim(min_label)",
                      src.replace("\n", " ").replace("  ", " "),
                      "the label shown in the report must be the one written "
                      "to the project")


if __name__ == "__main__":
    unittest.main()
