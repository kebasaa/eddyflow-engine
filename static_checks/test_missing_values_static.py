"""What a file can say for "no reading", and where it has to be recognised.

Four things mean the same thing and used to be handled in four places or not at
all: an empty field, an unparseable token, a NaN, and a numeric fill. The last
two were the dangerous ones. A NaN reads as a perfectly good real and passes
every `/= error` test downstream. A fill was caught only by CleanUpE2Set's -300
floor, which runs on E2Set - after the unit conversion - so a -9999 in a
nmol mol-1 column arrived there as -9.999 and walked straight through.

And an unparseable token cost the whole record. Measured on the CH-LAE dataset:
1202 of 36001 rows in one file carry a quoted "NAN" in a MIRO gas column, and
every one of those rows was discarded entire - wind data included. The sonic
recovered 526 rows per averaging period when that stopped.
"""

import re
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def read(path):
    return (ROOT / path).read_text(encoding="utf-8", errors="replace")


class MissingValueStaticTests(unittest.TestCase):
    def test_missing_values_are_blanked_before_any_conversion(self):
        """A fill value is only recognisable while it still looks like itself."""
        body = read("src/src_common/define_all_var_set.f90")
        assert "call BlankMissingValues(" in body
        fn = body[body.index("subroutine DefineAllVarSet"):]
        fn = fn[:fn.index("end subroutine DefineAllVarSet")]
        assert (fn.index("call BlankMissingValues(")
                < fn.index("call ConvertTraceGasUnits(")), (
            "the missing-value pass runs after a conversion has already scaled "
            "the fill value out of recognition"
        )

    def test_nan_is_caught_by_something_other_than_equality(self):
        """A NaN compares equal to nothing, itself included, so an equality test
        against a sentinel never sees one."""
        body = read("src/src_common/blank_missing_values.f90")
        assert "Vec(1:N) /= Vec(1:N)" in body, (
            "NaN must be caught by self-inequality, not by comparing to a value"
        )
        assert "huge(" in body, "infinity is not caught"

    def test_the_declared_value_only_adds(self):
        """An absent declaration must mean "the built-ins alone", not "zero" -
        blanking every genuine zero would be catastrophic and silent."""
        body = read("src/src_common/blank_missing_values.f90")
        assert "if (err_value /= error) then" in body, (
            "an undeclared column falls through to comparing against 0"
        )
        reader = read("src/src_common/read_metadata_file.f90")
        assert "LocCol(i)%err_value = error" in reader, (
            "the reader must leave an absent key at the sentinel"
        )

    def test_one_bad_field_costs_one_field(self):
        """The fast whole-record read stays, because it carries millions of
        rows; only a record it cannot parse is walked field by field."""
        body = read("src/src_common/import_ascii.f90")
        assert "read(dataline, *, iostat = parse_status)" in body, (
            "the whole-record read is gone, and with it the fast path"
        )
        assert "if (parse_status /= 0) &" in body
        assert body.count("call ParseDataRecord(") == 2, (
            "both record loops - plain and data-label - need the fallback"
        )

    def test_both_ascii_importers_share_one_parser(self):
        """Two copies of this loop would answer differently for the same record
        the first time either was touched."""
        for path in ("src/src_common/import_ascii.f90",
                     "src/src_common/import_ascii_with_text.f90"):
            assert "call ParseDataRecord(" in read(path), path
        #> The inline copy must be gone from the text importer, not merely
        #> bypassed.
        body = read("src/src_common/import_ascii_with_text.f90")
        assert "il: do j = 1, NumCol" not in body, (
            "the text importer still carries its own copy of the field loop"
        )

    def test_the_parser_answers_per_raw_column(self):
        """ImportAscii's buffer is indexed by raw file column and compacted
        afterwards. Filling hot slots instead shifts every column past the first
        ignored one - which it did, and the flux numbers moved before the
        counts gave it away."""
        body = read("src/src_common/parse_data_record.f90")
        assert "Vec(j)" in body and "nhot" not in body, (
            "the parser is writing hot slots again"
        )
        text = read("src/src_common/import_ascii_with_text.f90")
        assert "RawRec" in text, (
            "the text importer must compact the raw-indexed record into its "
            "hot-indexed array"
        )

    def test_an_untested_gas_is_reported_once(self):
        body = read("src/src_rp/test_absolute_limits.f90")
        assert "call ExceptionHandler(109)" in body
        assert "if (.not. AlLimitsWarned) then" in body, (
            "the test runs every averaging period; unlatched, the warning "
            "would bury the log"
        )
        assert "case(109)" in read("src/src_common/exception_handler.f90")


if __name__ == "__main__":
    unittest.main()
