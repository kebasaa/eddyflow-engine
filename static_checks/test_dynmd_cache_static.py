"""The dynamic metadata file is read once, and each period still gets its row.

RetrieveDynamicMetadata used to open the file and scan it from the top for
every period - twice per period in the main pass. The rows are now held in
memory (m_dynmd_rows). What must not change:

* the row chosen: the scan stops at the first row dated after the period, and
  applies the row before it - a carried-on scan may resume only while periods
  do not go backwards;
* rows are parsed only as far as the old scan read, so a malformed row past
  every period the run needs is still never touched;
* applying the row still merges into what earlier rows left - Read, Fix and
  Extract run on every call;
* Warning(105) once per overflowing row scanned, as before;
* the old path stays for a file without a date or time column.

Part of the EddyFlow engine's static checks.
"""

from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[1]


def code(rel):
    text = (ROOT / rel).read_text(encoding="utf-8", errors="replace")
    return "\n".join(ln for ln in text.splitlines() if not ln.lstrip().startswith("!"))


SRC = code("src/src_rp/retrieve_dynamic_metadata.f90")


def body(kind, name):
    i = SRC.index(kind + " " + name + "(")
    return SRC[i:SRC.index("end " + kind + " " + name, i)]


class TheRowIsTheOneTheScanFound(unittest.TestCase):

    def test_a_scan_resumes_only_while_periods_do_not_go_back(self):
        pick = body("subroutine", "PickDynMDRow")
        self.assertIn("if (FinalTimestamp < LastFinal) HaveLast = .false.", pick)
        self.assertIn("StopRow = 1", pick)

    def test_it_stops_at_the_first_row_after_the_period_and_skips_blanks(self):
        pick = body("subroutine", "PickDynMDRow")
        self.assertIn("if (Rows(k)%blank) then", pick)
        self.assertIn("if (Rows(k)%ts <= FinalTimestamp) then", pick)

    def test_rows_are_parsed_only_when_a_scan_reaches_them(self):
        pick = body("subroutine", "PickDynMDRow")
        self.assertIn("if (.not. Rows(k)%parsed) call ParseRow(k, nvars)", pick)
        self.assertNotIn("call ParseRow", body("subroutine", "LoadDynMDRows"))

    def test_the_warning_is_raised_for_every_overflowing_row_scanned(self):
        pick = body("subroutine", "PickDynMDRow")
        self.assertIn("call ExceptionHandler(105)", pick)
        self.assertIn("if (Rows(StopRow)%overwide) w = w + 1", pick)


class ApplyingItIsUnchanged(unittest.TestCase):

    def test_read_fix_and_extract_run_on_every_call(self):
        retrieve = body("subroutine", "RetrieveDynamicMetadata")
        cached = retrieve[retrieve.index("if (DynRowsUsable) then"):retrieve.index("return")]
        for call in ("call ReadMetadataFromTextVars(", "call FixDynamicMetadata()",
                     "call ExtractUsableMetadataFromDynamic("):
            self.assertIn(call, cached)

    def test_the_old_reading_stays_for_files_it_cannot_hold(self):
        load = body("subroutine", "LoadDynMDRows")
        self.assertIn("DynamicMetadataOrder(dynmd_date) == nint(error)", load)
        self.assertIn("if (open_status /= 0) return", load)
        self.assertIn("open(udf, file = AuxFile%DynMD", body("subroutine", "RetrieveDynamicMetadata"))


if __name__ == "__main__":
    unittest.main()
