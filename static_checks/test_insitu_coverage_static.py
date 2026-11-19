"""A gas that borrows a tube-mate's transfer function is reported, not warned.

Warning(110) said "the spectral assessment succeeded for some gases and not for
others" whenever a configured gas had no fit of its own. A gas without a fit
that takes another gas's transfer function on the same analyser - named in
gas_<i>_sa_source, or chosen by automatic selection - is corrected in situ all
the same, so it was warned about for nothing; and when the assessment came from
a file that was read, "succeeded" and "the file was written" were not true of
anything. The reader also raised Alert(65), "switched to Moncrieff et al.
(1997)", for any gas without a block, which nothing switched.

Now ReportInsituCoverage says it per gas, in both places: an Info> line naming
the donor for a gas that borrows, and a warning only for a gas left with no
in-situ correction at all - Warning(110) for an assessment made in the run,
Warning(111) for a file that was read. Which fit corrects a gas is asked of
SpectralSourceSlot, the same rule the flux correction follows.
"""

import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def code(rel):
    """Source without comment-only lines."""
    text = (ROOT / rel).read_text(encoding="utf-8", errors="replace")
    return "\n".join(l for l in text.splitlines() if not l.strip().startswith("!"))


HELPER = code("src/src_common/gas_slot_resolution.f90")
DIAG = code("src/src_fcc/spectral_assessment_diagnostics.f90")
READER = code("src/src_fcc/read_spectral_assessment_file.f90")
WRITER = code("src/src_fcc/output_spectral_assessment_results.f90")


def body(src, start, end):
    a = src.index(start)
    return src[a:src.index(end, a)]


class TheDonorIsTheCorrectionsOwn(unittest.TestCase):

    def test_asked_of_the_rule_the_flux_correction_follows(self):
        fn = body(HELPER, "integer function SpectralInsituDonor",
                  "end function SpectralInsituDonor")
        self.assertIn("src = SpectralSourceSlot(gas_slot, month, SA)", fn)
        self.assertIn("if (.not. GasHasSpectralFit(src)) cycle", fn)

    def test_water_never_borrows(self):
        fn = body(HELPER, "integer function SpectralInsituDonor",
                  "end function SpectralInsituDonor")
        self.assertIn("if (GasSlotIsWater(gas_slot)) return", fn)


class OneRoutineSaysItPerGas(unittest.TestCase):

    def setUp(self):
        self.fn = body(DIAG, "subroutine ReportInsituCoverage",
                       "end subroutine ReportInsituCoverage")

    def test_a_borrowing_gas_gets_an_info_line_and_is_not_counted(self):
        borrow = self.fn[self.fn.index("if (donor > 0) then"):]
        borrow = borrow[:borrow.index("else")]
        self.assertIn("call LogSay(' Info> '", borrow)
        self.assertIn("n_covered = n_covered + 1", borrow)
        self.assertIn("has no spectral transfer function of its own; corrected with", borrow)
        self.assertIn("'set by user'", self.fn)

    def test_the_warning_depends_on_where_the_assessment_came_from(self):
        self.assertIn("if (n_covered == 0 .or. n_none == 0) return", self.fn)
        branch = self.fn[self.fn.index("if (from_file) then"):]
        self.assertLess(branch.index("call ExceptionHandler(111)"),
                        branch.index("call ExceptionHandler(110)"))

    def test_the_gases_left_out_are_named(self):
        self.assertIn("Corrected with the analytical method instead: ", self.fn)


class BothPathsUseIt(unittest.TestCase):

    def test_the_reader(self):
        self.assertIn("call ReportInsituCoverage(.true., own_fit)", READER)
        self.assertNotIn("ExceptionHandler(110)", READER)

    def test_the_writer_for_one_rate_and_for_several(self):
        self.assertEqual(WRITER.count("call ReportInsituCoverage(.false., own_fit)"), 2)
        self.assertNotIn("ExceptionHandler(110)", WRITER)

    def test_alert_65_is_for_the_malformed_file_only(self):
        self.assertEqual(READER.count("call ExceptionHandler(65)"), 1)
        malformed = READER[READER.index("if (read_status /= 0) then"):]
        malformed = malformed[:malformed.index("end if")]
        self.assertIn("call ExceptionHandler(65)", malformed)
        self.assertIn("EddyFlowProj%hf_meth = 'moncrieff_97'", malformed)


if __name__ == "__main__":
    unittest.main()
