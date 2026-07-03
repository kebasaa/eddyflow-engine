from pathlib import Path
import re
import unittest


ROOT = Path(__file__).resolve().parents[1]


def read(path):
    return (ROOT / path).read_text(encoding="utf-8")


class PwbTimelagStaticTests(unittest.TestCase):
    def test_timelag_handle_uses_shared_covmax_default_fallback(self):
        source = read("src/src_rp/timelag_handle.f90")

        self.assertIn("subroutine ApplyCovMaxDefaultFallback", source)
        self.assertEqual(source.count("call ApplyCovMaxDefaultFallback"), 2)
        self.assertIn("used_tlag = actual_tlag", source)
        self.assertIn("def_tlag_used = .true.", source)
        self.assertIn("used_tlag = dble(def_rl) / Metadata%ac_freq", source)
        self.assertNotIn("ActTLag(j) = TLag(j)\n                            RowLags(j) = def_rl(j)", source)

    def test_pwb_success_and_fallback_paths_preserve_status(self):
        source = read("src/src_rp/timelag_handle.f90")
        pwb_block = source[source.index("case ('pwb')") : source.index("case ('none')")]

        self.assertIn("if (pwb_success) then", pwb_block)
        self.assertIn("TLag(j) = lPwbResult%selected_lag", pwb_block)
        self.assertIn("DefTlagUsed(j) = .false.", pwb_block)
        self.assertIn("call ApplyCovMaxDefaultFallback", pwb_block)
        self.assertIn("PWBResult(j)%fallback_used = .true.", pwb_block)

    def test_pwb_run_summary_is_printed_and_saved(self):
        module_source = read("src/src_rp/pwb_timelag_handle.f90")
        main_source = read("src/src_rp/eddyflow-rp_main.f90")
        globals_source = read("src/src_common/m_common_global_var.f90")

        self.assertIn("PwbSummary_FilePadding", globals_source)
        self.assertIn("public :: PwbDetectGas, ResetPwbDiagnostics, ReportPwbDiagnostics", module_source)
        self.assertIn("WARNING: all PWB detections fell back", module_source)
        self.assertIn("gas,attempts,native_pwb,fallback,maxcov_default_fallback_only", module_source)
        self.assertIn("call ResetPwbDiagnostics()", main_source)
        self.assertIn("if (Meth%tlag == 'pwb') call ReportPwbDiagnostics()", main_source)

    def test_fluxnet_output_exposes_pwb_fallback_flags(self):
        header_source = read("src/src_rp/init_fluxnet_file_rp.f90")
        writer_source = read("src/src_rp/write_out_fluxnet.f90")

        for gas in ("CO2", "H2O", "CH4", "GS4"):
            self.assertIn(f"{gas}_TLAG_PWB_FALLBACK", header_source)

        self.assertIn("0 = native PWB, 1 = fallback", writer_source)
        self.assertIn("Meth%tlag == 'pwb' .and. E2Col(gas)%present", writer_source)
        self.assertIn("PWBResult(gas)%fallback_used", writer_source)


if __name__ == "__main__":
    unittest.main()
