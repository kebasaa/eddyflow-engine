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
        self.assertGreaterEqual(source.count("call ApplyCovMaxDefaultFallback"), 2)
        self.assertIn("used_tlag = actual_tlag", source)
        self.assertIn("def_tlag_used = .true.", source)
        self.assertIn("used_tlag = dble(def_rl) / Metadata%ac_freq", source)
        self.assertNotIn("ActTLag(j) = TLag(j)\n                            RowLags(j) = def_rl(j)", source)

    def test_pwb_success_and_fallback_paths_preserve_status(self):
        source = read("src/src_rp/timelag_handle.f90")
        pwb_block = source[source.index("case ('pwb')") : source.index("case ('none')")]

        self.assertNotIn("call GetPwbFinalResult", pwb_block)
        self.assertIn("call PwbDetectGas", pwb_block)
        self.assertIn("lPwbResult%reliability_class = 'S1_optimal'", pwb_block)
        self.assertIn("lPwbResult%reliability_class = 'S2_optimal'", pwb_block)
        self.assertIn("lPwbResult%reliability_class = 'S3_carryforward'", pwb_block)
        self.assertIn("ActTLag(j) = lPwbResult%selected_lag", pwb_block)
        self.assertIn("DefTlagUsed(j) = .false.", pwb_block)
        self.assertIn("call ApplyCovMaxDefaultFallback", pwb_block)
        self.assertIn("lPwbResult%fallback_source = 'maxcov_default'", pwb_block)
        self.assertIn("lPwbResult%fallback_source = 'S3_carryforward'", pwb_block)

    def test_pwb_run_summary_is_printed_and_saved(self):
        module_source = read("src/src_rp/pwb_timelag_handle.f90")
        main_source = read("src/src_rp/eddyflow-rp_main.f90")
        globals_source = read("src/src_common/m_common_global_var.f90")

        self.assertIn("PwbSummary_FilePadding", globals_source)
        for token in ("PreparePwbBatch", "FinalizePwbBatch", "StorePwbRawResult", "GetPwbFinalResult"):
            self.assertNotIn(token, module_source)
            self.assertNotIn(token, main_source)
        self.assertIn("WARNING: all PWB detections fell back", module_source)
        self.assertIn("raw_selected_lag_s,raw_row_lag,applied_lag_s,applied_row_lag", module_source)
        self.assertIn("effective_block_length_s", module_source)
        self.assertIn("maxcov_default,nominal_default,other_fallback", module_source)
        self.assertIn("call ResetPwbDiagnostics()", main_source)
        self.assertIn("if (Meth%tlag == 'pwb') call ReportPwbDiagnostics()", main_source)

    def test_fluxnet_output_exposes_pwb_source_flags(self):
        header_source = read("src/src_rp/init_fluxnet_file_rp.f90")
        writer_source = read("src/src_rp/write_out_fluxnet.f90")

        for gas in ("CO2", "H2O", "CH4", "GS4"):
            self.assertIn(f"{gas}_TLAG_PWB_SOURCE", header_source)

        self.assertIn("0=native, 1=S3 carry-forward", writer_source)
        self.assertIn("Meth%tlag == 'pwb' .and. E2Col(gas)%present", writer_source)
        self.assertIn("PWBResult(gas)%fallback_source", writer_source)
        self.assertNotIn("median_raw", writer_source)


if __name__ == "__main__":
    unittest.main()
