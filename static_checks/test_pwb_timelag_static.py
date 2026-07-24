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
        self.assertIn("PWBResult(j)%reliability_class = 'S3_carryforward'", pwb_block)
        self.assertIn("DefTlagUsed(j) = .false.", pwb_block)
        self.assertIn("call ApplyCovMaxDefaultFallback", pwb_block)
        self.assertIn("PWBResult(j)%reliability_class = 'S4_instrument_shared'", pwb_block)
        self.assertIn("PWBResult(j)%fallback_source = 'instrument_shared'", pwb_block)
        self.assertIn("PWBResult(j)%donor_gas = GasLabel(k)", pwb_block)

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
        self.assertIn("donor_gas", module_source)
        self.assertIn("maxcov_default,nominal_default,other_fallback", module_source)
        self.assertIn("call ResetPwbDiagnostics()", main_source)
        self.assertIn("if (Meth%tlag == 'pwb') call ReportPwbDiagnostics()", main_source)

    def test_fluxnet_output_exposes_pwb_source_flags(self):
        header_source = read("src/src_rp/init_fluxnet_file_rp.f90")
        writer_source = read("src/src_rp/write_out_fluxnet.f90")

        for gas in ("CO2", "H2O", "CH4", "GS4"):
            self.assertIn(f"{gas}_TLAG_PWB_SOURCE", header_source)

        self.assertIn("0=native, 1=S3 carry-forward, 2=instrument_shared", writer_source)
        self.assertIn("Meth%tlag == 'pwb' .and. E2Col(gas)%present", writer_source)
        self.assertIn("PWBResult(gas)%fallback_source", writer_source)
        self.assertNotIn("median_raw", writer_source)

    def test_pwb_per_period_cache_is_versioned_and_cache_aware(self):
        module_source = read("src/src_rp/pwb_timelag_handle.f90")
        handle_source = read("src/src_rp/timelag_handle.f90")

        self.assertIn("PWB_TIMELAG_CACHE_VERSION=2", module_source)
        self.assertIn("PWB_TIMELAG_CACHE_VERSION=1", module_source)
        self.assertIn("fingerprint=", module_source)
        self.assertIn("period_seconds=", module_source)
        self.assertIn("date,time,gas,stage,actual_lag_s,used_lag_s", module_source)
        self.assertIn("subroutine ReadPwbTimelagCache", module_source)
        self.assertIn("subroutine WritePwbTimelagCache", module_source)
        self.assertIn("call LookupPwbTimelagCache", handle_source)
        self.assertIn("call StorePwbTimelagCache", handle_source)
        self.assertIn("cache_stage = 'pre_wpl'", handle_source)
        self.assertIn("cache_stage = 'post_wpl'", handle_source)
        self.assertIn("PWBResult = pwb_raw_Result", handle_source)
        self.assertIn("SetPwbPeriodTimestamp", module_source)
        self.assertIn("PwbPeriodDate", module_source)
        self.assertIn("PwbPeriodTime", module_source)
        self.assertIn("ValidPwbPeriodTimestamp", module_source)
        self.assertNotIn("PwbTimelagCache(i)%date == Stats%date", module_source)

    def test_pwb_aggregate_summary_uses_native_lags_and_donor_provenance(self):
        module_source = read("src/src_rp/pwb_timelag_handle.f90")
        handle_source = read("src/src_rp/timelag_handle.f90")
        writer_source = read("src/src_rp/writeout_timelag_optimization.f90")

        self.assertIn("AddPwbTimelagSummaryDataset", module_source)
        self.assertIn("S1_optimal", module_source)
        self.assertIn("S2_optimal", module_source)
        self.assertIn("PwbSummaryDonorCount", module_source)
        self.assertIn("PWBResult(j)%origin_gas", handle_source)
        self.assertIn("pwb_last_optimal_origin", handle_source)
        self.assertIn("PWB_aggregate_summary: true", writer_source)
        self.assertIn("PWB_summary_source_for_", writer_source)

    def test_pwb_supports_live_cache_and_aggregate_reuse_modes(self):
        parser_source = read("src/src_rp/read_ini_rp.f90")
        main_source = read("src/src_rp/eddyflow-rp_main.f90")
        writer_source = read("src/src_rp/writeout_timelag_optimization.f90")

        self.assertIn("PwbCacheUpdateRequested = len_trim(AuxFile%to) > 0", parser_source)
        self.assertIn("PWB mode: live detection during production processing.", main_source)
        self.assertIn("aggregate/RH-class time-lag reuse; PWB detection is disabled.", main_source)
        self.assertIn("PWB mode: exact per-period cache reuse", main_source)
        self.assertIn("'_pwb_timelag_opt'", writer_source)
        self.assertIn("PwbAggregateSummary", writer_source)

    def test_pwb_cache_generation_precedes_production_and_supports_assessment_only(self):
        parser_source = read("src/src_rp/read_ini_rp.f90")
        main_source = read("src/src_rp/eddyflow-rp_main.f90")

        self.assertIn("PwbCacheGenerate = .true.", parser_source)
        self.assertIn("PwbCacheUpdateRequested = len_trim(AuxFile%to) > 0", parser_source)
        self.assertIn("call WritePwbTimelagCache()", main_source)
        self.assertIn("PWB time-lag cache generation session terminated.", main_source)
        self.assertIn("PwbCacheGenerate .and. PWBSetup%detect_prewpl", main_source)
        self.assertIn("Meth%tlag == 'pwb' .and. PwbCacheDirty", main_source)
        self.assertIn("PwbTimelagN", main_source)
        self.assertIn("WriteOutTimelagOptimization(tlagn, E2NumVar, toH2On", main_source)

    def test_timelag_option_buffers_avoid_allocatable_size_warning_patterns(self):
        main_source = read("src/src_rp/eddyflow-rp_main.f90")

        self.assertIn("integer :: TimelagOptSize = 0", main_source)
        self.assertIn("integer :: PwbTimelagOptSize = 0", main_source)
        self.assertIn("PwbTimelagOptSize = size(RawTimeSeries) - 1", main_source)
        self.assertIn("allocate(PwbTimelagOpt(PwbTimelagOptSize))", main_source)
        self.assertIn("TimelagOptSize = toEndTimestampIndx - toStartTimestampIndx", main_source)
        self.assertIn("allocate(TimelagOpt(TimelagOptSize))", main_source)
        self.assertNotIn("AddPwbTimelagSummaryDataset(PwbTimelagOpt, size(PwbTimelagOpt)", main_source)
        self.assertNotIn("AddToTimelagOptDataset(TimelagOpt, size(TimelagOpt)", main_source)
        self.assertNotIn("FixTimelagOptDataset(TimelagOpt, size(TimelagOpt)", main_source)
        self.assertNotIn("PwbTimelagOpt(1:PwbTimelagN)", main_source)
        self.assertIn("FixTimelagOptDataset(PwbTimelagOpt, PwbTimelagOptSize", main_source)
        self.assertIn("PwbTimelagN > PwbTimelagOptSize", main_source)
        self.assertIn("ton > TimelagOptSize", main_source)

    def test_pwb_provenance_source_buffer_holds_inferred_labels(self):
        writer_source = read("src/src_rp/writeout_timelag_optimization.f90")

        self.assertIn("character(32) :: source", writer_source)
        self.assertNotIn("character(16) :: source", writer_source)
        for label in ("inferred_from_co2", "inferred_from_h2o", "inferred_from_ch4"):
            self.assertLessEqual(len(label), 32)
            self.assertIn(f"source = '{label}'", writer_source)


if __name__ == "__main__":
    unittest.main()
