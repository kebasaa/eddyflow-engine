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

        #> The tallies are a run-log report now, not a file of their own.
        self.assertNotIn("PwbSummary_FilePadding", globals_source)
        for token in ("PreparePwbBatch", "FinalizePwbBatch", "StorePwbRawResult", "GetPwbFinalResult"):
            self.assertNotIn(token, module_source)
            self.assertNotIn(token, main_source)
        self.assertIn("WARNING: all PWB detections fell back", module_source)
        self.assertIn("raw_lag_s,raw_row_lag,hdi_low_s", module_source)
        self.assertIn("effective_block_length_s", module_source)
        self.assertIn("donor_gas", module_source)
        self.assertIn("call ResetPwbDiagnostics()", main_source)
        self.assertIn("if (Meth%tlag == 'pwb') call ReportPwbDiagnostics()", main_source)

    def test_fluxnet_output_exposes_pwb_source_flags(self):
        header_source = read("src/src_rp/init_fluxnet_file_rp.f90")
        writer_source = read("src/src_rp/write_out_fluxnet.f90")

        # The column is now generated per configured gas rather than spelled
        # out for the four historical slots, so the literal names are gone.
        # What must survive is that every gas still gets the column.
        self.assertIn(
            "trim(FluxnetLayoutTags(j)) // '_TLAG_PWB_SOURCE'", header_source
        )
        for gas in ("CO2", "H2O", "CH4", "GS4"):
            self.assertNotIn(f"{gas}_TLAG_PWB_SOURCE", header_source)

        self.assertIn("0=native, 1=S3 carry-forward, 2=instrument_shared", writer_source)
        self.assertIn("Meth%tlag == 'pwb' .and. E2Col(gas)%present", writer_source)
        self.assertIn("PWBResult(gas)%fallback_source", writer_source)
        self.assertNotIn("median_raw", writer_source)

    def test_pwb_per_period_cache_is_versioned_and_cache_aware(self):
        module_source = read("src/src_rp/pwb_timelag_handle.f90")
        handle_source = read("src/src_rp/timelag_handle.f90")

        #> Versions 1 and 2 carried a pre_wpl/post_wpl stage column for a
        #> choice that no longer exists, and predate the retired speed
        #> settings, so their fingerprint could not match this build anyway.
        self.assertIn("PWB_TIMELAG_CACHE_VERSION=4", module_source)
        for stale in (1, 2, 3):
            self.assertNotIn("PWB_TIMELAG_CACHE_VERSION=%d" % stale, module_source)
        self.assertIn("fingerprint=", module_source)
        self.assertIn("period_seconds=", module_source)
        self.assertIn("date,time,gas,", module_source)
        self.assertIn("fill_method", module_source)
        self.assertIn("subroutine ReadPwbTimelagCache", module_source)
        self.assertIn("subroutine WritePwbTimelagCache", module_source)
        self.assertIn("call LookupPwbTimelagCache", handle_source)
        self.assertIn("call StorePwbTimelagCache", handle_source)
        self.assertNotIn("cache_stage", handle_source)
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
        self.assertIn("PWB mode: exact per-period reuse", main_source)
        self.assertIn("'_pwb_timelag_opt'", writer_source)
        self.assertIn("PwbAggregateSummary", writer_source)

    def test_pwb_cache_generation_precedes_production_and_supports_assessment_only(self):
        parser_source = read("src/src_rp/read_ini_rp.f90")
        main_source = read("src/src_rp/eddyflow-rp_main.f90")

        self.assertIn("PwbCacheGenerate = .true.", parser_source)
        self.assertIn("PwbCacheUpdateRequested = len_trim(AuxFile%to) > 0", parser_source)
        self.assertIn("call WritePwbTimelagCache()", main_source)
        self.assertIn("PWB time-lag pre-pass finished.", main_source)
        #> The whole record is settled before it is written, which is the
        #> point of having a pre-pass at all.
        self.assertIn("call PostProcessPwbTimelagCache()", main_source)
        self.assertLess(main_source.index("call PostProcessPwbTimelagCache()"),
                        main_source.index("call WritePwbTimelagCache()"))
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

    def test_pwb_provenance_names_its_donor_from_the_record(self):
        """The donor is whichever gas lent the summary, not one of three.

        This was a `select case` over co2/h2o/ch4 with every other donor
        falling to a bare 'inferred', so a summary borrowed from a COS or from
        a second CO2 did not say where it came from - which is the one thing
        a provenance line exists to record.

        The buffer has to grow with the label. 'inferred_from_' is fourteen
        characters and GasLabel returns up to thirty-two, so character(32)
        would truncate a record-derived donor name, and a truncated provenance
        string still parses as a valid one.
        """
        writer_source = read("src/src_rp/writeout_timelag_optimization.f90")

        self.assertIn("character(64) :: source", writer_source)
        self.assertNotIn("character(16) :: source", writer_source)
        self.assertIn(
            "source = 'inferred_from_' // trim(GasLabel(PwbSummarySource(gas)))",
            writer_source,
            "the donor must be named from its own record",
        )
        for label in ("inferred_from_co2", "inferred_from_h2o", "inferred_from_ch4"):
            self.assertNotIn(
                f"source = '{label}'",
                writer_source,
                f"{label} is back as a literal; the three-case chain cannot "
                "name a donor outside it",
            )


class TheOptimisationSummaryNamesSpeciesNotSlots(unittest.TestCase):
    """One naming for a gas, and it comes from the record.

    `TimelagOptGasLabel` spelled slots five to eight co2/h2o/ch4/4th_gas by
    position and deferred to `GasLabel` only past the fourth. It was kept so an
    existing optimisation file would still match, and that cost the file its
    meaning: on a project whose records are ordered COS, CO2, H2O it headed the
    COS block 'co2', the CO2 block 'h2o', and wrote a 'ch4' block for a project
    measuring no methane.

    Writer and reader shared it, so a single project round-tripped and no flux
    moved. The damage is that the name reads as a species and means a slot: a
    user reordering their records between two runs gets a cached window
    restored onto the wrong gas, silently. `GasLabel` directly above it already
    carries this exact note for the PWB cache - fixing one file and not the
    other left the two disagreeing about what a gas is called.
    """

    HANDLE = "src/src_rp/pwb_timelag_handle.f90"
    WRITER = "src/src_rp/writeout_timelag_optimization.f90"
    READER = "src/src_rp/read_timelag_opt_file.f90"

    def code(self, path):
        """Source with comment-only lines dropped.

        The removal is documented in prose where the function stood, so a
        naive search matches the explanation rather than a live call.
        """
        return "\n".join(
            ln for ln in read(path).splitlines() if not ln.lstrip().startswith("!")
        )

    def test_the_positional_helper_is_gone(self):
        self.assertNotIn(
            "TimelagOptGasLabel",
            self.code(self.HANDLE),
            "a second naming for the same slots is what this removed",
        )

    def test_no_historical_spelling_is_pinned_to_a_slot(self):
        body = self.code(self.HANDLE)
        for slot in ("histGas1", "histGas2", "histGas3", "histGas4"):
            self.assertNotIn(
                f"case ({slot})",
                body,
                f"{slot} is being given a name again; a label selected by "
                "slot position says species and means position",
            )
        self.assertNotIn(
            "'4th_gas'",
            body,
            "'4th_gas' is the spelling that existed nowhere else in the tree",
        )

    def test_writer_and_reader_both_name_from_the_record(self):
        for path in (self.WRITER, self.READER):
            src = self.code(path)
            self.assertIn(
                "use m_pwb_timelag, only: GasLabel",
                src,
                f"{path} must take the record-derived label",
            )
            self.assertIn(
                "GasLabel(gas)",
                src,
                f"{path} must name each block from its own record",
            )

    def test_gaslabel_still_resolves_through_the_records(self):
        """The helper everything now depends on. If this stops asking the
        records, every file named from it goes back to meaning positions."""
        body = self.code(self.HANDLE)
        self.assertIn("call SpectralVarTags(tags)", body)
        self.assertIn("if (len_trim(tags(gas)) > 0) GasLabel = tags(gas)", body)


if __name__ == "__main__":
    unittest.main()
