from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]


def read(path):
    return (ROOT / path).read_text(encoding="utf-8")


def code(rel):
    """Source with full-line comments removed.

    A retirement leaves a note behind saying what was retired and why, so a
    plain `assertNotIn` against the raw source matches the explanation rather
    than the code and reports the removal as a failure.
    """
    return (chr(10)).join(ln for ln in read(rel).splitlines()
                          if not ln.lstrip().startswith("!"))


class PwbStaticIntegrationTests(unittest.TestCase):
    def test_method_id_five_maps_to_pwb_without_moving_existing_ids(self):
        source = read("src/src_rp/read_ini_rp.f90")
        expected_order = [
            "case ('0')",
            "case ('1')",
            "case ('2')",
            "case ('3')",
            "case ('4')",
            "case ('5')",
        ]
        positions = [source.index(token) for token in expected_order]
        self.assertEqual(positions, sorted(positions))
        self.assertIn("Meth%tlag = 'pwb'", source)

    def test_pwb_setup_tags_defaults_and_typedefs_exist(self):
        typedefs = read("src/src_common/m_typedef.f90")
        globals_ = read("src/src_rp/m_rp_global_var.f90")
        reader = read("src/src_rp/read_ini_rp.f90")
        self.assertIn("type :: PWBSetupType", typedefs)
        self.assertIn("type :: PWBResultType", typedefs)
        #> 406..413 are the flat per-gas lag bounds - co2, h2o, ch4 and the
        #> fourth gas, min and max. They are retired with the rest of the flat
        #> layer: the window comes from gas_<i>_pwb_min_lag now. The labels
        #> stay in the table because it is positional and a retired key is
        #> blanked, not removed, but nothing reads them.
        for tag in range(406, 414):
            self.assertIn(f"SNTags({tag})", globals_)
            self.assertNotIn(
                f"SNTags({tag})", reader,
                f"SNTags({tag}) is a flat per-gas lag bound; the record "
                "carries it now")
        #> 414..421 are whole-run scalars, not per-gas, and stay.
        for tag in range(414, 422):
            self.assertIn(f"SNTags({tag})", globals_)
            self.assertIn(f"SNTags({tag})", reader)
        for default in (
            "PWBSetup%n_bootstrap = 99",
            "PWBSetup%min_valid_frac = 0.3d0",
            "PWBSetup%hdi_thresh_s = 0.5d0",
            "PWBSetup%dev_thresh_s = 0.5d0",
            "PWBSetup%hdi_prefilter_s = 1.0d0",
            "PWBSetup%smoothing_width = 5",
            "PWBSetup%random_seed = 2024",
        ):
            self.assertIn(default, reader)

    def test_native_detector_and_diagnostics_are_wired_without_python_runtime(self):
        source = read("src/src_rp/pwb_timelag_handle.f90")
        #> The arithmetic lives in m_pwb_core, which the reference test drives
        #> directly; this module is the engine glue around it.
        core = read("src/src_common/m_pwb_core.f90")
        self.assertIn("module m_pwb_timelag", source)
        self.assertIn("use m_pwb_core", source)
        self.assertIn("subroutine PwbDetectGas", source)
        self.assertIn("FitArAic", core)
        self.assertIn("RunPwbCombination", source)
        self.assertIn("MapLagEstimate", core)
        self.assertIn("Hdi95", core)
        self.assertIn("edge_pinned", source)
        self.assertIn("fallback_used", source)
        self.assertIn("block_length_clamped", source)
        self.assertIn("effective_block_length_s", source)
        self.assertIn("EddyFlowProj%id(1:len_trim(EddyFlowProj%id))", source)
        self.assertIn("PwbTimelag_FilePadding", source)
        self.assertIn("Timestamp_FilePadding", source)
        globals_ = read("src/src_common/m_common_global_var.f90")
        self.assertIn("PwbTimelag_FilePadding  = '_pwb_timelag'", globals_)
        self.assertNotIn("import scipy", source.lower())

    def test_pwb_writes_two_files_and_not_four(self):
        """The half-hourly table and the aggregate summary, nothing else.

        _pwb_diagnostics repeated the per-period table one column apart, and
        _pwb_summary held per-gas tallies the run log already prints. Both are
        folded away, so a reader has one file to open per question.
        """
        globals_ = code("src/src_common/m_common_global_var.f90")
        self.assertNotIn("_pwb_diagnostics", globals_)
        self.assertNotIn("_pwb_summary", globals_)
        for rel in ("src/src_rp/pwb_timelag_handle.f90",
                    "src/src_rp/timelag_handle.f90",
                    "src/src_rp/eddyflow-rp_main.f90"):
            source = code(rel)
            self.assertNotIn("PwbTimelagDiag_FilePadding", source, rel)
            self.assertNotIn("PwbSummary_FilePadding", source, rel)
            self.assertNotIn("WritePwbDiagnostic", source, rel)

    def test_the_two_speed_options_are_gone(self):
        """Both cost accuracy for well under a percent of runtime.

        Skipping the CCF normalisation left the four pre-whitening
        combinations compared on unnormalised covariances carrying different
        physical units, so the winner was decided by unit scale rather than by
        peak prominence. Capping the AR order under-fits the pre-whitener,
        which is the step that sharpens the peak in the first place.
        """
        for rel in ("src/src_common/m_typedef.f90",
                    "src/src_rp/read_ini_rp.f90",
                    "src/src_rp/pwb_timelag_handle.f90"):
            source = code(rel)
            self.assertNotIn("PWBSetup%approx_ccf", source, rel)
            self.assertNotIn("PWBSetup%max_ar_order", source, rel)
        #> Blanked in the positional tag table, not removed from it - and 422
        #> has since been reused for pwb_max_carry_h, which is safe for the
        #> same reason retiring it was: the parser matches a tag by its label,
        #> so a project still stating pwb_approx_ccf finds nothing.
        globals_ = read("src/src_rp/m_rp_global_var.f90")
        for tag in (422, 423, 424):
            self.assertIn("SNTags(%d)" % tag, globals_)
        for tag in (423, 424):
            self.assertNotIn("SNTags(%d)%%Label / 'pwb_" % tag, globals_)
        self.assertIn("SNTags(422)%Label / 'pwb_max_carry_h' /", globals_)

    def test_detection_runs_pre_wpl_without_a_choice(self):
        """Both alternatives ran on rotated 20 Hz data.

        The only difference was whether the gas series had been through the
        pointwise mixing-ratio conversion, and that conversion runs before
        time-lag compensation - so after it, cell temperature and water sit in
        the gas series at the wrong relative lag, and the gas series is what
        is being cross-correlated.
        """
        for rel in ("src/src_common/m_typedef.f90",
                    "src/src_rp/read_ini_rp.f90",
                    "src/src_rp/eddyflow-rp_main.f90"):
            self.assertNotIn("detect_prewpl", code(rel), rel)
        main_source = code("src/src_rp/eddyflow-rp_main.f90")
        #> Detection sits between the rotation and the WPL conversion.
        rot = main_source.index("call TiltCorrection(Meth%rot")
        det = main_source.index("if (Meth%tlag == 'pwb') then", rot)
        wpl = main_source.index("call PointByPointToMixingRatio", rot)
        self.assertLess(rot, det)
        self.assertLess(det, wpl)

    def test_the_default_lag_window_covers_every_gas(self):
        """A gas with no window searches nothing and returns the default lag.

        Written as four scalar assignments to co2/h2o/ch4/gas4, every slot
        past the fourth kept whatever the loader left in PWBSetup, so a fifth
        gas entered the block-bootstrap with a zero-width bound. It is a
        whole-array assignment now, which also makes it independent of how
        many slots the legacy names happen to cover.
        """
        source = read("src/src_rp/read_ini_rp.f90")
        self.assertIn("PWBSetup%min_lag = -10d0", source)
        self.assertIn("PWBSetup%max_lag =  10d0", source)
        for slot in ("co2", "h2o", "ch4", "gas4"):
            self.assertNotIn(
                "PWBSetup%%min_lag(%s) = -10d0" % slot, source,
                "the default window is back to naming slots")

    def test_bounds_sensitive_loops_do_not_rely_on_short_circuiting(self):
        #> Fortran does not promise to stop evaluating .and., so an index test
        #> guarding a subscript has to be its own statement. The insertion
        #> sorts that need this moved into the core with everything else.
        core = read("src/src_common/m_pwb_core.f90")
        self.assertNotIn("j >= 1 .and. x(j)", core)
        self.assertIn("do while (j >= 1)", core)
        self.assertIn("if (x(j) <= tmp) exit", core)
        #> FillMissingLinear stayed behind: it is the one loop here that knows
        #> the engine's missing-value code, so it is engine glue rather than
        #> arithmetic.
        handler = read("src/src_rp/pwb_timelag_handle.f90")
        self.assertNotIn("k <= n .and. x(k)", handler)
        self.assertIn("if (k > n) exit", handler)
        self.assertIn("if (x(k) /= error) exit", handler)

    def test_timelag_handle_falls_back_and_makefile_references_source(self):
        handler = read("src/src_rp/timelag_handle.f90")
        main = read("src/src_rp/eddyflow-rp_main.f90")
        makefile = read("prj/Makefile")
        self.assertIn("case ('pwb')", handler)
        self.assertIn("call PwbDetectGas", handler)
        self.assertIn("call CovMax", handler)
        # fallback_used was not renamed, but its carrier was: lPwbResult is now
        # only the Pass-1 scratch copy, and the bookkeeping writes PWBResult(j).
        self.assertIn("PWBResult(j)%fallback_used = .true.", handler)

        # The single failure path became three labelled outcomes across the
        # multi-pass logic. Pin the labels, not just the boolean - a pass that
        # sets the flag without a source is what fallback_source was added for.
        self.assertIn("PWBResult(j)%fallback_used = .false.", handler)
        for label in ("'instrument_shared'", "'S3_carryforward'",
                      "'maxcov_default'", "'native'"):
            self.assertIn(f"fallback_source = {label}", handler)
        self.assertNotIn("call GetPwbFinalResult", handler)
        self.assertNotIn("pwb_prepass_loop", main)
        self.assertNotIn("PreparePwbBatch", main)
        self.assertIn("pwb_timelag_handle.o", makefile)


if __name__ == "__main__":
    unittest.main()
