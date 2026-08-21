"""Static regression checks for FCC spectral-assessment diagnostics."""

from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]


def read(path):
    return (ROOT / path).read_text(encoding="utf-8")


class SpectralAssessmentDiagnosticsStaticTests(unittest.TestCase):
    def test_main_collects_and_reports_diagnostics_before_output(self):
        source = read("src/src_fcc/eddyflow-fcc_main.f90")
        self.assertIn("call ResetSpectralAssessmentDiagnostics()", source)
        self.assertIn("SADiagSelectedFiles = SADiagSelectedFiles + 1", source)
        report = source.index("call ReportSpectralAssessmentDiagnostics")
        output = source.index("call OutputSpectralAssessmentResults")
        self.assertLess(report, output)

    def test_legacy_assessment_errors_and_fallback_are_preserved(self):
        source = read("src/src_fcc/eddyflow-fcc_main.f90")
        report = source.index("call ReportSpectralAssessmentDiagnostics")
        output = source.index("call OutputSpectralAssessmentResults(nbins)", report)
        self.assertLess(report, output)
        assessment_output = read("src/src_fcc/output_spectral_assessment_results.f90")
        self.assertIn("call ExceptionHandler(76)", assessment_output)
        self.assertIn("call ExceptionHandler(77)", assessment_output)
        correction = read("src/src_common/bpcf_bandpass_spectral_corrections.f90")
        self.assertIn("call ExceptionHandler(69)", correction)
        self.assertIn("actual_hf_method = 'moncrieff_97'", correction)

    def test_diagnostics_cover_required_prerequisites_and_filters(self):
        source = read("src/src_fcc/spectral_assessment_diagnostics.f90")
        for phrase in (
            "Full cospectra (Fratini et al. 2012)",
            "VM filtering",
            "Foken filtering",
            "H2O RH classes",
            "Valid degraded wT covariance",
            "Spectral assessment: SUCCESS",
            "Spectral assessment: FAILED",
        ):
            self.assertIn(phrase, source)

    def test_the_report_states_what_each_gas_will_actually_get(self):
        """PASS/FAIL says whether the assessment covered a gas; it does not say
        what that means for the flux. A gas marked FAIL is still corrected -
        analytically - and a run where half the gases are is not the same run
        as one where none are."""
        source = read("src/src_fcc/spectral_assessment_diagnostics.f90")
        self.assertIn("Spectral assessment: PARTIAL", source,
                      "a file covering some gases must not report as SUCCESS")
        self.assertIn("OutcomeLabel", source,
                      "each gas line must say in situ or analytically")
        self.assertIn("assessment_ready = n_insitu > 0", source,
                      "readiness is 'any gas fitted', not 'every gas fitted'")

    def test_the_suggestion_says_when_a_flux_limit_cannot_help(self):
        """A gas rejected for small fluxes is fixed by lowering the floor. A
        gas whose classes are too thinly populated is not, and the suggested
        floor - a percentile of what is already there - removes records when it
        lands above the current one. The two used to print identically."""
        source = read("src/src_fcc/spectral_assessment_diagnostics.f90")
        fn = source[source.index("subroutine ReportFluxLimitSuggestions"):]
        fn = fn[:fn.index("end subroutine ReportFluxLimitSuggestions")]
        self.assertIn("Not a flux-limit problem", fn)
        self.assertIn("sa_min_smpl", fn)
        self.assertIn("maxval(class_counts)", fn,
                      "say how far the best class actually got")
        self.assertIn("suggested_min > current_min", fn,
                      "a suggestion above the current floor discards records "
                      "and must be flagged as such")
        #> The auto-apply guard is what kept water's unusable 33.2 out of the
        #> project file. It is not part of the wording change and must stay.
        self.assertIn("suggested_min < current_min .and. valid_classes >= 1", fn)

    def test_qaqc_tracks_flux_vm_foken_and_accepted_records(self):
        source = read("src/src_fcc/cospectra_qaqc.f90")
        for counter in (
            "SADiagRejectedFlux",
            "SADiagRejectedUstar",
            "SADiagRejectedVM",
            "SADiagRejectedFoken",
            "SADiagAccepted",
        ):
            self.assertIn(counter, source)

    def test_legacy_assessment_errors_are_emitted(self):
        source = read("src/src_fcc/output_spectral_assessment_results.f90")
        self.assertIn("call ExceptionHandler(76)", source)
        self.assertIn("call ExceptionHandler(77)", source)

    def test_automatic_configuration_is_optional_and_updates_only_output_project(self):
        globals_source = read("src/src_fcc/m_fx_global_var_mod.f90")
        parser = read("src/src_fcc/read_ini_fcc.f90")
        diagnostics = read("src/src_fcc/spectral_assessment_diagnostics.f90")
        main = read("src/src_fcc/eddyflow-fcc_main.f90")

        self.assertIn("automatic_spectra_config", globals_source)
        self.assertIn("if (SCTagFound(27))", parser)
        self.assertIn("FCCsetup%SA%automatic_config = .false.", parser)
        self.assertIn("subroutine ApplyAutomaticSpectralConfiguration(output_project)", diagnostics)
        self.assertIn("call EditIniFile(trim(output_project)", diagnostics)
        self.assertIn("'automatic_spectra_config', '0'", diagnostics)
        copy = main.index("call CopyFile(trim(adjustl(PrjPath))")
        automatic = main.index("call ApplyAutomaticSpectralConfiguration")
        self.assertLess(copy, automatic)

    def test_automatic_flux_recommendations_require_usable_class_coverage(self):
        source = read("src/src_fcc/spectral_assessment_diagnostics.f90")
        self.assertIn("suggested_min < current_min .and. valid_classes >= 1", source)
        self.assertIn("valid_classes > current_valid_classes", source)
        self.assertIn("'(f0.6)'", source)

    def test_assessment_only_mode_runs_requested_auxiliary_work_then_exits(self):
        tags = read("src/src_rp/m_rp_global_var.f90")
        parser = read("src/src_rp/read_ini_rp.f90")
        main = read("src/src_rp/eddyflow-rp_main.f90")

        self.assertIn("rot_pf_assessment_only", tags)
        self.assertIn("tlag_assessment_only", tags)
        self.assertIn("SCTagFound(100) .and.", parser)
        self.assertIn("SCTagFound(101) .and.", parser)
        self.assertIn("RPsetup%pf_assessment_only = RPsetup%pf_assessment_only .and.", parser)
        self.assertIn("RPsetup%tlag_assessment_only = RPsetup%tlag_assessment_only .and.", parser)
        self.assertLess(main.index("TIME LAG OPTIMIZATION IF REQUESTED"),
                        main.index("PLANAR FIT IF REQUESTED"))
        self.assertIn(".not. AssessmentOnly .or. RPsetup%tlag_assessment_only", main)
        self.assertIn(".not. AssessmentOnly .or. RPsetup%pf_assessment_only", main)
        self.assertLess(main.index("Auxiliary assessment-only session completed."),
                        main.index("Create TimeSeries for actual raw data processing"))


if __name__ == "__main__":
    unittest.main()


class WaterIsResolvedInTheAssessment(unittest.TestCase):
    """The RH-class assessment and the flux-limit thresholds are water's, and
    both read the sixth slot - water only when record two holds it."""

    def test_the_diagnostics_resolve_the_water_record(self):
        src = read("src/src_fcc/spectral_assessment_diagnostics.f90")
        self.assertIn("PrimaryWaterOutSlot()", src)
        self.assertNotIn("MeanBinSpecAvailable(cls, h2o)", src)
        self.assertNotIn("RegPar(h2o, cls)", src)

    def test_the_results_writer_resolves_it_too(self):
        src = read("src/src_fcc/output_spectral_assessment_results.f90")
        self.assertIn("wsl = PrimaryWaterOutSlot()", src)
        self.assertNotIn("%fn(h2o)", src)
        self.assertNotIn("RegPar(h2o,", src)

    def test_the_rh_class_table_is_one_loop(self):
        """Nine copies of the same three-line write, each naming its class in
        a literal that could disagree with the class it reported."""
        src = read("src/src_fcc/output_spectral_assessment_results.f90")
        self.assertNotIn("RH class   5 - 15% = ", src)
        self.assertIn("10 * cls - 5", src)

    def test_every_gas_reads_its_own_flux_limits(self):
        """The four arms ended in a `case default` reading gas4's thresholds,
        so a fifth gas was filtered on the fourth slot's limits even though
        ReadIniFCC had loaded its own from the record. A water record past the
        fourth slot was filtered on trace-gas limits rather than on LE."""
        src = read("src/src_fcc/spectral_assessment_diagnostics.f90")
        self.assertIn("FCCsetup%SA%min_un_gas(gas)", src)
        self.assertIn("FCCsetup%SA%max_gas(gas)", src)
        self.assertNotIn("FCCsetup%SA%min_un_gas(gas4)", src)
        self.assertIn("GasSlotIsWater(gas)", src,
                      "the latent-heat arm is a species question")
