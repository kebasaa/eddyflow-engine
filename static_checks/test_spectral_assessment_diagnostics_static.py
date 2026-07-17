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
