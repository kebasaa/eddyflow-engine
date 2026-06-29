from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]


def read(path):
    return (ROOT / path).read_text(encoding="utf-8")


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
        for tag in range(406, 422):
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
        self.assertIn("module m_pwb_timelag", source)
        self.assertIn("subroutine PwbDetectGas", source)
        self.assertIn("FitArAic", source)
        self.assertIn("RunPwbCombination", source)
        self.assertIn("Hdi95", source)
        self.assertIn("edge_pinned", source)
        self.assertIn("fallback_used", source)
        self.assertIn("eddyflow_pwb_timelag_diagnostics.csv", source)
        self.assertNotIn("import scipy", source.lower())

    def test_bounds_sensitive_loops_do_not_rely_on_short_circuiting(self):
        source = read("src/src_rp/pwb_timelag_handle.f90")
        self.assertNotIn("j >= 1 .and. x(j)", source)
        self.assertNotIn("k <= n .and. x(k)", source)
        self.assertIn("do while (j >= 1)", source)
        self.assertIn("if (x(j) <= tmp) exit", source)
        self.assertIn("if (k > n) exit", source)
        self.assertIn("if (x(k) /= error) exit", source)

    def test_timelag_handle_falls_back_and_makefile_references_source(self):
        handler = read("src/src_rp/timelag_handle.f90")
        makefile = read("prj/Makefile")
        self.assertIn("case ('pwb')", handler)
        self.assertIn("call PwbDetectGas", handler)
        self.assertIn("call CovMax", handler)
        self.assertIn("PWBResult(j)%fallback_used = .true.", handler)
        self.assertIn("pwb_timelag_handle.o", makefile)


if __name__ == "__main__":
    unittest.main()
