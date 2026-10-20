"""Limit of detection after Langford et al. (2015), wired once and read twice.

Langford, B., Acton, W., Ammann, C., Valach, A., and Nemitz, E.: Eddy-covariance
data with low signal-to-noise ratio: time-lag determination, uncertainties and
limit of detection, Atmos. Meas. Tech., 8, 4197-4213, 2015, Sect. 2.2.

The random flux error is the scatter of the w-scalar cross-covariance far from
the lag: RE_RMSE = sqrt(sigma^2 + mu^2) (their Eq. 9) or RE_SD = sigma, and
LoD = alpha RE. One routine, ``LangfordNoise``, measures it; two callers use it:

- ``ru_meth = 6`` (``langford_15``), the random error, on the detrended and
  aligned series the flux comes from, centred on lag zero;
- ``detlim_meth = 2`` (``langford_15``), the detection limit time-lag borrowing
  tests against, on the raw series before the shift, centred on each gas's lag.

Two call points rather than one shared result, deliberately: the borrowing test
has to run before the shift, and the random error has to see the detrending the
flux saw - Langford's Fig. 2c shows the LoD moving with it.

Known answers for the arithmetic are in ``src/src_tools/langford_check_main.f90``
(``mingw32-make langfordcheck``), run below when built.

Part of the EddyFlow engine's static checks.
"""

from pathlib import Path
import subprocess
import unittest


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
DRIVER = ROOT / "obj" / "win" / "langford_check.exe"


def read(*parts):
    return SRC.joinpath(*parts).read_text(encoding="utf-8", errors="replace")


def code(*parts):
    return "\n".join(ln for ln in read(*parts).splitlines()
                     if not ln.lstrip().startswith("!"))


DECODE = code("src_common", "write_processing_project_variables.f90")
READER = code("src_rp", "read_ini_rp.f90")
NOISE = code("src_rp", "langford_noise.f90")
RU = code("src_rp", "random_error_handle.f90")
DETLIM = code("src_rp", "flux_detection_limit.f90")
BORROW = code("src_rp", "borrow_timelag.f90")
CORE = code("src_common", "m_covmax_core.f90")
FULL_HDR_RP = code("src_rp", "init_outfiles_rp.f90")
FULL_RP = code("src_rp", "write_out_full.f90")
FULL_HDR_FCC = code("src_fcc", "init_out_files.f90")
FULL_FCC = code("src_fcc", "write_out_full_fcc.f90")
FLUXNET_HDR = code("src_rp", "init_fluxnet_file_rp.f90")
FLUXNET_RP = code("src_rp", "write_out_fluxnet.f90")
FLUXNET_FCC = code("src_fcc", "write_out_fluxnet_fcc.f90")
EX = code("src_common", "read_ex_record.f90")
GEN = (ROOT / "prj" / "gen_project_tags.py").read_text(encoding="utf-8")


class TheSettingsDecode(unittest.TestCase):

    def test_ru_meth_six_is_langford(self):
        self.assertIn("case(6)", DECODE)
        self.assertIn("RUsetup%meth = 'langford_15'", DECODE)

    def test_detlim_meth_two_is_langford(self):
        self.assertIn("RPSetup%detlim_meth = 'langford_15'", READER)

    def test_keys_sit_below_the_record_origins(self):
        for idx, key in ((65, "lod_window"), (66, "lod_stat"), (67, "lod_its_mult"),
                         (11, "lod_alpha")):
            self.assertIn('%d: "%s"' % (idx, key), GEN)

    def test_every_lod_read_is_guarded_with_a_literal_default(self):
        for n in (65, 66, 67):
            self.assertIn("if (SNTagFound(%d))" % n, READER)
        self.assertIn("RUsetup%lod_window = 0", READER)
        self.assertIn("RUsetup%lod_stat = 'rmse'", READER)
        self.assertIn("RUsetup%lod_its_mult = 10d0", READER)
        self.assertIn("RUsetup%lod_alpha = 1.96d0", DECODE)
        self.assertIn("if (EPPrjNTagFound(11))", DECODE)

    def test_the_scale_settings_are_decoded_without_a_random_error_method(self):
        # The automatic window of the detection limit needs them even when
        # ru_meth is none.
        head = DECODE[DECODE.index("RUsetup%its_meth = 'cross_0'") - 200:]
        self.assertNotIn("if (RUsetup%meth /= 'none') then\n        select case (nint(EPPrjNTags(23)",
                         DECODE)
        self.assertIn("RUsetup%tlag_max = 10", head)


class OneRoutineTwoCallers(unittest.TestCase):

    def test_the_random_error_is_measured_about_lag_zero(self):
        self.assertIn("case('langford_15')", RU)
        body = RU[RU.index("subroutine RU_Langford_15"):]
        self.assertIn("centre = 0", body)
        self.assertIn("call LangfordNoise(Set, N, M, centre, want, re)", body)
        for var in ("want(u) = .true.", "want(ts) = .true.",
                    "want(firstGas:lastGas) = .true."):
            self.assertIn(var, body)

    def test_the_detection_limit_is_measured_about_each_gas_lag(self):
        body = DETLIM[DETLIM.index("subroutine LangfordDetectionLimit"):]
        self.assertIn("centre(firstGas:lastGas) = RowLags(firstGas:lastGas)", body)
        self.assertIn("call LangfordNoise(Set, nrow, ncol, centre, want, re)", body)

    def test_langford_does_not_fall_through_to_wienhold(self):
        lang = DETLIM.index("if (RPSetup%detlim_meth == 'langford_15') then")
        wien = DETLIM.index("if (RPSetup%detlim_meth /= 'wienhold_94') return")
        self.assertLess(lang, wien)

    def test_borrowing_accepts_either_detection_limit(self):
        self.assertIn("RPSetup%detlim_meth /= 'langford_15') return", BORROW)

    def test_the_windows(self):
        for literal in ("150d0", "160d0", "180d0", "RUsetup%lod_its_mult",
                        "nint(10d0 * f)", "nint(30d0 * f)"):
            self.assertIn(literal, NOISE)
        self.assertIn("call FarLagCovarianceNoise(Set(:, w), Set(:, var), nrow, centre(var)",
                      NOISE)

    def test_the_core_is_pure_and_offers_both_statistics(self):
        self.assertIn("pure subroutine FarLagCovarianceNoise", CORE)
        self.assertIn("if (stat == 'sd') then", CORE)
        self.assertIn("re = sqrt(s2 / dble(nlags))", CORE)


class TheOutputs(unittest.TestCase):

    def test_full_output_columns_follow_each_random_error(self):
        for src in (FULL_HDR_RP, FULL_HDR_FCC):
            for col in ("'lod_Tau'", "'lod_H'", "'lod_LE'", "'lod_' // e2sg(gas)"):
                self.assertIn(col, src)
        for src in (FULL_RP, FULL_FCC):
            self.assertEqual(src.count("LimitOfDetection("), 4)

    def test_full_output_columns_only_with_langford(self):
        # Default runs keep their layout.
        for src in (FULL_HDR_RP, FULL_HDR_FCC, FULL_RP, FULL_FCC):
            self.assertEqual(src.count("if (RUsetup%meth == 'langford_15') then"), 4)

    def test_fluxnet_family_follows_the_random_errors_and_is_skipped_on_read(self):
        rand = FLUXNET_HDR.index("call AddFluxFamily('_RANDUNC_HF')")
        lod = FLUXNET_HDR.index("call AddFluxFamily('_LOD')")
        self.assertLess(rand, lod)
        for src in (FLUXNET_RP, FLUXNET_FCC):
            self.assertIn("lod_on = RUsetup%meth == 'langford_15'", src)
            self.assertEqual(src.count("LodField("), 6)
        self.assertIn("+ 4 + nExGas                &  !< skipped: Langford LoD", EX)
        self.assertIn("aux(1 : 4 + n_layout_gas), & !< Skip the LoD", EX)


class KnownAnswers(unittest.TestCase):

    @unittest.skipUnless(DRIVER.exists(), "langford_check not built (mingw32-make langfordcheck)")
    def test_the_driver_passes(self):
        r = subprocess.run([str(DRIVER)], capture_output=True, text=True, timeout=600)
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertIn("langford_check: all passed", r.stdout)


if __name__ == "__main__":
    unittest.main()
