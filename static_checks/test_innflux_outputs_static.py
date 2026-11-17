"""Optional FLUXNET columns after innFLUX, each off by default.

innFLUX (Striednig et al. 2020, Atmos. Meas. Tech. 13, 1447-1465) reports the
potential and virtual potential temperature fluxes (MET.wtheta, MET.wtheta_v)
and, for every flux, the mean and SD of the cross-covariance far from the lag
and a signal-to-noise ratio (qaqc.flux_noise_mean/_std, flux_SNR; far-lag
noise after Wienhold et al. 1994 and Spirig et al. 2005).

  out_theta_flux       RP.SNTags(427)  W_THETA, W_THETA_V
  out_lod_noise_stats  RP.SNTags(428)  <VAR>_LOD_NOISE_MEAN/_SD, <VAR>_FLUX_SNR

Both go in the FLUXNET tail after the RFlux block - the part read_ex_record.f90
captures whole and FCC re-emits verbatim - and the header and the row must
add them in the same order, under the same switch.

Part of the EddyFlow engine's static checks.
"""

from pathlib import Path
import re
import unittest

ROOT = Path(__file__).resolve().parents[1]


def read(*parts):
    return ROOT.joinpath(*parts).read_text(encoding="utf-8", errors="replace")


def code(*parts):
    return "\n".join(ln for ln in read(*parts).splitlines()
                     if not ln.lstrip().startswith("!"))


TAGS = read("src", "src_rp", "m_rp_global_var.f90")
READER = code("src", "src_rp", "read_ini_rp.f90")
HDR = code("src", "src_rp", "init_fluxnet_file_rp.f90")
ROW = code("src", "src_rp", "write_out_fluxnet.f90")
THETA = read("src", "src_rp", "theta_fluxes.f90")
RU = code("src", "src_rp", "random_error_handle.f90")
NOISE = code("src", "src_rp", "langford_noise.f90")
MAIN = code("src", "src_rp", "eddyflow-rp_main.f90")


class TheKeys(unittest.TestCase):
    def test_slots(self):
        self.assertRegex(TAGS, r"SNTags\(427\)%Label\s*/\s*'out_theta_flux'")
        self.assertRegex(TAGS, r"SNTags\(428\)%Label\s*/\s*'out_lod_noise_stats'")

    def test_off_unless_stated(self):
        self.assertIn("RPSetup%out_theta_flux = .false.", READER)
        self.assertIn("RPSetup%out_lod_noise_stats = .false.", READER)
        self.assertIn("if (SNTagFound(427)) RPSetup%out_theta_flux", READER)
        self.assertIn("if (SNTagFound(428)) RPSetup%out_lod_noise_stats", READER)


class TheTail(unittest.TestCase):
    def test_after_the_rflux_block_on_both_sides(self):
        for src in (HDR, ROW):
            rf = src.index("if (Test%rf) then")
            th = src.index("if (RPSetup%out_theta_flux) then")
            ns = src.index("if (RPSetup%out_lod_noise_stats) then")
            self.assertLess(rf, th)
            self.assertLess(th, ns)

    def test_header_and_row_agree(self):
        self.assertIn("call AddDatum(csv_row, 'W_THETA', separator)", HDR)
        self.assertIn("call AddDatum(csv_row, 'W_THETA_V', separator)", HDR)
        self.assertIn("Essentials%wtheta, csv_row", ROW)
        self.assertIn("Essentials%wthetav, csv_row", ROW)
        fams = re.findall(r"call AddVariableFamily\('(_LOD_NOISE_MEAN|_LOD_NOISE_SD|_FLUX_SNR)'\)",
                          HDR)
        self.assertEqual(fams, ["_LOD_NOISE_MEAN", "_LOD_NOISE_SD", "_FLUX_SNR"])
        fields = re.findall(r"Essentials%(lod_noise_mean|lod_noise_sd|lod_noise_snr)\(var\)",
                            ROW)
        self.assertEqual(fields, ["lod_noise_mean", "lod_noise_sd", "lod_noise_snr"])
        # Same loop as the RFlux families: every row variable, which is the
        # set AddVariableFamily names.
        block = ROW[ROW.index("if (RPSetup%out_lod_noise_stats) then"):]
        self.assertEqual(block[:1200].count("do jg = 1, nRowVar"), 3)


class ThetaFluxes(unittest.TestCase):
    def test_called_after_level_zero(self):
        a = MAIN.index("call Fluxes0_rp(.true.)")
        b = MAIN.index("call ThetaFluxes()")
        self.assertLess(a, b)
        self.assertLess(b - a, 300)

    def test_innflux_constants(self):
        self.assertIn("p0_hpa = 1000d0", THETA)
        self.assertIn("kappa = 0.286d0", THETA)
        self.assertIn("Stats%Pr * 1d-2", THETA)

    def test_off_and_reset_every_period(self):
        body = THETA[THETA.index("subroutine ThetaFluxes"):]
        first = body.index("Essentials%wtheta = error")
        gate = body.index("if (.not. RPSetup%out_theta_flux) return")
        self.assertLess(first, gate)

    def test_credited(self):
        self.assertIn("Striednig et al. 2020", THETA)
        self.assertIn("Stull 1988", THETA)


class NoiseStats(unittest.TestCase):
    def test_reset_every_period_before_the_lod_gate(self):
        body = RU[RU.index("subroutine LodLangford15"):]
        reset = body.index("Essentials%lod_noise_snr = error")
        gate = body.index("if (RUsetup%lod_meth /= 'langford_15') return")
        self.assertLess(reset, gate)

    def test_snr_is_innflux_s(self):
        self.assertIn("abs(LangfordNoisePeak(j) - LangfordNoiseMean(j))", RU)
        self.assertIn("/ LangfordNoiseSd(j)", RU)

    def test_noise_kept_only_when_asked(self):
        self.assertIn("if (RPSetup%out_lod_noise_stats .and. var <= E2NumVar) then", NOISE)


if __name__ == "__main__":
    unittest.main()
