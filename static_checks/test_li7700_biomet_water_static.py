"""A gas that takes its water from the biomet is corrected in every period.

Started with the LI-7700 multipliers; the same gap was in the dilution and the
drift correction's band broadening, covered at the end.

A gas whose moisture reference is the biomet (`gas_<i>_moist=-1`) has two
consumers of that humidity in RP: MoistTerms, which hands the WPL terms the
site's `Ambient%sigma` / `RHO%w`, and the LI-7700 spectroscopic multipliers
`Mul7700(gas)%A/B/C`. The multipliers read `Ambient%chi_biomet`.

The two agree whenever the period has a usable biomet RH: FluxParams then sets
RHO%w from that RH and forms chi_biomet from RHO%w. They part ways when it does
not - a gap in the biomet file, an RH out of range, or `biom_rh` naming a column
that is not an RH (embedded biomet selections count DATE and TIME as columns 1
and 2, so the 15th variable is `biom_rh=17`; 15 lands on P_RAIN). chi_biomet is
then error while RHO%w has fallen back to the hygrometer. The WPL went ahead,
the multipliers were skipped, SPEC_CORR_LI7700_*_CH4 read -9999, and
Level2GasFlux - whose guard exists because multiplying by the sentinel once
gave FCH4 around -8e8 - reported the methane flux missing.

On the GHG fixture with `gas_3_moist=-1, biom_rh=15` FCH4 was -9999 for both
half-hours; with the fix it is -18.537 / 75.827, the same as when CH4 names
the LI-7500A's water, and -18.555 / 75.862 with the biomet RH column selected
correctly (`biom_rh=17`, byte-identical before and after the fix).
"""

from pathlib import Path
import re
import unittest


ROOT = Path(__file__).resolve().parents[1]

RP_MAIN = "src/src_rp/eddyflow-rp_main.f90"
MOLEFRAC = "src/src_rp/molefractions_and_mixingratios.f90"
DRIFT = "src/src_rp/drift_correction.f90"
FLUX_PARAMS = "src/src_rp/flux_params.f90"
FLUXES23 = "src/src_rp/fluxes23_rp.f90"


def code(path):
    text = (ROOT / path).read_text(encoding="utf-8")
    return "\n".join(ln for ln in text.splitlines()
                     if not ln.lstrip().startswith("!"))


def multiplier_block():
    """From the per-period reset of Mul7700 to the call that fills it."""
    body = code(RP_MAIN)
    start = body.index("Mul7700(:)%A = error")
    end = body.index("call Multipliers7700(", start)
    return body, start, body[start:end]


def biomet_arm(block):
    m = re.search(r"if \(msl == biometMoistRef\) then\n(.*?)\n\s*elseif",
                  block, re.S)
    if m is None:
        raise AssertionError("the LI-7700 loop has no biomet arm")
    return m.group(1)


class TheMultipliersFollowMoistTerms(unittest.TestCase):
    FORMULA = "RHO%w * Ambient%Va / MW_H2O * 1d3"

    def test_the_biomet_arm_reads_the_site_density(self):
        arm = biomet_arm(multiplier_block()[2])
        self.assertIn("chi_moist = " + self.FORMULA, arm)
        self.assertNotIn("chi_biomet", arm,
                         "chi_biomet is error whenever the period has no "
                         "usable biomet RH, while the WPL of the same gas "
                         "goes ahead on RHO%w - the multipliers must read "
                         "what the WPL reads")

    def test_it_declines_only_when_the_site_has_no_water_either(self):
        arm = biomet_arm(multiplier_block()[2])
        self.assertIn("if (RHO%w == error .or. Ambient%Va == error) cycle",
                      arm)
        self.assertLess(arm.index("== error"), arm.index("chi_moist ="),
                        "the guard must come before the sentinel is "
                        "multiplied into a concentration")

    def test_it_is_the_same_arithmetic_as_chi_biomet(self):
        """So a period with a usable biomet RH is bit-identical to before."""
        self.assertIn("Ambient%chi_biomet = " + self.FORMULA,
                      code(FLUX_PARAMS))

    def test_moist_terms_hands_the_wpl_the_same_density(self):
        body = code(FLUXES23)
        start = body.index("subroutine MoistTerms(")
        terms = body[start:body.index("end subroutine MoistTerms", start)]
        self.assertIn("rhow_out  = RHO%w", terms)
        self.assertLess(terms.index("rhow_out  = RHO%w"),
                        terms.index("if (msl == biometMoistRef) return"))

    def test_rho_w_is_this_periods_when_the_multipliers_read_it(self):
        body, start, _ = multiplier_block()
        flux_params = body.rfind("call FluxParams(.true.)", 0, start)
        self.assertNotEqual(flux_params, -1,
                            "FluxParams sets RHO%w for the period and must "
                            "run before the multipliers read it")



class TheEarlierConsumersFallBackTheSameWay(unittest.TestCase):
    """The dilution and the drift broadening run before FluxParams, so RHO%w
    there is still the previous period's and cannot be what they read.
    FluxParams falls back to the primary hygrometer's chi - RHO%w is formed
    from Stats%chi(PrimaryWaterSlot()) - and these take that same chi
    directly. Without it the gas went undiluted, and a CO2 unbroadened, in
    exactly the periods its WPL and multipliers did go ahead."""

    def test_flux_params_falls_back_to_the_primary_hygrometer(self):
        self.assertIn("RHO%w = (Stats%chi(wsl) / Ambient%Va) * MW_H2O * 1d-3",
                      code(FLUX_PARAMS))

    def test_the_dilution_falls_back_to_it(self):
        body = code(MOLEFRAC)
        arm = biomet_arm(body)
        self.assertIn("waterChi = Ambient%chi_biomet", arm)
        self.assertIn("if (waterChi == error) then", arm)
        self.assertIn("wsl = PrimaryWaterSlot()", arm)
        self.assertIn("waterR   = Stats%r(wsl)", arm)
        self.assertIn("waterChi = Stats%chi(wsl)", arm)
        self.assertLess(arm.index("waterChi = Ambient%chi_biomet"),
                        arm.index("waterChi = Stats%chi(wsl)"),
                        "the biomet comes first; the hygrometer only when "
                        "the biomet has nothing this period")
        self.assertIn("integer, external :: PrimaryWaterSlot", body)

    def test_the_drift_broadening_falls_back_to_it(self):
        arm = biomet_arm(code(DRIFT))
        self.assertIn("chi_moist = Ambient%chi_biomet", arm)
        self.assertIn("if (chi_moist == error .and. wsl >= firstGas", arm)
        self.assertIn("chi_moist = Stats%chi(wsl)", arm)

    def test_the_drift_correction_reads_this_periods_biomet(self):
        """The main loop calls BiometWaterVapour after the drift correction,
        so the drift correction has to form it itself."""
        body = code(DRIFT)
        sub = body[body.index("subroutine DriftCorrection("):
                   body.index("end subroutine DriftCorrection")]
        call = sub.index("call BiometWaterVapour()")
        self.assertLess(sub.index("call AirAndCellParameters()"), call)
        self.assertLess(call, sub.index("call MoleFractionsAndMixingRatios()"))
        self.assertLess(call, sub.index("chi_moist = Ambient%chi_biomet"))


if __name__ == "__main__":
    unittest.main()
