"""A second hygrometer is a measurement, not a copy of the first.

A project may declare more than one H2O record. One of them is the site's -
PrimaryWaterSlot - and supplies the single latent heat flux, RH and
evapotranspiration a site has. The others are ordinary measurements that
happen to be water, and they are routed through the trace-gas path so that
their own flux is reported at all.

Two things went wrong with that routing, and neither is reachable with the
fixtures here. Both are pinned by construction and labelled as such.

**The molar-density factor.** MoleFractionsAndMixingRatios gives a trace gas
`d = chi/Va * 1d-3` and water `d = chi/Va`, because a trace gas's chi is on
the umol basis and water's is already on the mmol one. The open-path WPL arm
multiplied by a bare `1d3` to recover chi/Va - correct for a trace gas, and a
factor of a thousand too large for a second hygrometer, which reaches that arm
as a trace gas. No fixture has an open-path hygrometer, so this is arithmetic,
not observation.

**The flux the spectra are screened on.** Every water slot was judged on
`lEx%Flux0%LE`, the site's latent heat flux, so a second hygrometer's spectra
were accepted or rejected on the primary's flux. Each is now judged on its
own, converted to latent heat with lambda - which for the primary reproduces
Flux0%LE exactly, since that is how Flux0%LE is computed, so a
single-hygrometer site is unchanged.
"""

from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]


def code(rel):
    return "\n".join(ln for ln in (ROOT / rel).read_text(
        encoding="utf-8", errors="replace").splitlines()
        if not ln.lstrip().startswith("!"))


class TheDensityFactorFollowsTheSpecies(unittest.TestCase):
    def test_the_open_path_arm_does_not_hardcode_the_trace_gas_factor(self):
        src = code("src/src_rp/fluxes23_rp.f90")
        self.assertNotIn(
            "Stats%d(gas) * 1d3", src,
            "the bare 1d3 is the trace-gas conversion; water's molar density "
            "is already on the mmol basis, so a second hygrometer gets a WPL "
            "term a thousand times too large")
        self.assertIn("dens_to_chi", src)

    def test_water_takes_unity(self):
        src = code("src/src_rp/fluxes23_rp.f90")
        block = src[src.index("subroutine Level2GasFlux"):]
        block = block[:block.index("end subroutine Level2GasFlux")]
        self.assertIn("if (GasSlotIsWater(gas)) then", block)
        self.assertIn("dens_to_chi = 1d0", block)
        self.assertIn("dens_to_chi = 1d3", block)

    def test_the_two_conventions_are_still_what_this_assumes(self):
        """If MoleFractionsAndMixingRatios ever stops differing, this factor
        stops being needed - and would silently become wrong."""
        src = code("src/src_rp/molefractions_and_mixingratios.f90")
        self.assertIn("Stats%d(gas) = Stats%chi(gas) / LocVa(gas) * 1d-3", src,
                      "the trace-gas branch no longer applies 1d-3")
        self.assertIn("Stats%d(gas) = Stats%chi(gas) / LocVa(gas)\n", src + "\n",
                      "the water branch no longer omits it")


class EachHygrometerIsScreenedOnItsOwnFlux(unittest.TestCase):
    SITES = ("src/src_fcc/cospectra_qaqc.f90",
             "src/src_common/bpcf_fratini_12.f90")

    def test_the_site_flux_is_only_the_fallback(self):
        for rel in self.SITES:
            src = code(rel)
            self.assertIn("lEx%Flux0%gas(gas) * lEx%lambda", src,
                          "%s judges every water slot on the site's latent "
                          "heat flux" % rel)

    def test_the_primary_still_reduces_to_the_site_value(self):
        """Flux0%LE = Flux0%gas(wsl) * lambda * MW_H2O * 1d-3, so the
        conversion above is the same number for the primary. That is what
        keeps a single-hygrometer project unchanged."""
        src = code("src/src_rp/fluxes0_rp.f90")
        self.assertIn("Flux0%LE = Flux0%gas(wsl) * Ambient%lambda * MW_H2O * 1d-3",
                      src)


if __name__ == "__main__":
    unittest.main()
