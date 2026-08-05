"""A gas is corrected with the water it names, not with the site's.

`E2Col(gas)%moist_ref` is how a gas says which hygrometer corrects it, and
ResolveGasRef fills it in preferring one on the gas's own analyser. Several
corrections ignored it and read `PrimaryWaterSlot()` instead, which is the
site's humidity - one answer for a site that may have two analysers sampling
different air.

Worst of them was the water-flux term of the WPL correction. TimeLagHandle
computed Cov(w, water) only when the *primary* hygrometer was closed-path and
only for gases sharing the *primary's* instrument model, so a gas on a second
analyser got no covariance and therefore no water-flux term - while its sigma
and rho_w came from that second analyser's hygrometer. The two halves of one
term disagreed about which water they meant.

There is a real distinction underneath, and these checks encode it:

  time-series operations - the covariance, and point-by-point dilution - use
    the gas's reference AND require it to be on the same analyser. They work
    on raw columns at the gas's own row lag, and a hygrometer down another
    tube has a different lag, so the pairing means nothing.

  mean/WPL terms - sigma and rho_w - use the reference and accept
    ResolveGasRef's fallback to "the first hygrometer anywhere". Ambient
    humidity is a reasonable stand-in there; MoistTerms starts from the site
    value and overrides only where the per-gas one is real.

What must not come back is either of them reading the site's water directly.
"""

from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]


def code(rel):
    return "\n".join(ln for ln in (ROOT / rel).read_text(
        encoding="utf-8", errors="replace").splitlines()
        if not ln.lstrip().startswith("!"))


class EachGasUsesItsOwnReference(unittest.TestCase):
    SITES = {
        "src/src_rp/timelag_handle.f90": "the water-flux covariance",
        "src/src_rp/fluxes0_rp.f90": "the E_gas divisor",
        "src/src_rp/drift_correction.f90": "the CO2 band broadening",
        "src/src_rp/eddyflow-rp_main.f90": "the LI-7700 spectroscopic term",
    }

    def test_each_reads_moist_ref(self):
        for rel, what in self.SITES.items():
            self.assertIn("moist_ref", code(rel),
                          "%s must take %s from the gas's own record"
                          % (rel, what))

    def test_the_covariance_is_not_gated_on_the_primary(self):
        """The gate that cost a second analyser's gases their WPL term."""
        src = code("src/src_rp/timelag_handle.f90")
        self.assertNotIn(
            "E2Col(j)%instr%model /= E2Col(wsl)%instr%model", src,
            "gases are being matched against the primary hygrometer's "
            "analyser again, so a gas on a second one gets no covariance")
        self.assertIn("msl = E2Col(j)%moist_ref", src)

    def test_the_li7700_term_no_longer_reads_the_fallback_slot(self):
        """PrimaryWaterOutSlot names a trace gas when a project has no water."""
        src = code("src/src_rp/eddyflow-rp_main.f90")
        self.assertNotIn("Stats%chi(PrimaryWaterOutSlot())", src)


class TimeSeriesPairingStaysOnOneAnalyser(unittest.TestCase):
    """Cov(w, water) at this gas's lag, and point-by-point dilution, both
    need the water to have come down the same tube."""

    SITES = ("src/src_rp/timelag_handle.f90",
             "src/src_common/point_by_point_to_mixing_ratio.f90")

    def test_both_require_the_same_analyser(self):
        for rel in self.SITES:
            self.assertIn(
                "%instr%model /= E2Col(msl)%instr%model", code(rel),
                "%s pairs a gas with a hygrometer on another analyser; the "
                "row lag does not transfer between tubes" % rel)


class TheMeanTermsAcceptTheFallback(unittest.TestCase):
    """MoistTerms is the other half of the rule and must not gain the guard."""

    def test_moist_terms_starts_from_the_site_value(self):
        src = code("src/src_rp/fluxes23_rp.f90")
        body = src[src.index("subroutine MoistTerms"):]
        body = body[:body.index("end subroutine MoistTerms")]
        self.assertIn("sigma_out = Ambient%sigma", body,
                      "the per-gas terms must fall back to the site value, "
                      "not to `error`, when a gas names no usable hygrometer")
        self.assertIn("msl = E2Col(gas)%moist_ref", body)


if __name__ == "__main__":
    unittest.main()
