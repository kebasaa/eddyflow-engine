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

The distinction these checks used to encode has since been withdrawn, and the
reason is worth recording where it was asserted.

Time-series operations - the covariance and the point-by-point dilution - used
the gas's reference AND required it to be on the same analyser, on the ground
that a hygrometer down another tube has a different lag. Mean WPL terms - sigma
and rho_w - used the reference and accepted ResolveGasRef's fallback to "the
first hygrometer anywhere".

Splitting them that way reproduced the very fault described above. A gas whose
own analyser carries no hygrometer took sigma and rho_w from the borrowed water
and then got no flux term and no dilution: the correction ran half-built, from a
water the other two halves held it did not share. Declining is not the neutral
choice it looks like.

All three honour the reference now. The lag objection was real and is answered
by taking the covariance at the *hygrometer's* lag rather than the gas's, and
the compromise is announced through ExceptionHandler(106) instead of being
silently absorbed. `test_moisture_reference_static.py` owns that rule.

What must not come back is any of them reading the site's water directly.
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


class TimeSeriesPairingHonoursTheReference(unittest.TestCase):
    """Cov(w, water) and point-by-point dilution use the water the gas names,
    across analysers as well as on one.

    They refused it, and the refusal is what left a gas corrected by half a
    term. The lag objection behind the refusal is answered by choosing the lag;
    `test_moisture_reference_static.py` pins that, and this only guards against
    the flat refusal returning.
    """

    SITES = ("src/src_rp/timelag_handle.f90",
             "src/src_common/point_by_point_to_mixing_ratio.f90")

    def test_neither_refuses_the_pairing(self):
        for rel in self.SITES:
            self.assertNotIn(
                "%instr%model /= E2Col(msl)%instr%model) cycle", code(rel),
                "%s is declining the water its own record names, which leaves "
                "MoistTerms correcting the gas with a water this site says it "
                "does not share" % rel)


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
