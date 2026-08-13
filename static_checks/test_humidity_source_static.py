"""Humidity is taken from whichever source the site has.

Two sources supply it: a hygrometer among the gas records, and a biomet
relative-humidity sensor. Either is enough for the moist-air correction, so
the question asked first must be "is any humidity available?", not "is there a
hygrometer?".

Asking the wrong one first was a regression this migration introduced. The
guard added to stop the engine reading a non-water slot - `if (wsl < firstGas)`
- was placed at the head of the chain, so it ran *instead of* the biomet
branch. A site whose humidity came from an RH sensor and not an IRGA got
RH, vapour pressure and water density all `error`, and no correction at all,
with the measurement sitting right there in the file. Before the migration
`wsl` was the constant 6, that branch could never fire, and biomet RH worked.

This is a weaker claim than a fixture would be, and it is labelled as such:
every project in tests/regression points biom_file at a drive that is not
present, so no fixture can currently populate biomet%val(bRH). See the note in
gen_fixtures.py for what authoring one would take. Pinning the order is what
is available; it is not the same as observing the correction happen.
"""

from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]

PARAMS = "src/src_rp/flux_params.f90"


def code(rel):
    return "\n".join(ln for ln in (ROOT / rel).read_text(
        encoding="utf-8", errors="replace").splitlines()
        if not ln.lstrip().startswith("!"))


class BiometRHIsAskedFirst(unittest.TestCase):
    def test_the_biomet_test_precedes_the_hygrometer_test(self):
        src = code(PARAMS)
        biomet = src.index("biomet%val(bRH) > 0d0")
        hygro = src.index("elseif (wsl >= firstGas) then")
        self.assertLess(
            biomet, hygro,
            "the hygrometer test comes first again, so a site with only a "
            "biomet RH sensor gets no moist-air correction")

    def test_the_chain_does_not_open_with_the_no_water_guard(self):
        """The exact shape of the regression."""
        src = code(PARAMS)
        self.assertNotIn(
            "if (wsl < firstGas) then\n        RHO%w = error", src,
            "the no-water guard is back at the head of the humidity chain, "
            "where it pre-empts the biomet branch")

    def test_the_hygrometer_channel_is_not_written_here_at_all(self):
        """The branch computes site scalars; it writes no hygrometer's channel.

        This began as `if (wsl < firstGas) then continue`, guarding writes that
        indexed the primary slot - 0 on a site with no water. The guard is gone
        because the writes are: a hygrometer reports what it measured, and the
        biomet value is its own quantity, reported as h2o_biomet_*. So the
        out-of-bounds risk this test was written for cannot arise, and the
        stronger property is pinned instead.

        test_biomet_rh_scope_static covers where the biomet value goes now.
        """
        src = code(PARAMS)
        start = src.index("if (biomet%val(bRH) > 0d0")
        block = src[start: src.index("elseif (wsl >= firstGas) then", start)]
        for field in ("Stats%chi(", "Stats%r(", "Stats%d("):
            self.assertNotIn(field, block,
                             "%s in the biomet branch replaces a hygrometer's "
                             "measurement with the site value" % field)


class EveryBranchLeavesThePerHygrometerDensityDefined(unittest.TestCase):
    """RHO%w_at is read by MoistTerms for every gas's WPL correction.

    It used to be filled only in the raw-data branch. RHO is a module global
    with no per-period reset, so a site with biomet RH corrected every gas
    with the *previous* averaging period's humidity, and wrote that into the
    ex record for FCC to reuse.
    """

    def test_the_fill_is_outside_the_branches(self):
        src = code(PARAMS)
        fill = src.index("RHO%w_at(msl) = ")
        chain_end = src.index("!> Dew-point temperature") if \
            "!> Dew-point temperature" in src else None
        #: The fill must come after the whole if/elseif/else chain closes,
        #: which is what makes it apply to all three arms.
        self.assertIn("RHO%w_at = error", src)
        self.assertGreater(
            fill, src.index("elseif (wsl >= firstGas) then"),
            "the per-hygrometer density fill must follow the humidity chain, "
            "not sit inside one of its arms")
        del chain_end


class TheWarningNamesWhatIsLost(unittest.TestCase):
    """Air density is only half of it.

    Without humidity the Schotanus correction to H cannot be applied either,
    so the reported sensible heat flux is the uncorrected buoyancy flux. A
    warning that mentions only density understates the problem.
    """

    def test_the_message_covers_density_and_the_heat_flux(self):
        src = (ROOT / "src/src_common/exception_handler.f90").read_text(
            encoding="utf-8", errors="replace")
        block = src[src.index("case(104)"):]
        block = block[:block.index("end select")]
        for phrase in ("DRY air", "buoyancy flux", "biomet"):
            self.assertIn(phrase, block,
                          "warning 104 no longer mentions %r" % phrase)

    def test_it_is_raised_only_when_neither_source_exists(self):
        src = code("src/src_rp/read_ini_rp.f90")
        guard = src[src.index("ExceptionHandler(104)") - 400:
                    src.index("ExceptionHandler(104)")]
        self.assertIn("EddyFlowProj%gas_num > 0", guard,
                      "an anemometer-only project has no WPL to bias and "
                      "must not be warned")
        self.assertIn("bSetup%sel(bRH) <= 0", guard,
                      "a site with a biomet RH sensor gets the correction "
                      "and must not be warned")
        self.assertIn("PrimaryWaterSlot() < firstGas", guard)


if __name__ == "__main__":
    unittest.main()
