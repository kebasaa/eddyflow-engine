"""The (co)spectra are *computed* for every configured gas, not just four.

The binned file learned to declare a column per configured gas before the
maths learned to fill one. `DetectFeasibleSpectraAndCospectra` enables a slot
whenever its column is present, and the writers loop firstGas..lastGas - but
`AllCospectra`, `NormalizeCoSpectra`, `ExpAvrgCospectra` and the ogive binning
all stopped at the fourth gas. The binned arrays are intent(out), so a slot
past it was never assigned and the writer's only guard, `/= error`, happily
passed whatever the stack held.

Observed on the 8-gas fixture before the fix: `f_nat*spec(n2o)`,
`spec(co2_2)`, `spec(h2o_2)`, `spec(n2o_2)` and their four cospectra read
exactly `0.00000` in all 50 bins - a normalised cospectrum of zero, which is
a physically meaningful claim and a false one. The same slots in
`full_cospectra` were undefined for the same reason. Two runs of an unchanged
tree disagreed on those columns, which is how "undefined" rather than "wrong"
was established.

The cospectral index space is the same space as the gas slots - w_u = u = 1,
w_co2 = co2 = 5 - so the fix is a terminator, `w_lastGas`, and the four
`w_*` names stay for the slots that legitimately mean the anemometer.
"""

from pathlib import Path
import re
import unittest


ROOT = Path(__file__).resolve().parents[1]

ANALYSIS = "src/src_rp/spectral_analysis.f90"
TYPEDEF = "src/src_common/m_typedef.f90"


def read(path):
    return (ROOT / path).read_text(encoding="utf-8")


def code(path):
    """Source with comment-only lines dropped.

    The file documents the historical bounds in prose, so a naive search
    would match the explanation rather than the code.
    """
    return "\n".join(ln for ln in read(path).splitlines()
                     if not ln.lstrip().startswith("!"))


class TheCospectralRangeHasATerminator(unittest.TestCase):
    def test_w_lastgas_is_defined_and_tracks_the_gas_slots(self):
        src = code(TYPEDEF)
        self.assertRegex(
            src, r"w_firstGas\s*=\s*firstGas",
            "w_firstGas must be derived from firstGas, not written as a literal")
        self.assertRegex(
            src, r"w_lastGas\s*=\s*lastGas",
            "w_lastGas must be derived from lastGas, so the two index spaces "
            "cannot drift apart")

    def test_the_four_historical_w_names_are_kept(self):
        """They still name real things - w_u..w_ts are anemometer channels,
        and w_co2..w_gas4 are read by the legacy out_full_cosp_* tags."""
        src = code(TYPEDEF)
        for name in ("w_u", "w_v", "w_w", "w_ts", "w_gas4"):
            self.assertRegex(src, r"\b%s\s*=" % name,
                             "%s must keep its definition" % name)


class NothingInTheSpectralChainStopsAtTheFourthGas(unittest.TestCase):
    def test_no_four_gas_bound_survives(self):
        src = code(ANALYSIS)
        for pattern, label in (
            (r"w_u\s*,\s*w_gas4", "do j = w_u, w_gas4"),
            (r"w_u\s*:\s*w_gas4", "%of(w_u:w_gas4)"),
            (r"=\s*u\s*,\s*gas4", "do j = u, gas4"),
            (r"\bu\s*:\s*gas4", "%of(u:gas4)"),
        ):
            self.assertNotRegex(
                src, pattern,
                "%s bounds the compute chain at the fourth gas while the "
                "writer declares a column per configured gas" % label)

    def test_the_binning_covers_the_whole_range(self):
        """Both binned families must accumulate over the full slot range;
        it is the assignment, not the loop bound alone, that was missing."""
        src = code(ANALYSIS)
        for family in ("BinnedSpectrum", "BinnedCospectrum",
                       "BinnedOgive", "BinnedCoOgive"):
            self.assertIn(family, src, "%s should still exist" % family)
        self.assertGreaterEqual(
            len(re.findall(r"%of\(u:GHGNumVar\)", src)), 8,
            "the spectra and ogive binning must span u:GHGNumVar")
        self.assertGreaterEqual(
            len(re.findall(r"%of\(w_u:w_lastGas\)", src)), 8,
            "the cospectra and co-ogive binning must span w_u:w_lastGas")


class ASkippedSlotSaysNotPerformed(unittest.TestCase):
    """The rule this repository keeps paying for: widening a loop promotes an
    unconfigured slot from never-consulted to consulted at its default. Here
    the default was uninitialised storage, so the initialisation has to be
    explicit and it has to be `error`."""

    def test_allcospectra_initialises_both_outputs(self):
        src = code(ANALYSIS)
        body = src.split("subroutine AllCospectra")[1].split("end subroutine")[0]
        self.assertRegex(
            body, r"Spectrum\(1:N/2 \+ 1\)%of\(j\)\s*=\s*error",
            "AllCospectra must set every spectrum slot to error before "
            "filling the feasible ones")
        self.assertRegex(
            body, r"Cospectrum\(1:N/2 \+ 1\)%of\(j\)\s*=\s*error",
            "AllCospectra must set every cospectrum slot to error first")

    def test_allogives_initialises_both_outputs(self):
        src = code(ANALYSIS)
        body = src.split("subroutine AllOgives")[1].split("end subroutine")[0]
        self.assertRegex(body, r"Ogive\(1:N/2 \+ 1\)%of\(j\)\s*=\s*error",
                         "AllOgives must set every ogive slot to error first")
        self.assertRegex(body, r"CoOgive\(1:N/2 \+ 1\)%of\(j\)\s*=\s*error",
                         "AllOgives must set every co-ogive slot to error first")


if __name__ == "__main__":
    unittest.main()
