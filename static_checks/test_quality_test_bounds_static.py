"""KID, ZCD and the correlation difference are produced for every gas.

`write_out_fluxnet.f90` emits one LGD, one KID and one ZCD per *configured*
gas - `do var = u, ts + nFluxnetLayoutSlots`. LGD's producer was widened with
the rest of the flux chain; KID's and ZCD's were not, and `Essentials` is not
cleared between periods, so on the 8-gas fixture `N2O_KID`, `CO2_2_KID`,
`H2O_2_KID` and `N2O_2_KID` all read exactly `0.00000` and their ZCDs `0`.
A kurtosis index of zero is a claim about the data, not a missing value.

`Fisher` had the same bound and one more defect behind it: it was handed
`E2Primes(:, 1:GHGNumVar)` - 68 columns - together with `ncol =
size(E2Primes, 2)`, which is E2NumVar. The explicit-shape dummy therefore
described half again as many columns as were passed, and
`CorrelationMatrixNoError` read past them. Inert only because the loops
stopped at the fourth gas and never reached the overrun; widening them
without fixing the count would have walked straight into it.
"""

from pathlib import Path
import re
import unittest


ROOT = Path(__file__).resolve().parents[1]

KID = "src/src_rp/kid.f90"
FISHER = "src/src_rp/fisher.f90"
MAIN = "src/src_rp/eddyflow-rp_main.f90"


def code(path):
    return "\n".join(ln for ln in (ROOT / path).read_text(encoding="utf-8").splitlines()
                     if not ln.lstrip().startswith("!"))


class TheProducersCoverEveryConfiguredGas(unittest.TestCase):
    def test_kid_is_not_bounded_at_the_fourth_gas(self):
        src = code(KID)
        self.assertNotRegex(src, r"=\s*co2\s*,\s*gas4",
                            "KID must run firstGas..lastGas; the FLUXNET "
                            "writer emits a column per configured gas")
        self.assertRegex(src, r"=\s*firstGas\s*,\s*lastGas")

    def test_kid_still_has_its_sentinel_arm(self):
        """The widening is only safe because an absent slot is written as
        error/ierror rather than left at whatever Essentials carried over."""
        src = code(KID)
        self.assertRegex(src, r"Essentials%KID\(icol\)\s*=\s*error")
        self.assertRegex(src, r"Essentials%ZCD\(icol\)\s*=\s*ierror")

    def test_fisher_is_not_bounded_at_the_fourth_gas(self):
        src = code(FISHER)
        self.assertNotRegex(src, r"=\s*u\s*,\s*gas4",
                            "Fisher must run u..lastGas; CorrDiff is "
                            "dimensioned (GHGNumVar, GHGNumVar)")
        self.assertGreaterEqual(
            len(re.findall(r"=\s*u\s*,\s*lastGas", code(FISHER))), 3,
            "all three of Fisher's loops must span the full range")


class FisherIsToldHowManyColumnsItWasGiven(unittest.TestCase):
    def test_the_call_passes_the_section_width_not_the_array_width(self):
        src = code(MAIN)
        call = re.search(r"call Fisher\(([^)]*\)[^)]*)\)", src)
        self.assertIsNotNone(call, "the Fisher call site should still exist")
        text = call.group(1)
        self.assertIn("1:GHGNumVar", text)
        self.assertNotIn("size(E2Primes, 2)", text,
                         "ncol must describe the section passed (GHGNumVar), "
                         "not E2Primes's full second dimension (E2NumVar)")


if __name__ == "__main__":
    unittest.main()
