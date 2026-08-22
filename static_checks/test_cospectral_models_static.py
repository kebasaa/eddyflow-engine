"""The selectable analytic cospectrum, and the three things that make it safe.

``CospectralModel`` picks the shape every low-pass spectral correction is
integrated against. Three properties carry the whole design:

1. The default path is unreachable from the selector. ``CospectraMoncrieff97``
   runs unconditionally and the ``select case`` has no ``moncrieff_97`` branch,
   so a project that does not state the key gets bit-for-bit what it got before
   this routine existed. Asserted structurally here rather than trusted.

2. The Reynolds stress is never touched. Four of the five alternatives are
   scalar cospectra with no momentum form, so ``of(w_u)`` stays Moncrieff's for
   every option - which the regression confirms by holding ``Tau_scf`` fixed
   across all six models.

3. The constants are EddyUH's. Compared literal by literal against
   ``EC_Software_Spectral_Analysis/modelcospectra.m`` when that tree is beside
   this one. A transposed digit would move a flux by a percent and nothing
   else in the build would notice.
"""

import re
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MODELS = ROOT / "src" / "src_common" / "bpcf_cospectral_models.f90"
DECODER = ROOT / "src" / "src_common" / "write_processing_project_variables.f90"
EDDYUH = (ROOT.parent / "EddyUH_testing" / "EddyUH" / "EddyUH_1.7b_COS"
          / "EC_Software_Spectral_Analysis" / "modelcospectra.m")

#: The six callers of the dispatcher. Every analytic correction path has to go
#: through it, or the option would apply to some methods and not others with
#: nothing saying which.
CALLERS = (
    "bpcf_anemometric_fluxes.f90",
    "bpcf_Horst_97.f90",
    "bpcf_Ibrom_07.f90",
    "bpcf_li7550_analog_filters.f90",
    "bpcf_moncrieff_97.f90",
    "bpcf_only_lowfrequency_correction.f90",
)

#: ini value -> the name the decoder maps it to.
NAMES = {
    1: "kaimal_72",
    2: "sakai_01",
    3: "su_03",
    4: "moraes_08",
    5: "kristensen_97",
}


def read(path):
    return path.read_text(encoding="utf-8", errors="replace")


SRC = read(MODELS)
DISPATCH = SRC[SRC.index("subroutine CospectralModel("):]


def branches():
    """Each alternative's source, keyed by the name the decoder produces."""
    out = {}
    for name in NAMES.values():
        block = DISPATCH[DISPATCH.index("case ('%s')" % name):]
        end = block.find("case ('", 5)
        out[name] = block[:end if end > 0 else block.index("end select")]
    return out


def literals(text):
    """Every numeric literal in a Fortran fragment, as floats."""
    return sorted(set(float(x[:-2]) for x in re.findall(r"\d+\.?\d*d0", text)))


class TheDefaultPathIsUnchanged(unittest.TestCase):

    def test_moncrieff_runs_before_the_selector(self):
        body = DISPATCH[:DISPATCH.index("select case")]
        self.assertIn("call CospectraMoncrieff97(nf, kf, Cospectrum, zL, N)",
                      body,
                      "the unconditional Moncrieff call is what makes an "
                      "unstated project bit-identical")

    def test_the_selector_has_no_moncrieff_branch(self):
        #> If one were added it would overwrite what CospectraMoncrieff97 just
        #> wrote, and the default would stop being provably unchanged.
        cases = re.findall(r"case \('([a-z0-9_]+)'\)", DISPATCH)
        self.assertNotIn("moncrieff_97", cases)
        self.assertEqual(sorted(cases), sorted(NAMES.values()))

    def test_the_decoder_defaults_to_moncrieff_before_reading_the_tag(self):
        src = read(DECODER)
        i = src.index("EddyFlowProj%cosp_model = 'moncrieff_97'")
        j = src.index("EPPrjNTagFound(7)")
        self.assertLess(i, j,
                        "the literal default has to precede the guarded "
                        "override, or an absent key leaves the field holding "
                        "whatever the previous parse left")

    def test_an_unknown_value_falls_back_rather_than_being_kept(self):
        block = read(DECODER)
        block = block[block.index("EPPrjNTagFound(7)"):]
        block = block[:block.index("end select")]
        self.assertIn("case default", block)
        self.assertIn("EddyFlowProj%cosp_model = 'moncrieff_97'", block)


class TheReynoldsStressIsNeverTouched(unittest.TestCase):

    def test_no_branch_writes_the_momentum_slot(self):
        after = DISPATCH[DISPATCH.index("select case"):]
        self.assertNotIn("w_u", after,
                         "a model branch writes of(w_u); the four scalar-only "
                         "models have no momentum form to write")

    def test_the_scalar_helper_writes_temperature_and_every_gas(self):
        #> Naming slots one by one is how a fifth gas was once left holding the
        #> error value, which SpectralCorrectionFactors reads as "no correction
        #> available". The loop is the fix and has to stay a loop.
        helper = DISPATCH[DISPATCH.index("subroutine PutScalar"):]
        self.assertIn("Cospectrum(k)%of(w_ts) = value", helper)
        self.assertIn("do gas = firstGas, lastGas", helper)


class EveryAnalyticPathGoesThroughIt(unittest.TestCase):

    def test_no_caller_still_calls_moncrieff_directly(self):
        for name in CALLERS:
            src = read(ROOT / "src" / "src_common" / name)
            self.assertNotIn("call CospectraMoncrieff97(", src,
                             "%s bypasses the selector" % name)
            self.assertIn("call CospectralModel(nf, kf, Cospectrum, zL, nfreq)",
                          src, "%s does not call the selector" % name)

    def test_the_call_sites_match_the_signature(self):
        #> External subroutines with no interface block: an argument-list
        #> mismatch compiles clean and corrupts memory at run time. This is the
        #> only thing that catches it.
        sig = re.search(r"subroutine CospectralModel\(([^)]*)\)",
                        DISPATCH).group(1)
        formal = [a.strip() for a in sig.split(",")]
        self.assertEqual(formal, ["nf", "kf", "Cospectrum", "zL", "N"])
        for name in CALLERS:
            src = read(ROOT / "src" / "src_common" / name)
            actual = re.search(r"call CospectralModel\(([^)]*)\)",
                               src).group(1)
            self.assertEqual(len(actual.split(",")), len(formal),
                             "%s passes a different number of arguments" % name)


class TheConstantsAreEddyUHs(unittest.TestCase):

    #: The four single-form models, each one MATLAB line. Kaimal is excluded:
    #: it spans an if/else in both sources and is checked separately below.
    MARKERS = {
        "sakai_01": "Csakai            =",
        "su_03": "Csun              =",
        "moraes_08": "Cmoraes           =",
        "kristensen_97": "CKristensen       =",
    }

    @unittest.skipUnless(EDDYUH.is_file(), "EddyUH tree not beside this one")
    def test_each_branch_holds_exactly_eddyuhs_numbers(self):
        #> Set equality, not containment. This is why Kristensen's exponents
        #> are written 2*0.23 and 7/(6*0.23) rather than folded to 0.46 and
        #> 5.0725: folded, a transposed digit would be invisible here.
        m = read(EDDYUH)
        blocks = branches()
        for name, marker in self.MARKERS.items():
            line = m[m.index(marker):]
            line = line[:line.index(chr(10))]
            #> Everything after the first % is the citation comment, whose
            #> year would otherwise read as a coefficient.
            line = line.split("%")[0]
            want = sorted(set(float(x) for x in
                              re.findall(r"\d+\.?\d*", line)))
            self.assertEqual(literals(blocks[name]), want,
                             "%s: the Fortran and EddyUH's line do not hold "
                             "the same constants" % name)

    @unittest.skipUnless(EDDYUH.is_file(), "EddyUH tree not beside this one")
    def test_kaimals_two_branches_hold_eddyuhs_numbers(self):
        m = read(EDDYUH)
        block = branches()["kaimal_72"]
        #> Stable: H = 1+6.4*stab, n0 = 0.23*H^(3/4), then 0.88 and 1.5*x^2.1.
        #> The Fortran folds 3/4 to 0.75 and holds 2.1 as a literal, so the
        #> comparison is on the distinctive coefficients rather than the set.
        for c in ("6.4d0", "0.23d0", "0.88d0", "1.5d0", "2.1d0", "0.75d0"):
            self.assertIn(c, block, "the stable branch has lost %s" % c)
        #> Unstable: 11n/(1+13.3n)^(7/4) below n=1, 4n/(1+3.8n)^(7/3) above.
        for c in ("11d0", "13.3d0", "1.75d0", "4d0", "3.8d0"):
            self.assertIn(c, block, "the unstable branch has lost %s" % c)
        for c in ("6.4*stab", "0.23*H", "0.88", "13.3", "3.8"):
            self.assertIn(c, m, "EddyUH no longer spells Kaimal that way")

    def test_kristensen_drops_the_prefactor_and_says_why(self):
        block = branches()["kristensen_97"]
        self.assertNotIn("/ z", block)
        self.assertIn("cancels", block.lower(),
                      "dropping a factor from a published formula needs the "
                      "reason written where the formula is")


class TheChoiceIsDocumentedWhereItMatters(unittest.TestCase):

    def test_the_routine_says_the_library_is_not_eddyuhs_method(self):
        #> EddyUH uses these for diagnostic plots and corrects against measured
        #> cospectra. Selecting one here is not "running EddyUH's correction",
        #> and the header has to say so or the next reader will assume it is.
        head = SRC[SRC.index("PROVENANCE"):
                   SRC.index("subroutine CospectralModel(")]
        self.assertIn("not for its own correction", head)

    def test_the_moore_moncrieff_identity_is_recorded(self):
        self.assertIn("CMoore", SRC,
                      "EddyUH's CMoore IS this program's Moncrieff curve; "
                      "without that written down the two look like different "
                      "models and someone will add the second")

    def test_the_cancelling_normalisation_is_explained(self):
        head = SRC[SRC.index("WHAT CANCELS"):
                   SRC.index("subroutine CospectralModel(")]
        self.assertIn("IntCO / IntTFCO", head)


if __name__ == "__main__":
    unittest.main()
