"""A gas may name the biomet as its humidity, and eight places must agree.

`moist_ref` used to be a gas slot or nothing. It can now also be
`biometMoistRef`, and that value has to survive from the project file to every
correction that reads it. The failure mode if one site misses it is the worst
kind this codebase has: a gas corrected with a humidity nobody chose, in a run
that reports no error.

Two shapes of miss, and both have already happened once here:

  - **Coercion.** `DefineE2Set` defaults an unresolved reference to the primary
    hygrometer. The sentinel is outside `firstGas..lastGas` by construction, so
    left to the bounds test that default would re-point every gas the user sent
    to the biomet, and the interface's selection would do nothing at all.

  - **A bounds check that silently declines.** Several readers cycle on an
    out-of-range reference. For some that is exactly right - there is no
    high-frequency biomet series, so the water-flux covariance and the
    point-by-point dilution genuinely have nothing to do - and for others it
    would drop a correction the user asked for. The difference is not visible
    from the shape of the code, so each one is pinned with its reason.

The resolution order changed too. It was: explicit, the gas's own analyser,
then *the first H2O on any analyser*. That last arm handed a gas whose own
instrument carries no hygrometer another instrument's water, drawn through a
different cell at a different time lag - the compromise Warning(106) exists to
announce. It is now the biomet instead, which measures the air rather than the
inside of an unrelated analyser. Borrowing remains reachable by asking for it.
"""

from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]

TYPEDEF = "src/src_common/m_typedef.f90"
E2SET = "src/src_common/define_e2_set.f90"
GLOBALS = "src/src_common/m_common_global_var.f90"
READ_INI = "src/src_rp/read_ini_rp.f90"
FLUXES23 = "src/src_rp/fluxes23_rp.f90"
DILUTION = "src/src_common/point_by_point_to_mixing_ratio.f90"
MOLEFRAC = "src/src_rp/molefractions_and_mixingratios.f90"
DRIFT = "src/src_rp/drift_correction.f90"
RP_MAIN = "src/src_rp/eddyflow-rp_main.f90"
TIMELAG = "src/src_rp/timelag_handle.f90"
FLUXES0 = "src/src_rp/fluxes0_rp.f90"
EX_WRITER = "src/src_rp/write_out_fluxnet.f90"


def read(path):
    return (ROOT / path).read_text(encoding="utf-8")


def code(path):
    return "\n".join(ln for ln in read(path).splitlines()
                     if not ln.lstrip().startswith("!"))


class TheSentinelIsNamed(unittest.TestCase):
    def test_it_is_a_parameter(self):
        self.assertIn("integer, parameter :: biometMoistRef = -1", code(TYPEDEF))

    def test_it_cannot_collide_with_a_slot(self):
        """Negative, so no arithmetic on firstGas ever produces it."""
        body = code(TYPEDEF)
        at = body.index("biometMoistRef = ")
        self.assertIn("-1", body[at: at + 40])

    def test_the_interface_uses_the_same_number(self):
        """One value written down twice, joined only by the project file."""
        gui = ROOT.parent / "eddyflow-gui" / "src" / "measurement_record.h"
        if not gui.is_file():
            self.skipTest("interface tree not beside this one")
        self.assertIn("biometMoistureRef = -1",
                      gui.read_text(encoding="utf-8"))


class ResolutionPrefersTheGasOwnAnalyserThenTheBiomet(unittest.TestCase):
    def setUp(self):
        src = code(E2SET)
        start = src.index("integer function ResolveGasRef")
        self.block = src[start: src.index("end function ResolveGasRef")]

    def test_an_explicit_biomet_choice_is_honoured(self):
        self.assertIn("if (ref == biometMoistRef) then", self.block)
        self.assertIn("ResolveGasRef = biometMoistRef", self.block)

    def test_it_is_gated_on_the_project_having_one(self):
        """A reference to a sensor the project does not name is not a choice,
        it is a stale file."""
        self.assertIn("BiometRhConfigured", self.block)

    def test_the_borrow_arm_is_gone(self):
        """`first record of the wanted species, whichever instrument`."""
        self.assertNotIn("if (trim(adjustl(candidate)) == trim(wanted)) then\n"
                         "            ResolveGasRef = slot", self.block,
                         "borrowing another analyser's water must not be an "
                         "automatic outcome any more")

    def test_the_biomet_arm_is_water_only(self):
        """The function takes a species because it was written to be general.
        The biomet has a humidity and nothing else."""
        self.assertIn("trim(wanted) == 'H2O'", self.block)

    def test_the_flag_is_set_from_the_ini(self):
        """biomet%val is a per-period value and is not retrieved when
        DefineE2Set runs, so the question has to be asked of the project."""
        self.assertIn("logical :: BiometRhConfigured", code(GLOBALS))
        self.assertIn("BiometRhConfigured = bSetup%sel(bRH) > 0", code(READ_INI))


class TheDefaultDoesNotEatTheSentinel(unittest.TestCase):
    """DefineE2Set points an unresolved reference at the primary hygrometer.
    The sentinel is out of slot range, so without an explicit skip every gas
    the user sent to the biomet is quietly re-pointed."""

    def test_the_coercion_skips_it(self):
        body = code(E2SET)
        at = body.index("E2Col(j)%moist_ref = wsl")
        before = body[max(0, at - 400): at]
        self.assertIn("if (E2Col(j)%moist_ref == biometMoistRef) cycle", before)


class EveryConsumerDecidesDeliberately(unittest.TestCase):
    def test_the_wpl_terms_recognise_it(self):
        """MoistTerms already defaults to the site scalars, which are the
        biomet values - but "asked for the biomet" and "nothing resolved" are
        different statements that happen to agree, and the code should say
        which one it means."""
        self.assertIn("if (msl == biometMoistRef) return", code(FLUXES23))

    def test_the_mean_conversion_uses_the_biomet_concentrations(self):
        body = code(MOLEFRAC)
        self.assertIn("if (msl == biometMoistRef) then", body)
        self.assertIn("waterR   = Ambient%r_biomet", body)
        self.assertIn("waterChi = Ambient%chi_biomet", body)

    def test_the_drift_correction_uses_them(self):
        body = code(DRIFT)
        self.assertIn("if (msl == biometMoistRef) then", body)
        self.assertIn("chi_moist = Ambient%chi_biomet", body)
        self.assertNotIn("0.15d0 * Stats%chi(msl)", body,
                         "the broadening must read the resolved source, not "
                         "the referenced slot directly")

    def test_the_li7700_multipliers_use_them(self):
        body = code(RP_MAIN)
        self.assertIn("chi_moist = Ambient%chi_biomet", body)
        self.assertNotIn("Multipliers7700(Stats%Pr, Ambient%Ta, &\n"
                         "                    Stats%chi(msl), &", body)

    #> Prose assertions run on the comment text with its markers removed and
    #> whitespace collapsed. A sentence broken across two comment lines is
    #> still that sentence - and with `!>` left in it reads as "point by !>
    #> point", so an assertion on the words would fail on the wrap rather than
    #> on the meaning.
    @staticmethod
    def _prose(text):
        stripped = []
        for ln in text.splitlines():
            s = ln.strip()
            if s.startswith("!>"):
                s = s[2:]
            elif s.startswith("!"):
                s = s[1:]
            stripped.append(s)
        return " ".join(" ".join(stripped).split())

    def test_the_covariance_declines_and_says_why(self):
        """A half-hourly RH sensor has no high-frequency series. Declining is
        the honest answer, not a gap."""
        body = read(TIMELAG)
        at = body.index("Stats%h2ocov_tl = error")
        self.assertIn("half-hourly", self._prose(body[at: at + 1400]))

    def test_the_dilution_declines_and_says_why(self):
        """Same reason, plus FluxParams has not run at that point in the
        period, so the biomet mole fraction does not exist yet."""
        body = read(DILUTION)
        at = body.index("msl = E2Col(gas)%moist_ref")
        block = self._prose(body[max(0, at - 1600): at])
        self.assertIn("point by point", block)
        self.assertIn("ordinary WPL", block)

    def test_the_ex_record_carries_the_terms_for_it(self):
        """FCC would reach the same numbers by falling through a bounds check,
        but a reader of the file would see -9999 against a corrected gas."""
        self.assertIn("if (indx == biometMoistRef) then", code(EX_WRITER))


class TheWarningOnlyFollowsAnExplicitBorrow(unittest.TestCase):
    def test_a_biomet_reference_raises_no_cross_analyser_warning(self):
        """It is not another analyser's water."""
        body = code(E2SET)
        at = body.index("call ExceptionHandler(106)")
        block = body[max(0, at - 700): at]
        self.assertIn("if (msl < firstGas .or. msl > lastGas) cycle", block)


if __name__ == "__main__":
    unittest.main()
