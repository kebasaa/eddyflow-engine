"""The LI-7550 analog-filter correction is dormant, and correct if it wakes.

The block-averaging term describes what the LI-7550 does to the analog signals
it digitises. It is the interface box of the LI-7500 and LI-7200, so the term
belongs to the gases measured on such an analyser and to no others.

It used to be given to slots 5 and 6, on the reasoning that those analysers
measure CO2 and H2O. But a slot is not an analyser: on a site whose first
records are a QCL's, the term went to the QCL's channels and the LI-7200's
were left uncorrected, and a second LI-7200's CO2 sits well past slot 6
either way. The selection is made in BPCF_LI7550AnalogFilters now, from each
gas's own instrument record, and the arm in the transfer function takes what
it is given.

It is still unreachable. write_processing_project_variables hard-sets both
hf_correct_ghg_ba and hf_correct_ghg_zoh to .false., with the assignments that
would read them from the project file commented out, and
BPCF_LI7550AnalogFilters - the only caller that passes a gas - is guarded on
those two flags. The other three call sites pass u, w and ts.

So no fixture can observe any of this, which is exactly why it is pinned here:
the selection was fixed while it was cheap to reason about rather than left
for whoever re-enables the flags to discover. This check fails if the arm
returns to naming slots, and it fails if the flags come back without the
instrument test still being in place.
"""

from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]

PROJVARS = "src/src_common/write_processing_project_variables.f90"
TF = "src/src_common/bpcf_analytic_transfer_functions.f90"
CALLER = "src/src_common/bpcf_li7550_analog_filters.f90"


def code(path):
    return "\n".join(ln for ln in (ROOT / path).read_text(encoding="utf-8").splitlines()
                     if not ln.lstrip().startswith("!"))


class TheCorrectionIsStillDormant(unittest.TestCase):
    def test_both_flags_are_forced_false(self):
        src = code(PROJVARS)
        for flag in ("hf_correct_ghg_ba", "hf_correct_ghg_zoh"):
            self.assertRegex(
                src, r"EddyFlowProj%%%s\s*=\s*\.false\." % flag,
                "%s is no longer forced false. That is allowed, but the "
                "correction becomes observable, so the instrument selection "
                "below now needs a fixture behind it" % flag)

    def test_neither_flag_is_read_from_the_project_file(self):
        src = code(PROJVARS)
        for tag in ("EPPrjCTags(46)", "EPPrjCTags(47)"):
            self.assertNotIn(
                tag, src,
                "%s is being read again, so the correction is live and wants "
                "a fixture that reaches it" % tag)


class TheSelectionIsByInstrumentNotBySlot(unittest.TestCase):
    def test_the_caller_selects_on_the_analyser_model(self):
        """Which gases sit behind an LI-7550 is a question about instruments."""
        src = code(CALLER)
        self.assertIn("do gas = firstGas, lastGas", src)
        for model in ("li7500", "li7200"):
            self.assertIn("index(E2Col(gas)%%instr%%model, '%s')" % model, src,
                          "the term must follow the analyser, not the slot")
        for slot in ("histGas1", "histGas2", "histGas3", "histGas4"):
            self.assertNotIn(
                "size(nf), %s," % slot, src,
                "%s is back as the thing passed to the transfer function; "
                "a slot is not an analyser" % slot)

    def test_the_arm_takes_whatever_gas_it_is_given(self):
        src = code(TF)
        self.assertIn("case(firstGas:lastGas)", src)
        self.assertNotIn(
            "case(histGas1, histGas2)", src,
            "the arm names slots again; the caller decides which gases are "
            "on an LI-7550")

    def test_the_dormancy_is_documented_where_the_arm_is(self):
        """A reader who finds the arm must not have to discover for
        themselves that nothing reaches it."""
        text = (ROOT / TF).read_text(encoding="utf-8")
        arm = text.split("case(firstGas:lastGas)")[0][-2000:]
        self.assertIn("dormant", arm.lower())
        self.assertIn("hf_correct_ghg_ba", arm)


if __name__ == "__main__":
    unittest.main()
