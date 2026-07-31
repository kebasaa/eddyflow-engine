"""The LI-7550 analog-filter correction is dormant, and must not wake up
slot-bound.

`LI7550_AnalogSignalsTransferFunctions` has an arm `case(co2, h2o)` that gives
the block-averaging term to slots 5 and 6 only. On a two-analyser site those
slots need not be an LI-7500/7200's channels at all, and a second LI-7200's
CO2 sits well past them - so read as slot numbers the arm is wrong twice over.

It is currently unreachable. `write_processing_project_variables` hard-sets
both `hf_correct_ghg_ba` and `hf_correct_ghg_zoh` to `.false.`, with the
assignments that would read them from the project file commented out, and
`BPCF_LI7550AnalogFilters` - the only caller that ever passes a gas - is
guarded on those two flags. The other three call sites pass u, w and ts.

So it was left alone: widening it would be a change no fixture can observe,
and this repository has already paid for one correction that looked right and
was shipped without a perturbation that reached it.

What this check defends is the *next* person. If the flags become live again,
the arm has to be revisited at the same time - and the fix is to select by
instrument model in the caller, since the LI-7550 is the interface box of the
LI-7500/7200 and the term belongs to the gases on such an analyser rather
than to a slot range.
"""

from pathlib import Path
import re
import unittest


ROOT = Path(__file__).resolve().parents[1]

PROJVARS = "src/src_common/write_processing_project_variables.f90"
TF = "src/src_common/bpcf_analytic_transfer_functions.f90"


def code(path):
    return "\n".join(ln for ln in (ROOT / path).read_text(encoding="utf-8").splitlines()
                     if not ln.lstrip().startswith("!"))


class TheCorrectionIsStillDormant(unittest.TestCase):
    """If either assertion fails the feature has been re-enabled, and the
    slot-bound arm in bpcf_analytic_transfer_functions.f90 must be converted
    to an instrument-model test in the same change."""

    def test_both_flags_are_forced_false(self):
        src = code(PROJVARS)
        for flag in ("hf_correct_ghg_ba", "hf_correct_ghg_zoh"):
            self.assertRegex(
                src, r"EddyFlowProj%%%s\s*=\s*\.false\." % flag,
                "%s is no longer forced false - the LI-7550 arm is reachable "
                "again and is still bound to slots 5 and 6" % flag)

    def test_neither_flag_is_read_from_the_project_file(self):
        src = code(PROJVARS)
        for tag in ("EPPrjCTags(46)", "EPPrjCTags(47)"):
            self.assertNotIn(
                tag, src,
                "%s is being read again, so the correction is live; convert "
                "the case(co2, h2o) arm to select by instrument model" % tag)


class TheArmSaysWhyItIsLeftAlone(unittest.TestCase):
    def test_the_dormancy_is_documented_where_the_arm_is(self):
        """A reader who finds case(co2, h2o) must not have to discover the
        dormancy for themselves before deciding whether it is a defect."""
        text = (ROOT / TF).read_text(encoding="utf-8")
        arm = text.split("case(co2, h2o)")[0][-2000:]
        self.assertIn("DORMANT", arm)
        self.assertIn("hf_correct_ghg_ba", arm)


if __name__ == "__main__":
    unittest.main()
