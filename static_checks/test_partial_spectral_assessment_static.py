"""Guards the partial spectral assessment: some gases in situ, the rest not.

The assessment used to be all or nothing in three places at once, and fixing
any one of them alone left the run globally demoted:

  1. the writer decided whether the file existed by asking whether the primary
     hygrometer had a usable RH class, so a site whose water could not be
     fitted lost the CO2, N2O and COS blocks that had been fitted perfectly
     well;
  2. at correction time one unfitted gas - or an unfitted RH regression -
     demoted every gas to Moncrieff, which the code's own comment admitted;
  3. that per-period demotion was written back into EddyFlowProj%hf_meth, so
     the first bad period demoted every period after it for the rest of the
     run, whatever their own data looked like.

None of this is visible in a single-gas regression run, and the third is
invisible in any run whose assessment never succeeds - which is why it survived
so long. The checks below are on the mechanisms rather than on results.
"""

import re
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

WRITER = "src/src_fcc/output_spectral_assessment_results.f90"
READER = "src/src_fcc/read_spectral_assessment_file.f90"
CORRECT = "src/src_common/bpcf_bandpass_spectral_corrections.f90"
HELPER = "src/src_common/gas_slot_resolution.f90"


def read(rel):
    return (ROOT / rel).read_text(encoding="utf-8", errors="replace")


def code(rel):
    """Source with comment-only lines dropped, so a phrase in a comment cannot
    satisfy a check about what the code does."""
    out = []
    for line in read(rel).splitlines():
        stripped = line.strip()
        if stripped.startswith("!"):
            continue
        out.append(line)
    return "\n".join(out)


class OneQuestionAskedOnceTests(unittest.TestCase):
    def test_the_fit_test_is_shared_not_reimplemented(self):
        """The writer tested whether water's ensemble spectra existed while the
        report tested whether each gas had been fitted, so the file could be
        refused while the report called four gases PASS."""
        helper = code(HELPER)
        self.assertIn("logical function GasHasSpectralFit", helper)
        for user in (WRITER, READER, "src/src_fcc/spectral_assessment_diagnostics.f90"):
            self.assertIn("GasHasSpectralFit", code(user),
                          f"{user} answers 'is this gas fitted' its own way again")


class WriterTests(unittest.TestCase):
    def test_the_file_is_written_for_any_fitted_gas(self):
        body = code(WRITER)
        self.assertIn("if (n_fitted == 0) then", body)
        gate = body[body.index("if (FCCsetup%do_spectral_assessment) then"):]
        gate = gate[:gate.index("call LogSay(' Writing spectral assessment")]
        self.assertNotIn("goodj == ierror", gate,
                         "goodj is water's first usable RH class and belongs "
                         "to the H2O spectra file, not to this gate")

    def test_a_partial_file_says_so(self):
        self.assertIn("call ExceptionHandler(110)", code(WRITER))
        self.assertIn("call ExceptionHandler(110)", code(READER),
                      "a partial file read from disk is just as partial as one "
                      "written on the fly")

    def test_the_reader_does_not_key_the_warning_off_short_file(self):
        """short_file deliberately skips hygrometers, and water unfitted while
        the gases are fitted is the common case - keying off it would stay
        silent in exactly the situation the warning exists for."""
        body = code(READER)
        m = re.search(r"if \(n_fitted > 0[^\n]*\)\s*call ExceptionHandler\(110\)", body)
        self.assertIsNotNone(m, "the partial-file warning is gone")
        self.assertNotIn("short_file", m.group(0))


class CorrectionTests(unittest.TestCase):
    def test_the_method_is_chosen_per_gas(self):
        body = code(CORRECT)
        self.assertIn("insitu_ok", body)
        self.assertIn("analytic_only", body,
                      "the gases the in-situ method could not take must still "
                      "be corrected, by the analytic method")
        #> Order matters: the in-situ routines write across the whole gas range
        #> on their way past, so the analytic pass has to come second to leave
        #> the fitted gases standing and replace only the rest.
        self.assertLess(body.index("select case(trim(adjustl(actual_hf_method)))"),
                        body.index("analytic_only = loc_var_present"),
                        "the analytic pass must follow the in-situ dispatch")

    def test_the_fallback_does_not_outlive_its_period(self):
        self.assertNotIn(
            "EddyFlowProj%hf_meth = actual_hf_method", code(CORRECT),
            "writing the per-period fallback back into project state makes "
            "one bad period demote every later one, because the select case "
            "that re-evaluates it no longer matches an in-situ method")

    def test_switching_wholesale_still_reports_itself(self):
        """Error(69) says the method was switched. It should still fire when
        that is true of every gas, and not be lost in the per-gas rework."""
        body = code(CORRECT)
        self.assertIn("call ExceptionHandler(69)", body)
        self.assertIn("actual_hf_method = 'moncrieff_97'", body)


class FixtureTests(unittest.TestCase):
    def test_the_partial_fixture_exists_and_is_partial(self):
        sa = ROOT / "tests/regression/sa_n_gas_water_unfitted.txt"
        self.assertTrue(sa.exists(), "the partial assessment fixture is gone")
        text = sa.read_text(encoding="utf-8", errors="replace")
        water = [ln for ln in text.splitlines() if ln.startswith("RH class")]
        self.assertTrue(water, "the water block vanished from the fixture")
        self.assertTrue(all("-9999" in ln for ln in water),
                        "every RH class must be unfitted for this to be the "
                        "water-failed case")
        self.assertIn("January            =     0.20000     1.00000", text,
                      "the gas block must stay fitted, or the file is not "
                      "partial - it is empty")


if __name__ == "__main__":
    unittest.main()
