"""A half-hour's results do not depend on which half-hours were processed before it.

Splitting the main pass across processes needs each half-hour to come out the
same whether the run reached it from the start or a worker began just before
it. Several values were carried from one half-hour into later ones for no
physical reason, and each made results depend on where a run began:

* ET storage (Stor%ET) was set only when computed, never reset, so a period
  without its own kept the last consecutive pair's, however far back.
* OverrideSettings set bu_corr = 'none' for good the first time a period had
  no LI-7500: one file missing a record switched the Burba correction off for
  the rest of the run.
* A skipped period's FLUXNET row took the NIGHT flag of the last processed
  period, and with embedded biomet its biomet values too.
* The sonic output rate, when the project leaves it unset, was defaulted once
  from whichever period first needed it.
* The main pass began from whatever the pre-passes left: the dynamic metadata
  settings of the run's last record wherever a record left a field blank, and
  the PWB streaming classifier's last settled lag.
* A gas named to the biomet humidity was diluted by the previous period's:
  FluxParams computes it, after the mole fractions are formed from it.

Part of the EddyFlow engine's static checks.
"""

from pathlib import Path
import re
import unittest

ROOT = Path(__file__).resolve().parents[1]


def code(rel):
    text = (ROOT / rel).read_text(encoding="utf-8", errors="replace")
    return "\n".join(ln for ln in text.splitlines() if not ln.lstrip().startswith("!"))


MAIN = code("src/src_rp/eddyflow-rp_main.f90")
STORAGE = code("src/src_rp/storage.f90")
OVERRIDE = code("src/src_rp/override_settings.f90")
SKIP = code("src/src_rp/write_out_fluxnet_only_biomet.f90")
BPCF = code("src/src_common/bpcf_bandpass_spectral_corrections.f90")
FLUXPARAMS = code("src/src_rp/flux_params.f90")


#> The main pass's loop, not the pre-passes' (to_periods_loop, pf_periods_loop).
MAIN_LOOP = MAIN.index("\n    periods_loop: do")


def main_pass_start():
    return MAIN[MAIN.rindex("NumberOfOkPeriods = 0", 0, MAIN_LOOP):MAIN_LOOP]


class EtStorageIsNeverStale(unittest.TestCase):

    def test_storage_resets_et_before_any_early_return(self):
        self.assertLess(STORAGE.index("Stor%ET = error"),
                        STORAGE.index("if (tmp_date /= Stats%date .or. tmp_time /= Stats%time) then"))

    def test_the_first_period_resets_et(self):
        i = MAIN.index("if(InitializeStorage) then")
        self.assertIn("Stor%ET = error", MAIN[i:MAIN.index("InitializeStorage = .false.", i)])


class BurbaIsDecidedPerPeriod(unittest.TestCase):

    def test_each_call_starts_from_the_projects_choice(self):
        self.assertIn("RPsetup%bu_corr = project_bu_corr", OVERRIDE)
        self.assertLess(OVERRIDE.index("RPsetup%bu_corr = project_bu_corr"),
                        OVERRIDE.index("if (.not. has_li7500) RPsetup%bu_corr = 'none'"))
        self.assertIn("character(32), save :: project_bu_corr", OVERRIDE)


class SkippedRowsDescribeTheirOwnPeriod(unittest.TestCase):

    def test_the_skipped_row_assesses_its_own_daytime(self):
        body = SKIP[SKIP.index("subroutine WriteOutFluxnetOnlyBiomet("):
                    SKIP.index("end subroutine WriteOutFluxnetOnlyBiomet")]
        self.assertIn("call AssessDaytime(Stats%date, Stats%time)", body)
        self.assertLess(body.index("call AssessDaytime("), body.index("if (FluxnetFileOpen) then"))

    def test_embedded_biomet_is_cleared_every_period(self):
        loop_top = MAIN[MAIN_LOOP:MAIN_LOOP + 2500]
        self.assertIn("if (EddyFlowProj%biomet_data == 'embedded') then", loop_top)
        self.assertIn("if (allocated(bAggr)) bAggr = error", loop_top)


class TheSonicOutputRateIsRestoredPerPeriod(unittest.TestCase):

    def test_it_is_saved_once_and_restored(self):
        self.assertIn("ProjSonicOutputRate = EddyFlowProj%sonic_output_rate", BPCF)
        self.assertIn("EddyFlowProj%sonic_output_rate = ProjSonicOutputRate", BPCF)


class TheMainPassStartsFromStartUpNotFromThePrepasses(unittest.TestCase):

    def test_dynamic_metadata_settings_are_kept_before_the_prepasses(self):
        self.assertLess(MAIN.index("StartupMetadata = Metadata"),
                        MAIN.index("call PlanPrepassBatches(toEndTimestampIndx"))

    def test_and_restored_when_the_main_pass_starts(self):
        start = main_pass_start()
        for line in ("Metadata%lat = StartupMetadata%lat", "Metadata%canopy_height = StartupMetadata%canopy_height",
                     "RPsetup%wdf_num_secs = StartupWdfNumSecs"):
            self.assertIn(line, start)

    def test_the_pwb_streaming_state_is_reset_there(self):
        start = main_pass_start()
        for line in ("pwb_last_optimal_lag = error", "pwb_has_previous = .false."):
            self.assertIn(line, start)


class PrepassWorkersOpenNoOutputFiles(unittest.TestCase):

    def test_startup_output_opens_skip_prepass_workers(self):
        # A production worker ('pr', 'pd') writes into a folder of its own and
        # needs these open; a pre-pass worker ('to', 'pf') writes nothing.
        self.assertRegex(MAIN, r"nbVars > 0 &\s*\n\s*\.and\. \(BatchIndex == 0 \.or\. BatchKind == 'pr' \.or\. BatchKind == 'pd'\)\) &\s*\n\s*call InitBiometOut\(\)")
        self.assertRegex(MAIN, r"if \(NumUserVar > 0 \.and\. \(BatchIndex == 0 \.or\. BatchKind == 'pr' &\s*\n\s*\.or\. BatchKind == 'pd'\)\) call InitUserOutFiles\(\)")


class MoleFractionsUseThisPeriodsBiometHumidity(unittest.TestCase):

    def test_it_is_computed_before_each_mole_fraction_that_precedes_flux_params(self):
        calls = [m.start() for m in re.finditer(r"call MoleFractionsAndMixingRatios\(\)", MAIN)]
        self.assertEqual(len(calls), 2)
        for i in calls:
            before = MAIN[:i].rstrip().splitlines()[-1]
            self.assertIn("call BiometWaterVapour()", before)

    def test_by_the_same_arithmetic_as_flux_params(self):
        body = FLUXPARAMS[FLUXPARAMS.index("subroutine BiometWaterVapour()"):]
        for expr in ("(dexp(77.345d0 + 0.0057d0 * Stats%T", "Ambient%Va / MW_H2O * 1d3",
                     "/ (1.d0 - Ambient%chi_biomet * 1d-3)", "Ambient%chi_biomet / Ambient%Va"):
            self.assertIn(expr, body)
            self.assertIn(expr, FLUXPARAMS[:FLUXPARAMS.index("subroutine BiometWaterVapour()")])


if __name__ == "__main__":
    unittest.main()
