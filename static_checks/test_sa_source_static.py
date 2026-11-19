"""Where a gas's spectral transfer function comes from.

The in-situ correction fits a first-order transfer function per gas, from the
gas's spectrum against the sonic temperature's. A gas too noisy to fit - COS
and N2O on CH-LAE's MIRO came out at fc 0.002 and 0.035 Hz with Fn 0.04-0.07,
against 0.54 Hz and Fn 0.46 for the MIRO's own CO2 - is corrected with
nonsense. What the fit measures is the tube, the cell and the sensor's
response, which belong to the analyser, so such a gas can take a tube-mate's
fit instead: named per gas (gas_<i>_sa_source), or chosen automatically when
its own looks implausible (sa_auto_source).

Rules this file keeps true, because a regression run with neither key set -
which is every existing fixture - cannot see any of them:

  * **Same analyser, no water.** Another analyser has another tube; water's
    cut-off depends on humidity through sorption. Both in the reader, which
    refuses such a pairing, and in the rule, which never proposes one.
  * **One answer for the gate and the correction.** The hub's in-situ gate
    and RetrieveLPTFpars both ask SpectralSourceSlot, and both read the
    source's own class for the month. A gate that tested the gas's own fit
    while the correction read the donor's would send a gas with no fit of
    its own to the analytic method although a donor was there.
  * **The assessment file keeps each gas's own fit.** The fitter, the writer
    and the reader never ask for the source, so a bad fit stays visible.
  * **The constants are the interface's.** The GUI previews the same rule on
    a selected file; tests/test_sa_source_contract_static.py there holds the
    two to the same values.
  * **Suggestions come last.** Every suggestion the run makes is printed
    together as its very last message before the banner.

Part of the EddyFlow engine's static checks.
"""

from pathlib import Path
import re
import unittest

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"


def raw(*parts):
    return SRC.joinpath(*parts).read_text(encoding="utf-8", errors="replace")


def code(*parts):
    """Source without comment lines."""
    return "\n".join(ln for ln in raw(*parts).splitlines()
                     if not ln.lstrip().startswith("!"))


def body(src, start, end):
    i = src.index(start)
    return src[i:src.index(end, i)]


GEN = (ROOT / "prj" / "gen_project_tags.py").read_text(encoding="utf-8")
TAGS = raw("src_fcc", "m_fx_global_var_mod.f90")
COMMON = raw("src_common", "m_common_global_var.f90")
TYPEDEF = code("src_common", "m_typedef.f90")
READER = code("src_fcc", "read_ini_fcc.f90")
SLOTS = code("src_common", "gas_slot_resolution.f90")
HUB = code("src_common", "bpcf_bandpass_spectral_corrections.f90")
AUX = code("src_common", "bpcf_aux_subs.f90")
DIAG = code("src_fcc", "spectral_assessment_diagnostics.f90")
MAIN = code("src_fcc", "eddyflow-fcc_main.f90")
INI = code("src_common", "edit_ini_file.f90")


class TheKeysExist(unittest.TestCase):

    def test_the_per_gas_source_is_a_record_field(self):
        fields = GEN[GEN.index("FCC_GAS_NUMERIC = ["):]
        fields = fields[:fields.index("]")]
        self.assertIn('"sa_source"', fields)
        self.assertIn("'gas_1_sa_source'", TAGS)
        self.assertIn("'gas_64_sa_source'", TAGS)

    def test_the_switch_sits_in_a_blank_slot_below_the_records(self):
        #> Below fccGasOriginN, so no record moves.
        self.assertRegex(GEN, r'10:\s*"sa_auto_source"')
        self.assertIn("SNTags(10)%Label  / 'sa_auto_source' /", TAGS)

    def test_the_reader_uses_the_seventh_field(self):
        self.assertIn("fccGasLeapN   = 7", COMMON)
        self.assertIn("fccGasOriginN + (gas - 1) * fccGasLeapN + 6", READER)
        self.assertIn("SNTagFound(10)", READER)

    def test_the_setup_carries_them_to_the_workers(self):
        sa = body(TYPEDEF, "type :: SASetupType", "end type SASetupType")
        self.assertIn("integer :: source(E2NumVar)", sa)
        self.assertIn("logical :: auto_source", sa)


class TheReaderRefusesWhatCannotBeShared(unittest.TestCase):

    def test_absent_is_automatic(self):
        self.assertIn("FCCsetup%SA%source = SASourceAuto", READER)
        self.assertIn("FCCsetup%SA%auto_source = .false.", READER)

    def test_water_never_takes_part(self):
        block = READER[READER.index("FCCsetup%SA%source = SASourceAuto"):]
        block = block[:block.index("end do")]
        self.assertIn("if (GasSlotIsWater(slot)) cycle", block)
        self.assertIn("SpectralTubeMates(slot, firstGas + donor - 1)", block)
        self.assertIn("donor == gas", block)

    def test_tube_mates_are_one_named_analyser_without_water(self):
        f = body(SLOTS, "logical function SpectralTubeMates",
                 "end function SpectralTubeMates")
        self.assertIn("GasSlotIsWater(slot_a)", f)
        self.assertIn("GasSlotIsWater(slot_b)", f)
        self.assertIn("trim(ia) == 'other'", f)
        self.assertIn("SpectralTubeMates = trim(ia) == trim(ib)", f)
        #> Never by model, which cannot tell two LI-7200s apart
        self.assertNotIn("%model", f)


class TheRule(unittest.TestCase):

    def test_the_constants(self):
        self.assertRegex(COMMON, r"SaFnMin\s*=\s*0\.2d0")
        self.assertRegex(COMMON, r"SaFcRatioMin\s*=\s*0\.25d0")

    def test_a_usable_fit(self):
        f = body(SLOTS, "logical function SpectralFitUsable",
                 "end function SpectralFitUsable")
        self.assertIn("RegPar(gas_slot, cls)%fc <= 0d0", f)
        self.assertIn("RegPar(gas_slot, cls)%Fn >= SaFnMin", f)

    def test_the_reference_is_the_best_resolved_tube_mate(self):
        f = body(SLOTS, "integer function SpectralTubeReference",
                 "end function SpectralTubeReference")
        self.assertIn("SpectralTubeMates(gas_slot, cand)", f)
        self.assertIn("cls = SA%class(cand, month)", f)
        self.assertIn("RegPar(cand, cls)%Fn > best_fn", f)

    def test_implausible_is_unusable_or_a_factor_off(self):
        f = body(SLOTS, "integer function SpectralSuggestedSource",
                 "end function SpectralSuggestedSource")
        self.assertIn("if (ref == 0 .or. ref == gas_slot) return", f)
        self.assertIn("ratio < SaFcRatioMin .or. ratio > 1d0 / SaFcRatioMin", f)

    def test_who_decides(self):
        f = body(SLOTS, "integer function SpectralSourceSlot",
                 "end function SpectralSourceSlot")
        #> Water first, then the project's own word, then the switch
        order = [f.index("GasSlotIsWater(gas_slot)"),
                 f.index("SA%source(gas_slot) > 0"),
                 f.index("SA%source(gas_slot) == SASourceOwn"),
                 f.index("if (.not. SA%auto_source) return"),
                 f.index("SpectralSuggestedSource(gas_slot, month, SA)")]
        self.assertEqual(order, sorted(order))


class TheGateAndTheCorrectionAgree(unittest.TestCase):

    def test_the_hub_gate_asks_the_source(self):
        self.assertIn("src = SpectralSourceSlot(gas, month, LocSetup%SA)", HUB)
        self.assertIn("RegPar(src, LocSetup%SA%class(src, month))%fc == error", HUB)
        self.assertNotIn("RegPar(gas, LocSetup%SA%class(gas, month))", HUB)

    def test_the_cut_off_is_read_from_the_source(self):
        f = body(AUX, "subroutine RetrieveLPTFpars", "end subroutine RetrieveLPTFpars")
        self.assertEqual(f.count("src = SpectralSourceSlot(gas, month, LocSetup%SA)"), 2)
        self.assertIn("f_c(gas) = RegPar(src, LocSetup%SA%class(src, month))%fc", f)
        self.assertIn("f_2(gas) = RegPar(src, LocSetup%SA%class(src, month))%f2", f)
        self.assertNotIn("RegPar(gas, LocSetup%SA%class(gas, month))", f)

    def test_the_readiness_outcome_follows_the_source(self):
        self.assertIn("corrected_ok = GasInsituWithSource(gas, donor)", DIAG)
        self.assertIn("trim(OutcomeLabel(corrected_ok))", DIAG)


class TheFileKeepsEachGasItsOwnFit(unittest.TestCase):

    def test_fitting_writing_and_reading_never_ask(self):
        for parts in (("src_fcc", "fit_tf_models.f90"),
                      ("src_fcc", "output_spectral_assessment_results.f90"),
                      ("src_fcc", "read_spectral_assessment_file.f90"),
                      ("src_fcc", "sa_rates.f90"),
                      ("src_fcc", "fit_rh_to_cutoff.f90")):
            self.assertNotIn("SpectralSourceSlot", code(*parts), parts[-1])


class SuggestionsComeLast(unittest.TestCase):

    def test_the_sources_are_reported_once_by_the_parent(self):
        self.assertIn("if (BatchIndex == 0) call ReportSpectralSources()", MAIN)

    def test_nothing_is_printed_between_the_suggestions_and_the_banner(self):
        start = MAIN.index("call ApplyAutomaticSpectralConfiguration(")
        tail = MAIN[start:]
        between = tail[tail.index("call ReportRunSuggestions()"):
                       tail.index("' Program EddyFlow executed gracefully.'")]
        calls = re.findall(r"call (\w+)\(", between)
        self.assertEqual(set(calls) - {"ReportRunSuggestions", "LogSay"}, set(),
                         "only the banner may follow the suggestions: %s" % calls)
        #> The blank line, the row of stars, and the banner line's own call,
        #> which opens before the text the slice ends at.
        self.assertLessEqual(calls.count("LogSay"), 3)

    def test_the_flux_limit_suggestions_are_gathered(self):
        f = body(DIAG, "subroutine ReportFluxLimitSuggestions",
                 "end subroutine ReportFluxLimitSuggestions")
        self.assertIn("call AddRunSuggestion('flux'", f)

    def test_a_source_suggestion_is_gathered(self):
        f = body(DIAG, "subroutine ReportSpectralSources",
                 "end subroutine ReportSpectralSources")
        self.assertIn("call AddRunSuggestion('sources'", f)

    def test_what_automatic_selection_applied_is_saved_with_the_results(self):
        #> Not advice for the next run but what corrected this one, so it
        #> goes into the output project whatever automatic_spectra_config says.
        tail = MAIN[MAIN.index("call ApplyAutomaticSpectralConfiguration("):]
        self.assertLess(tail.index("call SaveAppliedSpectralSources("),
                        tail.index("call ReportRunSuggestions()"))
        f = body(DIAG, "subroutine SaveAppliedSpectralSources",
                 "end subroutine SaveAppliedSpectralSources")
        self.assertNotIn("automatic_config", f)
        self.assertIn("call WriteAutomaticSpectralSource(", f)
        #> Only a gas automatic selection switched in every month: the key has
        #> no months, and a user's own choice is already in the project.
        r = body(DIAG, "subroutine ReportSpectralSources",
                 "end subroutine ReportSpectralSources")
        self.assertIn("FCCsetup%SA%auto_source &\n            .and. all_same", r)
        self.assertIn("FCCsetup%SA%source(gas) == SASourceAuto", r)

    def test_the_automatic_configuration_adds_the_key(self):
        #> An automatic gas has no gas_<i>_sa_source line, and EditIniFile
        #> rewrites only lines that are there.
        f = body(DIAG, "subroutine WriteAutomaticSpectralSource",
                 "end subroutine WriteAutomaticSpectralSource")
        self.assertIn("call SetIniKeyInSection(", f)
        self.assertIn("'FluxCorrection_SpectralAnalysis_General'", f)
        self.assertNotIn("EditIniFile", f)
        self.assertIn("subroutine SetIniKeyInSection(", INI)


if __name__ == "__main__":
    unittest.main()
