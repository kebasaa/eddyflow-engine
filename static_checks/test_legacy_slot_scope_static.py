"""The legacy gas-slot names stay, but only where they earn it.

`co2 = 5, h2o = 6, ch4 = 7, n2o = 8, gas4 = 8` are not deleted yet. A handful
of uses are correct by design:

  on-disk aliases       LegacySpectralVarTag, HistoricGasSlot,
                        TimelagOptGasLabel and GasSlotFromDynMDTag accept the
                        four historical spellings so files written before the
                        records keep being readable. These only ever *widen*
                        what a reader accepts.

  external formats      The FLUXNET FP-In `GA_*` block has named CO2, H2O, CH4
                        and GS4 columns; the .metadata file defines its own
                        variable names. Neither is ours to renumber.

  compatibility modes   fix_out_format promises the fixed EddyPro 7.x column
                        set, which is four gas blocks whatever the project
                        holds. Widening it would break the compatibility the
                        flag exists to provide.

  the express guess    DefaultVarsSelection picks columns automatically in
                        embedded express mode, and is specified in terms of
                        the LI-COR trio: CO2 and H2O from a 7500 or 7200, CH4
                        from a 7700. "Which column did this site mean" has no
                        answer for an arbitrary species, so the guess names
                        those three and emits records for them; a site
                        measuring anything else declares its gases instead.

  fallbacks             PrimaryWaterOutSlot and PrimaryCarbonOutSlot fall back
                        to the historical slot when a project describes no
                        water or no CO2 - which is what the gates they replace
                        evaluated in that case.

The flat ini layer used to be on this list, and is not any more: read_ini_rp
and read_ini_fcc read every per-gas setting from its record, and a project
without records is refused rather than half-processed. That is what took
read_ini_rp from 76 occurrences to none.

What must not happen is a *new* `Set(:, co2)` or `E2Col(h2o)` appearing in a
processing path, which is how every defect this effort fixed began. So this
pins the count per file. A number going up fails; a number going down is the
work continuing and the expectation should be lowered to match.

Deliberately a count and not a pattern: the legitimate uses are too varied to
describe by regex, and a whitelist of exact strings would have to be rewritten
every time a comment moved.
"""

from pathlib import Path
import re
import subprocess
import unittest


ROOT = Path(__file__).resolve().parents[1]

SLOT_NAMES = re.compile(r"\b(co2|h2o|ch4|gas4|n2o)\b")

#: Occurrences outside comments, per file. Generated from the tree; see the
#: module docstring for why each file still has any.
ALLOWED = {
    "src/src_common/bpcf_analytic_transfer_functions.f90": 2,
    "src/src_common/bpcf_bandpass_spectral_corrections.f90": 1,
    "src/src_common/bpcf_li7550_analog_filters.f90": 4,
    "src/src_common/define_all_var_set.f90": 19,
    "src/src_common/define_e2_set.f90": 1,
    "src/src_common/define_used_variables.f90": 6,
    "src/src_common/gas4_output_units.f90": 18,
    "src/src_common/m_common_global_var.f90": 4,
    "src/src_common/m_typedef.f90": 12,
    "src/src_common/metadata_file_validation.f90": 6,
    "src/src_common/read_ex_record.f90": 3,
    "src/src_common/write_processing_project_variables.f90": 7,
    "src/src_fcc/cospectra_sorting_and_averaging.f90": 1,
    "src/src_fcc/init_out_files.f90": 30,
    "src/src_fcc/output_spectral_assessment_results.f90": 1,
    "src/src_fcc/read_ini_fcc.f90": 8,
    "src/src_fcc/read_spectral_assessment_file.f90": 2,
    "src/src_fcc/spectral_assessment_diagnostics.f90": 6,
    "src/src_fcc/write_out_fluxnet_fcc.f90": 20,
    "src/src_rp/add_to_timelag_opt_dataset.f90": 1,
    "src/src_rp/configure_for_express.f90": 5,
    "src/src_rp/default_vars_selection.f90": 21,
    "src/src_rp/drift_correction.f90": 1,
    "src/src_rp/init_fluxnet_file_rp.f90": 8,
    "src/src_rp/init_outfiles_rp.f90": 42,
    "src/src_rp/out_raw_data.f90": 3,
    "src/src_rp/pwb_timelag_handle.f90": 28,
    "src/src_rp/write_out_fluxnet.f90": 1,
}


def tracked_sources():
    out = subprocess.run(["git", "ls-files", "src"], cwd=ROOT,
                         capture_output=True, text=True, check=True).stdout
    return [p for p in out.split()
            if p.endswith(".f90") or p.endswith(".inc")]


def slot_uses(rel):
    """Occurrences outside full-line comments."""
    text = (ROOT / rel).read_text(encoding="utf-8", errors="replace")
    return sum(len(SLOT_NAMES.findall(ln)) for ln in text.splitlines()
               if not ln.lstrip().startswith("!"))


class LegacySlotUseIsPinned(unittest.TestCase):
    def test_no_file_gains_a_legacy_slot_reference(self):
        grew = []
        for rel in tracked_sources():
            n = slot_uses(rel)
            allowed = ALLOWED.get(rel, 0)
            if n > allowed:
                grew.append("%s: %d, expected at most %d" % (rel, n, allowed))
        self.assertFalse(
            grew,
            "these files gained a fixed gas-slot reference:\n  "
            + "\n  ".join(grew)
            + "\n\nA gas is addressed by its record. If the new use is one of "
              "the legitimate kinds - a legacy ini fallback, an on-disk alias, "
              "an external format, a compatibility mode or a no-water "
              "fallback - say which in a comment and raise the count here.")

    def test_the_expectations_are_not_stale(self):
        """A count that has dropped is the work continuing; lower it, so the
        file cannot silently regain what it gave up."""
        stale = []
        for rel, allowed in sorted(ALLOWED.items()):
            if not (ROOT / rel).is_file():
                stale.append("%s: listed but does not exist" % rel)
                continue
            n = slot_uses(rel)
            if n < allowed:
                stale.append("%s: %d now, expectation still %d" % (rel, n, allowed))
        self.assertFalse(
            stale,
            "lower these expectations to match the tree:\n  " + "\n  ".join(stale))


class TheSlotNamesThemselvesSurvive(unittest.TestCase):
    """Deleting them would break the layers that are supposed to use them."""

    def test_the_parameters_are_still_declared(self):
        src = (ROOT / "src/src_common/m_typedef.f90").read_text(
            encoding="utf-8", errors="replace")
        for name in ("co2", "h2o", "ch4", "n2o", "gas4"):
            self.assertRegex(
                src, r"integer, parameter :: %s\s*=" % name,
                "%s is still needed by the on-disk aliases and the fixed "
                "output format" % name)

    def test_the_dynamic_bounds_are_declared_from_the_capacity(self):
        """firstGas/lastGas are what everything else should iterate."""
        src = (ROOT / "src/src_common/m_typedef.f90").read_text(
            encoding="utf-8", errors="replace")
        self.assertIn("firstGas = NumAnemVar + 1", src)
        self.assertIn("lastGas  = GHGNumVar", src)


if __name__ == "__main__":
    unittest.main()
