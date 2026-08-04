"""The drift subsystem addresses gases by record, not by slot number.

The calibration-events reader packed one array with the gases' offset columns
at 5..8 and their reference columns at 9..12:

    integer, parameter :: co2_ref = 9, h2o_ref = 10, ch4_ref = 11, gas4_ref = 12
    ...
    do j = co2_ref, gas4_ref
        if (mdcol(j) /= 0) read(text_vars(mdcol(j)), *) Calib(i)%ref(j - 4)
    end do

The `- 4` is the whole problem in one expression: it only works because the
two families are exactly four wide and exactly four apart. A project with
eight gases had no name it could give the fifth one's reference column.

They are now two arrays, offcol and refcol, with no arithmetic relating them,
and the header names resolve through GasSlotFromDynMDTag - the record's own
label first, the four legacy spellings after, since dynamic metadata files in
the wild carry those and renaming them would break every one.

The same helper serves the raw-data scan in ReferenceCounts. That matters:
the signal-strength method takes its reference counts from the raw data and
its offsets from the metadata file, so two separate spellings of "which gas
is n2o_ref" could silently pair one gas's counts with another's calibration.

Coefficients are per-gas in the record schema (drift_dir_0..6,
drift_inv_0..6). A gas that supplies none keeps DriftCorr's `error` and is
skipped: the polynomials are per-channel instrument calibrations and there is
no general form for an arbitrary species on an arbitrary analyser, so drift
is opt-in per gas rather than defaulted.
"""

from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]

EVENTS = "src/src_rp/drift_retrieve_calibration_events.f90"
CORRECT = "src/src_rp/drift_correction.f90"
HELPER = "src/src_common/gas4_output_units.f90"
READINI = "src/src_rp/read_ini_rp.f90"
GENERATOR = "prj/gen_project_tags.py"

FIXTURES = ("tests/regression/base_drift.eddyflow",
            "tests/regression/base_drift_n_gas.eddyflow")


def code(path):
    return "\n".join(ln for ln in (ROOT / path).read_text(encoding="utf-8",
                                                          errors="replace").splitlines()
                     if not ln.lstrip().startswith("!"))


class TheReferenceColumnsAreResolved(unittest.TestCase):
    def test_no_hand_packed_reference_slots_remain(self):
        src = code(EVENTS)
        for token in ("co2_ref", "h2o_ref", "ch4_ref", "gas4_ref"):
            self.assertNotIn("%s   =" % token, src)
            self.assertNotIn("%s  =" % token, src)
        self.assertNotIn("- 4)", src,
                         "the offset relating the two index spaces is what "
                         "capped them at four gases")
        self.assertNotIn("mdcol", src,
                         "offsets and references must be separate arrays")

    def test_the_two_families_are_separate_arrays(self):
        src = code(EVENTS)
        self.assertIn("offcol(GHGNumVar)", src)
        self.assertIn("refcol(GHGNumVar)", src)

    def test_the_resolver_exists_and_accepts_the_legacy_names(self):
        src = code(HELPER)
        self.assertIn("integer function GasSlotFromDynMDTag", src)
        body = src.split("function GasSlotFromDynMDTag")[1].split("end function")[0]
        self.assertIn("SpectralVarTags", body,
                      "a record's own label is what lets a project name "
                      "n2o_ref or co2_2_ref")
        self.assertIn("HistoricGasSlot", body,
                      "co2/h2o/ch4/gas4 must keep working: every dynamic "
                      "metadata file in the wild spells them")

    def test_both_readers_go_through_it(self):
        """The metadata file and the raw-column scan must agree about which
        gas a name refers to, or signal-strength drift pairs one gas's counts
        with another's calibration."""
        for path in (EVENTS, CORRECT):
            self.assertIn("GasSlotFromDynMDTag", code(path))


class TheCorrectionIsPerGas(unittest.TestCase):
    def test_no_two_gas_slices_remain(self):
        src = code(CORRECT)
        self.assertNotIn("lDrift(co2:h2o)", src)
        self.assertNotIn("refCounts(co2:h2o)", src)
        self.assertNotIn("Set(:, co2)", src)
        self.assertNotIn("Set(:, h2o)", src)

    def test_water_is_resolved_not_indexed(self):
        self.assertNotIn("Stats%chi(h2o)", code(CORRECT),
                         "the broadening term reads the site's water, and the "
                         "sixth slot is water only by convention")
        self.assertIn("PrimaryWaterOutSlot()", code(CORRECT))

    def test_the_band_broadening_is_applied_to_carbon_dioxide_only(self):
        """0.15 is the water-vapour broadening of the CO2 band (LI-7200 manual
        Rev 5), not a general property of a trace gas. The h2o channel never
        had it; nor should N2O."""
        src = code(CORRECT)
        self.assertIn("0.15d0", src)
        self.assertIn("GasOutputLabel", src,
                      "the term must be chosen by species, not by slot")

    def test_the_debug_unit_is_gone(self):
        self.assertNotIn("write(987", code(CORRECT))

    def test_the_reference_counts_span_the_whole_gas_block(self):
        """%ri and %rf are filled in the main loop, not here, and they were
        the last (co2:h2o) slices in the drift chain.

        Everything either side spans firstGas:lastGas -
        DriftRetrieveCalibrationEvents fills %offset and %ref over it, and the
        `where` mask in drift_correction reads all four arrays over it. Calib
        is initialised to error, so a third gas's %ri and %rf stayed error,
        the mask excluded it, and a gas given its own reference and offset
        columns was silently never corrected.

        Only the signal_strength method reads %ri and %rf; the linear method
        uses %offset alone. Every drift fixture is linear, so no fixture
        covers this and the check is the only guard.
        """
        src = code("src/src_rp/eddyflow-rp_main.f90")
        for field in ("ri", "rf"):
            self.assertNotIn(
                "%%%s(co2:h2o)" % field, src,
                "Calib%%%s is assigned over two slots; the correction reads "
                "it over firstGas:lastGas" % field)
        self.assertIn("Calib(0)%ri(firstGas:lastGas)", src)
        self.assertNotIn("refCounts(co2:h2o)", src)


class TheCoefficientsArePerRecord(unittest.TestCase):
    def test_the_generator_emits_them(self):
        src = (ROOT / GENERATOR).read_text(encoding="utf-8", errors="replace")
        self.assertIn("drift_dir_", src)
        self.assertIn("drift_inv_", src)

    def test_read_ini_overrides_the_legacy_slots_from_the_records(self):
        src = code(READINI)
        self.assertIn("DriftCorr%dir_cal(j, firstGas + i - 1)", src,
                      "without a record override the legacy tags cap drift "
                      "at the four historical slots")
        self.assertIn("DriftCorr%inv_cal(j, firstGas + i - 1)", src)

    def test_the_polynomials_come_only_from_the_records(self):
        """The flat drift_dir_<species>_<k> block is gone with the rest of the
        legacy tag layer; a 5.0.0 project states its coefficients per record.

        This used to assert the opposite - that the flat tables at 301 and 329
        *were* read - because CH4's and the fourth gas's had been declared and
        then never consulted, so those two gases were silently never
        drift-corrected. The records cover every gas, so the fix outlived the
        tags it was made in.
        """
        src = code(READINI)
        self.assertNotIn("301 + (i - co2) * 7", src)
        self.assertNotIn("329 + (i - co2) * 7", src)
        self.assertIn("rpGasOriginN + (i - 1) * rpGasLeapN + 11 + j", src)
        self.assertIn("rpGasOriginN + (i - 1) * rpGasLeapN + 18 + j", src)
        self.assertIn("DriftCorr%dir_cal = error", src,
                      "a gas with no polynomial must stay at error, which is "
                      "what excludes it from the correction")


class FixturesExerciseBothDirections(unittest.TestCase):
    def test_the_fixtures_exist(self):
        """There was no drift fixture at all, which is why the whole
        subsystem could carry a debug write to unit 987 and two unread
        coefficient sets without anything noticing."""
        for rel in FIXTURES:
            path = ROOT / rel
            self.assertTrue(path.exists(), "%s must exist" % rel)
            text = path.read_text(encoding="utf-8", errors="replace")
            self.assertRegex(text, r"(?m)^drift_method=1$")
            self.assertRegex(text, r"(?m)^use_dyn_md_file=1$")

    def test_the_forward_fixture_drifts_a_gas_past_the_fourth_slot(self):
        """base_drift is the backward-compatibility proof - four gases, legacy
        column names, byte-identical output. base_drift_n_gas is the forward
        one: it must name a reference column no legacy spelling could reach."""
        text = (ROOT / "tests/regression/base_drift_n_gas.eddyflow").read_text(
            encoding="utf-8", errors="replace")
        self.assertRegex(text, r"(?m)^gas_5_drift_inv_1=")
        self.assertRegex(text, r"(?m)^gas_6_drift_inv_1=")
        md = (ROOT / "tests/regression/base_drift_n_gas.metadata").read_text(
            encoding="utf-8", errors="replace")
        self.assertIn("n2o_ref", md)
        self.assertIn("co2_2_ref", md,
                      "a second record of the same species is addressed by "
                      "the _2 suffix FullOutputGasTags gives it")


if __name__ == "__main__":
    unittest.main()
