"""Static checks for the closed-path spectroscopic correction.

Peltola et al. (2014): water vapour broadens a laser analyser's absorption
lines, so what it reports scales as ``1 + a*chi_q + b*chi_q^2`` and the
correction is a division by that. Applied point by point after Chen et al.
(2010), on the raw series.

Four things about it are worth pinning.

**It must stay off, and be the identity when on but undeclared.** Two gates:
``spectro_meth`` defaults to none, and both coefficients default to zero,
which makes the polynomial one.

**It must run before the mixing-ratio conversion and before the lag shift.**
The correction divides a gas by the water its own analyser read *at the same
sample*. That pairing exists only while the series are as measured; after
``TimeLagHandle`` shifts each column by its own lag, row i of the gas and row
i of the water are no longer one moment.

**It must not be gated on WPL.** The bias is in what the analyser reported,
whatever units it reported in and whether or not a density correction was
asked for. ``PointByPointToMixingRatio`` next door *is* gated on WPL, so the
two calls sit adjacent with different conditions and are easy to conflate.

**The metadata field must not collide with the calibration a and b.** Every
column already carries ``col_N_a_value`` and ``col_N_b_value``, which are the
linear gain and offset and have nothing to do with this.
"""

import re
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / "src"


def read(rel):
    return (SRC / rel).read_text(encoding="utf-8", errors="replace")


SPECTRO = read("src_rp/spectroscopic_closed_path.f90")
MAIN = read("src_rp/eddyflow-rp_main.f90")
READER = read("src_rp/read_ini_rp.f90")
TAGS = read("src_rp/m_rp_global_var.f90")
COMMON = read("src_common/m_common_global_var.f90")
METAREAD = read("src_common/read_metadata_file.f90")
TYPEDEF = read("src_common/m_typedef.f90")

#: project keys and their SNTags slots
KEYS = {58: "spectro_meth", 59: "spectro_water"}


class TheSettingsReachTheEngine(unittest.TestCase):

    def test_each_key_holds_its_slot(self):
        for index, label in KEYS.items():
            self.assertRegex(
                TAGS, r"SNTags\(%d\)%%Label\s*/\s*'%s'\s*/" % (index, label))

    def test_each_key_is_read_under_its_found_guard(self):
        for index in KEYS:
            self.assertIn("SNTagFound(%d)" % index, READER)

    def test_the_method_defaults_to_none(self):
        self.assertRegex(READER, r"RPSetup%spectro_meth\s*=\s*'none'")

    def test_the_water_channel_defaults_to_off(self):
        self.assertRegex(READER, r"RPSetup%spectro_water\s*=\s*\.false\.")


class TheCoefficientsAreDeclaredPerColumn(unittest.TestCase):

    def test_the_metadata_tags_exist(self):
        for suffix in ("spectro_a", "spectro_b"):
            self.assertRegex(
                COMMON, r"ANTags\(\d+\)%%Label / 'col_1_%s'" % suffix,
                "col_1_%s is missing from the metadata tag table" % suffix)

    def test_they_do_not_collide_with_the_calibration_gain_and_offset(self):
        #> col_N_a_value and col_N_b_value are the linear calibration. Both
        #> pairs live on the same column and neither name may be a substring
        #> of the other, which is also what test_ini_tag_collisions_static
        #> enforces across the whole table.
        for old, new in (("col_1_a_value", "col_1_spectro_a"),
                         ("col_1_b_value", "col_1_spectro_b")):
            self.assertNotIn(old, new)
            self.assertNotIn(new, old)

    def test_the_column_type_carries_them_separately_from_a_and_b(self):
        for field in ("spectro_a", "spectro_b"):
            self.assertRegex(
                TYPEDEF, r"real\(kind = dbl\) :: %s" % field,
                "ColType has no %s" % field)

    def test_the_reader_defaults_them_to_zero_when_absent(self):
        #> Absent must mean zero, not whatever the shared tag array last held:
        #> the identity, so a metadata file written before the key existed
        #> declares correctly that there is nothing to remove.
        self.assertRegex(METAREAD, r"LocCol\(i\)%spectro_a = 0d0")
        self.assertRegex(METAREAD, r"LocCol\(i\)%spectro_b = 0d0")
        self.assertIn("ANTagFound(init_an_col + i*leap_an_col + 9)", METAREAD)
        self.assertIn("ANTagFound(init_an_col + i*leap_an_col + 10)", METAREAD)

    def test_the_column_stride_grew_with_them(self):
        #> Two numeric fields were appended to each column's run, so the leap
        #> is 11 where it was 9. Getting this wrong does not raise - it reads
        #> every column's values from the wrong offsets.
        self.assertRegex(METAREAD, r"leap_an_col = 11")


class TheCorrectionIsInertUnlessAsked(unittest.TestCase):

    def test_it_returns_before_touching_anything(self):
        self.assertRegex(
            SPECTRO, r"if \(RPSetup%spectro_meth /= 'chen_10'\) return")

    def test_a_column_with_no_coefficients_is_skipped(self):
        self.assertRegex(
            SPECTRO, r"if \(coef_a == 0d0 \.and\. coef_b == 0d0\) cycle")

    def test_the_whole_pass_is_skipped_when_nothing_declares_one(self):
        self.assertIn("if (.not. anyCoefficients) return", SPECTRO)

    def test_the_water_channel_needs_its_own_switch(self):
        self.assertRegex(
            SPECTRO,
            r"if \(GasSlotIsWater\(gas\) \.and\. \.not\. RPSetup%spectro_water\) cycle",
            "the water channel is corrected without its own opt-in; that "
            "form is not part of the published result",
        )


class TheCoefficientConventionIsSpeltOut(unittest.TestCase):
    """EddyUH's coefficients are not interchangeable with these.

    Both packages divide by ``1 + a*chi_q + b*chi_q^2`` with the water in
    mol/mol - identical formulae - but EddyUH folds the dilution into the same
    polynomial: ``dilucorr.m`` uses ``a = -1, b = 0`` for the pure-dilution
    case. EddyFlow corrects the density separately, so its identity is
    ``a = b = 0``.

    A published value entered unconverted therefore counts the dilution twice,
    and nothing about the two formulae makes that visible. The mapping is
    ``a_here = a_EddyUH + 1``. This has to stay written down where someone
    typing a number will see it.
    """

    def test_the_docstring_states_the_mapping(self):
        self.assertIn("a_here = a_EddyUH + 1", SPECTRO)

    def test_it_says_the_coefficients_are_not_eddyuhs(self):
        self.assertIn("SPECTROSCOPIC ONLY", SPECTRO)
        self.assertIn("count the dilution a second time", SPECTRO)

    def test_it_records_that_the_crosstalk_routine_is_dead(self):
        #> Worth keeping: the file exists, is named as though it were a
        #> separate correction, and reads like one. The next person to compare
        #> the two packages will find it again.
        self.assertIn("EddyUH_ctc_CP.m", SPECTRO)
        self.assertIn("unreachable", SPECTRO)


class TheCallSitePairsGasWithItsOwnWater(unittest.TestCase):

    def test_it_runs_before_the_mixing_ratio_conversion(self):
        spectro = MAIN.index("call SpectroscopicClosedPath(")
        #> The conversion that follows it in the main loop. There is an
        #> earlier one, in the lag-optimisation pre-pass, which this
        #> deliberately does not pair with - see the test below.
        pbp = MAIN.index("call PointByPointToMixingRatio(", spectro)
        self.assertLess(
            spectro, pbp,
            "the correction runs after the conversion, so it would divide "
            "already-converted values by a water that has also moved",
        )

    def test_the_lag_optimisation_pre_pass_is_left_uncorrected(self):
        #> That pre-pass exists to find time lags, and it works on the rawest
        #> series it can: the CHANGELOG records PWB detection being made
        #> always pre-WPL for the same reason. A sub-percent multiplicative
        #> factor that varies slowly with humidity cannot move the lag at
        #> which a cross-covariance peaks, so correcting there would buy
        #> nothing and cost a second pass over every period twice over.
        #>
        #> One call only, and it is the one in the main loop.
        self.assertEqual(
            MAIN.count("call SpectroscopicClosedPath("), 1,
            "the correction is called more than once; the lag-optimisation "
            "pre-pass is meant to see the series as measured",
        )
        first_conversion = MAIN.index("call PointByPointToMixingRatio(")
        spectro = MAIN.index("call SpectroscopicClosedPath(")
        self.assertGreater(
            spectro, first_conversion,
            "the correction has moved into the lag-optimisation pre-pass",
        )

    def test_it_runs_before_the_time_lag_shift(self):
        spectro = MAIN.index("call SpectroscopicClosedPath(")
        tlag = MAIN.index("call TimeLagHandle(Meth%tlag")
        self.assertLess(
            spectro, tlag,
            "the correction runs after the lag shift, when row i of a gas "
            "and row i of its water are no longer the same moment",
        )

    def test_it_is_not_gated_on_wpl(self):
        #> The line before it is, which is exactly why this is worth pinning.
        window = MAIN[MAIN.index("call SpectroscopicClosedPath(") - 400:
                      MAIN.index("call SpectroscopicClosedPath(")]
        tail = window.rsplit("!>", 1)[-1]
        self.assertNotIn(
            "EddyFlowProj%wpl", tail,
            "the spectroscopic correction has been put behind the WPL "
            "switch; the bias is in what the analyser reported and is there "
            "whether or not a density correction was asked for",
        )

    def test_the_water_is_read_before_any_column_is_divided(self):
        #> Two passes. A hygrometer may itself be corrected, and every gas has
        #> to be divided by the water that was read, not one already moved.
        build = SPECTRO.index("call WaterMoleFraction(")
        divide = SPECTRO.index("call DivideBySensitivity(")
        self.assertLess(build, divide)

    def test_the_water_must_be_on_the_same_analyser(self):
        self.assertRegex(
            SPECTRO,
            r"if \(E2Col\(gas\)%instr%model /= E2Col\(msl\)%instr%model\) cycle",
            "a gas may now be corrected by a hygrometer on another "
            "instrument, which sampled different air",
        )

    def test_open_path_columns_are_declined(self):
        self.assertGreaterEqual(
            SPECTRO.count("path_type /= 'closed') cycle"), 2,
            "the correction is a property of the cell the sample passed "
            "through; an open path has none",
        )


class TheDivisionIsGuarded(unittest.TestCase):

    def test_a_non_positive_sensitivity_is_refused(self):
        self.assertRegex(
            SPECTRO, r"d\(:\) >= min_sensitivity",
            "dividing by a near-zero sensitivity would turn a "
            "sub-percent bias into a wild one",
        )

    def test_a_sample_with_no_water_is_left_alone_rather_than_voided(self):
        #> Leaving it uncorrected costs a fraction of a per cent on a handful
        #> of samples; voiding it puts gaps into a series that the despiking
        #> and the spectra have already been told is continuous.
        self.assertIn("col(:) /= error .and. chi(:) /= error", SPECTRO)
        self.assertNotRegex(
            SPECTRO, r"col\(:\) = error",
            "samples are being voided where the water is missing")


if __name__ == "__main__":
    unittest.main()
