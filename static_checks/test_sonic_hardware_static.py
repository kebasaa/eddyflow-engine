"""The two sonic hardware corrections, and the four things that make them safe.

Both act on the raw wind before any rotation - a Metek USA-1 head correction
(``Functions_Library/METEK_HC.m``) and an inclinometer tilt correction
(``EC_Software_Common/EddyUH_tiltangle.m``). Neither affects CH-LAE, which runs
a Gill HS with ``tiltcorr = 0``, so nothing in the regression fixtures
exercises them; these assertions and the synthetic driver they describe are the
whole of the coverage.

Four properties carry the design:

1. **Off is a bare return.** Both routines leave on ``'none'`` before touching
   anything, and both call sites sit behind that, so byte-identity with the
   keys absent is structural rather than incidental.
2. **The angles come through the custom-column machinery**, not a new record
   family. ``UserSet`` is a local of the main program - so it is *passed*,
   not reached for. That is what the first attempt got wrong, and it failed
   to compile rather than failing quietly, which is the good case.

   The cost is that ``DefineUserSet`` runs after both pre-passes, so the
   inclinometer correction reaches the flux loop only, while the head
   correction reaches all three. The routine's own note records why that is
   left standing.
3. **The Metek tables are not shipped.** They are Metek GmbH's measurements,
   which EddyUH redistributes under the University of Helsinki's agreement.
   Nothing in this repository may carry them, and a missing table declines the
   correction for the whole run instead of correcting some periods and not
   others.
4. **The swinging term is EddyUH's, dot product and all.** See below.

Worth recording, because it is a defect being deliberately reproduced:
**EddyUH's swinging correction is a scalar added to all three wind
components.** ``EddyUH_tiltangle.m:104`` writes ``V_true(i,:) + omega*T*L``
with ``omega`` 1x3, ``T`` 3x3 and ``L = [-1.5 -1.5 -1.5]'`` 3x1, so the product
is one number. The velocity of a point on a rotating body is ``omega`` CROSS
``(T L)``, a vector with three different components. The units survive - radians
per second times metres is metres per second - which is why it is easy to miss.
It is reproduced as written because the option exists to reproduce EddyUH's
numbers, and a silently corrected version would reproduce nothing.
"""

import re
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TILT = ROOT / "src" / "src_rp" / "inclinometer_tilt.f90"
HEAD = ROOT / "src" / "src_rp" / "metek_head_correction.f90"
RP_MAIN = ROOT / "src" / "src_rp" / "eddyflow-rp_main.f90"
READ_INI = ROOT / "src" / "src_rp" / "read_ini_rp.f90"
TYPEDEF = ROOT / "src" / "src_common" / "m_typedef.f90"
TAGS = ROOT / "src" / "src_rp" / "m_rp_global_var.f90"
GEN = ROOT / "prj" / "gen_project_tags.py"
EDDYUH = (ROOT.parent / "EddyUH_testing" / "EddyUH" / "EddyUH_1.7b_COS")
EDDYUH_TILT = EDDYUH / "EC_Software_Common" / "EddyUH_tiltangle.m"
EDDYUH_HEAD = EDDYUH / "Functions_Library" / "METEK_HC.m"


def read(path):
    return path.read_text(encoding="utf-8", errors="replace")


TILT_SRC = read(TILT)
HEAD_SRC = read(HEAD)
MAIN_SRC = read(RP_MAIN)
INI_SRC = read(READ_INI)


class OffIsABareReturn(unittest.TestCase):
    """The one property that keeps every existing project byte-identical."""

    def test_each_routine_leaves_before_it_touches_the_wind(self):
        for name, src, key in (("tilt", TILT_SRC, "tilt_sensor_meth"),
                               ("head", HEAD_SRC, "head_corr_meth")):
            guard = "if (RPSetup%%%s == 'none') return" % key
            self.assertIn(guard, src, name)
            #> Before the first write to Set, not merely somewhere in the file.
            first_write = re.search(r"^\s*Set\(i, u\) =", src, re.M)
            self.assertIsNotNone(first_write, name)
            self.assertLess(src.index(guard), first_write.start(), name)

    def test_the_defaults_are_off(self):
        self.assertIn("RPSetup%tilt_sensor_meth = 'none'", INI_SRC)
        self.assertIn("RPSetup%head_corr_meth = 'none'", INI_SRC)

    def test_every_key_is_read_only_when_the_project_carries_it(self):
        #> A key absent from the file must leave the default standing rather
        #> than read whatever the tag array happened to hold.
        for tag in ("SCTagFound(31)", "SCTagFound(32)", "SCTagFound(33)",
                    "SNTagFound(117)", "SNTagFound(118)", "SNTagFound(119)",
                    "SNTagFound(120)", "SNTagFound(121)"):
            self.assertIn(tag, INI_SRC, tag)


class TheCallSitesRunBeforeRotation(unittest.TestCase):
    """RP reads the raw data three times: once to optimise the time lags,
    once to assess the planar fit, and once to compute the fluxes. A
    correction applied to only one of them works a lag or a plane out from a
    wind the fluxes never see."""

    def test_the_head_correction_runs_in_all_three_passes(self):
        #> As many times as the angle-of-attack calibration beside it.
        self.assertEqual(MAIN_SRC.count("call MetekHeadCorrection("),
                         MAIN_SRC.count("call AoaCalibration("))
        self.assertEqual(MAIN_SRC.count("call MetekHeadCorrection("), 3)

    def test_each_one_sits_after_the_other_hardware_corrections(self):
        aoa = [m.start() for m in re.finditer(r"call AoaCalibration\(",
                                              MAIN_SRC)]
        head = [m.start() for m in re.finditer(r"call MetekHeadCorrection\(",
                                               MAIN_SRC)]
        self.assertEqual(len(aoa), len(head))
        for a, h in zip(aoa, head):
            self.assertLess(a, h)

    def test_every_rotation_is_preceded_by_a_head_correction(self):
        rot = [m.start() for m in re.finditer(r"call TiltCorrection\(",
                                              MAIN_SRC)]
        head = [m.start() for m in re.finditer(r"call MetekHeadCorrection\(",
                                               MAIN_SRC)]
        #> Two of the three passes rotate; the planar-fit assessment only
        #> collects means. Each rotation must have a head correction between
        #> it and the previous one.
        for r in rot:
            self.assertTrue(any(h < r for h in head))
            nearest = max(h for h in head if h < r)
            self.assertNotIn("call TiltCorrection(", MAIN_SRC[nearest:r])

    def test_the_inclinometer_runs_in_the_flux_loop_only(self):
        #> Not an oversight - see the note in the routine. UserSet does not
        #> exist during either pre-pass, and this asserts the limitation is
        #> the documented one rather than a fourth call site drifting in.
        self.assertEqual(MAIN_SRC.count("call InclinometerTilt("), 1)
        self.assertIn("THE FLUX LOOP ONLY", TILT_SRC)

    def test_the_head_correction_feeds_the_inclinometer_one(self):
        tilt = MAIN_SRC.index("call InclinometerTilt(")
        head = max(m.start() for m in
                   re.finditer(r"call MetekHeadCorrection\(", MAIN_SRC)
                   if m.start() < tilt)
        self.assertLess(head, tilt)

    def test_the_inclinometer_precedes_the_rotation_that_follows_it(self):
        tilt = MAIN_SRC.index("call InclinometerTilt(")
        rot = [m.start() for m in re.finditer(r"call TiltCorrection\(",
                                              MAIN_SRC) if m.start() > tilt]
        self.assertTrue(rot, "nothing rotates after the correction at all")

    def test_the_angles_are_defined_before_they_are_read(self):
        #> UserSet is filled by DefineUserSet. Calling the tilt correction
        #> ahead of it would read an unallocated array.
        self.assertLess(MAIN_SRC.index("call DefineUserSet("),
                        MAIN_SRC.index("call InclinometerTilt("))


class TheAnglesArePassedNotReachedFor(unittest.TestCase):
    """``UserSet`` is a local of the main program, ``UserCol`` is a global."""

    def test_userset_is_an_argument(self):
        sig = re.search(r"subroutine InclinometerTilt\((.*?)\)", TILT_SRC, re.S)
        self.assertIsNotNone(sig)
        self.assertIn("UserSet", sig.group(1))

    def test_the_call_site_guards_a_project_with_no_extra_columns(self):
        #> UserSet is unallocated when NumUserVar is zero, and passing an
        #> unallocated array to an explicit-shape dummy is undefined.
        #> Anchored on the tilt call, not on the head correction: that one
        #> appears three times and the first is in a pre-pass.
        at = MAIN_SRC.index("call InclinometerTilt(")
        block = MAIN_SRC[MAIN_SRC.rindex("\n", 0, at) - 400:]
        block = block[:block.index("call TiltCorrection(")]
        self.assertIn("NumUserVar > 0", block)
        self.assertIn("allocated(UserSet)", block)

    def test_the_channels_are_matched_by_name(self):
        for name in ("'theta'", "'phi'", "'psi'"):
            self.assertIn("call AngleChannel(%s" % name, TILT_SRC)

    def test_a_reading_past_full_scale_is_clamped(self):
        #> asin is undefined outside [-1, 1]; a NaN here would spread into
        #> every wind component of that sample and then into the covariances.
        self.assertIn("if (s > 1d0) s = 1d0", TILT_SRC)
        self.assertIn("if (s < -1d0) s = -1d0", TILT_SRC)


class PsiIsAlwaysZero(unittest.TestCase):
    """EddyUH reads a psi channel and then throws it away."""

    def test_the_angle_is_forced_to_zero_after_it_is_read(self):
        i = TILT_SRC.index("call AngleChannel('psi'")
        rest = TILT_SRC[i:]
        self.assertIn("psi = 0d0", rest[:rest.index("if (n_found == 0)")])

    @unittest.skipUnless(EDDYUH_TILT.is_file(), "EddyUH sources not present")
    def test_that_is_what_eddyuh_does(self):
        src = read(EDDYUH_TILT)
        self.assertIn("not measured", src)
        self.assertIn("psi = zeros(length(data.u),1)", src)


class TheSwingingTermIsEddyUHs(unittest.TestCase):
    """A dot product where the physics wants a cross product."""

    def test_ours_is_a_dot_product_added_to_all_three_components(self):
        self.assertIn("dot_product(matmul(omega, t), RPSetup%tilt_arm)",
                      TILT_SRC)
        #> The elementwise product this started as is a different answer
        #> whenever the arm's components differ - and even with EddyUH's own
        #> equal arm, it adds three numbers where EddyUH adds one.
        self.assertNotIn("matmul(omega, t) * RPSetup%tilt_arm", TILT_SRC)

    def test_the_note_says_it_is_wrong_as_physics(self):
        #> A reader who spots the dot product must find, in the file, that it
        #> was noticed and kept on purpose. Otherwise this reads as a bug.
        self.assertIn("CROSS", TILT_SRC)
        self.assertIn("EddyUH_tiltangle.m:104", TILT_SRC)

    @unittest.skipUnless(EDDYUH_TILT.is_file(), "EddyUH sources not present")
    def test_eddyuhs_arm_is_a_column_so_the_product_is_scalar(self):
        src = read(EDDYUH_TILT)
        self.assertRegex(src, r"L = \[-1\.5 -1\.5 -1\.5\]'")
        self.assertIn("omega*T*L", src)


class TheMetekTablesAreNotShipped(unittest.TestCase):

    def test_no_table_file_is_in_the_repository(self):
        for name in ("phicorr.dat", "ucorr.dat", "alphacorr.dat"):
            found = list(ROOT.rglob(name))
            self.assertEqual(found, [],
                             "%s is Metek GmbH data and cannot be "
                             "redistributed here" % name)

    def test_the_directory_is_a_setting(self):
        self.assertIn("RPSetup%head_corr_dir", HEAD_SRC)
        self.assertIn("RPSetup%head_corr_dir = ''", INI_SRC)

    def test_a_missing_table_declines_for_the_whole_run(self):
        #> Latched, not retried. Correcting some periods and not others would
        #> put two different winds in one output file.
        self.assertIn("logical, save :: refused", HEAD_SRC)
        self.assertIn("if (refused) return", HEAD_SRC)
        self.assertIn("refused = .true.", HEAD_SRC)

    def test_the_message_says_why_they_are_absent(self):
        self.assertIn("Metek GmbH data", HEAD_SRC)

    def test_the_tables_are_read_once_rather_than_once_a_period(self):
        self.assertIn("logical, save :: loaded", HEAD_SRC)
        self.assertIn("loaded = .true.", HEAD_SRC)


class TheMetekGridMatchesTheTables(unittest.TestCase):

    def test_twenty_rows_from_minus_fifty_in_fives(self):
        self.assertIn("integer, parameter :: nrows = 20", HEAD_SRC)
        self.assertIn("grd_first = -50d0", HEAD_SRC)
        self.assertIn("grd_step = 5d0", HEAD_SRC)

    def test_the_ends_are_held_rather_than_extrapolated(self):
        #> A ninth-harmonic fit extrapolated past the measured range produces
        #> confident nonsense.
        fn = re.search(r"subroutine Interpolate\(.*?end subroutine Interpolate",
                       HEAD_SRC, re.S)
        self.assertIsNotNone(fn)
        body = fn.group(0)
        self.assertIn("c = tab(1, :)", body)
        self.assertIn("c = tab(nrows, :)", body)
        #> EddyUH indexes Cf0(I+1) with I = fix((phi+50)/5)+1, which reaches
        #> row 21 at exactly +45 degrees. This guard is why ours does not.
        self.assertIn("if (idx >= nrows) then", body)

    @unittest.skipUnless(EDDYUH_HEAD.is_file(), "EddyUH sources not present")
    def test_the_harmonics_are_the_ones_eddyuh_evaluates(self):
        src = read(EDDYUH_HEAD)
        for h in ("3.*alfaM", "6.*alfaM", "9.*alfaM"):
            self.assertIn(h, src)
        for h in ("3d0 * az", "6d0 * az", "9d0 * az"):
            self.assertIn(h, HEAD_SRC)


class TheKeysAreRegistered(unittest.TestCase):

    def test_the_generator_owns_all_eight(self):
        gen = read(GEN)
        for key in ("tilt_sensor_meth", "tilt_sensor_v_g", "tilt_arm_x",
                    "tilt_arm_y", "tilt_arm_z", "tilt_lpf_s",
                    "head_corr_meth", "head_corr_dir"):
            self.assertIn(key, gen, key)

    def test_the_generated_table_holds_them_at_the_expected_slots(self):
        tags = read(TAGS)
        for slot, key in ((31, "tilt_sensor_meth"), (32, "head_corr_meth"),
                          (33, "head_corr_dir")):
            self.assertRegex(
                tags, r"SCTags\(%d\)%%Label\s*/\s*'%s'\s*/" % (slot, key))
        for slot, key in ((117, "tilt_sensor_v_g"), (118, "tilt_arm_x"),
                          (119, "tilt_arm_y"), (120, "tilt_arm_z"),
                          (121, "tilt_lpf_s")):
            self.assertRegex(
                tags, r"SNTags\(%d\)%%Label\s*/\s*'%s'\s*/" % (slot, key))

    def test_the_setup_type_carries_them(self):
        src = read(TYPEDEF)
        #> The declaration reads "type :: RPsetupType", so anchoring on
        #> "type RPsetupType" finds the END of the type, not its start - and
        #> the block then runs off the end of the file. And not the first
        #> "end type" either: RPsetupType has nested types of its own, and
        #> slicing at that one stops after five lines.
        end_at = src.index("end type RPsetupType")
        block = src[src.rindex("type :: RPsetupType", 0, end_at):end_at]
        for field in ("tilt_sensor_meth", "tilt_sensor_v_g", "tilt_arm",
                      "tilt_lpf_s", "head_corr_meth", "head_corr_dir"):
            self.assertIn(field, block, field)


class TheNumbersAreEddyUHsOwn(unittest.TestCase):
    """Switching one of these on must reproduce EddyUH, not something new."""

    def defaults_block(self):
        #> The DEFAULT assignment, not merely the literal somewhere in the
        #> file - the same literal appears again in the guard below that
        #> refuses a non-positive sensitivity, so an unscoped search passes
        #> even when the default itself has been changed. Found by injecting
        #> exactly that and watching this check not bite.
        start = INI_SRC.index("RPSetup%tilt_sensor_meth = 'none'")
        return INI_SRC[start:INI_SRC.index("if (SNTagFound(117))", start)]

    def test_the_sensitivity_and_the_arm(self):
        block = self.defaults_block()
        self.assertIn("RPSetup%tilt_sensor_v_g = 4d0", block)
        self.assertIn("RPSetup%tilt_arm = -1.5d0", block)
        self.assertIn("RPSetup%tilt_lpf_s = 0d0", block)

    @unittest.skipUnless(EDDYUH_TILT.is_file(), "EddyUH sources not present")
    def test_they_are_eddyuhs_literals(self):
        src = read(EDDYUH_TILT)
        self.assertRegex(src, r"sensitivity = 4;")
        self.assertRegex(src, r"L = \[-1\.5 -1\.5 -1\.5\]'")

    def test_a_sensitivity_of_zero_is_refused(self):
        #> It would divide the whole angle series by nothing.
        self.assertIn("if (RPSetup%tilt_sensor_v_g <= 0d0)", INI_SRC)

    def test_a_negative_filter_length_is_refused(self):
        self.assertIn("if (RPSetup%tilt_lpf_s < 0d0)", INI_SRC)


class TheMetekFrameIsFlippedBothWays(unittest.TestCase):
    """USA-1 reports a left-handed frame; the rotation needs a right-handed
    one, and the caller expects what it handed over."""

    def test_v_is_negated_before_and_after(self):
        self.assertEqual(TILT_SRC.count("if (metek) v_obs(2) = -v_obs(2)"), 1)
        self.assertEqual(TILT_SRC.count("if (metek) v_true(2) = -v_true(2)"), 1)

    def test_the_model_is_what_decides(self):
        self.assertIn("index(E2Col(u)%instr%model, 'usa1')", TILT_SRC)


class ErrorSamplesAreLeftAlone(unittest.TestCase):

    def test_neither_routine_rewrites_a_missing_wind_component(self):
        for name, src in (("tilt", TILT_SRC), ("head", HEAD_SRC)):
            self.assertRegex(
                src,
                r"if \(Set\(i, u\) == error \.or\. Set\(i, v\) == error\s*"
                r"&?\s*!?\s*\.or\. Set\(i, w\) == error\) cycle",
                name)


if __name__ == "__main__":
    unittest.main()
