"""Analysers are addressed by gas slot, not by a second role numbering.

ExType%instr was indexed by `sonic, ico2, ih2o, ich4, igas4` - a five-wide
numbering of the *same* analysers the gas slots already number. read_ex_record
bridged the two with `igas = ico2 + (gas - co2)` and mirrored the results into
the slot-indexed gas_instr afterwards.

Two numberings for one thing cost three ways:

  - anything reached by role stopped at four gases, which is why fluxes1's
    oxygen correction only ever corrected one hygrometer;
  - the ex record converted the first four analysers' units in a loop over the
    roles and the rest in a duplicate of that arithmetic at the read site;
  - the mirror's real job turned out to be *ordering*. The ex record carries
    analyser metadata twice - the GA_* columns in metadata units and a
    self-describing block in SI - and the mirror ran after both, quietly
    restoring the GA_* values for slots five to eight. Remove it without
    saying which block owns which slots and the SI block gets scaled a second
    time; measured on base_n_gas, every separation came out a hundredfold
    small.

The anemometer keeps its slot in instr, because it is not a gas.
"""

from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]

TYPEDEF = "src/src_common/m_typedef.f90"
EXREAD = "src/src_common/read_ex_record.f90"

ROLES = ("ico2", "ih2o", "ich4", "igas4")


def code(path):
    return "\n".join(ln for ln in (ROOT / path).read_text(encoding="utf-8",
                                                          errors="replace").splitlines()
                     if not ln.lstrip().startswith("!"))


class TheRoleNumberingIsGone(unittest.TestCase):
    def test_the_gas_roles_are_not_declared(self):
        src = code(TYPEDEF)
        for role in ROLES:
            self.assertNotIn("parameter :: %s" % role, src,
                             "%s numbers an analyser a second time" % role)

    def test_the_anemometer_keeps_its_slot(self):
        self.assertIn("parameter :: sonic", code(TYPEDEF),
                      "the anemometer is not a gas and stays in instr")

    def test_nothing_addresses_an_analyser_by_role(self):
        for path in sorted(ROOT.glob("src/src_*/*.f90")):
            src = code(path.relative_to(ROOT).as_posix())
            for role in ROLES:
                self.assertNotIn("instr(%s" % role, src,
                                 "%s reaches an analyser by role" % path.name)

    def test_the_bridge_and_the_mirror_are_gone(self):
        src = code(EXREAD)
        self.assertNotIn("ico2 + (gas - histCO2)", src,
                         "the two numberings no longer need bridging")
        self.assertNotIn("lEx%gas_instr(igas - ico2 + firstGas)", src,
                         "nothing to mirror when there is one array")


class TheTwoMetadataBlocksHaveAnOwner(unittest.TestCase):
    """The GA_* columns own co2..gas4 and the self-describing block owns the
    rest. The mirror used to enforce that by running last; stated explicitly
    now, because the two blocks are in different units."""

    def test_units_are_converted_at_the_read(self):
        src = code(EXREAD)
        self.assertIn("merge(instr_nsep * 1d-2, instr_nsep", src,
                      "converting after the SI block would scale it twice")

    def test_the_error_code_survives_conversion(self):
        """A gas with no analyser carries the error code in every field.
        Guarding the assignment instead of the arithmetic leaves the field at
        its initialised zero, and a zero separation reads as a measurement -
        base_rec's CH4, which has no analyser, reported 0.00000 separations."""
        src = code(EXREAD)
        self.assertIn("instr_nsep /= error)", src)
        self.assertNotIn("if (instr_nsep /= error) lEx%gas_instr", src)

    def test_the_si_block_covers_every_slot_it_names(self):
        """The self-describing analyser block carries its own slot number.

        It used to defer to the GA_* columns for the first four, which was
        right while GA_* was four fixed blocks and this one existed to carry
        the rest. GA_* is per gas now and both end in SI, so the boundary
        protected nothing.

        The writer's gate and this one must agree, because FCC re-emits the
        count of entries it *accepted*: a reader discarding four of eight
        would make FCC's header over-declare by four blocks.
        """
        src = code(EXREAD)
        self.assertNotIn("gas > histGas4 .and. gas >= firstGas", src)
        self.assertIn("if (gas >= firstGas .and. gas <= lastGas) then", src)
        writer = code("src/src_rp/init_fluxnet_file_rp.f90")
        self.assertNotIn("if (k > 4) then", writer,
                         "the writer still starts the block at the fifth "
                         "record; the reader now accepts all of them")


class CalibrationReferencesArePerGas(unittest.TestCase):
    """Gas4CalRefCol was one raw-column index, so the calibration reference
    could only ever rescale the fourth slot - a site calibrating any other gas
    had its fourth gas rescaled by an unrelated reference and the gas it meant
    to calibrate left alone."""

    def test_the_scalar_is_gone(self):
        for path in sorted(ROOT.glob("src/src_*/*.f90")):
            self.assertNotIn("Gas4CalRefCol",
                             code(path.relative_to(ROOT).as_posix()),
                             "%s still holds one calibration column" % path.name)

    def test_the_column_map_is_per_slot(self):
        self.assertIn("GasCalRefCol(GHGNumVar)",
                      code("src/src_common/m_common_global_var.f90"))

    def test_the_user_column_remembers_which_gas_it_calibrates(self):
        """The 'cal-ref' marker says only *that* a column is a reference."""
        self.assertIn("UserCalRefSlot(MaxUserVar)",
                      code("src/src_common/m_common_global_var.f90"))
        src = code("src/src_common/calibrate_gas.f90")
        self.assertIn("UserCalRefSlot(j)", src)
        self.assertNotIn("Set(:, gas4)", src)

    def test_the_routine_is_named_for_what_it_does(self):
        self.assertTrue((ROOT / "src/src_common/calibrate_gas.f90").exists())
        self.assertFalse((ROOT / "src/src_common/calibrate_gas4.f90").exists())
        self.assertIn("subroutine CalibrateGases",
                      code("src/src_common/calibrate_gas.f90"))


if __name__ == "__main__":
    unittest.main()
