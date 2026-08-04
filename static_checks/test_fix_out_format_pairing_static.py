"""The full output's header and rows walk one shared list of gas slots.

Same failure as the statistics files, in the other output file and one flag
away. The full output has two header branches. The dynamic branch names one
block per present gas. The `fix_out_format` branch is a literal naming exactly
co2, h2o, ch4 and the fourth slot - that is the point of the flag, which
promises the fixed EddyPro 7.x column set.

Both row writers, though, looped `firstGas, lastGas`, and their
`elseif (fix_out_format)` arms emit placeholder fields for slots holding *no
gas at all*. While lastGas was 8 the two agreed. At 68 the row carries sixty
phantom gas blocks. Measured on base_rec with fix_out_format=1, before the
fix: header 194 fields, rows 1094. Nine hundred extra - sixty absent gases
times fifteen fields - and every field after the gas block is shifted. FCC
parses the ex record by comma count, so it consumes the row without
complaining.

No fixture set the flag, which is why it went unseen. base_rec_fix and
base_5gas_fix do now.

Both sides call FullOutputGasSlots. What this defends is that they keep doing
so, in RP and in its FCC twin: a block added to one side and not the other
reintroduces the same silent shift.
"""

from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]

HELPER = "src/src_common/gas_slot_resolution.f90"
RP_HEADER = "src/src_rp/init_outfiles_rp.f90"
RP_WRITER = "src/src_rp/write_out_full.f90"
FCC_HEADER = "src/src_fcc/init_out_files.f90"
FCC_WRITER = "src/src_fcc/write_out_full_fcc.f90"

FIXTURES = ("tests/regression/base_rec_fix.eddyflow",
            "tests/regression/base_5gas_fix.eddyflow")


def code(path):
    return "\n".join(ln for ln in (ROOT / path).read_text(encoding="utf-8",
                                                          errors="replace").splitlines()
                     if not ln.lstrip().startswith("!"))


class OneListServesBothSides(unittest.TestCase):
    def test_the_helper_exists_and_knows_about_the_flag(self):
        src = code(HELPER)
        self.assertIn("subroutine FullOutputGasSlots", src)
        body = src.split("subroutine FullOutputGasSlots")[1].split("end subroutine")[0]
        self.assertIn("EddyFlowProj%fix_out_format", body,
                      "the fixed format's four-block promise is what the list "
                      "exists to keep, so the helper has to test the flag")

    def test_the_fixed_arm_is_exactly_four_slots(self):
        """Widening it would break the compatibility the flag provides."""
        src = code(HELPER)
        body = src.split("subroutine FullOutputGasSlots")[1].split("end subroutine")[0]
        fixed = body.split("fix_out_format")[1].split("return")[0]
        self.assertIn("do gas = histGas1, histGas4", fixed,
                      "the fixed format names co2, h2o, ch4 and the fourth "
                      "slot, present or not - the row fills absent ones with "
                      "the error label")

    def test_all_four_files_call_it(self):
        for path in (RP_HEADER, RP_WRITER, FCC_HEADER, FCC_WRITER):
            self.assertIn("call FullOutputGasSlots", code(path),
                          "%s must walk the shared list, not its own idea of "
                          "which slots the file carries" % path)

    def test_no_writer_still_loops_the_whole_gas_block(self):
        """`firstGas, lastGas` is 64 slots; the fixed header describes four."""
        for path in (RP_WRITER, FCC_WRITER):
            src = code(path)
            self.assertNotIn("do gas = firstGas, lastGas", src,
                             "%s must take its gas slots from the list" % path)
            self.assertNotIn("do var = firstGas, lastGas", src,
                             "%s must take its gas slots from the list" % path)

    def test_the_water_correction_factor_follows_the_resolved_slot(self):
        """LE's un_LE/LE_scf pair is gated on PrimaryWaterOutSlot, so its
        correction factor cannot come from the historical sixth slot - that is
        water only when record two happens to hold it."""
        for path in (RP_WRITER, FCC_WRITER):
            self.assertNotIn("BPCF%of(w_h2o)", code(path),
                             "%s reads a slot number and means water" % path)


class FixturesTurnTheFlagOn(unittest.TestCase):
    def test_the_fixed_format_is_exercised(self):
        """The defect survived because every fixture left fix_out_format at 0,
        so that whole header branch was never written and never compared."""
        for rel in FIXTURES:
            path = ROOT / rel
            self.assertTrue(path.exists(),
                            "%s must exist so the fixed format is exercised "
                            "by the harness" % rel)
            text = path.read_text(encoding="utf-8", errors="replace")
            self.assertRegex(text, r"(?m)^fix_out_format=1$",
                             "%s must switch the fixed format on" % rel)

    def test_one_fixture_configures_more_gases_than_the_format_holds(self):
        """Four gases prove the format still works; five prove the gases past
        it are dropped rather than shifting the row."""
        text = (ROOT / "tests/regression/base_5gas_fix.eddyflow").read_text(
            encoding="utf-8", errors="replace")
        self.assertRegex(text, r"(?m)^gas_num=5$")


if __name__ == "__main__":
    unittest.main()
