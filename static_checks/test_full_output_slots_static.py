"""The full output's header and rows walk one shared list of gas slots.

This is what remains of test_fix_out_format_pairing_static once the flag it
was named for is gone. Its subject died; its defence did not.

The defect it recorded: both row writers looped `firstGas, lastGas` while the
header described a smaller set. While lastGas was 8 the two agreed. At 68 the
row carried sixty phantom gas blocks - measured on base_rec, header 194 fields
against rows 1094. Nine hundred extra, sixty absent gases times fifteen
fields, and every field after the gas block shifted. FCC parses the ex record
by comma count, so it consumed the row without complaining.

Both sides call FullOutputGasSlots. What this defends is that they keep doing
so, in RP and in its FCC twin: a block added to one side and not the other
reintroduces the same silent shift. The helper no longer has a fixed arm - the
full output covers every configured gas now - so the list is one loop, and
what matters is that nobody goes around it.
"""

from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]

HELPER = "src/src_common/gas_slot_resolution.f90"
RP_HEADER = "src/src_rp/init_outfiles_rp.f90"
RP_WRITER = "src/src_rp/write_out_full.f90"
FCC_HEADER = "src/src_fcc/init_out_files.f90"
FCC_WRITER = "src/src_fcc/write_out_full_fcc.f90"


def code(path):
    return "\n".join(ln for ln in (ROOT / path).read_text(encoding="utf-8",
                                                          errors="replace").splitlines()
                     if not ln.lstrip().startswith("!"))


class OneListServesBothSides(unittest.TestCase):
    def test_the_helper_exists(self):
        self.assertIn("subroutine FullOutputGasSlots", code(HELPER))

    def test_all_four_files_call_it(self):
        for path in (RP_HEADER, RP_WRITER, FCC_HEADER, FCC_WRITER):
            self.assertIn("call FullOutputGasSlots", code(path),
                          "%s must walk the shared list, not its own idea of "
                          "which slots the file carries" % path)

    def test_no_writer_still_loops_the_whole_gas_block(self):
        """`firstGas, lastGas` is 64 slots; the header describes the ones the
        project configured and the row must agree."""
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


if __name__ == "__main__":
    unittest.main()
