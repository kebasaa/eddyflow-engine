"""Which columns an analyser's diagnostic word may invalidate.

The list used to be spelled `case (co2:gas4, pi:pe)`, and both halves stopped
meaning what they said once the slots widened:

  - `co2:gas4` is the first four gas slots, so a gas past the fourth kept every
    record its analyser's diagnostic rejected.
  - `pi:pe` was instrument 1's cell pressure through to air pressure - three
    slots. With one cell block per instrument it spans instruments 2..8
    entirely, so their cell *temperatures* started being filtered on a rule
    instrument 1's `tc` has never been subject to. On a two-analyser site that
    wiped the second analyser's cell temperature outright, and the physics then
    silently fell back to instrument 1's conditions for every gas.

The predicate names quantities instead of a slot span, which is what the
original list described back when there was only one instrument's cell block.
"""

from pathlib import Path
import re
import unittest


ROOT = Path(__file__).resolve().parents[1]
SRC = "src/src_rp/filter_dataset_for_diagnostics.f90"


def read(path):
    return (ROOT / path).read_text(encoding="utf-8")


class DiagnosticFilterScope(unittest.TestCase):
    def test_the_slot_span_is_gone(self):
        source = read(SRC)
        self.assertNotIn(
            "case (co2:gas4, pi:pe)", source,
            "the diagnostic filter is back to a slot span: it misses gases 5+ "
            "and swallows instruments 2..8's cell temperatures")

    def test_the_loop_covers_every_gas(self):
        source = read(SRC)
        self.assertIn("do var = firstGas, pe", source)

    def test_the_predicate_names_quantities(self):
        source = read(SRC)
        self.assertIn("logical function DiagFilterable(slot)", source)
        # Gas slots.
        self.assertIn("slot >= firstGas .and. slot <= lastGas", source)
        # One cell pressure per instrument - offset 3, where `pi` sits.
        self.assertIn("mod(slot - firstCell, NumCellPerInstr) == 3", source)
        # Ambient T/P.
        self.assertIn("slot == te .or. slot == pe", source)

    def test_cell_temperatures_are_not_filterable(self):
        """The asymmetry that caused the bug.

        Offsets 0, 1 and 2 of an instrument's block are cell_t, int_t_1 and
        int_t_2. Instrument 1's are tc/ti1/ti2, which sit below `pi` and were
        never in the historical list; making instrument 2's filterable while
        instrument 1's are not is what wiped the second analyser.
        """
        source = read(SRC)
        for offset in (0, 1, 2):
            self.assertNotIn(
                f"mod(slot - firstCell, NumCellPerInstr) == {offset}", source,
                f"cell block offset {offset} is a temperature; filtering it "
                f"would treat a second analyser differently from the first")

    def test_the_missing_gas_report_covers_every_gas(self):
        source = read(SRC)
        self.assertIn("E2Col(firstGas:lastGas)%present", source)
        self.assertNotIn("E2Col(co2:gas4)%present", source)

    def test_the_loop_index_is_declared(self):
        """`i` was only ever used inside the unrolled arms.

        Widening the loop moved the arms under a new guard; without the
        declaration the file relies on implicit typing, which `implicit none`
        would reject and which no test would otherwise notice.
        """
        source = read(SRC)
        decls = re.findall(r"^\s*integer :: (.+)$", source, re.M)
        names = {n.strip() for line in decls for n in line.split(",")}
        self.assertIn("i", names)
        self.assertIn("var", names)


if __name__ == "__main__":
    unittest.main()
