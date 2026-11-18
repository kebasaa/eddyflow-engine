"""PackFlagString is never told an array is longer than it is.

PackFlagString(digits, n, string) reads digits(1:n), up to the length of the
flag string - 1 + GHGNumVar. The three spike routines declare their flags
hflags(u:lastGas) and passed GHGNumVar as the count, so it read past the end of
a slice of about ten values: leftover memory as flag digits after the last gas,
and since the 128-gas ceiling a segmentation fault in the spike test, partway
through a production run on CH-LAE.

Pinned: no call pairs a slice ending at lastGas or pe with a ceiling count, and
the three spike routines pass their array's own size.
"""

import re
import unittest
from pathlib import Path

SRC = Path(__file__).resolve().parent.parent / "src"
CALL = re.compile(r"call PackFlagString\(\s*([^,]+?)\s*,\s*([^,]+?)\s*,")


def calls():
    for path in sorted(SRC.rglob("*.f90")):
        text = path.read_text(encoding="utf-8", errors="replace")
        for m in CALL.finditer(text):
            yield path.name, m.group(1), m.group(2)


class NoCountOutrunsItsArray(unittest.TestCase):

    def test_there_are_calls_to_check(self):
        self.assertGreaterEqual(len(list(calls())), 10)

    def test_no_narrow_slice_with_a_ceiling_count(self):
        for name, arg, count in calls():
            with self.subTest(file=name, call=arg):
                narrow = re.search(r":\s*(lastGas|pe)\s*\)$", arg)
                ceiling = count in ("GHGNumVar", "E2NumVar")
                self.assertFalse(narrow and ceiling,
                                 "%s passes %s with count %s" % (name, arg, count))

    def test_the_spike_routines_pass_their_own_size(self):
        for f in ("test_spike_detection_vickers_97.f90",
                  "test_spike_detection_mauder_13.f90",
                  "despike_consecutive_diff.f90"):
            text = (SRC / "src_rp" / f).read_text(encoding="utf-8", errors="replace")
            self.assertIn("integer :: hflags(u:lastGas)", text, f)
            self.assertIn("call PackFlagString(hflags, size(hflags), CharHF%sr)", text, f)


if __name__ == "__main__":
    unittest.main()
