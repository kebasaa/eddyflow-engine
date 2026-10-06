"""FCC's serial hot spots, each removed without changing a byte of output.

Measured on the Yatir day repeated over 30 days (1440 periods): the run took
196 s, the flux loop 2.7 s per day. Most of it was replace2. ReadExRecord
turns the error label into -9999 in every essentials row it reads, FCC reads
the file about four times, and replace2 rebuilt a buffer ten times the row -
640 KB - for every match: some 150 per row with the default label, where every
match is replaced by itself. After: 36 s, 0.8 s per day.
"""

import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def read(rel):
    return (ROOT / rel).read_text(encoding="utf-8", errors="replace").replace("\r\n", "\n")


def body(text, kind, name):
    i = text.index("%s %s(" % (kind, name))
    return text[i:text.index("end %s %s" % (kind, name), i)]


class Replace2ReturnsTheInputWhenNothingChanges(unittest.TestCase):

    B = body(read("src/src_common/string_sub.f90"), "function", "replace2")

    def test_a_match_replaced_by_itself_returns_at_once(self):
        i = self.B.index("if (len(what) == len(with)) then")
        self.assertLess(i, self.B.index("tstring = string"))
        self.assertIn("if (what == with) then", self.B[i:i + 120])

    def test_no_match_returns_at_once(self):
        self.assertLess(self.B.index("if (index(string, what) == 0) then"),
                        self.B.index("tstring = string"))

    def test_both_return_the_trimmed_input_as_the_loop_would(self):
        self.assertEqual(self.B.count("nstring = trim(string)"), 2)


class TheFluxDespikingCacheGrowsByDoubling(unittest.TestCase):

    def test_it_is_copied_only_when_full(self):
        b = body(read("src/src_fcc/pfd_handle.f90"), "subroutine", "StorePfdCache")
        self.assertIn("if (PfdCacheN >= size(PfdCache)) then", b)
        self.assertIn("allocate(tmp(2 * size(PfdCache)))", b)
        self.assertNotIn("allocate(tmp(PfdCacheN + 1))", b)


class TheCospectraFitDatasetGrowsAsItFills(unittest.TestCase):

    MAIN = read("src/src_fcc/eddyflow-fcc_main.f90")

    def test_it_is_not_sized_for_every_period_of_the_range(self):
        self.assertNotIn("(saEndTimestampIndx - saStartTimestampIndx + 1)))", self.MAIN)

    def test_it_grows_before_a_file_s_bins_are_added(self):
        i = self.MAIN.index("if (maxval(nfit) + nbins > size(FitStable)) then")
        self.assertLess(i, self.MAIN.index("call AddToCospectraFitDataset("))
        grow = self.MAIN[i:self.MAIN.index("end if", i)]
        self.assertEqual(grow.count("call move_alloc(FitGrown,"), 2)
        self.assertEqual(grow.count("FitGrown = NullFitCosp"), 2)


if __name__ == "__main__":
    unittest.main()
