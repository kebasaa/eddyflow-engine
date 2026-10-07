"""Fratini's full-cospectra file is found by a binary search, as the scan found it.

Every period used to walk the whole list of full-cospectra files, converting
the period's own date again at each step - quadratic in the length of the
run. The list is fixed for a run, so it is sorted once and each period
searches it. The answer has to be exactly the scan's: the lowest list index
whose timestamp equals the period's, by the same equality.

Part of the EddyFlow engine's static checks.
"""

from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[1]


def code(rel):
    text = (ROOT / rel).read_text(encoding="utf-8", errors="replace")
    return "\n".join(ln for ln in text.splitlines() if not ln.lstrip().startswith("!"))


FRATINI = code("src/src_common/bpcf_fratini_12.f90")
LIBDATE = code("src/src_common/m_libdate.f90")


def body(kind, name):
    i = FRATINI.index(kind + " " + name + "(")
    return FRATINI[i:FRATINI.index("end " + kind + " " + name, i)]


class TheFileIsFoundOnce(unittest.TestCase):

    def test_the_period_date_is_converted_once(self):
        top = FRATINI[:FRATINI.index("contains")]
        self.assertEqual(top.count("call DateTimeToDateType(lEx%end_date, lEx%end_time, Timestamp)"), 1)
        self.assertIn("indx = FirstFileAt(Timestamp)", top)
        self.assertNotIn("do i = 1, nfull", top)

    def test_the_index_is_kept_for_the_run(self):
        for decl in ("integer, allocatable, save :: Order(:)",
                     "type(DateType), allocatable, save :: SortedTs(:)",
                     "integer, save :: nIndexed = -1"):
            self.assertIn(decl, FRATINI)
        self.assertIn("if (.not. IndexCurrent()) call BuildIndex()", body("function", "FirstFileAt"))


class ItFindsWhatTheScanFound(unittest.TestCase):

    def test_order_compares_exactly_the_fields_equality_does(self):
        #> Equality in m_libdate is field by field; ordering there goes through
        #> Julian dates. The index orders by the same five fields equality
        #> reads, so "neither before the other" means "equal" here too.
        before = body("function", "Before")
        for field in ("Year", "Month", "Day", "Hour", "Minute"):
            self.assertIn("Date1%" + field, LIBDATE)
            self.assertIn("Before = d1%" + field + " < d2%" + field, before)

    def test_ties_go_to_the_lowest_index(self):
        self.assertIn("Precedes = i < j", body("function", "Precedes"))

    def test_the_search_lands_on_the_first_match(self):
        search = body("function", "FirstFileAt")
        self.assertIn("if (Before(SortedTs(mid), ts)) then", search)
        self.assertIn("if (SortedTs(lo) == ts) FirstFileAt = Order(lo)", search)
        self.assertIn("FirstFileAt = nint(error)", search)


if __name__ == "__main__":
    unittest.main()
