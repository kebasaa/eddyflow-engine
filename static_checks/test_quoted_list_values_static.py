"""A comma list written in quotes still parses.

The GUI writes a project value that holds a comma in quotes:
`cec_3_extra="5,4"`, `gas_1_sa_months="1-6,7-12"`. The INI reader strips quotes
only at the start of a whole line (`key=value` never starts with one), so the
value reaches the parsers with its quotes on.

ParseCecExtraList turned commas into blanks and read integers, so the first
read met `"5`, failed and dropped the whole list: the pairing ran with no extra
gases and said nothing. A single extra (`cec_2_extra=6`) is written without
quotes and worked, which is why it went unnoticed until a pairing named two.
ParseMonthGrouping is all-or-nothing and refused a quoted list outright.

Both now treat a quote as a blank before they read anything.
"""

from pathlib import Path
import re
import unittest


ROOT = Path(__file__).resolve().parents[1]
PROJECT = "src/src_common/write_processing_project_variables.f90"
MONTHS = "src/src_fcc/parse_month_grouping.f90"


def code(path):
    text = (ROOT / path).read_text(encoding="utf-8", errors="replace")
    return "\n".join(ln for ln in text.splitlines()
                     if not ln.lstrip().startswith("!"))


def routine(path, name):
    body = code(path)
    start = body.index("subroutine %s(" % name)
    return body[start:body.index("end subroutine %s" % name, start)]


def blanks_quotes(block, var):
    """The routine blanks both quote characters in `var`."""
    dq = re.search(r"%s\(\w+:\w+\) == '\"'" % var, block)
    sq = re.search(r"%s\(\w+:\w+\) == \"'\"" % var, block)
    return dq is not None and sq is not None


class TheCecExtraListIgnoresQuotes(unittest.TestCase):
    def setUp(self):
        self.block = routine(PROJECT, "ParseCecExtraList")

    def test_quotes_are_blanked(self):
        self.assertTrue(blanks_quotes(self.block, "buf"))

    def test_before_the_first_integer_read(self):
        blank = self.block.index("== '\"'")
        self.assertLess(blank, self.block.index("read(buf"))


class TheMonthGroupingIgnoresQuotes(unittest.TestCase):
    def setUp(self):
        self.block = routine(MONTHS, "ParseMonthGrouping")

    def test_quotes_are_blanked(self):
        self.assertTrue(blanks_quotes(self.block, "body"))

    def test_before_the_list_is_split(self):
        blank = self.block.index("== '\"'")
        self.assertLess(blank, self.block.index("index(body(pos:"))


if __name__ == "__main__":
    unittest.main()
