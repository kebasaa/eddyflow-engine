"""The PWB half-hourly time-lag table reads back what it wrote.

corr_pw and cv_99 went out as f12.8. Both are bounded by one, except where a
period had no detection of its own and kept the error code, which to eight
decimals is fourteen characters - so the table carried twelve asterisks there.
Its reader failed on the first such row and refused the whole table without a
word, and the caller stopped the run with "see the message above" over an
empty log: CH-LAE above canopy, step 2, against a step-1 table twelve hours in
the making.

Pinned: those two columns are written through PwbFixed8, which spells the
error code as every other column does; the reader turns a field of asterisks
back into the error code before it parses a row, so tables written before
this read too; and every refusal after the version check says why.
"""

import re
import unittest
from pathlib import Path

SRC = Path(__file__).resolve().parent.parent / "src" / "src_rp" / "pwb_timelag_handle.f90"
MOD = SRC.read_text(encoding="utf-8", errors="replace")


def between(start, end):
    a = MOD.index(start)
    return MOD[a:MOD.index(end, a)]


WRITER = between("subroutine WritePwbTimelagCache()", "end subroutine WritePwbTimelagCache")
READER = between("subroutine ReadPwbTimelagCache(", "end subroutine ReadPwbTimelagCache")
FIXED8 = between("character(12) function PwbFixed8(", "end function PwbFixed8")


class TheWriterFitsTheErrorCode(unittest.TestCase):

    def test_no_eight_decimal_field_is_left_in_the_row_format(self):
        self.assertNotIn("f12.8", WRITER)

    def test_both_columns_go_through_the_helper(self):
        self.assertIn("PwbFixed8(PwbTimelagCache(i)%result%corr_pw)", WRITER)
        self.assertIn("PwbFixed8(PwbTimelagCache(i)%result%cv_99)", WRITER)

    def test_the_error_code_is_written_as_the_other_columns_write_it(self):
        self.assertIn("if (x == error) then", FIXED8)
        self.assertIn("write(PwbFixed8, '(f12.6)') x", FIXED8)
        self.assertIn("write(PwbFixed8, '(f12.8)') x", FIXED8)


class TheReaderTakesOldTables(unittest.TestCase):

    def test_asterisks_are_cleared_before_the_row_is_parsed(self):
        self.assertLess(READER.index("call PwbUnstarFields(line)"),
                        READER.index("read(line, *, iostat=ios) date, time"))

    def test_they_become_the_error_code(self):
        fn = between("subroutine PwbUnstarFields(", "end subroutine PwbUnstarFields")
        self.assertIn("errText = '-9999.0'", fn)
        self.assertIn("verify(trim(adjustl(line(start:k))), '*') == 0", fn)


class EveryRefusalSaysWhy(unittest.TestCase):

    def test_no_silent_return_after_the_version_check(self):
        body = READER[READER.index("recognized = .true."):READER.index("\ncontains")]
        lines = body.splitlines()
        returns = 0
        for i, text in enumerate(lines):
            if text.strip() != "return":
                continue
            returns += 1
            #> Only the if-block this return closes: a message from the
            #> check before it must not count for this one.
            k = i - 1
            while k >= 0 and not lines[k].strip().startswith("if ("):
                k -= 1
            block = "\n".join(lines[k:i])
            self.assertTrue(re.search(r"call (LogSay|SayBadHeader)\(", block),
                            "a refusal without a message: %s" % lines[k].strip())
        self.assertGreaterEqual(returns, 7)

    def test_an_empty_table_says_so(self):
        tail = READER[READER.index("valid = PwbTimelagCacheN > 0"):]
        self.assertIn("if (.not. valid) then", tail)
        self.assertIn("has a header but no rows", tail)

    def test_an_unreadable_row_is_named_by_its_line(self):
        self.assertIn("lineno = lineno + 1", READER)
        self.assertIn("' could not be read:'", READER)


if __name__ == "__main__":
    unittest.main()
