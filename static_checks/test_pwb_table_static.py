"""The PWB table stores and finds rows without a scan from row 1 - same rows.

Three quadratics, all in the per-period PWB table and all growing with the
square of the record length:

* every row stored reallocated the table to exactly one row more and copied
  the lot (loading a cache file went the same way);
* every store first scanned the table from row 1 for an existing key;
* every lookup in the production pass scanned from row 1 too.

Measured with the check program (a year of half-hours, four gases, 70 080
rows): storing took 172 s before and 1.5 s after; looking every key up took
3.0 s before and 0.8 s after.

Now the table has spare capacity that doubles when full (PwbTimelagCacheN is
still the number of rows), and while its rows are in nondecreasing period
order - which the engine always builds - a period's rows are found by
bisection on a minutes column kept beside the table. Rows with the same date
and time have the same minutes, so in an ordered table they all lie in one
block, and the first match in the block is the first match in the table: the
same row a scan from row 1 returns, duplicate keys included. An out-of-order
table is scanned as before.

That rests on every append going through one place that keeps the order flag
and the minutes column true, and on no row's key ever changing - which these
checks pin - and is proven by ``pwb_table_check`` (``mingw32-make
pwbtablecheck``), which compares every lookup with a scan from row 1.

Part of the EddyFlow engine's static checks.
"""

from pathlib import Path
import re
import subprocess
import unittest

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
DRIVER = ROOT / "obj" / "win" / "pwb_table_check.exe"


def code(path):
    text = path.read_text(encoding="utf-8", errors="replace")
    return "\n".join(ln for ln in text.splitlines() if not ln.lstrip().startswith("!"))


PWB = code(SRC / "src_rp" / "pwb_timelag_handle.f90")
ALL = {p.name: code(p) for p in SRC.rglob("*.f90") if "src_tools" not in p.parts}


def statements(text):
    """Code lines that are not declarations, with the routine each sits in."""
    where = None
    for ln in text.splitlines():
        t = ln.strip()
        low = t.lower()
        for kw in ("subroutine ", "function "):
            i = low.find(kw)
            if i >= 0 and not low.startswith("end") and "call " not in low:
                where = t[i + len(kw):].split("(")[0].strip()
        if "::" in t:
            continue
        yield where, t


def body(name, kind="subroutine"):
    return PWB[PWB.index("%s %s(" % (kind, name)):PWB.index("end %s %s" % (kind, name))]


class RowsAreAddedInOnePlace(unittest.TestCase):

    def test_only_two_routines_grow_the_row_count(self):
        growers = set()
        for name, text in ALL.items():
            for where, t in statements(text):
                if re.match(r"PwbTimelagCacheN\s*=(?!=)", t):
                    growers.add(where)
        self.assertEqual(growers, {"InitPwbTimelagCache", "AppendPwbCacheRows",
                                   "StorePwbTimelagCacheAt"})

    def test_both_appenders_keep_the_flag_and_minutes(self):
        for name in ("AppendPwbCacheRows", "StorePwbTimelagCacheAt"):
            b = body(name)
            self.assertIn("call EnsurePwbCacheCapacity(", b)
            self.assertIn("call NotePwbRowAppended(", b)

    def test_the_note_records_minutes_and_order(self):
        b = body("NotePwbRowAppended")
        self.assertIn("PwbCacheMinutes(i) = m", b)
        self.assertIn("if (i > 1 .and. m < PwbCacheLastMinutes) PwbCacheInTimeOrder = .false.", b)

    def test_capacity_doubles_and_moves_the_minutes_with_the_rows(self):
        b = body("EnsurePwbCacheCapacity")
        self.assertIn("capacity = max(n, 2 * size(PwbTimelagCache))", b)
        self.assertIn("call move_alloc(grown, PwbTimelagCache)", b)
        self.assertIn("call move_alloc(grown_min, PwbCacheMinutes)", b)

    def test_the_merge_uses_the_appender(self):
        self.assertIn("call AppendPwbCacheRows(rows, nrec)", ALL["prepass_parallel.f90"])

    def test_init_resets_everything(self):
        b = body("InitPwbTimelagCache")
        for line in ("deallocate(PwbTimelagCache)", "deallocate(PwbCacheMinutes)",
                     "PwbCacheInTimeOrder = .true."):
            self.assertIn(line, b)


class NoKeyChangesAndNothingReadsTheSpareCapacity(unittest.TestCase):

    def test_no_row_key_is_assigned_outside_the_appender(self):
        for name, text in ALL.items():
            for where, t in statements(text):
                if re.match(r"PwbTimelagCache\([^)]*\)%(date|time|gas)\s*=(?!=)", t):
                    self.assertEqual(where, "StorePwbTimelagCacheAt", "%s in %s" % (t, name))

    def test_the_table_is_never_read_whole(self):
        """It is larger than its contents now; size() is for the capacity
        routine alone, and nothing may use the array without a 1:N range."""
        for name, text in ALL.items():
            for where, t in statements(text):
                if "size(PwbTimelagCache)" in t:
                    self.assertEqual(where, "EnsurePwbCacheCapacity", name)
                self.assertNotIn("PwbTimelagCache(:)", t, name)
                self.assertNotIn("PwbTimelagCache%", t, name)


class LookupIsTheFirstMatch(unittest.TestCase):

    def test_store_and_lookup_both_locate(self):
        self.assertIn("i = LocatePwbRow(PwbPeriodDate, PwbPeriodTime, gas)", body("StorePwbTimelagCache"))
        self.assertIn("i = LocatePwbRow(PwbPeriodDate, PwbPeriodTime, gas)", body("LookupPwbTimelagCache"))

    def test_bisection_only_while_ordered_otherwise_the_old_scan(self):
        b = body("LocatePwbRow", "function")
        self.assertIn("if (.not. PwbCacheInTimeOrder) then", b)
        self.assertLess(b.index("if (.not. PwbCacheInTimeOrder) then"), b.index("do while (lo < hi)"))

    def test_the_block_is_scanned_in_row_order_from_its_start(self):
        b = body("LocatePwbRow", "function")
        self.assertIn("if (PwbCacheMinutes(mid) < key) then", b)
        self.assertIn("do i = lo, PwbTimelagCacheN", b)
        self.assertIn("if (PwbCacheMinutes(i) /= key) exit", b)

    @unittest.skipUnless(DRIVER.exists(), "pwb_table_check not built (mingw32-make pwbtablecheck)")
    def test_the_shipped_routines_return_the_scanned_row(self):
        r = subprocess.run([str(DRIVER)], capture_output=True, text=True, timeout=600)
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertIn("SAME ROW AS A SCAN FROM ROW 1", r.stdout)
        self.assertIn("wrong rows: 0", r.stdout)


if __name__ == "__main__":
    unittest.main()
