"""A period skipped before the FLUXNET header exists still gets its row.

RP writes the FLUXNET header lazily, after the first period that imports
data, because the custom-variable columns are only known then - and with
them the width of every row, a skipped period's included. Periods skipped
before that point were written anyway, by WriteOutFluxnetOnlyBiomet, to
``uflxnt`` while nothing had opened it. gfortran connected the unit to
``fort.132`` in the working directory, the rows went there, and the FLUXNET
file began at the first period with data. Every period ahead of the first
raw file is such a period, so any window that opens before the data lost its
leading rows. They were malformed as well, padded to an ``nFluxnetFixedCols``
of zero.

The fix holds those periods back and has InitFluxnetFile_rp write their rows
straight after the header. What this pins:

* the skipped-period entry point never writes the unit itself - it writes
  only once the file is open, and defers otherwise;
* the row builder takes everything particular to its period as arguments.
  Read from ``Stats`` or the biomet aggregates instead, a deferred row would be
  built from whichever period happened to be current when it was flushed;
* InitFluxnetFile_rp flushes after the header, and nothing else claims the
  file is open;
* the data-row writer is still called only after the header exists.

The end-to-end gate is tests/regression/check_lead_gap.py on base_lead_gap.
"""

import re
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RP = ROOT / "src" / "src_rp"


def read(name):
    return (RP / name).read_text(encoding="utf-8", errors="replace")


def routine(text, name):
    """The body of subroutine `name`, up to its end statement."""
    start = re.search(r"^\s*subroutine\s+%s\s*\(" % name, text,
                      re.MULTILINE | re.IGNORECASE)
    assert start, "no subroutine %s" % name
    end = re.search(r"^\s*end\s+subroutine\s+%s\b" % name, text[start.end():],
                    re.MULTILINE | re.IGNORECASE)
    assert end, "no end of %s" % name
    return text[start.start():start.end() + end.end()]


def code(text):
    """Drop comments, so prose about a name does not count as using it."""
    return "\n".join(ln.split("!", 1)[0] for ln in text.splitlines())


WRITER = read("write_out_fluxnet_only_biomet.f90")
INIT = read("init_fluxnet_file_rp.f90")
MAIN = read("eddyflow-rp_main.f90")
UNIT_WRITE = re.compile(r"write\s*\(\s*uflxnt\b", re.IGNORECASE)


class TheSkippedPeriodEntryPointDefers(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.BODY = code(routine(WRITER, "WriteOutFluxnetOnlyBiomet"))

    def test_it_does_not_write_the_unit_itself(self):
        self.assertIsNone(
            UNIT_WRITE.search(self.BODY),
            "a write here runs before the header on a leading skipped period "
            "and lands in fort.132")

    def test_it_writes_only_once_the_file_is_open(self):
        gate = self.BODY.find("if (FluxnetFileOpen)")
        call = self.BODY.find("call WriteFluxnetOnlyBiometRow(")
        self.assertGreaterEqual(gate, 0)
        self.assertGreater(call, gate)

    def test_otherwise_it_keeps_the_period(self):
        self.assertIn("DeferredFluxnetRows(", self.BODY)
        self.assertIn("nDeferredFluxnetRows =", self.BODY)


class TheRowBuilderTakesItsPeriodAsArguments(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.BODY = code(routine(WRITER, "WriteFluxnetOnlyBiometRow"))

    def test_no_period_state_is_read_from_globals(self):
        for name in ("Stats%", "bAggr"):
            self.assertNotIn(
                name, self.BODY,
                "%s is the period current at flush time, not the one the row "
                "is for" % name)

    def test_it_is_the_only_skipped_row_writer(self):
        self.assertEqual(len(UNIT_WRITE.findall(code(WRITER))), 1)
        self.assertIsNotNone(UNIT_WRITE.search(self.BODY))


class TheHeaderFlushesWhatWasHeldBack(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.BODY = code(INIT[:INIT.index("\ncontains")])

    def test_flush_follows_the_header_write(self):
        header = [m.start() for m in UNIT_WRITE.finditer(self.BODY)]
        self.assertEqual(len(header), 1, "one header line")
        opened = self.BODY.find("FluxnetFileOpen = open_status == 0")
        flush = self.BODY.find("call FlushDeferredFluxnetRows()")
        self.assertGreater(opened, header[0])
        self.assertGreater(flush, opened)

    def test_only_the_header_routine_sets_the_flag(self):
        setters = []
        for path in sorted((ROOT / "src").rglob("*.f90")):
            text = code(path.read_text(encoding="utf-8", errors="replace"))
            if re.search(r"^\s*FluxnetFileOpen\s*=", text, re.MULTILINE):
                setters.append(path.name)
        self.assertEqual(["init_fluxnet_file_rp.f90"], setters)


class DataRowsStillFollowTheHeader(unittest.TestCase):

    def test_write_out_fluxnet_is_called_after_init(self):
        body = code(MAIN)
        init = body.find("call InitFluxnetFile_rp()")
        rows = [m.start() for m in
                re.finditer(r"call WriteOutFluxnet\(", body)]
        self.assertGreaterEqual(init, 0)
        self.assertTrue(rows)
        self.assertTrue(all(r > init for r in rows),
                        "a data row ahead of InitFluxnetFile_rp has no file")

    def test_the_unit_has_no_other_rp_writer(self):
        writers = sorted(p.name for p in RP.glob("*.f90")
                         if UNIT_WRITE.search(code(p.read_text(
                             encoding="utf-8", errors="replace"))))
        self.assertEqual(
            ["init_fluxnet_file_rp.f90", "write_out_fluxnet.f90",
             "write_out_fluxnet_only_biomet.f90"], writers)


if __name__ == "__main__":
    unittest.main()
