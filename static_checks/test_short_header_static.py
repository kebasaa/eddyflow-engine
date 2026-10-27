"""A file shorter than its stated header is skipped, never fatal.

A Jonkershoek run died with "Fortran runtime error: End of file" on its last
archive, written as the logger stopped: two minutes of raw data, and an
embedded biomet table holding its header and no record. The biomet reader
skipped the header with a bare `read(udf, *)`, so the read past the end of that
six-line table - whose metadata said seven - was a runtime error, and a
runtime error in one pre-pass worker ends the whole run. The raw ASCII reader
skipped its header the same way. The reader for files with text columns did
not crash but hung: gfortran answers every read after the end with a positive
code, and that loop took it for a bad line and read again.

The same archives declared one header row more than they had in every file.
Their biomet tables hold 30 one-minute rows, of which only the first carries
values and the other 29 are -9999, so the header skip took the one row that
mattered: every period reported 29 records imported and every biomet output
was -9999, with nothing naming why. The end of the run now says so once, when
no embedded biomet table gave a single value.

Part of the EddyFlow engine's static checks.
"""

from pathlib import Path
import re
import unittest


ROOT = Path(__file__).resolve().parents[1]


def read(path):
    return (ROOT / path).read_text(encoding="utf-8")


def code(path):
    return (chr(10)).join(ln for ln in read(path).splitlines()
                          if not ln.lstrip().startswith("!"))


def between(text, start, end):
    i = text.index(start)
    return text[i:text.index(end, i)]


def header_skip(path, routine):
    """The header-skip loop of routine, from its header_rows test onwards."""
    body = between(code(path), "subroutine " + routine,
                   "end subroutine " + routine)
    m = re.search(r"if \((bFileMetadata%nhead|FileInterpreter%header_rows) > 0\)",
                  body)
    return between(body[m.start():], "do i = 1", "end do")


class AShortFileIsSkipped(unittest.TestCase):

    def test_the_embedded_biomet_header_skip_cannot_run_off_the_end(self):
        loop = header_skip("src/src_rp/read_biomet_file.f90", "ReadBiometFile")
        self.assertIn("iostat = io_status", loop)
        self.assertIn("skip_file = .true.", loop)
        self.assertIn("return", loop)

    def test_the_raw_ascii_header_skip_cannot_run_off_the_end(self):
        loop = header_skip("src/src_common/import_ascii.f90", "ImportAscii")
        self.assertIn("iostat = io_status", loop)
        self.assertIn("header_short = .true.", loop)

    def test_a_short_raw_file_ends_as_one_read_to_its_end(self):
        #> FileEndReached sends the period on to its next file; without it a
        #> period with a LastRecord short of a full file stops at this one.
        body = code("src/src_common/import_ascii.f90")
        self.assertRegex(body, r"if \(header_short\) then\s+FileEndReached = \.true\.")
        #> And it still reaches the compaction and the close below.
        tail = body[body.index("if (header_short) then"):]
        self.assertIn("LocCol = TmpCol", tail)
        self.assertIn("close(unat)", tail)

    def test_no_header_skip_is_left_without_an_iostat(self):
        for path in ("src/src_rp/read_biomet_file.f90",
                     "src/src_common/import_ascii.f90"):
            bare = re.findall(r"^\s*read\s*\(\s*u\w+\s*,\s*\*\s*\)\s*$",
                              code(path), re.M)
            self.assertFalse(bare, "%s still has %r" % (path, bare))

    def test_the_text_reader_stops_on_a_read_after_the_end(self):
        body = code("src/src_common/import_ascii_with_text.f90")
        loop = between(body, "record_loop: do", "end do record_loop")
        first = loop[loop.index("read(unat"):loop.index("cycle record_loop")]
        self.assertIn("io_status == 5001", first)
        self.assertIn("io_status == 5008", first)
        self.assertLess(first.index("io_status == 5001"),
                        first.index("exit record_loop"))


class ARunWithoutAnyBiometValueSaysSo(unittest.TestCase):

    def test_each_embedded_table_is_noted_where_it_is_read(self):
        body = code("src/src_rp/read_licor_ghg_archive.f90")
        block = between(body, "call ReadBiometMetaFile", "if (MetaIsNeeded)")
        self.assertIn("EmbBiometTableSeen = .true.", block)
        #> A value, not a record: 29 rows of -9999 are records.
        seen = block[block.index("EmbBiometTableSeen = .true."):]
        self.assertIn("fbSet(1:fnbRecs, :) /= error", seen)
        self.assertIn("EmbBiometValueSeen = .true.", seen)

    def test_only_the_main_pass_counts(self):
        #> The pre-passes read the same archives; reset as the main pass starts.
        body = code("src/src_rp/eddyflow-rp_main.f90")
        head = body[:re.search(r"^\s*periods_loop: do", body, re.M).start()]
        head = head[head.rindex("call ResetRunColumns()"):]
        self.assertIn("EmbBiometTableSeen = .false.", head)
        self.assertIn("EmbBiometValueSeen = .false.", head)

    def test_the_warning_comes_after_the_workers_are_merged(self):
        body = code("src/src_rp/eddyflow-rp_main.f90")
        warn = body.index("if (EmbBiometTableSeen .and. .not. EmbBiometValueSeen)")
        self.assertLess(body.rindex("call MergeProdPieces", 0, warn), warn)
        self.assertLess(body.index("call ProgressEnd()"), warn)
        self.assertIn("biomet_header_rows", body[warn:warn + 800])

    def test_a_worker_hands_back_what_its_piece_met(self):
        body = code("src/src_rp/production_parallel.f90")
        finish = between(body, "subroutine FinishProdWorker",
                         "end subroutine FinishProdWorker")
        merge = between(body, "subroutine MergeProdPieces",
                        "end subroutine MergeProdPieces")
        #> Same place in the dump on both sides: right after the counts.
        w = finish.index("write(ud) EmbBiometTableSeen, EmbBiometValueSeen")
        self.assertLess(finish.index("write(ud) nOk, NumUserVar, nbVars"), w)
        self.assertLess(w, finish.index("write(ud) n" + chr(10)))
        r = merge.index("read(ud) wTableSeen, wValueSeen")
        self.assertLess(merge.index("read(ud) wOk, wUser, wBiomet"), r)
        self.assertLess(r, merge.index("read(ud) n" + chr(10)))
        self.assertIn("EmbBiometValueSeen = EmbBiometValueSeen .or. wValueSeen",
                      merge)


if __name__ == "__main__":
    unittest.main()
