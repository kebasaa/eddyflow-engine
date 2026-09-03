"""Decompressing the next GHG archive must not change which one is read.

A LI-COR archive costs about 170 ms more than the same data as a plain file,
and roughly 100 ms of that is 7-Zip running over the archive - time this
program spends waiting on another process, so it can be spent computing a flux
instead. Measured over `base_ghg` against `base_tlag_opt` on identical data,
prefetching recovers about 2.8 % of a run.

That is a small return, and it is only acceptable because the failure mode is
nothing at all. A prefetch is claimed **only** when the archive path matches
exactly and the extraction has said it finished; every other case - no request,
a request for a different file, one still running, one that failed, a directory
left over from a killed run - falls through to extracting synchronously exactly
as before. So the checks here are almost all about the ways a claim must NOT be
granted.

Two details carry the whole thing.

**The sentinel is a separate command after the extraction.** Testing for the
extracted files instead would race: they appear while 7-Zip is still writing
them, and the reader would get half a file.

**The request comes after the current archive's files are deleted.** Otherwise
the next extraction writes into a directory the caller is still reading out of.

There is also a Fortran constraint worth stating, because breaking it compiles
on one call and not the other: `UnZipArchive` and `ReadLicorGhgArchive` are
external procedures with no explicit interface, so their new arguments cannot
be `optional`. gfortran infers an interface from whichever call it sees first
and then rejects the other.
"""

import re
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / "src"


def read(rel):
    return (SRC / rel).read_text(encoding="utf-8", errors="replace")


PREFETCH = read("src_common/ghg_prefetch.f90")
UNZIP = read("src_common/unzip_archive.f90")
GHG = read("src_rp/read_licor_ghg_archive.f90")
IMPORT = read("src_rp/import_current_period.f90")
MAIN = read("src_rp/eddyflow-rp_main.f90")


def body(source, signature, end):
    start = source.index(signature)
    return source[start:source.index(end, start)]


class AClaimIsGrantedOnlyForTheRightFinishedArchive(unittest.TestCase):

    CLAIM = property(lambda self: body(
        PREFETCH, "subroutine GhgPrefetchClaim",
        "end subroutine GhgPrefetchClaim"))

    def test_nothing_pending_is_not_a_claim(self):
        self.assertIn("if (len_trim(Pending) == 0) return", self.CLAIM)

    def test_a_different_archive_is_not_a_claim(self):
        self.assertIn("if (trim(Pending) /= trim(ZipFile)) return", self.CLAIM)

    def test_an_unfinished_extraction_is_not_a_claim(self):
        self.assertIn("'ready.flag', exist = ex)", self.CLAIM)
        self.assertIn("if (.not. ex) return", self.CLAIM)

    def test_the_default_is_no(self):
        """Every path out before the checks must leave ready false."""
        i = self.CLAIM.index("ready = .false.")
        self.assertLess(i, self.CLAIM.index("if (len_trim(Pending) == 0) return"))

    def test_a_claim_is_consumed(self):
        """Or the same directory would be claimed twice, once after it has
        been overwritten by the following prefetch."""
        self.assertIn("Pending = ''", self.CLAIM)


class TheSentinelSaysTheExtractionFinished(unittest.TestCase):

    START = property(lambda self: body(
        PREFETCH, "subroutine GhgPrefetchStart",
        "end subroutine GhgPrefetchStart"))

    def test_it_is_written_after_the_extraction_not_by_it(self):
        """A separate command, so it cannot run until 7-Zip returned."""
        script = self.START
        seven = script.index("comm_7zip")
        flag = script.index("ready.flag")
        self.assertLess(seven, flag,
                        "the flag must be written after the extraction")

    def test_the_directory_is_cleared_before_extracting(self):
        """A killed run leaves files and possibly a flag behind."""
        self.assertRegex(self.START, r"rmdir|rm -rf")


class NothingIsExtractedIntoADirectoryBeingRead(unittest.TestCase):

    def test_the_request_follows_the_release(self):
        """Straight after, so nothing can slip between them and find the
        directory still marked in use."""
        start = GHG.index("call GhgPrefetchStart(")
        before = GHG[:start].rstrip()
        self.assertTrue(before.endswith("call GhgPrefetchRelease()"),
                        "the request must come straight after the release")


    def test_the_request_follows_the_delete(self):
        """The files have to be gone before the directory is reused."""
        delete = GHG.index("del_status = system(trim(comm))")
        self.assertLess(delete, GHG.index("call GhgPrefetchStart("))

    def test_a_request_while_in_use_is_refused(self):
        start = body(PREFETCH, "subroutine GhgPrefetchStart",
                     "end subroutine GhgPrefetchStart")
        self.assertIn("if (InUse) return", start)

    def test_only_one_request_is_outstanding(self):
        start = body(PREFETCH, "subroutine GhgPrefetchStart",
                     "end subroutine GhgPrefetchStart")
        self.assertIn("if (len_trim(Pending) /= 0) return", start)

    def test_every_early_return_releases(self):
        """A return that skips the release leaves InUse set, and every later
        prefetch in the run is then refused. Each early return must be
        immediately preceded by one - the final release is the normal exit and
        is paired with the request instead, so it is excluded."""
        section = GHG[GHG.index("call GhgPrefetchClaim("):
                      GHG.rindex("call GhgPrefetchRelease()")]
        rows = section.split(chr(10))
        unguarded = [n for n, line in enumerate(rows)
                     if line.strip() == "return"
                     and rows[n - 1].strip() != "call GhgPrefetchRelease()"]
        self.assertEqual([], unguarded,
                         "a return that does not release first, at line(s) %s "
                         "of the section" % unguarded)


class ItRespectsTheSerialSwitch(unittest.TestCase):

    def test_j_one_starts_nothing(self):
        start = body(PREFETCH, "subroutine GhgPrefetchStart",
                     "end subroutine GhgPrefetchStart")
        self.assertIn("if (NumJobs == 1) return", start)


class TheNewArgumentsAreRequired(unittest.TestCase):
    """Both are external procedures with no explicit interface, so `optional`
    is not something the standard defines for them."""

    def test_unzip_takes_them_positionally(self):
        self.assertIn(
            "BiometFile, BiometMetaFile, skip_file, WorkDir, prefetched)", UNZIP)
        decl = UNZIP[:UNZIP.index("call clearstr(MetaFile)")]
        #> The attribute, not the word - the comment beside them explains
        #> why they are not optional.
        self.assertNotIn(", optional", decl)

    def test_the_reader_takes_its_next_archive_positionally(self):
        decl = GHG[:GHG.index("skip_file = .false.")]
        self.assertNotIn(", optional", decl)

    def test_the_preamble_says_it_has_no_next(self):
        """It reads one file to learn the columns and stops."""
        i = MAIN.index("call ReadLicorGhgArchive(RawFileList(i)%path")
        self.assertIn("'')", MAIN[i:i + 500])

    def test_the_period_loop_names_the_next_file(self):
        i = IMPORT.index("call ReadLicorGhgArchive(")
        before = IMPORT[max(0, i - 500):i]
        self.assertIn("if (CurrentFile < NumFiles) then", before)
        self.assertIn("NextZip = FileList(CurrentFile + 1)%path", before)


class OneListingNotTwo(unittest.TestCase):
    """Both were only ever asking what the archive held, and each cost a shell."""

    def test_a_single_listing(self):
        self.assertEqual(UNZIP.count("dir_status = system(comm)"), 1)
        self.assertNotIn("meta_flist.tmp", UNZIP)
        self.assertNotIn("data_flist.tmp", UNZIP)

    def test_it_ignores_its_own_output_and_the_flag(self):
        self.assertIn("if (index(dataline, '.tmp') /= 0) cycle", UNZIP)
        self.assertIn("if (index(dataline, 'ready.flag') /= 0) cycle", UNZIP)

    def test_biomet_is_matched_before_the_plain_extension(self):
        """`-biomet.data` also ends in `.data`, so order decides."""
        block = UNZIP[UNZIP.index("arch_flist.tmp'"):]
        biomet = block.index("'-biomet.' // trim(adjustl(DataExt))")
        plain = block.index("'.' // trim(adjustl(DataExt))")
        self.assertLess(biomet, plain)


class TheDirectoryIsCleanedUp(unittest.TestCase):

    def test_the_run_removes_it(self):
        self.assertIn("call GhgPrefetchCleanup()", MAIN)

    def test_before_the_temporary_directory_goes(self):
        """Desktop mode deletes the whole tmp tree; embedded mode keeps it."""
        self.assertLess(MAIN.index("call GhgPrefetchCleanup()"),
                        MAIN.index("Delete tmp folder if running in embedded mode"))


if __name__ == "__main__":
    unittest.main()
