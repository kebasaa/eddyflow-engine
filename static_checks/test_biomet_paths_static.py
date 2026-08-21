"""A reader must never create the file it came to read.

The biomet path is the one user-named path the engine used to lose quietly.
`scanCsvFile` opened it without `status='old'`, so a path that did not exist was
CREATED, empty; the open then succeeded, the scan found no rows, and its own
"too few rows or columns" guard fired. The run reported a *content* problem for
a file that was simply not there, carried on with no biomet data, and left an
empty file behind in the working directory.

That is how a corrupted biomet path cost a real run its biomet data with
nothing in the log naming the actual fault. The checks here pin the two halves
of the fix: nothing in the biomet family creates its input, and a named file
that is absent stops the run the way every other named path already does.

Part of the EddyFlow engine's static checks.
"""

from pathlib import Path
import re
import unittest


ROOT = Path(__file__).resolve().parents[1]

#: The files that open a biomet path. dir_sub.f90 is deliberately absent:
#: it also holds directory helpers that write a listing to a temporary file
#: and so must be allowed to create. Its one relevant routine, scanCsvFile,
#: is checked on its own below.
BIOMET_SOURCES = (
    "src/src_rp/init_external_biomet.f90",
    "src/src_rp/read_biomet_file.f90",
    "src/src_rp/biomet_retrieve_external_data.f90",
)


def read(path):
    return (ROOT / path).read_text(encoding="utf-8")


def code(path):
    return (chr(10)).join(ln for ln in read(path).splitlines()
                          if not ln.lstrip().startswith("!"))


def opens_in(path):
    """Every open(...) statement, continuation lines joined.

    The argument list runs to the end of the statement rather than to
    the first closing bracket: `bFileList(nfl)%path` closes one of its
    own, and stopping there hides the `status=` that follows it.
    """
    joined = re.sub(r"&\s*\n\s*", " ", code(path))
    return [m.group(1) for m in
            re.finditer(r"^\s*open\s*\((.*)$", joined, re.M | re.I)]


class AReaderNeverCreatesItsInput(unittest.TestCase):

    def test_the_csv_scanner_opens_an_existing_file_only(self):
        #> The one that did the damage. A scanner has nothing to write.
        block = code("src/src_common/dir_sub.f90")
        block = block[block.index("subroutine scanCsvFile"):
                      block.index("end subroutine scanCsvFile")]
        self.assertIn("status='old'", block.replace(" ", ""))

    def test_no_biomet_open_omits_a_status(self):
        for path in BIOMET_SOURCES:
            for args in opens_in(path):
                self.assertIn(
                    "status", args.lower(),
                    "%s opens without a status, so a missing file is created "
                    "instead of reported: open(%s)" % (path, args.strip()))

    def test_the_scanner_reports_missing_apart_from_unusable(self):
        #> Without this the caller cannot tell a wrong path from a bad file,
        #> and every message it writes has to hedge between the two.
        block = code("src/src_common/dir_sub.f90")
        block = block[block.index("subroutine scanCsvFile"):
                      block.index("end subroutine scanCsvFile")]
        self.assertIn("logical, intent(out), optional :: missing", block)
        self.assertIn("if (present(missing)) missing = .true.", block)


class ANamedButMissingBiometFileStopsTheRun(unittest.TestCase):

    def test_it_aborts_through_the_shared_helper(self):
        #> pf_file, to_file, sa_file and proj_file all stop here already.
        #> Biomet was the one that did not.
        source = code("src/src_rp/init_external_biomet.f90")
        self.assertIn("call AbortOnMissingPath('biom_file'", source)

    def test_only_a_file_the_project_named(self):
        #> A file that came from a directory listing existed a moment ago, so
        #> its disappearance is a different and much rarer event; taking the
        #> whole run down for it would be wrong.
        source = code("src/src_rp/init_external_biomet.f90")
        self.assertIn("EddyFlowProj%biomet_data == 'ext_file'", source)
        guard = source[source.index("if (missing"):
                       source.index("call AbortOnMissingPath('biom_file'")]
        self.assertIn("ext_file", guard)

    def test_a_short_file_still_reports_the_old_error(self):
        #> Error(2) describes a file that scanned badly, which is right when
        #> the file is actually there.
        source = code("src/src_rp/init_external_biomet.f90")
        self.assertIn("call ExceptionHandler(2)", source)


class TheFixturesSayWhatTheyDo(unittest.TestCase):

    def test_a_fixture_naming_a_biomet_file_can_reach_it(self):
        #> Thirty-eight fixtures declared an external biomet file on a drive
        #> that does not exist, and got Error(2) plus a run with no biomet.
        #> That was invisible while a missing file was tolerated and is a hard
        #> stop now, so a fixture either names a file it can open or says it
        #> uses none.
        import os
        offenders = []
        for p in sorted((ROOT / "tests/regression").glob("*.eddyflow")):
            txt = p.read_text(encoding="utf-8", errors="replace")
            mode = re.search(r"^use_biom=(\d)", txt, re.M)
            if mode is None or mode.group(1) != "2":
                continue
            named = re.search(r"^biom_file=(.*)$", txt, re.M)
            path = named.group(1).strip() if named else ""
            if not path or not os.path.exists(path):
                offenders.append("%s -> %r" % (p.name, path))
        self.assertFalse(
            offenders,
            "these fixtures request an external biomet file they cannot "
            "open, which now aborts the run:\n  " + "\n  ".join(offenders))


if __name__ == "__main__":
    unittest.main()
