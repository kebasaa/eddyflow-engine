"""Inputs may come from a shared Google Drive or Dropbox link; outputs may not.

src/src_common/remote_source.f90 does the downloading. What keeps it correct
is not in that file but in where it is called from, so that is what is held
here:

**Every raw file is fetched before it is opened.** A raw file from a link has
a local path from the moment the folder is listed, but no file behind it until
RemoteEnsure downloads it. The two readers every raw file goes through, and the
acquisition frequency survey, which runs 7-Zip over files directly, must each
ask first. A reader added later that opens RawFileList paths itself would read
nothing and skip every period without an error that says why.

**The listing goes where `dir` went.** FileListByExt and NumberOfFilesInDir
both list the raw folder; both must hand it to the provider, or the count and
the list disagree.

**The output folder is refused, and first.** EddyFlow writes there and a
shared link is read only. The refusal comes before any download, so a run that
cannot finish does not fetch a season of data first.

**Only the module runs curl.** One place sets the options that make a failed
download a failure (--fail, the HTML check) and the temporary-file rename that
keeps a half-written file from ever being read.
"""

import re
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"


def code(rel):
    text = (SRC / rel).read_text(encoding="utf-8", errors="replace")
    return "\n".join(ln for ln in text.splitlines() if not ln.strip().startswith("!"))


class RawFilesAreFetchedBeforeReading(unittest.TestCase):
    def test_ghg_reader_ensures_before_unzipping(self):
        c = code("src_rp/read_licor_ghg_archive.f90")
        self.assertIn("call RemoteEnsure(ZipFile, .true.)", c)
        self.assertLess(c.index("call RemoteEnsure(ZipFile"),
                        c.index("call GhgPrefetchClaim(ZipFile"))

    def test_native_reader_ensures_before_opening(self):
        c = code("src_common/import_native_data.f90")
        self.assertIn("call RemoteEnsure(Filepath, .true.)", c)
        self.assertLess(c.index("call RemoteEnsure(Filepath"), c.index("select case"))

    def test_survey_peek_ensures_without_fetching_ahead(self):
        c = code("src_rp/survey_ghg_ac_freq.f90")
        self.assertIn("call RemoteEnsure(GhgFile%path, .false.)", c)
        self.assertLess(c.index("call RemoteEnsure(GhgFile%path"),
                        c.index("unzip_status = system("))

    def test_ghg_prefetch_waits_for_the_download(self):
        c = code("src_common/ghg_prefetch.f90")
        self.assertIn("if (.not. RemoteIsLocal(ZipFile)) return", c)

    def test_processing_order_is_adopted_after_sorting(self):
        c = code("src_rp/eddyflow-rp_main.f90")
        self.assertLess(c.index("call FilesInChronologicalOrder(RawFileList"),
                        c.index("call RemoteAdoptOrder(RawFileList"))


class ListingGoesWhereDirWent(unittest.TestCase):
    def test_both_listers_defer_to_the_provider(self):
        for rel in ("src_common/filelist_by_ext.f90", "src_common/dir_sub.f90"):
            c = code(rel)
            with self.subTest(rel=rel):
                self.assertIn("if (RemoteOwnsDir(DirIn)) then", c)
                self.assertIn("call RemoteWriteFileList(Ext", c)


class OutputFolderIsRefusedFirst(unittest.TestCase):
    def test_rp(self):
        c = code("src_rp/read_ini_rp.f90")
        refuse = c.index("call RemoteRefuseOutput('out_path', Dir%main_out)")
        self.assertLess(refuse, c.index("call RemoteFetchFile("))
        self.assertLess(refuse, c.index("call RemoteResolveInputs("))

    def test_fcc(self):
        c = code("src_fcc/read_ini_fcc.f90")
        refuse = c.index("call RemoteRefuseOutput('out_path', Dir%main_out)")
        self.assertLess(refuse, c.index("call RemoteFetchFile("))
        self.assertLess(refuse, c.index("call RemoteFetchFolder("))

    def test_fcc_inputs_are_fetched_before_they_are_checked(self):
        # A link is not a directory, so the existence checks below would abort
        # on it with "missing path" if the download came after them.
        c = code("src_fcc/read_ini_fcc.f90")
        fetch = c.index("call RemoteFetchFolder('sa_full_spectra'")
        self.assertLess(fetch, c.index("inquire(file = Dir%full"))
        self.assertLess(fetch, c.index("inquire(file = Dir%binned"))

    def test_error_numbers_exist(self):
        c = code("src_common/exception_handler.f90")
        for n in (120, 121, 122):
            with self.subTest(n=n):
                self.assertIn(f"case({n})", c)


class OnlyTheModuleRunsCurl(unittest.TestCase):
    def test_no_curl_elsewhere(self):
        pattern = re.compile(r"['\"]curl(\.exe)?\b", re.IGNORECASE)
        offenders = []
        for f in SRC.rglob("*.f90"):
            if f.name == "remote_source.f90":
                continue
            if pattern.search(code(f.relative_to(SRC).as_posix())):
                offenders.append(f.name)
        self.assertEqual(offenders, [])

    def test_downloads_fail_on_http_errors(self):
        c = code("src_common/remote_source.f90")
        self.assertIn("--fail", c)
        self.assertIn("LooksLikeHtml", c)
        self.assertIn(".part", c)


if __name__ == "__main__":
    unittest.main()
