"""A home directory with a quote in it stops the run before anything is made.

`-e "C:\\dir\\"` from a shell reaches the engine as `C:\\dir" -j 3 <project>`:
the trailing backslash escapes the closing quote, and the C runtime folds every
argument after it into the home directory. InitEnv then built TmpDir from that
and CreateDir handed cmd `mkdir "C:\\dir" -j 3\\tmp\\tmp_<ts>"`, which cmd reads
as several folders - so `-j`, `2`, `4` and `3\\tmp\\tmp_<ts>` turned up in the
repo root and in prj/, wherever the run had been started from.

The guard has to sit before the first path is built from the home directory
and before anything is created, and it has to stop: the project path has
already been swallowed, so falling back to the default home would quietly run
whatever project that default names.
"""

from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]
INIT_ENV = "src/src_common/init_env.f90"


def code(path):
    text = (ROOT / path).read_text(encoding="utf-8")
    return "\n".join(ln for ln in text.splitlines()
                     if not ln.lstrip().startswith("!"))


class AQuotedHomeDirectoryStops(unittest.TestCase):
    GUARD = """if (index(homedir, '"') > 0 .or. index(projPath, '"') > 0) then"""

    def setUp(self):
        self.body = code(INIT_ENV)

    def test_both_paths_are_tested_for_a_quote(self):
        self.assertIn(self.GUARD, self.body)

    def test_it_comes_before_any_path_is_built_or_created(self):
        guard = self.body.index(self.GUARD)
        self.assertLess(guard, self.body.index("call AdjDir(homedir"))
        self.assertLess(guard, self.body.index("CreateDir("))

    def test_it_says_so_on_stderr_and_stops(self):
        guard = self.body.index(self.GUARD)
        block = self.body[guard:self.body.index("end if", guard)]
        self.assertIn("write(error_unit", block,
                      "the run log is not open yet; a write to its unit "
                      "would land in a fort.NNN file")
        self.assertIn("stop 3", block)
        self.assertNotIn("homedir = ", block,
                         "falling back to a default home would run another "
                         "project, the real path having been swallowed")


if __name__ == "__main__":
    unittest.main()
