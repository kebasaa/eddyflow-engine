"""libgfortran's locale bracket reads freed memory under the UCRT.

Without newlocale, which Windows does not have, libgfortran wraps every
formatted READ and WRITE in

    old_locale = setlocale(LC_NUMERIC, NULL)
    setlocale(LC_NUMERIC, "C")
    ...
    setlocale(LC_NUMERIC, old_locale)

and keeps the pointer rather than a copy. Under the UCRT that string belongs
to the locale data the "C" call frees, so the restore reads freed memory.
Usually it still says "C"; about one run in eight of pwb_table_check the heap
had reused it and setlocale faulted inside ucrtbase.dll - MAIN__,
StorePwbTimelagCache, LocatePwbRow, PeriodMinutes's internal READ,
finalize_transfer, setlocale. GNU Fortran 15.2.0, MinGW-w64 UCRT.

src/src_os/win/setlocale_override.f90 replaces setlocale in every executable
that links libgfortran statically: a query answers "C" from storage that is
never freed, "C" is a no-op, anything else is refused. The program never
leaves "C", so that is all it needs.
"""

from pathlib import Path
import re
import unittest


ROOT = Path(__file__).resolve().parents[1]
OVERRIDE = ROOT / "src" / "src_os" / "win" / "setlocale_override.f90"
POSIX = ROOT / "src" / "src_os" / "posix"
MAKEFILE = ROOT / "prj" / "Makefile"


def code(path):
    return "\n".join(ln for ln in path.read_text(encoding="utf-8").splitlines()
                     if not ln.lstrip().startswith("!"))


class TheOverrideIsWhatTheLinkerBinds(unittest.TestCase):

    def test_it_is_the_c_symbol(self):
        self.assertRegex(code(OVERRIDE),
                         r"bind\(C, *name *= *'setlocale'\)")

    def test_it_is_a_windows_source(self):
        """On POSIX libgfortran has newlocale and never calls setlocale, and
        replacing libc's would affect everything else that calls it."""
        for f in POSIX.glob("*.f90"):
            self.assertNotIn("'setlocale'", f.read_text(encoding="utf-8"), f.name)

    def test_the_makefile_compiles_it(self):
        self.assertIn("setlocale_override.o:", MAKEFILE.read_text(encoding="utf-8"))


class AQueryNeverPointsAtMemoryThatCanBeFreed(unittest.TestCase):

    def setUp(self):
        self.body = code(OVERRIDE)

    def test_the_answer_lives_in_saved_storage(self):
        self.assertRegex(self.body, r"target, *save *:: *c_name\(2\)")
        self.assertIn("current = c_loc(c_name)", self.body)

    def test_nothing_calls_the_runtime(self):
        """Forwarding to the real setlocale would bring the freeing back."""
        for forbidden in ("GetProcAddress", "_wsetlocale", "__imp_setlocale"):
            self.assertNotIn(forbidden, self.body)

    def test_only_c_is_accepted(self):
        self.assertIn("if (wanted(1) /= 'C') return", self.body)
        self.assertIn("if (wanted(2) /= c_null_char) return", self.body)
        #> The refusal is the null pointer the result starts as.
        self.assertLess(self.body.index("current = c_null_ptr"),
                        self.body.index("if (wanted(1) /= 'C') return"))


if __name__ == "__main__":
    unittest.main()
