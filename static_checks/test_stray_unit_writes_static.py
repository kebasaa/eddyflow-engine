"""A write to a unit nobody opened does not fail - it goes to fort.NN.

gfortran connects an unconnected unit to ``fort.<unit>`` in whatever directory
the run started in, and the write succeeds. So this shape

    open(u, file = path, iostat = open_status)
    if (open_status /= 0) call ExceptionHandler(64)

    write(u, '(a)') ...

does not do what it reads as. When the open fails the handler logs, returns,
and every write below it lands in ``fort.<u>`` under a name carrying no
meaning - ``udf`` is the *millisecond the run started in*
(``init_env.f90``), so it is a different name every run.

That is not hypothetical. ``fort.615`` sat in the repository root holding 36 KB
of ensemble cospectra and Massman fits while ``Error(64)`` told the user
"Some spectral assessment results will not be written on output file". They
were written, in full, to a file nobody would look for.

This is the **third** time the repo has been bitten by the same class:

* unit 11 - ``edit_ini_file.f90``, fixed 2026-08-17 in ``7ed425c``; the failing
  path also deleted its input, so the three ``fort.11`` files it left behind
  were the only surviving copies of what those runs destroyed.
* unit 163 - the run log, guarded by ``test_run_log_static.py`` whose docstring
  names this exact mechanism.
* unit ``udf`` - the five sites here.

Two things are checked. That no ``fort.*`` file is in the tree, which is the
symptom and is now visible again because ``.gitignore`` no longer hides it. And
that an ``open`` whose failure is handled by a *non-fatal* exception cannot fall
through into writes.

The distinction matters: most sites in the codebase pair a failed open with
``ExceptionHandler(60)``, which ends in ``stop 1``. Those are already safe and
must not be flagged - banning the shape outright would be wrong.
"""

import re
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / "src"

HANDLER = (SRC / "src_common" / "exception_handler.f90").read_text(
    encoding="utf-8", errors="replace")


def fatal_codes():
    """Exception codes whose arm ends the program.

    Parsed rather than listed, so a handler that stops today and is changed to
    warn tomorrow moves into scope here on its own.
    """
    fatal = set()
    arms = re.split(r"\n\s*case\s*\((\d+)\)", HANDLER)
    for code, body in zip(arms[1::2], arms[2::2]):
        if re.search(r"^\s*(error )?stop\b", body, re.MULTILINE):
            fatal.add(int(code))
    return fatal


FATAL = fatal_codes()


class NoStrayUnitFilesAreInTheTree(unittest.TestCase):
    """The symptom. Visible again now that .gitignore does not hide it."""

    def test_no_fort_files(self):
        strays = [p for p in ROOT.rglob("fort.*") if ".git" not in p.parts]
        self.assertEqual(
            [], [str(p.relative_to(ROOT)) for p in strays],
            "a fort.NN file is gfortran connecting a unit nobody opened; find "
            "the unguarded write rather than deleting the file")

    def test_the_ignore_rule_was_not_put_back(self):
        text = (ROOT / ".gitignore").read_text(encoding="utf-8")
        self.assertNotIn(
            "\nfort.*", text,
            "ignoring fort.* turns a loud failure into a silent one")


class AFailedOpenCannotFallThroughIntoWrites(unittest.TestCase):
    """The cause. Only non-fatal handlers are in scope."""

    #: open, then a guard that calls a handler and does nothing else.
    PATTERN = re.compile(
        r"open\s*\([^\n]*iostat\s*=\s*(\w+)[^\n]*\)\s*\n"
        r"\s*if\s*\(\s*\1\s*/=\s*0\s*\)\s*call\s+ExceptionHandler\s*\(\s*(\d+)\s*\)"
        r"\s*\n")

    def test_every_non_fatal_guard_is_followed_by_a_stop_or_a_skip(self):
        offenders = []
        for path in sorted(SRC.rglob("*.f90")):
            text = path.read_text(encoding="utf-8", errors="replace")
            for m in self.PATTERN.finditer(text):
                code = int(m.group(2))
                if code in FATAL:
                    continue
                #> Non-fatal: the statement after the guard has to be a return,
                #> or the guard has to be a block that does something about it.
                after = text[m.end():m.end() + 200].lstrip()
                if after.startswith(("return", "cycle", "exit")):
                    continue
                line = text[:m.start()].count("\n") + 1
                offenders.append(
                    "%s:%d (ExceptionHandler(%d) does not stop)"
                    % (path.relative_to(ROOT), line, code))
        self.assertEqual(
            [], offenders,
            "a failed open followed by writes sends them to fort.NN")

    def test_the_fatal_set_was_actually_found(self):
        """Guards the parser above: an empty set would pass everything."""
        self.assertIn(60, FATAL, "ExceptionHandler(60) ends in stop 1")
        self.assertNotIn(64, FATAL, "ExceptionHandler(64) only logs")


class TheSpectralOutputsGoThroughTheGuardedOpen(unittest.TestCase):
    """The five sites that produced fort.615."""

    SOURCE = (SRC / "src_fcc" / "output_spectral_assessment_results.f90").read_text(
        encoding="utf-8", errors="replace")

    def test_no_site_opens_udf_directly(self):
        body = self.SOURCE[:self.SOURCE.index(
            "subroutine OpenSpectralOutputFile")]
        self.assertNotIn("open(udf", body,
                         "every output file goes through OpenSpectralOutputFile")

    def test_all_five_use_the_helper(self):
        self.assertEqual(self.SOURCE.count("call OpenSpectralOutputFile("), 5)

    def test_the_helper_sinks_a_failure_rather_than_returning(self):
        """Returning would skip the other four files, which may be fine."""
        body = self.SOURCE[self.SOURCE.index("subroutine OpenSpectralOutputFile"):]
        self.assertIn("open(udf, status = 'scratch')", body)
        self.assertIn("call ExceptionHandler(64)", body)


if __name__ == "__main__":
    unittest.main()
