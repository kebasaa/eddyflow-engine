"""Every object has a rule, and the rules are generated rather than typed.

A make target with no prerequisites is up to date the moment the file exists.
So a source file missing from the Makefile's dependency list is not a missing
optimisation - it is a source that, once compiled, is never compiled again.
The edit lands, the object does not change, and the previous compile is linked
in its place with nothing said.

That is not hypothetical here. It happened to gas_slot_resolution.o, was fixed
by hand, and the fix left a comment. It then happened to parse_month_grouping.o
and surfaced only as an undefined reference at link time - which is luck: the
edit happened to add a symbol. An edit that only changed the body of an
existing routine would have produced a binary quietly built from the old code,
and every regression run against it would have been meaningless.

So the list is generated now, and this is what stops it drifting again. It
also checks the two things a generated list is supposed to guarantee and a
hand-written one never did: that every source has an entry, and that the entry
names the source itself.
"""

from pathlib import Path
import re
import subprocess
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
MAKEFILE = ROOT / "prj" / "Makefile"
GENERATOR = ROOT / "prj" / "gen_makefile_deps.py"

SOURCE_GLOBS = ("src/src_rp/*.f90", "src/src_rp/fft4/*.F",
                "src/src_fcc/*.f90", "src/src_common/*.f90")


def rules():
    """object stem -> set of prerequisites, continuations joined."""
    text = MAKEFILE.read_text(encoding="utf-8", errors="replace")
    text = text.replace("\r\n", "\n")
    text = re.sub(r"\\\n\s*", " ", text)
    out = {}
    for line in text.split("\n"):
        m = re.match(r"^(\S+)\.o:(?!=)(.*)$", line)
        if m and not line.startswith("#"):
            out[m.group(1)] = set(m.group(2).split())
    return out


def sources():
    out = {}
    for pattern in SOURCE_GLOBS:
        for path in ROOT.glob(pattern):
            out[path.stem] = path
    return out


class EverySourceHasARule(unittest.TestCase):
    def test_no_object_is_missing(self):
        missing = sorted(set(sources()) - set(rules()))
        self.assertFalse(
            missing,
            "these sources have no dependency rule:\n  "
            + "\n  ".join(missing)
            + "\n\nMake treats a target with no prerequisites as up to date, "
              "so each of these would be compiled once and never again. Run "
              "prj/gen_makefile_deps.py.")

    def test_no_rule_names_a_source_that_is_gone(self):
        """A stale entry is harmless to the build and misleading to a reader.

        There was one: user_timelag_handle.o, whose source had been deleted.
        """
        stale = sorted(set(rules()) - set(sources()))
        self.assertFalse(
            stale,
            "these rules name sources that do not exist:\n  "
            + "\n  ".join(stale))

    def test_every_rule_depends_on_its_own_source(self):
        """The dependency that makes an edit rebuild anything at all."""
        src = sources()
        wrong = []
        for stem, deps in rules().items():
            if stem not in src:
                continue
            if src[stem].name not in deps:
                wrong.append("%s.o does not depend on %s" % (stem, src[stem].name))
        self.assertFalse(wrong, "\n  ".join(wrong))


class TheIncludesAreFollowed(unittest.TestCase):
    """interfaces.inc includes interfaces_1.inc, so the chain nests.

    The hand-written list had no include dependencies at all, so editing
    either file rebuilt nothing at all. Resolving only the first level would
    still miss every source that reaches interfaces_1.inc through
    interfaces.inc, which is half of them.
    """

    def test_a_source_including_the_outer_file_depends_on_the_inner_one(self):
        by_stem = rules()
        outer = [s for s, p in sources().items()
                 if "interfaces.inc" in p.read_text(encoding="utf-8",
                                                    errors="replace")]
        self.assertTrue(outer, "nothing includes interfaces.inc any more")
        for stem in outer:
            self.assertIn(
                "interfaces_1.inc", by_stem.get(stem, set()),
                "%s.o includes interfaces.inc, which includes "
                "interfaces_1.inc, but does not depend on it" % stem)


class TheBlockIsGenerated(unittest.TestCase):
    def test_the_generator_agrees_with_the_makefile(self):
        self.assertTrue(GENERATOR.exists(), GENERATOR)
        r = subprocess.run([sys.executable, str(GENERATOR), "--check"],
                           capture_output=True, text=True, cwd=str(ROOT / "prj"))
        self.assertEqual(
            r.returncode, 0,
            "the Makefile dependency block is stale; re-run "
            "prj/gen_makefile_deps.py\n%s%s" % (r.stdout, r.stderr))


if __name__ == "__main__":
    unittest.main()
