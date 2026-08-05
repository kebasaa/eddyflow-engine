"""The fixed full-output format is gone, and its tag slot is blanked.

`fix_out_format` selected a literal column set naming exactly four gas blocks
whatever the project held. It was a compatibility promise, and three things
were true of it by the end:

  * a fifth gas had no columns in that format and fell out of the full output
    altogether - the file the flag exists to stabilise was the one file that
    could not describe the project;
  * its two header branches had already drifted. FCC named both analysers from
    a single `co2_new_sw_ver` test where RP asked each one's own firmware, and
    both spelled `mean_value_RSSI_LI-7200` into the header text;
  * the format it reproduced was attributed to "EddyPro 7.x" by this fork's
    own comments and was never checked against a 7.x file in the tree.

So it is retired rather than repaired, and the full output covers every
configured gas.

Two things this pins. First, that no code reads the flag again - the field is
gone from the project type, so a reference would not compile, but the *tag*
survives and could be re-decoded. Second, and the reason this file exists at
all: that slot 37 of EPPrjCTags is BLANKED and not DELETED. These tables are
positional. Deleting an entry renumbers every tag after it and silently
rebinds hundreds of settings, with nothing to catch it - no compile error, no
crash, just a project file that means something different than it says.
"""

from pathlib import Path
import re
import subprocess
import unittest


ROOT = Path(__file__).resolve().parents[1]

TABLE = "src/src_common/m_common_global_var.f90"


def tracked_sources():
    out = subprocess.run(["git", "ls-files", "src"], cwd=ROOT,
                         capture_output=True, text=True, check=True).stdout
    return [p for p in out.split()
            if p.endswith(".f90") or p.endswith(".inc")]


def code(rel):
    """Source with full-line comments removed.

    A plain `in` against Fortran source almost always wants the comments gone
    first, or the assertion matches the note explaining the removal.
    """
    return "\n".join(ln for ln in (ROOT / rel).read_text(
        encoding="utf-8", errors="replace").splitlines()
        if not ln.lstrip().startswith("!"))


class TheFlagIsGone(unittest.TestCase):
    def test_no_source_reads_it(self):
        offenders = [rel for rel in tracked_sources()
                     if "fix_out_format" in code(rel)]
        self.assertFalse(
            offenders,
            "these read the retired fixed-format flag:\n  "
            + "\n  ".join(offenders)
            + "\n\nThe full output covers every configured gas. There is no "
              "second column set to select between.")

    def test_the_project_type_no_longer_carries_it(self):
        self.assertNotIn("logical :: fix_out_format",
                         code("src/src_common/m_typedef.f90"))

    def test_the_helper_has_one_loop(self):
        """FullOutputGasSlots chose between four slots and all of them."""
        body = code("src/src_common/gas_slot_resolution.f90") \
            .split("subroutine FullOutputGasSlots")[1].split("end subroutine")[0]
        self.assertEqual(body.count("do gas ="), 1,
                         "the fixed arm is back; the full output has one list")
        self.assertNotIn("histGas", body,
                         "the list names four historical slots again")


class TheTagSlotIsBlankedNotDeleted(unittest.TestCase):
    """The tables are positional. This is the whole point of the file."""

    def test_slot_37_is_present_and_empty(self):
        self.assertRegex(
            code(TABLE), r"EPPrjCTags\(37\)%Label\s*/\s*''\s*/",
            "slot 37 must keep its place with an empty label - deleting it "
            "renumbers every tag after it")

    def test_its_neighbours_did_not_move(self):
        """What a renumbering would look like, caught by name."""
        src = code(TABLE)
        for idx, label in ((36, "err_label"), (38, "qc_meth")):
            self.assertRegex(
                src, r"EPPrjCTags\(%d\)%%Label\s*/\s*'%s'\s*/" % (idx, label),
                "EPPrjCTags(%d) should still be '%s'; if it is not, the table "
                "has shifted and every project setting past slot 37 now means "
                "something else" % (idx, label))

    def test_the_table_did_not_shrink(self):
        src = code(TABLE)
        m = re.search(r"integer, parameter :: Npc = (\d+)", src)
        self.assertIsNotNone(m)
        self.assertEqual(int(m.group(1)), 274,
                         "Npc changed. Retiring a tag blanks its label and "
                         "leaves the size alone; the generator never shrinks "
                         "a table.")

    def test_the_generator_owns_the_retirement(self):
        """Hand-editing the block is how it drifts from the generator."""
        gen = (ROOT / "prj/gen_project_tags.py").read_text(encoding="utf-8")
        self.assertIn('"fix_out_format",', gen,
                      "the label must be retired through RETIRED_LABELS, so "
                      "re-running the generator keeps it blank")


class NoFixtureStillSetsIt(unittest.TestCase):
    def test_no_project_file_carries_the_key(self):
        offenders = [p.name for p in (ROOT / "tests/regression").glob("*.eddyflow")
                     if "fix_out_format" in p.read_text(encoding="utf-8",
                                                        errors="replace")]
        self.assertFalse(
            offenders,
            "these fixtures still set a key nothing reads:\n  "
            + "\n  ".join(sorted(offenders)))


if __name__ == "__main__":
    unittest.main()
