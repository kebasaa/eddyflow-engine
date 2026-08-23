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

One file may still name the flag, and naming it is how it gets thrown away.
An EddyPro project still carries `fix_out_format`, so the importer has to
recognise the key in order to DROP it; without that case the key would be
copied through into the converted EddyFlow project, where it would sit as a
setting nothing reads - exactly the orphan `col_co2` and friends are listed
there to avoid. That is a widening of what a reader accepts, not a use of the
flag, and it is the same distinction test_legacy_slot_scope_static.py draws
between an on-disk alias and a rename.

So the allowance is narrow in both directions: one occurrence, and it has to
be inside IsConsumedKey. A read anywhere - including anywhere else in that
same file - still fails.
"""

from pathlib import Path
import re
import subprocess
import unittest


ROOT = Path(__file__).resolve().parents[1]

TABLE = "src/src_common/m_common_global_var.f90"
IMPORT = "src/src_common/m_eddypro_import.f90"

#: Occurrences outside whole-line comments, per file, and why each is there.
#: A number going up is a new read and fails. A number going down is the
#: allowance outliving its reason - lower it, so the file cannot quietly
#: regain what it gave up.
ALLOWED = {
    #> IsConsumedKey, where naming a retired key is how the importer discards
    #> it rather than copying it into the file it writes. Pinned to that
    #> function by test_the_importer_only_names_it_to_discard_it below.
    IMPORT: 1,
}


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


def occurrences(rel):
    """How many times a file names the flag outside a whole-line comment."""
    return sum(ln.count("fix_out_format") for ln in code(rel).splitlines())


def function_body(rel, opening, closing):
    src = code(rel)
    start = src.index(opening)
    return src[start:src.index(closing, start) + len(closing)]


class TheFlagIsGone(unittest.TestCase):
    def test_no_source_reads_it(self):
        offenders = []
        for rel in tracked_sources():
            n = occurrences(rel)
            allowed = ALLOWED.get(rel, 0)
            if n > allowed:
                offenders.append("%s: %d, expected at most %d" % (rel, n, allowed))
        self.assertFalse(
            offenders,
            "these read the retired fixed-format flag:\n  "
            + "\n  ".join(offenders)
            + "\n\nThe full output covers every configured gas. There is no "
              "second column set to select between. The one legitimate reason "
              "to name the key is to discard it on import; if that is what "
              "this is, say so in a comment and raise the count in ALLOWED.")

    def test_the_allowance_is_not_stale(self):
        """An allowance that has stopped being needed is one that could be
        spent on a real read."""
        stale = []
        for rel, allowed in sorted(ALLOWED.items()):
            if not (ROOT / rel).is_file():
                stale.append("%s: listed but does not exist" % rel)
                continue
            n = occurrences(rel)
            if n < allowed:
                stale.append("%s: %d now, expectation still %d" % (rel, n, allowed))
        self.assertFalse(
            stale,
            "lower these expectations to match the tree:\n  " + "\n  ".join(stale))

    def test_the_importer_only_names_it_to_discard_it(self):
        """The allowance is for one function, not for the file.

        IsConsumedKey answers "is this a key the records replace, and so must
        not be copied". Naming the flag anywhere else in the importer would be
        reading it, and would still be wrong.
        """
        body = function_body(IMPORT,
                             "logical function IsConsumedKey(",
                             "end function IsConsumedKey")
        self.assertEqual(body.count("fix_out_format"), 1,
                         "the importer names the flag outside IsConsumedKey")
        self.assertIn("case ('fix_out_format')", body)
        #> And it is a discard, not a decode: the arm returns without reading
        #> a value, exactly as the col_* arm above it does.
        arm = body[body.index("case ('fix_out_format')"):]
        arm = arm[:arm.index("end select")]
        self.assertEqual(arm.split(), ["case", "('fix_out_format')", "return"],
                         "the fix_out_format arm does more than discard the key")

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
        #> A floor, not an equality: the rule being pinned is that the table
        #> never shrinks. Asserting the exact size made every new [Project]
        #> key edit this line, which is noise around the one thing that
        #> matters - a tag index moving under settings already written to
        #> files.
        self.assertGreaterEqual(
            int(m.group(1)), 274,
            "Npc shrank. Retiring a tag blanks its label and leaves the "
            "size alone; the generator never shrinks a table.")

    def test_the_generator_owns_the_retirement(self):
        """Hand-editing the block is how it drifts from the generator."""
        gen = (ROOT / "prj/gen_project_tags.py").read_text(encoding="utf-8")
        self.assertIn('"fix_out_format",', gen,
                      "the label must be retired through RETIRED_LABELS, so "
                      "re-running the generator keeps it blank")


class NoFixtureStillSetsIt(unittest.TestCase):
    def test_no_project_file_carries_the_key(self):
        #> Fixtures only. run.sh writes run_ref/run_chk copies of whichever
        #> fixture it is running into the same directory - they are gitignored
        #> scratch, and one of them carries whatever the fixture it was copied
        #> from carried. Counting them made this check pass or fail on whether
        #> anyone had run the harness lately, which is not a property of the
        #> repository.
        offenders = [p.name for p in (ROOT / "tests/regression").glob("*.eddyflow")
                     if not p.name.startswith("run_")
                     and "fix_out_format" in p.read_text(encoding="utf-8",
                                                        errors="replace")]
        self.assertFalse(
            offenders,
            "these fixtures still set a key nothing reads:\n  "
            + "\n  ".join(sorted(offenders)))


if __name__ == "__main__":
    unittest.main()
