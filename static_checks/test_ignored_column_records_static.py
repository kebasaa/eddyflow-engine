"""A project record naming an ignored column must not select it.

The metadata says what a column *is*; the project only says which columns it
wants. `ignore` and `not_numeric` mean the column holds nothing usable, and the
import drops such columns outright - so a record naming one could never resolve
to data.

DefineUsedVariables marked a column used purely because a record pointed at it,
never asking what the metadata declared it to be. MetadataFileValidation then
found a used column declared `ignore` and killed the run. The engine was arguing
with itself: marking a column, then dying because it was marked.

That is reachable from an ordinary edit. Moving a diagnostic from one column to
another leaves a record behind pointing at the old one; the old column is then
marked `ignore` in the metadata, and every subsequent run aborts on a column the
user had explicitly told the engine to ignore.
"""

from pathlib import Path
import re
import unittest


ROOT = Path(__file__).resolve().parents[1]
USED = "src/src_common/define_used_variables.f90"
VALIDATION = "src/src_common/metadata_file_validation.f90"
INFORM = "src/src_common/inform_of_metadata_problem.f90"


def read(path):
    return (ROOT / path).read_text(encoding="utf-8")


def body_of(text, name):
    """The text of a subroutine or function, by name."""
    start = text.index(name)
    return text[start:]


class RecordColumnsHonourTheMetadata(unittest.TestCase):
    def test_the_predicate_exists(self):
        source = read(USED)
        self.assertIn("logical function ColumnIsSelectable(col_num)", source)

    def test_it_rejects_both_spellings(self):
        """`ignore` and `not_numeric` are the two the validation tests for.

        Covering only one leaves the other reaching validation, which is the
        same fatal by a different route.
        """
        source = read(USED)
        body = body_of(source, "logical function ColumnIsSelectable")
        body = body[: body.index("end function ColumnIsSelectable")]
        self.assertIn("'ignore', 'not_numeric'", body)
        self.assertIn("call lowercase(var)", body)

    def test_every_record_loop_goes_through_it(self):
        """Gas, cell and diagnostic alike.

        These were three copies of the same unguarded pattern, which is how one
        of them being wrong meant all three were. A new record kind copied from
        its neighbours inherits the guard along with the shape.
        """
        source = read(USED)
        for kind, cap in (
            ("gas", "MaxNumGases"),
            ("cell", "MaxNumCellCols"),
            ("diag", "MaxNumDiagCols"),
        ):
            loop = re.search(
                r"do i = 1, min\(EddyFlowProj%%%s_num, %s\)\s*\n(.*?)\n\s*end do"
                % (kind, cap),
                source,
                re.S,
            )
            self.assertIsNotNone(loop, "the %s record loop is gone" % kind)
            self.assertIn(
                "ColumnIsSelectable(EddyFlowProj%%%s(i)%%col)" % kind,
                loop.group(1),
                "the %s record loop marks a column without asking whether the "
                "metadata allows it to be selected" % kind,
            )

    def test_the_slot_array_is_guarded_too(self):
        """The bare `where` blocks marked from the slot array unguarded.

        ApplyDiagnosticRecordColumns fills those slots before any metadata is
        read and is last-writer-wins, so a record naming an ignored column can
        take the slot from a usable one. Guarding only the record loops leaves
        that path marking the column and the run still dying on it.
        """
        source = read(USED)
        self.assertNotIn(
            "where (EddyFlowProj%Col(firstGas:E2NumVar) > 0)",
            source,
            "the slot array is marked without consulting the metadata",
        )
        self.assertIn("if (ColumnIsSelectable(EddyFlowProj%Col(i))) &", source)

    def test_the_diagnostic_slots_are_re_resolved(self):
        """And the winner has to be a record that can actually be read.

        Clearing first is what makes it correct: a slot whose only record names
        an ignored column must come back empty, not keep what the pre-metadata
        pass left in it, or the presence test reports a diagnostic the file does
        not carry.
        """
        source = read(USED)
        self.assertIn("if (EddyFlowProj%diag_num > 0) then", source)
        self.assertIn(
            "if (.not. ColumnIsSelectable(EddyFlowProj%Col(i))) &", source
        )
        self.assertIn(
            "if (.not. ColumnIsSelectable(EddyFlowProj%diag(i)%col)) cycle",
            source,
        )


class RecordsCompetingForOneSlotAreRefused(unittest.TestCase):
    """The engine holds one slot per diagnostic kind; two records overwrote.

    The loser vanished with nothing said, and which record lost turned on their
    order in the file - so a stale record was inert until the day it won.
    """

    def test_the_check_exists_and_has_its_own_flag(self):
        source = read(VALIDATION)
        self.assertIn("passed(27) = .false.", source)
        self.assertIn("logical function RecordIsLive", source)

    def test_inert_records_do_not_count(self):
        """A record on an ignored column is already inert.

        Counting it would refuse the very projects the ignore handling exists to
        keep running. So is a record whose column has since been re-declared as
        something else - see test_signal_strength_records_static.py, which is
        the case that made that matter.
        """
        source = read(VALIDATION)
        self.assertIn(
            "if (.not. RecordIsLive(LocCol, EddyFlowProj%diag(i))) cycle",
            source,
        )

    def test_cells_are_keyed_by_instrument_too(self):
        """Cell records repeat their variable by design - one cell_t per
        analyser - so the key is (variable, instrument), not the variable.

        Keying on the name alone refuses any two-analyser project, which is a
        normal configuration and not a duplicate at all.
        """
        source = read(VALIDATION)
        cells = source[source.index("EddyFlowProj%cell_num, MaxNumCellCols"):]
        self.assertIn("%cell(i)%instr", cells)
        self.assertIn("%cell(j)%instr", cells)

    def test_gases_are_exempt(self):
        """Their slot is the record index, so they cannot collide.

        Two CO2 records are a second analyser, not a duplicate.
        """
        source = read(VALIDATION)
        tail = source[source.index("passed(27)"):]
        self.assertNotIn("EddyFlowProj%gas_num", tail)

    def test_the_message_names_both_records(self):
        """Naming only the column leaves the user to find the other side.

        The interface shows one row; the duplicate is visible only in the
        project file, so the message has to say which records to look at.
        """
        source = read(INFORM)
        self.assertIn("if (.not. passed(27)) then", source)
        block = source[source.index("if (.not. passed(27)) then"):]
        self.assertIn("Diagnostic record ", block)
        self.assertIn("Cell record ", block)


if __name__ == "__main__":
    unittest.main()
