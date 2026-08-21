"""The CEC block in the essentials row is read forward, and nothing is
anchored to the end of that row any more.

RP computes the partition ratios from the raw high-frequency series and FCC
applies them to its own corrected totals, so the descriptors have to cross
between the two executables. They crossed in the FLUXNET/essentials row, and
they were found there by taking the **last eleven fields** of it.

That made "nothing may ever be appended after the descriptor" an invariant four
files had to keep by hand, with the compiler unable to help - and it could not
survive a block whose width depends on how many pairings a project declares.
So the block is self-describing now, like the hygrometer and analyser blocks
before it: a count, then that many fixed-width records, each announcing how many
targets follow it. Biomet is the genuine tail.

De-anchoring exposed something the anchor had been hiding. The skipped-period
writer emitted three fields per gas where the reader steps over seven, fourteen
per analyser where it steps over fifteen, and no hygrometer count at all. The
shortfall was silently absorbed into the biomet chunk, because the descriptor
was found by counting backwards and the misalignment never reached it. With
everything read forward, a short row is a misread row.
"""

from pathlib import Path
import re
import unittest


ROOT = Path(__file__).resolve().parents[1]

READER = "src/src_common/read_ex_record.f90"
HEADER = "src/src_rp/init_fluxnet_file_rp.f90"
ROW_RP = "src/src_rp/write_out_fluxnet.f90"
ROW_SKIPPED = "src/src_rp/write_out_fluxnet_only_biomet.f90"
ROW_FCC = "src/src_fcc/write_out_fluxnet_fcc.f90"


def read(relative):
    return (ROOT / relative).read_text(encoding="utf-8", errors="replace")


def code(relative):
    """Source with whole-line comments removed."""
    return "\n".join(ln for ln in read(relative).splitlines()
                     if not ln.lstrip().startswith("!"))


def param(source, name):
    m = re.search(r"integer, parameter :: %s\s*=\s*(\d+)" % name, source)
    assert m is not None, "%s not found" % name
    return int(m.group(1))


class NothingIsAnchoredToTheEnd(unittest.TestCase):
    def test_the_old_end_anchored_tail_is_gone(self):
        reader = code(READER)
        self.assertNotIn("nCecFields", reader)
        self.assertNotIn("remaining_fields", reader)
        #> And no writer still emits it.
        for path in (HEADER, ROW_RP, ROW_SKIPPED, ROW_FCC):
            self.assertNotIn("CEC_H2O_VALID", code(path))
            self.assertNotIn("lEx%cec%r_ET", code(path))

    def test_biomet_is_the_tail_now(self):
        """The CEC block sits ahead of biomet in the header and in every row,
        so the reader can consume it before the part whose width it does not
        know."""
        header = code(HEADER)
        self.assertLess(header.index("'NUM_CEC_PAIRS'"),
                        header.index("'NUM_BIOMET_VARS'"))

        rp = code(ROW_RP)
        self.assertLess(rp.index("AddIntDatumToDataline(n_cec_pairs"),
                        rp.index("AddIntDatumToDataline(nbVars"))

        skipped = code(ROW_SKIPPED)
        self.assertLess(skipped.index("AddIntDatumToDataline(n_cec_pairs"),
                        skipped.index("AddIntDatumToDataline(nbVars"))

        #> FCC re-emits the biomet chunk verbatim, so its CEC block has to come
        #> before that chunk rather than before a count it does not write.
        fcc = code(ROW_FCC)
        self.assertLess(fcc.index("AddIntDatumToDataline(n_cec_pairs"),
                        fcc.index("fluxnetChunks%s(6)"))


class TheBlockDescribesItsOwnWidth(unittest.TestCase):
    def setUp(self):
        self.reader = code(READER)

    def test_the_widths_are_named_and_stepped_over_by_name(self):
        self.assertEqual(param(self.reader, "nCecPairFixedFields"), 9)
        self.assertIn("strCharIndex(dataline, ',', nCecPairFixedFields)", self.reader)
        self.assertIn("strCharIndex(dataline, ',', nCecTargetFields)", self.reader)
        #> The per-target width is not pinned to a number here. It has grown
        #> once already, and the thing that matters is not what it is but that
        #> the reader, the writer and the header agree on it - which is what
        #> the two tests below check. A literal here only ever meant a second
        #> place to update.
        self.assertGreaterEqual(param(self.reader, "nCecTargetFields"), 6)

    def test_a_file_written_before_the_block_widened_is_refused(self):
        """Read one field short per target and everything after it is garbage.

        The guard is conditional on the block being there at all: a project
        with the partition switched off writes no CEC fields, and such a
        header is not old, just quiet.
        """
        vintage = code("src/src_fcc/eddyflow-fcc_main.f90")
        self.assertIn(
            "index(header, 'CEC_METH') > 0 .and. index(header, 'CEC_NS_') <= 0",
            vintage)
        self.assertIn("call ExceptionHandler(107)", vintage)

    def test_the_vintage_check_runs_before_anything_parses_a_record(self):
        """Otherwise it can never fire, and it did not.

        InitExVars reads the whole file by field position. An old file fails
        every record there, the run stops on "no valid data records found",
        and the message that would have explained why is never reached. Tested
        against a real pre-widening file, which reported error 61 instead of
        107 until the check moved ahead of it.
        """
        main = code("src/src_fcc/eddyflow-fcc_main.f90")
        checked = main.index("call CheckExFileVintageAt(AuxFile%ex)")
        parsed = main.index("call InitExVars(")
        self.assertLess(checked, parsed,
                        "the vintage check is back behind the first thing "
                        "that parses a record, where it cannot fire")
        #> Self-contained, because at that point no unit is open yet.
        self.assertIn("subroutine CheckExFileVintageAt(path)", main)
        self.assertIn("open(newunit = unt, file = path", main)
        #> One judgement shared by both entry points, so they cannot diverge.
        self.assertEqual(main.count("call JudgeExHeader(header)"), 2)

    def test_the_counts_are_bounded_before_they_are_used(self):
        #> A corrupt count is a record to reject, not a loop to run.
        self.assertIn("n_cec_pairs < 0 .or. n_cec_pairs > MaxNumCecPairs", self.reader)
        self.assertIn("n_cec_target < 0 .or. n_cec_target > MaxNumCecTargets",
                      self.reader)

    def test_every_pairing_slot_is_cleared_before_the_row_is_parsed(self):
        #> FCC walks the project's pairing list, which may be longer than what
        #> a given row carried. The surplus has to read as an empty descriptor,
        #> not as the previous record's.
        self.assertIn("call ResetCecDescriptor(lEx%cec(cec_p))", self.reader)
        self.assertIn("do cec_p = 1, MaxNumCecPairs", self.reader)

    def test_the_writer_emits_what_the_reader_steps_over(self):
        """One helper builds the row for all three writers, so the widths
        cannot drift between them."""
        resolver = code("src/src_common/gas_slot_resolution.f90")
        block = resolver[resolver.index("subroutine CecExRowValues"):]
        block = block[:block.index("end subroutine CecExRowValues")]
        self.assertEqual(block.count("call EmitInt(") + block.count("call EmitReal("),
                         param(self.reader, "nCecPairFixedFields")
                         + param(self.reader, "nCecTargetFields"),
                         "CecExRowValues and the reader disagree about how "
                         "wide a pairing is")
        for path in (ROW_RP, ROW_SKIPPED, ROW_FCC):
            self.assertIn("call CecExRowValues(", code(path))

    def test_the_header_names_as_many_fields_as_the_row_writes(self):
        header = code(HEADER)
        block = header[header.index("'NUM_CEC_PAIRS'"):header.index("'NUM_BIOMET_VARS'")]
        #> The same two numbers the reader steps over and the row helper
        #> emits. Three sources, one width, none of them a literal here.
        per_pair = block.count("call AddDatum(csv_row, 'CEC_")
        self.assertEqual(per_pair,
                         param(self.reader, "nCecPairFixedFields")
                         + param(self.reader, "nCecTargetFields"),
                         "the essentials header and the reader disagree about "
                         "how wide a pairing is")


class TheSkippedPeriodRowIsTheSameWidth(unittest.TestCase):
    """It was not, and the end anchor hid it."""

    def setUp(self):
        self.skipped = code(ROW_SKIPPED)
        self.reader = code(READER)

    def test_the_per_gas_moisture_block_is_full_width(self):
        want = param(self.reader, "nGasMoistFields")
        self.assertEqual(want, 7)
        #> The slot is real so the row stays parseable; the rest are error.
        self.assertIn("do indx = 1, %d" % (want - 1), self.skipped)

    def test_the_per_analyser_block_is_full_width(self):
        want = param(self.reader, "nGasInstrFields")
        self.assertEqual(want, 15)
        self.assertIn("do indx = 1, %d" % (want - 1), self.skipped)
        self.assertNotIn("do indx = 1, 13", self.skipped)

    def test_the_hygrometer_block_is_written_at_all(self):
        #> It was missing entirely: the writer went from the analyser loop
        #> straight to the biomet count.
        self.assertIn("call WaterOutSlots(w_slots, w_tags, n_w_slots)", self.skipped)
        self.assertIn("nWaterFluxFields = 7", self.reader)
        self.assertLess(self.skipped.index("WaterOutSlots"),
                        self.skipped.index("AddIntDatumToDataline(n_cec_pairs"))


if __name__ == "__main__":
    unittest.main()
