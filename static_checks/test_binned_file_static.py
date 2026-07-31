"""The binned (co)spectra file carries every configured gas.

It is the sole input to the on-the-fly spectral assessment, and it was fixed
at 18 columns - three frequencies, eight spectra, seven cospectra - on both
sides. So a gas past the fourth could never be assessed however wide the
loops downstream were: its data were never written, and the reader had no
column to put them in.

Writer and reader now name columns from the same helper, and the reader
matches by name rather than by position. That matters more here than
elsewhere: a project points sa_bin_spectra at a directory kept from an
earlier run, so reading a file written by a different build is the normal
case rather than the exception.
"""

from pathlib import Path
import re
import unittest


ROOT = Path(__file__).resolve().parents[1]

WRITER = "src/src_rp/spectral_analysis.f90"
READER = "src/src_fcc/read_binned_file.f90"


def read(path):
    return (ROOT / path).read_text(encoding="utf-8")


def code(path):
    return "\n".join(ln for ln in read(path).splitlines()
                     if not ln.lstrip().startswith("!"))


class WriterAndReaderShareTheirColumnNames(unittest.TestCase):
    def test_both_build_names_from_the_same_helper(self):
        for path in (WRITER, READER):
            self.assertIn("call SpectralVarTags", code(path),
                          "%s must name columns from the records" % path)

    def test_the_writer_names_no_gas_literally(self):
        """The header was a literal naming u..ch4 with only slot 8
        substituted - and substituted from SpecCol%label, a third spelling
        that is uppercase where every other file is lower."""
        source = code(WRITER)
        for literal in ("spec(co2)", "spec(h2o)", "spec(ch4)", "og(co2)",
                        "cospec(w_co2)", "SpecCol(gas4)%label"):
            self.assertNotIn(literal, source,
                             "%s still writes %s as a literal"
                             % (WRITER, literal))

    def test_the_reader_finds_the_header_rather_than_counting_lines(self):
        """Counting preamble lines breaks the moment the preamble changes
        length, and the reader this replaces discarded it only by letting the
        numeric read fail - which cannot tell a header from a corrupt row."""
        source = code(READER)
        self.assertIn("'#_freq'", source)
        self.assertNotIn("read(udf, *, iostat = read_status) BinSpec(i)%fnum",
                         source,
                         "the fixed 18-item positional read is back")

    def test_the_reader_matches_case_insensitively(self):
        """An older binned directory spells the fourth gas however the
        metadata did - 'COS' where the current writer says 'cos'."""
        source = code(READER)
        self.assertIn("call uppercase(probe)", source)
        self.assertIn("call uppercase(speclabs(j))", source)

    def test_an_unmatched_slot_declines(self):
        """ErrSpec, not a neighbour's column. A four-gas directory read by a
        project with eight must leave the last four unassessed."""
        source = code(READER)
        self.assertIn("BinSpec = ErrSpec", source)
        self.assertIn("if (spec_ord(j) > 0)", source)
        self.assertIn("if (cosp_ord(j) > 0)", source)

    def test_a_file_with_no_header_is_reported_not_silently_empty(self):
        """Returning nbins = 0 reads as "no data this period" and the
        assessment quietly shrinks."""
        source = code(READER)
        block = source[source.index("'#_freq'"):]
        self.assertIn("call ExceptionHandler(62)", block)


class TheAssessmentChainCoversEveryGas(unittest.TestCase):
    def test_the_cospectra_accumulator_is_not_four_bounded(self):
        source = code("src/src_fcc/cospectra_sorting_and_averaging.f90")
        self.assertNotIn("do gas = w_ts, w_gas4", source,
                         "ensemble cospectra a later gas contributed were "
                         "read and then discarded here")

    def test_the_flux_candidate_gate_covers_every_gas(self):
        source = code("src/src_fcc/spectral_assessment_diagnostics.f90")
        self.assertNotIn("if (gas < co2 .or. gas > gas4) return", source,
                         "the readiness report could only ever say flux=0 "
                         "for a gas past the fourth")

    def test_the_razor_blade_is_one_loop_with_a_sentinel_guard(self):
        """Four hand-written blocks per stability case, and thresholds that
        exist only for co2/ch4/gas4 plus whatever has an sa_* record.

        A threshold left at zero reads as "accept everything" for a minimum
        and "reject everything" for a maximum - neither a decision the
        project made - so unset means the test is skipped.
        """
        source = code("src/src_fcc/cospectra_qaqc.f90")
        self.assertEqual(source.count("do gas = firstGas, lastGas"), 2,
                         "one loop per stability case")
        self.assertIn("lo /= error .and. gas_flux < lo", source)
        self.assertIn("hi /= error .and. gas_flux > hi", source)
        self.assertIn("GasSlotIsWater(gas)", source,
                      "water is tested on its latent heat flux")

        init = " ".join(code("src/src_fcc/read_ini_fcc.f90").split())
        for name in ("min_un_gas", "min_st_gas", "max_gas"):
            self.assertIn("FCCsetup%%SA%%%s = error" % name, init,
                          "%s must start at the sentinel, before the legacy "
                          "slots are assigned" % name)


if __name__ == "__main__":
    unittest.main()
