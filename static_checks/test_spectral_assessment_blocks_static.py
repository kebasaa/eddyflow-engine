"""The spectral assessment file's per-gas blocks are matched by name.

Each block in that file opens with a header naming its gas -
"CO2_2            TFP            Fn          fc" - and the reader used to
ignore that name entirely. It walked the expected gas list and assigned the
Nth block to the Nth non-water slot, testing only that the header contained
the word TFP. A file whose block set differed from what the project expects,
by even one gas, had every block after the difference assigned to the wrong
species, silently and with plausible values.

That is not hypothetical. It is what blocks widening the water carve-out
here: every file written so far contains a block for a second hygrometer, so
a reader that skipped both waters would expect one block fewer than the file
holds and shift everything after it.

Proven by reordering: swapping the N2O and CO2_2 blocks in
sa_n_gas_fitted.txt leaves every correction factor unchanged (FN2O_SCF
4.05783, FCO2_2_SCF 4.06799 either way). Under the positional reader those
two would have traded places.
"""

from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]

READER = "src/src_fcc/read_spectral_assessment_file.f90"
WRITER = "src/src_fcc/output_spectral_assessment_results.f90"
GATE = "src/src_fcc/spectral_assessment_diagnostics.f90"


def code(path):
    return "\n".join(ln for ln in (ROOT / path).read_text(encoding="utf-8",
                                                          errors="replace").splitlines()
                     if not ln.lstrip().startswith("!"))


class BlocksAreResolvedByName(unittest.TestCase):
    def test_both_sides_name_blocks_from_the_same_helper(self):
        """The writer names each block from SpectralGasNames; the reader has
        to resolve them through the same one or the two can disagree."""
        for path in (READER, WRITER):
            self.assertIn("call SpectralGasNames", code(path),
                          "%s must use the shared per-slot names" % path)

    def test_the_reader_compares_the_header_name(self):
        src = code(READER)
        self.assertIn("blockname", src,
                      "the reader must read the gas name out of the block "
                      "header, not just test for the word TFP")
        self.assertRegex(
            src, r"trim\(adjustl\(sa_tags\(gas\)\)\) == trim\(blockname\)",
            "the header name must be compared against the per-slot tags")

    def test_the_reader_is_driven_by_the_file_not_by_the_gas_list(self):
        """Walking the expected gas list is what made it positional. The loop
        reads blocks until the headers stop."""
        src = code(READER)
        body = src[src.index("call SpectralGasNames"):]
        head = body[:body.index("do cls = JAN, DEC")]
        self.assertNotRegex(
            head, r"do\s+gas\s*=\s*firstGas,\s*lastGas\s*\n\s*if \(gas == h2o\) cycle\s*\n\s*if \(gas - firstGas",
            "the block loop must not be bounded by the expected gas list")

    def test_an_unwanted_block_is_consumed_not_skipped(self):
        """A block for a gas this project does not carry still occupies its
        lines; stepping over it without reading them desynchronises the file."""
        src = code(READER)
        self.assertIn("skipFn", src,
                      "a block with no matching slot must still be read")

    def test_shortness_is_decided_by_what_was_found(self):
        """It used to mean 'the loop stopped early', which a file-driven loop
        no longer expresses. It now means a wanted gas got no block."""
        src = code(READER)
        #: Classes, not months. Both ranges are 1..12 and both used to be
        #: spelled JAN:DEC, which is how a month index came to be stored where
        #: a class index belongs. The spelling is pinned so the two spaces
        #: cannot quietly merge again.
        self.assertRegex(src, r"all\(RegPar\(gas, 1:MaxGasClasses\)%fc == error\)")
        self.assertIn("short_file = .true.", src)


class SecondHygrometerRoundTrip(unittest.TestCase):
    """Writer, reader and readiness gate must agree about a second hygrometer.

    They did not. The writer skipped every water slot, so a project's second
    hygrometer was fitted on each run and thrown away; the reader skipped every
    water slot, so it could not have read one back; and the readiness gate
    demanded a fit for it before declaring the assessment usable. The gate was
    therefore asking for something the format had no way to carry.

    Three-sided, so a check on any one side alone would let the other two drift
    back out of agreement.
    """

    def test_the_writer_emits_a_block_for_every_non_primary_hygrometer(self):
        src = code(WRITER)
        self.assertIn("vapour TFP", src,
                      "a second hygrometer needs its own named RH table")
        self.assertRegex(
            src, r"if \(gas == wsl\) cycle",
            "only the PRIMARY is excluded from the named blocks - its table is "
            "the unnamed one at the fixed position")

    def test_the_reader_dispatches_on_numerosity(self):
        """A hygrometer's block is nine RH rows, a gas's is twelve months.

        The count column `numerosity` appears only on the RH tables and has
        since before the gas records, so the header says which shape follows.
        Reading twelve rows from a nine-row block would consume three lines of
        the next section and desynchronise everything after it - including the
        exponential coefficients that give water its cut-off.
        """
        src = code(READER)
        self.assertRegex(
            src, r"index\(dataline, 'numerosity'\) /= 0",
            "the block shape must be decided by the header, not assumed")

    def test_the_reader_excludes_only_the_primary(self):
        src = code(READER)
        #: The name-matching loop alone. The shortness loop further down does
        #: skip every water slot, and rightly: a hygrometer absent from the
        #: file is answered by the readiness gate and the analytic fallback,
        #: not by declaring the whole file short.
        body = src[src.index("blockname = adjustl"):src.index("'numerosity'")]
        self.assertNotIn(
            "GasSlotIsWater(gas)) cycle", body,
            "matching must exclude the primary alone, or a second hygrometer's "
            "named block can never resolve to a slot")
        self.assertIn("if (gas == wsl) cycle", body)

    def test_the_gate_looks_for_a_hygrometer_fit_in_the_rh_range(self):
        """The gate asked the month range of a slot binned by humidity."""
        src = code(GATE)
        body = src[src.index("assessment_ready = h2o_ok"):]
        arm = body[body.index("if (GasSlotIsWater(gas)) then"):]
        self.assertLess(
            arm.index("do cls = RH10, RH90"), arm.index("else"),
            "a hygrometer's fit lives in RH10..RH90, not 1..MaxGasClasses")


if __name__ == "__main__":
    unittest.main()
