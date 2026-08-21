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

    def test_the_reader_excludes_only_the_slot_the_unnamed_table_filled(self):
        src = code(READER)
        #: The name-matching loop alone. The shortness loop further down does
        #: skip every water slot, and rightly: a hygrometer absent from the
        #: file is answered by the readiness gate and the analytic fallback,
        #: not by declaring the whole file short.
        body = src[src.index("blockname = adjustl"):src.index("'numerosity'")]
        self.assertNotIn(
            "GasSlotIsWater(gas)) cycle", body,
            "matching must exclude one hygrometer alone, or a second "
            "hygrometer's named block can never resolve to a slot")
        #: `water_slot`, not `wsl`. The unnamed table at the top of the file
        #: goes to the hygrometer its stamp names, which is this project's
        #: primary only while the file agrees about which one that is.
        #: Excluding `wsl` instead would leave the stamped slot open to a
        #: second assignment and lock the real primary out of its own block.
        self.assertIn("if (gas == water_slot) cycle", body)

    def test_a_block_resolves_by_what_it_is_for_before_what_it_is_called(self):
        """The block name is an ordinal over repeats of a species, so it says
        nothing about which analyser a block belongs to. Two CO2 records and a
        re-ordered project is enough to hand one analyser's transfer function
        to the other - and on CH-LAE the two cells are 943 hPa and 70 hPa."""
        src = code(READER)
        body = src[src.index("blockname = adjustl"):src.index("'numerosity'")]
        self.assertIn("slot = SlotFromSpectralStamp(dataline)", body,
                      "the stamp the writer puts on every block must be tried")
        self.assertLess(
            body.index("SlotFromSpectralStamp"),
            body.index("trim(adjustl(sa_tags(gas))) == trim(blockname)"),
            "the stamp has to win over the name, or reordering still misassigns")
        self.assertIn("if (slot == 0) then", body,
                      "and the name match has to remain the fallback, since "
                      "every file written before the stamp relies on it")

    def test_old_files_naming_the_second_hygrometer_still_resolve(self):
        """`<TAG> vapour TFP` was what the writer emitted, and the reader takes
        everything before TFP as the name - so the block called itself
        `H2O_2 VAPOUR`, matched nothing, and was silently discarded."""
        src = code(READER)
        self.assertIn("' VAPOUR'", src,
                      "a file already written with the ' vapour' header must "
                      "still resolve, or the fix strands the files that "
                      "needed it")

    def test_every_hygrometer_gets_its_own_rh_relation(self):
        """The nine RH-class cut-offs are not what the iir correction uses.

        It evaluates exp(A*RH^2 + B*RH + C), and A/B/C used to be one
        project-wide set fitted from the primary alone. So a second hygrometer
        was fitted, written, read back - and then given the primary's curve,
        differing only in the humidity it was evaluated at. On CH-LAE the two
        read some fourteen points apart, more than an RH class.
        """
        fit = code("src/src_fcc/fit_rh_to_cutoff.f90")
        self.assertNotIn("RegPar(dum, dum)%e1 = EXPPar(1)", fit,
                         "the fit must land on the hygrometer it fitted")
        self.assertIn("RegPar(wsl, dum)%e1 = EXPPar(1)", fit)
        self.assertIn("if (.not. GasSlotIsWater(wsl)) cycle", fit,
                      "the fit has to run per water slot, not once")
        #> The primary's copy still has to reach the shared slot: it is what
        #> the file's standalone section carries and what an unfitted
        #> hygrometer falls back to.
        self.assertIn("RegPar(dum, dum)%e1 = RegPar(primary, dum)%e1", fit)
        self.assertIn("FCCMetadata%GasPathType(wsl)", fit,
                      "which arm a fit takes is a property of the analyser "
                      "being fitted, not of the primary")

        aux = code("src/src_common/bpcf_aux_subs.f90")
        arm = aux[aux.index("case('iir')"):aux.index("case('sigma')")]
        self.assertIn("A = RegPar(gas, dum)%e1", arm,
                      "the iir arm must take each hygrometer's own coefficients")
        self.assertIn("A = RegPar(dum, dum)%e1", arm,
                      "with the project-wide set as the fallback for a "
                      "hygrometer the assessment never fitted")

    def test_the_rh_relation_survives_the_file(self):
        """Fitting it per hygrometer is no use if the file cannot carry it."""
        self.assertIn("exp=", code(WRITER),
                      "a hygrometer block must state its own coefficients")
        self.assertIn("'exp='", code(READER),
                      "and the reader must take them back off it")

    def test_the_writer_stamps_every_block(self):
        src = code(WRITER)
        self.assertEqual(
            3, src.count("call SpectralBlockStamp("),
            "all three block kinds are stamped: the unnamed primary water "
            "table, each gas, and each further hygrometer")
        self.assertNotIn(
            "' vapour TFP", src,
            "the hygrometer block name must be the bare tag - the reader "
            "slices at TFP, so any word before it becomes part of the name")

    def test_the_gate_looks_for_a_hygrometer_fit_in_the_rh_range(self):
        """The gate asked the month range of a slot binned by humidity.

        The test now reads GasHasSpectralFit, which is where that question
        moved when the writer and the report stopped each answering it their
        own way. The requirement is unchanged.
        """
        src = code("src/src_common/gas_slot_resolution.f90")
        body = src[src.index("logical function GasHasSpectralFit"):]
        body = body[:body.index("end function GasHasSpectralFit")]
        arm = body[body.index("if (GasSlotIsWater(gas_slot)) then"):]
        self.assertLess(
            arm.index("do cls = RH10, RH90"), arm.index("else"),
            "a hygrometer's fit lives in RH10..RH90, not 1..MaxGasClasses")

    def test_the_file_is_written_for_any_fitted_gas_not_only_for_water(self):
        """Water failing used to discard the blocks of gases that had been
        fitted perfectly well, and every gas then fell back to analytic."""
        src = code(WRITER)
        self.assertIn("if (n_fitted == 0) then", src,
                      "the assessment file is gated on some gas being fitted")
        gate = src[src.index("if (FCCsetup%do_spectral_assessment) then"):]
        gate = gate[:gate.index("call LogSay(' Writing spectral assessment")]
        self.assertNotIn(
            "goodj == ierror", gate,
            "the file is gated on water's RH classes again; goodj is water's "
            "first usable class and belongs to the H2O spectra file below")
        self.assertIn("call ExceptionHandler(110)", src,
                      "a partial file must announce itself on the console")

    def test_the_partial_warning_exists_and_names_both_outcomes(self):
        src = code("src/src_common/exception_handler.f90")
        self.assertIn("case(110)", src)
        block = src[src.index("case(110)"):]
        block = block[:block.index("end select")]
        for phrase in ("analytical", "readiness report"):
            self.assertIn(phrase, block,
                          "Warning(110) must say what the unfitted gases got "
                          "and where to read which gas is which")

    def test_the_method_is_chosen_per_gas_and_does_not_stick(self):
        """One gas without a fit used to demote every gas, and the demotion
        was written back into project state so it lasted the whole run."""
        src = code("src/src_common/bpcf_bandpass_spectral_corrections.f90")
        self.assertNotIn(
            "EddyFlowProj%hf_meth = actual_hf_method", src,
            "the per-period fallback is written back into project state "
            "again, which makes one bad period demote every later one")
        self.assertIn("insitu_ok", src, "the in-situ mask is gone")
        self.assertIn("call BPCF_Moncrieff97", src)
        body = src[src.index("select case(trim(adjustl(actual_hf_method)))"):]
        self.assertIn(
            "analytic_only", body,
            "the gases the in-situ method could not take must still be "
            "corrected, by the analytic method, after the in-situ pass")


if __name__ == "__main__":
    unittest.main()
