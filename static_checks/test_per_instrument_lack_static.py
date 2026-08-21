"""Guards the per-instrument missing-samples allowance.

The completeness test measured every column against the file's row grid. An
instrument slower than that grid cannot fill it - a 1 Hz analyser in a 20 Hz
file writes one row in twenty - so its column read as 95 % missing and was
dropped outright, with no way to say "this instrument is slower".

Three things have to stay true for that to keep working, and none of them is
visible in a single-rate regression run, because every fixture has one rate and
the new arithmetic reduces exactly to the old one there:

  1. the drop test is measured per column, not against MaxPeriodNumRecords;
  2. the whole-period gates stay global - they count rows in which ANY column
     is valid, so a sparse column contributes nothing to them and a per-column
     allowance would say nothing about them;
  3. the allowance is read under SNTagFound, so an absent key falls back to the
     project-wide setting instead of to zero. A max_lack of 0 drops every
     column that is missing a single sample.
"""

import re
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[1]


def read(path):
    return (ROOT / path).read_text(encoding="utf-8", errors="replace")


class PerInstrumentLackStaticTests(unittest.TestCase):
    def test_the_drop_test_is_per_column(self):
        body = read("src/src_rp/eliminate_corrupted_variables.f90")
        assert "ColumnMaxLack(i)" in body, (
            "the allowance is back to a single project-wide number"
        )
        assert "ColumnAcFreq(i)" in body, (
            "the expected sample count is back on the file's rate"
        )
        assert "MaxPeriodNumRecords * RPsetup%max_lack" not in body, (
            "the row-grid count is back as the denominator; it is what read a "
            "complete 1 Hz column as 95 % missing"
        )

    def test_the_drop_test_counts_what_is_present(self):
        """`expected - count(present)` and `count(error)` differ the moment a
        column is slower than the rows: the structural gaps between its samples
        are error rows that were never going to hold data."""
        body = read("src/src_rp/eliminate_corrupted_variables.f90")
        assert "count(LocSet(:, i) /= error)" in body, body[body.index("do i = 1"):]

    def test_the_column_expectation_reduces_to_the_period_length(self):
        """A column at the file's own rate must give exactly the old numbers, or
        every existing project moves. The denominator has to be built from the
        same expression MaxPeriodNumRecords is."""
        body = read("src/src_rp/eliminate_corrupted_variables.f90")
        main = read("src/src_rp/eddyflow-rp_main.f90")
        assert "nint(RPsetup%avrg_len * 60d0 * freq)" in body
        assert "nint(RPsetup%avrg_len     * 60d0 * Metadata%ac_freq)" in main, (
            "MaxPeriodNumRecords is no longer built the way the per-column "
            "expectation assumes; the reduction to the old test is broken"
        )

    def test_the_spectra_feasibility_test_is_per_column(self):
        """A third of what the column should have produced, not a third of the
        rows. Against the rows, a 1 Hz column in a 10 Hz file is nine tenths
        error and was excluded from the spectra whatever its data said."""
        body = read("src/src_rp/fix_dataset_for_spectra.f90")
        assert "ColumnAcFreq(j)" in body
        assert "> nrow / 3" not in body, (
            "the arbitrary third is back to a third of the rows"
        )
        assert "expected / 3" in body

    def test_the_gap_duration_discounts_the_sampling_interval(self):
        """A slow column has F/f - 1 error rows between consecutive samples by
        construction. Counted as a gap, that is a permanent 0.9 s outage that
        never happened."""
        body = read("src/src_rp/longest_gap_duration.f90")
        assert "nint(Metadata%ac_freq / ColumnAcFreq(icol)) - 1" in body
        assert "max(0, LongestVariableGap" in body

    def test_a_column_reports_nothing_above_its_own_nyquist(self):
        """FixDatasetForSpectra interpolates a slow column's missing rows up
        onto the fast grid, so every bin to the station's Nyquist carries a
        number and the ones above the column's own are interpolation
        artefacts. Measured on base_slow: the MIRO's CO2 at 1 Hz is blanked
        from 0.515 Hz up while the LI-7200's continues to 4.58 Hz."""
        body = read("src/src_rp/spectral_analysis.f90")
        assert "subroutine CapSpectraAtColumnNyquist" in body
        #> Both axes: the unbinned one that goes to the full-cospectra file and
        #> the binned one. Blanking only one leaves the other fabricating.
        assert body.count("call CapSpectraAtColumnNyquist(") == 2, (
            "the cap is applied to only one of the two frequency axes"
        )
        cap = body[body.index("subroutine CapSpectraAtColumnNyquist"):]
        cap = cap[:cap.index("end subroutine CapSpectraAtColumnNyquist")]
        assert "if (ColumnAcFreq(j) >= Metadata%ac_freq) cycle" in cap, (
            "a column at the file's own rate must not be touched, or every "
            "existing project moves"
        )
        assert "= error" in cap and "= 0d0" not in cap, (
            "blank to the error code; a spectral density of zero is a claim "
            "about the data, not the absence of one"
        )

    def test_a_slower_column_is_rebuilt_on_its_own_grid(self):
        """Below its Nyquist the full-rate pass gives a slow column the
        spectrum of its INTERPOLATION - the linear fill between real samples
        is a low-pass filter, and it showed: on base_slow the MIRO's CO2 read
        0.242 at 0.43 Hz interpolated against 0.623 rebuilt."""
        body = read("src/src_rp/spectral_analysis.f90")
        assert "subroutine SlowColumnSpectra" in body
        #> Both passes: the unbinned one that feeds the full-cospectra file and
        #> the degraded covariances, and the binned one.
        assert body.count("call SlowColumnSpectra(") == 2, (
            "the rebuild runs on only one of the two passes"
        )
        #> From the untouched samples. Tapering and FourierTransform both work
        #> in place, so the binned pass transforms a copy and leaves the
        #> argument alone - which is also what lets the dummy stay intent(in).
        assert "real(kind = dbl) :: WorkSet(N, M)" in body
        assert "WorkSet = Set" in body
        assert "call SlowColumnSpectra(Set," in body

    def test_the_spectral_input_is_not_written_through(self):
        """SpectralAnalysis declares Set intent(in) and used to hand it
        straight to Tapering and FourierTransform, both intent(inout). It
        compiled only because those two have implicit interfaces, and it left
        the caller's array transformed in place."""
        body = read("src/src_rp/spectral_analysis.f90")
        fn = body[body.index("subroutine SpectralAnalysis"):]
        fn = fn[:fn.index("end subroutine SpectralAnalysis")]
        assert "real(kind = dbl), intent(in) :: Set(N, M)" in fn
        for call in ("call Tapering(RPsetup%tap_win, Set,",
                     "call FourierTransform(Set,"):
            assert call not in fn, (
                f"{call.strip()} writes through an intent(in) dummy"
            )

        fn = body[body.index("subroutine SlowColumnSpectra"):]
        fn = fn[:fn.index("end subroutine SlowColumnSpectra")]
        assert "if (freq >= Metadata%ac_freq) cycle" in fn, (
            "a column at the file's own rate must be left alone, or every "
            "existing project moves"
        )
        assert "modulo(SpecPhase(j) - SpecRowOffset, stride)" in fn, (
            "sampling at the wrong phase reads interpolated blends of two "
            "real samples and quietly attenuates what this exists to measure"
        )
        #> The normalisation has to come from the decimated series: the
        #> full-rate variance belongs to a different signal.
        assert "dspec / var_gas" in fn and "dcosp / cov_wgas" in fn

    def test_the_rebuild_runs_before_the_ogives(self):
        """An ogive is the integral of the spectrum beside it. Rebuilt after,
        the two files disagree about the same column - which they did, for
        every one of the slow fixture's six periods."""
        body = read("src/src_rp/spectral_analysis.f90")
        assert (body.index("call SlowColumnSpectra(Set, N, M, 'squared'")
                < body.index("call AllOgives(")), (
            "the ogives integrate the interpolated spectrum again"
        )

    def test_the_decimation_fabricates_no_sample(self):
        """The series is counted from the phase so the last sample lands on or
        before row N. Clamping instead repeated the last row - a made-up
        sample, in the one routine whose purpose is to use only real ones."""
        fn = read("src/src_rp/spectral_analysis.f90")
        fn = fn[fn.index("subroutine SlowColumnSpectra"):]
        fn = fn[:fn.index("end subroutine SlowColumnSpectra")]
        assert "nd = (N - phase) / stride" in fn
        assert "if (row > N) row = N" not in fn, (
            "the clamp is back, and with it a duplicated sample"
        )

    def test_the_sample_offset_is_the_commonest_not_the_first(self):
        """A column whose first sample is missing - an ordinary gap at the
        start of a period - reported phase zero, and then every rebuilt sample
        read an interpolated blend rather than a measurement."""
        body = read("src/src_rp/fix_dataset_for_spectra.f90")
        assert "PhaseIntervals" in body, (
            "the phase is taken from a single interval again"
        )
        assert "maxloc(tally(0:stride - 1), dim = 1) - 1" in body, (
            "the offset must be the commonest over many intervals"
        )

    def test_the_rebuild_is_normalised_once(self):
        """NormalizeCoSpectra divides by the full-rate statistics. A column
        rebuilt before it would be divided twice, by two different numbers."""
        body = read("src/src_rp/spectral_analysis.f90")
        binned = body[body.index("if (RPsetup%out_bin_sp .or. RPsetup%out_bin_og) then"):]
        assert (binned.index("call NormalizeCoSpectra(")
                < binned.index("call SlowColumnSpectra(")), (
            "the rebuild runs before the normalisation"
        )

    def test_w_is_paired_by_what_the_instrument_did(self):
        """A point-sampled gas against an averaged w biases the covariance, so
        instantaneous is the default and the average is opt-in. Measured on
        base_slow against base_slow_integr: 32 cospectral bins move and not one
        spectral bin does."""
        fn = read("src/src_rp/spectral_analysis.f90")
        fn = fn[fn.index("subroutine SlowColumnSpectra"):]
        fn = fn[:fn.index("end subroutine SlowColumnSpectra")]
        assert "if (E2Col(j)%instr%integrates) then" in fn, (
            "instr_<K>_integrates is read by nothing again"
        )
        assert "sum(Set(lo:row, w)) / dble(row - lo + 1)" in fn, (
            "the integrating branch must average w over the gas's own "
            "sampling interval"
        )
        assert "raw_w(i) = Set(row, w)" in fn, (
            "the default branch must point-sample w at the gas's instant"
        )

    def test_the_whole_period_gates_stay_global(self):
        """Essentials%n_in and friends count rows in which any column is valid,
        so these are the file's own completeness, not any instrument's."""
        main = read("src/src_rp/eddyflow-rp_main.f90")
        gates = re.findall(r"MissingRecords > RPsetup%max_lack", main)
        assert len(gates) >= 3, (
            f"only {len(gates)} whole-period gates use the global allowance; "
            "one has been made per-instrument, which it cannot be"
        )

    def test_the_resolvers_fall_back_rather_than_default_to_zero(self):
        body = read("src/src_rp/column_sampling.f90")
        assert "ColumnAcFreq = Metadata%ac_freq" in body
        assert "ColumnMaxLack = RPsetup%max_lack" in body
        assert "if (E2Col(icol)%instr%ac_freq > 0d0)" in body, (
            "an unset ac_freq is the error sentinel, not a rate"
        )
        assert "min(E2Col(icol)%instr%ac_freq, Metadata%ac_freq)" in body, (
            "a rate declared above the file's would expect more samples than "
            "the period has rows, and drop every column of that instrument"
        )
        assert "if (RPsetup%instr_max_lack(slot) >= 0d0)" in body, (
            "an unset allowance is the error sentinel; using it would drop "
            "every column"
        )

    def test_the_allowance_is_read_only_where_the_project_states_it(self):
        reader = read("src/src_rp/read_ini_rp.f90")
        assert "RPsetup%instr_max_lack = error" in reader, (
            "the allowance must start at the sentinel, or an instrument the "
            "project says nothing about inherits the previous run's value"
        )
        assert "SNTagFound(rpInstrMaxLackN + i - 1)" in reader, (
            "an absent key must fall back to the global, not to zero"
        )

    def test_a_column_can_name_its_instrument(self):
        """A column stores a COPY of its instrument, so without the slot number
        the project's per-instrument key cannot be resolved back to a column."""
        typedef = read("src/src_common/m_typedef.f90")
        instr = typedef[typedef.index("type :: InstrumentType"):]
        instr = instr[:instr.index("end type InstrumentType")]
        assert "integer :: slot" in instr

        meta = read("src/src_common/read_metadata_file.f90")
        assert "Instr(i)%slot = i" in meta, (
            "nothing fills the slot number, so every column resolves to 0 and "
            "the per-instrument allowance is unreachable"
        )


if __name__ == "__main__":
    unittest.main()
