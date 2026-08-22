"""Static checks for the flux detection limit (Wienhold et al. 1994).

Three things about this feature are easy to break silently and expensive to
notice, so they are pinned here.

**It must stay off unless asked for.** The whole point of adding it as an
option is that an existing project keeps its numbers. If the reader ever
stops defaulting `detlim_meth` to 'none', every project that never heard of
the key starts computing - and, worse, publishing - a column it did not ask
for.

**It must be computed before the time-lag shift.** The value is read off the
cross-covariance function, which only exists while the series are on their
raw alignment. `TimeLagHandle` shifts `Set` by each gas's lag near the end of
the routine; a call placed after that would measure the scatter of a function
that has already been collapsed, and would return plausible-looking numbers
that mean nothing. The call therefore lives inside `TimeLagHandle`, beside
the water covariances that have the same constraint.

**The ex record is positional.** `_DETLIM` is written by RP between the
verbatim LGD/KID/ZCD/CORRDIFF/NSR chunk and the Foken statistics, and parsed
back by `ReadExRecord` at exactly that point. If the writer and the reader
ever disagree about where it sits, nothing raises: every field after it
shifts by one column and the record is silently misread.
"""

import re
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / "src"

TAGS = (SRC / "src_rp" / "m_rp_global_var.f90").read_text(encoding="utf-8", errors="replace")
READER = (SRC / "src_rp" / "read_ini_rp.f90").read_text(encoding="utf-8", errors="replace")
DETLIM = (SRC / "src_rp" / "flux_detection_limit.f90").read_text(encoding="utf-8", errors="replace")
TLAG = (SRC / "src_rp" / "timelag_handle.f90").read_text(encoding="utf-8", errors="replace")
FLUXNET_HDR = (SRC / "src_rp" / "init_fluxnet_file_rp.f90").read_text(encoding="utf-8", errors="replace")
FLUXNET_OUT = (SRC / "src_rp" / "write_out_fluxnet.f90").read_text(encoding="utf-8", errors="replace")
EX = (SRC / "src_common" / "read_ex_record.f90").read_text(encoding="utf-8", errors="replace")
FCC_FLUXNET = (SRC / "src_fcc" / "write_out_fluxnet_fcc.f90").read_text(encoding="utf-8", errors="replace")
RP_FULL_HDR = (SRC / "src_rp" / "init_outfiles_rp.f90").read_text(encoding="utf-8", errors="replace")
FCC_FULL_HDR = (SRC / "src_fcc" / "init_out_files.f90").read_text(encoding="utf-8", errors="replace")
RP_FULL = (SRC / "src_rp" / "write_out_full.f90").read_text(encoding="utf-8", errors="replace")
FCC_FULL = (SRC / "src_fcc" / "write_out_full_fcc.f90").read_text(encoding="utf-8", errors="replace")

#: The three keys and the SNTags slots they were given, out of the blank run
#: at 55-68. All three are well below rpGasOriginN, so the per-gas records do
#: not move; test_ini_tag_collisions_static.py is what actually asserts that.
KEYS = {
    55: "detlim_meth",
    56: "detlim_offset_s",
    57: "detlim_window_s",
}

#: Wienhold's own window geometry. These are the values the feature carries so
#: that switching it on reproduces the published method rather than something
#: chosen here, and they are the values the GUI must offer as defaults.
WIENHOLD_OFFSET_S = "100d0"
WIENHOLD_WINDOW_S = "50d0"


class TheKeysReachTheEngine(unittest.TestCase):

    def test_each_key_holds_its_slot(self):
        for index, label in KEYS.items():
            self.assertRegex(
                TAGS,
                r"SNTags\(%d\)%%Label\s*/\s*'%s'\s*/" % (index, label),
                "SNTags(%d) is no longer %r" % (index, label),
            )

    def test_each_key_is_read_under_its_found_guard(self):
        for index in KEYS:
            self.assertIn(
                "SNTagFound(%d)" % index,
                READER,
                "SNTags(%d) is declared but never read - an absent key would "
                "then carry whatever the previous parse left in the saved "
                "array" % index,
            )


class TheFeatureIsOffUnlessAsked(unittest.TestCase):

    def test_the_method_defaults_to_none(self):
        self.assertRegex(
            READER,
            r"RPSetup%detlim_meth\s*=\s*'none'",
            "the literal default is gone; a project that never states "
            "detlim_meth would stop being inert",
        )

    def test_the_windows_default_to_wienholds_values(self):
        self.assertRegex(READER, r"RPSetup%detlim_offset_s\s*=\s*" + WIENHOLD_OFFSET_S)
        self.assertRegex(READER, r"RPSetup%detlim_window_s\s*=\s*" + WIENHOLD_WINDOW_S)

    def test_the_routine_returns_immediately_when_not_selected(self):
        self.assertRegex(
            DETLIM,
            r"if \(RPSetup%detlim_meth /= 'wienhold_94'\) return",
            "the early return is what makes the feature free when it is off",
        )

    def test_the_limit_is_error_filled_before_that_return(self):
        #> Order matters: the array has to be cleared before the early return,
        #> or a period that does not compute it publishes the previous
        #> period's numbers.
        clear = DETLIM.index("Essentials%detlim = error")
        bail = DETLIM.index("if (RPSetup%detlim_meth /= 'wienhold_94') return")
        self.assertLess(
            clear, bail,
            "Essentials%detlim must be cleared before the early return, or a "
            "run with the feature off carries stale values forward",
        )


class TheCallSiteSeesTheRawAlignment(unittest.TestCase):

    def test_it_is_called_from_timelag_handle(self):
        self.assertIn("call FluxDetectionLimit(", TLAG)

    def test_it_is_called_before_the_series_are_shifted(self):
        call = TLAG.index("call FluxDetectionLimit(")
        #> The shift is the block that fills TmpSet from Set at each gas's
        #> lag. Anchored on its comment rather than on a line number.
        shift = TLAG.index("Align data according to relevant time-lags")
        self.assertLess(
            call, shift,
            "FluxDetectionLimit runs after the time-lag shift, so the "
            "cross-covariance function it measures no longer exists",
        )

    def test_it_is_not_also_called_from_the_main_program(self):
        main = (SRC / "src_rp" / "eddyflow-rp_main.f90").read_text(
            encoding="utf-8", errors="replace")
        self.assertNotIn(
            "call FluxDetectionLimit(", main,
            "a second call from the main program would run against the "
            "shifted array",
        )


class TheRecordLayoutAgrees(unittest.TestCase):

    def test_rp_writes_the_family_and_fcc_reads_it_back(self):
        self.assertIn("'_DETLIM'", FLUXNET_HDR)
        self.assertIn("Essentials%detlim(FluxnetLayoutSlots(j))", FLUXNET_OUT)
        self.assertIn("lEx%detlim(exSlots(jx))", EX)
        self.assertIn("lEx%detlim(gas)", FCC_FLUXNET)

    def test_the_reader_counts_the_family_by_a_named_parameter(self):
        #> read_ex_record navigates by counted commas and forbids bare
        #> literals for exactly this reason.
        self.assertRegex(
            EX,
            r"nDetlimFields\s*=\s*n_layout_gas",
            "the detlim field count is not expressed as a named parameter "
            "over the gas count",
        )

    def test_the_family_is_parsed_between_the_chunk_and_the_foken_stats(self):
        chunk = EX.index("fluxnetChunks%s(2) = dataline")
        detlim = EX.index("lEx%detlim(exSlots(jx))")
        foken = EX.index("lEx%TAU_SS, lEx%H_SS")
        self.assertLess(chunk, detlim, "detlim is parsed before the chunk it follows")
        self.assertLess(detlim, foken, "detlim is parsed after the Foken statistics")

    def test_the_writer_puts_it_in_the_same_place(self):
        nsr = FLUXNET_OUT.index("Essentials%mahrt98_NR(FluxnetLayoutSlots(j))")
        detlim = FLUXNET_OUT.index("Essentials%detlim(FluxnetLayoutSlots(j))")
        ss = FLUXNET_OUT.index("STDiff%w_u")
        self.assertLess(nsr, detlim)
        self.assertLess(
            detlim, ss,
            "the writer emits detlim somewhere other than between the NSR "
            "family and the Foken statistics, where the reader looks for it",
        )


class TheFullOutputMirrorsItself(unittest.TestCase):
    """RP and FCC both write the full output, and must agree on its shape."""

    def test_both_headers_name_the_column(self):
        for name, text in (("src_rp", RP_FULL_HDR), ("src_fcc", FCC_FULL_HDR)):
            self.assertIn(
                "'detlim'", text,
                "%s does not name the detlim column in the full output header"
                % name,
            )

    def test_both_headers_widened_the_group_by_one(self):
        #> The per-gas concentration group grew from five members to six, so
        #> header1's comma run has to grow with it or every following group
        #> label lands one column early.
        for name, text in (("src_rp", RP_FULL_HDR), ("src_fcc", FCC_FULL_HDR)):
            self.assertIn(
                "',,,,,'", text,
                "%s still pads the concentration group for five members, not "
                "six" % name,
            )

    def test_both_writers_emit_it(self):
        self.assertIn("Essentials%detlim(gas)", RP_FULL)
        self.assertIn("lEx%detlim(gas)", FCC_FULL)

    def test_both_headers_label_the_units_as_a_covariance(self):
        #> Not a flux. The limit qualifies the covariance and is never scaled
        #> with it, so labelling it in flux units would invite a comparison
        #> that does not hold.
        for name, text in (("src_rp", RP_FULL_HDR), ("src_fcc", FCC_FULL_HDR)):
            self.assertIn(
                "[1=default],[cov]", text,
                "%s does not label the detlim column in covariance units"
                % name,
            )


class TheEstimatorIsTheOneWienholdDescribes(unittest.TestCase):

    def test_it_averages_both_sides_of_the_peak(self):
        self.assertRegex(
            DETLIM,
            r"do side = -1, 1, 2",
            "the two windows either side of the lag are what Wienhold "
            "averages; one side alone is a different estimator",
        )

    def test_a_window_needs_more_than_two_lags(self):
        self.assertRegex(
            DETLIM,
            r"min_lags\s*=\s*3",
            "the standard deviation of fewer than three covariances says "
            "nothing about the spread of the function",
        )

    def test_the_scatter_is_a_sample_standard_deviation(self):
        self.assertIn("dble(n - 1)", DETLIM)


class TheCovarianceHelperTakesEitherSign(unittest.TestCase):
    """CovarianceW had to learn negative lags for this to work at all."""

    def test_the_loop_bounds_admit_a_negative_lag(self):
        #> A gas at 16 s sampled 100 s earlier sits at -84 s. The old bounds
        #> `do i = 1, nrow - lag` with col2(i+lag) read off both ends there.
        self.assertRegex(
            TLAG,
            r"do i = max\(1, 1 - lag\), min\(nrow, nrow - lag\)",
            "CovarianceW no longer guards both ends for a negative lag",
        )

    def test_a_positive_lag_still_walks_the_same_terms(self):
        #> For lag >= 0 the bounds collapse to 1 and nrow-lag, which is what
        #> the routine always did - this is what keeps the existing callers
        #> byte-identical.
        for lag, nrow in ((0, 100), (5, 100), (37, 18000)):
            self.assertEqual(max(1, 1 - lag), 1)
            self.assertEqual(min(nrow, nrow - lag), nrow - lag)


if __name__ == "__main__":
    unittest.main()
