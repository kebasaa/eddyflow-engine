"""Where a PWB time lag is allowed to come from.

A time lag is a property of one sampling tube. Handing one analyser's lag to
another analyser's gases is not an approximation of the right answer, it is a
different measurement -- so the rules about where a lag may travel from are
worth stating separately from the rest of the detector, and worth pinning.

Two rules, and this file exists to keep them true:

  * never across instruments, and never on the strength of a model string,
    which cannot tell two LI-7200s at one site apart;
  * the gas's OWN lag first, in all three of its forms, before any borrowed
    one -- two gases down one tube still have measurably different delays, so
    borrowing trades a stale number for a biased one.

The second is dyco's precedence and its reasoning; see dyco/pwb.py's
fill_tlag_gaps.

Part of the EddyFlow engine's static checks.
"""

from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]

MODULE = "src/src_rp/pwb_timelag_handle.f90"
HANDLER = "src/src_rp/timelag_handle.f90"
CORE = "src/src_common/m_pwb_core.f90"


def read(path):
    return (ROOT / path).read_text(encoding="utf-8")


def code(path):
    """Source with full-line comments removed.

    Every rule below is explained in a comment that names the thing it
    forbids, so a plain scan of the raw text matches the explanation.
    """
    return (chr(10)).join(ln for ln in read(path).splitlines()
                          if not ln.lstrip().startswith("!"))


def body_of(source, opener, closer):
    return source[source.index(opener):source.index(closer)]


class ALagIsNeverBorrowedAcrossInstruments(unittest.TestCase):

    def test_the_identity_test_demands_a_named_instrument_on_both_sides(self):
        #> The model string is not an identity. A record naming no instrument
        #> cannot be shown to share a tube with anything, so it neither
        #> donates nor receives: absence of evidence is not identity.
        block = body_of(code(MODULE), "logical function SameAnalyser",
                        "end function SameAnalyser")
        self.assertIn("SameAnalyser = .false.", block)
        self.assertIn("if (len_trim(E2Col(a)%instr_name) == 0) return", block)
        self.assertIn("if (len_trim(E2Col(b)%instr_name) == 0) return", block)
        self.assertIn("if (E2Col(a)%instr_name /= E2Col(b)%instr_name) return",
                      block)

    def test_the_model_string_can_only_refuse_never_admit(self):
        #> It is compared, but the comparison only ever returns early: equal
        #> names with unequal models is a malformed project, not a shared
        #> analyser. There must be no branch that ACCEPTS on the model.
        block = body_of(code(MODULE), "logical function SameAnalyser",
                        "end function SameAnalyser")
        self.assertIn("if (E2Col(a)%instr%model /= E2Col(b)%instr%model) return",
                      block)
        self.assertNotIn("else", block)
        self.assertNotIn("SameAnalyser = E2Col(a)%instr%model", block)

    def test_every_place_that_could_share_a_lag_asks_the_same_question(self):
        #> Three of them: the settled per-period table, the live classifier,
        #> and the aggregate summary -- which had no instrument test of any
        #> kind and so could hand one tube's optimised window to another
        #> tube's gas, then feed it back into the next run.
        module = code(MODULE)
        handler = code(HANDLER)
        self.assertIn("if (.not. SameAnalyser(gas, g)) cycle", module)
        self.assertIn("if (.not. SameAnalyser(gas, candidate)) cycle", module)
        self.assertIn("SameAnalyser(j, k)", handler)

    def test_no_donor_is_chosen_by_comparing_models(self):
        for path in (MODULE, HANDLER):
            source = code(path)
            self.assertNotIn("%instr%model /= E2Col(j)%instr%model", source)
            self.assertNotIn("%instr%model == E2Col(j)%instr%model", source)

    def test_water_never_donates_on_any_path(self):
        #> Water's lag depends on humidity in a way the trace gases' does not,
        #> so borrowing from it is worse than not borrowing.
        #>
        #> Stated as a pairing rather than a count: every donor search - the
        #> two per-period ones and the aggregate summary - refuses water, so
        #> the two totals move together and a fourth path added later has to
        #> refuse it as well to keep this passing. Counting to a literal meant
        #> adding a path failed here for no reason.
        module = code(MODULE)
        searches = module.count("SameAnalyser(gas, ")
        self.assertGreaterEqual(searches, 3, "the donor searches moved")
        self.assertEqual(module.count("if (GasSlotIsWater("), searches,
                         "a donor search does not refuse water")
        self.assertIn("GasSlotIsWater(k)", code(HANDLER))

    def test_a_rejected_lag_is_the_last_resort_not_the_second(self):
        """A gas the rule rejected everywhere must try its tube-mate first.

        Carbonyl sulfide is the case: its HDI routinely spans the whole search
        window, so every detection is prefiltered, there is nothing to take a
        median of, and before this the only thing left was its own rejected
        covariance maximisation. On the run this was found in, that put COS at
        10.6 s in a tube whose delay is 16.2 s, while the CO2 beside it in the
        same tube carried an interpolated 16.5 s.

        So there is a second shared pass, and it sits between the median and
        the terminal arm. Order is the whole point: before the median it would
        outrank the gas's own evidence, after the terminal arm it would never
        run.
        """
        module = code(MODULE)
        median = module.index("reliability_class = 'S3_median'")
        filled = module.index("reliability_class = 'S4_instrument_filled'")
        terminal = module.index("fill_method = 'maxcov_default'")
        self.assertLess(median, filled,
                        "the tube-mate now outranks the gas's own median")
        self.assertLess(filled, terminal,
                        "the terminal arm runs before the tube-mate is tried")

    def test_the_second_pass_takes_a_filled_donor_but_not_a_borrowed_one(self):
        #> The point of the pass is to accept a donor the rule did not trust
        #> outright - interpolated, carried, backfilled. What it must not
        #> accept is a lag this same pass just borrowed, or the value would
        #> walk from tube-mate to tube-mate with nothing behind it.
        module = code(MODULE)
        block = module[module.index("reliability_class) == 'fallback') cycle"):]
        block = block[:block.index("reliability_class = 'S4_instrument_filled'")]
        self.assertIn("== 'S4_instrument_filled') cycle", block)

    def test_the_terminal_arm_reads_the_streaming_answer_by_row(self):
        #> fallback_lag used to be filled per gas, indexed by that gas's
        #> ordinal, and read in the same loop. Three passes now sit between
        #> the capture and the read, so it is captured by cache row instead -
        #> indexing it by the old ordinal would hand a row another row's lag.
        module = code(MODULE)
        self.assertIn("fallback_lag(i) = PwbTimelagCache(i)%used_lag", module)
        self.assertIn("PwbTimelagCache(i)%used_lag = fallback_lag(i)", module)
        self.assertNotIn("fallback_lag(n) = ", module)
        self.assertNotIn("= fallback_lag(j)", module)


class AStaleFileCannotReintroduceIt(unittest.TestCase):
    """The aggregate file is read back by a later run.

    One written before the borrowing rule was restricted can name a donor on
    another analyser. The reader used to ignore those provenance lines
    entirely and take the window regardless.
    """

    def test_a_foreign_window_is_refused_not_taken(self):
        source = code("src/src_rp/read_timelag_opt_file.f90")
        self.assertIn("PWB_summary_source_for_", source)
        self.assertIn("inferred_from_", source)
        self.assertIn("if (SameAnalyser(gas, donor)) exit", source)
        self.assertIn("foreign(gas) = .true.", source)
        #> The flag has to suppress the three values that follow it.
        self.assertIn("if (foreign(gas)) cycle", source)


class TheGasOwnLagComesFirst(unittest.TestCase):

    def setUp(self):
        self.block = body_of(code(MODULE),
                             "subroutine PostProcessPwbTimelagCache",
                             "end subroutine PostProcessPwbTimelagCache")

    def test_the_donor_search_runs_after_every_own_lag_fill(self):
        order = [self.block.index("'interpolated'"),
                 self.block.index("'carryforward'"),
                 self.block.index("'backfilled'"),
                 self.block.index("'instrument_shared'"),
                 self.block.index("'median'")]
        self.assertEqual(
            order, sorted(order),
            "all three forms of the gas's own lag must be exhausted before a "
            "donor is consulted")

    def test_the_carry_limit_bounds_all_three_directions(self):
        #> Bounding one would achieve nothing: with detections either side of
        #> a long unusable stretch, an unbounded backward fill covers exactly
        #> the span the forward carry was forbidden to cross.
        self.assertEqual(
            self.block.count("if (limit > 0d0 .and. dist > limit) cycle"), 3)
        self.assertIn("limit = PWBSetup%max_carry_h * 60d0", self.block)
        self.assertIn("'S3_expired'", self.block)

    def test_the_limit_is_elapsed_time_not_a_row_count(self):
        #> The table has a row only where a period was processed, so counting
        #> rows would reach straight across a gap in the raw files.
        self.assertIn("PeriodMinutes", code(MODULE))
        self.assertIn("carry_hours", self.block)

    def test_an_evidence_based_fill_clears_the_fallback_flag(self):
        #> fallback_used is raised at detection for a period that produced no
        #> usable lag, so an arm that later finds one has to lower it or the
        #> run summary reports a settled period as a fallback.
        self.assertGreaterEqual(
            self.block.count("%fallback_used = .false."), 4)


class TheSettingReachesTheEngine(unittest.TestCase):

    def test_the_carry_limit_is_read_and_fingerprinted(self):
        self.assertIn("PWBSetup%max_carry_h = 24d0", code("src/src_rp/read_ini_rp.f90"))
        self.assertIn("SNTagFound(422)", code("src/src_rp/read_ini_rp.f90"))
        self.assertIn("SNTags(422)%Label / 'pwb_max_carry_h' /",
                      read("src/src_rp/m_rp_global_var.f90"))
        #> In the fingerprint, or a table computed under a different limit is
        #> reused as though it were current.
        self.assertIn("_carry=", code(MODULE))


class TheNumbersFollowTheReference(unittest.TestCase):
    """The three places this engine had disagreed with RFlux.

    All three WERE pinned numerically, by the RFlux comparison that used to
    live in test_pwb_reference_static.py; that went with dyco, whose fixtures
    it needed. These assertions are what is left: they say where in the source
    each agreement comes from, so a later edit that reintroduces one fails
    here with a reason. Nothing now fails with a number.
    """

    def test_the_covariance_is_read_off_the_undifferenced_series(self):
        core = code(CORE)
        self.assertIn("ss_undiff = ss", core)
        self.assertIn("ww_undiff = ww", core)
        self.assertIn("call LinearDetrend(ww_undiff, n)", core)
        self.assertIn("call ComputeCcovWindow(ww_undiff, ss_undiff", core)

    def test_the_covariance_divides_by_n_not_by_the_overlap(self):
        #> Dividing by the overlap inflates it by n/(n-|lag|), which at lag
        #> 169 of 6000 is 2.9%.
        block = body_of(code(CORE), "subroutine ComputeCcovWindow",
                        "end subroutine ComputeCcovWindow")
        self.assertIn("ccov(lag) = cov / dble(n)", block)
        self.assertNotIn("/ dble(nn)", block)

    def test_differencing_returns_one_sample_fewer(self):
        block = body_of(code(CORE), "subroutine DifferenceSeries",
                        "end subroutine DifferenceSeries")
        self.assertIn("ne = n - 1", block)

    def test_the_full_data_ccf_drops_the_ar_initialisation_samples(self):
        #> R drops them with na.action = na.omit. Keeping them as zeros moved
        #> cor_pww in the fourth significant digit at AR order 67. The
        #> BOOTSTRAP keeps them, as R does.
        core = code(CORE)
        self.assertIn("w_fs(p_s+1:ne), s_fs(p_s+1:ne)", core)


class TheBlockLengthIsDerived(unittest.TestCase):

    def test_it_is_floored_at_twice_the_widest_bound(self):
        #> A block shorter than the lag range cannot contain the structure the
        #> bootstrap exists to preserve. This used to warn and then use the
        #> short block anyway.
        source = code(MODULE)
        self.assertIn("block_len = max(requested_block_len, 2 * widest)", source)
        self.assertIn("res%block_length_clamped = block_len > requested_block_len",
                      source)
        self.assertNotIn("is shorter than 2*lag_max", source)


class TheAppliedLagIsTheOneTheDataMovedBy(unittest.TestCase):

    def test_it_is_measured_back_off_the_record_shift(self):
        #> Interpolation makes fractional lags ordinary and the data can only
        #> move by whole records, so the requested and the applied lag now
        #> differ by up to half a sample as a matter of course.
        source = code(MODULE)
        self.assertIn("dble(PwbTimelagCache(i)%row_lag) / Metadata%ac_freq",
                      source)


class TheAggregateSummaryCannotSwitchTheMethod(unittest.TestCase):

    def test_a_failed_rh_class_table_leaves_pwb_running(self):
        #> Under PWB the aggregate table is a by-product; the half-hourly
        #> table is the method's output and is complete before this runs.
        #> Switching here discarded all of it and processed the whole run with
        #> covariance maximization, reporting every _TLAG_PWB_SOURCE as
        #> missing.
        source = code("src/src_rp/optimize_timelags.f90")
        block = source[source.index("if (toH2O(1)%def == error"):]
        self.assertIn("Meth%tlag == 'pwb'", block)
        self.assertLess(block.index("Meth%tlag == 'pwb'"),
                        block.index("Meth%tlag = 'maxcov'"))


if __name__ == "__main__":
    unittest.main()
