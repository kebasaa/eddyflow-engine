import math
import re
from pathlib import Path
import unittest


REJECTED = 0
NORMAL = 1
ALL_STOMATAL = 2
ALL_NONSTOMATAL = 3
SINGULAR = 4
WRONG_SIGN = 5
H2O_TO_ET = 0.0648
DEFAULT_SINGULAR_BAND = 0.2
ROOT = Path(__file__).resolve().parents[1]


def read(relative):
    return (ROOT / relative).read_text(encoding="utf-8")


def normalize_percent(value, default):
    """An occupancy limit is a percentage and nothing else.

    This used to accept a fraction too and pick between them by magnitude,
    which made every setting at or below one percent mean a hundred times what
    it said.
    """
    if 0 <= value <= 100:
        return value / 100
    return default


def normalize_band(value, default=DEFAULT_SINGULAR_BAND):
    if 0 <= value <= 1:
        return value
    return default


def normalize_signal_strength(value, default=70):
    if value <= 0:
        return 0
    if value <= 100:
        return value
    return default


def interpolate_short_gaps(values, max_gap=4):
    values = list(values)
    i = 1
    while i < len(values) - 1:
        if values[i] is not None:
            i += 1
            continue
        start = i
        while i < len(values) and values[i] is None:
            i += 1
        length = i - start
        if i < len(values) and length <= max_gap and values[start - 1] is not None:
            step = (values[i] - values[start - 1]) / (length + 1)
            for offset in range(length):
                values[start + offset] = values[start - 1] + step * (offset + 1)
    return values


def filter_signal(values, signal_strength, threshold):
    values = list(values)
    if threshold <= 0 or signal_strength is None:
        return values
    for index, signal in enumerate(signal_strength[: len(values)]):
        if signal is not None and signal < threshold:
            values[index] = None
    return values


def _sigma(values):
    finite = [v for v in values if v is not None]
    if len(finite) < 2:
        return 0.0
    mean = sum(finite) / len(finite)
    return math.sqrt(max(0.0, sum(v * v for v in finite) / len(finite) - mean * mean))


def standard_deviations(w_values, c_values, q_values):
    """Over the samples all three series have, which is what the octants use."""
    triples = [t for t in zip(w_values, c_values, q_values) if None not in t]
    if len(triples) < 2:
        return 0.0, 0.0, 0.0
    return (_sigma([t[0] for t in triples]),
            _sigma([t[1] for t in triples]),
            _sigma([t[2] for t in triples]))


def passes_hole(w, q, c, h, sigma_w, sigma_q, sigma_c):
    """The hyperbolic hole of Thomas et al. (2008), Eqs. 7a/7b.

    H is dimensionless: it scales the instantaneous flux against
    sigma_w * sigma_s, so one value means the same thing at every site. The
    engine used to compare against the raw product, which made H depend on
    whether the gas was carried in umol/mol or mmol/m3.

    The hole is a property of the PAIR. Zahn et al. define the indicator
    functions on (w', q', c') alone and reuse them unchanged for every
    partitioned species; testing each scalar against its own sigma would make
    the octants depend on which scalar was being summed and lose that.
    """
    if h <= 0:
        return True
    return (abs(w * q) >= h * sigma_w * sigma_q
            and abs(w * c) >= h * sigma_w * sigma_c)


def count_octants(w_values, c_values, q_values, h=0):
    sigma_w, sigma_c, sigma_q = standard_deviations(w_values, c_values, q_values)
    valid = o1 = o2 = 0
    for w_value, c_value, q_value in zip(w_values, c_values, q_values):
        if None in (w_value, c_value, q_value):
            continue
        #> A point the hole rejects still counts toward N - the paper
        #> normalises by every sample - but not toward the octant, so raising H
        #> also tightens the occupancy gates. That is Eq. 7, where H sits
        #> inside the indicator function.
        valid += 1
        if w_value > 0 and q_value > 0 and c_value > 0:
            if passes_hole(w_value, q_value, c_value, h, sigma_w, sigma_q, sigma_c):
                o1 += 1
        elif w_value > 0 and q_value > 0 and c_value < 0:
            if passes_hole(w_value, q_value, c_value, h, sigma_w, sigma_q, sigma_c):
                o2 += 1
    return valid, o1, o2


def target_sums(w_values, c_values, q_values, scalar, h=0):
    """One target's two sample fluxes, over the samples that target has.

    The octant mask comes from w, c and q alone - Zahn et al. Sec. 2.4, "the
    indicator functions I_R and I_P remain the same as I_E and I_T" - so a
    target with no value for a sample contributes nothing for that sample and
    does NOT take the sample out of the octant. That is what stops an extra
    species moving the water and carbon it was added beside.

    N is this target's own count rather than the pairing's. The ratio is
    invariant to that choice, so nothing computed downstream moves; it only
    makes f_O1 and f_O2 mean "per sample this target actually had".
    """
    sigma_w, sigma_c, sigma_q = standard_deviations(w_values, c_values, q_values)
    o1 = o2 = 0.0
    n = 0
    for w, c, q, s in zip(w_values, c_values, q_values, scalar):
        if None in (w, c, q) or s is None:
            continue
        n += 1
        if not (w > 0 and q > 0):
            continue
        if not passes_hole(w, q, c, h, sigma_w, sigma_q, sigma_c):
            continue
        if c > 0:
            o1 += w * s
        elif c < 0:
            o2 += w * s
    if n == 0:
        return None, None
    return o1 / n, o2 / n


def partition_stability(sub_c1, sub_c2, sub_n, whole_O1, whole_O2):
    """The corrected statistic: how far the octant split moves.

    Both splits are normalised by |f_O1| + |f_O2|, which cannot vanish while
    the octants hold anything - unlike Foken's denominator, which is the
    covariance under test and goes to zero with the night-time carbon flux.
    A trend that scales both octants together moves neither split and
    cancels, which is the robustness the ratio was meant to buy.

    Bounded by 200. None where fewer than three sixths could be centred.
    """
    denom = abs(whole_O1) + abs(whole_O2)
    if denom == 0:
        return None
    split_whole = whole_O1 / denom
    usable = [(a, b, n) for a, b, n in zip(sub_c1, sub_c2, sub_n) if n >= 2]
    if len(usable) < 3:
        return None
    a1 = sum(a / n for a, b, n in usable) / len(usable)
    a2 = sum(b / n for a, b, n in usable) / len(usable)
    denom = abs(a1) + abs(a2)
    if denom == 0:
        return None
    return abs(split_whole - a1 / denom) * 100.0


def flux_is_unresolved(total, err, k):
    """Is this flux distinguishable from zero?

    Finkelstein & Sims give |F|/RE ~ |r| * sqrt(N_indep / 2), so comparing a
    flux against a multiple of its own random error IS the significance of the
    w-scalar correlation, with the number of independent samples taken from
    the period's own integral timescale rather than assumed.

    An absent error is no opinion, never a failure: a run with the estimator
    switched off leaves every random error unset, and refusing every period on
    that basis would turn a switched-off diagnostic into a switched-off
    partition.
    """
    if k <= 0:
        return False
    if err is None or err <= 0:
        return False
    if total is None:
        return False
    return abs(total) < k * err


def apply_partition(status, ratio, total):
    if status == NORMAL:
        return total / (1 + 1 / ratio), total / (1 + ratio)
    if status == ALL_STOMATAL:
        return 0.0, total
    if status == ALL_NONSTOMATAL:
        return total, 0.0
    return math.nan, math.nan


def carbon_status(ratio, band=DEFAULT_SINGULAR_BAND):
    """R and P nearly cancelling puts 1 + r on top of zero.

    Zahn et al. reject -1.2 < r < -0.8 and say the width is dataset-dependent,
    so it is a setting. The engine hard-coded 0.05, which this very module
    already contradicted.

    Exactly -1 is rejected whatever the band, including none: it is a division
    by zero rather than a judgement about how near the singularity is too near.
    """
    if ratio == -1:
        return SINGULAR
    if band > 0 and abs(ratio + 1) < band:
        return SINGULAR
    return NORMAL


def sign_checked_status(descriptor_status, f_o1, f_o2, total):
    """The total has to point the way the two sample fluxes point together.

    Substituting the ratio into the partition and cancelling gives

        nonstomatal = total * f_O1 / (f_O1 + f_O2)
        stomatal    = total * f_O2 / (f_O1 + f_O2)

    so each component carries the sign of its own sample flux exactly when
    sign(total) == sign(f_O1 + f_O2). Otherwise the arithmetic still returns
    two numbers that sum to the total, and they are a negative respiration
    beside a positive photosynthesis.

    One rule for every species. Water is the case where both sample fluxes are
    positive by construction, so a downward total - dewfall - is always caught.
    Carbon is the case where they have opposite signs and can cancel, so it is
    caught only when they cancel the wrong way.

    Only the ratio branch is guarded. Where one octant holds too few points the
    paper hands the whole flux to the other component; that branch divides by
    nothing and leaves the total's sign on the single component it fills, so
    there is no inversion here to catch and vetoing it would only throw away
    periods - the nighttime ones above all, which is where CEC does most of its
    work.
    """
    if descriptor_status != NORMAL:
        return descriptor_status
    if total in (None, 0):
        return descriptor_status
    sampled = f_o1 + f_o2
    if sampled == 0:
        return descriptor_status
    if sampled * total < 0:
        return WRONG_SIGN
    return descriptor_status


class CecReferenceTests(unittest.TestCase):
    def test_interpolates_at_most_four_internal_samples(self):
        self.assertEqual(
            interpolate_short_gaps([0.0, None, None, None, None, 5.0]),
            [0.0, 1.0, 2.0, 3.0, 4.0, 5.0],
        )
        self.assertEqual(
            interpolate_short_gaps([0.0, None, None, None, None, None, 6.0]),
            [0.0, None, None, None, None, None, 6.0],
        )
        self.assertEqual(
            interpolate_short_gaps([0.0, None, 2.0], max_gap=0),
            [0.0, None, 2.0],
        )

    def test_completeness_and_stationarity_boundaries(self):
        self.assertTrue(10 * 90 >= 9 * 100)
        self.assertFalse(10 * 89 >= 9 * 100)
        self.assertEqual(normalize_percent(90, 0.5), 0.90)
        self.assertEqual(normalize_percent(-1, 0.5), 0.5)
        self.assertTrue(25 <= 25)
        self.assertFalse(26 <= 25)

    def test_signal_strength_filters_only_when_sample_data_are_available(self):
        values = [1.0, 2.0, 3.0]
        self.assertEqual(filter_signal(values, None, 70), values)
        self.assertEqual(filter_signal(values, [71, 69, 70], 70), [1.0, None, 3.0])
        self.assertEqual(filter_signal(values, [10, 20, 30], 0), values)
        self.assertEqual(normalize_signal_strength(110), 70)

    def test_h_threshold_reduces_eligible_octant_events(self):
        w_values = [1.0, 1.0, 1.0]
        q_values = [0.01, 0.5, 0.5]
        c_values = [0.5, 0.01, -0.5]
        self.assertEqual(count_octants(w_values, c_values, q_values, h=0), (3, 2, 1))
        #> The near-origin events go first, and they go at a threshold that
        #> means the same thing whatever units the gas is carried in.
        sigma_w, sigma_c, sigma_q = standard_deviations(w_values, c_values, q_values)
        self.assertEqual(sigma_w, 0.0)
        thinned = count_octants(w_values, c_values, q_values, h=1.0)
        self.assertEqual(thinned[0], 3)
        self.assertLessEqual(thinned[1] + thinned[2], 3)

    def test_the_hole_is_scaled_by_sigma_not_by_the_raw_product(self):
        #> The same turbulence described in two unit systems must give the same
        #> octants. Multiply the gas by a thousand and nothing may move.
        w_values = [0.4, -0.3, 0.9, 0.2, -0.6, 0.5]
        q_values = [0.02, 0.5, 0.4, 0.01, -0.2, 0.3]
        c_values = [0.3, 0.02, -0.4, 0.005, 0.1, -0.2]
        scaled_c = [v * 1000 for v in c_values]
        scaled_q = [v * 1000 for v in q_values]
        for h in (0.0, 0.25, 1.0):
            self.assertEqual(
                count_octants(w_values, c_values, q_values, h=h),
                count_octants(w_values, scaled_c, scaled_q, h=h),
                f"the hole is not unit-invariant at H={h}")

    def test_the_hole_is_one_test_for_the_pair_not_one_per_scalar(self):
        #> The invariant that lets an arbitrary species be partitioned: the
        #> indicator functions are built from (w', q', c') and nothing else, so
        #> every target shares one set of octants.
        source = read("src/src_common/m_cec.f90")
        self.assertIn("CecPassesHyperbolicThreshold(primes(i, 1), primes(i, iw), &",
                      source)
        self.assertIn("primes(i, ic), setup%h, sigma_w, sigma_q, sigma_c", source)
        #> The octant test names w, the water column and the carbon column and
        #> nothing else. An extra species contributes two sums inside those
        #> octants and has no say in which points are in them.
        self.assertIn("sum_O1(k) = sum_O1(k) + primes(i, 1) * primes(i, k + 1)",
                      source)
        self.assertIn("sum_O2(k) = sum_O2(k) + primes(i, 1) * primes(i, k + 1)",
                      source)

    def test_the_screen_runs_before_the_gap_filler(self):
        """Zahn et al. Sec. 3.2, in the order they state it.

        "CO2 or H2O measurements with a signal strength [below] 70% ... were
        deleted. Small gaps (up to four consecutive points) were filled by
        linear interpolation." The gaps worth filling are the ones the screen
        has just made, so a brief dropout costs four interpolated samples
        instead of costing the whole period at the 90% completeness gate.

        This ran the other way round and was described in the source as the
        paper's order. It was not. cec_max_gap_fill = 0 is how a project asks
        for the strict reading, where a condemned sample is never rebuilt.
        """
        source = read("src/src_common/m_cec.f90")
        screen = source.index("call FilterCecSignalStrength(work(:, k + 1)")
        gap_fill = source.index("call InterpolateShortCecGaps(work(:, k)")
        self.assertLess(screen, gap_fill,
                        "the signal-strength screen must precede the gap filler")
        #> And both must precede the detrending, which is the whole reason the
        #> partition builds its own compact copy rather than reading E2Primes.
        self.assertLess(gap_fill, source.index("call Fluctuations(work,"))

    def test_a_total_pointing_against_the_octants_is_flagged_not_split(self):
        #> Water: both sample fluxes positive, so a downward total is dewfall
        #> and cannot be split into evaporation and transpiration.
        self.assertEqual(sign_checked_status(NORMAL, 0.4, 0.6, 3.0), NORMAL)
        self.assertEqual(sign_checked_status(NORMAL, 0.4, 0.6, -3.0), WRONG_SIGN)

        #> But NOT where the octant was too sparse to divide by. Zahn et al.
        #> Sec. 2.4 hand the whole flux to the other component there, and that
        #> assignment cannot invert a sign because it performs no division -
        #> the component IS the total. Vetoing it would discard exactly the
        #> nighttime periods the fallback exists to carry.
        self.assertEqual(sign_checked_status(ALL_STOMATAL, 0.4, 0.6, -0.1),
                         ALL_STOMATAL)
        self.assertEqual(sign_checked_status(ALL_NONSTOMATAL, 0.4, 0.6, -0.1),
                         ALL_NONSTOMATAL)
        for status in (ALL_STOMATAL, ALL_NONSTOMATAL):
            first, second = apply_partition(status, 0.0, -0.1)
            self.assertAlmostEqual(first + second, -0.1)
            self.assertLessEqual(first, 0.0)
            self.assertLessEqual(second, 0.0)

        #> Carbon, daytime: a little respiration against a lot of uptake, and
        #> a downward total. Consistent.
        self.assertEqual(sign_checked_status(NORMAL, 0.2, -0.8, -9.0), NORMAL)
        #> Carbon, nighttime: respiration dominates and the total is upward.
        self.assertEqual(sign_checked_status(NORMAL, 0.8, -0.2, 4.0), NORMAL)
        #> Carbon, contradictory: the octants say net upward - r_Fc below -1 -
        #> while the corrected total is downward. Splitting it here returns a
        #> NEGATIVE respiration beside a positive photosynthesis, which is what
        #> the engine did on the first real period this was run against.
        self.assertEqual(sign_checked_status(NORMAL, 0.8, -0.536, -1.80),
                         WRONG_SIGN)
        nonstomatal, stomatal = apply_partition(NORMAL, 0.8 / -0.536, -1.80)
        self.assertLess(nonstomatal, 0.0)
        self.assertGreater(stomatal, 0.0)
        self.assertAlmostEqual(nonstomatal + stomatal, -1.80)
        source = read("src/src_common/m_cec.f90")
        self.assertIn("flux%comp(k)%status = cec_wrong_sign", source)
        #> Stated once and generally, over the sum of the two sample fluxes, so
        #> it holds for water, for carbon and for any species partitioned in
        #> the same octants.
        self.assertIn("sampled = target%f_O1 + target%f_O2", source)
        self.assertIn("CecTotalContradictsOctants = sampled * total < 0d0",
                      source)

        #> Inside the ratio arm, not ahead of the select. Asserted by slicing
        #> the arm out, because the call sitting anywhere in the subroutine is
        #> what it used to do and would still satisfy a bare `assertIn`.
        arm = source[source.index("            case (cec_normal)"):
                     source.index("            case (cec_all_stomatal)")]
        self.assertIn("if (CecTotalContradictsOctants(descriptor%target(k), total))",
                      arm,
                      "the sign guard no longer sits inside the ratio branch")
        head = source[source.index("subroutine ApplyCecDescriptor"):
                      source.index("            case (cec_normal)")]
        self.assertNotIn("CecTotalContradictsOctants", head,
                         "the sign guard is back ahead of the select, where it "
                         "also vetoes the sparse-octant fallbacks")

    def test_partitioning_an_extra_scalar_needs_no_new_octants(self):
        #> Zahn et al. Sec 2.4: "the indicator functions I_R and I_P remain the
        #> same as I_E and I_T". Any scalar can therefore be partitioned by
        #> summing it over the octants the pair already defines - which is what
        #> makes COS a drop-in rather than a second implementation.
        w_values = [0.4, -0.3, 0.9, 0.2, -0.6, 0.5]
        q_values = [0.2, 0.5, 0.4, 0.1, -0.2, 0.3]
        c_values = [0.3, 0.2, -0.4, 0.05, 0.1, -0.2]
        cos_values = [-0.1, 0.4, -0.3, 0.02, 0.3, -0.05]

        def sums(scalar):
            o1 = o2 = 0.0
            n = len(w_values)
            for w, c, q, s in zip(w_values, c_values, q_values, scalar):
                if w > 0 and q > 0 and c > 0:
                    o1 += w * s
                elif w > 0 and q > 0 and c < 0:
                    o2 += w * s
            return o1 / n, o2 / n

        for scalar, total in ((q_values, 4.0), (c_values, -9.0), (cos_values, -2.0)):
            f_o1, f_o2 = sums(scalar)
            ratio = f_o1 / f_o2
            nonstomatal, stomatal = apply_partition(carbon_status(ratio), ratio, total)
            self.assertAlmostEqual(nonstomatal + stomatal, total)

    def test_an_extra_species_cannot_move_the_pair_it_was_added_to(self):
        """Adding COS to a pairing must leave its water and carbon untouched.

        The sample gate names w, the water column and the carbon column and
        nothing else, so a target that is missing for a sample no longer takes
        that sample out of the period. Before this, one gappy extra species
        shrank n_valid for every target and silently revised the E/T/R/P a
        pairing had already published.
        """
        source = read("src/src_common/m_cec.f90")
        body = source[source.index("subroutine ExtractCecDescriptor"):
                      source.index("end subroutine ExtractCecDescriptor")]
        for column in ("primes(i, 1)", "primes(i, iw)", "primes(i, ic)"):
            self.assertIn("if (.not. CecValueIsValid(%s)) cycle" % column, body)
        #> The old all-targets veto. Named exactly, because it was the
        #> variable that has to be gone, not the word - the first version of
        #> this matched any prose containing "usable" and tripped on a comment.
        self.assertNotIn("logical :: usable", body)
        self.assertNotIn("if (.not. usable) cycle", body,
                         "a target other than w, q and c can veto a sample "
                         "again, so extras move the pair once more")
        #> And each target divides by its own count, not the pairing's.
        for line in ("descriptor%target(k)%f_O1 = sum_O1(k) / dble(n_target_valid(k))",
                     "descriptor%target(k)%f_O2 = sum_O2(k) / dble(n_target_valid(k))"):
            self.assertIn(line, body)

    def test_a_targets_own_sample_count_leaves_its_ratio_alone(self):
        """Which is why normalising per target moves no flux.

        f_O1 and f_O2 both scale with 1/N, so the ratio - the only thing
        ApplyCecDescriptor reads - is invariant. The magnitudes change, and
        nothing downstream reads them for magnitude.
        """
        w_values = [0.4, -0.3, 0.9, 0.2, -0.6, 0.5]
        q_values = [0.2, 0.5, 0.4, 0.1, -0.2, 0.3]
        c_values = [0.3, 0.2, -0.4, 0.05, 0.1, -0.2]
        cos_values = [-0.1, 0.4, -0.3, 0.02, 0.3, -0.05]

        #> Blank the extra species only where the sample is in no octant, so
        #> the two sums are identical and only N differs.
        gappy = list(cos_values)
        for i, (w, q) in enumerate(zip(w_values, q_values)):
            if not (w > 0 and q > 0):
                gappy[i] = None
        self.assertIn(None, gappy, "the fixture stopped exercising the gap")

        full_o1, full_o2 = target_sums(w_values, c_values, q_values, cos_values)
        gap_o1, gap_o2 = target_sums(w_values, c_values, q_values, gappy)
        self.assertNotAlmostEqual(full_o1, gap_o1)
        self.assertAlmostEqual(full_o1 / full_o2, gap_o1 / gap_o2)

        #> And the water and carbon of the same pairing are untouched by either.
        for scalar in (q_values, c_values):
            self.assertEqual(target_sums(w_values, c_values, q_values, scalar),
                             target_sums(w_values, c_values, q_values, scalar))

    def test_a_rejected_period_still_reports_its_octants(self):
        """Counts and fractions describe the period either way.

        The counts are accumulated in the sample loop whatever happens next,
        so computing the fractions after the completeness and stationarity
        gates left a rejected period publishing half its diagnostics: counts
        present, fractions at error, and no way to divide one by the other
        because n_valid is not written out. On the run this was found in, 30
        of 48 periods came out that way.

        A period that did not partition is exactly when someone wants to read
        these, so they are computed before the gates. The occupancy gate still
        sits after them, because it is the one that reads them.
        """
        source = read("src/src_common/m_cec.f90")
        body = source[source.index("subroutine ExtractCecDescriptor"):
                      source.index("end subroutine ExtractCecDescriptor")]
        frac = body.index("descriptor%frac_O1 = dble(descriptor%n_O1)")
        complete = body.index("< setup%min_valid * dble(nrow)) return")
        stationary = body.index("> setup%max_stationarity")
        occupancy = body.index("< setup%min_o1_o2) return")
        self.assertLess(frac, complete, "the completeness gate returns first")
        self.assertLess(frac, stationary, "the stationarity gate returns first")
        self.assertLess(complete, occupancy,
                        "the occupancy gate must still read them, so it stays last")
        #> And it cannot divide by zero on the way.
        self.assertIn("if (descriptor%n_valid > 0) then", body)

    def test_a_gas_with_no_diagnostic_skips_the_screen_not_the_partition(self):
        """No AGC or RSSI column is a normal site, not a broken one.

        Most analysers report no signal strength at all - a quantum cascade
        laser has nothing to report - and the partition does not depend on
        one. The lookup returning nothing must therefore skip that target's
        screen and carry on, never abandon the pairing.
        """
        source = read("src/src_common/m_cec.f90")
        body = source[source.index("subroutine BuildCecPrimes"):
                      source.index("end subroutine BuildCecPrimes")]
        block = body[body.index("if (setup%signal_strength > 0d0"):]
        block = block[:block.index("call Fluctuations(work,")]
        self.assertIn("if (sig_col <= 0 .or. sig_col > nuser) cycle", block)
        self.assertNotIn("return", block,
                         "a missing diagnostic column now abandons the pairing "
                         "instead of skipping the screen")

    def test_each_sixth_is_recentred_on_its_own_mean(self):
        """Without this the statistic measures nothing whatsoever.

        `primes` already holds whole-period fluctuations. Slice them without
        re-centring and the mean of the six sub-interval values equals the
        whole-period value identically, so the difference between them - which
        IS the trend signal, and is the whole of Foken's test - is exactly
        zero by construction. The first version of this did that, and read
        70000% on periods Foken called stationary; all that survived was the
        nonlinearity of a ratio of noisy sixths.

        StationarityTest re-centres because it slices the raw set and takes
        each slice's covariance about that slice's own mean. This has to do
        the same by hand.
        """
        engine = read("src/src_rp/stationarity_test.f90")
        self.assertIn("call CovarianceMatrixNoError(SubSet, subn, GHGNumVar, SubCov, error)",
                      engine)
        cec = read("src/src_common/m_cec.f90")
        block = cec[cec.index("subroutine CecPartitionStability"):
                    cec.index("end subroutine CecPartitionStability")]
        self.assertIn("(primes(i, 1) - wbar) * (primes(i, j + 1) - sbar)", block)
        self.assertIn("wbar = sw(sub, j) / dble(n_all(sub, j))", block)

        #> And the mean is taken over every sample the sixth can use, not only
        #> its octant members - centring on those would remove the very
        #> asymmetry the octant selects for.
        first = block[block.index("n_all = 0"):block.index("c1 = 0d0")]
        self.assertIn("if (octant(i) < 0) cycle", first)
        self.assertNotIn("if (octant(i) <= 0) cycle", first)

    def test_the_denominator_cannot_vanish(self):
        """The other half of the fix, and the half that answers the problem.

        Foken divides by the covariance whose stationarity he is testing, so
        his statistic runs away as that covariance approaches zero - which for
        carbon at night it does, and that is what was rejecting half a day of
        well-sampled partitions. Both splits here are normalised by
        |f_O1| + |f_O2|, which cannot vanish while the octants hold anything.
        """
        cec = read("src/src_common/m_cec.f90")
        block = cec[cec.index("subroutine CecPartitionStability"):
                    cec.index("end subroutine CecPartitionStability")]
        self.assertEqual(block.count("denom = dabs(f_O1(j)) + dabs(f_O2(j))"), 1)
        self.assertEqual(block.count("denom = dabs(a1) + dabs(a2)"), 1)
        self.assertIn("split_whole = f_O1(j) / denom", block)
        self.assertIn("split_sub = a1 / denom", block)
        #> Never by the quantity being tested, which is the mistake being
        #> corrected. The statistic is bounded by 200 as a result.
        self.assertNotIn("/ f_O1(j)", block)
        self.assertNotIn("/ r_whole", block)

    def test_it_divides_the_period_as_the_existing_test_does(self):
        #> Six, or the two numbers are not comparable on the same run.
        self.assertIn("integer, parameter :: ndiv = 6",
                      read("src/src_rp/stationarity_test.f90"))
        self.assertIn("integer, parameter :: ncecdiv = 6",
                      read("src/src_common/m_cec.f90"))

    def test_a_steady_split_reads_zero_however_the_flux_drifts(self):
        """The property the whole mode rests on.

        Both octants doubling through the period is a drifting flux and a
        perfectly steady partition. Foken's statistic on the covariance sees
        the drift; this one sees the partition, and the partition is what CEC
        multiplies into the total.
        """
        n = [10] * 6
        c1 = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
        c2 = [2.0, 4.0, 6.0, 8.0, 10.0, 12.0]
        whole1 = sum(c1) / 6.0
        whole2 = sum(c2) / 6.0
        self.assertAlmostEqual(
            partition_stability([a * 10 for a in c1], [b * 10 for b in c2],
                                n, whole1, whole2), 0.0)

    def test_a_split_that_moves_is_what_it_catches(self):
        #> Non-stomatal dominant early, stomatal dominant late: the same total
        #> activity throughout, a partition that inverts.
        n = [10] * 6
        c1 = [9.0, 8.0, 7.0, 3.0, 2.0, 1.0]
        c2 = [1.0, 2.0, 3.0, 7.0, 8.0, 9.0]
        whole1 = sum(c1) / 6.0
        whole2 = sum(c2) / 6.0
        got = partition_stability([a * 10 for a in c1], [b * 10 for b in c2],
                                  n, whole1, whole2)
        self.assertGreaterEqual(got, 0.0)
        self.assertLessEqual(got, 200.0, "the statistic is meant to be bounded")

    def test_it_is_bounded_where_the_old_one_ran_to_70000(self):
        """The failure that sent the first attempt back.

        One sixth whose stomatal octant nearly cancels used to send a ratio to
        infinity and drag the mean of ratios with it. Averaging the
        covariances and normalising by the partition's own size cannot do
        that: the result never leaves [0, 200].
        """
        n = [10] * 6
        c1 = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
        c2 = [1.0, 1.0, 1.0, 1.0, 1.0, 1e-12]
        got = partition_stability([a * 10 for a in c1], [b * 10 for b in c2],
                                  n, sum(c1) / 6.0, sum(c2) / 6.0)
        self.assertLessEqual(got, 200.0)

    def test_too_few_sub_intervals_is_undefined_not_failed(self):
        """Rejecting on a number that could not be computed is the fault this
        mode exists to correct, so the engine reads error as "no opinion"."""
        self.assertIsNone(partition_stability([1.0] * 6, [1.0] * 6,
                                              [10, 10, 1, 1, 1, 1], 1.0, 1.0))

        body = read("src/src_common/m_cec.f90")
        gate = body[body.index("if (setup%stationarity_mode == cec_stat_ratio) then"):]
        gate = gate[:gate.index("else")]
        #> Guarded on /= error before it is compared, so an undefined value
        #> cannot trip the threshold.
        self.assertIn("%ns_r /= error", gate)
        for target in ("cecTargetWater", "cecTargetCarbon"):
            self.assertIn("descriptor%%target(%s)%%ns_r" % target, gate)

    def test_the_published_criterion_is_what_a_project_gets_by_default(self):
        types = read("src/src_common/m_typedef.f90")
        self.assertIn("integer, parameter :: cec_stat_flux = 0", types)
        self.assertIn("integer, parameter :: cec_stat_ratio = 1", types)

        source = read("src/src_common/write_processing_project_variables.f90")
        #> Sliced to the block that sets the defaults, because the same line
        #> also appears in the parser's fallback below and an unsliced
        #> assertIn would pass on that alone - which would leave a project
        #> that states nothing silently taking the new mode.
        defaults = source[:source.index("if (EPPrjNTagFound(")]
        self.assertIn("EddyFlowProj%cec%stationarity_mode = cec_stat_flux", defaults)

        #> And an unrecognised value falls back to the paper rather than to
        #> the new mode. The safe direction to be wrong in.
        parser = source[source.index("if (EPPrjNTagFound(5)) then"):]
        parser = parser[:parser.index("end if", parser.index("else"))]
        self.assertIn("nint(EPPrjNTags(5)%value) == cec_stat_ratio", parser)
        self.assertIn("EddyFlowProj%cec%stationarity_mode = cec_stat_flux", parser)

    def test_the_mode_switches_the_gate_and_nothing_else(self):
        """It replaces one gate. Completeness and occupancy are unconditional,
        and a period that fails those fails in either mode."""
        body = read("src/src_common/m_cec.f90")
        extract = body[body.index("subroutine ExtractCecDescriptor"):
                       body.index("end subroutine ExtractCecDescriptor")]
        mode = extract.index("setup%stationarity_mode")
        complete = extract.index("< setup%min_valid * dble(nrow)) return")
        occupancy = extract.index("< setup%min_o1_o2) return")
        self.assertLess(complete, mode, "completeness moved inside the mode")
        self.assertLess(mode, occupancy, "occupancy moved inside the mode")
        gated = extract[extract.index("if (setup%max_stationarity > 0d0) then"):occupancy]
        self.assertNotIn("min_valid", gated)
        self.assertNotIn("min_o1_o2", gated)

    def test_the_statistic_is_reported_whichever_gate_is_chosen(self):
        #> So the two criteria can be compared on one run, and so a rejected
        #> period still says what it was judged on - the same reason the
        #> octant fractions moved above the gates.
        body = read("src/src_common/m_cec.f90")
        extract = body[body.index("subroutine ExtractCecDescriptor"):
                       body.index("end subroutine ExtractCecDescriptor")]
        assigned = extract.index("call CecPartitionStability(")
        self.assertLess(assigned,
                        extract.index("if (setup%max_stationarity > 0d0) then"),
                        "the statistic is computed after the gate that reads it")
        self.assertLess(assigned,
                        extract.index("< setup%min_valid * dble(nrow)) return"),
                        "a period rejected for completeness reports no statistic")

    def test_an_unresolved_flux_is_not_worth_partitioning(self):
        """The failure this exists to stop.

        The octants mean something only if the sign of c\' carries a surface
        signature. On the record this was written for, night |r(w,CO2)| ran at
        0.079 against 0.336 by day - below one sigma, so the moist ejections
        split near evenly and O2 never emptied. The partition would still
        return two numbers summing to the total; they would mean nothing.
        """
        #> A flux ten times its own error is real; one a third of it is not.
        self.assertFalse(flux_is_unresolved(10.0, 1.0, 2.0))
        self.assertTrue(flux_is_unresolved(0.3, 1.0, 2.0))
        #> Sign does not matter - a downward night flux is as unresolved as an
        #> upward one of the same size.
        self.assertTrue(flux_is_unresolved(-0.3, 1.0, 2.0))

    def test_the_test_is_off_by_default_and_off_means_off(self):
        self.assertFalse(flux_is_unresolved(0.001, 1.0, 0.0))
        source = read("src/src_common/write_processing_project_variables.f90")
        defaults = source[:source.index("if (EPPrjNTagFound(")]
        self.assertIn("EddyFlowProj%cec%min_flux_sigma = 0d0", defaults)
        body = read("src/src_common/m_cec.f90")
        fn = body[body.index("logical function CecFluxIsUnresolved"):
                  body.index("end function CecFluxIsUnresolved")]
        self.assertIn("if (setup%min_flux_sigma <= 0d0) return", fn)

    def test_a_missing_random_error_is_no_opinion_not_a_failure(self):
        """ru_meth = 0 leaves every random error unset. Refusing every period
        on that basis would turn a switched-off diagnostic into a switched-off
        partition, so the engine warns instead and the test abstains."""
        self.assertFalse(flux_is_unresolved(0.001, None, 2.0))
        body = read("src/src_common/m_cec.f90")
        fn = body[body.index("logical function CecFluxIsUnresolved"):
                  body.index("end function CecFluxIsUnresolved")]
        self.assertIn("if (.not. CecValueIsValid(err)) return", fn)
        self.assertIn("if (err <= 0d0) return", fn)
        #> And it is said out loud, once, where it can still be acted on.
        self.assertIn("call ExceptionHandler(115)",
                      read("src/src_rp/read_ini_rp.f90"))
        self.assertIn("Warning(115)", read("src/src_common/exception_handler.f90"))

    def test_the_pairing_is_judged_on_both_of_its_channels(self):
        """The pool is the moist ejections and the split is the sign of c\',
        so an unresolved water flux means "moist ejection" is not selecting
        surface-influenced air, and an unresolved carbon flux means the split
        is a coin toss. Either one sinks the pairing."""
        body = read("src/src_common/m_cec.f90")
        fn = body[body.index("logical function CecPairingIsUnresolved"):
                  body.index("end function CecPairingIsUnresolved")]
        self.assertIn("cecTargetWater", fn)
        self.assertIn("cecTargetCarbon", fn)
        self.assertIn(".or.", fn)

    def test_the_pairing_check_precedes_the_per_target_loop(self):
        #> It rejects the whole pairing, so it cannot live inside the loop
        #> over targets - and the extras' own check must live inside it.
        body = read("src/src_common/m_cec.f90")
        apply = body[body.index("subroutine ApplyCecDescriptor"):
                     body.index("end subroutine ApplyCecDescriptor")]
        pairing = apply.index("if (CecPairingIsUnresolved(")
        loop = apply.index("do k = 1, descriptor%n_target", pairing)
        self.assertLess(pairing, loop)
        #> The extras gate themselves, as they do for completeness.
        extra = apply.index("if (k /= cecTargetWater .and. k /= cecTargetCarbon) then")
        self.assertLess(loop, extra)

    def test_the_refusal_has_its_own_name(self):
        """"No resolvable flux" is a different statement from "not enough
        data", and only one of them is worth revisiting with a longer record."""
        types = read("src/src_common/m_typedef.f90")
        #> Read out of the source and compared to each other, so that reusing a
        #> code already spoken for fails here rather than silently relabelling
        #> some other verdict in qc_cec_*.
        codes = dict((m.group(1), int(m.group(2))) for m in re.finditer(
            r"integer, parameter :: (cec_[a-z_]+) = (\d+)", types))
        self.assertEqual(codes.pop("cec_insignificant", None), 6)
        others = sorted(v for k, v in codes.items()
                        if not k.startswith("cec_stat_"))
        self.assertNotIn(6, others)
        self.assertEqual(others, list(range(len(others))),
                         "the codes are a contiguous block; 6 was the next free")
        body = read("src/src_common/m_cec.f90")
        apply = body[body.index("subroutine ApplyCecDescriptor"):
                     body.index("end subroutine ApplyCecDescriptor")]
        self.assertEqual(apply.count("cec_insignificant"), 2,
                         "the pairing refusal and the extra refusal, no more")

    def test_the_papers_own_refusal_is_not_relabelled(self):
        """A period the occupancy gate or the singularity band already refused
        was not lost to this test. Overwriting its reason would inflate what
        the test appears to cost and hide the reason worth acting on."""
        body = read("src/src_common/m_cec.f90")
        apply = body[body.index("subroutine ApplyCecDescriptor"):
                     body.index("end subroutine ApplyCecDescriptor")]
        branch = apply[apply.index("if (CecPairingIsUnresolved("):
                       apply.index("do k = 1, descriptor%n_target",
                                   apply.index("return", apply.index(
                                       "if (CecPairingIsUnresolved(")))]
        self.assertIn("flux%comp(k)%status = descriptor%target(k)%status", branch)
        #> valid is false for cec_rejected and for cec_singular alike, so one
        #> test covers both of the paper's refusals.
        self.assertIn("if (descriptor%target(k)%valid) &", branch)

    def test_each_error_comes_from_the_slot_its_total_came_from(self):
        """A mismatched index would test one gas's flux against another gas's
        error and report a verdict that means nothing, silently."""
        rp = read("src/src_rp/eddyflow-rp_main.f90")
        block = rp[rp.index("cec_totals(cec_k) = Flux3%gas(cec_slots(cec_k))"):]
        block = block[:block.index("end do")]
        self.assertIn("cec_errors(cec_k) = Essentials%rand_uncer(cec_slots(cec_k))",
                      block)
        fcc = read("src/src_fcc/eddyflow-fcc_main.f90")
        block = fcc[fcc.index("cec_totals(cec_k) = Flux3%gas(cec_slot)"):]
        block = block[:block.index("end do")]
        self.assertIn("cec_errors(cec_k) = lEx%rand_uncer(cec_slot)", block)

    def test_it_applies_whichever_stationarity_mode_is_chosen(self):
        #> The failure is not mode-specific. Mode 0 merely happens to reject
        #> those periods today for a different reason, and would stop the
        #> moment someone set cec_max_stationarity to 0.
        body = read("src/src_common/m_cec.f90")
        apply = body[body.index("subroutine ApplyCecDescriptor"):
                     body.index("end subroutine ApplyCecDescriptor")]
        self.assertNotIn("stationarity_mode", apply)

    def test_normal_partition_conserves_totals(self):
        evaporation, transpiration = apply_partition(NORMAL, 0.5, 3.0)
        respiration, photosynthesis = apply_partition(NORMAL, -0.25, -12.0)
        self.assertAlmostEqual(evaporation + transpiration, 3.0)
        self.assertAlmostEqual(respiration + photosynthesis, -12.0)

    def test_water_depth_outputs_use_h2o_conversion(self):
        evaporation, transpiration = apply_partition(NORMAL, 0.5, 3.0)
        evaporation_et = evaporation * H2O_TO_ET
        transpiration_et = transpiration * H2O_TO_ET
        self.assertAlmostEqual(evaporation_et + transpiration_et, 3.0 * H2O_TO_ET)

    def test_sparse_octant_fallbacks_conserve_totals(self):
        for status in (ALL_STOMATAL, ALL_NONSTOMATAL):
            first, second = apply_partition(status, 0.0, 4.0)
            self.assertAlmostEqual(first + second, 4.0)

    def test_singular_carbon_ratio_is_rejected(self):
        for ratio in (-1.199, -1.02, -0.801):
            self.assertEqual(carbon_status(ratio), SINGULAR)
        #> Not -1.2/-0.8 exactly: the engine tests abs(r + 1) < band and so
        #> does this, and 1 - 0.8 is not 0.2 in binary. The knife edge is not a
        #> physical distinction, so the check stays off it.
        for ratio in (-1.21, -0.79, -0.5, -2.0):
            self.assertEqual(carbon_status(ratio), NORMAL)

    def test_a_degenerate_ratio_is_named_rather_than_divided_by(self):
        """Both ends of the ratio are a division by zero, and both have a
        right answer that does not need one.

        f_O2 = 0 is "nothing in the stomatal octant", which is the whole flux
        to the other component - the same answer 1/(1 + 1/inf) gives. f_O1 = 0
        is its mirror. And r = -1 exactly is a division by zero however wide
        the rejection band is, including when the band is switched off.
        """
        source = read("src/src_common/m_cec.f90")
        self.assertIn("else if (target%f_O2 == 0d0) then", source)
        self.assertIn("else if (target%f_O1 == 0d0) then", source)
        self.assertIn("if (target%r == -1d0) then", source)

        #> The reference model agrees on the answers those arms give.
        self.assertAlmostEqual(sum(apply_partition(ALL_NONSTOMATAL, 0.0, 4.0)), 4.0)
        self.assertEqual(apply_partition(ALL_NONSTOMATAL, 0.0, 4.0), (4.0, 0.0))
        self.assertEqual(apply_partition(ALL_STOMATAL, 0.0, 4.0), (0.0, 4.0))
        #> And that r = -1 is singular at any band, including none.
        self.assertEqual(carbon_status(-1.0, band=0.0), SINGULAR)
        self.assertEqual(carbon_status(-1.0, band=DEFAULT_SINGULAR_BAND), SINGULAR)

    def test_the_singular_band_is_a_setting(self):
        #> The engine used 0.05 while this module used 0.2, so the reference
        #> and the implementation disagreed by a factor of four on which
        #> periods survive. Now both read the setting.
        self.assertEqual(carbon_status(-0.9, band=0.05), NORMAL)
        self.assertEqual(carbon_status(-0.9, band=0.2), SINGULAR)
        #> Switching the band off leaves the near-singular periods in; only the
        #> exact division by zero is still refused.
        self.assertEqual(carbon_status(-0.99, band=0.0), NORMAL)
        self.assertEqual(normalize_band(0.35), 0.35)
        self.assertEqual(normalize_band(7.0), DEFAULT_SINGULAR_BAND)

    def test_occupancy_limits_are_percentages_only(self):
        self.assertEqual(normalize_percent(90, 0.5), 0.90)
        self.assertEqual(normalize_percent(20, 0.5), 0.20)
        #> The reading that used to be a fraction. Half a percent is half a
        #> percent, not half.
        self.assertEqual(normalize_percent(0.5, 0.9), 0.005)
        self.assertEqual(normalize_percent(-1, 0.5), 0.5)

    def test_fcc_scaling_preserves_partition(self):
        rp = apply_partition(NORMAL, 0.5, 3.0)
        fcc = apply_partition(NORMAL, 0.5, 4.5)
        self.assertAlmostEqual(fcc[0] / rp[0], 1.5)
        self.assertAlmostEqual(fcc[1] / rp[1], 1.5)

    def test_shared_fortran_applies_authoritative_totals(self):
        source = read("src/src_common/m_cec.f90")
        #> Every component is the pairing's own corrected total times a ratio,
        #> and the total is carried alongside so a reader can check the sum.
        self.assertIn("total = totals(k)", source)
        self.assertIn("flux%comp(k)%total = total", source)
        self.assertIn("total / (1d0 + 1d0 / descriptor%target(k)%r)", source)
        self.assertIn("total / (1d0 + descriptor%target(k)%r)", source)
        self.assertIn("flux%comp(cecTargetWater)%nonstomatal * h2o_to_ET", source)
        self.assertIn("flux%comp(cecTargetWater)%stomatal * h2o_to_ET", source)
        #> The band is a setting now, defaulting to the 0.2 Zahn et al. use.
        #> It was 0.05, which this file's own carbon_status() already
        #> contradicted - the reference model and the engine disagreed by
        #> a factor of four on which periods survive.
        self.assertIn("CecIsSingular(target%r, setup%singular_band)", source)
        self.assertNotIn("abs(descriptor%r_Fc + 1d0) < 0.05d0", source)
        self.assertIn("abs(r + 1d0) < band", source)

    def test_the_header_and_the_row_are_one_contract(self):
        """Both executables name the columns and emit the values from the same
        pair of helpers, walked in the same order.

        They used to be two hand-written lists in four files. A row that emits
        its values in a different order than its header names them does not
        fail to build; every column after the divergence is silently misread.
        """
        for header_path, writer_path in (
            ("src/src_rp/init_outfiles_rp.f90", "src/src_rp/write_out_full.f90"),
            ("src/src_fcc/init_out_files.f90", "src/src_fcc/write_out_full_fcc.f90"),
        ):
            header = read(header_path)
            writer = read(writer_path)
            self.assertIn("call CecPairs(cec_pairs, n_cec_pairs)", header)
            self.assertIn("call CecOutputColumns(", header)
            self.assertIn("call CecPairs(cec_pairs, n_cec_pairs)", writer)
            self.assertIn("call CecRowValues(", writer)
            #> Neither may spell a column name or a field of its own.
            self.assertNotIn("'E_cec,", header)
            self.assertNotIn("CECFlux%E_cec", writer)

        resolver = read("src/src_common/gas_slot_resolution.f90")
        #> Water first, carbon second, extras after - stated once, so the
        #> descriptor, the totals and the columns all agree on what target
        #> number two is.
        self.assertIn("slots(cecTargetWater) = pair%water_slot", resolver)
        self.assertIn("slots(cecTargetCarbon) = pair%carbon_slot", resolver)

    def test_fcc_cec_data_follow_covariances_like_the_header(self):
        writer = read("src/src_fcc/write_out_full_fcc.f90")
        self.assertGreater(writer.index("call CecRowValues("),
                           writer.index("lEx%cov_w(gas)"))


if __name__ == "__main__":
    unittest.main()
