import math
from pathlib import Path
import unittest


REJECTED = 0
NORMAL = 1
ALL_STOMATAL = 2
ALL_NONSTOMATAL = 3
SINGULAR = 4
H2O_TO_ET = 0.0648
ROOT = Path(__file__).resolve().parents[1]


def normalize_fraction(value, default):
    if 0 <= value <= 1:
        return value
    if 1 < value <= 100:
        return value / 100
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


def count_octants(w_values, c_values, q_values, h=0):
    valid = o1 = o2 = 0
    for w_value, c_value, q_value in zip(w_values, c_values, q_values):
        if None in (w_value, c_value, q_value):
            continue
        valid += 1
        if w_value > 0 and q_value > 0 and c_value > 0:
            if h <= 0 or (abs(w_value * q_value) >= h and abs(w_value * c_value) >= h):
                o1 += 1
        elif w_value > 0 and q_value > 0 and c_value < 0:
            if h <= 0 or (abs(w_value * q_value) >= h and abs(w_value * c_value) >= h):
                o2 += 1
    return valid, o1, o2


def apply_partition(status, ratio, total):
    if status == NORMAL:
        return total / (1 + 1 / ratio), total / (1 + ratio)
    if status == ALL_STOMATAL:
        return 0.0, total
    if status == ALL_NONSTOMATAL:
        return total, 0.0
    return math.nan, math.nan


def carbon_status(ratio):
    return SINGULAR if -1.2 < ratio < -0.8 else NORMAL


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
        self.assertEqual(normalize_fraction(0.90, 0.5), 0.90)
        self.assertEqual(normalize_fraction(90, 0.5), 0.90)
        self.assertEqual(normalize_fraction(-1, 0.5), 0.5)
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
        self.assertEqual(count_octants(w_values, c_values, q_values, h=0.1), (3, 0, 1))

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
        for ratio in (-1.2, -0.8, -0.79, -1.21):
            self.assertEqual(carbon_status(ratio), NORMAL)

    def test_fcc_scaling_preserves_partition(self):
        rp = apply_partition(NORMAL, 0.5, 3.0)
        fcc = apply_partition(NORMAL, 0.5, 4.5)
        self.assertAlmostEqual(fcc[0] / rp[0], 1.5)
        self.assertAlmostEqual(fcc[1] / rp[1], 1.5)

    def test_shared_fortran_applies_authoritative_totals(self):
        source = (ROOT / "src/src_common/m_cec.f90").read_text(encoding="utf-8")
        self.assertIn("flux%NEE_cec = Fc_total", source)
        self.assertIn("flux%E_cec_ET = flux%E_cec * h2o_to_ET", source)
        self.assertIn("flux%Tr_cec_ET = flux%Tr_cec * h2o_to_ET", source)
        self.assertIn("abs(descriptor%r_Fc + 1d0) < 0.05d0", source)

    def test_rp_and_fcc_use_h2o_flux_and_matching_output_order(self):
        expected_header = "E_cec,Tr_cec,E_cec_ET,Tr_cec_ET,r_ET_cec"
        expected_fields = [
            "CECFlux%E_cec",
            "CECFlux%Tr_cec",
            "CECFlux%E_cec_ET",
            "CECFlux%Tr_cec_ET",
            "CECFlux%r_ET_cec",
        ]
        programs = (
            ("src/src_rp/eddyflow-rp_main.f90", "src/src_rp/init_outfiles_rp.f90",
             "src/src_rp/write_out_full.f90"),
            ("src/src_fcc/eddyflow-fcc_main.f90", "src/src_fcc/init_out_files.f90",
             "src/src_fcc/write_out_full_fcc.f90"),
        )
        for main_path, header_path, writer_path in programs:
            main = (ROOT / main_path).read_text(encoding="utf-8")
            header = (ROOT / header_path).read_text(encoding="utf-8")
            writer = (ROOT / writer_path).read_text(encoding="utf-8")
            #> Water first, then CO2 - the order ApplyCecDescriptor expects.
            #> This read `Flux3%h2o, Flux3%co2` until the multi-gas refactor
            #> replaced the per-species scalars with Flux3%gas(slot), and then
            #> the literal slots with the resolved ones: CEC is defined on a
            #> CO2/water pair, but the pair is a species question, and slots
            #> five and six are CO2 and water by convention only.
            water = main.index("Flux3%gas(PrimaryWaterOutSlot())")
            carbon = main.index("Flux3%gas(PrimaryCarbonOutSlot())")
            self.assertLess(water, carbon,
                            "%s passes the CEC pair in the wrong order" % main_path)
            self.assertNotIn("Flux3%gas(h2o), Flux3%gas(co2)", main)
            self.assertIn(expected_header, header)
            positions = [writer.index(field) for field in expected_fields]
            self.assertEqual(positions, sorted(positions))

    def test_fcc_cec_data_follow_covariances_like_the_header(self):
        writer = (ROOT / "src/src_fcc/write_out_full_fcc.f90").read_text(encoding="utf-8")
        self.assertGreater(writer.index("CECFlux%E_cec"), writer.index("lEx%cov_w(gas)"))


if __name__ == "__main__":
    unittest.main()
