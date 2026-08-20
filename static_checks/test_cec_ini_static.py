from pathlib import Path
import csv
import unittest


ROOT = Path(__file__).resolve().parents[1]


def read(path):
    return (ROOT / path).read_text(encoding="utf-8")


def normalize_stationarity(value, default=25):
    value = float(value)
    if value >= 0:
        return value
    return default


def split_cec_tail(fields):
    if len(fields) == 11:
        return [], fields
    if len(fields) > 11:
        return fields[:-11], fields[-11:]
    return fields, []


class CecIniStaticTests(unittest.TestCase):
    def test_cec_ini_keys_are_registered_as_numeric_project_tags(self):
        source = read("src/src_common/m_common_global_var.f90")
        for key in (
            "cec_h",
            "cec_min_o1_o2",
            "cec_min_octant",
            "cec_min_valid",
            "cec_signal_strength",
            "cec_max_gap_fill",
            "cec_max_stationarity",
        ):
            assert f"/ '{key}'" in source


    def test_project_parser_uses_project_tag_found_arrays(self):
        rp = read("src/src_rp/read_ini_rp.f90")
        fcc = read("src/src_fcc/read_ini_fcc.f90")
        assert "EPPrjNTagFound, EPPrjCTagFound" in rp
        assert "EPPrjNTagFound, EPPrjCTagFound" in fcc


    def test_cec_defaults_and_normalizers_are_project_level(self):
        source = read("src/src_common/write_processing_project_variables.f90")
        expected_defaults = (
            "EddyFlowProj%cec%h = 0d0",
            "EddyFlowProj%cec%singular_band = 0.2d0",
            "EddyFlowProj%cec%min_o1_o2 = 0.20d0",
            "EddyFlowProj%cec%min_octant = 0.05d0",
            "EddyFlowProj%cec%min_valid = 0.90d0",
            "EddyFlowProj%cec%signal_strength = 70d0",
            "EddyFlowProj%cec%max_stationarity = 25d0",
            "EddyFlowProj%cec%max_gap_fill = 4",
        )
        for expected in expected_defaults:
            assert expected in source
        assert "NormalizeCecPercent" in source
        assert "NormalizeCecBand" in source
        #> An occupancy limit is a percentage and nothing else. The old
        #> normalizer read anything in (0, 1] as a fraction instead, so an
        #> interface spin box labelled [%] and set to 0.5 came back as fifty.
        assert "NormalizeCecFraction" not in source
        assert "NormalizeCecSignalStrength" in source
        assert "NormalizeCecMaxGapFill" in source
        assert "NormalizeCecStationarity" in source


    def test_cec_descriptor_uses_settings_instead_of_magic_thresholds(self):
        source = read("src/src_common/m_cec.f90")
        #> `setup` rather than `active_setup`: the defaults are the parser's
        #> now and the caller always states them, so there is no second copy
        #> for the two to drift apart in.
        assert "setup%max_gap_fill" in source
        assert "setup%max_stationarity" in source
        assert "setup%min_valid" in source
        assert "setup%min_o1_o2" in source
        assert "setup%min_octant" in source
        assert "setup%signal_strength" in source
        assert "CecPassesHyperbolicThreshold" in source
        assert "setup%singular_band" in source
        assert "call InterpolateShortCecGaps(work(:, k), 4)" not in source
        assert "10 * descriptor%n_valid < 9 * nrow" not in source
        assert "descriptor%frac_O1 + descriptor%frac_O2 < 0.20d0" not in source
        assert "frac_O1 < 0.05d0" not in source
        assert "stationarity_carbon > 25" not in source
        assert "stationarity_water > 25" not in source


    def test_cec_stationarity_threshold_accepts_default_relaxed_and_disabled_modes(self):
        source = read("src/src_common/m_cec.f90")

        assert "if (setup%max_stationarity > 0d0) then" in source
        assert "dble(stationarity_carbon) > setup%max_stationarity" in source
        assert "dble(stationarity_water) > setup%max_stationarity" in source
        #> The default lives in the parser, and only there. The module used
        #> to carry a second copy for an optional argument nobody passed,
        #> which is a default that can drift from the one in force.
        assert "EddyFlowProj%cec%max_stationarity = 25d0" in read(
            "src/src_common/write_processing_project_variables.f90")
        assert "setup%max_stationarity = 25d0" not in source

        assert normalize_stationarity("25") == 25
        assert normalize_stationarity("50") == 50
        assert normalize_stationarity("0") == 0
        assert normalize_stationarity("-1") == 25

        def accepted(st_co2, st_h2o, threshold):
            return threshold <= 0 or (st_co2 <= threshold and st_h2o <= threshold)

        assert not accepted(30, 20, 25)
        assert accepted(30, 20, 50)
        assert accepted(999, 999, 0)


    def test_the_partition_screens_the_raw_series_before_detrending(self):
        """The screen has to run before the trend is fitted, not after.

        It used to edit the finished E2Primes, so the trend those fluctuations
        came from had been fitted straight through the samples being rejected.
        And it found one AGC or RSSI column by model substring and used it for
        both gases, which on a two-analyser site screened one analyser's gas
        with the other's diagnostic and on anything that is not a LI-7200 or a
        LI-7500 - a quantum cascade laser measuring carbonyl sulfide, say -
        screened nothing at all.
        """
        source = read("src/src_rp/eddyflow-rp_main.f90")
        assert "SetLicorDiagnostics" in read("src/src_rp/set_licor_diagnostics.f90")

        cec_block = source[
            source.index("!> ===== 6.3 CONDITIONAL EDDY COVARIANCE") :
            source.index("!> ===== 7. DETRENDING")
        ]
        assert "call BuildCecPrimes(" in cec_block
        assert "E2Set" in cec_block
        assert "UserSet" in cec_block
        #> Nothing here may reach for a diagnostic by model name any more.
        assert "li7200" not in cec_block
        assert "li7500" not in cec_block
        assert "Essentials%AGC" not in cec_block
        assert "Diag7200%AGC" not in cec_block

        #> Before the run's own detrending, and therefore before E2Set is
        #> deallocated a few lines below it.
        assert source.index("!> ===== 6.3 CONDITIONAL EDDY COVARIANCE") < \
            source.index("!> ===== 7. DETRENDING")

        #> Each gas is screened against its own analyser, and the two
        #> conventions are read in opposite directions.
        resolver = read("src/src_common/gas_slot_resolution.f90")
        assert "function CecSignalColumnFor" in resolver
        assert "UserCol(j)%instr%slot /= E2Col(gas_slot)%instr%slot" in resolver
        assert "function CecSignalIsRssi" in resolver
        assert "SwVerFromString('5.3.0')" in resolver
        cec = read("src/src_common/m_cec.f90")
        assert "if (signal_strength(i) < threshold) values(i) = error" in cec
        assert "if (signal_strength(i) > threshold) values(i) = error" in cec


    def test_read_ex_record_consumes_the_cec_block_before_biomet(self):
        source = read("src/src_common/read_ex_record.f90")

        #> The block is read forward off a count of pairings, like the
        #> hygrometer and analyser blocks before it, instead of being found by
        #> counting a fixed eleven fields back from the end of the row. That
        #> anchor made "nothing may be appended after the descriptor" a rule
        #> four files had to keep by hand, and it could not survive a block
        #> whose width depends on how many pairings a project declares.
        assert "nCecPairFixedFields = 9" in source
        assert "nCecTargetFields = 6" in source
        assert "n_cec_pairs > MaxNumCecPairs" in source
        assert "strCharIndex(dataline, ',', nCecPairFixedFields)" in source
        assert "strCharIndex(dataline, ',', nCecTargetFields)" in source
        assert "lEx%n_cec = n_cec_pairs" in source
        #> Nothing is anchored to the end of the row any more.
        assert "nCecFields" not in source
        assert "remaining_fields" not in source


    # test_ch_lae_project_cec_percent_defaults_normalize_to_expected_values was
    # removed here: it read data/CH-LAE_COS.eddyflow (deleted in 84fbb7a) and
    # asserted against a Python re-implementation of the rule, so it never
    # exercised the engine. The Fortran normalizers are covered by
    # test_cec_percent_style_project_values_are_normalized_in_fortran in
    # test_fcc_full_header_static.py.

    def test_ch_lae_full_output_regression_has_no_biomet_section_when_cec_tail_needed(self):
        sample = ROOT / "data" / "eddyflow_CH-LAE_COS_full_output_2026-07-02T093117_adv.csv"
        if not sample.exists():
            return

        with sample.open(encoding="utf-8", newline="") as handle:
            rows = list(csv.reader(handle))

        assert "biomet" not in ",".join(rows[0]).lower()
        header = rows[1]
        for name in (
            "E_cec",
            "Tr_cec",
            "E_cec_ET",
            "Tr_cec_ET",
            "r_ET_cec",
            "Reco_cec",
            "P_cec",
            "NEE_cec",
            "r_Fc_cec",
        ):
            assert name in header


if __name__ == "__main__":
    unittest.main()
