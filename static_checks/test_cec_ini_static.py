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
            "EddyFlowProj%cec%min_o1_o2 = 0.20d0",
            "EddyFlowProj%cec%min_octant = 0.05d0",
            "EddyFlowProj%cec%min_valid = 0.90d0",
            "EddyFlowProj%cec%signal_strength = 70d0",
            "EddyFlowProj%cec%max_stationarity = 25d0",
            "EddyFlowProj%cec%max_gap_fill = 4",
        )
        for expected in expected_defaults:
            assert expected in source
        assert "NormalizeCecFraction" in source
        assert "NormalizeCecSignalStrength" in source
        assert "NormalizeCecMaxGapFill" in source
        assert "NormalizeCecStationarity" in source


    def test_cec_descriptor_uses_settings_instead_of_magic_thresholds(self):
        source = read("src/src_common/m_cec.f90")
        assert "active_setup%max_gap_fill" in source
        assert "active_setup%max_stationarity" in source
        assert "active_setup%min_valid" in source
        assert "active_setup%min_o1_o2" in source
        assert "active_setup%min_octant" in source
        assert "active_setup%signal_strength" in source
        assert "CecPassesHyperbolicThreshold" in source
        assert "call InterpolateShortCecGaps(w_prime, 4)" not in source
        assert "10 * descriptor%n_valid < 9 * nrow" not in source
        assert "descriptor%frac_O1 + descriptor%frac_O2 < 0.20d0" not in source
        assert "descriptor%frac_O1 < 0.05d0" not in source
        assert "stationarity_co2 > 25" not in source
        assert "stationarity_h2o > 25" not in source


    def test_cec_stationarity_threshold_accepts_default_relaxed_and_disabled_modes(self):
        source = read("src/src_common/m_cec.f90")

        assert "if (active_setup%max_stationarity > 0d0) then" in source
        assert "dble(stationarity_co2) > active_setup%max_stationarity" in source
        assert "dble(stationarity_h2o) > active_setup%max_stationarity" in source
        assert "setup%max_stationarity = 25d0" in source

        assert normalize_stationarity("25") == 25
        assert normalize_stationarity("50") == 50
        assert normalize_stationarity("0") == 0
        assert normalize_stationarity("-1") == 25

        def accepted(st_co2, st_h2o, threshold):
            return threshold <= 0 or (st_co2 <= threshold and st_h2o <= threshold)

        assert not accepted(30, 20, 25)
        assert accepted(30, 20, 50)
        assert accepted(999, 999, 0)


    def test_rp_passes_sample_level_signal_strength_when_available(self):
        source = read("src/src_rp/eddyflow-rp_main.f90")
        assert "UserSet(:, cec_co2_signal_col)" in source
        assert "UserSet(:, cec_h2o_signal_col)" in source
        assert "UserCol(j)%var == 'AGC'" in source
        assert "UserCol(j)%var == 'RSSI'" in source
        assert "SetLicorDiagnostics" in read("src/src_rp/set_licor_diagnostics.f90")

        cec_block = source[
            source.index("!> Extract CEC before spectral processing interpolates E2Primes.") :
            source.index("CECFlux%r_Fc_cec = CECDescriptor%r_Fc")
        ]
        assert "UserSet(:, cec_co2_signal_col)" in cec_block
        assert "UserSet(:, cec_h2o_signal_col)" in cec_block
        assert "Essentials%AGC" not in cec_block
        assert "Diag7200%AGC" not in cec_block


    def test_read_ex_record_parses_cec_tail_without_biomet_fields(self):
        source = read("src/src_common/read_ex_record.f90")

        # The CEC descriptor's own width is the named nCecFields (still 11); see
        # test_ex_record_layout_static.py, which pins the layout counts.
        assert "if (remaining_fields == nCecFields) then" in source
        assert "cec_line = dataline(1:len_trim(dataline))" in source
        assert "dataline = ''" in source
        assert "elseif (remaining_fields > nCecFields) then" in source
        assert "if (len_trim(cec_line) > 0) then" in source
        assert "strCharIndex(dataline, ',', remaining_fields - nCecFields)" in source

        cec_fields = [
            "0.5", "-0.25", "17900", "1400", "1000", "0.078", "0.056",
            "1", "1", "1", "1",
        ]
        biomet, cec = split_cec_tail(cec_fields)
        assert biomet == []
        assert cec == cec_fields

        biomet, cec = split_cec_tail(["TA", "PA"] + cec_fields)
        assert biomet == ["TA", "PA"]
        assert cec == cec_fields

        biomet, cec = split_cec_tail(cec_fields[:-1])
        assert biomet == cec_fields[:-1]
        assert cec == []


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
