from pathlib import Path
import re
import unittest


ROOT = Path(__file__).resolve().parents[1]


def read(path):
    return (ROOT / path).read_text(encoding="utf-8")


def normalize_custom_labels(raw_labels, ncustom):
    labels = []
    for i in range(ncustom):
        label = raw_labels[i].replace("CUSTOM_", "").lower()
        if label == "flowrate":
            label = "flowrate_li7200_1_mean"
        elif not re.search(r"[a-z]", label):
            label = f"custom_{i + 1}_mean"
        elif "_mean" not in label:
            label = f"{label}_mean"
        labels.append(label)
    return labels


class FccFullHeaderStaticTests(unittest.TestCase):
    def test_fcc_custom_header_is_not_written_as_unbounded_csv_blob(self):
        source = read("src/src_fcc/init_out_files.f90")
        globals_source = read("src/src_fcc/m_fx_global_var_mod.f90")

        self.assertIn("character(32) :: UserVarHeader(MaxUserVar)", globals_source)
        self.assertNotIn("AddDatum(header2, UserVarHeader(1:len_trim(UserVarHeader))", source)
        self.assertEqual(source.count("custom_label = UserVarHeader(i)"), 2)
        self.assertEqual(source.count('write(custom_label, \'("custom_", i0, "_mean")\') i'), 2)
        self.assertEqual(source.count("if (i > 1) call AddDatum(header1, '', separator)"), 2)

    def test_fcc_flowrate_custom_headers_get_flowrate_units(self):
        source = read("src/src_fcc/init_out_files.f90")

        self.assertEqual(source.count("index(custom_label, 'flowrate_') == 1"), 2)
        self.assertEqual(source.count("call AddDatum(header3, '[m+3s-1]', separator)"), 2)

    def test_raw_flowrate_override_is_gas_scoped_and_instrument_specific(self):
        source = read("src/src_rp/eddyflow-rp_main.f90")

        override_block = source[source.index("replace instrument"):]
        self.assertIn("do i = co2, gas4", override_block)
        self.assertIn("UserCol(j)%var == 'flowrate'", override_block)
        self.assertIn("UserCol(j)%instr_name == E2Col(i)%instr_name", override_block)
        self.assertIn("E2Col(i)%instr%tube_f = UserStats%Mean(j)", override_block)
        self.assertNotIn("do i = 1, E2NumVar", override_block)

    def test_rp_flowrate_custom_headers_are_model_numbered_and_unitful(self):
        source = read("src/src_rp/init_outfiles_rp.f90")

        self.assertIn("user_header(NumUserVar)", source)
        self.assertIn("function FullOutputCustomLabel", source)
        self.assertIn("UserCol(j)%var == 'flowrate'", source)
        self.assertIn("case ('flowrate', 'co2', 'h2o', 'cell_t', 'int_p')", source)
        self.assertIn("clean_label = trim(var_token) // '_' // trim(model_token)", source)
        self.assertNotIn('write(user_header(j), \'("flowrate_", a, "_", i0, "_mean")\')', source)
        self.assertIn("user_unit(j) = '[m+3s-1]'", source)
        self.assertEqual(source.count("call AddDatum(header3, user_unit(var)"), 2)
        self.assertNotIn("usg(var)(1:len_trim(usg(var))) // 'mean'", source)

    def test_rp_multi_irga_custom_column_selection_is_consistent(self):
        for source_path, marker in (
            ("src/src_common/define_used_variables.f90", "NumUserVar = NumUserVar + 1"),
            ("src/src_rp/define_vars.f90", "NumUserVar = usr_cnt"),
            ("src/src_common/define_user_set.f90", "NumUserVar = jj"),
        ):
            source = read(source_path)
            self.assertIn("logical function IsCustomOutputColumn(col)", source)
            self.assertIn("if (col%useit) return", source)
            self.assertIn("if (len_trim(var) == 0) return", source)
            self.assertIn("'ignore', 'not_numeric', 'none'", source)
            self.assertIn("'agc', 'rssi'", source)
            self.assertIn(marker, source)

        define_vars = read("src/src_rp/define_vars.f90")
        pass3 = define_vars[
            define_vars.index("! Pass 3: collect user columns") :
            define_vars.index("\ncontains")
        ]
        self.assertIn("UserCol = NullCol", define_vars)
        self.assertIn("if (.not. IsCustomOutputColumn(LocCol(idx))) cycle", pass3)
        self.assertNotIn("case ('u','v','w','ts','sos','co2','h2o','ch4','n2o'", pass3)

    def test_rp_multi_irga_custom_headers_are_model_qualified(self):
        full_source = read("src/src_rp/init_outfiles_rp.f90")
        fluxnet_source = read("src/src_rp/init_fluxnet_file_rp.f90")

        for source in (full_source, fluxnet_source):
            self.assertIn("case ('flowrate', 'co2', 'h2o', 'cell_t', 'int_p')", source)
            self.assertIn("clean_label = trim(var_token) // '_' // trim(model_token)", source)
            self.assertIn("case ('a':'z', '0':'9', '_', '-')", source)
            self.assertNotIn("flow_ordinal", source)

        self.assertNotIn("previous_flow_model", full_source)
        self.assertNotIn("flowrate_, a", full_source)

    def test_metadata_blank_variable_guard_is_column_local(self):
        source = read("src/src_common/read_metadata_file.f90")

        self.assertIn("if (len_trim(LocCol(i)%var) == 0) LocCol(i)%var = 'ignore'", source)
        self.assertNotIn("if (len(LocCol(i)%var) == 0) LocCol%var = 'ignore'", source)

    def test_rp_full_headers_are_not_latin1_reencoded_before_utf8_write(self):
        source = read("src/src_rp/init_outfiles_rp.f90")
        full_output_block = source[: source.index("!>==========================================================================")]

        self.assertNotIn("latin1_to_utf8", full_output_block)
        self.assertIn("write(uflx, '(a)') header1(1:len_trim(header1) - 1)", full_output_block)
        self.assertIn("write(uflx, '(a)') header2(1:len_trim(header2) - 1)", full_output_block)
        self.assertIn("write(uflx, '(a)') header3(1:len_trim(header3) - 1)", full_output_block)
        self.assertIn("'[' // char(181) // 'mol", full_output_block)

    def test_fcc_custom_header_parser_bounds_and_sanitizes_labels(self):
        source = read("src/src_fcc/init_ex_vars.f90")

        self.assertIn("marker_custom = index(fluxnet_header, 'NUM_CUSTOM_VARS')", source)
        self.assertIn("marker_biomet = index(fluxnet_header, 'NUM_BIOMET_VARS')", source)
        self.assertIn("field_count < MaxUserVar", source)
        self.assertIn("UserVarHeader(field_count) = custom_label", source)
        self.assertIn("label_has_alpha", source)
        self.assertIn("index(custom_label, '_mean') == 0", source)

    def test_fcc_full_headers_are_not_latin1_reencoded_before_utf8_write(self):
        source = read("src/src_fcc/init_out_files.f90")
        full_output_block = source[: source.index("!> METADATA file")]

        self.assertNotIn("latin1_to_utf8", full_output_block)
        self.assertIn("write(uflx, '(a)') header1(1:len_trim(header1) - 1)", full_output_block)
        self.assertIn("write(uflx, '(a)') header2(1:len_trim(header2) - 1)", full_output_block)
        self.assertIn("write(uflx, '(a)') header3(1:len_trim(header3) - 1)", full_output_block)

    def test_bad_sample_tail_repairs_to_matching_header_and_data_counts(self):
        sample = ROOT / "eddyflow_CH-LAE_COS_full_output_2026-06-30T185123_adv.csv"
        if not sample.exists():
            self.skipTest("legacy full-output sample fixture is not present")
        rows = [line.split(",") for line in sample.read_text(encoding="utf-8").splitlines()[:4]]
        header1, header2, header3, data = rows

        self.assertEqual([len(row) for row in rows], [144, 183, 149, 149])

        custom_start = header1.index("custom_variables")
        ncustom = len(data) - custom_start - 9
        raw_custom_labels = header2[custom_start: header2.index("E_cec")]
        custom_labels = normalize_custom_labels(raw_custom_labels, ncustom)

        repaired_header1 = (
            header1[:custom_start]
            + ["custom_variables"]
            + [""] * (ncustom - 1)
            + ["conditional_eddy_covariance_(H2O)"]
            + [""] * 4
            + ["conditional_eddy_covariance_(CO2)"]
            + [""] * 3
        )
        repaired_header2 = (
            header2[:custom_start]
            + custom_labels
            + ["E_cec", "Tr_cec", "E_cec_ET", "Tr_cec_ET", "r_ET_cec"]
            + ["Reco_cec", "P_cec", "NEE_cec", "r_Fc_cec"]
        )
        repaired_header3 = header3[:custom_start] + ["--"] * ncustom + header3[-9:]

        self.assertEqual(len(repaired_header1), len(data))
        self.assertEqual(len(repaired_header2), len(data))
        self.assertEqual(len(repaired_header3), len(data))
        self.assertEqual(repaired_header2[custom_start], "flowrate_li7200_1_mean")
        self.assertEqual(repaired_header2[custom_start + 1], "custom_2_mean")
        self.assertEqual(repaired_header2[custom_start + ncustom], "E_cec")
        self.assertFalse(any(re.fullmatch(r"[-+]?(?:\d+|\d*\.\d+)", label) for label in custom_labels))


if __name__ == "__main__":
    unittest.main()
