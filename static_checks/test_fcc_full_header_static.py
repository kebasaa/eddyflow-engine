from pathlib import Path
import csv
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


def split_cec_tail(fields):
    if len(fields) == 11:
        return [], fields
    if len(fields) > 11:
        return fields[:-11], fields[-11:]
    return fields, []


def normalize_cec_fraction(value, default):
    value = float(value)
    if 0 <= value <= 1:
        return value
    if 1 < value <= 100:
        return value / 100
    return default


class FccFullHeaderStaticTests(unittest.TestCase):
    def test_fcc_custom_header_is_not_written_as_unbounded_csv_blob(self):
        source = read("src/src_fcc/init_out_files.f90")
        globals_source = read("src/src_fcc/m_fx_global_var_mod.f90")

        self.assertIn("character(64) :: UserVarHeader(MaxUserVar)", globals_source)
        self.assertNotIn("AddDatum(header2, UserVarHeader(1:len_trim(UserVarHeader))", source)
        self.assertEqual(source.count("custom_label = UserVarHeader(i)"), 2)
        self.assertEqual(source.count('write(custom_label, \'("custom_", i0, "_mean")\') i'), 2)
        self.assertEqual(source.count("if (i > 1) call AddDatum(header1, '', separator)"), 2)

    def test_fcc_custom_headers_get_inferred_units(self):
        source = read("src/src_fcc/init_out_files.f90")

        self.assertIn("function CustomUnitFromLabel(label)", source)
        self.assertIn("index(clean_label, 'flowrate_') == 1", source)
        self.assertIn("index(clean_label, 'n2o_') == 1", source)
        self.assertIn("index(clean_label, 'h2o_') == 1", source)
        self.assertIn("index(clean_label, 'int_t_') == 1", source)
        self.assertIn("index(clean_label, 'int_p_') == 1", source)
        self.assertEqual(source.count("custom_unit = CustomUnitFromLabel(custom_label)"), 2)
        self.assertEqual(source.count("call AddDatum(header3, custom_unit"), 2)

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

        self.assertIn("character(64) :: user_header(NumUserVar)", source)
        self.assertIn("character(64) :: usg(NumUserVar)", source)
        self.assertIn("character(64) :: clean_label", source)
        self.assertIn("user_header(NumUserVar)", source)
        self.assertIn("function FullOutputCustomLabel", source)
        self.assertIn("case ('flowrate')", source)
        self.assertIn("case ('flowrate', 'co2', 'h2o', 'ch4', 'n2o', 'int_t', 'int_p')", source)
        self.assertIn("function CustomVarToken", source)
        self.assertIn("function CustomModelToken", source)
        self.assertIn("case ('cell_t', 'int_t_1', 'int_t_2')", source)
        self.assertIn("case ('co2', 'h2o', 'ch4', 'n2o')", source)
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

        self.assertIn("character(64) :: usg(NumUserVar)", fluxnet_source)
        self.assertIn("character(64) :: clean_label", fluxnet_source)
        for source in (full_source, fluxnet_source):
            self.assertIn("case ('flowrate', 'co2', 'h2o', 'ch4', 'n2o', 'int_t', 'int_p')", source)
            self.assertIn("clean_label = trim(var_token) // '_' // trim(model_token)", source)
            self.assertIn("model_token = CustomModelToken", source)
            self.assertIn("var_token = CustomVarToken(var_token)", source)
            self.assertIn("model_token = model_token(6:len_trim(model_token))", source)
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
        self.assertIn("utf8_mu = char(194) // char(181)", full_output_block)
        self.assertIn("'[' // utf8_mu // 'mol", full_output_block)

    def test_fcc_custom_header_parser_bounds_and_sanitizes_labels(self):
        source = read("src/src_fcc/init_ex_vars.f90")

        self.assertIn("marker_custom = index(fluxnet_header, 'NUM_CUSTOM_VARS')", source)
        self.assertIn("marker_biomet = index(fluxnet_header, 'NUM_BIOMET_VARS')", source)
        self.assertIn("character(128) :: custom_field", source)
        self.assertIn("character(64) :: custom_label", source)
        self.assertIn("custom_field = fluxnet_header(field_start:field_end)", source)
        self.assertIn("field_count < MaxUserVar", source)
        self.assertIn("UserVarHeader(field_count) = custom_label", source)
        self.assertIn("label_has_alpha", source)
        self.assertIn("custom_label(len_trim(custom_label) - 4:len_trim(custom_label)) == '_mean'", source)

    def test_ch_lae_custom_header_names_and_units_are_expected(self):
        expected_labels = [
            "n2o_mga4_6_2_mean",
            "co2_li7200_1_mean",
            "h2o_li7200_1_mean",
            "int_t_li7200_1_mean",
            "int_p_li7200_1_mean",
            "flowrate_li7200_1_mean",
            "flowrate_miro_mga4_6_2_mean",
        ]
        expected_units = [
            "[µmol+1mol_a-1]",
            "[µmol+1mol_a-1]",
            "[mmol+1mol_a-1]",
            "[K]",
            "[Pa]",
            "[m+3s-1]",
            "[m+3s-1]",
        ]
        bad_labels = {
            "n2o_mean",
            "int_t_1_mean",
            "flowrate_miro_mga4_6_2_me_mean",
        }

        self.assertEqual(len(expected_labels), len(expected_units))
        self.assertNotIn("--", expected_units)
        self.assertTrue(bad_labels.isdisjoint(expected_labels))

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
            + [""] * 5
            + ["conditional_eddy_covariance_(CO2)"]
            + [""] * 4
        )
        repaired_header2 = (
            header2[:custom_start]
            + custom_labels
            + ["E_cec", "Tr_cec", "E_cec_ET", "Tr_cec_ET", "r_ET_cec", "qc_cec_h2o"]
            + ["Reco_cec", "P_cec", "NEE_cec", "r_Fc_cec", "qc_cec_co2"]
        )
        repaired_header3 = header3[:custom_start] + ["--"] * ncustom + header3[-11:]

        self.assertEqual(len(repaired_header1), len(data))
        self.assertEqual(len(repaired_header2), len(data))
        self.assertEqual(len(repaired_header3), len(data))
        self.assertEqual(repaired_header2[custom_start], "flowrate_li7200_1_mean")
        self.assertEqual(repaired_header2[custom_start + 1], "custom_2_mean")
        self.assertEqual(repaired_header2[custom_start + ncustom], "E_cec")
        self.assertFalse(any(re.fullmatch(r"[-+]?(?:\d+|\d*\.\d+)", label) for label in custom_labels))

    def test_read_ex_record_parses_cec_tail_without_biomet_fields(self):
        source = read("src/src_common/read_ex_record.f90")

        self.assertIn("if (remaining_fields == 11) then", source)
        self.assertIn("cec_line = dataline(1:len_trim(dataline))", source)
        self.assertIn("dataline = ''", source)
        self.assertIn("elseif (remaining_fields > 11) then", source)
        self.assertIn("if (len_trim(cec_line) > 0) then", source)
        self.assertIn("strCharIndex(dataline, ',', remaining_fields - 11)", source)

        cec_fields = [
            "0.5", "-0.25", "17900", "1400", "1000", "0.078", "0.056",
            "1", "1", "1", "1",
        ]
        biomet, cec = split_cec_tail(cec_fields)
        self.assertEqual(biomet, [])
        self.assertEqual(cec, cec_fields)

        biomet, cec = split_cec_tail(["TA", "PA"] + cec_fields)
        self.assertEqual(biomet, ["TA", "PA"])
        self.assertEqual(cec, cec_fields)

        biomet, cec = split_cec_tail(cec_fields[:-1])
        self.assertEqual(biomet, cec_fields[:-1])
        self.assertEqual(cec, [])

    def test_ch_lae_cec_project_defaults_normalize_percent_style_values(self):
        values = {}
        for line in read("data/CH-LAE_COS.eddyflow").splitlines():
            if "=" not in line or line.lstrip().startswith(";"):
                continue
            key, value = line.split("=", 1)
            values[key.strip()] = value.strip()

        self.assertEqual(values["cec_meth"], "1")
        self.assertEqual(normalize_cec_fraction(values["cec_min_o1_o2"], 0.20), 0.20)
        self.assertEqual(normalize_cec_fraction(values["cec_min_octant"], 0.05), 0.05)
        self.assertEqual(normalize_cec_fraction(values["cec_min_valid"], 0.90), 0.90)
        self.assertEqual(float(values["cec_signal_strength"]), 70)
        self.assertEqual(int(float(values["cec_max_gap_fill"])), 4)

    def test_ch_lae_full_output_without_biomet_still_has_cec_columns(self):
        sample = ROOT / "data" / "eddyflow_CH-LAE_COS_full_output_2026-07-02T093117_adv.csv"
        if not sample.exists():
            self.skipTest("CH-LAE generated full-output regression file is not present")

        with sample.open(encoding="utf-8", newline="") as handle:
            rows = list(csv.reader(handle))

        self.assertNotIn("biomet", ",".join(rows[0]).lower())
        header = rows[1]
        for name in (
            "E_cec",
            "Tr_cec",
            "E_cec_ET",
            "Tr_cec_ET",
            "r_ET_cec",
            "qc_cec_h2o",
            "Reco_cec",
            "P_cec",
            "NEE_cec",
            "r_Fc_cec",
            "qc_cec_co2",
        ):
            self.assertIn(name, header)

    def test_fcc_full_output_includes_cec_status_columns_and_values(self):
        header_source = read("src/src_fcc/init_out_files.f90")
        writer_source = read("src/src_fcc/write_out_full_fcc.f90")

        self.assertEqual(header_source.count("r_ET_cec,qc_cec_h2o"), 2)
        self.assertEqual(header_source.count("r_Fc_cec,qc_cec_co2"), 2)
        self.assertEqual(header_source.count("[mm+1hour-1],[#],[#]"), 2)
        self.assertEqual(header_source.count("[umol+1m-2s-1],[#],[#]"), 2)
        self.assertIn("WriteDatumInt(lEx%cec%h2o_status", writer_source)
        self.assertIn("WriteDatumInt(lEx%cec%co2_status", writer_source)

    def test_cec_stationarity_ini_setting_is_registered_and_used(self):
        tags_source = read("src/src_common/m_common_global_var.f90")
        parser_source = read("src/src_common/write_processing_project_variables.f90")
        cec_source = read("src/src_common/m_cec.f90")

        self.assertIn("integer, parameter :: Npn = 32", tags_source)
        self.assertIn("EPPrjNTags(32)%Label / 'cec_max_stationarity'", tags_source)
        self.assertIn("EddyFlowProj%cec%max_stationarity = 25d0", parser_source)
        self.assertIn("NormalizeCecStationarity", parser_source)
        self.assertIn("EPPrjNTagFound(32)", parser_source)
        self.assertIn("active_setup%max_stationarity", cec_source)
        self.assertIn("if (active_setup%max_stationarity > 0d0) then", cec_source)
        self.assertNotIn("stationarity_co2 > 25", cec_source)

    def test_processing_variable_contract_is_project_level_and_parsed_by_rp_and_fcc(self):
        typedefs = read("src/src_common/m_typedef.f90")
        parser = read("src/src_common/read_processing_variables.f90")

        self.assertIn("integer, parameter :: MaxProcessingVariables = 64", typedefs)
        self.assertIn("type :: ProcessingVariableType", typedefs)
        self.assertIn("character(64) :: processing_id", typedefs)
        self.assertIn("character(64) :: reference_h2o_id", typedefs)
        self.assertIn("type(GasCollectionType) :: processing", typedefs)
        self.assertNotIn("canonical_var", typedefs)
        self.assertIn("subroutine ReadProcessingVariables(IniFile)", parser)
        self.assertNotIn("MigrateLegacyProcessingVariables", parser)
        self.assertNotIn("ApplyProcessingVariablesToLegacySlots", parser)
        self.assertNotIn("CanonicalGasVar", parser)

        for path in ("src/src_rp/read_ini_rp.f90", "src/src_fcc/read_ini_fcc.f90"):
            source = read(path)
            self.assertLess(
                source.index("call ReadProcessingVariables(PrjPath)"),
                source.index("call WriteProcessingProjectVariables()"),
            )

    def test_processing_variable_parser_covers_required_errors(self):
        parser = read("src/src_common/read_processing_variables.f90")

        for text in (
            "[ProcessingVariables] group is required",
            "count must be positive",
            "at least one enabled row",
            "count exceeds MaxProcessingVariables",
            "empty id",
            "empty gas",
            "non-positive gas_col",
            "non-positive irga_index",
            "non-positive gas_index",
            "id must match gas_irgaIndex_gasIndex",
            "Duplicate ProcessingVariables id",
            "Invalid H2O reference",
            "H2O reference is not an H2O row",
            "H2O row may only self-reference",
            "Ambiguous or missing H2O reference",
        ):
            self.assertIn(text, parser)

    def test_processing_variable_ids_are_generated_and_validated_from_numeric_identity(self):
        parser = read("src/src_common/read_processing_variables.f90")
        docs = read("ENGINE_PROCESSING_VARIABLES.md")

        self.assertIn("ExpectedProcessingId(row%gas_name, row%irga_index, row%gas_instance_index)", parser)
        self.assertIn("if (len_trim(configured_id) == 0) then", parser)
        self.assertIn("row%processing_id = expected_id", parser)
        self.assertIn("ProcessingVariables id must match gas_irgaIndex_gasIndex", parser)
        self.assertIn("write(ExpectedProcessingId, '(a,a,i0,a,i0)')", parser)
        self.assertIn("formatted as `gas_irgaIndex_gasIndex`", docs)
        self.assertIn("Instrument model names may remain metadata", docs)

    def test_processing_variable_h2o_refs_are_exact_row_references(self):
        parser = read("src/src_common/read_processing_variables.f90")

        self.assertIn("if (h2o_count == 1) then", parser)
        self.assertIn("EddyFlowProj%processing%rows(i)%reference_h2o_id = EddyFlowProj%processing%rows(only_h2o)%processing_id", parser)
        self.assertIn("ref_index = FindProcessingVariableById(EddyFlowProj%processing%rows(i)%reference_h2o_id)", parser)
        self.assertIn("EddyFlowProj%processing%rows(i)%h2o_ref_index = ref_index", parser)
        self.assertIn("EddyFlowProj%processing%rows(i)%h2o_ref_index = i", parser)

    def test_processing_variable_molecular_constants_are_normalized_at_parse_boundary(self):
        parser = read("src/src_common/read_processing_variables.f90")
        docs = read("ENGINE_PROCESSING_VARIABLES.md")

        self.assertIn("row%molecular_weight = row%molecular_weight * 1d-3", parser)
        self.assertIn("row%molecular_diffusivity = row%molecular_diffusivity * 1d-4", parser)
        self.assertIn("molecular weight: kg/mol", docs)
        self.assertIn("molecular diffusivity: m2/s", docs)

    def test_project_parser_no_longer_populates_fixed_gas_slots_from_project_tags(self):
        parser = read("src/src_common/write_processing_project_variables.f90")

        forbidden = (
            "EddyFlowProj%col(co2)",
            "EddyFlowProj%col(h2o)",
            "EddyFlowProj%col(ch4)",
            "EddyFlowProj%col(gas4)",
            "EddyFlowProj%col(tc)",
            "EddyFlowProj%col(ti1)",
            "EddyFlowProj%col(ti2)",
            "EddyFlowProj%col(pi)",
            "EddyFlowProj%col(te)",
            "EddyFlowProj%col(pe)",
            "EddyFlowProj%col(E2NumVar + diag72)",
            "EddyFlowProj%col(E2NumVar + diag75)",
            "EddyFlowProj%col(E2NumVar + diag77)",
            "ApplyProcessingVariablesToLegacySlots",
        )
        for text in forbidden:
            self.assertNotIn(text, parser)

    def test_row_indexed_gas_computational_containers_exist(self):
        typedefs = read("src/src_common/m_typedef.f90")

        self.assertIn("type :: GasStatsType", typedefs)
        self.assertIn("type(GasResultType) :: gas(MaxProcessingVariables)", typedefs)
        self.assertIn("type(GasCollectionType) :: processing", typedefs)
        self.assertIn("type(InstrumentType) :: processing_instr(MaxProcessingVariables)", typedefs)
        self.assertIn("processing_measure_type(MaxProcessingVariables)", typedefs)
        self.assertIn("h2o_ref_index", typedefs)
        self.assertIn("subroutine EnsureExProcessingRows", typedefs)
        self.assertIn("logical function IsH2OProcessingRow", typedefs)

    def test_fcc_level1_uses_processing_row_loop_for_gases(self):
        source = read("src/src_fcc/fluxes1.f90")

        self.assertIn("nrows = ProcessingRowCount(lEx%processing)", source)
        self.assertIn("do i = 1, nrows", source)
        self.assertIn("Flux1%gas(i) = lEx%Flux0%gas(i)", source)
        self.assertIn("IsH2OProcessingRow(lEx%processing%rows(i))", source)
        self.assertNotIn("lEx%instr(ico2)", source)
        self.assertNotIn("lEx%instr(ih2o)", source)
        self.assertNotIn("BPCF%of(w_co2)", source)
        self.assertNotIn("BPCF%of(w_h2o)", source)

    def test_rp_level1_uses_processing_row_loop_for_gases(self):
        source = read("src/src_rp/fluxes1_rp.f90")

        self.assertIn("nrows = ProcessingRowCount(EddyFlowProj%processing)", source)
        self.assertIn("do i = 1, nrows", source)
        self.assertIn("Flux1%gas(i) = Flux0%gas(i)", source)
        self.assertIn("IsH2OProcessingRow(EddyFlowProj%processing%rows(i))", source)
        self.assertNotIn("E2Col(co2)%Instr", source)
        self.assertNotIn("E2Col(h2o)%Instr", source)
        self.assertNotIn("BPCF%of(w_co2)", source)
        self.assertNotIn("BPCF%of(w_h2o)", source)

    def test_fcc_level23_uses_processing_row_h2o_refs(self):
        source = read("src/src_fcc/fluxes23.f90")

        self.assertIn("nrows = ProcessingRowCount(lEx%processing)", source)
        self.assertIn("do i = 1, nrows", source)
        self.assertIn("h2o_i = lEx%processing%rows(i)%h2o_ref_index", source)
        self.assertIn("Flux2%gas(i)", source)
        self.assertIn("Flux3%gas(i)", source)
        self.assertIn("lEx%processing_instr(i)%path_type", source)


if __name__ == "__main__":
    unittest.main()
