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
        if not re.search(r"[a-z]", label):
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
        self.assertEqual(repaired_header2[custom_start], "flowrate_mean")
        self.assertEqual(repaired_header2[custom_start + 1], "custom_2_mean")
        self.assertEqual(repaired_header2[custom_start + ncustom], "E_cec")
        self.assertFalse(any(re.fullmatch(r"[-+]?(?:\d+|\d*\.\d+)", label) for label in custom_labels))


if __name__ == "__main__":
    unittest.main()
