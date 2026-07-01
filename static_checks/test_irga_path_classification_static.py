from pathlib import Path
import re
import unittest


ROOT = Path(__file__).resolve().parents[1]


def read(path):
    return (ROOT / path).read_text(encoding="utf-8")


def function_body(source, name):
    match = re.search(
        rf"{name}\(.*?end function {name}",
        source,
        flags=re.DOTALL | re.IGNORECASE,
    )
    if match is None:
        raise AssertionError(f"{name} not found")
    return match.group(0)


class IrgaPathClassificationStaticTests(unittest.TestCase):
    def test_shared_helper_classifies_csi_irgas(self):
        typedefs = read("src/src_common/m_typedef.f90")
        open_path = function_body(typedefs, "IsOpenPathIrgaModel")
        closed_path = function_body(typedefs, "IsClosedPathIrgaModel")

        for model in ("csi_ec150", "csi_irgason_irga"):
            self.assertIn(model, open_path)
            self.assertNotIn(model, closed_path)

        self.assertIn("csi_ec155", closed_path)
        self.assertNotIn("csi_ec155", open_path)

    def test_standard_and_dynamic_metadata_use_shared_path_classifier(self):
        for path in (
            "src/src_common/read_metadata_file.f90",
            "src/src_common/read_ex_record.f90",
            "src/src_rp/retrieve_dynamic_metadata.f90",
        ):
            self.assertIn("IrgaPathTypeFromModel", read(path))

    def test_dynamic_metadata_recognizes_non_licor_irga_firms(self):
        source = read("src/src_rp/retrieve_dynamic_metadata.f90")
        for firm in ("'csi_irga'", "'miro'", "'aerodyne'"):
            self.assertIn(firm, source)
        self.assertIn("case('csi_ec150', 'csi_ec155', 'csi_tga200a'", source)
        self.assertIn("case('licor', 'other_irga', 'csi_irga', 'miro', 'aerodyne')", source)

    def test_metadata_validation_treats_ec155_as_closed_path(self):
        source = read("src/src_common/metadata_file_validation.f90")
        self.assertNotIn(
            "case ('generic_open_path', 'csi_ec150', 'csi_ec155', "
            "'csi_irgason_irga')",
            source,
        )
        self.assertIn("case ('generic_closed_path', 'csi_ec155'", source)

    def test_legacy_model_prefixes_are_not_active_identifiers(self):
        legacy_prefix = "camp" + "bell_"
        for directory in ("src", "static_checks", "tests"):
            for path in (ROOT / directory).rglob("*"):
                if path.is_file() and path.suffix in {".f90", ".F", ".py", ".inc"}:
                    self.assertNotIn(legacy_prefix, path.read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
