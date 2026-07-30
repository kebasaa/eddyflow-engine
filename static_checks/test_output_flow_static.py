from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]


def read(path):
    return (ROOT / path).read_text(encoding="utf-8")


class OutputFlowStaticTests(unittest.TestCase):
    def test_fcc_following_rp_run_does_not_publish_parent_full_output(self):
        source = read("src/src_rp/eddyflow-rp_main.f90")
        fcc_branch = source[
            source.index("EddyFlowProj%fcc_follows     = .true."):
            source.index("make_dataset_common         = .false.")
        ]
        self.assertIn("EddyFlowProj%out_full        = .false.", fcc_branch)
        self.assertNotIn("if (EddyFlowProj%col(gas4) == 0) EddyFlowProj%out_full", fcc_branch)
        self.assertIn("parent FLUXNET essentials file", fcc_branch)

    def test_fcc_spectral_assessment_quotes_main_output_directory(self):
        source = read("src/src_fcc/output_spectral_assessment_results.f90")
        self.assertIn(
            'CreateDir(\'"\' // Dir%main_out(1:len_trim(Dir%main_out)) // \'"\')',
            source,
        )

if __name__ == "__main__":
    unittest.main()
