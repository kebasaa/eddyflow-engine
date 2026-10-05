"""A skipped period's FLUXNET row is as wide as the header it sits under.

The header names the custom variables the project declared when it was
written. A skipped period's row padded its custom block to NumUserVar + 1 -
the count of the period being skipped. Where every raw file brings its own
metadata, as LI-COR GHG archives do, that count can differ from the header's:
base_ghg_mixed_60's 02:00 period declared 13 custom variables against the
header's 16, and its row came out three columns short, every column after the
custom block under the wrong name. Nothing saw it until run.sh began keeping
RP's own FLUXNET file for check_columns.py (FCC's file used to replace it).

ReadExRecord stops at a not_enough_data row before the custom block, so FCC
never misread it; the file was simply not what its header says.

Part of the EddyFlow engine's static checks.
"""

from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[1]


def code(rel):
    text = (ROOT / rel).read_text(encoding="utf-8", errors="replace")
    return "\n".join(ln for ln in text.splitlines() if not ln.lstrip().startswith("!"))


INIT = code("src/src_rp/init_fluxnet_file_rp.f90")
SKIP = code("src/src_rp/write_out_fluxnet_only_biomet.f90")


class TheSkippedRowTakesTheHeadersCount(unittest.TestCase):

    def test_the_header_records_how_many_custom_variables_it_names(self):
        block = INIT[INIT.index("call AddDatum(csv_row, 'NUM_CUSTOM_VARS', separator)"):]
        block = block[:block.index("end if")]
        self.assertIn("nFluxnetCustomVars = max(0, NumUserVar)", block)

    def test_the_skipped_row_pads_to_that_count(self):
        row = SKIP[SKIP.index("subroutine WriteFluxnetOnlyBiometRow("):
                   SKIP.index("end subroutine WriteFluxnetOnlyBiometRow")]
        self.assertIn("do i = 1, nFluxnetCustomVars + 1", row)
        self.assertNotIn("do i = 1, NumUserVar + 1", row)


if __name__ == "__main__":
    unittest.main()
