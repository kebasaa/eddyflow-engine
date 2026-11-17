"""The gas ceiling of 128, and what kept it from costing a small project.

MaxNumGases went from 64 to 128 (MaxNumCol 200 -> 320, MaxUserVar 30 -> 64)
for PTR-MS datasets, innFLUX's use case (Striednig et al. 2020). The first
build with it was 3.5 times slower on a 4-gas project and overflowed the
Windows stack. Three changes made it free, each pinned here because undoing
one would bring that back without failing anything else:

- CovarianceMatrixNoError skips slots with no data and copies the symmetric
  half, so its cost follows the variables a project has, not the ceiling.
- The (co)spectra and ogive writers build rows in a row-sized buffer
  (SpecLineLen): AddDatum costs the buffer's whole length per call, and those
  files write a row per frequency.
- The engine reserves 64 MB of stack on Windows.

LongOutstringLen is derived from the ceilings because the FLUXNET header
grows with the square of the gas count (one covariance per pair).

Part of the EddyFlow engine's static checks.
"""

from pathlib import Path
import re
import unittest

ROOT = Path(__file__).resolve().parents[1]


def read(*parts):
    return ROOT.joinpath(*parts).read_text(encoding="utf-8", errors="replace")


TYPEDEF = read("src", "src_common", "m_typedef.f90")
STATS = read("src", "src_common", "stats_operator_no_error.f90")
SPEC = read("src", "src_rp", "spectral_analysis.f90")
CMAKE = read("CMakeLists.txt")
GEN_MD = read("prj", "gen_metadata_tags.py")
GUI_DEFS = ROOT.parent / "eddyflow-gui" / "src" / "defs.h"


class TheCeilings(unittest.TestCase):
    def test_values(self):
        self.assertIn("integer, parameter :: MaxNumGases = 128", TYPEDEF)
        self.assertIn("integer, parameter :: MaxNumCol = 320", TYPEDEF)
        self.assertIn("integer, parameter :: MaxUserVar = 64", TYPEDEF)
        self.assertIn("NUM_COLS = 320", GEN_MD)

    def test_output_line_grows_with_the_pairs(self):
        self.assertRegex(TYPEDEF, r"LongOutstringLen = 65536 \+ 2048 \* GHGNumVar &\s*\n\s*"
                                  r"\+ 32 \* \(MaxNumGases \* \(MaxNumGases - 1\)\) / 2")

    @unittest.skipUnless(GUI_DEFS.is_file(), "GUI checkout not beside the engine")
    def test_gui_agrees(self):
        self.assertIn("const int MAX_GASES = 128;", GUI_DEFS.read_text(encoding="utf-8"))


class WhatKeptItFree(unittest.TestCase):
    def test_covariance_skips_empty_slots_and_mirrors(self):
        body = STATS[STATS.index("subroutine CovarianceMatrixNoError"):]
        body = body[:body.index("end subroutine CovarianceMatrixNoError")]
        self.assertIn("hasData(i) = any(Set(:, i) /= err_float)", body)
        self.assertIn("Cov(i, j) = Cov(j, i)", body)
        self.assertIn("if (.not. (hasData(i) .and. hasData(j))) then", body)

    def test_spectra_rows_are_row_sized(self):
        self.assertIn("SpecLineLen = DatumLen * (2 * GHGNumVar + 8)", TYPEDEF)
        self.assertNotIn("character(LongOutstringLen) :: dataline", SPEC)
        self.assertGreaterEqual(SPEC.count("character(SpecLineLen) :: dataline"), 6)

    def test_windows_stack(self):
        self.assertIn('target_link_options(${target} PRIVATE "-Wl,--stack,67108864")', CMAKE)
        self.assertRegex(CMAKE, r"if\(WIN32\)\s*\n\s*target_link_options\(\$\{target\} PRIVATE \"-Wl,--stack")


if __name__ == "__main__":
    unittest.main()
