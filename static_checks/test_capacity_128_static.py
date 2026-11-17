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

The ceiling then cost memory: a 4-gas run peaked at 246 MB committed where the
64-gas build took 132 MB. Two run-time sizes give it back, pinned below:

- The slot layout follows the project (SetSlotLayout): the cell slots, te and
  pe sit right after the last gas, so the per-period datasets (E2Set,
  E2Primes, SpecSet and every copy made of them) are pe wide, not E2NumVar.
  Every loop that used to stop at GHGNumVar now stops at lastGas, or treats
  the slots past it as absent - past the last gas they are cells now.
- The output buffers are LineLen long - LongOutstringLen's formula at the
  project's gas count - and on the heap. AddDatum scans its buffer's length,
  so that is a speed-up as much as a saving.

Part of the EddyFlow engine's static checks.
"""

from pathlib import Path
import re
import unittest

ROOT = Path(__file__).resolve().parents[1]


def read(*parts):
    return ROOT.joinpath(*parts).read_text(encoding="utf-8", errors="replace")


TYPEDEF = read("src", "src_common", "m_typedef.f90")
PRJVARS = read("src", "src_common", "write_processing_project_variables.f90")
BASIC = read("src", "src_common", "basic_stats.f90")
RPMAIN = read("src", "src_rp", "eddyflow-rp_main.f90")
SRC = ROOT / "src"
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



def body(text, name):
    start = text.index("subroutine " + name)
    return text[start:text.index("end subroutine " + name, start)]


class TheLayoutFollowsTheProject(unittest.TestCase):
    def test_layout_set_from_the_gas_count(self):
        layout = body(TYPEDEF, "SetSlotLayout")
        self.assertIn("lastGas = firstGas + max(4, min(ngas, MaxNumGases)) - 1", layout)
        self.assertIn("firstCell = lastGas + 1", layout)
        self.assertIn("pe  = lastCell + 2", layout)
        records = body(PRJVARS, "ReadMeasurementRecords")
        call = records.index("call SetSlotLayout(EddyFlowProj%gas_num)")
        #> Before the first record addresses a slot.
        self.assertLess(call, records.index("do i = 1, MaxNumGases"))
        self.assertGreater(call, records.index("EddyFlowProj%gas_num  = min(max("))

    def test_datasets_are_project_wide(self):
        for decl in ("allocate(E2Set(PeriodRecords, pe))",
                     "allocate(E2Primes(PeriodRecords, pe))",
                     "allocate(SpecSet(SpecRow, pe))"):
            self.assertIn(decl, RPMAIN)
        self.assertNotRegex(RPMAIN, r"allocate\((E2Set|E2Primes|SpecSet)\([^)]*E2NumVar")
        self.assertNotIn("E2Set(:, 1:GHGNumVar)", RPMAIN)
        self.assertNotIn("E2Primes(:, 1:GHGNumVar)", RPMAIN)

    def test_no_explicit_shape_dataset_dummy_at_the_ceiling(self):
        #> Sequence association: a pe-wide set handed to Set(N, E2NumVar) is
        #> read past its end. No compiler or runtime check sees it.
        for f in (SRC / "src_rp").glob("*.f90"):
            text = f.read_text(encoding="utf-8", errors="replace")
            self.assertNotRegex(text, r"::\s*Set\(\s*N\s*,\s*E2NumVar\s*\)", f.name)

    def test_covariance_into_the_leading_block(self):
        self.assertIn("Stats%Cov(1:ncol, 1:ncol), error)", BASIC)

    def test_cells_are_not_gases(self):
        spec = read("src", "src_rp", "spectral_analysis.f90")
        self.assertIn("where (SpecCol(u:lastGas)%present)", spec)
        self.assertNotIn("where (SpecCol(u:GHGNumVar)%present)", spec)
        self.assertIn("GasSlotPresent(lastGas + 1:) = .false.", RPMAIN)
        self.assertEqual(RPMAIN.count("Metadata%d, GasSlotPresent,"), 2)


class TheBuffersFollowTheProject(unittest.TestCase):
    def test_line_length_from_the_gas_count(self):
        self.assertIn("integer :: LineLen = LongOutstringLen", TYPEDEF)
        self.assertRegex(body(TYPEDEF, "SetSlotLayout"),
                         r"LineLen = 65536 \+ 2048 \* lastGas &\s*\n\s*"
                         r"\+ 32 \* \(\(lastGas - firstGas \+ 1\) \* \(lastGas - firstGas\)\) / 2")

    def test_buffers_on_the_heap(self):
        #> Fixed length, so assignment pads as before; allocatable, so the
        #> storage is the heap's - a stack copy at 128 gases would be megabytes
        #> on a Linux main thread that has eight.
        n = 0
        for f in SRC.rglob("*.f90"):
            text = f.read_text(encoding="utf-8", errors="replace")
            if f.name in ("m_typedef.f90", "m_fx_global_var_mod.f90"):
                continue
            self.assertNotIn("character(LongOutstringLen) ::", text, f.name)
            self.assertNotIn("character(LineLen) ::", text, f.name)
            n += text.count("character(LineLen), allocatable ::")
        self.assertGreaterEqual(n, 45)


if __name__ == "__main__":
    unittest.main()
