"""FCC's binned (co)spectra are read by workers; everything else stays put.

The import loop depends on the order its files come in: it looks each up in
the essentials file reading forward, screens it, and adds it to sums whose
floating-point result depends on that order. So only the reading is split. A
worker parses its files and hands back their bins; the parent runs the loop
as it always did and takes each file's bins in file order instead of reading
the file. Byte-identical by construction - and modest: on 30 days of the Yatir
record the import took 3.1 s instead of 4.1 s, the parent's share setting the
pace.
"""

import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def read(rel):
    return (ROOT / rel).read_text(encoding="utf-8", errors="replace").replace("\r\n", "\n")


MAIN = read("src/src_fcc/eddyflow-fcc_main.f90")
PAR = read("src/src_fcc/fcc_parallel.f90")
POOL = read("src/src_common/batch_pool.f90")


def body(text, kind, name):
    i = text.index("%s %s(" % (kind, name))
    return text[i:text.index("end %s %s" % (kind, name), i)]


class OnlyTheReadingMoves(unittest.TestCase):

    def test_the_loop_takes_each_file_from_get_binned_file(self):
        i = MAIN.index("binned_loop: do")
        loop = MAIN[i:MAIN.index("end do binned_loop", i)]
        self.assertIn("call GetBinnedFile(fcount, BinnedFileList(fcount), BinSpec, BinCosp,", loop)
        self.assertNotIn("call ReadBinnedFile(", loop)
        self.assertIn("call RetrieveExVarsByTimestamp(", loop)

    def test_the_workers_start_before_the_loop_and_are_done_after_it(self):
        self.assertLess(MAIN.index("call StartBinnedSplit("), MAIN.index("binned_loop: do"))
        i = MAIN.index("end do binned_loop")
        self.assertIn("call FinishBinnedSplit()", MAIN[i:i + 60])

    def test_a_worker_branches_before_the_essentials_file_is_read(self):
        i = MAIN.index("if (BatchKind == 'sb') call RunBinnedReadWorker()")
        self.assertLess(i, MAIN.index("call CheckExFileVintageAt(AuxFile%ex)"))
        self.assertLess(i, MAIN.index("call InitExVars("))


class TheParentGetsWhatReadBinnedFileWouldHaveGiven(unittest.TestCase):

    G = body(PAR, "subroutine", "GetBinnedFile")

    def test_the_first_piece_is_read_here(self):
        self.assertIn("call PrepassChunk(1, s, e)", self.G)
        self.assertEqual(self.G.count("call ReadBinnedFile(InFile, BinSpec, BinCosp, nrow, nbins, skip)"), 2)

    def test_files_come_in_order_and_are_checked(self):
        self.assertIn("if (.not. more .or. fc /= fcount) &", self.G)

    def test_a_file_that_could_not_be_read_leaves_nbins_and_says_so(self):
        i = self.G.index("if (skip) then")
        block = self.G[i:self.G.index("end if", i)]
        self.assertIn("call ExceptionHandler(62)", block)
        self.assertLess(self.G.index("if (skip) then"), self.G.index("nbins = nb"))

    def test_bins_are_rebuilt_from_errspec(self):
        self.assertIn("Bins = ErrSpec", body(PAR, "subroutine", "ReadBins"))

    def test_a_dump_for_other_gases_is_refused(self):
        self.assertIn("if (any(mine /= theirs)) &", body(PAR, "subroutine", "OpenBinDump"))


class ThePoolCanBeConsumedInOrder(unittest.TestCase):

    def test_waiting_for_one_piece_keeps_the_others_busy(self):
        b = body(POOL, "subroutine", "WaitForBatchPiece")
        self.assertIn("call PollPrepassPool(PoolEff - 1, .true., progressed)", b)
        self.assertIn("if (k < NextChunk .and. .not. Running(k)) exit", b)

    def test_an_early_end_waits_for_the_running_and_starts_no_more(self):
        b = body(POOL, "subroutine", "CancelBatchPool")
        self.assertIn("NextChunk = NumChunks + 1", b)
        self.assertIn("PoolActive = .false.", b)
        f = body(PAR, "subroutine", "FinishBinnedSplit")
        self.assertIn("call CancelBatchPool()", f)
        self.assertIn("close(u, status = 'delete')", f)


if __name__ == "__main__":
    unittest.main()
