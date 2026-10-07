"""FCC's flux computation, split across worker processes, writes what one
process would have.

After the first record that opens the output files - the head - the rest of
the essentials file is cut into pieces and handed to workers (m_batch_pool).
A worker takes the parent's spectral assessment from a context file instead
of repeating it, processes the head again into a folder of its own so its
files open as the parent's did, skips to its piece unread, and hands back
where its rows begin; the parent appends them in record order.
"""

import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def read(rel):
    return (ROOT / rel).read_text(encoding="utf-8", errors="replace").replace("\r\n", "\n")


MAIN = read("src/src_fcc/eddyflow-fcc_main.f90")
PAR = read("src/src_fcc/fcc_parallel.f90")


def body(text, kind, name):
    i = text.index("%s %s(" % (kind, name))
    return text[i:text.index("end %s %s" % (kind, name), i)]


class AWorkerTakesTheParentsAssessment(unittest.TestCase):

    def test_it_branches_before_the_assessment_and_skips_it(self):
        i = MAIN.index("if (BatchKind == 'fx') then\n        call FccWorkerStart(FccHead)\n        goto 100")
        self.assertLess(i, MAIN.index("SPECTRAL ASSESSMENT SECTION"))
        self.assertGreater(i, MAIN.index("call FullCospectraLength(FullFilelist(1)%path, nrow_full)"))

    def test_the_parent_keeps_the_state_where_the_flux_loop_starts(self):
        i = MAIN.index("100 continue")
        self.assertIn("if (BatchIndex == 0) call CaptureFccContext()", MAIN[i:i + 200])

    def test_the_context_carries_what_the_assessment_sets(self):
        w = body(PAR, "subroutine", "WriteFccContext")
        r = body(PAR, "subroutine", "FccWorkerStart")
        for name in ("Timestamp_FilePadding", "RegParR", "UnParR", "StParR",
                     "GasRateUse", "FileRateUse", "MultiRateSA"):
            self.assertIn(name, w, name)
            self.assertIn(name, r, name)
        self.assertIn("CtxProj, CtxFCCsetup, CtxRegPar, CtxStPar, CtxUnPar, CtxTFShape, CtxMassPar", w)
        self.assertIn("EddyFlowProj, FCCsetup, RegPar, StPar, UnPar, TFShape, MassPar", r)

    def test_a_worker_reports_as_itself_and_writes_into_its_own_folder(self):
        r = body(PAR, "subroutine", "FccWorkerStart")
        self.assertIn("EddyFlowProj%caller = ownCaller", r)
        self.assertIn("Dir%main_out = WorkerRoot(BatchParentPid, 'fx', BatchIndex)", r)


class ThePiecesAreCutAfterTheHead(unittest.TestCase):

    def test_the_split_is_tried_once_the_output_files_are_open(self):
        i = MAIN.index("call InitOutFiles(lEx)")
        self.assertIn("call TryFccFluxSplit(i, NumExRecords, FccParentEnd, FccWorkers)",
                      MAIN[i:i + 400])

    def test_embedded_mode_and_shared_links_are_refused(self):
        t = body(PAR, "subroutine", "TryFccFluxSplit")
        self.assertIn("if (EddyFlowProj%run_env == 'embedded') then", t)
        self.assertIn("if (RemoteFetchedAny()) then", t)

    def test_a_worker_reads_only_its_head_and_its_piece(self):
        i = MAIN.index("ex_loop: do i = 1, NumExRecords")
        loop = MAIN[i:MAIN.index("call ReadExRecord('', uex, -1,", i)]
        self.assertIn("if (i >= BatchSliceEnd) exit ex_loop", loop)
        self.assertIn("if (i /= FccHead .and. i < BatchSliceStart) then", loop)
        self.assertIn("read(uex, *, iostat = skip_status)", loop)
        self.assertIn("if (i == BatchSliceStart) call FccPieceBegins()", loop)
        self.assertIn("call StopIfParentGone()", loop)

    def test_the_parent_stops_at_the_first_cut_and_keeps_the_workers_busy(self):
        i = MAIN.index("ex_loop: do i = 1, NumExRecords")
        loop = MAIN[i:MAIN.index("call ReadExRecord('', uex, -1,", i)]
        self.assertIn("if (i >= FccParentEnd) exit ex_loop", loop)
        self.assertIn("call TopUpPrepassBatches()", loop)


class TheRowsAreAppendedInOrder(unittest.TestCase):

    def test_only_the_three_output_files_and_never_the_essentials(self):
        self.assertIn("integer, parameter :: OutUnits(NOut) = [uflx, umd, uflxnt]", PAR)
        self.assertNotIn("uex", PAR.replace("essentials file being read, uex,", ""))

    def test_a_worker_stops_before_the_end_of_the_run(self):
        i = MAIN.index("end do ex_loop")
        tail = MAIN[i:]
        w = tail.index("if (BatchKind == 'fx') then")
        self.assertLess(w, tail.index("PostProcessFluxDespiking"))
        self.assertLess(w, tail.index("keep_parent"))
        self.assertIn("call FinishFccWorker()\n        call FinishBatchWorker()\n        stop", tail[w:w + 200])

    def test_the_parent_waits_then_merges_before_closing_and_despiking(self):
        i = MAIN.index("end do ex_loop")
        tail = MAIN[i:]
        self.assertLess(tail.index("call WaitPrepassBatches('fx', FccWorkers)"),
                        tail.index("call MergeFccPieces()"))
        self.assertLess(tail.index("call MergeFccPieces()"),
                        tail.index("PostProcessFluxDespiking"))

    def test_despiking_periods_come_back_in_piece_order(self):
        m = body(PAR, "subroutine", "MergeFccPieces")
        self.assertIn("do k = 2, PrepassChunkCount()", m)
        self.assertIn("call StorePfdCache(rows(i)%date", m)
        self.assertIn("call AppendBytes(path, offset, parentPath(j))", m)
        self.assertIn("call RemoveWorkerRoots('fx')", m)


if __name__ == "__main__":
    unittest.main()
