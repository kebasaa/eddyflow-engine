"""The time-lag summaries list every gas the run had, not only those of its
last half-hour; and a split run counts every file a worker failed to download.

What the end of a pass writes - the time-lag optimisation and PWB summaries,
the windows the optimiser fits, the donors, the storage cleaning - goes by
E2Col, which describes the LAST period. A gas with no data in the final
half-hour was left out of the summary, as if never measured. Each processed
period now notes the gases it had (m_run_columns), and before the end of the
pre-pass and of the run a gas the last period lacked is put back as the last
period that had it left it. Split passes merge the notes piece by piece.

The shared-link summary - "N raw file(s) could not be downloaded" - was the
parent's own count only. Each process now marks a failed file in the staging
directory the run's processes share, and the parent counts the marks, so a
file failed by two processes counts once.
"""

import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def read(rel):
    return (ROOT / rel).read_text(encoding="utf-8", errors="replace").replace("\r\n", "\n")


COLS = read("src/src_rp/run_columns.f90")
MAIN = read("src/src_rp/eddyflow-rp_main.f90")
PROD = read("src/src_rp/production_parallel.f90")
PRE = read("src/src_rp/prepass_parallel.f90")
REMOTE = read("src/src_common/remote_source.f90")


def body(text, kind, name):
    i = text.index("%s %s(" % (kind, name))
    return text[i:text.index("end %s %s" % (kind, name), i)]


class EachPeriodNotesItsGases(unittest.TestCase):

    def test_a_noted_gas_keeps_the_last_period_that_had_it(self):
        b = body(COLS, "subroutine", "NoteRunColumns")
        self.assertIn("if (.not. E2Col(j)%present) cycle", b)
        self.assertIn("RunCol(j) = E2Col(j)", b)

    def test_only_gases_the_last_period_lacked_are_put_back(self):
        b = body(COLS, "subroutine", "UseRunColumns")
        self.assertIn("if (E2Col(j)%present .or. .not. RunColSeen(j)) cycle", b)

    def test_both_passes_note_every_processed_period(self):
        i = MAIN.index("call AddToTimelagOptDataset(TimelagOpt, TimelagOptSize, ton)")
        self.assertLess(MAIN.index("call NoteRunColumns()", i), MAIN.index("end do to_periods_loop"))
        i = MAIN.index("call WriteOutFluxnet(StDiff, DtDiff, STFlg, DTFlg)")
        self.assertIn("call NoteRunColumns()", MAIN[i:i + 120])

    def test_both_passes_start_afresh(self):
        self.assertLess(MAIN.index("call ResetRunColumns()"), MAIN.index("to_periods_loop: do"))
        i = MAIN.index("\n    periods_loop: do")
        self.assertIn("call ResetRunColumns()", MAIN[i - 300:i])

    def test_the_whole_run_s_gases_reach_the_end_of_each_pass(self):
        i = MAIN.index("!**** NOW STARTS TIME LAG OPT CALCULATIONS")
        self.assertLess(MAIN.index("call UseRunColumns()", i), MAIN.index("call PostProcessPwbTimelagCache()", i))
        j = MAIN.index("if (Test%stor_clean .and. StorCacheN > 0) call PostProcessStorClean()")
        self.assertLess(MAIN.rindex("call UseRunColumns()", 0, j), j)
        self.assertGreater(MAIN.rindex("call UseRunColumns()", 0, j), MAIN.index("call MergeProdPieces("))


class SplitPassesMergeTheNotesInOrder(unittest.TestCase):

    def test_a_production_worker_notes_only_its_piece(self):
        self.assertIn("call ResetRunColumns()", body(PROD, "subroutine", "ProdPieceBegins"))

    def test_production_dumps_carry_and_merge_them(self):
        self.assertIn("call WriteRunColumns(ud)", body(PROD, "subroutine", "FinishProdWorker"))
        self.assertIn("call MergeRunColumns(ud)", body(PROD, "subroutine", "MergeProdPieces"))

    def test_pre_pass_dumps_carry_and_merge_them(self):
        for w in ("WriteTlagBatchDump", "WritePwbBatchDump"):
            self.assertIn("call WriteRunColumns(u)", body(PRE, "subroutine", w), w)
        for m in ("MergeTlagBatchDumps", "MergePwbBatchDumps"):
            self.assertIn("call MergeRunColumns(u)", body(PRE, "subroutine", m), m)
        self.assertIn("BatchMagic = 'EDDYFLOW_PREPASS_08 '", PRE)

    def test_a_later_piece_overrides_an_earlier_one(self):
        b = body(COLS, "subroutine", "MergeRunColumns")
        self.assertIn("if (.not. seen(j)) cycle", b)
        self.assertIn("RunCol(j) = cols(j)", b)


class EveryFailedDownloadIsCountedOnce(unittest.TestCase):

    def test_a_failed_file_is_marked_where_every_process_sees_it(self):
        self.assertIn("call MarkFailed(Raw(e)%local)", REMOTE)
        self.assertIn("file = trim(local) // '.failed'", body(REMOTE, "subroutine", "MarkFailed"))

    def test_the_parent_counts_the_marks_and_a_worker_reports_no_total(self):
        b = body(REMOTE, "subroutine", "RemoteCleanup")
        self.assertIn("inquire(file = trim(Raw(e)%local) // '.failed', exist = failedMark)", b)
        self.assertIn("if (BatchIndex > 0) then\n            NumFailed = 0", b)
        self.assertLess(b.index("failedMark"), b.index("trim(comm_rmdir) // ' \"' // trim(StagingDir)"))


if __name__ == "__main__":
    unittest.main()
