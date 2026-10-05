"""The production pass split across workers gives what a single pass gives.

The end-to-end gate is tests/regression/check_parallel.sh with PAR_KIND=pr and
EDDYFLOW_PROD_PIECE_PERIODS=1, which cuts nearly every fixture at every
half-hour and diffs every output against -j 1. These checks pin the parts of
the design that a run only exercises when its data happen to need them:

* Every continuous output file is in the unit range the merge walks, so none
  is left behind in a worker's folder.
* Every per-period accumulator the loop updates is handed back and merged:
  the ok-period count, the storage cache, the PWB summary rows and donor
  tally.
* The PWB classifier lives in one place, m_pwb_stream, so the single pass and
  the parent's replay of a split one cannot drift apart.
* What cannot be split is refused: the Billesbach random stream, embedded
  mode, raw data from a shared link, an output file still to be opened.
* A cut needs the half-hour before it to begin with a raw file of its own,
  which is what makes a worker's file search land where a single pass's did.

Part of the EddyFlow engine's static checks.
"""

from pathlib import Path
import re
import unittest

ROOT = Path(__file__).resolve().parents[1]


def code(rel):
    text = (ROOT / rel).read_text(encoding="utf-8", errors="replace")
    return "\n".join(ln for ln in text.splitlines() if not ln.lstrip().startswith("!"))


MAIN = code("src/src_rp/eddyflow-rp_main.f90")
PROD = code("src/src_rp/production_parallel.f90")
STREAM = code("src/src_rp/pwb_stream.f90")
TLAG = code("src/src_rp/timelag_handle.f90")
UNITS = code("src/src_common/m_index_parameters.f90")


def unit_number(name):
    m = re.search(r"integer, parameter ::\s*" + name + r"\s*=\s*(\d+)", UNITS)
    return int(m.group(1))


def body(src, kind, name):
    i = src.index(kind + " " + name + "(")
    return src[i:src.index("end " + kind + " " + name, i)]


class EveryContinuousFileIsMerged(unittest.TestCase):

    def test_the_merge_walks_uqc_to_uflxnt(self):
        self.assertIn("integer, parameter :: FirstOutUnit = uqc", PROD)
        self.assertIn("integer, parameter :: LastOutUnit = uflxnt", PROD)

    def test_every_continuous_output_unit_is_in_that_range(self):
        lo, hi = unit_number("uqc"), unit_number("uflxnt")
        opened = set()
        for rel in ("init_outfiles_rp.f90", "init_user_outfiles.f90",
                    "init_fluxnet_file_rp.f90", "init_biomet_out.f90"):
            opened |= set(re.findall(r"open\((\w+), file", code("src/src_rp/" + rel)))
        self.assertTrue(opened)
        for u in opened:
            self.assertTrue(lo <= unit_number(u) <= hi, u)


class EveryAccumulatorIsHandedBack(unittest.TestCase):

    def test_the_worker_starts_them_afresh(self):
        begins = body(PROD, "subroutine", "ProdPieceBegins")
        for line in ("nOk = 0", "StorCacheN = 0", "PwbSummaryDonorCount = 0"):
            self.assertIn(line, begins)
        self.assertRegex(MAIN, r"call ProdPieceBegins\(NumberOfOkPeriods\)\s*\n\s*PwbTimelagN = 0")

    def test_and_the_parent_adds_them_up(self):
        merge = body(PROD, "subroutine", "MergeProdPieces")
        for line in ("nOk = nOk + wOk", "StorCacheN = StorCacheN + wStorN",
                     "pwbN = pwbN + wPwbN",
                     "PwbSummaryDonorCount = PwbSummaryDonorCount + wDonors"):
            self.assertIn(line, merge)

    def test_the_last_pieces_counts_reach_the_end_of_the_run(self):
        merge = body(PROD, "subroutine", "MergeProdPieces")
        self.assertIn("NumUserVar = wUser", merge)
        self.assertIn("nbVars = wBiomet", merge)


class ThePwbClassifierIsInOnePlace(unittest.TestCase):

    def test_timelag_handle_gathers_and_classifies_through_the_stream(self):
        self.assertIn("call PwbGatherEvidence(", TLAG)
        self.assertIn("call PwbClassifyPeriod(", TLAG)
        for verdict in ("'S1_optimal'", "'S2_optimal'", "'S3_carryforward'"):
            self.assertNotIn("reliability_class = " + verdict, TLAG)
            self.assertIn("reliability_class = " + verdict, STREAM)

    def test_evidence_does_not_read_the_streams_memory(self):
        gather = body(STREAM, "subroutine", "PwbGatherEvidence")
        for state in ("pwb_last_optimal", "pwb_has_previous"):
            self.assertNotIn(state, gather)

    def test_the_parent_replays_from_zeroed_row_lags_into_the_raw_arrays(self):
        replay = body(STREAM, "subroutine", "PwbReplayEvidence")
        self.assertLess(replay.index("RowLags = 0"), replay.index("call PwbClassifyPeriod("))
        self.assertIn("pwb_raw_ActTLag, pwb_raw_TLag, pwb_raw_DefTlagUsed", replay)

    def test_a_production_worker_takes_verdicts_from_its_lead_in_on(self):
        self.assertIn("if (BatchKind == 'pr' .and. pcount >= BatchSliceStart - 1) then", MAIN)


class WhatCannotBeSplitIsRefused(unittest.TestCase):

    def test_refusals(self):
        split = body(MAIN, "subroutine", "TryProdSplit")
        for cond in ("RUsetup%meth == 'billesbach_11'",
                     "EddyFlowProj%run_env == 'embedded'",
                     "RemoteActive()",
                     "EddyFlowProj%biomet_data == 'embedded' .and. initializeBiometOut",
                     "AddUserStatsHeader .and. userStatsOpen"):
            self.assertIn(cond, split)

    def test_a_worker_that_opens_a_file_mid_piece_stops(self):
        self.assertIn("if (op .and. .not. PieceOpen(u)) then",
                      body(PROD, "subroutine", "FinishProdWorker"))

    def test_only_a_serial_parent_decides(self):
        self.assertIn("ProdSplitArmed = BatchIndex == 0 .and. NumJobs /= 1", MAIN)


class CutsNeedACleanLeadIn(unittest.TestCase):

    def test_the_lead_in_starts_with_a_file_of_its_own(self):
        cut = body(PROD, "function", "ProdCutAllowed")
        self.assertIn("Series(q), Series(q + 1))) cycle", cut)
        self.assertIn("Series(q - 1), Series(q))) return", cut)

    def test_workers_skip_the_prepasses(self):
        self.assertIn("BatchKind /= 'pr' .and. BatchKind /= 'pd') then", MAIN)


class AStoppedRunTakesItsWorkersWithIt(unittest.TestCase):
    """The interface stops a run by terminating its job object, which takes
    every worker at once; a parent killed on its own - Task Manager, a
    command-line run - is noticed by each worker within a period."""

    def test_a_worker_checks_on_its_parent_even_while_skipping_ahead(self):
        start = MAIN.index("if (pcount > ProdHeadEnd .and. pcount < BatchSliceStart - 1) then")
        skip = MAIN[start:MAIN.index("cycle periods_loop", start)]
        self.assertIn("call StopIfParentGone()", skip)

    def test_and_while_repeating_the_drift_walk(self):
        start = MAIN.index("drift_loop: do")
        self.assertIn("call StopIfParentGone()", MAIN[start:start + 400])

    def test_an_orphaned_worker_removes_its_own_folder(self):
        pool = code("src/src_rp/prepass_parallel.f90")
        stop = body(pool, "subroutine", "StopIfParentGone")
        self.assertIn("BatchOwnOutDir", stop)
        self.assertLess(stop.index("if (op) close(u)"), stop.index("comm_rmdir"))
        self.assertIn("BatchOwnOutDir = WorkerMainOut",
                      body(PROD, "subroutine", "AdoptProdWorkerOutput"))

    def test_folders_of_killed_runs_are_removed_before_workers_start(self):
        split = body(MAIN, "subroutine", "TryProdSplit")
        self.assertLess(split.index("call RemoveStaleWorkerRoots()"),
                        split.index("call WriteProdContext(pcount)"))
        stale = body(PROD, "subroutine", "RemoveStaleWorkerRoots")
        self.assertIn("ProcessAlive(pid)", stale)
        self.assertIn("if (pid /= ProcessSelfId()) then", stale)

    def test_a_worker_never_writes_into_a_leftover_folder(self):
        adopt = body(PROD, "subroutine", "AdoptProdWorkerOutput")
        self.assertIn("comm_rmdir", adopt)

    def test_liveness_is_asked_the_same_way_on_every_system(self):
        for rel in ("src/src_os/win/m_process_os.f90", "src/src_os/posix/m_process_os.f90"):
            self.assertIn("logical function ProcessAlive(pid)", code(rel))


if __name__ == "__main__":
    unittest.main()
