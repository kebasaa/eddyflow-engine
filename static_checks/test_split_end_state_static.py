"""The end of a split run reads the run's last period, not the parent's.

After a production pass split across workers, the parent runs the end of the
run once - the time-lag optimiser, the PWB summary, the storage cleaning, the
datasets. Some of what they read is per-period state the last period leaves
behind: which gases were present, their instruments (E2Col). The parent held
its own piece's last period, so a Yatir run split four ways left a gas out of
its PWB summary that the serial run listed. The last piece now hands E2Col
back with its custom and biomet variable counts, and the parent adopts it.

One count the end reports was also lost: peaks outside the lag window are
counted while detecting, which in a split PWB run happens in the detection
workers. Each counts its own piece and hands the counts back with its
evidence; the parent adds them.
"""

import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def read(rel):
    return (ROOT / rel).read_text(encoding="utf-8", errors="replace").replace("\r\n", "\n")


PROD = read("src/src_rp/production_parallel.f90")
MAIN = read("src/src_rp/eddyflow-rp_main.f90")
PWB = read("src/src_rp/pwb_timelag_handle.f90")


def body(text, kind, name):
    i = text.index("%s %s(" % (kind, name))
    return text[i:text.index("end %s %s" % (kind, name), i)]


class TheLastPieceHandsBackItsVariables(unittest.TestCase):

    def test_a_production_worker_writes_e2col(self):
        b = body(PROD, "subroutine", "FinishProdWorker")
        self.assertLess(b.index("write(ud) PwbSummaryDonorCount"), b.index("write(ud) E2Col"))
        self.assertLess(b.index("write(ud) E2Col"), b.rindex("write(ud) DumpMagic"))

    def test_the_parent_takes_the_last_piece_s(self):
        b = body(PROD, "subroutine", "MergeProdPieces")
        self.assertIn("read(ud) wE2Col", b)
        i = b.index("if (k == nPieces) then")
        self.assertIn("E2Col = wE2Col", b[i:b.index("end if", i)])


class OutsideWindowPeaksAreCountedOnce(unittest.TestCase):

    def test_a_detection_worker_counts_only_its_piece(self):
        self.assertIn("if (pcount == BatchSliceStart .and. BatchKind == 'pd') &\n"
                      "                call ResetPwbOutsideWindow()", MAIN)

    def test_it_hands_the_counts_back_and_the_parent_adds_them(self):
        fin = body(PROD, "subroutine", "FinishEvidenceWorker")
        self.assertLess(fin.index("write(EvidenceUnit) .false."), fin.index("write(EvidenceUnit) counts"))
        rep = body(PROD, "subroutine", "ReplayProdEvidence")
        self.assertIn("call AddPwbOutsideWindow(wOutside)", rep)

    def test_the_helpers_touch_only_that_counter(self):
        self.assertIn("pwb_outside_window = 0", body(PWB, "subroutine", "ResetPwbOutsideWindow"))
        self.assertIn("pwb_outside_window = pwb_outside_window + counts",
                      body(PWB, "subroutine", "AddPwbOutsideWindow"))


if __name__ == "__main__":
    unittest.main()
