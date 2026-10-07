"""A pre-pass worker runs only the pre-pass it was launched for.

A worker is a copy of the program started from the top with ``--batch
<kind>:...``, and the time-lag pre-pass (time-lag optimisation or PWB cache
generation) comes before the planar fit. Nothing told a planar-fit worker to
stand aside, so on a project with both on the fly it walked its planar-fit
slice as time-lag periods, wrote time-lag records under the planar-fit dump
name, and stopped. The parent read those records as wind means, filtered the
garbage out, and every sector of the fit failed (-9999), where a serial run
fitted. No fixture had both pre-passes on, so nothing saw it;
``base_gappy_tlag`` (``gen_gappy.py``) has both, and ``check_parallel.sh`` on
it is the end-to-end gate.

Two guards now: the time-lag pre-pass is skipped by a worker of another kind,
and each dump carries its kind, which the writer and the merge both check -
so a mix-up stops the run instead of becoming data.

Part of the EddyFlow engine's static checks.
"""

from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[1]


def code(rel):
    text = (ROOT / rel).read_text(encoding="utf-8", errors="replace")
    return "\n".join(ln for ln in text.splitlines() if not ln.lstrip().startswith("!"))


MAIN = code("src/src_rp/eddyflow-rp_main.f90")
# The worker pool itself moved to src_common/batch_pool.f90, shared with FCC;
# what is checked here is read from both.
PARALLEL = code("src/src_rp/prepass_parallel.f90") + code("src/src_common/batch_pool.f90")


def body(name):
    return PARALLEL[PARALLEL.index("subroutine %s(" % name):PARALLEL.index("end subroutine %s" % name)]


class APlanarFitWorkerSkipsTheTimeLagPrepass(unittest.TestCase):

    def test_the_time_lag_block_admits_only_the_parent_and_to_workers(self):
        i = MAIN.index("if ((trim(adjustl(Meth%tlag)) == 'tlag_opt' .or. PwbCacheGenerate) .and.")
        head = MAIN[i:MAIN.index("then", i) + 4]
        self.assertIn("(BatchIndex == 0 .or. BatchKind == 'to')", head)

    def test_the_time_lag_block_precedes_the_planar_fit(self):
        """Which is why the planar-fit worker meets it at all."""
        self.assertLess(MAIN.index("call StartPrepassBatches('to'"),
                        MAIN.index("call StartPrepassBatches('pf'"))


class EveryDumpSaysWhichPrepassItBelongsTo(unittest.TestCase):

    def test_writers_refuse_the_wrong_kind_and_record_their_own(self):
        for name, kind in (("WriteTlagBatchDump", "to"), ("WritePwbBatchDump", "to"),
                           ("WritePfBatchDump", "pf")):
            b = body(name)
            self.assertIn("call RequireBatchKind('%s')" % kind, b)
            self.assertLess(b.index("call RequireBatchKind("), b.index("open(newunit"))
            self.assertIn("write(u) BatchKind", b)

    def test_merges_check_the_kind_they_read(self):
        for name, expect in (("MergeTlagBatchDumps", "kind"), ("MergePwbBatchDumps", "kind"),
                             ("MergePfBatchDumps", "'pf'")):
            b = body(name)
            self.assertIn("read(u) dumpKind", b)
            self.assertIn("if (dumpKind /= %s) &" % expect, b)
            self.assertLess(b.index("read(u) idx, idxCount"), b.index("read(u) dumpKind"))

    def test_the_guard_stops_the_worker(self):
        b = body("RequireBatchKind")
        self.assertIn("if (BatchKind == kind) return", b)
        self.assertIn("error stop", b)


if __name__ == "__main__":
    unittest.main()
