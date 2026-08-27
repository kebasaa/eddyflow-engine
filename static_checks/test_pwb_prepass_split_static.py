"""The PWB cache pre-pass may be split, and what a worker hands back is evidence.

`prepass_parallel.f90` used to say a PWB cache pre-pass could not be split at
all: the streaming classifier decides a period partly from the last settled
detection before it, and a slice starting cold has none. That was true of the
classifier, and it was the wrong conclusion, because in a cache-generation run
the classifier does not decide the output. `PostProcessPwbTimelagCache` does,
once the whole run has been read, and it overrules every row.

What blocked the split was not the classifier but four values that carried its
verdict into the table anyway. Each is now taken from the settled table
instead, and each has its own check:

  * the terminal fallback's lag - `test_pwb_terminal_fallback_static`
  * the aggregate dataset's membership - `test_pwb_aggregate_dataset_static`
  * the donor tally - `test_pwb_donor_tally_static`
  * every field saying how a period was settled -
    `test_pwb_postpass_owns_its_fields_static`

The last was found by the split rather than before it: a serial run cannot show
a field holding a stale value, because the staleness is identical every time.

Two properties make the transport work, and they are what this file pins.

**A worker settles nothing.** It runs detection and hands back rows. The
post-pass, the aggregate rebuild and the resolve all happen once, in the
parent, over every slice concatenated.

**Order is the answer.** The post-pass sorts by timestamp with a stable
insertion sort, so rows appended in slice order come out exactly as one loop
would have left them, and rows sharing a timestamp keep their gas order.

Verified end to end by `check_parallel.sh base_pwb_par.eddyflow`: 24 periods
across 6 worker processes, 85 files byte-identical to the same run under -j 1,
the half-hourly cache among them.

Part of the EddyFlow engine's static checks.
"""

from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]

PARALLEL = "src/src_rp/prepass_parallel.f90"
RP_MAIN = "src/src_rp/eddyflow-rp_main.f90"
GLOBALS = "src/src_rp/m_rp_global_var.f90"


def read(path):
    return (ROOT / path).read_text(encoding="utf-8")


def code(path):
    return (chr(10)).join(ln for ln in read(path).splitlines()
                          if not ln.lstrip().startswith("!"))


PAR = code(PARALLEL)
MAIN = code(RP_MAIN)
GLOB = code(GLOBALS)


class TheSplitIsOffered(unittest.TestCase):

    def test_the_pre_pass_is_no_longer_refused_for_generating_a_cache(self):
        self.assertNotIn(".not. PwbCacheGenerate, toWorkers", MAIN)

    def test_the_time_lag_pre_pass_is_allowed_outright(self):
        self.assertIn(".true., toWorkers", MAIN)


class AWorkerHandsBackEvidence(unittest.TestCase):

    def test_it_dumps_the_cache_rows(self):
        self.assertIn("write(u) PwbTimelagCache(1:PwbTimelagCacheN)", PAR)

    def test_it_dumps_the_period_axis_with_them(self):
        """The aggregate dataset has no timestamp of its own, so without this
        the parent cannot tell which period a row belongs to."""
        self.assertIn("write(u) PwbOptDate(1:nOpt)", PAR)
        self.assertIn("write(u) PwbOptTime(1:nOpt)", PAR)

    def test_the_period_axis_is_reachable_from_both_sides(self):
        """It was module state private to pwb_timelag_handle, which the dump
        cannot see."""
        self.assertIn("character(10), allocatable :: PwbOptDate(:)", GLOB)
        self.assertIn("character(5), allocatable :: PwbOptTime(:)", GLOB)

    def test_the_worker_picks_the_dump_that_matches_the_run(self):
        self.assertIn("call WritePwbBatchDump", MAIN)
        self.assertIn("call WriteTlagBatchDump", MAIN)


class TheParentDoesTheSettling(unittest.TestCase):

    def test_it_merges_before_it_post_processes(self):
        self.assertLess(MAIN.index("call MergePwbBatchDumps"),
                        MAIN.index("call PostProcessPwbTimelagCache()"))

    def test_the_merge_appends_rather_than_interleaves(self):
        """Slice order is period order, which is what makes the stable sort in
        the post-pass reproduce a single loop."""
        self.assertIn("do k = 2, nEff", PAR)

    def test_the_cache_grows_once_per_slice(self):
        """StorePwbTimelagCacheAt reallocates and copies the whole table for
        every row, which is quadratic; the merge must not inherit that."""
        self.assertIn("allocate(grown(PwbTimelagCacheN + nrec))", PAR)
        self.assertIn("call move_alloc(grown, PwbTimelagCache)", PAR)


class AnOlderDumpIsRefusedRatherThanMisread(unittest.TestCase):

    def test_the_format_moved(self):
        self.assertIn("EDDYFLOW_PREPASS_04", PAR)
        self.assertNotIn("EDDYFLOW_PREPASS_03", PAR)

    def test_the_magic_is_checked_on_the_way_in(self):
        self.assertIn("if (magic /= BatchMagic) &", PAR)


if __name__ == "__main__":
    unittest.main()
