"""The PWB period labels exist before a parallel pre-pass merges into them.

A parallel PWB pre-pass appends each worker's periods after the parent's own
slice, writing each period's date and time into PwbOptDate and PwbOptTime.
Those were allocated in RecordPwbTimelagOptPeriod, on the first period the
parent recorded itself. On CH-LAE below canopy the parent's slice - the first
two days of the range - had no usable period, so the merge wrote the first
worker period through an unallocated array: a segmentation fault after all 96
pieces had finished.

Pinned: the labels are allocated with the dataset they label and at its size;
the record routine no longer allocates them; and the merge checks for them
before its first write rather than faulting.
"""

import unittest
from pathlib import Path

SRC = Path(__file__).resolve().parent.parent / "src" / "src_rp"


def read(name):
    return (SRC / name).read_text(encoding="utf-8", errors="replace")


MAIN = read("eddyflow-rp_main.f90")
HANDLE = read("pwb_timelag_handle.f90")
PARALLEL = read("prepass_parallel.f90")


def between(src, start, end):
    a = src.index(start)
    return src[a:src.index(end, a)]


class TheLabelsComeWithTheDataset(unittest.TestCase):

    def test_allocated_beside_the_dataset_at_its_size(self):
        block = between(MAIN, "if (Meth%tlag == 'pwb') then\n", "end if")
        self.assertIn("allocate(PwbTimelagOpt(PwbTimelagOptSize))", block)
        self.assertIn("allocate(PwbOptDate(PwbTimelagOptSize), "
                      "PwbOptTime(PwbTimelagOptSize))", block)

    def test_before_the_pre_pass(self):
        self.assertLess(MAIN.index("allocate(PwbOptDate(PwbTimelagOptSize)"),
                        MAIN.index("to_periods_loop: do"))

    def test_the_record_routine_does_not_allocate(self):
        fn = between(HANDLE, "subroutine RecordPwbTimelagOptPeriod(",
                     "end subroutine RecordPwbTimelagOptPeriod")
        self.assertNotIn("allocate(PwbOptDate", fn)
        self.assertIn("error stop 'The PWB period labels were not allocated "
                      "before the pass.'", fn)


class TheMergeChecksFirst(unittest.TestCase):

    def test_the_guard_comes_before_the_first_write(self):
        fn = between(PARALLEL, "subroutine MergePwbBatchDumps(",
                     "end subroutine MergePwbBatchDumps")
        guard = fn.index("error stop 'The PWB period labels were not "
                         "allocated before the merge.'")
        self.assertLess(guard, fn.index("PwbOptDate(nOpt) = sdate(i)"))
        self.assertIn("size(PwbOptDate) < nmax", fn)


if __name__ == "__main__":
    unittest.main()
