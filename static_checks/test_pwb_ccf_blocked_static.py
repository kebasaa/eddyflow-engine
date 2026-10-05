"""The PWB cross-correlation sums a block of lags at once, and moves no bit.

``ComputeCcfWindow`` is the inner loop of the whole PWB pre-pass: 4
pre-whitening combinations x 99 bootstrap replicates x every gas x every
period. It used to be one accumulator per lag, summed a lag at a time::

    cov = cov + xc(i) * yc(i + lag)

Every add there waits for the previous one, so the loop ran at the latency of
a floating-point add, not its throughput. On the Yatir run that was most of the
~25 s each period cost.

Now a block of 32 lags shares one pass over ``i``. Each lag still has its own
accumulator and still receives the same products in the same order - ``i``
ascending, first over the range every lag in the block shares, then each lag's
own remaining products - so every result is bitwise the same. Measured at the
production shape (36 000 records, 581 lags, 396 calls): 5.7x, against 2.9x for
blocks of 8, 3.1x for 16, 4.4x for 64 and 3.4x for 128.

Two things make the bitwise claim true, and both are pinned here:

* **No reassociation.** gfortran does not reorder floating-point sums without
  -ffast-math, and nothing here asks it to.
* **No fused multiply-add.** A fused product is rounded once instead of twice.
  The default target has no FMA instruction, so nothing fuses today; the build
  says -ffp-contract=off so that a future -march cannot change that silently.

The proof is ``ccf_bitexact`` (``mingw32-make ccfcheck``): the pre-rewrite
routine frozen verbatim, compared with the shipped one on the raw bits of every
element over 6331 cases - random shapes, every window edge for small n, zero
variance, the production shape. It ran clean.

Part of the EddyFlow engine's static checks.
"""

from pathlib import Path
import re
import subprocess
import unittest


ROOT = Path(__file__).resolve().parents[1]

CORE = ROOT / "src" / "src_common" / "m_pwb_core.f90"
HARNESS = ROOT / "src" / "src_tools" / "ccf_bitexact_main.f90"
MAKEFILE = ROOT / "prj" / "Makefile"
DRIVER = ROOT / "obj" / "win" / "ccf_bitexact.exe"


def read(path):
    return path.read_text(encoding="utf-8", errors="replace")


def code(path):
    return (chr(10)).join(ln for ln in read(path).splitlines()
                          if not ln.lstrip().startswith("!"))


CCF = code(CORE)
CCF = CCF[CCF.index("subroutine ComputeCcfWindow"):CCF.index("end subroutine ComputeCcfWindow")]


class TheLagsAreSummedInBlocks(unittest.TestCase):

    def test_the_block_size_is_the_measured_best(self):
        self.assertIn("integer, parameter :: NB = 32", CCF)

    def test_one_accumulator_per_lag(self):
        self.assertIn("real(kind = dbl) :: acc(NB)", CCF)

    def test_both_directions_share_a_pass_over_i(self):
        self.assertIn("acc = acc + xc(i) * yc(i + lag:i + lag + NB - 1)", CCF)
        self.assertIn("acc = acc + xc(i + lag:i + lag + NB - 1) * yc(i)", CCF)

    def test_each_lag_finishes_its_own_products_in_ascending_order(self):
        """The shared range first, then this lag's remaining i - which is the
        order the single-lag loop used, split in two."""
        self.assertIn("do i = shared + 1, n - (lag + k - 1)", CCF)
        self.assertLess(CCF.index("do i = 1, shared"),
                        CCF.index("do i = shared + 1, n - (lag + k - 1)"))

    def test_a_lag_with_under_two_records_is_still_zero(self):
        self.assertIn("if (n - abs(lag) <= 1) ccf(lag) = 0d0", CCF)

    def test_nothing_asks_for_reassociation(self):
        mk = read(MAKEFILE)
        for flag in ("-ffast-math", "-Ofast", "-fassociative-math"):
            self.assertNotIn(flag, mk)


class NoMultiplyAndAddIsFused(unittest.TestCase):

    def test_the_build_forbids_contraction(self):
        flags = re.search(r"(?m)^CFLAGS = (.*)$", read(MAKEFILE)).group(1)
        self.assertIn("-ffp-contract=off", flags)


class TheProofIsKeptAndRuns(unittest.TestCase):

    def test_the_harness_freezes_the_old_routine(self):
        h = read(HARNESS)
        self.assertIn("subroutine RefCcf(", h)
        #> The single-accumulator loop, exactly as it was.
        self.assertIn("cov = cov + xc(i) * yc(i + lag)", h)
        self.assertIn("cov = cov + xc(i - lag) * yc(i)", h)

    def test_it_compares_raw_bits_not_values(self):
        self.assertIn("transfer(a(k), 0_int64) /= transfer(b(k), 0_int64)", read(HARNESS))

    def test_it_is_linked_against_the_shipped_object(self):
        self.assertRegex(read(MAKEFILE), r"ccfcheck : m_numeric_kinds\.o m_pwb_core\.o")

    @unittest.skipUnless(DRIVER.exists(), "ccf_bitexact not built (mingw32-make ccfcheck)")
    def test_the_shipped_routine_is_bit_identical(self):
        r = subprocess.run([str(DRIVER)], capture_output=True, text=True, timeout=600)
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertIn("BIT-IDENTICAL", r.stdout)
        self.assertIn("mismatched elements: 0", r.stdout)


if __name__ == "__main__":
    unittest.main()
