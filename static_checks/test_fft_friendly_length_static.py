"""power_of_two=2: the (co)spectra on the largest fast FFT length.

A power of two leaves out 9 % of 30 minutes at 20 Hz and 43 % at 32 Hz. The
FFT here is FFTPACK, fast for any length without a large prime factor -
measured per transform: 32768 rows 0.12 ms, 36000 0.14 ms, 35750 (with 11 and
13) 0.21 ms, 35760 (with 149) 0.33 ms, 35761 (with 3251) 11 ms, 35759 (prime)
865 ms. FftFriendlyLength takes the largest length with no prime factor above
13, which loses at most 1 % of the rows at any rate.

Pinned: the reader takes 2 as this choice and 1 as the power of two; the
production loop asks FftFriendlyLength only when the power of two was not
asked for; and the helper's cut-off is 13. The helper itself is also run,
when gfortran is on the PATH, against row counts at 5, 10, 20 and 32 Hz.
"""

import os
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / "src"
HELPER = (SRC / "src_common" / "power_of_two.f90").read_text(encoding="utf-8", errors="replace")
INI = (SRC / "src_rp" / "read_ini_rp.f90").read_text(encoding="utf-8", errors="replace")
MAIN = (SRC / "src_rp" / "eddyflow-rp_main.f90").read_text(encoding="utf-8", errors="replace")


class TheChoiceIsRead(unittest.TestCase):

    def test_two_is_the_fast_length_and_one_the_power_of_two(self):
        self.assertIn("RPsetup%power_of_two = SCTags(87)%value(1:1) == '1'", INI)
        self.assertIn("RPsetup%fft_friendly_len = SCTags(87)%value(1:1) == '2'", INI)

    def test_the_loop_asks_for_it_after_the_power_of_two(self):
        a = MAIN.index("call PowerOfTwo(max_nsmpl, SpecRow)")
        b = MAIN.index("elseif (RPsetup%fft_friendly_len) then")
        c = MAIN.index("call FftFriendlyLength(max_nsmpl, SpecRow)")
        self.assertLess(a, b)
        self.assertLess(b, c)

    def test_the_cut_off_is_thirteen(self):
        self.assertIn("integer, parameter :: primes(6) = [2, 3, 5, 7, 11, 13]", HELPER)


@unittest.skipUnless(shutil.which("gfortran"), "gfortran not on PATH")
class TheHelperRuns(unittest.TestCase):

    def test_lengths_at_the_usual_rates(self):
        driver = (
            "program t\n"
            "    implicit none\n"
            "    integer :: n(10) = [9000, 18000, 36000, 57600, 72000, 115200, "
            "35760, 35759, 57120, 1]\n"
            "    integer :: i, m\n"
            "    do i = 1, size(n)\n"
            "        call FftFriendlyLength(n(i), m)\n"
            "        print '(i0,1x,i0)', n(i), m\n"
            "    end do\n"
            "end program t\n")
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "t.f90")
            with open(p, "w") as f:
                f.write(driver)
            exe = os.path.join(d, "t.exe")
            subprocess.run(["gfortran", "-o", exe, p,
                            str(SRC / "src_common" / "power_of_two.f90")],
                           check=True, capture_output=True)
            out = subprocess.run([exe], check=True, capture_output=True, text=True).stdout
        got = dict(tuple(map(int, line.split())) for line in out.split("\n") if line.strip())

        def smooth(k):
            for p in (2, 3, 5, 7, 11, 13):
                while k % p == 0:
                    k //= p
            return k == 1

        for n, m in got.items():
            self.assertLessEqual(m, n)
            self.assertTrue(smooth(m), "%d has a prime factor above 13" % m)
            self.assertFalse(any(smooth(k) for k in range(m + 1, n + 1)),
                             "%d is not the largest such length up to %d" % (m, n))
            if n >= 1000:
                self.assertLess((n - m) / n, 0.01)
        #> The whole half-hour or hour at every usual rate, nothing left out.
        for n in (9000, 18000, 36000, 57600, 72000, 115200):
            self.assertEqual(got[n], n)


if __name__ == "__main__":
    unittest.main()
