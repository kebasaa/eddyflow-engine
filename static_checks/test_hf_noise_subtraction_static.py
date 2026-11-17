"""Each gas keeps its own denoised ensemble spectrum.

SubtractHighFreqNoise loops over every gas slot. A slot with no noise range
took the else branch, which copied the raw ensemble over the denoised one for
the WHOLE array - and every unused slot has no range, so the first one after a
denoised gas put that gas's raw spectrum back. The "denoised ensemble averaged
spectrum" columns of the spectral assessment were the raw ones.

Pinned: the raw copy is made once, before the loop, and nothing inside the
loop copies the whole array.
"""

import re
import unittest
from pathlib import Path

SRC = (Path(__file__).resolve().parent.parent / "src" / "src_fcc"
       / "subtract_high_freq_noise.f90").read_text(encoding="utf-8", errors="replace")


class TheRawCopyIsMadeOnce(unittest.TestCase):

    def test_the_copy_comes_before_the_gas_loop(self):
        self.assertLess(SRC.index("dMeanBinSpec = MeanBinSpec"),
                        SRC.index("do gas = firstGas, lastGas"))

    def test_nothing_in_the_loop_copies_the_whole_array(self):
        loop = SRC[SRC.index("do gas = firstGas, lastGas"):]
        self.assertIsNone(re.search(r"dMeanBinSpec\s*=\s*MeanBinSpec", loop))
        self.assertIn("dMeanBinSpec(i, cls)%of(gas) =", loop)

    def test_there_is_no_else_branch_to_reset_it(self):
        loop = SRC[SRC.index("do gas = firstGas, lastGas"):]
        self.assertNotRegex(loop, r"\n\s*else\s*\n")


if __name__ == "__main__":
    unittest.main()
