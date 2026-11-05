"""Defects found by the 8.1.8 audit of the RP and FCC chains, one class each.

Each of these gave a wrong number or flag without saying so, in a
configuration people use. The classes pin the corrected code shapes, so a
later edit that brings one back fails here rather than in someone's data.

Part of the EddyFlow engine's static checks.
"""

from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"


def code(*parts):
    text = SRC.joinpath(*parts).read_text(encoding="utf-8", errors="replace")
    return "\n".join(ln for ln in text.splitlines()
                     if not ln.lstrip().startswith("!"))


def body(src, start, end):
    return src[src.index(start):src.index(end)]


class MauderDespikingDrawsFromGoodSamples(unittest.TestCase):
    """Mauder (2013) drew its interpolation from the row before the run,
    which could be the error code: a spike beside one filtered sample became
    about (x - 9999) / 2."""

    MAUDER = code("src_rp", "test_spike_detection_mauder_13.f90")

    def test_no_line_from_the_row_before_the_run(self):
        self.assertNotIn("Set(i - (cnt + 1), j)", self.MAUDER)
        self.assertNotIn("IsSpike(i-k, j) = .true.", self.MAUDER)

    def test_the_line_starts_at_the_last_good_sample(self):
        self.assertIn("m = (Set(i, j) - Set(prev_good, j)) / dble(i - prev_good)", self.MAUDER)
        self.assertIn("if (prev_good /= 0) then", self.MAUDER)
        self.assertIn("Set(k, j) = m * dble(k - prev_good) + q", self.MAUDER)
        # Both arms track the good samples and flag samples only.
        self.assertEqual(self.MAUDER.count("prev_good = i"), 2)
        self.assertEqual(self.MAUDER.count("if (Set(k, j) /= error) IsSpike(k, j) = .true."), 2)
        self.assertEqual(self.MAUDER.count("new_spike = .not. IsSpike(run_last, j)"), 2)


if __name__ == "__main__":
    unittest.main()
