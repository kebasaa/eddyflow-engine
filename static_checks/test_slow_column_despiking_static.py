"""Guards two places where a row was mistaken for a sample.

Both were invisible while every column sampled at the file's rate, and both
became wrong the moment an instrument could declare a slower one. They are
guarded here rather than by a fixture because the arithmetic reduces exactly to
the old arithmetic on a single-rate dataset - a regression run cannot tell the
fixed code from the broken code unless the fixture itself is multi-rate.

1. Vickers-Moncrieff despiking counted consecutive ROWS against sr%num_spk,
   which counts consecutive SAMPLES. A slower column carries the error code
   between its samples, so a single-sample spike counted as `stride` and, with
   the default allowance of three, every run in a column at a quarter of the
   file rate or slower was discarded as a physical excursion. Worse, an
   accepted run's interpolation wrote the whole row span, turning the error
   rows inside it into fabricated samples - which is how a period with too
   little real data passed the CEC 90% gate.

2. ExpWeightAvrg computed the trend for EVERY column inside a loop over one
   column, `Trend(i, :) = ... Set(i, :)`, guarded on column k. At one rate the
   assignment was identical on each pass and the mistake was idempotent; at
   mixed rates a fast column's pass folds -9999 from a slow column's padding
   into that column's trend.
"""

import re
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

SPIKE = "src/src_rp/test_spike_detection_vickers_97.f90"
FLUCT = "src/src_common/fluctuations.f90"


def read(rel):
    return (ROOT / rel).read_text(encoding="utf-8", errors="replace")


def code_lines(body):
    """Source lines with comments and blanks dropped."""
    out = []
    for line in body.splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("!"):
            continue
        out.append(line)
    return out


class SpikeRunCountsSamplesTests(unittest.TestCase):
    def test_an_error_row_no_longer_extends_the_run(self):
        """The exact statement that made a gap look like a spike run."""
        for line in code_lines(read(SPIKE)):
            assert "cnt = cnt + 1" not in line or "/= 0" not in line, (
                "an error row advances the spike counter again: a gap between "
                "a slow column's samples will be counted against sr%num_spk "
                "and the column will stop being despiked\n  " + line.strip()
            )

    def test_the_run_is_located_by_its_own_bounds(self):
        """`i - k` arithmetic assumes the run's rows and samples coincide."""
        body = read(SPIKE)
        assert "run_start" in body and "run_last" in body, (
            "the run's row extent is no longer tracked, so the flagging and "
            "the interpolation are back to assuming one sample per row"
        )
        for line in code_lines(body):
            assert not re.search(r"IsSpike\(\s*i\s*-\s*k\s*,", line), (
                "spike flags are addressed by counting rows back from the run "
                "end again\n  " + line.strip()
            )

    def test_the_interpolation_never_writes_an_error_row(self):
        """Filling padding manufactures samples the instrument never reported,
        and the valid-sample count is what CEC and the completeness test read."""
        body = read(SPIKE)
        assert re.search(
            r"if \(Set\(k, j\) /= error\)\s*&\s*Set\(k, j\) = m \* dble",
            body), (
            "the interpolation writes every row of the run's span, so a slower "
            "column's error rows become fabricated samples"
        )

    def test_the_spike_percentage_is_of_the_columns_own_samples(self):
        body = read(SPIKE)
        assert "ColumnAcFreq(j)" in body, (
            "the hf_lim percentage is back on the file's row count, which "
            "divides a slower column's spike count by its stride"
        )
        assert "dble(tot_spikes(j)) / dble(N)" not in body, (
            "the spike percentage divides by the row count again"
        )


class ExpWeightAvrgWritesOneColumnTests(unittest.TestCase):
    def test_the_trend_assignment_names_its_own_column(self):
        body = read(FLUCT)
        start = body.index("subroutine ExpWeightAvrg")
        fn = body[start:body.index("end subroutine ExpWeightAvrg", start)]
        for line in code_lines(fn):
            assert not re.search(r"\b(Trend|Set|Primes)\([^)]*,\s*:\s*\)", line), (
                "ExpWeightAvrg assigns across every column from inside its "
                "loop over one column; at mixed sampling rates that writes "
                "one column's error code into another column's trend\n  "
                + line.strip()
            )

    def test_the_sibling_routines_did_not_acquire_the_same_shape(self):
        """RunningMean is the same algorithm and had it right; keep it that
        way, and keep LinDetrend that way too."""
        body = read(FLUCT)
        for name in ("RunningMean", "LinDetrend"):
            start = body.index("subroutine %s" % name)
            fn = body[start:body.index("end subroutine %s" % name, start)]
            for line in code_lines(fn):
                assert not re.search(
                    r"\b(Trend|Primes)\(\s*i\s*,\s*:\s*\)", line), (
                    "%s now writes every column from a single column's "
                    "iteration\n  %s" % (name, line.strip())
                )


if __name__ == "__main__":
    unittest.main()
