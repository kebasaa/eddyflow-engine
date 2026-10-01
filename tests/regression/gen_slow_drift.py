#!/usr/bin/env python3
"""Build the drifting slow-instrument fixtures: a known lag between two rows.

The decisive test for detecting a slow gas's time lag. The 1 Hz analyser at
Yatir is not on a fixed phase of the 20 Hz rows: it runs at about 0.987 Hz,
so its spacing is mostly 20 rows, often 21, and a missed sample leaves 40.
And the lag it has to be detected at is not a whole number of its own
intervals. A method that only ever looks at the gas's rate returns a multiple
of about one second; the true lag is somewhere between.

What it builds, from the same three hours of CH-LAE the other slow fixtures
use (10 Hz rows):

  * the MIRO's six columns (7-12) left real only on a drifting sample grid -
    one sample every 10.135 rows (0.987 Hz), about 1 % of samples missed -
    and -9999 elsewhere, the shape a slow instrument writes into a faster
    file;

  * its COS column (10) replaced at those samples by a signal MADE from the
    sonic's w, delayed by exactly 3.35 s - 33.5 rows, between two rows - plus
    noise. So the lag to recover is known, and it is not on the gas's grid
    and not even on the rows' grid: the best a row-resolution method can say
    is 3.3 or 3.4 s, either within 0.05 s.

Three variants, each with its own data folder and project:

  base_slow_drift         point-sampled, the delayed w smoothed over a centred
                          1 s window first - a peak broader than one gas
                          interval, as tube smoothing gives.
  base_slow_drift_sharp   point-sampled, no smoothing: a peak much narrower
                          than one gas interval, which tests the assumption
                          the refinement rests on (the true peak lies within
                          one interval of the coarse one).
  base_slow_drift_integr  the smoothed signal averaged over each sample's
                          interval, with the instrument declared integrating
                          (base_slow_integr.metadata), so w is paired the same
                          way.

The gate is the COS time lag the run reports per period: 3.35 +/- 0.05 s
(check_slow_lag.sh).

Usage:  python gen_slow_drift.py

Writes the projects beside this script and the raw files to
data_slow_drift_{smooth,sharp,integr}/ beside it, which are gitignored. The
source is the Lagern repository's own data folder, read and never written.
"""

import math
import random
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent

HEADER_ROWS = 1
W_COL = 3            # 1-based: the sonic's w
MIRO_COLS = [7, 8, 9, 10, 11, 12]
COS_COL = 10
FILL = "-9999"
PERIOD_ROWS = 10.135  # 0.987 Hz against 10 Hz rows
MISS_PROB = 0.01
LAG_ROWS = 33.5       # 3.35 s at 10 Hz
SMOOTH_ROWS = 11      # centred, 1.1 s
COS_MEAN, COS_GAIN, COS_NOISE = 0.50, 0.05, 0.005   # ppb, ppb per m/s, ppb
SEED = 20261001

FILES = [
    f"CH-LAE_ec_preproc_10hz_20250601-{h:02d}{m:02d}.csv"
    for h in range(3)
    for m in (0, 30)
]

VARIANTS = {
    "base_slow_drift": ("smooth", "base_slow.eddyflow", "base_slow.metadata"),
    "base_slow_drift_sharp": ("sharp", "base_slow.eddyflow", "base_slow.metadata"),
    "base_slow_drift_integr": ("integr", "base_slow_integr.eddyflow", "base_slow_integr.metadata"),
}


def ini_value(text, key):
    m = re.search(rf"^{re.escape(key)}=(.*)$", text, re.M)
    if not m:
        raise SystemExit(f"{key} not found")
    return m.group(1).strip()


def set_key(text, key, value):
    pat = re.compile(rf"^{re.escape(key)}=.*$", re.M)
    if not pat.search(text):
        raise SystemExit(f"{key} not found")
    return pat.sub(lambda m: f"{key}={value}", text)


def read_files(src_dir):
    """Every row of the six files, as lists of fields, plus their headers."""
    headers, rows, counts = [], [], []
    for name in FILES:
        src = src_dir / name
        if not src.exists():
            raise SystemExit(f"missing source file {src}")
        with open(src, encoding="utf-8", errors="replace") as f:
            lines = f.read().splitlines()
        headers.append(lines[:HEADER_ROWS])
        body = [ln.split(",") for ln in lines[HEADER_ROWS:]]
        rows.extend(body)
        counts.append(len(body))
    return headers, rows, counts


def w_series(rows):
    """w at every row, gaps bridged linearly - it only drives the signal."""
    w = []
    for r in rows:
        try:
            v = float(r[W_COL - 1])
        except (ValueError, IndexError):
            v = float("nan")
        w.append(v if v > -300 else float("nan"))
    last = next((v for v in w if not math.isnan(v)), 0.0)
    for i, v in enumerate(w):
        if math.isnan(v):
            w[i] = last
        else:
            last = v
    return w


def at(series, x):
    """series at fractional row x, linearly interpolated, clamped."""
    n = len(series)
    if x <= 0:
        return series[0]
    if x >= n - 1:
        return series[-1]
    i = int(math.floor(x))
    f = x - i
    return series[i] * (1 - f) + series[i + 1] * f


def main():
    template = (HERE / "base_slow.eddyflow").read_text(encoding="utf-8")
    src_dir = Path(ini_value(template, "data_path"))
    if "data_slow" in str(src_dir):
        # base_slow points at its own decimated copy; the source is the
        # Lagern folder base_n_gas reads.
        src_dir = Path(ini_value((HERE / "base_n_gas.eddyflow").read_text(encoding="utf-8"), "data_path"))
    headers, rows, counts = read_files(src_dir)
    n = len(rows)
    w = w_series(rows)

    half = SMOOTH_ROWS // 2
    smooth = [sum(w[max(0, i - half):min(n, i + half + 1)]) / (min(n, i + half + 1) - max(0, i - half))
              for i in range(n)]

    rng = random.Random(SEED)
    samples = []
    t = rng.uniform(0, PERIOD_ROWS)
    while t < n - 1:
        if rng.random() >= MISS_PROB:
            samples.append(int(round(t)))
        t += PERIOD_ROWS
    noise = {r: rng.gauss(0, COS_NOISE) for r in samples}
    stride = int(round(PERIOD_ROWS))

    for name, (kind, project, meta) in VARIANTS.items():
        out_dir = HERE / f"data_slow_drift_{kind}"
        out_dir.mkdir(parents=True, exist_ok=True)
        base = smooth if kind in ("smooth", "integr") else w
        sample_set = set(samples)
        cos = {}
        for r in samples:
            if kind == "integr":
                lo = max(0, r - stride + 1)
                v = sum(at(base, k - LAG_ROWS) for k in range(lo, r + 1)) / (r + 1 - lo)
            else:
                v = at(base, r - LAG_ROWS)
            cos[r] = COS_MEAN + COS_GAIN * v + noise[r]

        start = 0
        for fname, header, count in zip(FILES, headers, counts):
            with open(out_dir / fname, "w", encoding="utf-8", newline="") as fo:
                for h in header:
                    fo.write(h + "\n")
                for i in range(start, start + count):
                    f = list(rows[i])
                    if i in sample_set:
                        f[COS_COL - 1] = f"{cos[i]:.6f}"
                    else:
                        for c in MIRO_COLS:
                            if c - 1 < len(f):
                                f[c - 1] = FILL
                    fo.write(",".join(f) + "\n")
            start += count

        p = (HERE / project).read_text(encoding="utf-8")
        p = set_key(p, "data_path", out_dir.as_posix())
        p = set_key(p, "proj_file", (HERE / meta).as_posix())
        (HERE / f"{name}.eddyflow").write_text(p, encoding="utf-8")
        print(f"wrote {name}.eddyflow  ({kind}, {len(samples)} samples of {n} rows, {meta})")


if __name__ == "__main__":
    main()
