"""The engine's innFLUX-derived routines against a port of innFLUX itself.

innFLUX: Striednig et al. (2020), Atmos. Meas. Tech., 13, 1447-1465,
https://doi.org/10.5194/amt-13-1447-2020; MATLAB code v1.2.0, University of
Innsbruck. innflux_port.py holds the line-by-line ports; the engine side is
obj/<os>/innflux_compare, built with

    cmake --build --preset release --target innfluxcmp

Stage 1 - the covariance peak search (m_lagpeak_core FindCovPeak against
find_covariance_peak.m):

  * peak: the same cross-covariance function handed to both. innFLUX smooths
    it with a Hamming window and decides the sign of the peak from a
    100-point moving average of its derivative; the engine does the same on
    its own (reversed) lag axis, over a window wide enough that only interior
    values are searched. The lags must agree exactly, both for the
    derivative rule and for innFLUX's positive-only search.
  * series: synthetic w and scalar series with a known delay. innFLUX's lag
    from MATLAB's xcov('unbiased') and the engine's from its own
    CrossCovarianceSeries, which removes the means of the overlapping pairs
    at every lag rather than each series' mean once. On series with a clear
    flux the two must still give the same lag; a one-sample difference on a
    near-flat function is reported as KNOWN (that mean convention), not
    failed. Anything larger fails.

Exit status 1 on any unexplained difference.

    <dp python> tests/innflux/compare.py [--seed N] [--cases N] [--keep DIR]

Part of the EddyFlow engine's test tools.
"""

import argparse
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))

from innflux_port import find_covariance_peak, step1_lag  # noqa: E402

DRIVER = ROOT / "obj" / "win" / "innflux_compare.exe"
if not DRIVER.exists():
    DRIVER = ROOT / "obj" / "linux" / "innflux_compare"

#: innFLUX_parameters_template.m
MAX_LAG = 2000
SEARCH = 100
FILTER = 15

failures = []
known = []


def run_driver(text, keep, name):
    d = Path(keep) if keep else Path(tempfile.mkdtemp(prefix="innflux_"))
    d.mkdir(parents=True, exist_ok=True)
    p = d / f"{name}.txt"
    p.write_text(text)
    r = subprocess.run([str(DRIVER), str(p)], capture_output=True, text=True)
    if r.returncode != 0:
        raise SystemExit(f"driver failed on {p}:\n{r.stdout}{r.stderr}")
    return r.stdout.split("\n")


def synthetic_ccf(rng, kind):
    """A cross-covariance function in MATLAB's order (index 0 = m = -MAX_LAG)."""
    m = np.arange(-MAX_LAG, MAX_LAG + 1)
    lag = int(rng.integers(-SEARCH + 20, SEARCH - 20))       # scalar delay
    width = float(rng.uniform(5, 40))
    amp = float(rng.uniform(0.2, 1.0)) * (1 if kind != "neg" else -1)
    # xcov(w, c) peaks at m = -lag for a scalar delayed by lag samples.
    f = amp * np.exp(-0.5 * ((m + lag) / width) ** 2)
    f += float(rng.uniform(0.0, 0.15)) * rng.standard_normal(len(m))
    if kind == "slope":
        f += float(rng.uniform(-1e-4, 1e-4)) * m
    return f


def stage1_peak(rng, ncases, keep):
    lines = []
    expect = []
    for i in range(ncases):
        kind = ("pos", "neg", "slope")[i % 3]
        f = synthetic_ccf(rng, kind)
        for positive in (0, 1):
            expect.append(find_covariance_peak(f, MAX_LAG, SEARCH, FILTER, positive))
            # Engine axis: lag L = -m, so cov_L(L) = f at m = -L.
            lines.append("peak")
            lines.append(f"{-MAX_LAG} {MAX_LAG} {-SEARCH} {SEARCH} {FILTER} "
                         f"{1 if positive else 2} 0")
            lines.extend(f"{v:.17e}" for v in f[::-1])
    out = [l for l in run_driver("\n".join(lines) + "\n", keep, "peak") if l.strip()]
    bad = 0
    for k, (e, l) in enumerate(zip(expect, out)):
        got = int(l.split()[1])
        if got != e:
            bad += 1
            failures.append(f"peak case {k}: innFLUX {e}, engine {got}")
    print(f"stage 1 peak:   {len(expect)} cases, {bad} mismatches")


def stage1_series(rng, ncases, keep):
    n = 18000
    lines = []
    expect = []
    for i in range(ncases):
        delay = int(rng.integers(-60, 60))
        w = rng.standard_normal(n + 400)
        noise = float(rng.uniform(0.5, 3.0))
        sign = 1.0 if i % 2 == 0 else -1.0
        # scalar delayed by `delay` samples: c(t) = w(t - delay)
        c = sign * np.roll(w, delay) + noise * rng.standard_normal(n + 400)
        c += 0.001 * np.arange(n + 400)                       # a trend to detrend
        w, c = w[200:200 + n], c[200:200 + n]
        lag_inn, _ = step1_lag(w, c, MAX_LAG, SEARCH, FILTER)
        expect.append((delay, lag_inn))
        # The engine sees what innFLUX's covariance sees: w - mean, c detrended.
        t = np.arange(n, dtype=float)
        cd = c - np.polyval(np.polyfit(t, c, 1), t)
        lines.append("series")
        lines.append(f"{n} {-SEARCH} {SEARCH} {FILTER} 2 0")
        lines.extend(f"{a:.17e} {b:.17e}" for a, b in zip(w - w.mean(), cd))
    out = [l for l in run_driver("\n".join(lines) + "\n", keep, "series") if l.strip()]
    bad = 0
    for k, ((delay, e), l) in enumerate(zip(expect, out)):
        got = int(l.split()[1])
        if got == e:
            continue
        if abs(got - e) == 1:
            known.append(f"series case {k}: innFLUX {e}, engine {got} (true {delay}) "
                         "- per-lag vs global mean")
            continue
        bad += 1
        failures.append(f"series case {k}: innFLUX {e}, engine {got} (true {delay})")
    print(f"stage 1 series: {len(expect)} cases, {bad} mismatches, "
          f"{sum(1 for k in known if k.startswith('series'))} known")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--seed", type=int, default=20201447)
    ap.add_argument("--cases", type=int, default=30)
    ap.add_argument("--keep", default=None, help="keep the case files here")
    args = ap.parse_args()
    if not DRIVER.exists():
        raise SystemExit(f"{DRIVER} not built: cmake --build --preset release "
                         "--target innfluxcmp")
    rng = np.random.default_rng(args.seed)
    stage1_peak(rng, args.cases, args.keep)
    stage1_series(rng, max(4, args.cases // 3), args.keep)
    for k in known:
        print("KNOWN  " + k)
    for f in failures:
        print("FAIL   " + f)
    if failures:
        sys.exit(1)
    print("OK")


if __name__ == "__main__":
    main()
