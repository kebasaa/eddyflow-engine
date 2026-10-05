# PWB against dyco

`compare.py` runs the engine's PWB time-lag detector and dyco's
(github.com/holukas/dyco, the reference implementation of Vitale et al. 2024)
on the **same bootstrap resamples** and compares every number both produce.

`static_checks/test_pwb_reference_static.py` already pins the deterministic half
of the chain (unit-root test, AR fits, pre-whitened CCF) to RFlux's frozen
output. It stops at the bootstrap because the engine and R share no random
stream. This goes the rest of the way: a stub generator records the block
starts dyco draws, and the engine's chain runs on exactly those.

## Running it

```
cd prj && mingw32-make pwbcmp           # builds obj/win/pwb_compare.exe
```

dyco must be checked out as `dyco-main` beside the repository root (it is
gitignored), and the Python needs numpy, scipy, pandas and matplotlib. On the
development machine that is the conda env `dp`, whose DLL folders have to be on
PATH because it is not activated:

```
E=/c/Users/jonmuell/AppData/Local/miniconda3/envs/dp
export PATH="$E:$E/Library/mingw-w64/bin:$E/Library/usr/bin:$E/Library/bin:$E/Scripts:$PATH"
$E/python.exe tests/pwb_dyco/compare.py --files 6
$E/python.exe tests/pwb_dyco/compare.py --no-stage1 --table <run>/..._pwb_timelag_adv.csv
```

Exit status 1 on any unexplained difference.

## Stage 1: one period at a time

Cases: dyco's two synthetic reference fixtures (stationary and differencing
branch of the unit-root test, 20 Hz, symmetric and positive windows) and real
30-minute CH-LAE periods from `tests/regression/data_gappy` (10 Hz; LI-7200
and MIRO gases, the MIRO columns with their gaps), wind double-rotated with
dyco's own `WindDoubleRotation`. Data are only read.

`obj/<os>/pwb_compare` (`src/src_tools/pwb_compare_main.f90`) runs the same
core routines `PwbDetectGas` does - `FillMissingLinear`, `PwbPreWhiten`,
`PwbBootstrapCombination` four times, `PwbSummariseBootstrap`,
`PwbBestCombination` - linked against `m_pwb_core` alone.

Compared per period: the unit-root decision; the three AR orders and first
coefficients; the full-data pre-whitened CCF peak; the block length; and per
combination (cw, wc, ct, tc) the 99 replicate peak lags (exactly), the mean
smoothed CCF over the window (to 1e-11 absolute - the engine sums directly,
dyco by FFT), the 95 % HDI (exactly), the mode (exactly) and the |CCF| at the
mode; then the chosen combination and the final lag.

dyco jitters the replicate lags by N(0, 1e-4) before its kernel density; the
engine has no random stream there, so dyco is handed zero jitter.

**Reported, not failed:** `PINNED` - the engine prefers a combination whose mode
is not on the window edge before comparing magnitudes (m_pwb_core departure 2);
dyco compares magnitudes only and then throws the period away as edge-pinned.

Result on 2026-10-02: every number identical in all periods but one, the one
being a `PINNED` case. Before the engine took dyco's mode estimator
(`MapLagEstimate`), this found the final lag a record apart in 6 of 40 periods.

## Stage 2: the reliability classes across periods

`--table` takes an engine half-hourly PWB table (`*_pwb_timelag_*.csv`, e.g.
from `base_gappy_pwb`) and runs each gas's raw detections through dyco's
`apply_hdi_prefilter` and `apply_pwbopt` with the table's own thresholds. Every
period dyco settles as S1 or S2 must be S1 or S2 in the engine with the same
lag, and where both carry a lag forward it must be the same lag. The engine's
further fills - interpolation between settled periods, the carry limit in
hours, sharing within an analyser, the terminal fallback - are its documented
departures and are counted, not failed.
