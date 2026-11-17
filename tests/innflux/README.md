# The engine against innFLUX

`compare.py` runs the engine's innFLUX-derived routines and a line-by-line
Python port of innFLUX's own MATLAB code (`innflux_port.py`) on the same inputs
and compares the results.

innFLUX is the University of Innsbruck's eddy-covariance code:
Striednig, M., Graus, M., Märk, T. D., and Karl, T. G.: InnFLUX – an
open-source code for conventional and disjunct eddy covariance analysis of
trace gas measurements: an urban test case, Atmos. Meas. Tech., 13, 1447–1465,
https://doi.org/10.5194/amt-13-1447-2020, 2020. The port follows v1.2.0
(2023-01-14). No MATLAB is needed.

Every innFLUX-derived option in the engine is **off by default**.

## Running it

```
cmake --build --preset release --target innfluxcmp   # builds obj/win/innflux_compare.exe
```

The Python needs numpy and scipy; on the development machine that is the conda
env `dp`, whose DLL folders have to be on PATH because it is not activated:

```
E=/c/Users/jonmuell/AppData/Local/miniconda3/envs/dp
export PATH="$E:$E/Library/mingw-w64/bin:$E/Library/usr/bin:$E/Library/bin:$E/Scripts:$PATH"
$E/python.exe tests/innflux/compare.py [--cases N] [--seed N] [--keep DIR]
```

Exit status 1 on any unexplained difference.

## Stage 1: the covariance peak search

`m_lagpeak_core` `FindCovPeak` against `find_covariance_peak.m`. innFLUX
smooths the cross-covariance with a normalised Hamming window through MATLAB
`filtfilt` (smoothing before the peak search after Taipale et al. 2010), and
decides whether the peak is a maximum or a minimum from a 100-point moving
average of the function's derivative.

- **peak**: the same function handed to both. The lags must agree exactly, for
  the derivative rule and for innFLUX's positive-only search (the one its step
  2 uses on summed absolute functions). The engine works on its own lag axis
  (positive = scalar delayed), which is MATLAB's `xcov` axis reversed.
- **series**: synthetic `w` and scalar series with a known delay, noise and a
  trend. innFLUX takes its covariance from `xcov(..., 'unbiased')`, which
  removes each series' mean once; the engine's `CrossCovarianceSeries` removes
  the means of the overlapping pairs at every lag. A one-sample difference
  from that is reported as KNOWN; anything larger fails.

Result 2026-11-17, seeds 1–3 and the default, 60 cases each: 0 mismatches,
0 known.

## Stage 2: the clock-drift estimate

`m_driftfix_core` `SegmentedDrift` against a port of `tools/driftfix/innFLUX_driftfix.m`
(innFLUX v1.2.0, 2023; based on a routine by M. Graus, adapted by M. Striednig; not
described in Striednig et al. 2020). A synthetic run of 20 days of 30-minute
periods whose lag drifts by 0.15 s per day and jumps by 3 s on day 12, seen
through three reference gases with 0.3 s of noise and 5 % outliers.

- The engine and the port must agree to 1e-9 s.
- Both must find exactly one jump.
- The recovered drift must be within 0.3 s RMS of the injected one.

The robust LOWESS follows MATLAB's `smooth(..., 'rlowess')`: the span nearest
points in index order, tricube weights, five bisquare passes on residuals
scaled by six times their median (Cleveland 1979). innFLUX places each jump
with MATLAB's `risetime`. Both sides here use the midpoint of the first 10 %
and 90 % crossings within two days instead, which is a documented departure.

Result 2026-11-17, seeds 1–5: one jump each; engine and port within 1.6e-12 s;
RMS against the injected drift 0.015–0.167 s.
