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
