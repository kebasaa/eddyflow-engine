"""Line-by-line Python ports of the innFLUX routines the engine reimplements.

innFLUX is the University of Innsbruck's MATLAB eddy-covariance code,
Striednig, M., Graus, M., Maerk, T. D., and Karl, T. G.: InnFLUX - an
open-source code for conventional and disjunct eddy covariance analysis of
trace gas measurements: an urban test case, Atmos. Meas. Tech., 13,
1447-1465, https://doi.org/10.5194/amt-13-1447-2020, 2020. Version ported:
v1.2.0 (2023-01-14).

Each function keeps the MATLAB code's structure and its 1-based index
arithmetic, written out in the comments, so a reader can lay the two side by
side. Nothing here is used by the engine; tests/innflux/compare.py runs these
against the engine's own routines (obj/<os>/innflux_compare, built by the
CMake target innfluxcmp).

Needs numpy and scipy (on the development machine, the conda env `dp`).

Part of the EddyFlow engine's test tools.
"""

import numpy as np
from scipy.signal import filtfilt
from scipy.signal.windows import hamming


def _filtfilt(b, x):
    """MATLAB filtfilt(b, 1, x) for an FIR filter.

    MATLAB pads by 3*(nfilt-1) samples of odd reflection and starts both
    passes from the steady-state initial conditions; scipy's method='pad'
    does the same, but its default pad length is 3*max(len(a), len(b)), one
    filter length longer - hence the explicit padlen.
    """
    b = np.asarray(b, dtype=float)
    return filtfilt(b, [1.0], np.asarray(x, dtype=float), padlen=3 * (len(b) - 1))


def xcov_unbiased(x, y, maxlag):
    """MATLAB xcov(x, y, maxlag, 'unbiased'), lags -maxlag..maxlag.

    xcov removes each series' own mean, then c(m) = sum_n x(n+m) y(n) / (N-|m|).
    Returned in MATLAB's order, index 0 = m = -maxlag.
    """
    x = np.asarray(x, dtype=float) - np.mean(x)
    y = np.asarray(y, dtype=float) - np.mean(y)
    n = len(x)
    out = np.empty(2 * maxlag + 1)
    for i, m in enumerate(range(-maxlag, maxlag + 1)):
        if m >= 0:
            s = np.dot(x[m:], y[:n - m])
        else:
            s = np.dot(x[:n + m], y[-m:])
        out[i] = s / (n - abs(m))
    return out


def find_covariance_peak(covariance, max_lag, search_range, filter_length,
                         positive_peak_only):
    """subroutines/find_covariance_peak.m.

    For two shifted (anti)correlated signals s1 and s2 this estimates the lag
    of s2 with respect to s1 from the maximum or minimum of xcov(s1, s2).
    """
    covariance = np.asarray(covariance, dtype=float)
    h = hamming(filter_length, sym=True)
    # 1-based (max_lag + 1 - search_range):(max_lag + 1 + search_range)
    lo = max_lag - search_range
    hi = max_lag + search_range + 1
    if positive_peak_only:
        filtered_cov = _filtfilt(h / h.sum(), covariance)
        idx = int(np.argmax(filtered_cov[lo:hi])) + 1     # MATLAB 1-based idx
        return -(idx - search_range - 1)
    # find maximum or minimum of filtered derivative
    dfdt = _filtfilt(np.ones(100) / 100, np.diff(covariance))
    dfdt_min_idx = int(np.argmin(dfdt[lo:hi])) + 1
    dfdt_max_idx = int(np.argmax(dfdt[lo:hi])) + 1
    filtered_cov = _filtfilt(h / h.sum(), covariance)
    if dfdt_max_idx < dfdt_min_idx:
        # maximum of derivative is to the left -> covariance is positive
        idx = int(np.argmax(filtered_cov[lo:hi])) + 1
    else:
        # maximum of derivative is to the right -> covariance is negative
        idx = int(np.argmin(filtered_cov[lo:hi])) + 1
    return -(idx - search_range - 1)


def nandetrend(x):
    """subroutines/nandetrend.m: remove a linear trend, NaN allowed."""
    x = np.asarray(x, dtype=float)
    t = np.arange(len(x), dtype=float)
    ok = np.isfinite(x)
    coeff = np.polyfit(t[ok], x[ok], 1)
    return x - np.polyval(coeff, t)


def step1_lag(w, c, max_lag, search_range, filter_length, detrend=True):
    """The per-interval lag of innFLUX_step1.m (lagtime1), in samples.

    w_prime = w - mean(w); c_prime = nandetrend(c) (or c - mean(c));
    NaN replaced by 0; cov_wc = xcov(w_prime, c_prime, MAX_LAG, 'unbiased');
    lag_samples = find_covariance_peak(cov_wc, MAX_LAG, LAG_SEARCH_RANGE,
                                       COVPEAK_FILTER_LENGTH, 0).
    Returns (lag_samples, cov_wc).
    """
    w = np.asarray(w, dtype=float)
    c = np.asarray(c, dtype=float)
    w_prime = w - np.nanmean(w)
    c_prime = nandetrend(c) if detrend else c - np.nanmean(c)
    w_prime = np.where(np.isfinite(w_prime), w_prime, 0.0)
    c_prime = np.where(np.isfinite(c_prime), c_prime, 0.0)
    cov_wc = xcov_unbiased(w_prime, c_prime, max_lag)
    return find_covariance_peak(cov_wc, max_lag, search_range, filter_length, 0), cov_wc
