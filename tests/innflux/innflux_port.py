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


# --------------------------------------------------------------------------
# tools/driftfix/innFLUX_driftfix.m (innFLUX v1.2.0, 2023; based on a routine
# by M. Graus, adapted by M. Striednig; not in Striednig et al. 2020).
# --------------------------------------------------------------------------

def rlowess(x, y, span, niter=5):
    """MATLAB smooth(x, y, span, 'rlowess'): local linear fits over the span
    nearest points in index order, tricube weights on the distance in x
    scaled by the window's largest, then niter bisquare robustness passes on
    residuals scaled by 6 * median |residual| (Cleveland 1979)."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    n = len(x)
    if n < 3:
        return y.copy()
    k = max(2, min(int(span), n))
    rw = np.ones(n)
    ys = y.copy()
    for it in range(niter + 1):
        for i in range(n):
            lo = i - (k - 1) // 2
            hi = lo + k - 1
            if lo < 0:
                lo, hi = 0, k - 1
            if hi > n - 1:
                hi, lo = n - 1, n - k
            xs = x[lo:hi + 1]
            dmax = max(abs(x[i] - x[lo]), abs(x[hi] - x[i]))
            if dmax > 0:
                wt = (1 - (np.abs(xs - x[i]) / (dmax * 1.0000001)) ** 3) ** 3
            else:
                wt = np.ones(len(xs))
            wt = wt * rw[lo:hi + 1]
            sw = wt.sum()
            if sw <= 0:
                ys[i] = y[i]
                continue
            sx = (wt * xs).sum()
            sy = (wt * y[lo:hi + 1]).sum()
            sxx = (wt * xs * xs).sum()
            sxy = (wt * xs * y[lo:hi + 1]).sum()
            den = sw * sxx - sx * sx
            if abs(den) <= 1e-12 * max(1.0, abs(sw * sxx)):
                ys[i] = sy / sw
            else:
                ys[i] = (sy * sxx - sx * sxy) / den + (sw * sxy - sx * sy) / den * x[i]
        if it == niter:
            break
        res = y - ys
        mad = np.median(np.abs(res))
        if mad <= 0:
            break
        r = np.abs(res) / (6 * mad)
        rw = np.where(r < 1, (1 - r ** 2) ** 2, 0.0)
    return ys


def driftfix(px, py, t, span1d, span, thr):
    """The drift at times t from reference lags (px, py), innFLUX's steps:
    one-day rlowess, jumps where y(t) - y(t + 1 day) < thr placed at the
    midpoint of the 10 %/90 % rise within two days (MATLAB risetime stands
    there), then rlowess with the long span inside each stretch."""
    px = np.asarray(px, float)
    py = np.asarray(py, float)
    t = np.asarray(t, float)
    n = len(t)
    s1 = rlowess(px, py, span1d)
    y1 = np.interp(t, px, s1)
    rising = np.zeros(n, bool)
    for i in range(n):
        if t[i] + 1 <= t[-1]:
            rising[i] = y1[i] - np.interp(t[i] + 1, t, y1) < thr
    jumps = []
    if not (rising[0] or rising[-1]):
        i = 0
        while i < n - 1:
            if not rising[i] and rising[i + 1]:
                s = i + 1
                e = s
                while e < n - 1 and t[e + 1] <= t[s] + 2:
                    e += 1
                seg = y1[s:e + 1]
                l10 = seg.min() + 0.1 * (seg.max() - seg.min())
                l90 = seg.min() + 0.9 * (seg.max() - seg.min())
                lt = s + int(np.argmin(seg))
                while lt < e and y1[lt] < l10:
                    lt += 1
                ut = lt
                while ut < e and y1[ut] < l90:
                    ut += 1
                jumps.append((lt + ut + 3) // 2 - 1)      # the engine's 1-based (lt+ut+1)/2
                while i < n - 1 and rising[i + 1]:
                    i += 1
            i += 1
    drift = np.full(n, np.nan)
    bounds = [0] + jumps + [n]
    for k in range(len(bounds) - 1):
        a, b = bounds[k], bounds[k + 1]
        if b <= a:
            continue
        t_lo = -np.inf if k == 0 else t[a]
        t_hi = np.inf if k == len(bounds) - 2 else t[b]
        m = (px >= t_lo) & (px < t_hi)
        if not m.any():
            continue
        sw = rlowess(px[m], py[m], span) if m.sum() >= 3 else py[m]
        drift[a:b] = np.interp(t[a:b], px[m], sw)
    return drift, len(jumps)
