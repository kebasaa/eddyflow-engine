# innFLUX 1.2.0 vs EddyFlow 8.2

This page compares EddyFlow with innFLUX stage by stage. It answers three questions:

- Which of innFLUX's capabilities did EddyFlow lack?
- Which of them were added?
- Where do the two packages differ although they intend the same thing?

The innFLUX-derived options in EddyFlow are all **off by default**. With them off, a project's
output is byte-identical to what it was before they existed.

- **innFLUX** is MATLAB code from the University of Innsbruck. It is described in Striednig, M.,
  Graus, M., Märk, T. D., and Karl, T. G.: InnFLUX – an open-source code for conventional and
  disjunct eddy covariance analysis of trace gas measurements: an urban test case, Atmos. Meas.
  Tech., 13, 1447–1465, https://doi.org/10.5194/amt-13-1447-2020, 2020.
  - The version examined is v1.2.0 (2023-01-14): `innFLUX_step0.m` (planar fit), `_step1.m` (per
    interval), `_step2.m` (cumulated lag), `subroutines/`, `tools/qclass` and `tools/driftfix`,
    about 2800 lines.
  - v1.2.0 added the `driftfix` tool (a routine by M. Graus, adapted by M. Striednig), which the
    paper does not describe.
- **EddyFlow** is Fortran, a maintained fork of LI-COR EddyPro 6.2.2.

innFLUX's terms ask that work using it cite Striednig et al. (2020). EddyFlow reimplements the
algorithms and copies no code. Every EddyFlow control, documentation page and module that derives
from innFLUX names its source.

## Stage by stage

| Stage | innFLUX | EddyFlow |
|---|---|---|
| Raw input | `.mat` files: daily sonic `w<yyyymmdd>.mat`, tracer files separate | ASCII, binary, GHG, TOB1, SLT; one file stream |
| Slow / disjunct tracer | step-function upsampling; DEC by zero-fill × factor | per-column rate (`instr_<K>_ac_freq`); covariance over valid pairs, which equals zero-fill × factor |
| Tilt correction | none, double rotation, or planar fit in ±15° windows at 1° steps | none, double, triple, planar fit (≤ 36 sectors), **moving windows (new, opt-in)** |
| Tilt uncertainty | bootstrap of α, β (Fig. 1) | **bootstrap (new, opt-in)** |
| Detrending | linear for T and tracers, mean for wind | block, linear, running mean, exponential |
| Lag, per period | Hamming-smoothed `xcov`, derivative sign rule (`lagtime1`) | covariance maximisation, optimizer, PWB, **smoothed peak search (new, opt-in)** |
| Lag, cumulated | sum of \|`xcov`\| over filtered periods per segment (`lagtime2`; Park et al. 2013) | **method 6 (new, opt-in)** |
| Clock drift | `driftfix`: rLOWESS and jump detection (`lagtime3`) | **`tlag_drift_meth` (new, opt-in)** |
| CCF on disk | daily `cov` files | **per-period `eddyflow_ccf` (new, opt-in)** |
| Unit conversion | constant molar volume at the period's T and p | full air-property chain, per-gas mole fraction or mixing ratio |
| WPL | dilution term from the H2O flux | Webb et al. 1980 / Ibrom et al. 2007, point-by-point mixing ratio |
| Spectral corrections | none (cospectra reported) | analytic and in-situ, many models |
| Random error | Finkelstein & Sims, Billesbach shuffle, Lenschow white noise | the same, plus Mann & Lenschow and Mahrt |
| Detection limit | far-lag noise ±160–180 s: mean, SD, RMSE, SNR (Spirig et al. 2005) | Wienhold, Langford windows; **mean, SD and SNR as output (new, opt-in)** |
| Quality | SST, ITC, flow-direction class, Foken 1–9 (`qclass`) | Foken 1–9, Mauder-Foken 0/1/2, Göckede, Vitale; wind-direction filter |
| θ fluxes | `wtheta`, `wtheta_v` | **`W_THETA`, `W_THETA_V` (new, opt-in)** |

## Gaps and what became of them

| innFLUX capability | Status | EddyFlow |
|---|---|---|
| Smoothed peak search with derivative sign rule | **implemented** | `covmax_smooth_len`, `covmax_peak_sign`; `m_lagpeak_core` |
| Cumulated-covariance lag per segment | **implemented** | `tlag_meth=6`, `tlag_cum_*`; `m_cumcov` |
| Clock-drift correction | **implemented** | `tlag_drift_*`, `gas_<i>_tlag_drift_ref`; `m_driftfix_core` |
| Per-period cross-covariance function | **implemented** | `gas_<i>_out_ccf`, `out_ccf_*` |
| w'θ', w'θv' | **implemented** | `out_theta_flux` |
| Far-lag noise mean, SD, SNR | **implemented** | `out_lod_noise_stats` |
| Moving-window planar fit | **implemented** | `pf_mw_num`, `pf_mw_width` |
| Bootstrap of tilt angles | **implemented** | `pf_boot_n`, `pf_boot_seed`; `m_pfboot` |
| Hundreds of PTR masses in one run | planned | MaxNumGases 64 today |
| Flow-direction class, combined class, dc/dt, completeness fraction | not adopted | judged of little value beside existing outputs |
| `.mat` raw input, separate tracer file stream | not adopted | convert to a format EddyFlow reads |

## Same intent, different implementation

- **Cross-covariance.**
  - innFLUX: `xcov(…, 'unbiased')` removes each series' mean once.
  - EddyFlow: `CrossCovarianceSeries` removes the means of the overlapping pairs at each lag.
  - On series with a flux, both give the same lag in every test case (`tests/innflux`, stage 1).
- **Smoothing at the ends of the function.**
  - innFLUX: MATLAB `filtfilt` pads the ends.
  - EddyFlow: computes the function past the search window, so only interior values, which are
    identical, are searched.
- **Clock drift relative to its median.**
  - innFLUX: shifts every function by the full smoothed lag, so the summed peak lands near zero,
    then adds that lag back.
  - EddyFlow: shifts by the drift relative to its median, so the peak stays inside the gas's
    search window. The result is the same.
- **Jump placement.**
  - innFLUX: MATLAB `risetime`.
  - EddyFlow: the midpoint of the first 10 % and 90 % crossings within two days.
- **Segments.**
  - Both place a period by its start time.
  - innFLUX's segments are MATLAB datenums; EddyFlow's are `yyyy-mm-ddTHH:MM`.
- **θ fluxes.** Both come from uncorrected covariances.
  - EddyFlow humidity-corrects the sonic temperature as level 2 H is (van Dijk et al. 2004) and
    writes θv with the specific humidity (0.61 q).
  - innFLUX uses 0.38 on the water mole fraction, the same factor.
- **Unit conversion and WPL.**
  - innFLUX converts with a constant molar volume and adds only the dilution term.
  - EddyFlow keeps its full treatment; not adopted.
- **Thresholds.** innFLUX hard-codes its filters (u* > 0.3, completeness 0.95/0.5, SST < 0.3). In
  EddyFlow they are settings, with innFLUX's values as defaults.

## Reproducing an innFLUX run in EddyFlow

| innFLUX parameter | EddyFlow |
|---|---|
| `TILT_CORRECTION_MODE = 2` | planar fit, `pf_mw_num=360`, `pf_mw_width=30` |
| bootstrap in step 0 | `pf_boot_n=1000` |
| `DETREND_*_SIGNAL = 1` | linear detrending |
| `LAG_SEARCH_RANGE` | the gases' time-lag optimization windows |
| `COVPEAK_FILTER_LENGTH = 15` | `covmax_smooth_len=15`, `covmax_peak_sign=2` (per period), or `tlag_cum_filter_len=15` |
| step 2 (`lagtime2`), `data_segments` | `tlag_meth=6`, `tlag_cum_segments` |
| driftfix (`lagtime3`) | `tlag_drift_meth=1`, reference gases `gas_<i>_tlag_drift_ref=1` |
| `cov/` files | `gas_<i>_out_ccf=1`, `out_ccf_max_lag_s` |
| `MET.wtheta`, `wtheta_v` | `out_theta_flux=1` |
| `flux_noise_mean/std`, `flux_SNR` | `lod_meth=1`, `lod_window=1` (±160–180 s), `out_lod_noise_stats=1` |
| `APPLY_WPL_CORRECTION` | `wpl_meth` (EddyFlow's own treatment) |
| `DISJUNCT_EC = 1` | the analyser's `instr_<K>_ac_freq` at its true rate |

## Verification

- `tests/innflux/compare.py` compares the engine with a line-by-line Python port of the MATLAB
  routines.
  - Peak search: 360 cases, 0 mismatches.
  - Drift: the engine agrees with the port to 1e-12 s and recovers an injected drift and jump.
- Fixtures: `base_covmax_innflux`, `base_innflux_outputs`, `base_cumcov`, `base_cumcov_drift`,
  `base_cumcov_par`, `base_ccf_out`, `base_pf_mw` and `base_pf_boot`.
- With every option off, the existing fixtures stay byte-identical.

## References

- Billesbach, D. P. (2011), Agric. For. Meteorol., 151, 394–405.
- Cleveland, W. S. (1979), J. Am. Stat. Assoc., 74, 829–836.
- Finkelstein, P. L. and Sims, P. F. (2001), J. Geophys. Res., 106, 3503–3509.
- Foken, T. and Wichura, B. (1996), Agric. For. Meteorol., 78, 83–105.
- Langford, B. et al. (2015), Atmos. Meas. Tech., 8, 4197–4213.
- Lenschow, D. H., Mann, J. and Kristensen, L. (1994), J. Atmos. Ocean. Tech., 11, 661–673.
- Park, J.-H. et al. (2013), Atmos. Chem. Phys., 13, 1439–1456.
- Spirig, C. et al. (2005), Atmos. Chem. Phys., 5, 465–481.
- Striednig, M., Graus, M., Märk, T. D. and Karl, T. G. (2020), Atmos. Meas. Tech., 13, 1447–1465.
- Stull, R. B. (1988), An Introduction to Boundary Layer Meteorology, Kluwer.
- Taipale, R., Ruuskanen, T. M. and Rinne, J. (2010), Atmos. Meas. Tech., 3, 853–862.
- Wienhold, F. G., Frahm, H. and Harris, G. W. (1994), J. Geophys. Res., 99, 16557–16567.
- Wilczak, J. M., Oncley, S. P. and Stage, S. A. (2001), Boundary-Layer Meteorol., 99, 127–150.
