# EddyUH 1.7b_COS vs EddyFlow 8.0.0

A stage-by-stage comparison of the two packages, written to answer three
questions: what does EddyUH do differently, which of its capabilities are
missing from EddyFlow, and which of those are worth adding.

Reference document. Tables first; prose only where an algorithm needs
explaining. File paths are given on both sides so every claim can be checked.

- **EddyUH** — MATLAB, University of Helsinki, © 2011, contacts I. Mammarella
  and O. Peltola. The tree examined is `EddyUH_1.7b_COS`, a carbonyl-sulfide
  fork: 145 `.m` files in eight module folders plus a function library.
- **EddyFlow** — Fortran, a maintained fork of LI-COR **EddyPro 6.2.2**
  (`README.md:10`; the second commit in the history is
  `c74a11b "Initial commit based on EddyPro v6.2.2"`).

The two lineages are independent. Where both implement "the same" method, the
implementations still differ often enough that a numerical comparison needs the
["Same intent, different implementation"](#same-intent-different-implementation)
section below to be interpretable.

## Read the code, not the menu

EddyUH's shipped state does not always match its own interface. This is not a
criticism of the package — it is a working scientific fork — but it does mean
the GUI label is not evidence of what ran:

- `EC_Software_FluxCalc/EddyUH_SC_Flux.m:281` implements the `lagwindow == 1`
  "standard: maximum of cross-covariance" branch as **the baseline-subtracted
  criterion**; the plain `max|cov|` line is commented out immediately above at
  line 283. Selecting "standard" in the GUI therefore runs the *alternative*
  method.
- `EddyUH_SC_Flux.m:323-331` unconditionally copies the CO2 time lag onto COS
  and CO. `EddyUH_SC_Flux2.m:324-341` instead applies the conditional
  Nemitz et al. (2018) rule. Which file the flux loop calls decides the COS lags.
- The LagOptimizer 3σ outlier filter is guarded by a comment reading
  "Comment if Time lag optimizer is not used", not by a setting.
- Large alternative blocks (a running-mean cross-covariance for the COS lag; the
  Nemitz rule inside `SC_Flux.m`) sit commented out in place.

## The CH-LAE configuration, decoded

Decoded from the supplied project files. Three are MATLAB v5 MAT-files
(`preproc_2603231128_CH-`, `lag_2603231128_CH-.10cl`,
`planar_fit_2603231128_CH-.1cl`); `resptime_2603231128_CH-.8cl` is plain text.
The `Ncl` suffix is the class count — 10 RH classes for the lag fit, 1 unbinned
planar fit, 8 RH classes for response time.

| Setting | EddyUH value |
|---|---|
| Sonic | Gill HS, 10 Hz, path 0.15 m, τ 0.016 s; AoA **off**, crosswind **off**, `headcorr=1`, `tiltcorr=0` |
| Analyser 1 | "Aerodyne cw-QCL" — CO2/H2O/N2O/**COS**; tube 60 m × 6 mm, 11.5 lpm, τ 0.1 s, horizontal separation 0.14 m, cell path 76 cm |
| Analyser 2 | LI-7200 — CO2/H2O; tube 1 m × 5.3 mm, 13.5 lpm, τ 0.06 s, separation 0.11 m |
| Site | z 47 m, altitude 689 m, lat 47.478, h_c 37 m, d 24.79 m, z0 5.55 m, boom 209° |
| Averaging | 30 min; FFT block m = 16384; Hamming window; 50 logarithmic bins |
| Rotation | `Dim_coordrot = 2` → double rotation (a `planar_fit_*.1cl` exists but is not selected) |
| Detrending | `DetrendType = 0` → block average |
| Despiking | `spi_method = 1` (consecutive-difference limit), per-variable `dlim`, `MaxNoSpikes = 6000` |
| Lag windows | QCL: CO2 16 ± 1 s, H2O 17 ± 5 s, N2O 16 ± 1 s, COS 16 ± 1 s. LI-7200: CO2 0.6 ± 0.5 s, H2O 4 ± 4 s |
| Lag climatology | CO2 16.015 ± 0.412, COS 16.018 ± 0.411, N2O 15.984 ± 0.425, H2O 18.42 ± 3.63 s; LI-7200 CO2 0.651 ± 0.127, H2O 4.523 ± 3.44 s. H2O additionally binned into 10 RH classes |
| Response times | `resptime_*.8cl` — **all NaN; never determined** |
| Corrections | `WPL_flag = n`, `diluflag = n`, `crosstalkflag = n`, `Fluxcalib_flag = n`, `detlim_flag = 0`, `unc_flag = 0` |
| Meteorology | external CSV supplying P, T_air, RH |

Three of those entries govern everything else:

1. **The response-time file is empty.** EddyUH's in-situ high-frequency
   correction (`flag_HF = 'yem'`) had no τ to work with, so the run used the
   theoretical closed-path branch `'ytc'` — Moncrieff et al. (1997) with
   Massman (1991) tube attenuation.
2. **WPL and dilution were off on both analysers**, as were the detection limit
   and the uncertainty estimates.
3. The physical instrument is a **MIRO MGA**. EddyUH has a fixed list of 17
   analysers (`inst_names`) with no MIRO, so it was configured as an Aerodyne
   cw-QCL; the June 2025 output header (`s2506_CH-c_all.flx`) reads
   `Fluxes.COS_cwQCL`. Tube geometry and τ were entered correctly, so the
   substitution is cosmetic — but EddyUH cannot name the instrument that took
   the data, and EddyFlow can.

**Consequence for this site:** any EddyFlow-vs-EddyUH difference is dominated by
the spectral correction and the time lag, not by the density corrections.

## Stage-by-stage

| Stage | EddyUH | EddyFlow | Verdict |
|---|---|---|---|
| Raw formats | Generic binary/ASCII, Edisol `.slt`, Campbell `.TOA`, Campbell `.TOB1` | LI-COR `.ghg` archive, generic ASCII/TOA, TOB1 (IEEE4/FP2), EddySoft SLT, EdiSol SLT, generic binary | **EddyFlow-only** (GHG) |
| Despiking | 3 methods: consecutive-difference limit, absolute limits, modified Vickers & Mahrt (1997) | Vickers & Mahrt (1997) or Mauder et al. (2013) MAD | **differs** — see below |
| Detrending | Block average, linear, McMillen (1988) autoregressive running mean | Block average, linear, running mean, exponential running mean | match |
| Rotation | Yaw / pitch / roll (double, triple); planar fit binned by sector, unbinned, or by stability class — the last of which its own rotation cannot apply | none / double / triple / planar fit / planar fit no-bias; sector-wise, and unbinned as a single 360° sector (`MaxNumWSect = 36`) | match in practice — see G |
| Time lag | Windowed cross-covariance; `max\|cov\|` or **baseline-subtracted**; 3σ rejection against a monthly climatology; RH-classed H2O window | constant / maxcov&default / maxcov / automatic optimisation (RH classes for H2O) / **PWB block-bootstrap** | **gap (E)** + EddyFlow-only (PWB) |
| Spectral, high-pass | Analytic per detrending mode (Aubinet et al. 2000) | Analytic (Moncrieff et al. 2004) | match |
| Spectral, low-pass | Moncrieff 97 closed/open, sonic-only, first-order in-situ τ, Horst 97 analytic | Moncrieff 97, Massman 00, Horst 97, Ibrom 07, Fratini 12, custom | **EddyFlow-only** (Massman, Ibrom, Fratini) |
| Correction application | **Iterative** with WPL and spectroscopy to < 1 % change, buoyancy flux recomputed each pass | Single pass | **differs** |
| Response time | In-situ ensemble-cospectral ratio, fitted per RH/WD/wind-speed class; H2O fitted `τ = a + b(RH/100)^c` | Ibrom-style in-situ cut-off frequency, RH-regressed (`src/src_fcc/fit_rh_to_cutoff.f90`) | match in purpose |
| WPL | Webb et al. (1980) open path; dilution-only closed path | Webb et al. (1980), switchable | match |
| Spectroscopic, closed path | **Peltola et al. (2014)** flux-level; **Chen et al. (2010)** point-by-point | — | **gap (B)** |
| Spectroscopic, open path | LI-7700 after McDermitt et al. (2010) | LI-7700 after McDermitt et al. (2011), `src/src_common/multipliers_7700.f90` | match |
| Crosstalk | Not a separate correction — see below | covered by the spectroscopic correction | **not a gap** |
| Sonic T, crosswind | Liu et al. (2001); coefficients for **Gill R2 and Metek USA-1 only** - any other sonic prints "No crosswind correction for this sonic anemometer model available", so the Gill HS at CH-LAE could not have been corrected even if asked | Liu et al. (2001), per-model coefficients (`src/src_common/cross_wind_corr.f90`) | EddyFlow broader |
| Sonic T, humidity | Schotanus (1983) / van Dijk et al. (2004) | van Dijk et al. (2004) eq. 3.53; FCC also offers Kaimal & Gaynor (1991) | match |
| Angle of attack | Nakai et al. (2006) | Nakai et al. (2006) **and** Nakai & Shimoyama (2012), auto-inferred; Gill w-boost handling | **EddyFlow-only** |
| Sonic head correction | Metek USA-1 3-D flow-distortion lookup tables | — | gap (I), niche |
| Raw QC tests | Vickers & Mahrt (1997) suite | Same suite, plus **KID** (Vitale et al. 2020), Fisher correlation-matrix, Qn scale, repeated-values R², wind-direction sector filter | **EddyFlow-only** |
| Flux QC flags | Foken et al. (2004) 1–9, ITC after Göckede et al. (2004) | Mauder & Foken (2004) 0-1-2, Foken (2003) 1–9, Göckede et al. (2006) 1–5 | match |
| Stationarity | Foken & Wichura (1996); Mahrt (1998) intermittency; Vickers & Mahrt relative flux error; Haar mean/variance | Foken et al. (2004); Mahrt (1998) via `RU_Mahrt_98`; Haar inside the discontinuities test | match in substance |
| Random uncertainty | **Billesbach (2011)** shuffle (20 realisations), Finkelstein & Sims (2001), Lenschow et al. (2000) noise | Finkelstein & Sims (2001), Mann & Lenschow (1994), Mahrt (1998). Billesbach present but **dead code** | **gap (D)** |
| Detection limit | **Wienhold et al. (1994)** | — | **gap (A)** |
| Footprint | Kormann & Meixner (2001); **urban morphometric + land-use extension** | Kljun et al. (2004), Kormann & Meixner (2001), Hsieh et al. (2000) | **gap (I)** + EddyFlow-only |
| Storage | — (eddy flux only) | One-point storage for H, LE and every gas | **EddyFlow-only** |
| Spectral output | Per-period spectra; ensemble mean (co)spectra; fittable cospectral models; peak frequency vs stability; ogives | Per-period and binned spectra/cospectra, binned ogives, ensemble by time and stability, fitted models, and six selectable analytic cospectra | closed (H); peak-frequency parameterisation still EddyUH-only |
| Partitioning | — | **Conditional Eddy Covariance** (Zahn et al. 2022), COS-capable | **EddyFlow-only** |
| Gas capacity | 17 named analysers, no MIRO | 64 gas slots, arbitrary analysers, multi-hygrometer, per-gas humidity reference | **EddyFlow-only** |

## Gaps — EddyUH has, EddyFlow does not

### A. Flux detection limit — Wienhold et al. (1994) — *implemented*

The standard deviation of the cross-covariance function measured **far from the
peak**, where there is no flux signal, taken as the noise floor of the
covariance. `EC_Software_Preproc/EddyUH_detlim_Preproc.m` evaluates it in two
50 s windows placed ±100 s from the gas's own lag and averages them.

Nothing equivalent exists in EddyFlow (`grep -ri "wienhold\|detlim" src/` returns
nothing).

**Why it matters here.** At CH-LAE the COS covariance is order 1e-8 in its raw
units and the flux median is −0.04 nmol m⁻² s⁻¹. Whether a given half-hour is
a measurement or noise is the first question anyone asks of a COS dataset, and
an EddyFlow run currently provides no way to answer it. It is also the
precondition for gap F.

### B. Closed-path spectroscopic and dilution correction — Peltola et al. (2014) / Chen et al. (2010) — *implemented*

Water vapour broadens the absorption lines of a laser analyser, so the reported
mixing ratio depends on humidity beyond simple dilution. EddyFlow implements
this for the **open-path LI-7700 only** (`src/src_common/multipliers_7700.f90`
with the lookup tables in `m_methane_tables.f90`). There is no equivalent for a
**closed-path** laser — which is what the MIRO MGA at CH-LAE is.

EddyUH offers both forms:

- **flux-level**, `EC_Software_FluxCalc/EddyUH_spect_CP.m`, after Peltola et al.
  (2014):

  ```
  cov ← (1−χ_q)/(1 + a·χ_q + b·χ_q²)
        · [ cov − (a + 2b·χ_q − b·χ_q² + 1)/((1 + a·χ_q + b·χ_q²)(1−χ_q)) · x̄ · w'q' ]
  ```

  with `w'q'` evaluated **at the time lag of the gas being corrected**, not at
  the water's own lag;
- **point-by-point**, `Functions_Library/dilucorr.m`, after Chen et al. (2010):
  `x ← x / (1 + a·χ_v + b·χ_v²)` applied to the raw series, with `a = −1, b = 0`
  meaning pure dilution and `a = b = 0` meaning none.

At CH-LAE both were switched off (`diluflag = n`), so this does not explain any
current EddyUH-EddyFlow difference — but it is the correction a COS/N2O QCL
dataset will eventually need, and EddyFlow cannot apply it at all.

#### The coefficients are not interchangeable with EddyUH's

Both packages divide by the same polynomial, `1 + a·χq + b·χq²`, with `χq` the
water mole fraction in mol/mol. They differ in what the polynomial is asked to
carry, and the difference is not visible in the formula.

**EddyUH folds the dilution in.** `Functions_Library/dilucorr.m` sets
`a = −1, b = 0` when no spectroscopic correction is wanted, giving a divisor of
`(1 − χq)` — the moist-to-dry conversion. Every coefficient in
`crosstalk_coeff_CP.m` therefore sits on a baseline of `a = −1`.

**EddyFlow does not.** Its WPL corrects the density separately, and more fully
than EddyUH's closed-path WPL, which carries no cell temperature or pressure
term. So `spectro_a = spectro_b = 0` is the identity and the coefficient carries
only the spectroscopic excess.

Mapping one onto the other, `D_EddyFlow = D_EddyUH / (1 − χq)`, which to second
order in `χq` is

```
spectro_a = a_EddyUH + 1
spectro_b = b_EddyUH + a_EddyUH + 1
```

and the pure-dilution pair `(−1, 0)` maps to `(0, 0)` exactly, as it must.

**Entering an EddyUH or Rella value unconverted counts the dilution twice.**
The published pairs, with the EddyFlow equivalents alongside:

| Analyser | Gas | EddyUH `a`, `b` | EddyFlow `spectro_a`, `spectro_b` |
|---|---|---|---|
| Aerodyne cw-QCL / dual-QCL | N2O | −1.39, 0 | **−0.39, −0.39** |
| Los Gatos RMT-200 | CH4 | −1.219, 1.678 | **−0.219, 1.459** |
| Los Gatos RMT-200 | CO2 | −1.219, 1.229 | **−0.219, 1.010** |
| Picarro G1301-f / G2301-f / G2311-f | CH4 | −0.9823, −2.39 | **0.0177, −2.3723** |
| Picarro G1301-f / G2301-f / G2311-f | CO2 | −1.2, −2.674 | **−0.2, −2.874** |
| Picarro G1301-f / G2301-f / G2311-f | H2O | 0.772, 2.525 | **1.772, 4.297** |
| anything else in EddyUH's table | — | −1, 0 | **0, 0** (no spectroscopic term) |

The right-hand column is derived here, not published: the `spectro_b` values in
particular come from a second-order expansion, worth about 7 × 10⁻⁴ in the
divisor at 15 mmol/mol. Treat them as a starting point, not as authority.

Note also that most of EddyUH's table is `−1, 0` — that is, **no spectroscopic
correction at all**, only dilution. Only the Los Gatos, the Picarro and the
Aerodyne N2O entries carry a real spectroscopic term.

### C. Analyser crosstalk — *withdrawn; not a separate capability*

An earlier draft of this document listed crosstalk as a gap, on the strength of
`EC_Software_FluxCalc/EddyUH_ctc_CP.m` and `crosstalk_coeff_CP.m` existing.
Following the call graph shows there is no such correction:

- **`EddyUH_ctc_CP.m` is never called.** Nothing in the package references it;
  the only other mention of the name is a stale copy-paste header at the top of
  `EddyUH_spect_CP.m`. It would not work if it were called — it computes
  `data + coeff .* datq` where `coeff` is a two-element vector and `datq` a
  scalar, which does not assign into a scalar slot.
- **`crosstalk_coeff_CP.m` is a table of spectroscopic coefficients.** Its only
  live callers are the two setup dialogs, which use it to seed
  `set_Gan.crosstalkcoeff` — and `set_Gan.crosstalkcoeff.(gas)(1)` and `(2)`
  are what `EddyUH_spect_CP.m` and `dilucorr.m` read as their `a` and `b`.
- **`crosstalkflag` gates the spectroscopic correction**, not a crosstalk one:
  `EddyUH.m:745` selects analysers on `crosstalkflag == 'y'` and calls
  `EddyUH_spect_CP`.

"Crosstalk" is a historical name for the spectroscopic coefficients in this
package. The capability is gap **B**, and is implemented.

What the table is genuinely worth is its **published coefficient values** —
see [the convention note](#the-coefficients-are-not-interchangeable-with-eddyuhs)
before using any of them.

### D. Billesbach (2011) random-shuffle uncertainty — *repaired and wired*

The idea: randomly permute the scalar series, recompute its covariance with w,
and repeat. Any covariance that survives shuffling is noise, so the spread over
realisations is the flux **noise floor**. EddyUH averages 20 realisations
(`EC_Software_Preproc/EddyUH_unc_Preproc.m`).

EddyFlow has `RIN_Billesbach_11` at `src/src_rp/random_error_handle.f90:305`,
marked `\todo Under development`. It is **never called** — the only other
occurrence of the name in the tree is a Billesbach (2012) citation at
`integral_turbulence_scale.f90:134`. As written it would also not work: `Set` is
declared `Set(N, M)` but indexed `Set(w, 1:N)` — row `w`, columns `1:N` — the
transpose of the convention used throughout the rest of that file, so it reads
out of bounds whenever `M < N`. And `ntimes = 1`.

Finishing it is a bug fix plus wiring.

### E. Baseline-subtracted cross-covariance lag detection — *implemented*

Instead of the largest `|cross-covariance|` in the search window, take the
largest **deviation of the cross-covariance from the straight line joining the
two window end points**.

The problem it solves: for a weak flux the cross-covariance function often has a
sloping baseline — from a trend, from a nearby stronger correlation, or from
insufficient averaging — and `max|cov|` then lands on whichever window edge the
baseline is highest at, not on the peak. Subtracting the chord removes the slope
and leaves the peak.

EddyUH: `EddyUH_SC_Flux.m:281`. **This is what CH-LAE actually ran**, because the
"standard" branch has been edited to use it (see
[Read the code, not the menu](#read-the-code-not-the-menu)).

EddyFlow's `CovMax` (`src/src_rp/timelag_handle.f90:496`) is plain `max|cov|`.
Its `covmax_stocdet` option detrends the *time series* before maximising — a
different remedy aimed at the same failure, applied one level up.

This is the single largest reason EddyFlow and EddyUH lags will differ at this
site.

### F. Detection-limit-conditional lag borrowing — Nemitz et al. (2018) — *implemented*

`EC_Software_FluxCalc/EddyUH_SC_Flux2.m:324` — use the tube-mate's lag when
`|max cov| < 3 × detection limit`, **or** when the maximum lands on either edge
of the search window. The reasoning: gases sharing a tube share a transport
delay, so a species whose own peak is not distinguishable from noise should
inherit the delay from one whose peak is.

EddyFlow already borrows across tube-mates, but only inside the PWB detector —
`S4_instrument_shared` and `S4_instrument_filled` in the reliability hierarchy at
`src/src_rp/pwb_timelag_handle.f90:550` — and it triggers on **HDI width**, a
measure of the bootstrap's own indecision, rather than on signal-to-noise. The
plain covariance-maximisation methods never borrow at all.

The two rules are complementary rather than redundant. Implementing F requires A.

### G. Stability-binned planar fit — *withdrawn: unreachable in EddyUH*

`EC_Software_PlanarFit/planar_fit_coeff.m` fits `w = b0 + b1·u + b2·v` per bin,
where the bin may be a **stability class** (`datcls = 0`), the whole dataset
(`datcls = 1`, written to the file as `nobin`), or a wind sector
(`datcls = 2`). Only the stability dimension is absent from EddyFlow, so that
is what this item was going to add.

**It cannot be applied.** `EddyUH_coordrot.m` takes the stability parameter as
its sixth argument and selects the matching `stabbin` row from it. Every one of
the four call sites in the package passes `NaN`:

| Caller | Line | Sixth argument |
|---|---|---|
| `EC_Software_FluxCalc/EddyUH.m` — the flux calculation | 579 | `NaN` |
| `EC_Software_Preproc/EddyUH_preproc.m` | 457 | `NaN` |
| `EC_Software_Common/EddyUH_plot_rawdata.m` | 487 | `NaN` |

With `stab = NaN`, `find(llim<=stab & ulim>=stab)` is empty and neither
`stab<llim(1)` nor `stab>ulim(end)` is true, so the plane coefficients `b` are
never assigned — and the next statement, `k(3) = 1./(1+b(2)^2+b(3)^2)`, reads
them. A project selecting stability bins would not produce different fluxes; it
would stop with an undefined variable. `nobin` and `WDbin` both assign `b`
unconditionally and are unaffected.

So there is no reference behaviour to reproduce. This is the same shape of
finding as the crosstalk item (C) and the cospectral-model selector (H): the
file exists, the mathematics is written down, and nothing reaches it.

**And it was not used here.** The site's planar-fit file,
`planar_fit_2603231128_CH-.1cl`, holds a single `nobin` fit for 2025-11-01 to
2026-03-15 — `w = 0.0845 − 0.3624·u + 0.0559·v` — with no `stabbin` or `WDbin`
field at all. The preprocessor setup then sets `Dim_coordrot = 2`, so CH-LAE
ran **double rotation** and that plane was never applied to a flux.

**Nothing is missing on the EddyFlow side either.** EddyUH's `nobin` is one
plane fitted to everything, which is what EddyFlow computes for a single wind
sector of 360° — `read_ini_rp.f90` forces `wsect_end(num_sec) = 360`, and
`WindSector` puts every direction in sector 1 when `num_sec = 1`. Sector and
unbinned modes are therefore both already available; only the stability
dimension is not, and it is the one EddyUH cannot use.

**If it is ever wanted anyway**, the obstacle is architectural rather than
arithmetic, and worth writing down. EddyFlow rotates at
`eddyflow-rp_main.f90:2171` and first computes `Ambient%zL` in `Fluxes0_rp` at
`:2536`, so the stability that would select the bin does not exist yet when the
bin must be chosen. The planar-fit assessment pass is worse off still: it keeps
only `Stats4%Mean(u:w)` per period (`:1281-1283`) and never forms a heat flux.
Both would need a provisional double-rotated stability computed for the purpose
— which is what EddyUH does, its `stab` coming from a preprocessor run under a
different rotation — and that is a design decision with no working precedent to
copy.

### H. Cospectral model library — *implemented*

EddyUH's `EC_Software_Spectral_Analysis/modelcospectra.m` returns six scalar
cospectra and two momentum ones. Five of the six are now selectable in EddyFlow
through `cosp_model`, defaulting to the curve it has always used:

| Model | Weighted form `n·Co/cov` | Stability branches |
|---|---|---|
| Moncrieff et al. (1997) — **default** | `12.92n/(1+26.7n)^{1.375}`, `4.378n/(1+3.8n)^{2.4}` above n = 0.54 (unstable); `n/(A+Bn^{2.1})` (stable) | yes |
| Kaimal et al. (1972) | `11n/(1+13.3n)^{7/4}`, `4n/(1+3.8n)^{7/3}` above n = 1 (unstable); `0.88(n/n₀)/(1+1.5(n/n₀)^{2.1})`, `n₀ = 0.23(1+6.4ζ)^{3/4}` (stable) | yes |
| Sakai et al. (2001) | `8n/(1+(11n)^{7/3})` | no |
| Su et al. (2003) | `10.4n/(1+10.6n)²` | no |
| Moraes et al. (2008) | `4.1n/(1+(6.6n)^{7/3})` | no |
| Kristensen et al. (1997) | `(2π/z)·n/(1+(2.4πn)^{2μ})^{7/6μ}`, μ = 0.23 | no |

Three things about this table are easy to get wrong and were, in the first draft
of the implementation plan:

- **EddyUH's `CMoore`, after Moore (1986), *is* EddyFlow's "Moncrieff et al.
  (1997)"** — term for term, scalar and momentum both. One cospectrum, two
  names. There is no seventh model hiding behind the second name.
- **There is no Mammarella et al. (2010) cospectral model in EddyUH.** The name
  appears in the file headers as an author contact and nowhere else.
- **`kaimal()` in `kaimal_models.f90` does not compute Kaimal's curve.** It
  computes Moore's, and feeds the `kaimal_cosp` column of the spectral
  assessment output. The column name is a misnomer that predates this work; the
  Kaimal (1972) row above is a genuinely different curve.

Also worth recording: `FCCsetup%SA%cosp` and `SA%cosp_type`, set from
`sa_cosp_type` in `read_ini_fcc.f90`, are **write-only** — nothing in the tree
reads either. The plan's original route for this item was to extend that
selector, which would have produced an option that did nothing.

The correction factor is `∫Co dn / ∫TF·Co dn`, so only the *shape* matters: any
constant scaling the cospectrum cancels. Kristensen's `2π/z` prefactor is
therefore dropped, and so is the numerically integrated normalisation EddyUH
applies to it — which depends on hard-coded indices into whatever frequency
vector it is handed, and could not change a correction factor either way.

The Reynolds stress keeps Moncrieff's momentum cospectrum whichever model is
selected; four of the five alternatives are scalar-only. On the analytic
regression fixture `Tau_scf` is identical across all six models, and the gas
correction factors run 1.0170 (Moraes) to 1.0557 (Kristensen) against
Moncrieff's 1.0335 — about four percent of the correction. Kaimal comes out at
1.0330, which is the expected answer: Moore's curve is a fit to Kaimal's data.

**Not carried across:** EddyUH's peak-frequency-versus-stability
parameterisation (`Functions_Library/peakFreq.m`) and the fitted-from-data forms
in `create_model_Cospectra.m`. The first shifts a cospectrum's peak rather than
choosing its form, and the four neutral models it would most help are precisely
the ones whose published fits carry no stability term to override. The second
overlaps EddyFlow's own `FitCospectralModel`.

### I. Urban footprint extension

`EC_Software_Footprint/urbanfpr_read_h.m`, `urbanfpr_process_h.m`,
`urbanfpr_read_lu.m`, `urbanfpr_process_lu.m` with
`Functions_Library/fpr_parameters.m`: compute **footprint-weighted displacement
height and roughness length from a building-height raster** via Macdonald et al.
(1998) (α = 4.43, β = 0.55, C_D = 1.2), recomputed for every averaging period and
fed back into the stability calculation (`EddyUH.m` ≈ lines 835–930); and
compute footprint-weighted land-use class fractions.

Genuinely absent from EddyFlow and genuinely novel. It needs raster I/O and a
per-period z0/d feedback loop. Out of scope for a forest site.

## EddyFlow-only capabilities

So the comparison is not one-directional. Relative to EddyUH, EddyFlow adds:

- **Formats** — LI-COR `.ghg` archives, Campbell TOB1 (IEEE4 and FP2), EddySoft
  SLT.
- **Spectral corrections** — Massman (2000, 2001) analytic; Ibrom et al. (2007)
  fully in-situ; Fratini et al. (2012) in-situ using measured w/Ts cospectra;
  Horst & Lenschow (2009) sensor-separation correction. The in-situ choice is
  made **per gas**, with individual fallback to Moncrieff for gases lacking a
  usable fit (`src/src_common/bpcf_bandpass_spectral_corrections.f90:92-151`).
- **Footprints** — Kljun et al. (2004) and Hsieh et al. (2000) in addition to
  Kormann & Meixner, with automatic demotion outside Kljun's validity envelope.
- **Corrections** — Burba et al. (2008) LI-7500 surface-heating terms; oxygen
  correction for krypton/Lyman-α hygrometers (van Dijk et al. 2003); analyser
  drift correction from calibration events.
- **Angle of attack** — Nakai & Shimoyama (2012) as well as Nakai et al. (2006),
  inferred automatically from anemometer model and firmware; Gill WindMaster
  w-boost firmware-bug handling.
- **QC** — KID (Vitale et al. 2020), Fisher correlation-matrix difference, Qn
  robust scale (Rousseeuw & Croux 1993), repeated-values R² test, and
  instantaneous wind-direction sector filtering (up to 16 excluded sectors).
- **Time lag** — the pre-whitening block-bootstrap detector (Vitale et al. 2024)
  with a documented reliability hierarchy and a re-readable half-hourly cache.
- **Partitioning** — Conditional Eddy Covariance (Zahn et al. 2022), with
  arbitrary extra species riding the w'/q'/c' octants, so COS partitions as a
  drop-in.
- **Storage** — one-point storage terms for H, LE and every gas.
- **Scale** — 64 gas slots against EddyUH's 17 fixed analyser models, multiple
  analysers, multiple hygrometers, per-gas humidity reference.
- **Biomet** — embedded, external file, or recursive external directory, with
  FLUXNET label and unit standardisation.

## Same intent, different implementation

The section to consult when the two packages disagree numerically.

### Despiking

EddyUH's Vickers & Mahrt implementation (`EC_Software_Common/EddyUH_despike.m`)
documents four deliberate deviations from the paper:

| | Vickers & Mahrt (1997) | EddyUH |
|---|---|---|
| Window advance | 1 point | **400 points** |
| σ multiplier increment | 0.1 | **0.3** ("0.3 is better") |
| Runs of ≥4 outliers | spikes | **not spikes** |
| Replacement | hold previous | **linear interpolation** |

Spike counts will not match `src/src_rp/test_spike_detection_vickers_97.f90`, and
neither will the despiked series. At CH-LAE the point is moot in a different way:
`spi_method = 1` selects the consecutive-difference limit, which has **no
EddyFlow equivalent at all**.

### Iterative vs single-pass correction

`EC_Software_FluxCalc/EddyUH_Spectralcorr.m` runs the spectral correction, WPL
and the spectroscopic correction in a loop until the flux changes by ≤ 1 %,
recomputing the buoyancy flux — and hence z/L, and hence the model cospectrum's
peak frequency — on every pass. EddyFlow applies each once
(`src/src_rp/fluxes0_rp.f90`, `fluxes1_rp.f90`, `fluxes23_rp.f90`).

The difference is largest in strongly non-neutral conditions, where the
correction factor is most sensitive to z/L.

### Tube attenuation has no COS coefficient

`EC_Software_FluxCalc/EddyUH_TF.m:155-194` hard-codes Massman (1991) attenuation
coefficients at Re = 2300:

| Gas | Λ |
|---|---|
| H2O | 14.03 |
| CO2, O3, N2O | 24.39 |
| CH4 | 14.73 |

There is **no branch for COS**, so `exist('lx','var')` fails and the tube term
`ttx` is set to unity. EddyUH's theoretical closed-path correction for COS
therefore omits tube attenuation entirely — through 60 m of tubing.

EddyFlow has no such hole, because it does not use a hard-coded table. In
`src/src_common/bpcf_analytic_transfer_functions.f90:95-120` it computes the tube
velocity from the stated flow rate and diameter, the air viscosity from air
temperature, and the actual Reynolds number, then branches:

- **laminar** (Re < 2300) — `exp(−t_tube·(π·d/2·f)² / (6·D_c))`, using the gas's
  own molecular diffusivity;
- **turbulent** (Re ≥ 2300) — Lenschow & Raupach (1991) eq. 7, which is
  gas-independent.

At CH-LAE the QCL tube carries 11.5 lpm through 6 mm, giving 6.8 m/s and
**Re ≈ 2700 — turbulent**. EddyFlow therefore applies the gas-independent branch
and needs no COS coefficient at all, while EddyUH, which assumes Re = 2300 and
looks COS up in a table that does not contain it, applies no tube term whatever.

This difference is not a matter of taste: EddyFlow is right, and EddyUH
under-corrects the COS flux at this site.

### Quality-flag composition

`EC_Software_FluxCalc/EddyUH_QualityFlags.m` maps both the flux stationarity and
the integral turbulence characteristics onto 1–9 using the standard
15/30/50/75/100/250/500/1000 % thresholds, but the **overall flag combines the
stationarity with the ITC of w only**. EddyFlow's `qc_flags_subs.f90` follows the
Handbook of Micrometeorology table. Flags will differ for periods where the
scalar ITC is the binding constraint.

### Sensor separation

EddyUH stores a single scalar `Hsensor_separ` per analyser. EddyFlow stores
northward, eastward and vertical separations. A vertical offset that EddyFlow
accounts for cannot be expressed in EddyUH at all.

### Spectral estimation

Both taper with a Hamming window over a power-of-two block and bin
logarithmically. EddyUH uses the normalisation factor 2.52 from Kaimal & Finnigan
(1994, p. 267) (`EC_Software_Preproc/calculate_save_spectra.m`). At 10 Hz over
30 min both use 16384 samples and 50 bins, so this stage matches closely.

## Recommendations

| Gap | Verdict | Reasoning |
|---|---|---|
| **A** — Wienhold detection limit | **done** | `detlim_meth`, off by default. Reported in covariance units in the full output and FLUXNET files. |
| **B** — closed-path spectroscopic | **done** | `spectro_meth`, off by default, point-by-point on the raw series. Per-gas `col_N_spectro_a/b` in the metadata. |
| **C** — Rella crosstalk | **withdrawn** | Not a separate capability. EddyUH's "crosstalk" coefficients *are* the spectroscopic ones, and its `EddyUH_ctc_CP.m` is unreachable dead code. Covered by B. |
| **D** — Billesbach | **done** | `ru_meth = 4`, one more choice on the existing random-uncertainty method. 20 realisations, EddyUH’s number. |
| **E** — baseline-subtracted lag | **done** | `covmax_debaseline`, off by default. A modifier on covariance maximisation, not a sixth method. |
| **F** — conditional lag borrowing | **done** | `tlag_borrow_meth`, off by default, and refused without A. Complements, does not duplicate, the PWB borrowing. |
| **G** — stability-binned planar fit | **withdrawn** | EddyUH never passes a stability to the rotation, so its own `stabbin` branch leaves the plane undefined. CH-LAE used `nobin` and double rotation. Sector and unbinned modes are already available here. |
| **H** — cospectral models | **done** | `cosp_model`, five alternatives, Moncrieff stays the default. |
| **I** — urban footprint | document, **skip** | The only item needing both a new raster input format and a per-period feedback of z0/d into the stability calculation. Out of scope for a forest site. |

Six of the nine were ported and two withdrawn on evidence; the ninth is out of
scope. **Every ported capability is optional and off by default**, and
EddyFlow's own behaviour and defaults are unchanged - each one was checked
byte-identical against the previous engine with the new setting absent. When one
is enabled, its parameters default to EddyUH's own values, and every new control
carries a tooltip naming the method and its citation.

The two withdrawals were not judgement calls about worth. In both cases the
capability exists in EddyUH's source and **cannot be reached at run time**: the
crosstalk routine (C) has no caller, and the stability-binned planar fit (G) is
selected by a stability that every call site passes as `NaN`. The
spectral-assessment cospectral-model selector met the same fate on the way into
H. Following the call graph rather than the file names is what separated the
three of them from real gaps.

The house rules for adding a setting - which generator owns which table, which
INI group a key must sit in, and what fails silently if either is wrong - are in
[`adding-a-setting.md`](adding-a-setting.md).

## Reproducing an EddyUH run in EddyFlow

`EF_730b_COS_LAE_2025_EddyUH_comparison.eddyflow` and
`COS_LAE_2025_EddyUH_comparison.metadata` configure EddyFlow to match the
decoded CH-LAE settings as closely as the current code allows, with no source
changes. The choices that matter:

| EddyFlow key | Value | EddyUH basis |
|---|---|---|
| `rot_meth` | 1 — double rotation | `Dim_coordrot = 2` |
| `detrend_meth` | 0 — block average | `DetrendType = 0` |
| `avrg_len`, `power_of_two`, `tap_win`, `nbins` | 30, 1, 3 (Hamming), 50 | `ave_t = 30`, `m = 16384`, `calculate_save_spectra.m` |
| `tlag_meth` | 2 — maxcov & default | windowed covariance maximisation with fallback |
| `hf_meth` | 1 — Moncrieff et al. (1997) | `'ytc'`; the in-situ path was unusable (empty `resptime`) |
| `lf_meth` | 1 — analytic | block-average high-pass transfer function |
| `wpl_meth` | 0 | `WPL_flag = n`, `diluflag = n` |
| `qc_meth` | 2 — Foken 1–9 | `EddyUH_QualityFlags.m` |
| `foot_meth` | 2 — Kormann & Meixner (2001) | EddyUH's only footprint model |
| `cec_meth`, `ru_meth`, `cross_wind`, `bu_corr` | 0 | not present or switched off in EddyUH |
| `flow_distortion` | 0 — no AoA correction | `Aoa_corr = 'n'` |
| `gill_wm_wboost` | 0 | no EddyUH equivalent (inert on an HS-50 regardless) |
| `despike_vm` | 1 — Vickers & Mahrt | nearest relative to `spi_method = 1`; **not a match** |
| `out_mean_cosp`, `out_mean_spec` | 0 | either one sets `fcc_follows` (`eddyflow-rp_main.f90:266`), which hands `full_output` to FCC. Off keeps the run a single RP pass. Per-period binned cospectra and ogives stay on |
| `col_*_{nom,min,max}_timelag` | per the table above | `lags ± dlags` from `preproc_*` |

Two values in the metadata were wrong independently of this comparison and are
corrected in the new file: the Gill HS path length was `1.0` where metadata path
lengths are in **cm** and the instrument is 0.15 m (`read_metadata_file.f90:242`),
and the LI-7200 tube length was `1.0` where tube lengths are in **cm** and the
tube is 1 m (`read_metadata_file.f90:393`).

A one-day trial of that project reproduces EddyUH's lag population closely — COS
time lags span 15.1–16.9 s with a median of 16.0 s, against the EddyUH
climatology of 16.018 ± 0.411 s.

**Differences that will remain, and should be quantified rather than chased:**

- despiking (`spi_method = 1` has no equivalent);
- the lag criterion, until gap E is implemented;
- single-pass versus iterative correction;
- `MaxNoSpikes = 6000`, which has no EddyFlow analogue;
- the RH-classed H2O lag window — EddyFlow's only RH-classed lag path is
  `tlag_meth = 4`, a different algorithm requiring a prior optimisation run;
- tube attenuation for COS, which EddyUH omits and EddyFlow applies.

## References cited by the two packages

Aubinet et al. (2000) · Billesbach (2011) · Burba et al. (2008) · Chen et al.
(2010) · Finkelstein & Sims (2001) · Foken & Wichura (1996) · Foken et al. (2004)
· Fratini et al. (2012) · Göckede et al. (2004, 2006) · Horst (1997) · Horst &
Lenschow (2009) · Hsieh et al. (2000) · Ibrom et al. (2007) · Kaimal et al.
(1972) · Kaimal & Finnigan (1994) · Kaimal & Gaynor (1991) · Kaimal & Kristensen
(1991) · Kljun et al. (2004) · Kormann & Meixner (2001) · Kristensen et al.
(1997) · Lenschow & Raupach (1991) · Lenschow et al. (2000) · Liu et al. (2001) · Macdonald et al. (1998) ·
Mahrt et al. (1998) · Mammarella et al. (2010) · Mann & Lenschow (1994) · Massman
(1991, 1998, 2000, 2001) · Mauder et al. (2013) · Mauder & Foken (2004) ·
McDermitt et al. (2010, 2011) · McMillen (1988) · Moncrieff et al. (1997, 2004) ·
Moore (1986) · Moraes et al. (2008) · Nakai et al. (2006) · Nakai & Shimoyama
(2012) · Nemitz et al. (2018) · Peltola et al. (2014) · Rannik & Vesala (1999) ·
Rella (2010) · Rousseeuw & Croux (1993) · Sakai et al. (2001) · Schotanus et al.
(1983) · Su et al. (2003) · van Dijk et al. (2003, 2004) · Vickers & Mahrt (1997)
· Vitale et al. (2020, 2024) · Webb, Pearman & Leuning (1980) · Wienhold et al.
(1994) · Wilczak et al. (2001) · Zahn et al. (2022)
