"""The engine's PWB against dyco's, number for number, on the same resamples.

static_checks/test_pwb_reference_static.py pins the deterministic half of the
chain - the unit-root test, the AR fits, the pre-whitened CCF - against RFlux's
frozen output. It deliberately stops at the bootstrap, because the engine and
R share no random stream. This goes the rest of the way: dyco is run with a
stub generator that records the block starts it draws, the engine's chain
(obj/<os>/pwb_compare, built with `mingw32-make pwbcmp`) is run on exactly those
starts, and every quantity both produce is compared:

  * the unit-root decision, the three AR orders and first coefficients, the
    full-data pre-whitened CCF peak;
  * per combination (cw, wc, ct, tc): the 99 replicate peak lags, exactly; the
    mean smoothed CCF over the search window; the 95 % HDI of the replicate
    lags, exactly; the mode; the |mean smoothed CCF| at the mode;
  * the chosen combination and the final lag and HDI.

The mode is dyco's estimator in both (m_pwb_core MapLagEstimate), except that
dyco jitters the lags with N(0, 1e-4) from its random stream and the engine
has no stream there - so dyco is handed zero jitter and the modes must agree
exactly. One difference is known and reported rather than failed, because it
is a documented departure in m_pwb_core.f90:

  * PINNED - the engine prefers a combination whose mode is not on the window
    edge before comparing magnitudes; dyco compares magnitudes only.

Any other mismatch fails (exit status 1). Before the engine took dyco's mode
estimator this found the final lag a record apart in 6 of 40 real periods.

Inputs: dyco's two synthetic reference fixtures (stationary and differencing
branch, 20 Hz) and real 30-min CH-LAE periods from tests/regression/data_gappy
(10 Hz; LI-7200 and MIRO gases, the MIRO columns with their gaps), wind
double-rotated with dyco's own WindDoubleRotation. Read only.

Needs a Python with numpy, scipy, pandas and matplotlib (dyco imports it) - on
this machine the conda env `dp` - and dyco checked out as dyco-main beside
this repository's root.

    <dp python> tests/pwb_dyco/compare.py [--files N] [--keep DIR]

Part of the EddyFlow engine's test tools.
"""

import argparse
import gzip
import io
import math
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
DYCO = ROOT / "dyco-main"
sys.path.insert(0, str(DYCO))

import dyco.pwb  # noqa: E402
from dyco.pwb import PreWhiteningBootstrap, PwbBatchDetection  # noqa: E402
from dyco.rotation import WindDoubleRotation  # noqa: E402

#: dyco warns when the block is shorter than twice its search half-width. Here
#: that half-width is widened to cover the engine's guard band, so the warning
#: says nothing about the block the two actually use - which is compared.
dyco.pwb.warn = lambda *args, **kwargs: None

DRIVER = ROOT / "obj" / "win" / "pwb_compare.exe"
if not DRIVER.exists():
    DRIVER = ROOT / "obj" / "linux" / "pwb_compare"

COMBOS = ("cw", "wc", "ct", "tc")
#: The engine's defaults (read_ini_rp.f90), which are also dyco's.
N_BOOT = 99
SMOOTH = 5
BLOCK_S = 20.0
REL = 1e-9
#: The normalised CCF lies in [-1, 1]. The engine sums it directly, dyco by
#: FFT; over 18000 records the two agree to a few 1e-13 absolute, which on a
#: CCF value near 1e-4 is a relative difference of 1e-9 and more.
CCF_ABS = 1e-11

GAPPY = ROOT.parent / "eddyflow-engine" / "tests" / "regression" / "data_gappy"
if not GAPPY.exists():
    GAPPY = ROOT / "tests" / "regression" / "data_gappy"

#: Real gases and their windows [s]: the physical, positive tube delays the
#: CH-LAE fixtures search, plus one symmetric window to exercise negative lags.
REAL_GASES = (
    ("LI72_CO2_DRY", 0.0, 10.0),
    ("LI72_H2O_DRY", 0.0, 15.0),
    ("CO2_DRY", 0.0, 15.0),
    ("H2O", 0.0, 25.0),
    ("N2O_DRY", 0.0, 15.0),
    ("COS_DRY", -10.0, 10.0),
)


class RecordingRng:
    """Stands in for dyco's generator: draws as numpy would, and keeps the
    block starts so the engine can be handed the same ones."""

    def __init__(self, seed):
        self.g = np.random.default_rng(seed)
        self.starts = []

    def integers(self, low, high, size):
        a = self.g.integers(low, high, size=size)
        self.starts.append(np.array(a))
        return a

    def normal(self, loc, scale, size):
        #> dyco's only use: the MAP jitter. The engine has none.
        return np.zeros(size)


class RecordingPwb(PreWhiteningBootstrap):
    """dyco's class, keeping what its run() computes per combination."""

    def _run_combination_bootstrap(self, x_pw, y_pw):
        out = super()._run_combination_bootstrap(x_pw, y_pw)
        self.combos_seen = getattr(self, "combos_seen", []) + [out]
        return out

    def _fit_ar_model(self, x):
        phi, p = super()._fit_ar_model(x)
        self.phis = getattr(self, "phis", []) + [phi]
        return phi, p


def engine_geometry(n, hz, min_rl, max_rl):
    """PwbDetectGas's evaluated range and RunPwbCombination's block length."""
    trail = max(1, SMOOTH) // 2
    margin = max(trail, int(math.floor(2.0 * hz + 0.5)))
    eval_lo = max(min_rl - margin, -(n - 3))
    eval_hi = min(max_rl + margin, n - 3)
    widest = max(abs(min_rl), abs(max_rl))
    block_req = int(math.floor(BLOCK_S * hz + 0.5))
    block = max(block_req, 2 * widest)
    return eval_lo, eval_hi, block_req, block


def run_case(name, w, t, s, hz, lws, uws, seed, keep):
    n = len(w)
    min_rl = int(round(lws * hz))
    max_rl = int(round(uws * hz))
    eval_lo, eval_hi, block_req, block = engine_geometry(n, hz, min_rl, max_rl)
    lag_max_rec = max(abs(eval_lo), abs(eval_hi))

    df = pd.DataFrame({"w": w, "s": s, "t": t})
    pwb = RecordingPwb(df, var_w="w", var_scalar="s", var_tsonic="t", hz=hz,
                       lag_max_s=lag_max_rec / hz, n_bootstrap=N_BOOT,
                       block_length_s=block / hz, wdt=SMOOTH,
                       lws=min_rl / hz, uws=max_rl / hz, segment_name=name)
    rng = RecordingRng(seed)
    pwb._rng = rng
    pwb.run()
    assert len(rng.starts) == 4, "dyco drew block starts %d times" % len(rng.starts)
    n_eff = pwb._n_eff
    nblocks = rng.starts[0].shape[1]

    buf = io.StringIO()
    buf.write("%d %d %d %d %d %d %d %d\n" % (n, hz, min_rl, max_rl, SMOOTH,
                                            block_req, nblocks, N_BOOT))
    for a, b, c in zip(w, t, s):
        buf.write("%s %s %s\n" % tuple("-9999" if not np.isfinite(v) else repr(float(v))
                                       for v in (a, b, c)))
    for st in rng.starts:
        for row in st:
            buf.write(" ".join(str(int(v) + 1) for v in row) + "\n")
    with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False,
                                     dir=keep if keep else None) as f:
        f.write(buf.getvalue())
        inpath = f.name
    run = subprocess.run([str(DRIVER), inpath], capture_output=True, text=True)
    if not keep:
        Path(inpath).unlink()
    if run.returncode != 0:
        return [("FAIL", "driver: " + run.stdout.strip() + run.stderr.strip())], None
    eng = {}
    for line in run.stdout.splitlines():
        k, _, v = line.partition("=")
        eng[k.strip()] = v.strip()

    out = []

    def fail(msg):
        out.append(("FAIL", msg))

    def note(kind, msg):
        out.append((kind, msg))

    def close(a, b):
        return abs(a - b) <= REL * max(abs(a), abs(b), 1e-300)

    # The deterministic half.
    if int(eng["n_eff"]) != n_eff:
        fail("n_eff %s vs %d" % (eng["n_eff"], n_eff))
    if (eng["differenced"] == "T") != (n_eff == n - 1):
        fail("unit-root decision differs")
    for i, key in enumerate(("scalar", "w", "t")):
        dkey = {"scalar": "scalar", "w": "w", "t": "tsonic"}[key]
        if int(eng["ar_order_" + key]) != pwb._ar_orders[dkey]:
            fail("AR order %s: %s vs %d" % (key, eng["ar_order_" + key], pwb._ar_orders[dkey]))
        phi = pwb.phis[i]
        if len(phi) and not close(float(eng["phi1_" + key]), float(phi[0])):
            fail("phi1 %s: %s vs %r" % (key, eng["phi1_" + key], float(phi[0])))
    if int(eng["tlag_pw"]) != pwb._tlag_pw_records:
        fail("full-data PW peak %s vs %d" % (eng["tlag_pw"], pwb._tlag_pw_records))
    if not close(float(eng["corr_pw"]), pwb._corr_pw):
        fail("corr_pw %s vs %r" % (eng["corr_pw"], pwb._corr_pw))
    if int(eng["block_len"]) != pwb._block_length_records:
        fail("block length %s vs %d" % (eng["block_len"], pwb._block_length_records))

    # The bootstrap, combination by combination.
    lr = pwb._lag_max_records
    modes_equal = True
    for c, combo in zip(COMBOS, pwb.combos_seen):
        elags = [int(v) for v in eng["lags_" + c].split()]
        dlags = [int(v) for v in combo["lags"]]
        if elags != dlags:
            bad = sum(1 for a, b in zip(elags, dlags) if a != b)
            fail("%s replicate lags differ in %d of %d" % (c, bad, len(dlags)))
            continue
        ems = [float(v) for v in eng["mean_smooth_" + c].split()]
        for lag in range(min_rl, max_rl + 1):
            a = ems[lag - eval_lo]
            b = float(combo["mean_smooth_ccf"][lag + lr])
            if abs(a - b) > CCF_ABS:
                fail("%s mean smoothed CCF at lag %d: %r vs %r" % (c, lag, a, b))
                break
        dlo, dhi = PreWhiteningBootstrap._hdi(np.array(dlags) / hz, 0.95)
        if float(eng["hdi_lo_" + c]) != dlo or float(eng["hdi_hi_" + c]) != dhi:
            fail("%s HDI [%s, %s] vs [%r, %r]" % (c, eng["hdi_lo_" + c],
                                                  eng["hdi_hi_" + c], dlo, dhi))
        emode = int(eng["mode_" + c])
        if emode != int(combo["mode_lag"]):
            modes_equal = False
            fail("%s mode %d vs dyco %d" % (c, emode, combo["mode_lag"]))
        at = abs(float(combo["mean_smooth_ccf"][emode + lr]))
        if abs(float(eng["ccf_at_mode_" + c]) - at) > CCF_ABS:
            fail("%s |CCF| at mode %s vs %r" % (c, eng["ccf_at_mode_" + c], at))

    any_pinned = any(eng["ok_" + c] == "F" for c in COMBOS)
    if eng["best"] != pwb._best_combination:
        why = []
        if not modes_equal:
            why.append("MODE")
        if any_pinned:
            why.append("PINNED")
        if why:
            note("+".join(why), "chose %s vs dyco %s" % (eng["best"], pwb._best_combination))
        else:
            fail("chose %s vs dyco %s with equal modes and nothing pinned"
                 % (eng["best"], pwb._best_combination))
    elif int(eng["mode_" + eng["best"]]) != pwb._tlag_records and modes_equal:
        fail("final lag %s vs %d" % (eng["mode_" + eng["best"]], pwb._tlag_records))

    summary = (eng["best"], int(eng["mode_" + eng["best"]]),
               float(eng["hdi_hi_" + eng["best"]]) - float(eng["hdi_lo_" + eng["best"]]),
               pwb._best_combination, pwb._tlag_records, pwb._hdi_hi_s - pwb._hdi_lo_s,
               eng["differenced"])
    return out, summary


def ladder_stage(table_path):
    """Stage 2: the engine's settled classes against dyco's PWBOPT ladder.

    The engine's half-hourly table carries, per period and gas, the raw
    detection (lag, HDI range, edge-pinned, pre-filtered) and what
    PostProcessPwbTimelagCache settled. The raw columns go through dyco's
    apply_hdi_prefilter and apply_pwbopt with the table's own thresholds; then
    every period dyco calls S1 or S2 must be S1 or S2 in the engine with the
    same lag, and vice versa, and where both carry a lag forward it must be the
    same lag. The engine's further fills - interpolation between settled
    periods, the carry limit in hours, sharing within an analyser, the terminal
    fallback - are its documented departures (m_pwb_core.f90 header,
    PostProcessPwbTimelagCache) and are counted, not failed.
    """
    lines = Path(table_path).read_text(encoding="utf-8").splitlines()
    fp = next(l for l in lines if l.startswith("fingerprint="))
    def setting(key):
        part = fp.split(key + "=")[1]
        return float(part.split("_")[0].split(":")[0])
    hdi_t, dev_t, pre_t = setting("hdi"), setting("dev"), setting("prefilter")
    start = lines.index("data") + 1
    df = pd.read_csv(io.StringIO(chr(10).join(lines[start:])), skipinitialspace=True)
    df["stamp"] = df["date"].str.strip() + " " + df["time"].str.strip()
    fails, counts = [], {}
    for gas, g in df.groupby("gas", sort=False):
        g = g.sort_values("stamp", kind="stable")
        raw = g["raw_lag_s"].astype(float).where(g["edge_pinned"].str.strip() != "T")
        raw = raw.where(raw > -9998)
        hdi = g["hdi_range_s"].astype(float)
        hdi = hdi.where(hdi > -9998)
        tl = PwbBatchDetection.apply_hdi_prefilter(raw.to_numpy(), hdi.to_numpy(),
                                                   threshold=pre_t)
        lad = PwbBatchDetection.apply_pwbopt(np.asarray(tl, float), hdi.to_numpy(),
                                             hdi_thresh=hdi_t, dev_thresh=dev_t)
        for (_, row), (_, d) in zip(g.iterrows(), lad.iterrows()):
            ecls = row["reliability_class"].strip()
            eown = ecls in ("S1_optimal", "S2_optimal")
            down = d["flag"] in ("S1_optimal", "S2_optimal")
            used = float(row["used_lag_s"])
            key = "%s/%s" % (d["flag"], ecls)
            counts[key] = counts.get(key, 0) + 1
            where = "%s %s" % (row["stamp"], gas.strip())
            if eown != down or (eown and ecls != d["flag"]):
                fails.append("%s: engine %s, dyco %s" % (where, ecls, d["flag"]))
            elif eown and abs(used - d["pwbopt_s"]) > 1e-6:
                fails.append("%s: settled lag %r vs dyco %r" % (where, used, d["pwbopt_s"]))
            elif (ecls == "S3_carryforward" and np.isfinite(d["pwbopt_s"])
                    and abs(used - d["pwbopt_s"]) > 1e-6):
                fails.append("%s: carried %r vs dyco %r" % (where, used, d["pwbopt_s"]))
    print("stage 2 (%s): hdi<%g, dev<=%g, prefilter %g" % (Path(table_path).name, hdi_t, dev_t, pre_t))
    for key in sorted(counts):
        print("  dyco/engine %-40s %d" % (key, counts[key]))
    for f in fails:
        print("  FAIL", f)
    return len(fails)


def synthetic_cases():
    for label in ("stationary", "differencing"):
        path = DYCO / "tests" / "data" / ("pwb_reference_%s.csv.gz" % label)
        with gzip.open(path, "rt") as f:
            df = pd.read_csv(f)
        yield ("dyco_" + label, df["w"].to_numpy(float), df["tsonic"].to_numpy(float),
               df["scalar"].to_numpy(float), 20, -10.0, 10.0)
        yield ("dyco_" + label + "_window", df["w"].to_numpy(float),
               df["tsonic"].to_numpy(float), df["scalar"].to_numpy(float), 20, 0.0, 9.0)


def real_cases(nfiles):
    files = sorted(GAPPY.glob("CH-LAE_ec_preproc_10hz_*.csv"))
    if not files:
        return
    pick = [files[i] for i in np.linspace(0, len(files) - 1, min(nfiles, len(files))).astype(int)]
    for path in pick:
        df = pd.read_csv(path)
        rot = WindDoubleRotation(df["U"], df["V"], df["W"])
        w = np.asarray(rot.w2, dtype=float)
        t = df["T_SONIC"].to_numpy(float)
        for gas, lws, uws in REAL_GASES:
            s = df[gas].to_numpy(float)
            if np.isfinite(s).mean() < 0.3:
                continue
            yield ("%s:%s" % (path.stem[-13:], gas), w, t, s, 10, lws, uws)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--files", type=int, default=4, help="real CH-LAE periods to use")
    ap.add_argument("--keep", help="keep the driver inputs in this directory")
    ap.add_argument("--seed", type=int, default=20261002)
    ap.add_argument("--table", action="append", default=[],
                    help="an engine *_pwb_timelag_*.csv to run stage 2 on (repeatable)")
    ap.add_argument("--no-stage1", action="store_true", help="only stage 2")
    args = ap.parse_args()
    if not DRIVER.exists():
        print("pwb_compare driver not built: run `mingw32-make pwbcmp` in prj/")
        return 2

    cases = [] if args.no_stage1 else list(synthetic_cases()) + list(real_cases(args.files))
    counts = {"FAIL": 0}
    print("%-34s %-26s %-26s %s" % ("case", "engine (combo lag hdi)", "dyco (combo lag hdi)", "notes"))
    for i, (name, w, t, s, hz, lws, uws) in enumerate(cases):
        results, summ = run_case(name, w, t, s, hz, lws, uws, args.seed + i, args.keep)
        notes = []
        for kind, msg in results:
            counts[kind] = counts.get(kind, 0) + 1
            notes.append("%s: %s" % (kind, msg))
        if summ:
            eng = "%s %5d %6.2f s" % (summ[0], summ[1], summ[2])
            dy = "%s %5d %6.2f s" % (summ[3], summ[4], summ[5])
        else:
            eng = dy = "-"
        print("%-34s %-26s %-26s %s" % (name, eng, dy, "; ".join(notes) if notes else "identical"))
    print()
    print("cases: %d   %s" % (len(cases), "   ".join("%s: %d" % kv for kv in sorted(counts.items()))))
    for table in args.table:
        print()
        counts["FAIL"] += ladder_stage(table)
    return 1 if counts["FAIL"] else 0


if __name__ == "__main__":
    sys.exit(main())
