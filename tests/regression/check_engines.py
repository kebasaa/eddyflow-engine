"""Compare two full_output tables column by NAME, not by position.

Written for the extended-.ghg compatibility test, which has to answer two
different questions with the same machinery:

  * EddyPro on a classic archive against EddyPro on the extended one, where
    the answer must be "no column moved at all" - the whole claim of the format
    is that the extension is free.
  * EddyFlow against EddyPro on the same archive, where some columns must agree
    exactly, some may differ by a known amount, and one is expected to differ
    because we fixed something EddyPro gets wrong.

By name because the two programs do not write the same columns: EddyPro's
express output has ~191 and EddyFlow's ~164, and the shared ones are not in the
same order. Comparing by position would line `H` up against whatever EddyFlow
happens to write in EddyPro's eleventh slot and report a catastrophe.

Usage:
    check_engines.py A.csv B.csv --mode identical
    check_engines.py A.csv B.csv --mode engines

Exit status is 0 when every gated column passed.
"""
import argparse
import io
import sys


#> Header is row 2 and units row 3, so data starts at row 4. Both programs
#> write that shape; it is the LI-COR ancestry they share.
HEADER_ROW = 1
FIRST_DATA_ROW = 3

MISSING = {'-9999', '-9999.0', 'NaN', 'nan', ''}

#> Every gated column, its allowed relative difference, and the reason for it.
#> Prefix match, longest wins, so `co2_molar_density` and its siblings come
#> along without being listed one by one.
#>
#> Nothing here is exact equality, and that is deliberate. Both programs write
#> about six significant digits, so requiring two independently-built binaries
#> to print the same text is a claim about formatting rather than about the
#> physics: measured, the wind statistics agree to about 5e-6, which is the
#> last printed digit moving. Exact equality is still the right gate for
#> EddyPro against EddyPro (`--mode identical`), where the same binary produces
#> both sides.
PRECISION = 2e-5     #: last printed digit
CHAIN = 5e-4         #: the humidity chain, below
CANCELLATION = 5e-3  #: quantities formed by subtracting two similar numbers

GATES = (
    #> Quantities neither program has changed. Anything here moving by more
    #> than the printed digit means something in the shared physics moved.
    (PRECISION, 'wind and raw concentrations', (
        'Tau', 'u*', 'TKE', 'L', 'wind_speed', 'max_wind_speed', 'wind_dir',
        'u_unrot', 'v_unrot', 'w_unrot', 'u_rot', 'v_rot', 'w_rot',
        'co2_molar_density', 'co2_mixing_ratio', 'co2_mole_fraction',
        'h2o_molar_density', 'h2o_mixing_ratio', 'h2o_mole_fraction',
        'co2_time_lag', 'h2o_time_lag',
        'qc_Tau', 'qc_H', 'qc_LE', 'qc_co2_flux', 'qc_h2o_flux',
        #> Second-moment and footprint quantities. Measured on the first run of
        #> this test at 1.2e-6 to 8.6e-6 - the same last-digit level as the
        #> wind statistics they are built from. Listed rather than left
        #> unexplained so that a real change in them fails instead of being
        #> read as noise.
        'u_var', 'v_var', 'w_var', 'ts_var', 'co2_var', 'h2o_var', 'ch4_var',
        'x_peak', 'x_offset', 'x_10%', 'x_30%', 'x_50%', 'x_70%', 'x_90%',
        '(z-d)/L',
    )),
    #> EddyFlow refined the molecular weights used through the humidity chain.
    #> The changelog predicts ~0.026 %, and the 48-archive comparison recorded
    #> in sweep.sh measured H within 2e-4 and LE within 2.7e-4.
    (CHAIN, 'humidity chain, refined molecular weights', (
        'H', 'LE', 'ET', 'h2o_flux', 'co2_flux', 'bowen_ratio',
        'water_vapor_density', 'e', 'es', 'specific_humidity', 'RH', 'Tdew',
        'sonic_temperature', 'air_temperature', 'air_pressure', 'air_density',
        'air_molar_volume', 'air_heat_capacity',
        #> The random uncertainties and the vertical advection terms are
        #> computed FROM the fluxes above, so they carry the same difference.
        #> Measured on the first run at 2.6e-4 for un_LE, the rest below 1e-5.
        'un_', 'co2_v-adv', 'h2o_v-adv', 'ch4_v-adv', 'none_v-adv',
    )),
    #> VPD is es - e, two numbers within a couple of percent of each other, so
    #> the 0.026 % on each component surfaces here roughly a hundredfold. The
    #> amplification is arithmetic, not a second disagreement.
    (CANCELLATION, 'difference of two similar numbers', ('VPD',)),
)

#> Reported and never gated, each with the reason it differs. A column here is
#> one somebody has looked at and understood; anything NOT here and not gated
#> comes out as "unlisted", which is the report's way of saying nobody has.
UNGATED = (
    #> EddyPro leaves the LI-7500's signal-strength slot NaN and mislabels the
    #> analyser as an LI-7200; EddyFlow fills it. Everything downstream differs
    #> BECAUSE we fixed it, so gating it would be gating the bug.
    ('EddyPro leaves the AGC slot NaN',
     ('ch4_', 'qc_ch4_flux', 'rand_err_ch4', 'none_', 'qc_none_')),
    #> The packed per-variable flag strings: a leading 8 then one digit per
    #> variable. EddyPro always writes its four fixed gas slots, EddyFlow only
    #> the gases the project actually has - three here - so the strings are
    #> legitimately a digit shorter. Comparing them as numbers is meaningless
    #> either way; they are digit strings, not quantities.
    ('flag string, one digit per gas - EddyFlow has 3, EddyPro always 4', (
        'spikes_hf', 'spikes_sf', 'amplitude_resolution_hf',
        'amplitude_resolution_sf', 'drop_out_hf', 'drop_out_sf',
        'absolute_limits_hf', 'absolute_limits_sf', 'skewness_kurtosis_hf',
        'skewness_kurtosis_sf', 'discontinuities_hf', 'discontinuities_sf',
        'timelag_hf', 'timelag_sf', 'attack_angle_hf', 'non_steady_wind_hf',
    )),
    #> Whether the nominal lag was used rather than a detected one. The two
    #> programs choose differently; a boolean's relative difference from zero
    #> is infinite and says nothing.
    ('boolean - the two programs choose the lag differently',
     ('co2_def_timelag', 'h2o_def_timelag', 'ch4_def_timelag',
      'none_def_timelag')),
    #> MEASURED: identical to six digits and opposite in sign - EddyPro
    #> +0.0966893 against EddyFlow -0.0966890 on the same period. That is a
    #> sign convention, not a disagreement about the physics, and it is not
    #> about the archive format either way. Recorded here so the report does
    #> not present it as an unexplained 200 % error every run.
    ('sign convention differs; magnitudes agree to 6 digits', ('T*',)),
)


def read_table(path):
    """-> (header, [row, ...]), rows as lists of strings."""
    with io.open(path, encoding='cp1252', errors='replace') as fh:
        rows = [ln.rstrip('\n').rstrip('\r').split(',') for ln in fh]
    if len(rows) <= FIRST_DATA_ROW:
        sys.exit('%s has no data rows' % path)
    return rows[HEADER_ROW], rows[FIRST_DATA_ROW:]


def classify(name):
    """-> (tolerance, reason). None tolerance means ungated.

    Longest prefix wins, which matters: 'ch4_flux' is ungated while 'co2_flux'
    is on the humidity chain, and a plain startswith over an unordered list
    would let 'H' claim 'H_strg' or 'e' claim 'es'. Length settles it the same
    way every time rather than by however the tuples happen to be written.
    """
    best, hit = '', (None, None)
    for reason, prefixes in UNGATED:
        for p in prefixes:
            if (name == p or name.startswith(p)) and len(p) > len(best):
                best, hit = p, (None, reason)
    for tol, reason, prefixes in GATES:
        for p in prefixes:
            if (name == p or name.startswith(p)) and len(p) > len(best):
                best, hit = p, (tol, reason)
    return hit


def worst_delta(a_rows, b_rows, ia, ib):
    """Largest relative difference over the paired rows, or None if none differ.

    A pair where either side is missing is skipped rather than counted as a
    difference: the two programs do not always decline to report the same
    period, and a -9999 against a number says nothing about agreement.
    """
    worst = None
    for ra, rb in zip(a_rows, b_rows):
        if ia >= len(ra) or ib >= len(rb):
            continue
        x, y = ra[ia].strip(), rb[ib].strip()
        if x == y:
            continue
        if x in MISSING or y in MISSING:
            continue
        try:
            fx, fy = float(x), float(y)
        except ValueError:
            #> Text that differs - a filename, a flag string. Reported as an
            #> infinite delta so it cannot pass a numeric gate silently.
            return float('inf')
        if fx == fy:
            continue
        rel = abs(fy - fx) / abs(fx) if fx else float('inf')
        if worst is None or rel > worst:
            worst = rel
    return worst


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('a')
    ap.add_argument('b')
    ap.add_argument('--mode', choices=('identical', 'engines'), required=True)
    ap.add_argument('--label-a', default='A')
    ap.add_argument('--label-b', default='B')
    args = ap.parse_args()

    ha, ra = read_table(args.a)
    hb, rb = read_table(args.b)

    #> Columns that carry the run rather than the result. `filename` names the
    #> archive and legitimately differs between two output directories.
    skip = {'filename', 'date', 'time', 'DOY', 'daytime', 'file_records',
            'used_records'}
    shared = [n for n in ha if n in hb and n not in skip]

    #> A name-keyed comparison that matches nothing also passes. Say how many
    #> columns were actually compared, and refuse to call it a pass on none.
    print('  %s vs %s: %d shared columns, %d rows'
          % (args.label_a, args.label_b, len(shared), min(len(ra), len(rb))))
    if not shared:
        print('  FAIL: no shared columns - the two headers have nothing in '
              'common, so nothing was compared')
        return 1
    if len(ra) != len(rb):
        print('  FAIL: %d rows against %d - the two runs did not cover the '
              'same periods' % (len(ra), len(rb)))
        return 1

    ia = {n: i for i, n in enumerate(ha)}
    ib = {n: i for i, n in enumerate(hb)}

    failures, moved, reported = [], 0, []
    for name in shared:
        d = worst_delta(ra, rb, ia[name], ib[name])
        if d is None:
            continue
        moved += 1
        if args.mode == 'identical':
            #> Same binary on both sides, so nothing may move at all.
            failures.append((name, d, 'must be identical'))
            continue
        tol, reason = classify(name)
        if tol is None:
            reported.append((name, d, reason or 'unlisted'))
        elif d > tol:
            failures.append((name, d, 'over %.0e - %s' % (tol, reason)))
        else:
            reported.append((name, d, 'within %.0e' % tol))

    print('  %d of %d shared columns differ' % (moved, len(shared)))
    #> Split, because the two mean different things. An explained column is one
    #> somebody has looked at; an unlisted one is a difference nobody has
    #> accounted for yet, and burying those under a dozen known ones is how a
    #> new disagreement goes unnoticed.
    unlisted = [r for r in reported if r[2] == 'unlisted']
    explained = [r for r in reported if r[2] != 'unlisted']
    if unlisted:
        print('  NOT YET EXPLAINED - %d column(s):' % len(unlisted))
        for name, d, _ in sorted(unlisted, key=lambda r: -r[1]):
            print('    %-28s %10.3e' % (name, d))
    for name, d, why in sorted(explained, key=lambda r: -r[1])[:8]:
        print('    %-28s %10.3e  %s' % (name, d, why))
    if failures:
        print('  FAIL:')
        for name, d, why in sorted(failures, key=lambda r: -r[1]):
            print('    %-28s %10.3e  %s' % (name, d, why))
        return 1
    print('  ok')
    return 0


if __name__ == '__main__':
    sys.exit(main())
