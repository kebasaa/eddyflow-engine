"""Gate for base_h2o_late: displacing water from slot 6 must change nothing.

`base_h2o_late` is `base_n_gas` with records 2 and 5 swapped, so slot 6 holds
N2O and slot 9 holds H2O. The same physical columns feed the same physical
quantities, so every quantity derived from *the site's water* must come out
identical. Anything that moves is a place where the engine read the slot
instead of resolving the species.

Two column sets, deliberately separated:

SCALARS - one per site, derived from the water column and the anemometer
alone. These must match exactly. They are the gate.

PER_GAS_PREFIXES - each gas's own flux and mixing ratio. These must also
match once dilution resolves through `moist_ref`, but they additionally
depend on per-record settings that the swap legitimately moves (the PWB
timelag cache is keyed by species label, so H2O's entry is looked up
differently at slot 9). Reported, not asserted, until Stage B lands.

Usage:  check_h2o_late.py out_<ngas_run> out_<late_run>
"""
import csv
import glob
import os
import sys

#> Quantities that are one-per-site and water-derived. If the engine resolves
#> its water rather than assuming slot 6, every one of these is unchanged.
SCALARS = (
    'TAU', 'H', 'LE', 'ET', 'FH2O',
    'MO_LENGTH', 'ZL', 'BOWEN', 'TSTAR',
    'TA_EP', 'RH_EP', 'VPD_EP', 'TDEW',
    'AIR_DENSITY', 'AIR_RHO_CP', 'AIR_CP',
    'VAPOR_DENSITY', 'VAPOR_PARTIAL_PRESSURE', 'VAPOR_PARTIAL_PRESSURE_SAT',
    'SPECIFIC_HUMIDITY', 'SPECIFIC_HEAT_EVAP', 'VAPOR_DRYAIR_RATIO',
    'DRYAIR_PARTIAL_PRESSURE', 'DRYAIR_DENSITY', 'DRYAIR_MV',
)

PER_GAS_PREFIXES = ('FC', 'FCOS', 'FN2O', 'FCO2_2', 'FH2O_2', 'FN2O_2')
MIXING = '_MIXING_RATIO'


def load(directory):
    matches = glob.glob(os.path.join(directory, '*fluxnet_adv.csv'))
    if not matches:
        raise SystemExit('no fluxnet file in %s' % directory)
    rows = list(csv.reader(open(matches[0])))
    return rows[0], rows[1:]


def main():
    if len(sys.argv) != 3:
        raise SystemExit(__doc__)
    ref_dir, late_dir = sys.argv[1], sys.argv[2]
    href, dref = load(ref_dir)
    hlate, dlate = load(late_dir)

    if len(dref) != len(dlate):
        raise SystemExit('row counts differ: %d vs %d' % (len(dref), len(dlate)))

    #> By name, never by position: the record order changed, so the per-gas
    #> column groups sit at different offsets in the two files.
    failures = []
    for name in SCALARS:
        if name not in href or name not in hlate:
            failures.append('%s: absent from one of the two files' % name)
            continue
        i, j = href.index(name), hlate.index(name)
        for n, (ra, rb) in enumerate(zip(dref, dlate)):
            if ra[i] != rb[j]:
                failures.append('%s row %d: %s vs %s' % (name, n + 1, ra[i], rb[j]))
                break

    print('scalars checked: %d' % len(SCALARS))
    if failures:
        print('FAIL - water is still being read from its slot:')
        for f in failures:
            print('  ' + f)
    else:
        print('PASS - every water-derived scalar is unchanged')

    #> Informational: the per-gas family.
    print()
    print('per-gas (informational):')
    for name in PER_GAS_PREFIXES + tuple(
            c for c in href if c.endswith(MIXING)):
        if name not in href or name not in hlate:
            continue
        i, j = href.index(name), hlate.index(name)
        a, b = dref[0][i], dlate[0][j]
        print('  %-24s %-16s %-16s %s' % (name, a, b, '' if a == b else '<-- moved'))

    sys.exit(1 if failures else 0)


if __name__ == '__main__':
    main()
