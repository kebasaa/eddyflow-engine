"""Give every fixture its per-gas settings as records, not as flat tags.

The engine used to read a flat `sr_lim_co2`-style tag into each of the four
legacy slots and then let an `SNTagFound`-guarded record override it. The
fixtures were written for that: they carry the flat tags and almost no
`gas_<i>_<setting>` keys at all. Removing the flat reader without rewriting
them first would not fail loudly - it would silently drop every per-gas
threshold in every fixture and leave the outputs looking plausible.

So this runs first, on its own, and must change nothing: a record key holding
the same value the flat tag held overrides it with itself. Any fixture whose
output moves is this script getting a value wrong, which is exactly what it
is for - the migration is self-verifying.

Records one to four take the value of the flat tag for their slot, in slot
order co2, h2o, ch4, gas4. That mapping is positional, not by species: the
engine assigns `sr%lim_gas(co2)` to slot five whatever species record one
names, so a fixture whose first record is CH4 still takes `sr_lim_co2`. This
reproduces that.

Records five and up are left alone, because the flat layer never reached them
either - they fall to a whole-array default, and after the flat reader goes
that default is stated in Fortran instead of being whatever the loader left.
`sr_lim` is the exception: its flat read assigns the CO2 value to the *whole*
array before naming the four slots, so every later gas inherits it today and
has to be written out to keep doing so.

Run with --check to confirm the fixtures still match, the way the tag-table
generators are checked.
"""
import glob
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))

#> Slot order. Record i takes the flat tag of SPECIES[i - 1].
SPECIES = ('co2', 'h2o', 'ch4', 'gas4')

#> (record field, flat key template). The template is formatted with the
#> species slug; a group whose flat keys are absent from a fixture is skipped,
#> so a project that never set them does not gain them.
#>
#> Several flat groups are missing an h2o member - to_min_flux and the three
#> sa thresholds are asked of CO2, CH4 and the fourth gas only, because water
#> is judged by LE instead. Record two gets no key for those, which leaves it
#> on the same whole-array default it has today.
GROUPS = [
    ('sr_lim',          'sr_lim_{s}'),
    ('al_min',          'al_{s}_min'),
    ('al_max',          'al_{s}_max'),
    ('ds_hf',           'ds_hf_{s}'),
    ('ds_sf',           'ds_sf_{s}'),
    ('tl_def',          'tl_def_{s}'),
    ('to_min_flux',     'to_{s}_min_flux'),
    ('to_min_lag',      'to_{s}_min_lag'),
    ('to_max_lag',      'to_{s}_max_lag'),
    ('pwb_min_lag',     'pwb_{s}_min_lag'),
    ('pwb_max_lag',     'pwb_{s}_max_lag'),
    ('out_full_sp',     'out_full_sp_{s}'),
    ('out_full_cosp_w', 'out_full_cosp_w_{s}'),
    ('out_raw',         'out_raw_{s}'),
    ('sa_fmin',         'sa_fmin_{s}'),
    ('sa_fmax',         'sa_fmax_{s}'),
    ('sa_hfn_fmin',     'sa_hfn_{s}_fmin'),
    ('sa_min_st',       'sa_min_st_{s}'),
    ('sa_min_un',       'sa_min_un_{s}'),
    ('sa_max',          'sa_max_{s}'),
]
GROUPS += [('drift_dir_%d' % k, 'drift_dir_{s}_%d' % k) for k in range(7)]
GROUPS += [('drift_inv_%d' % k, 'drift_inv_{s}_%d' % k) for k in range(7)]

#> The one setting whose flat read fills the whole array before naming slots,
#> so gases past the fourth inherit the CO2 value and must be written out.
WHOLE_ARRAY_FROM_CO2 = 'sr_lim'


def parse(lines):
    """Flat values by key, and the index of the last line of each key."""
    values, where = {}, {}
    for n, ln in enumerate(lines):
        if '=' not in ln or ln.lstrip().startswith(('[', '#', ';')):
            continue
        key, _, value = ln.partition('=')
        values[key.strip()] = value
        where[key.strip()] = n
    return values, where


def gas_count(values):
    try:
        return int(values.get('gas_num', '0') or 0)
    except ValueError:
        return 0


def build(lines):
    values, where = parse(lines)
    count = gas_count(values)
    if count <= 0:
        return lines, 0

    #> Insertions keyed by the line they follow, so each new key lands in the
    #> same section as the flat key it mirrors. Section membership is what
    #> ParseIniFile matches on, so a key in the wrong one is simply not found.
    additions = {}
    added = 0
    for field, template in GROUPS:
        flat = {s: template.format(s=s) for s in SPECIES}
        present = [s for s in SPECIES if flat[s] in values]
        if not present:
            continue
        anchor = max(where[flat[s]] for s in present)

        emit = []
        for i in range(1, count + 1):
            key = 'gas_%d_%s' % (i, field)
            if key in values:
                continue          #< already a record; never overwrite
            if i <= len(SPECIES):
                slug = SPECIES[i - 1]
                if flat[slug] not in values:
                    continue      #< no flat member for this slot
                emit.append('%s=%s' % (key, values[flat[slug]]))
            elif field == WHOLE_ARRAY_FROM_CO2 and 'co2' in present:
                emit.append('%s=%s' % (key, values[flat['co2']]))
        if emit:
            additions.setdefault(anchor, []).extend(emit)
            added += len(emit)

    out = []
    for n, ln in enumerate(lines):
        out.append(ln)
        out.extend(additions.get(n, ()))
    return out, added


def main():
    check = '--check' in sys.argv
    failed = False
    for path in sorted(glob.glob(os.path.join(HERE, '*.eddyflow'))):
        name = os.path.basename(path)
        if re.match(r'^run_', name):
            continue              #< harness scratch, rewritten every run
        with open(path) as fh:
            lines = fh.read().splitlines()
        out, added = build(lines)
        text = '\n'.join(out) + '\n'
        if check:
            if added:
                print('%s is missing %d record keys' % (name, added))
                failed = True
        elif added:
            with open(path, 'w') as fh:
                fh.write(text)
            print('%s: added %d record keys' % (name, added))
    if failed:
        sys.exit(1)


if __name__ == '__main__':
    main()
