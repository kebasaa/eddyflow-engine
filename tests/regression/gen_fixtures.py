"""Derive fixtures that displace water from its historical slot.

Gas slots are assigned by record order - `slot = firstGas + i - 1` - so the
constants co2=5, h2o=6, ch4=7 are aliases for "record 1..4" and nothing more.
Every fixture so far happens to put water in record 2, which means every
"slot 6 is water" assumption in the engine has been untestable: the assumption
was always true, so no output could ever contradict it.

`base_h2o_late` is a straight swap of records 2 and 5 of `base_n_gas` -
nothing else changes, not even the molecular weights. Slot 6 then holds N2O
and slot 9 holds H2O. Both are on the MIRO, so every MIRO gas resolves its
moisture to slot 9 and every LI-7200 gas to slot 11, and the same physical
columns feed the same physical quantities. That makes the gate exact rather
than approximate:

    every LE / ET / RH / VPD / air-density / specific-humidity / H2O-flux cell
    must equal base_n_gas's, compared BY HEADER NAME - the record order
    changed, so the per-gas column groups move and a positional diff is
    meaningless.

Run `gen_fixtures.py --check` to confirm the emitted files still match what
this script would write, the same way the tag-table generators are checked.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))

#> The record fields, in the order the interface writes them.
FIELDS = ('var', 'instr', 'col', 'moist', 'cell', 'mw', 'diff')


def read_records(lines, count):
    """Pull gas_<i>_<field> into a list of dicts, indexed from 0."""
    recs = [dict() for _ in range(count)]
    for ln in lines:
        if not ln.startswith('gas_'):
            continue
        key, _, value = ln.partition('=')
        parts = key.split('_')
        if len(parts) < 3 or not parts[1].isdigit():
            continue
        idx = int(parts[1]) - 1
        field = '_'.join(parts[2:])
        if 0 <= idx < count and field in FIELDS:
            recs[idx][field] = value
    return recs


def write_records(lines, recs):
    """Substitute the record block back, leaving every other line untouched."""
    out = []
    for ln in lines:
        if not ln.startswith('gas_'):
            out.append(ln)
            continue
        key, _, _ = ln.partition('=')
        parts = key.split('_')
        if len(parts) < 3 or not parts[1].isdigit():
            out.append(ln)
            continue
        idx = int(parts[1]) - 1
        field = '_'.join(parts[2:])
        if 0 <= idx < len(recs) and field in FIELDS:
            out.append('%s=%s' % (key, recs[idx][field]))
        else:
            out.append(ln)
    return out


def build_h2o_late(src_lines):
    count = 0
    for ln in src_lines:
        if ln.startswith('gas_num='):
            count = int(ln.split('=', 1)[1])
    recs = read_records(src_lines, count)
    #> Records 2 and 5 (1-based) -> slot 6 gets N2O, slot 9 gets H2O.
    recs[1], recs[4] = recs[4], recs[1]
    return write_records(src_lines, recs)


def build_mw(src_lines, mw1='', mw6=''):
    """base_n_gas pointed at the local metadata, with two CO2 records given
    distinct molecular weights.

    Molecular weight is only consulted on the g_m3 / mg_m3 / ug_m3 arms of the
    unit conversion, so the metadata copy sets both CO2 columns to ug_m3 -
    physically meaningless for a mole fraction, but it is the only way to make
    MW observable at all, and the test is a ratio.

    Two records, not one, because the defect has two halves: the *fallback*
    was gated on `slot > gas4`, and the *lookup* mapped a column's species
    name to a fixed slot. Giving the second CO2 its own weight catches the
    second half - both CO2 columns used to convert with record one's.
    """
    count = 0
    for ln in src_lines:
        if ln.startswith('gas_num='):
            count = int(ln.split('=', 1)[1])
    recs = read_records(src_lines, count)
    recs[0]['mw'] = mw1
    recs[5]['mw'] = mw6
    out = write_records(src_lines, recs)
    local = os.path.join(HERE, 'base_mw.metadata').replace(os.sep, '/')
    return [('proj_file=' + local) if ln.startswith('proj_file=') else ln
            for ln in out]


def build_own_binned(src_lines):
    """base_n_gas reading the binned (co)spectra it wrote itself.

    Every other fixture points sa_bin_spectra at a shared directory that
    predates the N-gas binned format and carries four gases - which is what
    makes those runs the backward-compatibility case, and why they are right
    to leave gases 5+ unassessed. SELF is the forward case: run.sh rewrites
    the token to this run's own output between RP and FCC.
    """
    return [('sa_bin_spectra=SELF') if ln.startswith('sa_bin_spectra=') else ln
            for ln in src_lines]


def build_tlag_opt(src_lines):
    """base_n_gas running the time-lag optimiser on the fly.

    Every other fixture sets tlag_meth=5 (block-bootstrap) and to_mode=0, which
    together mean AdjustTimelagOptSettings is never called - so the whole
    "derive a search window from the instrument geometry" path had no coverage
    at all, in a project format where a gas states its own window per record.

    tlag_meth=4 selects the optimiser and to_mode=1 runs it on the fly. The gas
    records are left alone deliberately: base_n_gas states a window for records
    one to four and none for five to eight, so one run exercises both arms -
    the declared window is taken verbatim, and the undeclared one is derived
    rather than read as a request for [0, 0].

    The optimisation period is the processing window, not the fortnight
    inherited from base_n_gas. That fortnight had the prepass walk 672
    half-hours to inform a run covering six of them, which took minutes against
    the seconds every other fixture takes - enough to keep a fixture out of the
    loop anyone actually runs. This is 38 seconds; a day is 113 and buys
    nothing.

    What that costs, stated plainly: six periods leave the optimiser too few
    determinations to fit a window, so toPasGas stays at the error value and
    SetTimelags falls back on the declared one. The fixture therefore covers
    the derivation path and the empty-optimiser guard - every gas takes it -
    but not the optimiser's own statistics, which need the fortnight. It still
    separates the defect from the fix: before, a gas declaring no window got
    the garbage +-0.222387 range built from a median of nothing; now it gets
    its declared window back.
    """
    period = {
        'to_start_date': _value(src_lines, 'pr_start_date'),
        'to_start_time': _value(src_lines, 'pr_start_time'),
        'to_end_date': _value(src_lines, 'pr_end_date'),
        'to_end_time': _value(src_lines, 'pr_end_time'),
    }
    out = []
    for ln in src_lines:
        key = ln.split('=', 1)[0]
        if key == 'tlag_meth':
            out.append('tlag_meth=4')
        elif key == 'to_mode':
            out.append('to_mode=1')
        elif key in period:
            out.append('%s=%s' % (key, period[key]))
        else:
            out.append(ln)
    return out


def _value(lines, key):
    for ln in lines:
        if ln.startswith(key + '='):
            return ln.split('=', 1)[1]
    raise SystemExit('no %s in the source fixture' % key)


def build_cell_ref(src_lines):
    """base_n_gas_cell with one gas pointed at the other analyser's cell.

    gas_<i>_cell names the cell record whose analyser holds a gas. It was
    parsed and discarded for as long as the cell slots were a single global
    set, and every fixture leaves it at 0 - so the field has only ever been
    exercised on its "auto" arm, where the instrument name is matched instead.

    base_n_gas_cell carries two analysers: cell records 1 and 2 are the MIRO's,
    3 and 4 the LI-7200's. Record one's gas is on the MIRO, so pointing it at
    cell record 3 makes the explicit reference and the name match disagree,
    which is the only way to tell whether the field is read at all. Its cell
    temperature and pressure - and every quantity computed from them - then
    come from the LI-7200's block.
    """
    return [('gas_1_cell=3' if ln.startswith('gas_1_cell=') else ln)
            for ln in src_lines]


TARGETS = {
    'base_cell_ref.eddyflow': ('base_n_gas_cell.eddyflow', build_cell_ref),
    'base_tlag_opt.eddyflow': ('base_n_gas.eddyflow', build_tlag_opt),
    'base_h2o_late.eddyflow': ('base_n_gas.eddyflow', build_h2o_late),
    'base_n_gas_bin.eddyflow': ('base_n_gas.eddyflow', build_own_binned),
    'base_mw_ref.eddyflow': ('base_n_gas.eddyflow',
                             lambda ls: build_mw(ls, '', '')),
    'base_mw.eddyflow': ('base_n_gas.eddyflow',
                         lambda ls: build_mw(ls, '90.0000', '30.0000')),
}


def main():
    check = '--check' in sys.argv
    failed = False
    for name, (source, build) in TARGETS.items():
        with open(os.path.join(HERE, source)) as fh:
            src_lines = fh.read().splitlines()
        text = '\n'.join(build(src_lines)) + '\n'
        path = os.path.join(HERE, name)
        if check:
            existing = open(path).read() if os.path.isfile(path) else None
            if existing != text:
                print('%s is out of date' % name)
                failed = True
            else:
                print('%s up to date' % name)
        else:
            with open(path, 'w') as fh:
                fh.write(text)
            print('wrote %s from %s' % (name, source))
    if failed:
        sys.exit(1)


if __name__ == '__main__':
    main()
