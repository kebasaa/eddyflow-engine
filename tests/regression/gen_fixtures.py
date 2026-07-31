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


TARGETS = {
    'base_h2o_late.eddyflow': ('base_n_gas.eddyflow', build_h2o_late),
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
