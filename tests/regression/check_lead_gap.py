"""The FLUXNET file must have a row for every period of the window, in order.

Gate for base_lead_gap, whose window opens two hours before the first raw
file. The FLUXNET header is written lazily, after the first period that
imports data, because its custom-variable columns - and so the width of every
row - are only known then. The skipped periods ahead of that were written
anyway, to the unopened unit: gfortran connected it to fort.132 in the working
directory, the rows went there, and the FLUXNET file started at the first
period with data. check_columns.py could not see it, since every row that did
reach the file matched its header.

So this checks what that fault removed: the first row starts at the window's
start, every half-hour follows without a gap or a repeat up to the window's
end, and the periods before the first raw file are marked not_enough_data.

Usage:  check_lead_gap.py out_<dir> [start end first_file]
        timestamps as yyyymmddHHMM; the defaults are base_lead_gap's
"""
import csv
import datetime
import glob
import os
import sys

FMT = '%Y%m%d%H%M'
STEP = datetime.timedelta(minutes=30)


def fluxnet_file(out_dir):
    """RP's FLUXNET file, which run.sh keeps as *_rp.csv.

    Not the plain *fluxnet*.csv: that is FCC's, which starts at the first
    valid record and so would hide exactly the rows this is about.
    """
    found = glob.glob(os.path.join(out_dir, '*fluxnet*_rp.csv'))
    if len(found) != 1:
        raise SystemExit('expected one RP FLUXNET file in %s, found %s'
                         % (out_dir, found))
    return found[0]


def main():
    args = sys.argv[1:]
    if not args:
        raise SystemExit(__doc__)
    out_dir = args[0]
    start, end, first_file = (args[1:4] if len(args) >= 4
                              else ('202505312200', '202506010300',
                                    '202506010000'))
    start = datetime.datetime.strptime(start, FMT)
    end = datetime.datetime.strptime(end, FMT)
    first_file = datetime.datetime.strptime(first_file, FMT)

    path = fluxnet_file(out_dir)
    with open(path, newline='', encoding='utf-8') as fh:
        rows = list(csv.reader(fh))
    header, data = rows[0], rows[1:]
    i_start = header.index('TIMESTAMP_START')
    i_end = header.index('TIMESTAMP_END')
    i_file = header.index('FILENAME_HF')

    problems = []
    expected = start
    for n, row in enumerate(data, start=2):
        got = datetime.datetime.strptime(row[i_start], FMT)
        if got != expected:
            problems.append('line %d starts %s, expected %s'
                            % (n, got.strftime(FMT), expected.strftime(FMT)))
            break
        if datetime.datetime.strptime(row[i_end], FMT) != got + STEP:
            problems.append('line %d ends %s' % (n, row[i_end]))
        if got + STEP <= first_file and row[i_file] != 'not_enough_data':
            problems.append('line %d is before the first raw file but names %s'
                            % (n, row[i_file]))
        expected = got + STEP
    if not problems and expected != end:
        problems.append('rows stop at %s, the window ends at %s'
                        % (expected.strftime(FMT), end.strftime(FMT)))

    name = os.path.basename(path)
    if problems:
        for p in problems:
            print('%s: %s' % (name, p))
        print('FAIL: the FLUXNET file does not cover the window')
        sys.exit(1)
    print('%s: %d rows, %s to %s, the leading %d not_enough_data'
          % (name, len(data), start.strftime(FMT), end.strftime(FMT),
             (first_file - start) // STEP))
    print('every period of the window has its row, in order')


if __name__ == '__main__':
    main()
