"""Every CSV a run produces must have as many values as it has column names.

Nothing else catches this. `check_h2o_late.py` compares by header name, so a
uniform shift leaves every name intact and every comparison passing while the
values underneath belong to the column to their left. `diff -r` against a
reference catches it only if the reference was correct, and the reference is
regenerated from the same writer.

The fault this was written for: WriteOutFluxnetFcc passed `separator` where
AddCharDatumToDataline expects `err_label`, so a cell WriteDatumChar judged
missing was replaced by a comma - a field boundary rather than a value. The row
gained a column and everything after it shifted by one. It stayed hidden
because WriteDatumChar's missing-value test is the literal '899999999', which
the time-lag cell only produces at exactly eight gases: one filler digit plus
one per gas, all nines when the test was not performed. Four- and five-gas
projects made a shorter cell and never tripped it.

Usage:  check_columns.py out_<dir> [out_<dir> ...]
"""
import csv
import glob
import os
import sys


#> Paths whose rows are deliberately ragged. The binned and ensemble
#> (co)spectra and the spectral assessment carry a preamble of their own
#> before the table starts, so a plain field-count comparison does not
#> describe them. Matched against the whole path, because these live in
#> subdirectories named for the product rather than the file.
SKIP = ('binned_', 'cospectra', 'spectral_analysis', 'processing_')


def rows_of(path):
    with open(path, newline='') as fh:
        return list(csv.reader(fh))


def check(path):
    """(header_count, first bad row number, its count) or None if clean.

    The header is the first row carrying more than one field, not simply the
    first row. The statistics files open with a single-field title -
    `first_statistics:_on_raw_data` - and the names are on the line below, so
    taking row one flagged all four of them as one-column headers against
    hundred-and-twenty-value rows. That is the check crying wolf about the
    only files it was never able to describe, and a checker with four standing
    false positives is one nobody reads.
    """
    rows = rows_of(path)
    start = next((i for i, r in enumerate(rows) if len(r) > 1), None)
    if start is None or len(rows) < start + 2:
        return None
    width = len(rows[start])
    for n, row in enumerate(rows[start + 1:], start=start + 2):
        if not any(field.strip() for field in row):
            continue          #< trailing blank line
        if len(row) != width:
            return width, n, len(row)
    return None


def main():
    targets = sys.argv[1:]
    if not targets:
        print(__doc__.strip().splitlines()[-1])
        return 2

    failed = checked = 0
    for target in targets:
        for path in sorted(glob.glob(os.path.join(target, '**', '*.csv'),
                                     recursive=True)):
            if any(tag in path.replace(os.sep, '/') for tag in SKIP):
                continue
            checked += 1
            bad = check(path)
            if bad is None:
                continue
            width, n, got = bad
            print('%s: header names %d columns, row %d has %d values'
                  % (os.path.relpath(path), width, n, got))
            failed += 1

    if failed:
        print('\n%d of %d files are misaligned - a value is sitting under the '
              'wrong name.' % (failed, checked))
        return 1
    print('%d files checked, every row matches its header' % checked)
    return 0


if __name__ == '__main__':
    sys.exit(main())
