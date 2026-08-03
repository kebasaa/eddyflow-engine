#!/usr/bin/env bash
# Rebuilt regression harness for the eddyflow engine.
#
# Usage: run.sh ref|chk
#
# Runs RP then FCC over a 3-hour subset of the CH-LAE dataset and normalises
# the run timestamps out of the output, so two runs can be diffed.
#
# Three traps this guards against, all of which produced false negatives
# before:
#   1. RP only. FCC recomputes the fluxes under fcc_follows, so an RP-only run
#      compares files nothing wrote.
#   2. A stale output dir. `find | head -1` happily picks up a previous run's
#      CSV, so the dir is emptied and the match count asserted.
#   3. `find` returning nothing makes diff compare two empty streams and print
#      IDENTICAL. Files are asserted to exist before anything is compared.
set -euo pipefail

WHICH="${1:?usage: run.sh ref|chk}"
HERE="$(cd "$(dirname "$0")" && pwd)"
# Defaults to the release build the build script produces. Override with BIN=
# to compare a working tree against itself: both ref and chk must come from the
# same binaries, or the diff reports build differences as regressions.
BIN="${BIN:-/c/Users/jonmuell/Documents/GitHub/build/eddyflow-engine-win-release/bin}"
# The engine links the gfortran runtime dynamically and the build does not
# copy it next to the binaries.
export PATH="/c/Users/jonmuell/mingw64/bin:$PATH"
OUT="$HERE/out_$WHICH"
HOME_DIR="$HERE/home"

rm -rf "$OUT"; mkdir -p "$OUT"
rm -rf "$HOME_DIR/tmp"; mkdir -p "$HOME_DIR/tmp" "$HOME_DIR/ini"

# Point the project at this run's output dir. The ex_file is what FCC reads
# back, so it has to name the file RP writes in *this* run.
WIN_OUT="$(cygpath -w "$OUT" | sed 's|\\|/|g')"
sed -e "s|^out_path=.*|out_path=$WIN_OUT|" "$HERE/${BASE:-base.eddyflow}" > "$HERE/run_$WHICH.eddyflow"

echo "== RP =="
"$BIN/eddyflow_rp.exe" "$(cygpath -w "$HERE/run_$WHICH.eddyflow")" -e "$(cygpath -w "$HOME_DIR")/" > "$OUT/_rp.log" 2>&1 \
    || { echo "RP FAILED"; tail -30 "$OUT/_rp.log"; exit 1; }

# FCC reads the ex (FLUXNET) file RP just wrote.
EXFILE="$(find "$OUT" -name "*fluxnet*.csv" | head -1)"
[ -n "$EXFILE" ] || { echo "no fluxnet file from RP"; tail -30 "$OUT/_rp.log"; exit 1; }
sed -i "s|^ex_file=.*|ex_file=$(cygpath -w "$EXFILE" | sed 's|\\|/|g')|" "$HERE/run_$WHICH.eddyflow"

# A fixture whose sa_bin_spectra is the token SELF reads the binned files this
# run just wrote, rather than the shared directory the other fixtures point
# at. That shared directory predates the N-gas binned format and has four
# gases, which is what makes it the backward-compatibility case; SELF is the
# forward one, and there is no other way to get it without a second run.
if grep -q '^sa_bin_spectra=SELF' "$HERE/run_$WHICH.eddyflow"; then
    sed -i "s|^sa_bin_spectra=SELF|sa_bin_spectra=$WIN_OUT/eddyflow_binned_cospectra|"         "$HERE/run_$WHICH.eddyflow"
fi

echo "== FCC =="
"$BIN/eddyflow_fcc.exe" "$(cygpath -w "$HERE/run_$WHICH.eddyflow")" -e "$(cygpath -w "$HOME_DIR")/" > "$OUT/_fcc.log" 2>&1 \
    || { echo "FCC FAILED"; tail -30 "$OUT/_fcc.log"; exit 1; }

# Normalise: strip the run timestamp from names and from file contents.
#
# This MUST recurse. Most of the output lives in the per-period subdirectories
# (binned_cospectra, binned_ogives, full_cospectra, spectral_analysis) and every
# one of those files carries the run timestamp in its name; normalising only the
# top level leaves 21 of 25 files reading as added/removed under `diff -r`, which
# buries whatever really changed. The processing_*.eddyflow copy of the project
# needs it too - it is not a .csv, so the old glob skipped it.
cd "$OUT"
rm -f _rp.log _fcc.log

# Rename depth-first so a renamed parent cannot invalidate a child's path.
find . -depth \( -name '*.csv' -o -name '*.txt' -o -name '*.eddyflow' \) -print |
while IFS= read -r f; do
    n="$(echo "$f" | sed -E 's/_?[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{6}//g')"
    [ "$f" = "$n" ] || mv "$f" "$n"
done

find . -type f \( -name '*.csv' -o -name '*.txt' -o -name '*.eddyflow' \) -print |
while IFS= read -r f; do
    sed -i -E 's/[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{6}/TIMESTAMP/g' "$f"
    sed -i -E 's#(out_path|ex_file|file_name|proj_file)=.*#\1=PATH#' "$f"
done

echo "== $WHICH: $(find . -type f | wc -l) output files =="
find . -type f | sort
