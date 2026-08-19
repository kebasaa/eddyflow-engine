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

# Keep RP's own log. Both binaries write one, and normalisation strips the
# timestamp that distinguishes them - so FCC's overwrote RP's and every
# RP-side message, errors included, was absent from the compared artefacts.
# The README's "the run log is compared too" only ever held for FCC.
for rplog in "$OUT"/*_log_*.log; do
    [ -e "$rplog" ] || continue
    mv "$rplog" "${rplog%.log}_rp.log"
done

# FCC reads the ex (FLUXNET) file RP just wrote.
EXFILE="$(find "$OUT" -name "*fluxnet*.csv" | head -1)"
[ -n "$EXFILE" ] || { echo "no fluxnet file from RP"; tail -30 "$OUT/_rp.log"; exit 1; }
sed -i "s|^ex_file=.*|ex_file=$(cygpath -w "$EXFILE" | sed 's|\\|/|g')|" "$HERE/run_$WHICH.eddyflow"

# SELF means "the (co)spectra this run just wrote", which is the only way to
# get them without a second run. Every fixture uses it now. They used to name a
# shared directory outside the repo instead - which did not exist on any
# machine, so the engine hit Error(87), demoted the in-situ method to Moncrieff
# and finished, and the sweep passed it because a degraded run still matches its
# own headers. 29 of the fixtures configure Fratini and none of them was
# running it. A missing path is fatal now, so that cannot recur silently.
if grep -q '^sa_bin_spectra=SELF' "$HERE/run_$WHICH.eddyflow"; then
    sed -i "s|^sa_bin_spectra=SELF|sa_bin_spectra=$WIN_OUT/eddyflow_binned_cospectra|"         "$HERE/run_$WHICH.eddyflow"
fi
# The full cospectra directory needs the same treatment, and did not get it.
# Fratini reads it, 29 of the fixtures ask for Fratini, and a path that is not
# there is now fatal rather than a silent demotion to Moncrieff - so this line
# is the difference between the suite running the method it configures and not
# running at all.
if grep -q '^sa_full_spectra=SELF' "$HERE/run_$WHICH.eddyflow"; then
    sed -i "s|^sa_full_spectra=SELF|sa_full_spectra=$WIN_OUT/eddyflow_full_cospectra|"         "$HERE/run_$WHICH.eddyflow"
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
find . -depth \( -name '*.csv' -o -name '*.txt' -o -name '*.eddyflow' -o -name '*.log' \) -print |
while IFS= read -r f; do
    n="$(echo "$f" | sed -E 's/_?[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{6}//g')"
    [ "$f" = "$n" ] || mv "$f" "$n"
done

find . -type f \( -name '*.csv' -o -name '*.txt' -o -name '*.eddyflow' -o -name '*.log' \) -print |
while IFS= read -r f; do
    sed -i -E 's/[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{6}/TIMESTAMP/g' "$f"
    sed -i -E 's#(out_path|ex_file|file_name|proj_file)=.*#\1=PATH#' "$f"
    # The run log is the engine's console, captured, so it is compared like
    # everything else. Three things in it differ between two runs of an
    # unchanged tree and say nothing about the results: the wall clock, how
    # long each period took, and which of the two output directories this run
    # was pointed at. The last one also settles a long-standing false positive
    # in base_n_gas_bin, whose SELF token records the same path.
    sed -i -E 's/[0-9]{4}-[0-9]{2}-[0-9]{2} [0-9]{2}:[0-9]{2}:[0-9]{2}\.[0-9]{3}/CLOCK/g' "$f"
    sed -i -E 's/[0-9]+:[0-9]{2}:[0-9]{2}\.[0-9]{3}/ELAPSED/g' "$f"
    sed -i -E 's#out_(ref|chk)#OUT#g' "$f"
    sed -i -E 's#run_(ref|chk)\.eddyflow#run_WHICH.eddyflow#g' "$f"
done

echo "== $WHICH: $(find . -type f | wc -l) output files =="
find . -type f | sort
