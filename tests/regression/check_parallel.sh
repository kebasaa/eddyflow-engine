#!/usr/bin/env bash
# Does splitting a pre-pass across worker processes change the answer?
#
# Usage: check_parallel.sh [fixture.eddyflow]   (default: base_tlag_par.eddyflow)
#
# Runs the fixture twice through run.sh - once with -j 1, once with -j 0 -
# and diffs the two normalised output trees. Every file must match.
#
# Why this is a separate script and not just a sweep fixture: sweep.sh does
# not diff against a reference at all, and run.sh passes no -j, so the stored
# reference for base_tlag_par is ITSELF a parallel run. There was no gate on
# serial-versus-parallel equivalence anywhere in the suite - the claim was
# checked by hand once and then only asserted. This makes it repeatable.
#
# The run log is excluded, and only the run log. A parallel run concatenates
# each worker's own log into the parent's, so it legitimately differs; every
# other artefact - fluxes, spectra, the project copy - must be byte-identical.
#
# The fixture has to be one whose pre-pass is long enough for the engine to
# bother splitting. PlanPrepassBatches refuses a range too short to pay for
# the processes, so on most fixtures this would compare two serial runs and
# pass without testing anything. That is what the "did it actually split"
# assertion below is for.
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
FIXTURE="${1:-base_tlag_par.eddyflow}"

echo "== serial (-j 1) =="
RP_EXTRA="-j 1" BASE="$FIXTURE" "$HERE/run.sh" ref
echo "== parallel (-j 0) =="
RP_EXTRA="-j 0" BASE="$FIXTURE" "$HERE/run.sh" chk

# A pass means nothing if the parallel run never split. The parent says so in
# its log, which run.sh keeps as *_rp.log.
if ! grep -rqi "Splitting the pre-pass across" "$HERE/out_chk"; then
    echo "FAIL: the -j 0 run did not split - this fixture proves nothing."
    echo "      Use one whose pre-pass window spans enough averaging periods."
    exit 1
fi

status=0
while IFS= read -r f; do
    rel="${f#"$HERE/out_ref/"}"
    case "$rel" in *_rp.log) continue ;; esac
    if ! cmp -s "$f" "$HERE/out_chk/$rel"; then
        echo "DIFFERS: $rel"
        status=1
    fi
done < <(find "$HERE/out_ref" -type f)

# Catch a file that exists on only one side, which cmp above cannot see.
a="$(cd "$HERE/out_ref" && find . -type f ! -name "*_rp.log" | sort)"
b="$(cd "$HERE/out_chk" && find . -type f ! -name "*_rp.log" | sort)"
if [ "$a" != "$b" ]; then
    echo "FAIL: the two runs did not write the same set of files"
    diff <(echo "$a") <(echo "$b") || true
    status=1
fi

n="$(echo "$a" | grep -c . || true)"
if [ "$status" -eq 0 ]; then
    echo "IDENTICAL across $n files (run log excluded)"
else
    echo "PARALLEL RUN DIFFERS FROM SERIAL"
fi
exit "$status"
