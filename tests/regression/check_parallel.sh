#!/usr/bin/env bash
# Does splitting a pre-pass across worker processes change the answer?
#
# Usage: [PAR_JOBS=N] [PAR_KIND=pre|pr|fx] check_parallel.sh [fixture.eddyflow]
#        (defaults: PAR_JOBS=0 - one worker per core; PAR_KIND=pre - a
#         pre-pass must have been split, on base_tlag_par.eddyflow; with
#         PAR_KIND=pr, base_prod_par.eddyflow, two days cut as a real run is)
#
# PAR_KIND=pr checks the production pass instead - the main period loop that
# computes the fluxes - and asserts that it was split. Most fixtures span a
# day or less, which the engine would cut into a handful of pieces at most;
# set EDDYFLOW_PROD_PIECE_PERIODS=1 (or 2) to cut at nearly every half-hour,
# so that whatever a half-hour inherits from the one before is exercised at
# every cut. The variable only affects a run that splits, so the -j 1
# reference is unchanged by it.
#
# PAR_KIND=fx checks FCC's flux computation, split by essentials record. RP
# runs serially both times; only FCC is given -j. Set
# EDDYFLOW_FCC_PIECE_PERIODS=1 (or 2, 3) to cut at nearly every record. Both
# run logs are left out of the comparison: FCC's parent appends its workers'
# logs to its own.
#
# Runs the fixture twice through run.sh - once with -j 1, once with
# -j $PAR_JOBS - and diffs the two normalised output trees. Every file must
# match. The job count sets how many pieces the range is cut into (about four
# per worker), so running it at a few values checks more than one cut.
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
PAR_JOBS="${PAR_JOBS:-0}"
PAR_KIND="${PAR_KIND:-pre}"
case "$PAR_KIND" in
    pre) MARKER="Splitting the pre-pass across"; DEFAULT=base_tlag_par.eddyflow ;;
    pr)  MARKER="Splitting the production pass across"; DEFAULT=base_prod_par.eddyflow ;;
    fx)  MARKER="Splitting the flux computation across"; DEFAULT=base_prod_par.eddyflow ;;
    *)   echo "PAR_KIND must be pre, pr or fx"; exit 2 ;;
esac
FIXTURE="${1:-$DEFAULT}"

if [ "$PAR_KIND" = fx ]; then
    echo "== serial (FCC -j 1) =="
    RP_EXTRA="-j 1" FCC_EXTRA="-j 1" BASE="$FIXTURE" "$HERE/run.sh" ref
    echo "== parallel (FCC -j $PAR_JOBS) =="
    RP_EXTRA="-j 1" FCC_EXTRA="-j $PAR_JOBS" BASE="$FIXTURE" "$HERE/run.sh" chk
    LOGS="*_log_*"
else
    echo "== serial (-j 1) =="
    RP_EXTRA="-j 1" BASE="$FIXTURE" "$HERE/run.sh" ref
    echo "== parallel (-j $PAR_JOBS) =="
    RP_EXTRA="-j $PAR_JOBS" BASE="$FIXTURE" "$HERE/run.sh" chk
    LOGS="*_rp.log"
fi

# A pass means nothing if the parallel run never split. The parent says so in
# its log, which run.sh keeps as *_rp.log.
if ! grep -rqi "$MARKER" "$HERE/out_chk"; then
    echo "FAIL: the -j $PAR_JOBS run did not split - this fixture proves nothing."
    echo "      Use one whose window spans enough averaging periods."
    exit 1
fi

status=0
while IFS= read -r f; do
    rel="${f#"$HERE/out_ref/"}"
    case "$(basename "$rel")" in $LOGS) continue ;; esac
    if ! cmp -s "$f" "$HERE/out_chk/$rel"; then
        echo "DIFFERS: $rel"
        status=1
    fi
done < <(find "$HERE/out_ref" -type f)

# Catch a file that exists on only one side, which cmp above cannot see.
a="$(cd "$HERE/out_ref" && find . -type f ! -name "$LOGS" | sort)"
b="$(cd "$HERE/out_chk" && find . -type f ! -name "$LOGS" | sort)"
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
