#!/usr/bin/env bash
# Run every fixture and report which ones broke.
#
# Usage: sweep.sh [ref|chk]        (default: chk)
#        BIN=/path/to/bin sweep.sh chk
#
# run.sh takes one fixture per invocation, which is right for diffing a single
# case but meant nothing ever ran the whole set. base_no_gas sat failing at the
# FCC stage for an unknown length of time because of it: RP wrote a
# column-aligned file, FCC rejected every record, and no one was looking. The
# fixture was not even in the README table.
#
# Two gates per fixture, both of which that failure would have tripped:
#
#   1. run.sh exits zero - RP *and* FCC completed. RP alone is not enough; FCC
#      recomputes the fluxes under fcc_follows, so an RP-only pass compares
#      files nothing wrote.
#   2. check_columns.py passes - every row matches its header.
#
# Gate 2 alone would not have caught base_no_gas: header and row were short by
# the same block and so agreed with each other. That is why gate 1 exists.
#
# This does not diff against a reference. Use run.sh ref/chk for that, on the
# one fixture whose behaviour you are changing.

HERE="$(cd "$(dirname "$0")" && pwd)"
WHICH="${1:-chk}"
shift 2>/dev/null || true
#> Any names after ref|chk narrow the run to those fixtures, which is what you
#> want while fixing one of them - the full set is a good many minutes.
NAMED="$*"
BIN="${BIN:-/c/Users/jonmuell/Documents/GitHub/build/eddyflow-engine-win-release/bin}"
PY="${PY:-/c/Users/jonmuell/AppData/Local/miniconda3/python.exe}"

# Every fixture that is meant to run to completion. The ones deliberately left
# out and why:
#   base.eddyflow            the pre-record legacy project, kept for migration
#   base_dup                 the same species twice, diffed against base_5gas
#   base_neg                 expects one gas's columns to move; diff it
#   base_n_gas_sa_bad        expects the assessment to be rejected
#   base_drift/base_dynmd    need their own .metadata alongside
#
# The base_slow* set was built for the per-instrument rate work and then never
# registered here, so the one path where a column samples under the row rate had
# no gate at all - which is how a spike counter that counted rows instead of
# samples went unnoticed. They read their data from the '_slow' sibling
# directory gen_slow.py writes.
#
# base_ghg is the same three hours as LI-COR .ghg archives - the compressed
# format, which had no fixture at all. It is the only path that unzips (five
# shell invocations and a decompress per file, about 350 ms) and the only one
# that rewrites Metadata%ac_freq and NumCol per file, since each archive
# carries its own metadata. Built by gen_ghg.py; needs 7-Zip, and is skipped
# rather than failed without it.
#
# base_tlag_par is base_tlag_opt with a two-day time-lag optimisation window
# instead of a three-hour one. That is the only fixture here whose pre-pass is
# long enough for the engine to split it across worker processes: every other
# one covers too few averaging periods to be worth starting a process for, so
# without this the parallel path had no gate at all. It costs about a minute.
#
# base_pwb_cache is base_rec with to_mode=1, which is PWB cache generation:
# walk every period first, then decide every time lag at once from the
# finished table. 39 fixtures configure PWB and not one of them set to_mode=1,
# so PostProcessPwbTimelagCache - the routine that actually settles every lag,
# and the S1/S2/share/interpolate/back-fill/carry/median ladder inside it - had
# no coverage whatsoever. Three hours is enough to exercise the ladder: co2 and
# h2o carry forward, cos borrows across the analyser.
#
# base_pwb_prefilt is base_pwb_cache with the HDI pre-filter tightened to
# 0.10 s, which discards every detection every gas made. That drives all 21
# rows to the terminal arm of PostProcessPwbTimelagCache - the one labelled
# maxcov_default - which nothing else in this suite reaches at all: on
# base_pwb_cache every gas reports fallback=0.
#
# It is the case that separates a per-period covariance maximum from a
# carried lag wearing its label. Before that was fixed, consecutive periods
# came back with the SAME "maxcov_default" lag - 18.8 s three times running
# for h2o - which a per-period maximum cannot do.
#
# base_pwb_par is base_pwb_cache over twelve hours instead of three, which is
# the only PWB fixture long enough for the engine to split its pre-pass:
# PlanPrepassBatches wants nPeriods/4 workers, so 24 periods gets 6 and 7 gets
# none. sweep.sh passes no -j, so the engine takes its default of every core -
# which means this fixture runs the PARALLEL path here, and base_pwb_cache the
# serial one, on the same code.
#
# That they agree is not gated here. check_parallel.sh does that.
#
# Keep this list free of comments: it is a word-split string, not shell source,
# so a '#' line inside it becomes four or five bogus fixture names.
FIXTURES="
base_rec
base_cec_cos
base_no_gas
base_no_water
base_no_ch4
base_5gas
base_n_gas
base_n_gas_cell
base_n_gas_ru
base_h2o_late
base_cell_ref
base_biomet_water
base_biomet_rh
base_tlag_opt
base_tlag_par
base_pwb_cache
base_pwb_prefilt
base_pwb_par
base_ghg
base_auto_sa
base_mw
base_mw_ref
base_slow
base_slow_naive
base_slow_lack
base_slow_integr
base_n_gas_sa_partial
base_ep
base_ep_native
"

[ -n "$NAMED" ] && FIXTURES="$NAMED"

pass=0
fail=0
failed=""

for f in $FIXTURES; do
    #> A fixture is normally a .eddyflow. base_ep is an EddyPro project,
    #> which run.sh hands to the engine as it stands so the import runs for
    #> real; the extension is how run.sh knows to look for the imported file
    #> afterwards, and this is the only place that path runs end to end.
    if [ -f "$HERE/$f.eddyflow" ]; then
        fixture="$f.eddyflow"
    elif [ -f "$HERE/$f.eddypro" ]; then
        fixture="$f.eddypro"
    else
        printf '%-22s SKIP  (no such fixture)\n' "$f"; continue
    fi
    #> The GHG fixture is the only one that needs a tool the engine shells
    #> out to. Without 7-Zip its archives cannot be opened and the run
    #> produces nothing, which would read as a code failure - so it is
    #> skipped rather than failed, and said out loud.
    if [ "$f" = base_ghg ] && ! command -v 7z >/dev/null 2>&1 \
            && ! command -v 7za >/dev/null 2>&1; then
        printf '%-22s SKIP  (7-Zip not on PATH)\n' "$f"; continue
    fi
    rm -rf "$HERE/out_$WHICH"
    #> Log outside the output directory: run.sh clears that directory as its
    #> first act, so a redirect into it has nowhere to land.
    log="${TMPDIR:-/tmp}/sweep_$f.log"
    if ! BIN="$BIN" BASE="$fixture" bash "$HERE/run.sh" "$WHICH" \
            > "$log" 2>&1; then
        printf '%-22s FAIL  (run.sh exited non-zero)\n' "$f"
        tail -3 "$HERE/out_$WHICH/_fcc.log" 2>/dev/null \
            || tail -3 "$log"
        fail=$((fail + 1)); failed="$failed $f"
        continue
    fi
    cols="$("$PY" "$HERE/check_columns.py" "$HERE/out_$WHICH" 2>&1 | tail -1)"
    case "$cols" in
        *"every row matches"*)
            printf '%-22s ok    %s\n' "$f" "$cols"
            pass=$((pass + 1)) ;;
        *)
            printf '%-22s FAIL  %s\n' "$f" "$cols"
            fail=$((fail + 1)); failed="$failed $f" ;;
    esac
done

echo
echo "$pass passed, $fail failed"
if [ "$fail" -gt 0 ]; then
    echo "failed:$failed"
    exit 1
fi
