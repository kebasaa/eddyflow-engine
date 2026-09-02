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
# base_ghg_licor is two GENUINE LI-COR SmartFlux archives, committed under
# data_ghg/. base_ghg does not cover what its name suggests: gen_ghg.py builds
# it by copying base_site.metadata into a synthesized archive, and it runs
# use_biom=0 - so LI-COR's own embedded metadata was never read by any test,
# and embedded biomet was never exercised at all. That is exactly how a
# missing optional argument in ReadBiometFile's scanCsvFile call - undefined
# behaviour that segfaulted on every embedded-biomet run - survived until it
# was found by hand.
#
# So this one is use_pfile=0 (metadata read from inside each archive, written
# by the LI-7550, not by our interface) and use_biom=1 (30 biomet records per
# period, LI-COR's own tab-separated 6-header-row format). Open-path
# LI-7500A + LI-7700 on a Metek sonic, which is also the suite's only
# open-path gas analyser. Two archives rather than one so it spans two
# averaging periods and re-reads the metadata per file. Needs 7-Zip, same as
# base_ghg.
#
# base_ghg_ext is those same two archives with their embedded metadata rewritten
# into the EXTENDED .ghg format: the LI-7500A declares itself a
# generic_open_path stand-in with real geometry, states its true identity in
# instr_2_ef_model, and every col_N_instrument that named it is repointed at the
# stand-in - which is what an archive describing an analyser EddyPro cannot name
# has to look like.
#
# The gate is an EQUALITY, not just a clean exit: every DATA file must be
# byte-identical to base_ghg_licor's, because ef_model puts the same LI-7500A
# back. The run log legitimately differs by two lines - the data directory, and
# the one announcing the format version - and processing_adv.eddyflow by the
# data_path it echoes.
#
# What it catches is the override being read at all. The stand-in states the
# LI-7500A's TRUE geometry, so a run that ignored ef_model would compute every
# flux to the digit and differ only in what it calls the analyser - two column
# labels, co2_li7500a_1_mean against co2_generic_open_path_1_mean. Verified by
# stripping ef_model back out: four output files move, none of them by a value.
# That is the whole margin, and nothing weaker than a full diff sees it.
#
# sweep.sh runs the structural gates only; the equality is run.sh ref/chk across
# the two fixtures.
#
# Built by gen_ghg_ext.py into data_ghg_ext/, which is not committed - the
# archives are 3.5 MB of binary and the ten lines that differ belong in a script
# you can read. Skipped, like the 7-Zip cases, when it has not been generated.
#
# base_ghg_campbell is the other half of the extended format: a RENAME rather
# than a stand-in. The Metek is swapped for a Campbell CSAT3B, which EddyPro
# spells `csat3b` and EddyFlow spells `csi_csat3b` - the archive states EddyPro's
# spelling, because it is the only one EddyPro takes, and names the real one in
# instr_1_ef_model.
#
# The LI-COR archives cannot show this on their own: every instrument in them is
# spelt the same in both programs, so nothing here reached the canonicalising
# path at all. What lived behind it: col_N_instrument is never canonicalised, so
# an instrument written csat3b_1 and canonicalised to csi_csat3b_1 left all five
# sonic columns pointing at a name the instrument list no longer held. No column
# bound to the sonic, and the run died with "exactly one selected u, v, w and one
# selected ts or sos are required" - a message that does not mention instruments.
#
# So the gate is that it RUNS at all, and that the metadata output names
# csi_csat3b_1. A regression re-breaks column binding and fails outright rather
# than shifting a number.
#
# Its project pins master_sonic=csi_csat3b_1, the resolved name - the third place
# an instrument name propagates, after the instrument block and the columns.
#
# Built by gen_ghg_ext.py into data_ghg_campbell/, same as base_ghg_ext.
#
# base_ep_licor is the same two archives handed over as an EDDYPRO project, so
# run.sh routes it through the importer the way base_ep does. base_ep is the
# only other importer fixture and it is run_mode=0 (advanced), file_type=1
# (generic ASCII), use_pfile=1 (external metadata file). This one is
# run_mode=1 (express), file_type=0 (GHG) and use_pfile=0 (metadata embedded
# in each archive) - three import paths that had no coverage between them,
# including the use_pfile=0 branch where there is no metadata file to resolve
# gas records against and they come out naming no analyser.
#
# The gas columns are stated explicitly rather than left at EddyPro's -1.
# -1 means "resolve automatically", which only happens in embedded mode: a
# SmartFlux project run on the desktop leaves every gas column unresolved, and
# EddyPro answers that with silent NaN gas fluxes. Pinning that would make the
# fixture agree with a broken run, so it names the columns the metadata
# describes - co2 at 10, h2o at 12, ch4 at 41.
#
# Verified against EddyPro 7.0.9 on the full 48-archive sample set before being
# committed: Tau and u* bit-identical, concentrations and time lags
# bit-identical, H within 2e-4, LE within 2.7e-4 - the last matching a
# documented molecular-weight refinement whose changelog predicts ~0.026%.
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
base_ghg_licor
base_ghg_ext
base_ghg_campbell
base_auto_sa
base_mw
base_mw_ref
base_slow
base_slow_naive
base_slow_lack
base_slow_integr
base_n_gas_sa_partial
base_ep
base_ep_licor
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
    #> The GHG fixtures are the only ones that need a tool the engine shells
    #> out to. Without 7-Zip their archives cannot be opened and the run
    #> produces nothing, which would read as a code failure - so they are
    #> skipped rather than failed, and said out loud.
    case "$f" in
        base_ghg|base_ghg_licor|base_ghg_ext|base_ghg_campbell|base_ep_licor)
            if ! command -v 7z >/dev/null 2>&1 \
                    && ! command -v 7za >/dev/null 2>&1; then
                printf '%-22s SKIP  (7-Zip not on PATH)\n' "$f"; continue
            fi ;;
    esac
    #> Generated, not committed. Absent means gen_ghg_ext.py has not been run,
    #> which is a missing input rather than a broken engine - the same reason
    #> the 7-Zip cases skip rather than fail.
    #> Named outright rather than derived from $f. Both directories are
    #> generated by gen_ghg_ext.py and neither is committed, so an absent one
    #> means the generator has not been run - a missing input rather than a
    #> broken engine, which is why these skip rather than fail.
    case "$f" in
        base_ghg_ext)      gen_dir=data_ghg_ext ;;
        base_ghg_campbell) gen_dir=data_ghg_campbell ;;
        *)                 gen_dir= ;;
    esac
    if [ -n "$gen_dir" ] && [ ! -d "$HERE/$gen_dir" ]; then
        printf '%-22s SKIP  (run gen_ghg_ext.py first)
' "$f"; continue
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
