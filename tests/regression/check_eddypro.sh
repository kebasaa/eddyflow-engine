#!/usr/bin/env bash
# Does an extended .ghg still process in EDDYPRO, and do both engines agree?
#
# Usage: check_eddypro.sh
#        EDDYPRO_BIN=/path/to/eddypro_rp.exe check_eddypro.sh
#
# The whole point of the extended format is that one archive processes in both
# programs. EddyFlow's half of that is gated by sweep.sh - base_ghg_ext and
# base_ghg_campbell. EddyPro's half was gated by nothing: every result claimed
# for it, including "189 of 189 columns bit-identical", came from hand-runs in a
# scratch directory. That made the compatibility claim an anecdote, and nothing
# would have noticed if a change to gen_ghg_ext.py started producing archives
# EddyPro refuses.
#
# Not a sweep fixture, for the same reason check_parallel.sh is not one: this
# needs a second program, installed at a machine-specific path, and each of its
# runs takes about half a minute. It skips rather than fails wherever EddyPro,
# 7-Zip or the generated archives are missing.
#
# Three things about invoking EddyPro, each learnt from an actual error:
#
#   * `-e <home>/` must come BEFORE the project path. EddyPro's argument parser
#     lets the path swallow the next token, so the flag after it is eaten.
#   * That home has to be writable, and the App-V install root is not.
#   * A real 7z.exe must be on PATH. EddyPro shells out to cmd.exe to unzip,
#     and that child cannot see EddyPro's own virtualised copy.
set -uo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
PY="${PY:-/c/Users/jonmuell/AppData/Local/miniconda3/python.exe}"
EDDYPRO_BIN="${EDDYPRO_BIN:-/c/Users/jonmuell/AppData/Local/Microsoft/AppV/Client/Integration/86FA278F-C789-4657-A2BD-8BF9A9AA45C1/Root/EddyPro7/bin/eddypro_rp.exe}"

skip() { echo "SKIP: $*"; exit 0; }
fail() { echo "FAIL: $*"; exit 1; }

[ -x "$EDDYPRO_BIN" ] || skip "no EddyPro at $EDDYPRO_BIN (set EDDYPRO_BIN)"
command -v 7z >/dev/null 2>&1 || command -v 7za >/dev/null 2>&1 \
    || skip "7-Zip not on PATH - EddyPro shells out to it to unzip"
for d in data_ghg data_ghg_ext data_ghg_campbell; do
    [ -d "$HERE/$d" ] || skip "$d/ absent - run gen_ghg_ext.py first"
done

WORK="$(mktemp -d)"
#> KEEP=1 leaves the work directory behind - both engines' outputs, the
#> generated projects and every log - which is what you want the first time a
#> column shows up in the report and you need to see the two numbers.
if [ -z "${KEEP:-}" ]; then
    trap 'rm -rf "$WORK"' EXIT
else
    trap 'echo; echo "kept: $WORK"' EXIT
fi
mkdir -p "$WORK/home"

status=0

# --------------------------------------------------------------------------
# 0. Do the archives actually carry what makes them interesting?
#
# Everything below compares EddyPro against EddyPro. If the generator produced
# archives that are simply copies, step 1 compares a file with itself and
# passes while proving nothing. That is not hypothetical: the ef_model negative
# control written earlier in this work silently tested nothing - a relative
# path handed to `7z u` with cwd elsewhere created a new archive instead of
# updating the intended one, and the "control" re-ran identical data. It looked
# like a clean pass.
# --------------------------------------------------------------------------
echo "== 0. the archives carry the extension =="
zbin=$(command -v 7z || command -v 7za)
check_keys() {
    local dir="$1"; shift
    local ghg tmp md
    ghg=$(ls "$HERE/$dir"/*.ghg 2>/dev/null | head -1)
    [ -n "$ghg" ] || fail "$dir/ holds no archive"
    tmp="$WORK/peek_$dir"; mkdir -p "$tmp"
    md="$(basename "$ghg" .ghg).metadata"
    "$zbin" e "$ghg" -o"$tmp" "$md" -y >/dev/null 2>&1 \
        || fail "cannot read $md out of $dir/"
    for want in "$@"; do
        grep -q "$want" "$tmp/$md" \
            || fail "$dir/ metadata has no '$want' - the fixture is not extended, so nothing below would test anything"
    done
    echo "   $dir: $* present"
}
check_keys data_ghg_ext      ef_model generic_open_path
check_keys data_ghg_campbell ef_model csat3b
#> And the control must NOT be extended, or step 1 compares like with like.
if "$zbin" e "$(ls "$HERE"/data_ghg/*.ghg | head -1)" \
        -o"$WORK/peek_ctrl" '*.metadata' -y >/dev/null 2>&1 \
   && grep -rqi "ef_model" "$WORK/peek_ctrl" 2>/dev/null; then
    fail "data_ghg/ is itself extended - there is no control to compare against"
fi
echo "   data_ghg: plain, as the control must be"

# --------------------------------------------------------------------------
# EddyPro runs. The project is base_ep_licor.eddypro, which is a real EddyPro
# project already committed for the importer fixture - run_mode=1, file_type=0,
# use_pfile=0 - with only the paths rewritten. One committed file describes the
# run for both programs, which is the property this whole test is about.
# --------------------------------------------------------------------------
run_eddypro() {
    local name="$1" data="$2" sonic="${3:-}"
    local prj="$WORK/$name.eddypro" out="$WORK/out_$name"
    mkdir -p "$out"
    sed -e "s|^project_id=.*|project_id=$name|" \
        -e "s|^data_path=.*|data_path=$(cygpath -m "$HERE/$data")|" \
        -e "s|^out_path=.*|out_path=$(cygpath -m "$out")|" \
        "$HERE/base_ep_licor.eddypro" > "$prj"
    #> master_sonic names the instrument by MODEL STRING, so for a stand-in it
    #> has to be the stand-in's name - the only spelling EddyPro resolves.
    #> Wrong here, EddyPro reports missing wind and never mentions instruments.
    [ -n "$sonic" ] && sed -i "s|^master_sonic=.*|master_sonic=$sonic|" "$prj"
    rm -rf "${WORK:?}/home/tmp"
    ( cd "$(dirname "$EDDYPRO_BIN")" \
      && "./$(basename "$EDDYPRO_BIN")" -e "$(cygpath -w "$WORK/home")/" \
             "$(cygpath -w "$prj")" ) > "$WORK/$name.log" 2>&1
    local code=$?
    local csv
    csv=$(ls "$out"/*full_output*.csv 2>/dev/null | head -1)
    if [ $code -ne 0 ] || [ -z "$csv" ]; then
        echo "   $name: REFUSED (exit $code)"
        grep -iE "warning\(|error\(" "$WORK/$name.log" | head -4 | sed 's/^/      /'
        return 1
    fi
    #> Periods, not just an exit code: EddyPro can finish cleanly having
    #> processed nothing at all.
    local n=$(( $(grep -c "" "$csv") - 3 ))
    [ "$n" -gt 0 ] || { echo "   $name: ran but produced 0 periods"; return 1; }
    echo "   $name: ok, $n period(s)"
    echo "$csv" > "$WORK/$name.csvpath"
    return 0
}

echo
echo "== 1. EddyPro: the extension must cost nothing =="
run_eddypro ctrl data_ghg || status=1
run_eddypro ext  data_ghg_ext || status=1
if [ -f "$WORK/ctrl.csvpath" ] && [ -f "$WORK/ext.csvpath" ]; then
    "$PY" "$HERE/check_engines.py" \
        "$(cat "$WORK/ctrl.csvpath")" "$(cat "$WORK/ext.csvpath")" \
        --mode identical --label-a "EddyPro classic" --label-b "EddyPro extended" \
        || status=1
else
    status=1
fi

echo
echo "== 2. EddyPro: the renamed Campbell must process =="
#> Not compared against the control: the sonic genuinely changed from a Metek
#> to a CSAT3B, so the numbers legitimately move. The claim here is only that
#> EddyPro reads the file and produces fluxes from it.
if run_eddypro campbell data_ghg_campbell csat3b_1; then
    "$PY" - "$(cat "$WORK/campbell.csvpath")" <<'PY' || status=1
import io, sys
rows = [l.rstrip('\n').split(',') for l in
        io.open(sys.argv[1], encoding='cp1252', errors='replace')]
head, data = rows[1], rows[3:]
bad = []
for name in ('Tau', 'H', 'u*'):
    if name not in head:
        bad.append('%s: absent from the output' % name); continue
    vals = [r[head.index(name)] for r in data if head.index(name) < len(r)]
    if all(v.strip() in ('-9999', '-9999.0', 'NaN', '') for v in vals):
        bad.append('%s: every period missing' % name)
if bad:
    print('  FAIL: the sonic produced nothing -')
    for b in bad:
        print('    ' + b)
    sys.exit(1)
print('  ok: Tau, H and u* all present')
PY
else
    status=1
fi

# --------------------------------------------------------------------------
# 3. Both engines on the same archive.
#
# This is what "processes correctly in both" means beyond "neither crashes".
# check_engines.py gates the columns that must agree, allows a stated
# tolerance on the humidity chain, and merely reports the methane columns -
# which differ because EddyPro leaves the LI-7500's signal-strength slot NaN
# and EddyFlow fills it. Gating those would be gating the bug.
# --------------------------------------------------------------------------
echo
echo "== 3. EddyFlow against EddyPro, same archives AND same project =="
#> The SAME project file goes to both programs: EddyPro runs it natively,
#> EddyFlow imports it. Anything else is not a comparison. Handing EddyPro its
#> express project and EddyFlow the base_ghg_ext .eddyflow compares two
#> different processing configurations - different detrending, different
#> corrections, different tests - and reports differences of 15 % that say
#> nothing whatever about the archive format.
#>
#> The generated projects are named _cmp_*.eddypro, which .gitignore's
#> tests/regression whitelist does not cover, so they stay untracked.
for ep in ext campbell; do
    [ -f "$WORK/$ep.csvpath" ] || { echo "   $ep: no EddyPro run to compare"; status=1; continue; }
    cmp_prj="_cmp_$ep.eddypro"
    cp "$WORK/$ep.eddypro" "$HERE/$cmp_prj"
    rm -rf "$HERE/out_chk"
    if ! BASE="$cmp_prj" bash "$HERE/run.sh" chk > "$WORK/ef_$ep.log" 2>&1; then
        echo "   $ep: EddyFlow run failed"
        grep -iE "warning\(|error\(" "$WORK/ef_$ep.log" | head -3 | sed 's/^/      /'
        rm -f "$HERE/$cmp_prj"; status=1; continue
    fi
    rm -f "$HERE/$cmp_prj"
    efcsv=$(ls "$HERE/out_chk"/*full_output*.csv 2>/dev/null | head -1)
    [ -n "$efcsv" ] || { echo "   $ep: EddyFlow wrote no full_output"; status=1; continue; }
    "$PY" "$HERE/check_engines.py" "$(cat "$WORK/$ep.csvpath")" "$efcsv" \
        --mode engines --label-a "EddyPro $ep" --label-b "EddyFlow $ep" \
        || status=1
done

echo
if [ "$status" -eq 0 ]; then
    echo "PASS: the extended archives process in both engines"
else
    echo "FAIL: see above"
fi
exit "$status"
