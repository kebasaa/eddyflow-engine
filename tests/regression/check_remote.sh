#!/usr/bin/env bash
# Reading a fixture from a shared link must give the same results as reading
# it from the folder.
#
# Usage: BIN=<engine bin dir> check_remote.sh BASE gdrive|dropbox [server args]
#
# The fixture's raw data folder, and its metadata file if it names one, are
# copied into remote_tmp/data and served by remote_server.py the way Google
# Drive or Dropbox serve a shared link. The fixture then runs twice through
# run.sh: `ref` reading the copy from disk, `chk` with data_path and proj_file
# replaced by links and EDDYFLOW_REMOTE_BASE pointing the engine at the server.
# Every output file except the run log and the copy of the project must be
# identical; those two legitimately differ, since they name the link.
#
# Every raw file must cross the network once, however many passes read it: the
# server logs each file it serves whole, and a file served twice fails the
# check. EXTEND=N replaces the data with N consecutive half-hours, copies of
# the fixture's own archives under new names (the engine places a file by its
# name), with the processing and time-lag subsets switched off - enough periods
# for `RP_EXTRA="-j 4"` to split a pre-pass over workers, which need four
# periods each.
#
# DEDUPE=1 also names the shared folder as head_corr_dir, which the engine
# downloads whole at startup and these fixtures otherwise leave unused. Every
# raw file is then named by two settings and must still be transferred once.
#
# NEST=1 moves the files into two subfolders, for the recursive listing.
# EXTRACT_META=1 takes the .metadata out of the first GHG archive and uses it
# as the alternative metadata file (use_pfile=1), for a remote proj_file.
#
# Server args inject failures, e.g. `--missing <file>` or `--html <file>`,
# and `--page 1` makes the Dropbox listing page through its entries one by one.
set -euo pipefail

BASE="${1:?usage: check_remote.sh BASE gdrive|dropbox [server args]}"
PROVIDER="${2:?usage: check_remote.sh BASE gdrive|dropbox [server args]}"
shift 2
HERE="$(cd "$(dirname "$0")" && pwd)"
: "${BIN:?set BIN to the engine bin directory}"
PY="${PY:-/c/Users/jonmuell/AppData/Local/miniconda3/python.exe}"
# 7-Zip for the GHG fixtures; it is not on PATH in Git Bash
export PATH="/c/Users/jonmuell/Documents/GitHub/eddyflow-portable/bin:/c/Users/jonmuell/mingw64/bin:$PATH"

TMP="$HERE/remote_tmp"
rm -rf "$TMP"; mkdir -p "$TMP/data"

value() { grep -m1 "^$1=" "$HERE/$BASE" | cut -d= -f2- | tr -d '\r'; }
DATA="$(value data_path)"
cp -r "$(cygpath -u "$DATA")"/. "$TMP/data/"
EXTRA_SED=()
if [ "${EXTEND:-0}" -gt 0 ]; then
    mapfile -t templates < <(find "$TMP/data" -maxdepth 1 -name '*.ghg' | sort)
    [ "${#templates[@]}" -gt 0 ] || { echo "EXTEND needs .ghg files"; exit 2; }
    day="$(basename "${templates[0]}" | cut -c1-10)"
    suffix="$(basename "${templates[0]}" | cut -d_ -f2-)"
    mkdir -p "$TMP/extended"
    for i in $(seq 0 $((EXTEND - 1))); do
        name="$(printf '%sT%02d%02d00_%s' "$day" $((i / 2)) $(((i % 2) * 30)) "$suffix")"
        cp "${templates[$((i % ${#templates[@]}))]}" "$TMP/extended/$name"
    done
    rm -f "$TMP"/data/*.ghg
    mv "$TMP"/extended/*.ghg "$TMP/data/"
    EXTRA_SED=(-e "s|^pr_subset=.*|pr_subset=0|" -e "s|^to_subset=.*|to_subset=0|"
               -e "s|^pf_subset=.*|pf_subset=0|")
fi
if [ "${NEST:-0}" = "1" ]; then
    mkdir -p "$TMP/data/a" "$TMP/data/b/c"
    i=0
    for f in "$TMP"/data/*.*; do
        if [ $((i % 2)) = 0 ]; then mv "$f" "$TMP/data/a/"; else mv "$f" "$TMP/data/b/c/"; fi
        i=$((i + 1))
    done
fi
META=""
if [ "${EXTRACT_META:-0}" = "1" ]; then
    META=site.metadata
    mkdir -p "$TMP/data/meta"
    first="$(find "$TMP/data" -name '*.ghg' | sort | head -1)"
    7z e -so "$(cygpath -w "$first")" '*.metadata' '-x!*-biomet.metadata' > "$TMP/data/meta/$META"
elif [ "$(value use_pfile)" = "1" ]; then
    META="$(basename "$(value proj_file)")"
    mkdir -p "$TMP/data/meta"
    cp "$(cygpath -u "$(value proj_file)")" "$TMP/data/meta/$META"
fi
WDATA="$(cygpath -m "$TMP/data")"

"$PY" "$HERE/remote_server.py" --root "$(cygpath -w "$TMP/data")" \
    --port-file "$(cygpath -w "$TMP/port")" --lifetime 1800 \
    --log "$(cygpath -w "$TMP/served.log")" "$@" &
SERVER=$!
trap 'kill $SERVER 2>/dev/null || true' EXIT
for _ in $(seq 50); do [ -s "$TMP/port" ] && break; sleep 0.2; done
[ -s "$TMP/port" ] || { echo "server did not start"; exit 1; }
PORT="$(cat "$TMP/port")"

link() { "$PY" -c "import json,sys; d=json.load(open(sys.argv[1])); print(d[sys.argv[2]] if len(sys.argv)<4 else d[sys.argv[2]][sys.argv[3]])" "$(cygpath -w "$TMP/port.links")" "$@"; }
case "$PROVIDER" in
    gdrive)  ROOT_LINK="$(link gdrive_root)"; META_LINK="${META:+$(link gdrive_files "meta/$META")}";;
    dropbox) ROOT_LINK="$(link dropbox_root)"; META_LINK="${META:+$(link dropbox_files "meta/$META")}";;
    *) echo "provider must be gdrive or dropbox"; exit 2;;
esac

# sed replacement text: & and | are special
esc() { printf '%s' "$1" | sed -e 's/[&|\\]/\\&/g'; }
sed -e "s|^data_path=.*|data_path=$WDATA|" "${EXTRA_SED[@]}" \
    ${META:+-e "s|^proj_file=.*|proj_file=$WDATA/meta/$META|" -e "s|^use_pfile=.*|use_pfile=1|"} \
    "$HERE/$BASE" > "$HERE/remote_ref.eddyflow"
CHK_SED=()
if [ "${DEDUPE:-0}" = "1" ]; then
    CHK_SED=(-e "/^data_path=/a head_corr_dir=\"$(esc "$ROOT_LINK")\"")
fi
# Quoted, as the interface's QSettings writes any value containing = ; or ,
sed -e "s|^data_path=.*|data_path=\"$(esc "$ROOT_LINK")\"|" "${EXTRA_SED[@]}" "${CHK_SED[@]}" \
    ${META:+-e "s|^proj_file=.*|proj_file=\"$(esc "$META_LINK")#path=meta/$META\"|" -e "s|^use_pfile=.*|use_pfile=1|"} \
    "$HERE/$BASE" > "$HERE/remote_chk.eddyflow"

BASE=remote_ref.eddyflow bash "$HERE/run.sh" ref > "$TMP/ref.out" 2>&1 \
    || { echo "ref run failed"; tail -20 "$TMP/ref.out"; exit 1; }
EDDYFLOW_REMOTE_BASE="http://127.0.0.1:$PORT" BASE=remote_chk.eddyflow \
    bash "$HERE/run.sh" chk > "$TMP/chk.out" 2>&1 \
    || { echo "chk run failed"; tail -20 "$TMP/chk.out"; exit 1; }

echo "== $BASE via $PROVIDER =="
echo "output files: ref $(find "$HERE/out_ref" -type f | wc -l), chk $(find "$HERE/out_chk" -type f | wc -l)"
if diff -r -q -x '*.log' -x '*.eddyflow' "$HERE/out_ref" "$HERE/out_chk"; then
    echo "IDENTICAL (except the run log and project copy)"
    STATUS=0
else
    echo "DIFFERENT"
    STATUS=1
fi
served="$(sort "$TMP/served.log" 2>/dev/null | uniq -c)"
twice="$(printf '%s\n' "$served" | awk '$1 > 1' | sed '/^$/d')"
echo "downloads: $(printf '%s\n' "$served" | sed '/^$/d' | wc -l) files, $(cat "$TMP/served.log" 2>/dev/null | wc -l) transfers"
if [ -n "$twice" ]; then
    echo "DOWNLOADED MORE THAN ONCE:"; printf '%s\n' "$twice"
    STATUS=1
else
    echo "each file downloaded once"
fi
if grep -h "Warning(121)" "$HERE"/out_chk/*_log*_rp.log >/dev/null 2>&1 && [ "$#" -eq 0 ]; then
    echo "Warning(121) without an injected failure"; STATUS=1
fi
echo "-- run log lines only in chk --"
diff <(cat "$HERE"/out_ref/*_log*_rp.log 2>/dev/null) \
     <(cat "$HERE"/out_chk/*_log*_rp.log 2>/dev/null) | grep '^>' || true
rm -f "$HERE/remote_ref.eddyflow" "$HERE/remote_chk.eddyflow"
exit $STATUS
