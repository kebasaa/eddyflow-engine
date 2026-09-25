#!/usr/bin/env bash
# FCC reading RP's results from a shared link must give the same results as
# reading them from the folder.
#
# Usage: BIN=<engine bin dir> check_remote_fcc.sh BASE
#
# Google Drive links only; the Dropbox listing is covered by check_remote.sh.
#
# RP runs once, from disk, into remote_tmp/fcc/served. FCC then runs twice
# over it: once with ex_file, sa_bin_spectra and sa_full_spectra naming the
# folder, once naming links to it served by remote_server.py. The FCC outputs
# must be identical once the run timestamp is taken out of them.
#
# RP's output is served as written. The harness's normalised out_ref will not
# do: its file names have lost the timestamps FCC's file template expects.
set -euo pipefail

BASE="${1:?usage: check_remote_fcc.sh BASE}"
HERE="$(cd "$(dirname "$0")" && pwd)"
: "${BIN:?set BIN to the engine bin directory}"
PY="${PY:-/c/Users/jonmuell/AppData/Local/miniconda3/python.exe}"
export PATH="/c/Users/jonmuell/Documents/GitHub/eddyflow-portable/bin:/c/Users/jonmuell/mingw64/bin:$PATH"

T="$HERE/remote_tmp/fcc"
HOME_DIR="$HERE/remote_tmp/home"
rm -rf "$T" "$HOME_DIR"
mkdir -p "$T/served" "$T/outA" "$T/outB" "$HOME_DIR/tmp" "$HOME_DIR/ini"
M="$(cygpath -m "$T")"
esc() { printf '%s' "$1" | sed -e 's/[&|\\]/\\&/g'; }

# RP, from disk
sed -e "s|^out_path=.*|out_path=$M/served|" "$HERE/$BASE" > "$T/rp.eddyflow"
"$BIN/eddyflow_rp.exe" "$(cygpath -w "$T/rp.eddyflow")" -e "$(cygpath -w "$HOME_DIR")/" \
    > "$T/rp.log" 2>&1 || { echo "RP failed"; tail -20 "$T/rp.log"; exit 1; }
EX="$(cd "$T/served" && ls ./*fluxnet*.csv | head -1)"
EX="${EX#./}"
[ -n "$EX" ] || { echo "RP wrote no fluxnet file"; exit 1; }

"$PY" "$HERE/remote_server.py" --root "$(cygpath -w "$T/served")" \
    --port-file "$(cygpath -w "$T/port")" --lifetime 900 &
SERVER=$!
trap 'kill $SERVER 2>/dev/null || true' EXIT
for _ in $(seq 50); do [ -s "$T/port" ] && break; sleep 0.2; done
[ -s "$T/port" ] || { echo "server did not start"; exit 1; }
PORT="$(cat "$T/port")"
lk() { "$PY" -c "import json,sys; print(json.load(open(sys.argv[1]))[sys.argv[2]][sys.argv[3]])" \
    "$(cygpath -w "$T/port.links")" "$@"; }
EXL="$(lk gdrive_files "$EX")"
BINL="$(lk gdrive_folders eddyflow_binned_cospectra)"
FULL="$(lk gdrive_folders eddyflow_full_cospectra)"

# The run from disk reads its own copy: FCC deletes the ex file it read
# unless keep_parent is set, and must not delete the one being served.
cp -r "$T/served" "$T/local"
sed -e "s|^out_path=.*|out_path=$M/outA|" \
    -e "s|^ex_file=.*|ex_file=$M/local/$EX|" \
    -e "s|^sa_bin_spectra=.*|sa_bin_spectra=$M/local/eddyflow_binned_cospectra|" \
    -e "s|^sa_full_spectra=.*|sa_full_spectra=$M/local/eddyflow_full_cospectra|" \
    "$HERE/$BASE" > "$T/a.eddyflow"
sed -e "s|^out_path=.*|out_path=$M/outB|" \
    -e "s|^ex_file=.*|ex_file=$(esc "$EXL")|" \
    -e "s|^sa_bin_spectra=.*|sa_bin_spectra=$(esc "$BINL")|" \
    -e "s|^sa_full_spectra=.*|sa_full_spectra=$(esc "$FULL")|" \
    "$HERE/$BASE" > "$T/b.eddyflow"

"$BIN/eddyflow_fcc.exe" "$(cygpath -w "$T/a.eddyflow")" -e "$(cygpath -w "$HOME_DIR")/" \
    > "$T/a.log" 2>&1 || { echo "FCC from disk failed"; tail -20 "$T/a.log"; exit 1; }
EDDYFLOW_REMOTE_BASE="http://127.0.0.1:$PORT" \
    "$BIN/eddyflow_fcc.exe" "$(cygpath -w "$T/b.eddyflow")" -e "$(cygpath -w "$HOME_DIR")/" \
    > "$T/b.log" 2>&1 || { echo "FCC from links failed"; tail -20 "$T/b.log"; exit 1; }

norm() {
    find "$1" -type f \( -name '*.csv' -o -name '*.txt' \) | sort | while IFS= read -r f; do
        echo "== ${f#$1/}" | sed -E 's/_?[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{6}//g'
        sed -E 's/[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{6}/TIMESTAMP/g' "$f"
    done
}
echo "== $BASE, FCC inputs via gdrive =="
grep 'shared link' "$T/b.log" || true
echo "output files: disk $(find "$T/outA" -type f | wc -l), links $(find "$T/outB" -type f | wc -l)"
if diff <(norm "$T/outA") <(norm "$T/outB") > "$T/diff.txt"; then
    echo "IDENTICAL"
else
    echo "DIFFERENT"; head -20 "$T/diff.txt"; exit 1
fi
