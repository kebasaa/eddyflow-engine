#!/usr/bin/env bash
# Is a slow gas's time lag detected between its own samples?
#
# Usage: check_slow_lag.sh            (needs gen_slow_drift.py run first)
#
# Runs base_slow_drift, base_slow_drift_sharp and base_slow_drift_integr through
# run.sh and reads the COS time lag of every period from the full output. The
# COS there is made from w delayed by exactly 3.35 s and sampled on a drifting
# 0.987 Hz grid, so every period must report 3.35 +/- 0.05 s - 3.3 or 3.4 at
# 10 Hz rows. A method that only looks at the gas's own rate reports a multiple
# of about one second; the build before two-stage detection fell back on every
# period and reported 8.1 or 10.0 s.
set -uo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
PY="${PYTHON:-/c/Users/jonmuell/AppData/Local/miniconda3/python.exe}"
status=0
for fx in base_slow_drift base_slow_drift_sharp base_slow_drift_integr; do
    if ! BASE="$fx.eddyflow" "$HERE/run.sh" chk > /dev/null 2>&1; then
        echo "FAIL: $fx did not run"; status=1; continue
    fi
    f="$(ls "$HERE"/out_chk/*full_output*.csv | head -1)"
    if "$PY" - "$f" "$fx" <<'PYEOF'
import csv, sys
rows = list(csv.reader(open(sys.argv[1], encoding='utf-8', errors='replace')))
i = rows[1].index('cos_time_lag')
lags = [float(r[i]) for r in rows[3:]]
bad = [x for x in lags if abs(x - 3.35) > 0.05 + 1e-9]
print('%-24s cos lags %s  %s' % (sys.argv[2], ' '.join('%.2f' % x for x in lags), 'OK' if lags and not bad else 'WRONG'))
sys.exit(0 if lags and not bad else 1)
PYEOF
    then :; else status=1; fi
done
[ "$status" -eq 0 ] && echo "SLOW-GAS LAG RECOVERED IN EVERY PERIOD" || echo "SLOW-GAS LAG NOT RECOVERED"
exit "$status"
