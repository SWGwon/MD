#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MACRO_FILE="$SCRIPT_DIR/plot_npe_subtracted.C"

ok() {
    echo "[OK] $*"
}

warn() {
    echo "[WARN] $*" >&2
}

fail() {
    echo "[FAIL] $*" >&2
    exit 1
}

command -v bash >/dev/null 2>&1 || fail "bash not found in PATH"
ok "bash: $(bash --version | sed -n '1p')"

command -v python3 >/dev/null 2>&1 || fail "python3 not found in PATH"
ok "python3: $(python3 --version)"

python3 - <<'PY'
import importlib.util
import sys

if sys.version_info < (3, 10):
    raise SystemExit("Python 3.10 or newer is required.")

if importlib.util.find_spec("tkinter") is None:
    print("[WARN] tkinter: unavailable; run_npe_gui.sh will use the browser GUI fallback.")
else:
    print("[OK] tkinter: available")
PY

command -v root >/dev/null 2>&1 || fail "ROOT command 'root' not found in PATH"
ok "root: $(command -v root)"

if command -v root-config >/dev/null 2>&1; then
    root_version="$(root-config --version)"
    ok "root-config: $root_version"
    if [[ "$root_version" =~ ^6\.39 ]]; then
        warn "ROOT $root_version looks newer than the stable 6.36 series; verify reproducibility before publishing results."
    fi
else
    warn "root-config not found; ROOT is available but version detection is limited."
fi

[[ -f "$MACRO_FILE" ]] || fail "NPE ROOT macro not found: $MACRO_FILE"

root -l -b -q "$MACRO_FILE(\"/tmp/md_missing_source.root\",\"/tmp/md_missing_bg.root\",1.0e7,-1,1.0,400,0,-1,125.0e6,\"env_check\")" >/tmp/md_npe_root_check.log 2>&1 || {
    cat /tmp/md_npe_root_check.log >&2
    fail "ROOT failed to load the NPE macro."
}
ok "ROOT macro loads: $MACRO_FILE"

echo
echo "Environment check complete."
