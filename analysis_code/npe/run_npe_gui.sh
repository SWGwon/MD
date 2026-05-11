#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VENV_PY="$SCRIPT_DIR/.venv/bin/python"
PYTHON_CMD="python3"

if [[ ! -x "$VENV_PY" ]]; then
    echo "Virtual environment not found. Creating it first..."
    "$SCRIPT_DIR/setup_env.sh" || true
fi

if [[ -x "$VENV_PY" ]]; then
    PYTHON_CMD="$VENV_PY"
else
    echo "Using system Python because no virtual environment is available."
fi

if "$PYTHON_CMD" - <<'PY'
import importlib.util
raise SystemExit(0 if importlib.util.find_spec("tkinter") else 1)
PY
then
    exec "$PYTHON_CMD" "$SCRIPT_DIR/run_npe_gui.py"
else
    echo "Tkinter is not available. Starting the browser-based GUI instead."
    echo "Install python3-tk later if you want the desktop GUI."
    exec "$PYTHON_CMD" "$SCRIPT_DIR/run_npe_web.py"
fi
