#!/usr/bin/env bash
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VENV_DIR="$SCRIPT_DIR/.venv"

command -v python3 >/dev/null 2>&1 || {
    echo "Error: python3 not found in PATH" >&2
    exit 1
}

PYTHON_CMD="python3"

if python3 -m venv "$VENV_DIR"; then
    PYTHON_CMD="$VENV_DIR/bin/python"
    echo "Created virtual environment: $VENV_DIR"
else
    echo
    echo "Warning: could not create a virtual environment."
    echo "This usually means python3-venv is not installed."
    echo
    echo "Ubuntu/Debian:"
    echo "  sudo apt install python3-venv"
    echo
    echo "Continuing with system Python for the dependency check."
fi

"$PYTHON_CMD" - <<'PY'
import importlib.util
import sys

print(f"Python: {sys.executable}")

if importlib.util.find_spec("tkinter") is None:
    print()
    print("Tkinter is not available in this Python.")
    print("The web GUI can still run without extra Python packages.")
    print()
    print("To enable the desktop Tkinter GUI, install the OS package:")
    print("  Ubuntu/Debian: sudo apt install python3-tk")
    print("  Fedora:        sudo dnf install python3-tkinter")
    print("  Arch:          sudo pacman -S tk")
else:
    print("Tkinter: available")
PY

if command -v root >/dev/null 2>&1; then
    echo "ROOT: available ($(command -v root))"
else
    echo "Warning: ROOT command 'root' was not found in PATH."
fi

echo
echo "Run GUI with:"
echo "  ./run_npe_gui.sh"
