#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
START_DIR="$(pwd)"

usage() {
    cat <<'EOF'
Usage:
  ./run_compton_edge_scan.sh -i HISTOGRAMS.root -x MIN -X MAX -o PREFIX [options]

Options:
  -i FILE      Histogram ROOT file from run_compton_edge_analysis.sh
  -x MIN       Base fit x-axis minimum in NPE
  -X MAX       Base fit x-axis maximum in NPE
  -o PREFIX    Output prefix for scan plot, ROOT file, and text summary
  -m FILE      ROOT macro path (default: fit_compton_edge.C next to this script)
  -H NAME      Histogram name (default: hTotalSub, falling back to hTotalSourceOnly)
  -M MODEL     Fit model: erfc_linear or erfc_gaussian (default: erfc_linear)
  -F FRACTION  Fraction of fit width used to shift xmin/xmax (default: 0.10)
  -h           Show this help
EOF
}

die() {
    echo "Error: $*" >&2
    exit 1
}

abs_path() {
    local path="$1"
    if [[ "$path" = /* ]]; then
        printf "%s" "$path"
    else
        printf "%s/%s" "$START_DIR" "$path"
    fi
}

root_string_escape() {
    local value="$1"
    value="${value//\\/\\\\}"
    value="${value//\"/\\\"}"
    printf "%s" "$value"
}

root_string_arg() {
    printf '"%s"' "$(root_string_escape "$1")"
}

hist_file=""
fit_xmin=""
fit_xmax=""
out_prefix=""
macro_file="$SCRIPT_DIR/fit_compton_edge.C"
hist_name=""
model="erfc_linear"
scan_fraction="0.10"

while getopts ":i:x:X:o:m:H:M:F:h" opt; do
    case "$opt" in
        i) hist_file="$OPTARG" ;;
        x) fit_xmin="$OPTARG" ;;
        X) fit_xmax="$OPTARG" ;;
        o) out_prefix="$OPTARG" ;;
        m) macro_file="$OPTARG" ;;
        H) hist_name="$OPTARG" ;;
        M) model="$OPTARG" ;;
        F) scan_fraction="$OPTARG" ;;
        h) usage; exit 0 ;;
        :) die "Option -$OPTARG requires an argument." ;;
        \?) die "Unknown option: -$OPTARG" ;;
    esac
done

[[ -n "$hist_file" ]] || die "Histogram ROOT file is required."
[[ -n "$fit_xmin" ]] || die "Fit x minimum is required."
[[ -n "$fit_xmax" ]] || die "Fit x maximum is required."
[[ -n "$out_prefix" ]] || die "Output prefix is required."

hist_file="$(abs_path "$hist_file")"
macro_file="$(abs_path "$macro_file")"

[[ -f "$hist_file" ]] || die "Histogram ROOT file not found: $hist_file"
[[ -f "$macro_file" ]] || die "ROOT macro not found: $macro_file"
command -v root >/dev/null 2>&1 || die "ROOT command 'root' not found in PATH"

echo
echo "Running Compton edge range scan"
echo "  Macro:         $macro_file"
echo "  Hist file:     $hist_file"
echo "  Hist name:     ${hist_name:-auto}"
echo "  Model:         $model"
echo "  Base fit range: [$fit_xmin, $fit_xmax]"
echo "  Scan fraction: $scan_fraction"
echo "  Prefix:        $out_prefix"
echo

root_call="$(printf 'scan_compton_edge_fit(%s,%s,%s,%s,%s,%s,%s)' \
    "$(root_string_arg "$hist_file")" \
    "$fit_xmin" \
    "$fit_xmax" \
    "$(root_string_arg "$out_prefix")" \
    "$(root_string_arg "$hist_name")" \
    "$(root_string_arg "$model")" \
    "$scan_fraction")"
root -l -b -q -e ".L $(root_string_escape "$macro_file")" -e "$root_call"
