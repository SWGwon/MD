#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
START_DIR="$(pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"

usage() {
    cat <<'EOF'
Usage:
  ./run_npe_analysis.sh
  ./run_npe_analysis.sh -s SOURCE.root -b BACKGROUND.root -o PREFIX [options]

Options:
  -s FILE      Source ROOT file
  -b FILE      Background ROOT file
  -d DIR       Data directory used for interactive file selection
  -o PREFIX    Output PNG prefix
  -O DIR       Output directory for generated PNG files
  -L LABEL     Source label shown in plot legends (default: source file name)
  -m FILE      ROOT macro path (default: plot_npe_subtracted.C next to this script)
  -g GAIN      PMT gain (default: 1.0e7)
  -q QUANTILE  X-axis quantile, 1.0 means full range (default: 1.0)
  -n BINS      Histogram bins (default: 400)
  -x MIN       Manual x-axis minimum (use with -X)
  -X MAX       Manual x-axis maximum (use with -x)
  -c HZ        SyncTime_TTT clock Hz (default: 125.0e6)
  -h           Show this help

Examples:
  ./run_npe_analysis.sh
  ./run_npe_analysis.sh -s source.root -b background.root -o source_run
  ./run_npe_analysis.sh -s source.root -b background.root -o zoom -x 0 -X 200000
  ~/MD/analysis_code/run_npe_analysis.sh -d ~/MD/data/run1 -s source.root -b bg.root -O ~/MD/plots -o source_run
EOF
}

die() {
    echo "Error: $*" >&2
    exit 1
}

pick_root_file() {
    local label="$1"
    shift
    local files=("$@")

    echo
    echo "Select $label file:"
    local i
    for i in "${!files[@]}"; do
        printf "  %2d) %s\n" "$((i + 1))" "${files[$i]}"
    done

    local choice
    while true; do
        read -r -p "$label number: " choice
        if [[ "$choice" =~ ^[0-9]+$ ]] && (( choice >= 1 && choice <= ${#files[@]} )); then
            printf "%s" "${files[$((choice - 1))]}"
            return
        fi
        echo "Enter a number from 1 to ${#files[@]}."
    done
}

prompt_default() {
    local prompt="$1"
    local default="$2"
    local value
    read -r -p "$prompt [$default]: " value
    printf "%s" "${value:-$default}"
}

abs_path() {
    local path="$1"
    if [[ "$path" = /* ]]; then
        printf "%s" "$path"
    else
        printf "%s/%s" "$START_DIR" "$path"
    fi
}

abs_dir() {
    local path="$1"
    mkdir -p "$path"
    (cd "$path" && pwd)
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

build_root_macro_call() {
    local macro="$1"
    local src="$2"
    local bg="$3"
    local gain_value="$4"
    local bg_scale="$5"
    local quantile="$6"
    local bins="$7"
    local xmin="$8"
    local xmax="$9"
    local clock="${10}"
    local prefix="${11}"
    local dynamic_range="${12}"
    local sampling_ns="${13}"
    local resistance="${14}"
    local adc_bits_value="${15}"
    local label="${16}"

    printf '%s(%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)' \
        "$macro" \
        "$(root_string_arg "$src")" \
        "$(root_string_arg "$bg")" \
        "$gain_value" \
        "$bg_scale" \
        "$quantile" \
        "$bins" \
        "$xmin" \
        "$xmax" \
        "$clock" \
        "$(root_string_arg "$prefix")" \
        "$dynamic_range" \
        "$sampling_ns" \
        "$resistance" \
        "$adc_bits_value" \
        "$(root_string_arg "$label")"
}

source_file=""
bg_file=""
data_dir="$START_DIR"
out_prefix=""
output_dir="$REPO_DIR/results"
macro_file="$SCRIPT_DIR/plot_npe_subtracted.C"
gain="1.0e7"
bg_scale="-1"
x_quantile="1.0"
n_bins="400"
x_min="0"
x_max="-1"
ttt_clock="125.0e6"
source_label=""
dynamic_range_v="2.0"
sampling_time_ns="2.0"
resistance_ohm="50.0"
adc_bits="14"

while getopts ":s:b:d:o:O:L:m:g:q:n:x:X:c:h" opt; do
    case "$opt" in
        s) source_file="$OPTARG" ;;
        b) bg_file="$OPTARG" ;;
        d) data_dir="$OPTARG" ;;
        o) out_prefix="$OPTARG" ;;
        O) output_dir="$OPTARG" ;;
        L) source_label="$OPTARG" ;;
        m) macro_file="$OPTARG" ;;
        g) gain="$OPTARG" ;;
        q) x_quantile="$OPTARG" ;;
        n) n_bins="$OPTARG" ;;
        x) x_min="$OPTARG" ;;
        X) x_max="$OPTARG" ;;
        c) ttt_clock="$OPTARG" ;;
        h) usage; exit 0 ;;
        :) die "Option -$OPTARG requires an argument." ;;
        \?) die "Unknown option: -$OPTARG" ;;
    esac
done

macro_file="$(abs_path "$macro_file")"
data_dir="$(abs_dir "$data_dir")"
output_dir="$(abs_dir "$output_dir")"

[[ -f "$macro_file" ]] || die "ROOT macro not found: $macro_file"
command -v root >/dev/null 2>&1 || die "ROOT command 'root' not found in PATH"

if [[ -z "$source_file" || -z "$bg_file" ]]; then
    mapfile -t root_files < <(find "$data_dir" -maxdepth 1 -type f -name '*.root' -printf '%f\n' | sort)
    (( ${#root_files[@]} > 0 )) || die "No ROOT files found in $data_dir"

    [[ -n "$source_file" ]] || source_file="$(pick_root_file "source" "${root_files[@]}")"
    [[ -n "$bg_file" ]] || bg_file="$(pick_root_file "background" "${root_files[@]}")"
fi

if [[ "$source_file" != /* ]]; then
    if [[ -f "$data_dir/$source_file" ]]; then
        source_file="$data_dir/$source_file"
    else
        source_file="$(abs_path "$source_file")"
    fi
fi
if [[ "$bg_file" != /* ]]; then
    if [[ -f "$data_dir/$bg_file" ]]; then
        bg_file="$data_dir/$bg_file"
    else
        bg_file="$(abs_path "$bg_file")"
    fi
fi

[[ -f "$source_file" ]] || die "Source file not found: $source_file"
[[ -f "$bg_file" ]] || die "Background file not found: $bg_file"

if [[ -z "$out_prefix" ]]; then
    base="$(basename "${source_file%.root}")"
    base="${base//[^A-Za-z0-9_.-]/_}"
    out_prefix="$(prompt_default "Output prefix" "$base")"
fi

if [[ -t 0 && "$#" -eq 0 ]]; then
    gain="$(prompt_default "PMT gain" "$gain")"
    x_quantile="$(prompt_default "X quantile, 1.0 = full range" "$x_quantile")"
    n_bins="$(prompt_default "Histogram bins" "$n_bins")"
    echo "Manual x range: press Enter twice to use automatic/full range."
    x_min="$(prompt_default "X min" "$x_min")"
    x_max="$(prompt_default "X max (-1 = automatic)" "$x_max")"
fi

echo
echo "Running NPE analysis"
echo "  Macro:      $macro_file"
echo "  Data dir:   $data_dir"
echo "  Output dir: $output_dir"
echo "  Source:     $source_file"
echo "  Background: $bg_file"
echo "  Prefix:     $out_prefix"
echo "  Label:      ${source_label:-source file name}"
echo "  Gain:       $gain"
echo "  X quantile: $x_quantile"
echo "  Bins:       $n_bins"
echo "  X range:    [$x_min, $x_max]"
echo

cd "$output_dir"
macro_call="$(build_root_macro_call \
    "$macro_file" \
    "$source_file" \
    "$bg_file" \
    "$gain" \
    "$bg_scale" \
    "$x_quantile" \
    "$n_bins" \
    "$x_min" \
    "$x_max" \
    "$ttt_clock" \
    "$out_prefix" \
    "$dynamic_range_v" \
    "$sampling_time_ns" \
    "$resistance_ohm" \
    "$adc_bits" \
    "$source_label")"
root -l -b -q "$macro_call"

echo
echo "Key outputs:"
echo "  $output_dir/${out_prefix}_overlay_total_rate_log.png"
echo "  $output_dir/${out_prefix}_overlay_total_log.png"
echo "  $output_dir/${out_prefix}_subtracted_total.png"
