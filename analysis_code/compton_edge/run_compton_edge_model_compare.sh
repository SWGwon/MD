#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
START_DIR="$(pwd)"

usage() {
    cat <<'EOF'
Usage:
  ./run_compton_edge_model_compare.sh -i HISTOGRAMS.root -x MIN -X MAX -o PREFIX [options]

Options:
  -i FILE      Histogram ROOT file from run_compton_edge_analysis.sh
  -x MIN       Fit x-axis minimum in NPE
  -X MAX       Fit x-axis maximum in NPE
  -o PREFIX    Output prefix for model-comparison outputs
  -f FILE      Fit wrapper path (default: run_compton_edge_fit.sh next to this script)
  -H NAME      Histogram name (default: hTotalSub, falling back to hTotalSourceOnly)
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

get_result_value() {
    local line="$1"
    local key="$2"
    for item in $line; do
        if [[ "$item" == "$key="* ]]; then
            printf "%s" "${item#*=}"
            return 0
        fi
    done
    printf ""
}

hist_file=""
fit_xmin=""
fit_xmax=""
out_prefix=""
fit_wrapper="$SCRIPT_DIR/run_compton_edge_fit.sh"
hist_name=""

while getopts ":i:x:X:o:f:H:h" opt; do
    case "$opt" in
        i) hist_file="$OPTARG" ;;
        x) fit_xmin="$OPTARG" ;;
        X) fit_xmax="$OPTARG" ;;
        o) out_prefix="$OPTARG" ;;
        f) fit_wrapper="$OPTARG" ;;
        H) hist_name="$OPTARG" ;;
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
fit_wrapper="$(abs_path "$fit_wrapper")"

[[ -f "$hist_file" ]] || die "Histogram ROOT file not found: $hist_file"
[[ -x "$fit_wrapper" ]] || die "Fit wrapper not found or not executable: $fit_wrapper"

echo
echo "Running Compton edge model comparison"
echo "  Fit wrapper: $fit_wrapper"
echo "  Hist file:   $hist_file"
echo "  Hist name:   ${hist_name:-auto}"
echo "  Fit range:   [$fit_xmin, $fit_xmax]"
echo "  Prefix:      $out_prefix"
echo

fit_args=(-i "$hist_file" -x "$fit_xmin" -X "$fit_xmax")
if [[ -n "$hist_name" ]]; then
    fit_args+=(-H "$hist_name")
fi

linear_output="$("$fit_wrapper" "${fit_args[@]}" -o "${out_prefix}_erfc_linear" -M erfc_linear)"
printf "%s\n" "$linear_output"
gaussian_output="$("$fit_wrapper" "${fit_args[@]}" -o "${out_prefix}_erfc_gaussian" -M erfc_gaussian)"
printf "%s\n" "$gaussian_output"

linear_line="$(printf "%s\n" "$linear_output" | sed -n '/^FIT_RESULT /p' | tail -n 1)"
gaussian_line="$(printf "%s\n" "$gaussian_output" | sed -n '/^FIT_RESULT /p' | tail -n 1)"

[[ -n "$linear_line" ]] || die "Could not parse erfc_linear FIT_RESULT."
[[ -n "$gaussian_line" ]] || die "Could not parse erfc_gaussian FIT_RESULT."

linear_edge="$(get_result_value "$linear_line" edge)"
linear_error="$(get_result_value "$linear_line" edge_error)"
linear_chi2="$(get_result_value "$linear_line" chi2_ndf)"
linear_warnings="$(get_result_value "$linear_line" warnings)"
gaussian_edge="$(get_result_value "$gaussian_line" edge)"
gaussian_error="$(get_result_value "$gaussian_line" edge_error)"
gaussian_chi2="$(get_result_value "$gaussian_line" chi2_ndf)"
gaussian_warnings="$(get_result_value "$gaussian_line" warnings)"

read -r edge_diff edge_diff_abs combined_error edge_pull <<EOF
$(awk -v a="$linear_edge" -v b="$gaussian_edge" -v ea="$linear_error" -v eb="$gaussian_error" 'BEGIN {
    diff = b - a;
    absdiff = diff < 0 ? -diff : diff;
    err = sqrt(ea * ea + eb * eb);
    pull = err > 0 ? diff / err : 0;
    printf "%.12g %.12g %.12g %.12g", diff, absdiff, err, pull;
}')
EOF

summary_file="${out_prefix}_compton_edge_model_compare.txt"
{
    echo "hist_file $hist_file"
    echo "hist_name ${hist_name:-auto}"
    echo "fit_xmin $fit_xmin"
    echo "fit_xmax $fit_xmax"
    echo "linear_edge $linear_edge"
    echo "linear_edge_error $linear_error"
    echo "linear_chi2_ndf $linear_chi2"
    echo "linear_warnings ${linear_warnings:-none}"
    echo "gaussian_edge $gaussian_edge"
    echo "gaussian_edge_error $gaussian_error"
    echo "gaussian_chi2_ndf $gaussian_chi2"
    echo "gaussian_warnings ${gaussian_warnings:-none}"
    echo "edge_diff_gaussian_minus_linear $edge_diff"
    echo "edge_diff_abs $edge_diff_abs"
    echo "combined_stat_error $combined_error"
    echo "edge_diff_pull $edge_pull"
    echo "linear_png ${out_prefix}_erfc_linear_compton_edge_fit.png"
    echo "gaussian_png ${out_prefix}_erfc_gaussian_compton_edge_fit.png"
} > "$summary_file"

echo
echo "--- Compton Edge Model Comparison Complete ---"
echo "erfc_linear edge [NPE]: $linear_edge +/- $linear_error"
echo "erfc_gaussian edge [NPE]: $gaussian_edge +/- $gaussian_error"
echo "Gaussian - linear edge difference [NPE]: $edge_diff"
echo "Absolute model difference [NPE]: $edge_diff_abs"
echo "Combined statistical error [NPE]: $combined_error"
echo "Difference pull: $edge_pull"
echo "Model comparison summary: $summary_file"
echo "MODEL_COMPARE_RESULT linear_edge=$linear_edge linear_error=$linear_error gaussian_edge=$gaussian_edge gaussian_error=$gaussian_error edge_diff=$edge_diff edge_diff_abs=$edge_diff_abs combined_error=$combined_error edge_diff_pull=$edge_pull linear_chi2_ndf=$linear_chi2 gaussian_chi2_ndf=$gaussian_chi2 linear_warnings=${linear_warnings:-none} gaussian_warnings=${gaussian_warnings:-none}"
