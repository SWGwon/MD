#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_ROOT="${1:-/home/sgwon/MD/data/optical_test}"
RUN_TAG="$(date +%Y%m%d_%H%M%S)"
OUT_DIR="${2:-$SCRIPT_DIR/results/optical_gel_comparison_repro_$RUN_TAG}"

WITHOUT_DIR="$DATA_ROOT/without_optical"
WITH_DIR="$DATA_ROOT/with_optical_gel"

for file in \
    "$WITHOUT_DIR/Co60_1hr_prod.root" \
    "$WITHOUT_DIR/Cs137_1hr_prod.root" \
    "$WITHOUT_DIR/background_1hr_prod.root" \
    "$WITH_DIR/Co60_1hr_prod.root" \
    "$WITH_DIR/Cs137_1hr_prod.root" \
    "$WITH_DIR/background_1hr_prod.root"; do
    if [[ ! -f "$file" ]]; then
        echo "Error: required input file not found: $file" >&2
        exit 1
    fi
done

if [[ ! -f "$SCRIPT_DIR/compare_optical_gel_histograms.C" ]]; then
    echo "Error: optical-gel comparison macro not found: $SCRIPT_DIR/compare_optical_gel_histograms.C" >&2
    exit 1
fi
if [[ ! -f "$SCRIPT_DIR/compare_source_edge_linearity.C" ]]; then
    echo "Error: source-linearity comparison macro not found: $SCRIPT_DIR/compare_source_edge_linearity.C" >&2
    exit 1
fi

if [[ -d "$OUT_DIR" ]] && find "$OUT_DIR" -mindepth 1 -maxdepth 1 -print -quit | grep -q .; then
    echo "Error: output directory is not empty: $OUT_DIR" >&2
    echo "Choose a new directory so stale files cannot be mistaken for fresh results." >&2
    exit 1
fi
mkdir -p "$OUT_DIR"

"$SCRIPT_DIR/run_compton_edge_analysis.sh" \
    -d "$WITHOUT_DIR" \
    -s Co60_1hr_prod.root \
    -b background_1hr_prod.root \
    -O "$OUT_DIR" \
    -o Co60_without_gel_wide \
    -L "Co60 without optical gel" \
    -g 1.0e7 -B -1 -q 1.0 -c 125.0e6 \
    -x 0 -X 1200 -n 1200 -C 0,1

"$SCRIPT_DIR/run_compton_edge_analysis.sh" \
    -d "$WITH_DIR" \
    -s Co60_1hr_prod.root \
    -b background_1hr_prod.root \
    -O "$OUT_DIR" \
    -o Co60_with_gel_wide \
    -L "Co60 with optical gel" \
    -g 1.0e7 -B -1 -q 1.0 -c 125.0e6 \
    -x 0 -X 1200 -n 1200 -C 0,1

CO_WITHOUT_HIST="$OUT_DIR/Co60_without_gel_wide_histograms.root"
CO_WITH_HIST="$OUT_DIR/Co60_with_gel_wide_histograms.root"

# Condition-tailored windows used for the primary effective-edge comparison.
"$SCRIPT_DIR/run_compton_edge_fit.sh" \
    -i "$CO_WITHOUT_HIST" -x 100 -X 400 \
    -o "$OUT_DIR/Co60_without_gel_total_100_400" \
    -M erfc_gaussian -S Co60

"$SCRIPT_DIR/run_compton_edge_fit.sh" \
    -i "$CO_WITH_HIST" -x 250 -X 1000 \
    -o "$OUT_DIR/Co60_with_gel_total_250_1000" \
    -M erfc_gaussian -S Co60

"$SCRIPT_DIR/run_compton_edge_scan.sh" \
    -i "$CO_WITHOUT_HIST" -x 100 -X 400 \
    -o "$OUT_DIR/Co60_without_gel_total_100_400" \
    -M erfc_gaussian -S Co60 -F 0.10

"$SCRIPT_DIR/run_compton_edge_scan.sh" \
    -i "$CO_WITH_HIST" -x 250 -X 1000 \
    -o "$OUT_DIR/Co60_with_gel_total_250_1000" \
    -M erfc_gaussian -S Co60 -F 0.10

# A window scaled from the no-gel range is retained as a fit-window stress test.
"$SCRIPT_DIR/run_compton_edge_fit.sh" \
    -i "$CO_WITH_HIST" -x 214 -X 855 \
    -o "$OUT_DIR/Co60_with_gel_total_214_855" \
    -M erfc_gaussian -S Co60

"$SCRIPT_DIR/run_compton_edge_scan.sh" \
    -i "$CO_WITH_HIST" -x 214 -X 855 \
    -o "$OUT_DIR/Co60_with_gel_total_214_855" \
    -M erfc_gaussian -S Co60 -F 0.10

# Binning and empirical-model stress tests cited in the report.
"$SCRIPT_DIR/run_compton_edge_analysis.sh" \
    -d "$WITHOUT_DIR" \
    -s Co60_1hr_prod.root \
    -b background_1hr_prod.root \
    -O "$OUT_DIR" \
    -o Co60_without_gel_wide_2npe \
    -L "Co60 without optical gel" \
    -g 1.0e7 -B -1 -q 1.0 -c 125.0e6 \
    -x 0 -X 1200 -n 600 -C 0,1

"$SCRIPT_DIR/run_compton_edge_analysis.sh" \
    -d "$WITH_DIR" \
    -s Co60_1hr_prod.root \
    -b background_1hr_prod.root \
    -O "$OUT_DIR" \
    -o Co60_with_gel_wide_2npe \
    -L "Co60 with optical gel" \
    -g 1.0e7 -B -1 -q 1.0 -c 125.0e6 \
    -x 0 -X 1200 -n 600 -C 0,1

"$SCRIPT_DIR/run_compton_edge_fit.sh" \
    -i "$OUT_DIR/Co60_without_gel_wide_2npe_histograms.root" \
    -x 100 -X 400 \
    -o "$OUT_DIR/Co60_without_gel_2npe_total_100_400" \
    -M erfc_gaussian -S Co60

"$SCRIPT_DIR/run_compton_edge_fit.sh" \
    -i "$OUT_DIR/Co60_with_gel_wide_2npe_histograms.root" \
    -x 250 -X 1000 \
    -o "$OUT_DIR/Co60_with_gel_2npe_total_250_1000" \
    -M erfc_gaussian -S Co60

"$SCRIPT_DIR/run_compton_edge_model_compare.sh" \
    -i "$CO_WITHOUT_HIST" -x 100 -X 400 \
    -o "$OUT_DIR/Co60_without_gel_total_100_400_model"

"$SCRIPT_DIR/run_compton_edge_model_compare.sh" \
    -i "$CO_WITH_HIST" -x 250 -X 1000 \
    -o "$OUT_DIR/Co60_with_gel_total_250_1000_model"

root_escape() {
    local value="$1"
    value="${value//\\/\\\\}"
    value="${value//\"/\\\"}"
    printf '"%s"' "$value"
}

tree_ttt_span_seconds() {
    local file="$1"
    local root_expression
    local span
    root_expression="$(printf 'TFile f(%s); auto *t=(TTree*)f.Get("phys_tree"); if (!t) { gSystem->Exit(2); } printf("TTT_SPAN_SECONDS %%.12f\\n",(t->GetMaximum("SyncTime_TTT")-t->GetMinimum("SyncTime_TTT"))/125.0e6);' \
        "$(root_escape "$file")")"
    span="$(root -l -b -q -e "$root_expression" | awk '$1 == "TTT_SPAN_SECONDS" { print $2; exit }')"
    if [[ -z "$span" ]]; then
        echo "Error: could not determine SyncTime_TTT span for $file" >&2
        exit 1
    fi
    printf '%s' "$span"
}

CO_WITHOUT_TTT_SPAN_SECONDS="$(tree_ttt_span_seconds "$WITHOUT_DIR/Co60_1hr_prod.root")"
CO_WITH_TTT_SPAN_SECONDS="$(tree_ttt_span_seconds "$WITH_DIR/Co60_1hr_prod.root")"
CS_WITHOUT_TTT_SPAN_SECONDS="$(tree_ttt_span_seconds "$WITHOUT_DIR/Cs137_1hr_prod.root")"
CS_WITH_TTT_SPAN_SECONDS="$(tree_ttt_span_seconds "$WITH_DIR/Cs137_1hr_prod.root")"

compare_call="$(printf '%s(%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)' \
    "$SCRIPT_DIR/compare_optical_gel_histograms.C" \
    "$(root_escape "$CO_WITHOUT_HIST")" \
    "$(root_escape "$CO_WITH_HIST")" \
    "$(root_escape "without gel")" \
    "$(root_escape "with gel")" \
    "$(root_escape "$OUT_DIR/Co60_gel_comparison")" \
    "$CO_WITHOUT_TTT_SPAN_SECONDS" \
    "$CO_WITH_TTT_SPAN_SECONDS" \
    50 \
    1000 \
    50 \
    600)"
root -l -b -q "$compare_call"

"$SCRIPT_DIR/run_compton_edge_analysis.sh" \
    -d "$WITHOUT_DIR" \
    -s Cs137_1hr_prod.root \
    -b background_1hr_prod.root \
    -O "$OUT_DIR" \
    -o Cs137_without_gel_wide \
    -L "Cs137 without optical gel" \
    -g 1.0e7 -B -1 -q 1.0 -c 125.0e6 \
    -x 0 -X 1200 -n 1200 -C 0,1

"$SCRIPT_DIR/run_compton_edge_analysis.sh" \
    -d "$WITH_DIR" \
    -s Cs137_1hr_prod.root \
    -b background_1hr_prod.root \
    -O "$OUT_DIR" \
    -o Cs137_with_gel_wide \
    -L "Cs137 with optical gel" \
    -g 1.0e7 -B -1 -q 1.0 -c 125.0e6 \
    -x 0 -X 1200 -n 1200 -C 0,1

CS_WITHOUT_HIST="$OUT_DIR/Cs137_without_gel_wide_histograms.root"
CS_WITH_HIST="$OUT_DIR/Cs137_with_gel_wide_histograms.root"

"$SCRIPT_DIR/run_compton_edge_fit.sh" \
    -i "$CS_WITHOUT_HIST" -x 40 -X 180 \
    -o "$OUT_DIR/Cs137_without_gel_total_40_180" \
    -M erfc_gaussian -S Cs137

"$SCRIPT_DIR/run_compton_edge_fit.sh" \
    -i "$CS_WITH_HIST" -x 85 -X 385 \
    -o "$OUT_DIR/Cs137_with_gel_total_85_385" \
    -M erfc_gaussian -S Cs137

"$SCRIPT_DIR/run_compton_edge_scan.sh" \
    -i "$CS_WITHOUT_HIST" -x 40 -X 180 \
    -o "$OUT_DIR/Cs137_without_gel_total_40_180" \
    -M erfc_gaussian -S Cs137 -F 0.10

"$SCRIPT_DIR/run_compton_edge_scan.sh" \
    -i "$CS_WITH_HIST" -x 85 -X 385 \
    -o "$OUT_DIR/Cs137_with_gel_total_85_385" \
    -M erfc_gaussian -S Cs137 -F 0.10

"$SCRIPT_DIR/run_compton_edge_analysis.sh" \
    -d "$WITHOUT_DIR" \
    -s Cs137_1hr_prod.root \
    -b background_1hr_prod.root \
    -O "$OUT_DIR" \
    -o Cs137_without_gel_wide_2npe \
    -L "Cs137 without optical gel" \
    -g 1.0e7 -B -1 -q 1.0 -c 125.0e6 \
    -x 0 -X 1200 -n 600 -C 0,1

"$SCRIPT_DIR/run_compton_edge_analysis.sh" \
    -d "$WITH_DIR" \
    -s Cs137_1hr_prod.root \
    -b background_1hr_prod.root \
    -O "$OUT_DIR" \
    -o Cs137_with_gel_wide_2npe \
    -L "Cs137 with optical gel" \
    -g 1.0e7 -B -1 -q 1.0 -c 125.0e6 \
    -x 0 -X 1200 -n 600 -C 0,1

"$SCRIPT_DIR/run_compton_edge_fit.sh" \
    -i "$OUT_DIR/Cs137_without_gel_wide_2npe_histograms.root" \
    -x 40 -X 180 \
    -o "$OUT_DIR/Cs137_without_gel_2npe_total_40_180" \
    -M erfc_gaussian -S Cs137

"$SCRIPT_DIR/run_compton_edge_fit.sh" \
    -i "$OUT_DIR/Cs137_with_gel_wide_2npe_histograms.root" \
    -x 85 -X 385 \
    -o "$OUT_DIR/Cs137_with_gel_2npe_total_85_385" \
    -M erfc_gaussian -S Cs137

"$SCRIPT_DIR/run_compton_edge_model_compare.sh" \
    -i "$CS_WITHOUT_HIST" -x 40 -X 180 \
    -o "$OUT_DIR/Cs137_without_gel_total_40_180_model"

"$SCRIPT_DIR/run_compton_edge_model_compare.sh" \
    -i "$CS_WITH_HIST" -x 85 -X 385 \
    -o "$OUT_DIR/Cs137_with_gel_total_85_385_model"

cs_compare_call="$(printf '%s(%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)' \
    "$SCRIPT_DIR/compare_optical_gel_histograms.C" \
    "$(root_escape "$CS_WITHOUT_HIST")" \
    "$(root_escape "$CS_WITH_HIST")" \
    "$(root_escape "without gel")" \
    "$(root_escape "with gel")" \
    "$(root_escape "$OUT_DIR/Cs137_gel_comparison")" \
    "$CS_WITHOUT_TTT_SPAN_SECONDS" \
    "$CS_WITH_TTT_SPAN_SECONDS" \
    30 \
    500 \
    30 \
    350 \
    "$(root_escape "Cs-137")")"
root -l -b -q "$cs_compare_call"

source_compare_call="$(printf '%s(%s,%s,%s,%s,%s)' \
    "$SCRIPT_DIR/compare_source_edge_linearity.C" \
    "$(root_escape "$OUT_DIR/Cs137_without_gel_total_40_180_compton_edge_fit.txt")" \
    "$(root_escape "$OUT_DIR/Cs137_with_gel_total_85_385_compton_edge_fit.txt")" \
    "$(root_escape "$OUT_DIR/Co60_without_gel_total_100_400_compton_edge_fit.txt")" \
    "$(root_escape "$OUT_DIR/Co60_with_gel_total_250_1000_compton_edge_fit.txt")" \
    "$(root_escape "$OUT_DIR/Cs137_Co60_edge_linearity")")"
root -l -b -q "$source_compare_call"

echo
echo "Optical-gel comparison complete."
echo "Results: $OUT_DIR"
