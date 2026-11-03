#!/usr/bin/env bash
#Run from SMFNet_paper_analysis_code github repo: bash ./Processing_code/motif_annotation/run_ChromBPNet_finemo.sh
set -euo pipefail
OUTPUT=output/motif_annotation
mkdir -p "$OUTPUT/finemo"

finemo extract-regions-modisco-fmt -s "$OUTPUT/shap/ohe.npz" -a "$OUTPUT/shap/attr_hypo.npz" \
    -o "$OUTPUT/finemo/finemo_in" -p "$OUTPUT/shap/region_ids.narrowPeak" -w 2114
finemo call-hits -r "$OUTPUT/finemo/finemo_in.npz" -m "$OUTPUT/modisco/modisco_results.h5" -o "$OUTPUT/finemo"
finemo report -r "$OUTPUT/finemo/finemo_in.npz" -H "$OUTPUT/finemo" \
    -o "$OUTPUT/finemo/report" -m "$OUTPUT/modisco/modisco_results.h5"
