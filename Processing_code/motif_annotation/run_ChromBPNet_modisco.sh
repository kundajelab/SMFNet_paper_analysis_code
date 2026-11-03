#!/usr/bin/env bash
#Run from SMFNet_paper_analysis_code github repo: bash ./Processing_code/motif_annotation/run_ChromBPNet_modisco.sh
set -euo pipefail
INPUT=input_files/motif_annotation
OUTPUT=output/motif_annotation
mkdir -p "$OUTPUT/modisco/report"

modisco motifs -s "$OUTPUT/shap/ohe.npz" -a "$OUTPUT/shap/attr_hypo.npz" \
    -n 50000 -o "$OUTPUT/modisco/modisco_results.h5"
modisco report -i "$OUTPUT/modisco/modisco_results.h5" -o "$OUTPUT/modisco/report" \
    -m "$INPUT/JASPAR2026_vertebrates.meme" -t
modisco meme -i "$OUTPUT/modisco/modisco_results.h5" -o "$OUTPUT/modisco/chrombpnet_motifs_cwm.meme" -t CWM
modisco meme -i "$OUTPUT/modisco/modisco_results.h5" -o "$OUTPUT/modisco/chrombpnet_motifs_pfm.meme" -t PFM
python Processing_code/motif_annotation/name_motifs.py
