#!/usr/bin/env bash
#Run from SMFNet_paper_analysis_code github repo: bash ./Processing_code/motif_annotation/train_ChromBPNet.sh
set -euo pipefail
INPUT=input_files/motif_annotation
OUTPUT=output/motif_annotation
MODEL_FILE=trained_models/model_ChromBPNet_mESC_nobias.h5
ATAC_BAM="$OUTPUT/atac_preprocessing/mESC_ATAC_merged.bam"
mkdir -p "$OUTPUT" trained_models
test -s "$ATAC_BAM" || { echo "Missing training BAM: $ATAC_BAM; run prepare_ATAC_data.sh first" >&2; exit 1; }
for name in mm10.fa mm10.chrom.sizes fold0.json mESC_accessible_regions.narrowPeak mESC_background_regions.narrowPeak mESC_bias_nonpeaks_fold0.bed; do
    test -s "$INPUT/$name" || { echo "Missing input: $INPUT/$name" >&2; exit 1; }
done

#Train the bias model and the sequence model with the same chromosome split.
chrombpnet bias pipeline \
    -ibam "$ATAC_BAM" -d ATAC \
    -g "$INPUT/mm10.fa" -c "$INPUT/mm10.chrom.sizes" \
    -p "$INPUT/mESC_accessible_regions.narrowPeak" -n "$INPUT/mESC_bias_nonpeaks_fold0.bed" \
    -fl "$INPUT/fold0.json" -b 0.2 -o "$OUTPUT/bias_fold0" \
    -il 2114 -ol 1000 -s 1234 -e 50 -es 5 -l 0.001 -bs 64 -fil 128 -dil 4 -j 0

chrombpnet pipeline \
    -ibam "$ATAC_BAM" -d ATAC \
    -g "$INPUT/mm10.fa" -c "$INPUT/mm10.chrom.sizes" \
    -p "$INPUT/mESC_accessible_regions.narrowPeak" -n "$INPUT/mESC_background_regions.narrowPeak" \
    -fl "$INPUT/fold0.json" -b "$OUTPUT/bias_fold0/models/bias.h5" -o "$OUTPUT/fold0" \
    -il 2114 -ol 1000 -s 1234 -e 50 -es 5 -l 0.001 -bs 64 -fil 512 -dil 8 -j 500 -sr 0.1

#Store the model used for motif annotation under the Zenodo filename.
cp "$OUTPUT/fold0/models/chrombpnet_nobias.h5" "$MODEL_FILE"
