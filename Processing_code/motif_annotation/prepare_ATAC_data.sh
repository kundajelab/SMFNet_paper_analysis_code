#!/usr/bin/env bash
#Run from SMFNet_paper_analysis_code github repo: bash ./Processing_code/motif_annotation/prepare_ATAC_data.sh
set -euo pipefail
INPUT=input_files/motif_annotation
OUTPUT=output/motif_annotation
WORK="$OUTPUT/atac_preprocessing"
THREADS=8
TRIM_CORES=8
CALL_PEAKS=1
mkdir -p "$WORK"/{sra,fastq,trimmed,bam,logs,tmp,reference,peaks}

for tool in prefetch fasterq-dump trim_galore bowtie2 bowtie2-build samtools picard bedtools; do
    command -v "$tool" >/dev/null || { echo "Required command not found: $tool" >&2; exit 1; }
done
if [[ "$CALL_PEAKS" == 1 ]]; then command -v macs2 >/dev/null; fi
for name in mm10.fa mm10-blacklist.v2.bed; do
    test -s "$INPUT/$name" || { echo "Missing input: $INPUT/$name" >&2; exit 1; }
done

INDEX="$WORK/reference/mm10"
if [[ ! -s "$INDEX.1.bt2" && ! -s "$INDEX.1.bt2l" ]]; then
    bowtie2-build --threads "$THREADS" "$INPUT/mm10.fa" "$INDEX"
fi

#ESC-TKO biological replicates from PMID36471082.
RUNS=(SRR15974966 SRR15974967)
BAMS=()
for run in "${RUNS[@]}"; do
    mkdir -p "$WORK/tmp/$run"
    prefetch --max-size 100G --output-directory "$WORK/sra" "$run"
    fasterq-dump --split-files --threads "$THREADS" --temp "$WORK/tmp/$run" \
        --outdir "$WORK/fastq" "$WORK/sra/$run/$run.sra"
    trim_galore --cores "$TRIM_CORES" --paired --nextera \
        "$WORK/fastq/${run}_1.fastq" "$WORK/fastq/${run}_2.fastq" --output_dir "$WORK/trimmed"

    bowtie2 -p "$THREADS" --local --very-sensitive-local -X 2000 -x "$INDEX" \
        -1 "$WORK/trimmed/${run}_1_val_1.fq" -2 "$WORK/trimmed/${run}_2_val_2.fq" \
        2> "$WORK/logs/${run}_alignment.log" \
      | samtools view -h -f 2 -q 20 - \
      | awk '/^@/ || ($3 != "chrM" && $7 != "chrM")' \
      | samtools sort -@ "$THREADS" -m 1G -T "$WORK/tmp/$run/sort" -o "$WORK/bam/${run}.sorted.bam" -

    picard MarkDuplicates INPUT="$WORK/bam/${run}.sorted.bam" \
        OUTPUT="$WORK/bam/${run}.dedup.bam" METRICS_FILE="$WORK/logs/${run}_duplicates.txt" \
        VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=true
    bedtools intersect -v -a "$WORK/bam/${run}.dedup.bam" -b "$INPUT/mm10-blacklist.v2.bed" \
        > "$WORK/bam/${run}.final.bam"
    samtools index -@ "$THREADS" "$WORK/bam/${run}.final.bam"
    BAMS+=("$WORK/bam/${run}.final.bam")

    if [[ "$CALL_PEAKS" == 1 ]]; then
        macs2 callpeak -t "$WORK/bam/${run}.final.bam" -f BAMPE -g mm -q 0.05 \
            -n "$run" --outdir "$WORK/peaks"
    fi
done

#Keep the full chromosome set for ChromBPNet.
samtools merge -@ "$THREADS" -o "$WORK/mESC_ATAC_merged.bam" "${BAMS[@]}"
samtools index -@ "$THREADS" "$WORK/mESC_ATAC_merged.bam"
samtools quickcheck -v "$WORK/mESC_ATAC_merged.bam"
if [[ "$CALL_PEAKS" == 1 ]]; then
    bedtools intersect -a "$WORK/peaks/${RUNS[0]}_peaks.narrowPeak" \
        -b "$WORK/peaks/${RUNS[1]}_peaks.narrowPeak" > "$WORK/peaks/mESC_ATAC_replicate_overlap.bed"
fi
printf 'Training BAM: %s\n' "$WORK/mESC_ATAC_merged.bam"
