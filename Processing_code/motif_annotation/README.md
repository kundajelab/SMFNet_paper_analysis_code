# ChromBPNet motif annotation

## Overview

Train an mESC ATAC ChromBPNet model with the fold0 chromosome split, discover motifs with TF-MoDISco, and locate motif instances with Fi-NeMo.

- `prepare_ATAC_data.sh`: Example workflow to download and preprocess the published ESC-TKO ATAC replicates.
- `train_ChromBPNet.sh`: Trains a fold0 bias model, then the fold0 ChromBPNet model.
- `compute_ChromBPNet_DeepSHAP.py`: Computes count-head contributions from the model without the bias branch.
- `run_ChromBPNet_modisco.sh`: Discovers motifs, compares them with JASPAR, and exports CWM and PFM files.
- `name_motifs.py`: Creates an editable table of proposed TF names from the best Tomtom matches.
- `trim_motif_CWMs.py`: Trims motif flanks and writes individual motif logos.
- `run_ChromBPNet_finemo.sh`: Calls motif instances using the contributions and discovered motifs.
- `annotate_finemo_hits.py`: Adds TF names and the corresponding 2048 bp region IDs to the hits.

## Instructions

Run all commands from the repository root. Paths and settings are defined at the top of each script; the example scripts take no command-line arguments.

### Input files supplied with the repository

Keep these seven small files in `input_files/motif_annotation/`:

| File | Purpose |
|------|---------|
| `fold0.json` | Training, validation and test chromosome split |
| `mESC_accessible_regions.narrowPeak` | Accessible bait regions for ChromBPNet training |
| `mESC_background_regions.narrowPeak` | Background bait regions for ChromBPNet training |
| `mESC_bias_nonpeaks_fold0.bed` | Nonpeak regions for fold0 bias-model training |
| `mESC_annotation_regions.bed` | Regions for contributions and motif annotation |
| `mESC_annotation_regions_2048bp.bed` | Region labels for the final motif hits |
| `JASPAR2026_vertebrates.meme` | Reference motifs for TF identification |

Together these files are about 22 MiB uncompressed and can be included directly in GitHub. Background regions used for the sequence model and nonpeak regions used for bias training are separate inputs.

The full workflow uses these files plus public downloads: SRA runs `SRR15974966` and `SRR15974967`, the UCSC mm10 reference, and the mm10 v2 blacklist. Download commands are provided below. It generates the merged BAM, bias model, ChromBPNet model, contribution arrays, MoDISco results and Fi-NeMo hits. No additional Zenodo-only input is required to run from the beginning.

### Precomputed files available on Zenodo

These archives are already listed in the Zenodo upload inventory:

| Zenodo archive | Relevant file inside | Use |
|----------------|----------------------|-----|
| `model_ChromBPNet_mESC.tar.gz` | `model_ChromBPNet_mESC_nobias.h5` | Skip training and compute contributions with the supplied model |
| `mESC_ChromBPnet_annotated_hits.tar.gz` | `SMFNet_annotated_motif_hits.tsv` | Use the completed ChromBPNet/Fi-NeMo annotation directly |
| `F1_preprocessing_inputs.tar.gz` | `motif_annotation/modisco_results.h5` and `motif_annotation/chrombpnet_motifs_cwm_trim.meme` | Original discovered motifs and reviewed motif names |

The file `SMFNet_annotated_motif_hits.tsv` retains its existing archive filename; it contains the ChromBPNet-derived Fi-NeMo hits with TF names and 2048 bp region labels. It is the precomputed counterpart of this example's `output/motif_annotation/finemo/annotated_hits.tsv`.

To use the trained model and completed annotation, extract their archives from the repository root:

```bash
mkdir -p trained_models input_files/motif_annotation
tar -xzf model_ChromBPNet_mESC.tar.gz -C trained_models
tar -xzf mESC_ChromBPnet_annotated_hits.tar.gz -C input_files/motif_annotation
```

The resulting paths are `trained_models/model_ChromBPNet_mESC_nobias.h5` and `input_files/motif_annotation/SMFNet_annotated_motif_hits.tsv`. The completed annotation can be used without rerunning training, attribution or motif discovery.

To extract only the original motif files from the shared F1 archive:

```bash
tar -xzf F1_preprocessing_inputs.tar.gz -C input_files \
  motif_annotation/modisco_results.h5 \
  motif_annotation/chrombpnet_motifs_cwm_trim.meme
```

These optional files preserve the original motif patterns and reviewed labels. Running Fi-NeMo again still requires contribution arrays, which the DeepSHAP script generates. Newly trained models and newly discovered patterns can differ from the supplied results; use the archived annotated table when the original annotation is needed.

### Reference genome

Also place an indexed UCSC mm10 reference at `input_files/motif_annotation/mm10.fa`. To obtain it:

```bash
mkdir -p input_files/motif_annotation
curl -fL https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz \
  -o input_files/motif_annotation/mm10.fa.gz
gunzip input_files/motif_annotation/mm10.fa.gz
samtools faidx input_files/motif_annotation/mm10.fa
cut -f1,2 input_files/motif_annotation/mm10.fa.fai \
  > input_files/motif_annotation/mm10.chrom.sizes
```

Analysis outputs go to `output/motif_annotation/`. The model used for motif annotation is stored directly at `trained_models/model_ChromBPNet_mESC_nobias.h5`, matching the filename inside the Zenodo archive. Each script creates its output directory.

### Public ATAC data preparation

`prepare_ATAC_data.sh` is an executable example reconstructed from the original preprocessing scripts and the BAM processing history. It downloads and processes two ESC-TKO biological replicates from PMID36471082:

| Replicate | SRA run |
|-----------|---------|
| R1 | [SRR15974966](https://www.ncbi.nlm.nih.gov/sra/?term=SRR15974966) |
| R2 | [SRR15974967](https://www.ncbi.nlm.nih.gov/sra/?term=SRR15974967) |

After preparing mm10 above, obtain the mm10 blacklist and run:

```bash
curl -fL https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/mm10-blacklist.v2.bed.gz \
  -o input_files/motif_annotation/mm10-blacklist.v2.bed.gz
gunzip input_files/motif_annotation/mm10-blacklist.v2.bed.gz
bash Processing_code/motif_annotation/prepare_ATAC_data.sh
```

The example downloads paired FASTQs, trims Nextera adapters, maps to mm10 with Bowtie2 (`--local --very-sensitive-local -X 2000`), keeps properly paired alignments with MAPQ ≥20, removes mitochondrial alignments, removes duplicates per replicate, and excludes blacklist-overlapping alignments. Replicates are then merged and indexed. ChrX and chrY are retained. The script requires a `picard` command on PATH.

The training BAM is written to `output/motif_annotation/atac_preprocessing/mESC_ATAC_merged.bam`; the training script uses this path by default. `THREADS` and `TRIM_CORES` are set to 8 at the top of the script. FASTQs and intermediate BAMs are retained, so allow sufficient working disk space. No Tn5 shift is applied here; ChromBPNet handles its own signal preprocessing.

MACS2 calls each replicate with `-f BAMPE -g mm -q 0.05`. The script also writes their overlap as `peaks/mESC_ATAC_replicate_overlap.bed`, using the original `bedtools intersect` behavior. Edit `CALL_PEAKS=0` in the script to skip peak calling. These peaks document the ATAC preprocessing; model training uses the supplied bait-region-based accessible/background files, which are not regenerated by this mock example.

This template was syntax-checked and checked against the recorded commands; the full public-data download and preprocessing have not been rerun here.

### Training

Use the ChromBPNet environment:

```bash
bash Processing_code/motif_annotation/train_ChromBPNet.sh
```

To use another already processed BAM, edit `ATAC_BAM` at the top of `train_ChromBPNet.sh`.

Both models use fold0: test chr1/chr3/chr6; validation chr8/chr16; remaining listed chromosomes for training. Inputs are 2114 bp and predictions 1000 bp. Training uses learning rate 0.001, batch size 64, seed 1234, up to 50 epochs and patience 5. The bias model uses 128 filters, four dilation layers, no jitter and bias threshold 0.2. ChromBPNet uses 512 filters, eight dilation layers, up to 500 bp jitter and a negative sampling ratio of 0.1.

### Contributions

Use the attribution environment:

```bash
python Processing_code/motif_annotation/compute_ChromBPNet_DeepSHAP.py
```

This attributes the count head using 20 whole-sequence dinucleotide shuffles and seed 0. The script centers 2114 bp windows on the annotation BED, skips incomplete or non-ACGT sequences, and writes `ohe.npz`, `attr.npz`, `attr_hypo.npz`, and `region_ids.npz`, each with the `arr_0` key. It creates the Fi-NeMo narrowPeak file from the same retained rows, so peak IDs remain aligned. Annotation covers the supplied accessible regions across all chromosome splits; this is not a held-out performance evaluation.

To skip training, download `model_ChromBPNet_mESC.tar.gz` from Zenodo and extract it directly into `trained_models/`:

```bash
mkdir -p trained_models
tar -xzf model_ChromBPNet_mESC.tar.gz -C trained_models
python Processing_code/motif_annotation/compute_ChromBPNet_DeepSHAP.py
```

This creates `trained_models/model_ChromBPNet_mESC_nobias.h5`. The training script copies its newly trained model to the same path, and the contribution script always reads that file.

### Motif discovery and names

Use the MoDISco environment, with MEME Suite's `tomtom` available on PATH:

```bash
bash Processing_code/motif_annotation/run_ChromBPNet_modisco.sh
```

MoDISco uses hypothetical contributions and a maximum of 50,000 seqlets per metacluster within the central 400 bp (the original discovery settings). Fi-NeMo subsequently scans the complete 2114 bp attribution windows. Open `output/motif_annotation/modisco/report/motifs.html` and review `modisco/motif_names.tsv`. The table retains the best JASPAR match and q-value; its `motif` column is editable. A best match is a proposed label, including when the match is weak or shared by several TFs. Pattern numbers are specific to each discovery run; assign labels from that run's report.

```bash
python Processing_code/motif_annotation/trim_motif_CWMs.py
```

Trimming uses a five-base moving mean of absolute CWM importance and retains the span exceeding 30% of its maximum. It writes the trimmed MEME file and individual logos. Fi-NeMo reads the original MoDISco H5; the separately trimmed CWM is for motif inspection and downstream uses, and does not change Fi-NeMo's hit coordinates.

### Motif instances and annotation

Use the Fi-NeMo environment:

```bash
bash Processing_code/motif_annotation/run_ChromBPNet_finemo.sh
```

Then use the attribution environment:

```bash
python Processing_code/motif_annotation/annotate_finemo_hits.py
```

The final file is `output/motif_annotation/finemo/annotated_hits.tsv`. Genomic hit coordinates and scores are retained. Each 2114 bp attribution window is assigned to the 2048 bp BED interval with greatest overlap; ties follow BED order. Unmatched regions remain empty. Assigning a region label does not clip hits to that region. Positive and negative pattern IDs remain distinct when mapping TF names.

For an existing, reviewed named MEME file, set `NAMED_MEME_FILE = Path("path/to/named_motifs.meme")` at the top of `annotate_finemo_hits.py` to read its names instead of the TSV. Legacy names such as `pattern_0_CTCF` refer to positive patterns. Use that file only with the MoDISco results from which those names were assigned.

## Dependencies

Recorded local versions are listed below. Keep TensorFlow training separate from PyTorch attribution and Fi-NeMo. Training and attribution benefit from a GPU; the complete attribution and motif-discovery calculations are substantial jobs. The scripts contain no cluster-specific Slurm or environment activation commands.

### ATAC preprocessing

Versions recorded in the original scripts, BAM headers and trimming logs:

| Package | Version |
|---------|---------|
| SRA Toolkit | 3.0.0 |
| Trim Galore | 0.6.10 |
| Cutadapt | 5.0 |
| Bowtie2 | 2.5.4 |
| samtools | 1.21 |
| Picard | 3.1.0 |
| Java | 17 |
| bedtools | 2.30.0 |
| MACS2 | 2.2.7.1 |

Use separate environments if needed: the recorded Cutadapt run used Python 3.12.9, whereas the original MACS2 module supplied its own environment.

### Training: Python 3.8.20

| Package | Version |
|---------|---------|
| chrombpnet | 0.1.7 |
| tensorflow | 2.8.0 |
| numpy | 1.23.4 |
| pandas | 2.0.3 |
| h5py | 3.11.0 |
| pyfaidx | 0.6.1 |
| scipy | 1.10.1 |
| samtools | 1.21 |
| bedtools | 2.31.1 |
| ucsc-bedgraphtobigwig | 481 |
| meme | 5.5.8 |

### Attribution and annotation: Python 3.10.19

| Package | Version |
|---------|---------|
| bpnet-lite | 0.9.1 |
| torch | 2.5.1+cu121 |
| numpy | 1.26.4 |
| pandas | 2.3.3 |
| h5py | 3.15.1 |
| pyfaidx | 0.8.1.4 |

### Motif discovery and trimming: Python 3.9.0

| Package | Version |
|---------|---------|
| modisco-lite | 2.3.2 |
| meme | 5.5.2 |
| numpy | 2.0.2 |
| pandas | 2.2.3 |
| h5py | 3.13.0 |
| matplotlib | 3.9.4 |
| logomaker | 0.8.6 |
| numba | 0.60.0 |
| scipy | 1.13.1 |
| hdf5plugin | 5.0.0 |

### Fi-NeMo: Python 3.10.19

| Package | Version |
|---------|---------|
| FiNeMo | 0.30 |
| polars | 1.38.1 |
| torch | 2.5.1+cu121 |
| numpy | 1.26.4 |
| h5py | 3.15.1 |
| pyfaidx | 0.8.1.4 |
| scipy | 1.15.3 |
| matplotlib | 3.10.7 |
| hdf5plugin | 6.0.0 |

## Validation

The cleaned trimming code reproduces all 26 existing trimmed CWM matrices. Hit coordinates, scores and region assignments match all 1,102,841 existing hits; TF labels match the supplied reviewed name mapping. A one-region CPU test exercised model loading, both DeepLIFT-SHAP outputs and Fi-NeMo input preparation. The training commands were checked against the installed ChromBPNet parser. Full training and motif discovery were not rerun while preparing this example.
