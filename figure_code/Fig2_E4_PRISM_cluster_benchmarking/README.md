# SMFNet PRISM and cluster benchmarking

## Overview

- **`plot_CTCF_PRISM.R`**: Plots raw PRISM scores for individual CTCF-example reads for Fig. 2b.
- **`plot_PRISM_validation.R`**: Plots PRISM mutation scores against methylation agreement for Fig. E4h-i.
- **`plot_enhancer_clusters.R`**: Plots endogenous motif-cluster profiles, single molecules, cluster fractions and experiment-versus-SMFNet comparisons.

## Instructions

Run all commands from the repository root. Download the required archives from Zenodo and extract their contents into:

```text
input_files/Fig2_E4_PRISM_cluster_benchmarking/
```

Extract the listed archives to `input_files/Fig2_E4_PRISM_cluster_benchmarking/`. Output is written to `output/Fig2_E4_PRISM_cluster_benchmarking/`.

### CTCF PRISM single-read examples

**Required archive:**

- `CTCF_PRISM_results.tar.gz`

```bash
Rscript figure_code/Fig2_E4_PRISM_cluster_benchmarking/plot_CTCF_PRISM.R
```

### PRISM validation

**Required archive:**

- `PRISM_validation.tar.gz`

```bash
Rscript figure_code/Fig2_E4_PRISM_cluster_benchmarking/plot_PRISM_validation.R
```

Uses the supplied validation results for real reads and the fraction-agreement kernel. Regenerating the CTCF example is independent of this validation analysis.

### Endogenous motif-cluster analysis

**Required archives:**

- `SMFNet_annotated_motif_hits_for_codependency_analysis.tar.gz`
- `mESC_SMF-ONT_testRegions_reads.tar.gz`
- `mESC_SMF-ONT_testRegions_states.tar.gz`

Run the following commands in order:

```bash
python figure_code/Fig2_E4_PRISM_cluster_benchmarking/prepare_cluster_TFBS.py
python figure_code/Fig2_E4_PRISM_cluster_benchmarking/prepare_cluster_regions.py
python figure_code/Fig2_E4_PRISM_cluster_benchmarking/quantify_enhancer_clusters.py
Rscript figure_code/Fig2_E4_PRISM_cluster_benchmarking/plot_enhancer_clusters.R
```

## Dependencies

### Python 3.12

| Package | Version |
|---------|---------|
| numpy | 2.2.6 |
| pandas | 3.0.5 |
| h5py | 3.16.0 |
| tqdm | 4.70.1 |
| scikit-learn | 1.9.1 |
| matplotlib | 3.10.3 |
| intervaltree | 3.1.0 |

The optional PRISM generation script additionally requires PyTorch and SMFNet with its checkpoint dependencies. GPU inference has not been rerun as part of this cleanup.

### R 4.4

| Package | Version |
| --------- | --------- |
| dplyr | 1.2.1 |
| tidyr | 1.3.2 |
| ggplot2 | 4.0.3 |
| egg | 0.4.5 |
| gridExtra | 2.3.1 |
| ggseqlogo | 0.2.2 |
| Biostrings | 2.74.0 |
