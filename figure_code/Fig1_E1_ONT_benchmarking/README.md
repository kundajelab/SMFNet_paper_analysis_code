# mESC SMF-ONT model benchmarking

## Overview

- **`plot_sox2_example.R`**: Plots Sox2 locus examples for Fig. 1c-d.
- **`plot_SMF-ONT_qc.R`**: Plots experimental QC for mESCH DNMT-TKO SMF-ONT data for Fig. E1b-d.
- **`plot_SMFNet_pos_comp.R`**: Compares SMFNet-predicted and experimental occupancy at individual positions for Fig. 1e.
- **`plot_SMFNet_metrics.R`**: Plots performance metrics for the mESC DNMT-TKO SMF-ONT model for Fig. E1e, i-j.
- **`count_SMFNet_motif_states.py`, `compute_SMFNet_motif_avg.py`, and `plot_SMFNet_motif_states.R`**: Calculate state percentages at putative annotated TFBSs and average enrichment around annotated motifs, and generate composite plots and motif state correlations for Fig. E1f-h, k.

## Instructions

Run all commands from the repository root. Download the required archives from Zenodo and extract their contents into:

```text
input_files/Fig1_E1_ONT_benchmarking/
```

### Sox2 locus examples

**Required archive:** `Example_Sox2_locus_data.tar.gz`

```bash
Rscript figure_code/Fig1_E1_ONT_benchmarking/plot_sox2_example.R
```

### Experimental mESC DNMT-TKO SMF-ONT QC

**Required archive:** `WT_SMF-ONT_qc_data.tar.gz`

```bash
Rscript figure_code/Fig1_E1_ONT_benchmarking/plot_SMF-ONT_qc.R
```

### QC at individual SMF positions

**Required archive:** `SMFNet_mESC_basic_qc.tar.gz`

```bash
Rscript figure_code/Fig1_E1_ONT_benchmarking/plot_SMFNet_pos_comp.R
```

### SMFNet performance metrics for mESC DNMT-TKO SMF-ONT

**Required archive:** `SMFNet_mESC_basic_qc.tar.gz`

```bash
Rscript figure_code/Fig1_E1_ONT_benchmarking/plot_SMFNet_metrics.R
```

### SMFNet motif QC for mESC DNMT-TKO SMF-ONT

**Required archives:**

- `SMFNet_mESC_basic_qc.tar.gz`
- `mESC_SMF-ONT_testRegions_reads.tar.gz`
- `mESC_SMF-ONT_testRegions_states.tar.gz`

Run the following commands in order:

```bash
python figure_code/Fig1_E1_ONT_benchmarking/count_SMFNet_motif_states.py
python figure_code/Fig1_E1_ONT_benchmarking/compute_SMFNet_motif_avg.py
Rscript figure_code/Fig1_E1_ONT_benchmarking/plot_SMFNet_motif_states.R
```

## Dependencies

### Python 3.12

| Package | Version |
|---------|---------|
| numpy | 2.2.6 |
| pandas | 3.0.5 |
| h5py | 3.15.1 |
| tqdm | 4.67.1 |

### R 4.4.1

| Package | Version |
|---------|---------|
| tidyverse | 2.0.0 |
| ggplot2 | 4.0.2 |
| egg | 0.4.5 |
| gridExtra | 2.3 |
| hdf5r | 1.3.16 |
| data.table | 1.18.6.1 |
| viridis | 0.6.5 |
| ggrastr | 1.0.2 |
| dplyr | 1.2.0 |
| tidyr | 1.3.2 |
