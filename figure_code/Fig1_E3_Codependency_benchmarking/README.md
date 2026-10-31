# SMFNet co-dependency benchmarking

## Overview

- **`plot_coaccessibility_matrices.R`**: Plots smoothed experiment/model co-accessibility matrices for Fig. E3d.
- **`plot_coaccessibility_distance.R`**: Plots the distance vs co-accessibility score GAM fit for Fig. E3e.
- **`plot_CTCF_coaccessibility.R`**: Plots the CTCF-aligned co-accessibility summary for bound/unbound molecules for Fig. 1f.

## Instructions

Run all commands from the repository root. Download the required archives from Zenodo and extract their contents into:

```text
input_files/Fig1_E3_Codependency_benchmarking/
```

Intermediate files are generated automatically under `output/Fig1_E3_Codependency_benchmarking/`.

### General co-accessibility matrices

**Required archive:**
- `mESC_SMF-ONT_testRegions_states.tar.gz`

**Required files:**
- `coaccessibility_region_coverage.tsv`
- `coaccessibility_background_bias.txt`

Run the following commands in order:
```bash
python figure_code/Fig1_E3_Codependency_benchmarking/compute_coaccessibility.py
python figure_code/Fig1_E3_Codependency_benchmarking/smooth_coaccessibility_matrices.py
Rscript figure_code/Fig1_E3_Codependency_benchmarking/plot_coaccessibility_matrices.R
```

### Distance vs co-accessibility score GAM fit

**Required archive:**
- `mESC_SMF-ONT_testRegions_states.tar.gz`

**Required file:**
- `coaccessibility_region_coverage.tsv`

Run the following commands in order:
```bash
python figure_code/Fig1_E3_Codependency_benchmarking/compute_coaccessibility.py
python figure_code/Fig1_E3_Codependency_benchmarking/extract_coaccessibility_distances.py
Rscript figure_code/Fig1_E3_Codependency_benchmarking/plot_coaccessibility_distance.R
```

### CTCF-bound vs unbound co-accessibility summary

**Required archives:**
- `mESC_SMF-ONT_testRegions_states.tar.gz`
- `mESC_ChromBPnet_annotated_hits.tar.gz`
- `SMFNet_mESC_basic_qc.tar.gz`

**Required file:**
- `coaccessibility_background_bias.txt`

Run the following commands in order:
```bash
python figure_code/Fig1_E3_Codependency_benchmarking/count_CTCF_motif_states.py
python figure_code/Fig1_E3_Codependency_benchmarking/compute_CTCF_coaccessibility.py
Rscript figure_code/Fig1_E3_Codependency_benchmarking/plot_CTCF_coaccessibility.R
```

## Dependencies

### Python 3.12

| Package | Version |
|---------|---------|
| numpy | 2.2.6 |
| pandas | 3.0.5 |
| h5py | 3.16.0 |
| tqdm | 4.70.1 |
| scipy | 1.15.3 |

### R 4.3

| Package | Version |
|---------|---------|
| ggplot2 | 4.0.3 |
| egg | 0.4.5 |
| gridExtra | 2.3.1 |
| dplyr | 1.2.1 |
| stringr | 1.6.0 |
| mgcv | 1.9-1 |
