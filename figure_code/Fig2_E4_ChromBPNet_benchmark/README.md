# SMFNet ChromBPNet benchmarking

## Overview

- **`plot_marginalizations.R`**: Plots ChromBPNet and SMFNet motif marginalization profiles for Fig. E4a-c.
- **`plot_GpC_marginalizations.R`**: Plots motif marginalizations from the mESC, NPC, C2C12 and MEL GpC-only models for Fig. E4e.
- **`plot_SMFNet_ISM_heatmaps.R`**: Plots endogenous SMFNet ISM heatmaps, average profiles and sampled molecules for Fig. 2a + E4g.
- **`plot_ISM_correlations.R`**: Plots ChromBPNet vs SMFNet ISM correlations across TFBSs for Fig. E4d.

## Instructions

Run all commands from the repository root. Place the figure-specific input files in:

```text
input_files/Fig2_E4_ChromBPNet_benchmark/
```

Outputs are written to `output/Fig2_E4_ChromBPNet_benchmark/`. The scripts take no command-line arguments. 

### ChromBPNet and SMFNet motif marginalization

**Required files:**

- `ChromBPNet_marginalization_profiles.tsv`
- `SMFNet_marginalization_profiles.tsv`
- `SMFNet_marginalization_states.h5`

```bash
Rscript figure_code/Fig2_E4_ChromBPNet_benchmark/plot_marginalizations.R
```

### GpC-only model motif marginalization

**Required archive:**

- `GpC_marginalizations.tar.gz`

```bash
Rscript figure_code/Fig2_E4_ChromBPNet_benchmark/plot_GpC_marginalizations.R
```

### ChromBPNet and SMFNet total marginalization correlation

**Required archive:**

- `ChromBPNet_SMFNet_marginalizations.tar.gz`

```bash
Rscript figure_code/Fig2_E4_ChromBPNet_benchmark/plot_marginalizations.R
```

### GpC-only model motif marginalization

**Required archive:**

- `GpC_only_marginalizations.tar.gz`

```bash
Rscript figure_code/Fig2_E4_ChromBPNet_benchmark/plot_GpC_marginalizations.R
```

### Endogenous SMFNet ISM heatmaps

**Required archives:**

- `Endogenous_ISM_inputs.tar.gz`
- `mESC_SMF-ONT_testRegions_reads.tar.gz`
- `SMFNet_mESC_basic_qc.tar.gz`
- `SMFNet_annotated_motif_hits_for_codependency_analysis.tar.gz`

```bash
Rscript figure_code/Fig2_E4_ChromBPNet_benchmark/plot_SMFNet_ISM_heatmaps.R
```

### ChromBPNet vs SMFNet endogenous ISM

**Required archives:**

- `Endogenous_ISM_inputs.tar.gz`
- `model_ChromBPNet_mESC.tar.gz`
- `SMFNet_annotated_motif_hits_for_codependency_analysis.tar.gz`


Run the following commands in order:

```bash
python figure_code/Fig2_E4_ChromBPNet_benchmark/compute_ChromBPNet_ISM.py
Rscript figure_code/Fig2_E4_ChromBPNet_benchmark/plot_ISM_correlations.R
```

## Dependencies

### Python 3.12

| Package | Version |
|---------|---------|
| numpy | 1.26.4 |
| torch | 2.5.1+cu121 |
| h5py | 3.15.1 |
| matplotlib | 3.10.7 |
| pyfaidx | 0.8.1.4 |
| bpnet-lite | 0.9.1 |
| tangermeme | 0.5.1 |

### R 4.4

| Package | Version |
|---------|---------|
| ggplot2 | 4.0.2 |
| dplyr | 1.2.0 |
| tidyr | 1.3.2 |
| stringr | 1.6.0 |
| egg | 0.4.5 |
| gridExtra | 2.3 |
| hdf5r | 1.3.16 |
| zoo | 1.8.15 |
| GenomicRanges | 1.58.0 |
