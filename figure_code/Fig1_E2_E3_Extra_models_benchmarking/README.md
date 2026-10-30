# Extra SMFNet models Benchmarking

## Overview

- **`plot_Sonmezer_example.R`**: Plots example region for short-read SMFNet model for Fig. E2c.
- **`compute_Sonmezer_ONT_metrics.py`, and `plot_Sonmezer_qc.R`**: Calculate quality control metrics between ONT and short-read mESC SMFNet models. Plots quality control metrics for the short-read mESC SMFNet model for Fig. E2a-b,d-f.

## Instructions

Run all commands from the repository root. Download the required archives from Zenodo and extract their contents into:

```text
input_files/Fig1_E2_E3_Extra_models_benchmarking/
```

### Short-read mESC SMFnet example

**Required archive:** `mESC_Sonmezer_testRegions_data.tar.gz`

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
| h5py | 3.16.0 |
| tqdm | 4.70.1 |

### R 4.4

| Package | Version |
|---------|---------|
| ggplot2 | 4.0.3 |
| egg | 0.4.5 |
| gridExtra | 2.3.1 |
| hdf5r | 1.3.16 |
| data.table | 1.18.6.1 |
| viridis | 0.6.5 |
| ggrastr | 1.0.2 |
| dplyr | 1.2.1 |
| tidyr | 1.3.2 |

