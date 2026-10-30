# SMFNet co-dependency benchmarking

## Overview

- **`plot_Sonmezer_example.R`**: Plots example region for short-read SMFNet model for Fig. E2c.

## Instructions

Run all commands from the repository root. Download the required archives from Zenodo and extract their contents into:

```text
input_files/Fig1_E2_E3_Extra_models_benchmarking/
```

### Short-read mESC SMFnet example

**Required archives:** 
- `mESC_Sonmezer_testRegions_data.tar.gz`
- `mESC_ChromBPnet_annotated_hits.tar.gz`

```bash
Rscript figure_code/Fig1_E2_E3_Extra_models_benchmarking/plot_Sonmezer_example.R
```

### mESC short-read SMFNet QC

**Required archives:** 
- `mESC_Sonmezer_testRegions_data.tar.gz`
- `mESC_SMF-ONT_testRegions_reads.tar.gz`

Run the following commands in order:
```bash
python figure_code/Fig1_E2_E3_Extra_models_benchmarking/compute_Sonmezer_ONT_metrics.py
Rscript figure_code/Fig1_E2_E3_Extra_models_benchmarking/plot_Sonmezer_qc.R
```

### GpC-only SMFNet models QC

**Required archive:** 
`Extra_models_and_Fiber-seq_inputs.tar.gz`

```bash
Rscript figure_code/Fig1_E2_E3_Extra_models_benchmarking/plot_Kreibich_metrics.R
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

### R 4.4

| Package | Version |
|---------|---------|
| ggplot2 | 4.0.3 |
| egg | 0.4.5 |
| gridExtra | 2.3.1 |
| hdf5r | 1.3.16 |
| tibble | 3.3.1 |
| ggrastr | 1.0.2 |
| dplyr | 1.2.1 |
| tidyr | 1.3.2 |

