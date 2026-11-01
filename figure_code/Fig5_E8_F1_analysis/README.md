# SMFNet F1 variant benchmarking

## Overview

- **`plot_F1_qc.R`**: Read lengths, replicate correlations and F1 versus wildtype correlations for Fig. E8b–d. The counts and coverage in Fig. E8a are supplied in `read_length_N50.tsv`.
- **`plot_variant_metrics.R`**: Median and interquartile changes around disrupted motifs for Fig. 5e.
- **`plot_accessibility_prediction.R`**: Accessible-fraction changes across variant windows for Fig. 5d.
- **`plot_microstate_prediction.R`**: Microstate-fraction changes for Fig. E8g.
- **`plot_microstate_examples.R`**: Average profiles, single-molecule microstates and fraction changes at chr6:124399844–124400244 for Fig. 5b–c and E8e–f.
- **`plot_F1_vs_single_motif.R`**: Frequency, width and phasing correlations for Fig. 5f–h. The insilico single-motif summary is supplied as `insilico_metrics_sum.txt`.
- **`plot_disrupted_examples.R`**: Single-molecule states and average profiles at chr1:46905014–46905414 and chr7:27622166–27622566 for Fig. 5i–j and E8h–i.

## Instructions

Run all commands from the repository root. Download the required archives from Zenodo and extract their contents into:

```text
input_files/Fig5_E8_F1_analysis/
```

Outputs are written to `output/Fig5_E8_F1_analysis/`. Completed prerequisite calculations can be reused.

### F1 sequencing and replicate QC

**Required archives:**

- `F1_QC_inputs.tar.gz`

Run the following commands in order:

```bash
Rscript figure_code/Fig5_E8_F1_analysis/plot_F1_qc.R
```

### Variant effects on accessibility, width and phasing

**Required archives:**

- `F1_states.tar.gz`
- `F1_variant_inputs.tar.gz`

Run the following commands in order:

```bash
python figure_code/Fig5_E8_F1_analysis/compute_variant_metrics.py
Rscript figure_code/Fig5_E8_F1_analysis/plot_variant_metrics.R
```

### Experimental versus predicted accessibility changes

**Required archives:**

- `F1_states.tar.gz`
- `F1_variant_inputs.tar.gz`

Run the following commands in order:

```bash
python figure_code/Fig5_E8_F1_analysis/compute_variant_metrics.py
Rscript figure_code/Fig5_E8_F1_analysis/plot_variant_metrics.R
python figure_code/Fig5_E8_F1_analysis/make_variant_windows.py
python figure_code/Fig5_E8_F1_analysis/compute_accessibility_windows.py
Rscript figure_code/Fig5_E8_F1_analysis/plot_accessibility_prediction.R
```

### Experimental versus predicted microstate changes

**Required archives:**

- `F1_states.tar.gz`
- `F1_variant_inputs.tar.gz`

Run the following commands in order:

```bash
python figure_code/Fig5_E8_F1_analysis/compute_variant_metrics.py
Rscript figure_code/Fig5_E8_F1_analysis/plot_variant_metrics.R
python figure_code/Fig5_E8_F1_analysis/make_variant_windows.py
python figure_code/Fig5_E8_F1_analysis/compute_accessibility_windows.py
python figure_code/Fig5_E8_F1_analysis/cluster_variant_windows.py
Rscript figure_code/Fig5_E8_F1_analysis/plot_microstate_prediction.R
```

### KLF microstate example

**Required archives:**

- `F1_states.tar.gz`
- `F1_variant_inputs.tar.gz`

Run the following commands in order:

```bash
python figure_code/Fig5_E8_F1_analysis/compute_variant_metrics.py
Rscript figure_code/Fig5_E8_F1_analysis/plot_variant_metrics.R
python figure_code/Fig5_E8_F1_analysis/make_variant_windows.py
python figure_code/Fig5_E8_F1_analysis/compute_accessibility_windows.py
python figure_code/Fig5_E8_F1_analysis/cluster_variant_windows.py
Rscript figure_code/Fig5_E8_F1_analysis/plot_microstate_examples.R
```

### F1 disruption versus single-motif insertion

**Required archives:**

- `F1_states.tar.gz`
- `F1_variant_inputs.tar.gz`

Run the following commands in order:

```bash
python figure_code/Fig5_E8_F1_analysis/compute_variant_metrics.py
Rscript figure_code/Fig5_E8_F1_analysis/plot_variant_metrics.R
Rscript figure_code/Fig5_E8_F1_analysis/plot_F1_vs_single_motif.R
```

### KLF and OCT4::SOX2 disruption examples

**Required archives:**

- `F1_states.tar.gz`
- `F1_variant_inputs.tar.gz`

Run the following commands in order:

```bash
python figure_code/Fig5_E8_F1_analysis/compute_variant_metrics.py
python figure_code/Fig5_E8_F1_analysis/extract_disrupted_example_windows.py
Rscript figure_code/Fig5_E8_F1_analysis/plot_disrupted_examples.R
```

## Dependencies

### Python 3.12

| Package | Version |
|---------|---------|
| numpy | 2.2.6 |
| pandas | 3.0.5 |
| scipy | 1.15.3 |
| numba | 0.62.1 |
| h5py | 3.15.1 |
| scikit-learn | 1.7.2 |
| tqdm | 4.67.1 |

**NumPy exception:** For `cluster_variant_windows.py`, use NumPy **1.26.4** with scikit-learn **1.7.2** in a separate Python 3.12 environment. Keep the supplied random seed and input-window ordering.

### R 4.4.1

| Package | Version |
|---------|---------|
| tidyverse | 2.0.0 |
| ggplot2 | 4.0.2 |
| dplyr | 1.2.0 |
| tidyr | 1.3.2 |
| tibble | 3.3.1 |
| stringr | 1.6.0 |
| forcats | 1.0.1 |
| egg | 0.4.5 |
| gridExtra | 2.3 |
| data.table | 1.18.6.1 |
| viridis | 0.6.5 |
| hdf5r | 1.3.16 |
| ggrastr | 1.0.2 |
| hexbin | 1.28.5 |
| scales | 1.4.0 |
| grid | 4.4.1 |
