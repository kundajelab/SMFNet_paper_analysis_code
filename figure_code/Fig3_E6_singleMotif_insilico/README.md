# Single-motif, GC-content and homotypic in silico analysis

## Overview

- **`quantify_single_motifs.py`** and **`compute_single_motif_profiles.py`**: Calculate opening frequency, median accessible width, boundary dispersion and profiles for central single motifs and matched no-motif backgrounds.
- **`quantify_gc_motifs.py`**: Calculates the same metrics across the GC-content screen, including no-motif controls.
- **`quantify_homotypic_motifs.py`** and **`compute_homotypic_profiles.py`**: Quantify the homotypic libraries at 20–80-bp spacing and calculate the 50-bp profiles.
- **`plot_single_motifs.R`**: Plot in-silico single motif insertion metrics and bin plots for Fig. 3b-c + E6d.
- **`plot_gc_motifs.R`**: Plot GC-context heatmaps for Fig. E6c.
- **`plot_homotypic_motifs.R`**: Plot homotypic number, spacing heatmaps for Fig. 3e + E6e-g
- **`plot_homotypic_profiles.R`**: Plot homotypic number bin plots for Fig. 3f.
- **`plot_chip_gc.R`**: Calculate GC content around chip Peaks and plots violin for Fig. E6b.
- **`plot_motif_logos.py`**: Plots the selected CWMs using a maximum trim length of 16 bp for Fig. E6a.

## Instructions

Run all commands from the repository root. Download the required archives from Zenodo and extract their contents into:

```text
input_files/Fig3_E6_singleMotif_insilico/
```

Calculations and plots are written to `output/Fig3_E6_singleMotif_insilico/`.


### Single-motif metrics and profiles

**Required archive:**

- `SMFNet-ONT_insilico_single_motif_states.tar.gz`

Run the following commands in order:

```bash
python figure_code/Fig3_E6_singleMotif_insilico/quantify_single_motifs.py
python figure_code/Fig3_E6_singleMotif_insilico/compute_single_motif_profiles.py
Rscript figure_code/Fig3_E6_singleMotif_insilico/plot_single_motifs.R
```

### Homotypic motif metrics

**Required archive:**

- `SMFNet-ONT_insilico_homotypic_states.tar.gz`

Run the following commands in order:

```bash
python figure_code/Fig3_E6_singleMotif_insilico/quantify_homotypic_motifs.py
Rscript figure_code/Fig3_E6_singleMotif_insilico/plot_homotypic_motifs.R
```

### Homotypic motif bin plots

**Required archive:**

- `SMFNet-ONT_insilico_homotypic_states.tar.gz`

Run the following commands in order:

```bash
python figure_code/Fig3_E6_singleMotif_insilico/compute_homotypic_profiles.py
Rscript figure_code/Fig3_E6_singleMotif_insilico/plot_homotypic_profiles.R
```

### GC-content motif analysis

**Required archive:**

- `SMFNet-ONT_insilico_gc_inputs.tar.gz`

Run the following commands in order:

```bash
python figure_code/Fig3_E6_singleMotif_insilico/quantify_gc_motifs.py
Rscript figure_code/Fig3_E6_singleMotif_insilico/plot_gc_motifs.R
```

### ChIP-peak GC content

**Required archive:**

- `mESC_ChIP_narrowPeaks.tar.gz`

```bash
Rscript figure_code/Fig3_E6_singleMotif_insilico/plot_chip_gc.R
```

### Motif logos

**Required archive:**

- `insilico_motif_CWMs.tar.gz`

```bash
python figure_code/Fig3_E6_singleMotif_insilico/plot_motif_logos.py
```

## Dependencies

### Python 3.12

| Package | Version |
|---------|---------|
| numpy | 2.2.6 |
| pandas | 3.0.5 |
| h5py | 3.15.1 |
| numba | 0.62.1 |
| tqdm | 4.67.1 |
| matplotlib | 3.10.7 |
| tangermeme | 0.5.1 |
| torch | 2.5.1+cu121 |

**NumPy exception:** For `compute_single_motif_profiles.py` and `compute_homotypic_profiles.py`, use NumPy **1.26.4** in a separate Python 3.12 environment to preserve the original ordering of tied reads between bins.

### R 4.4.1

| Package | Version |
|---------|---------|
| tidyverse | 2.0.0 |
| dplyr | 1.2.0 |
| tidyr | 1.3.2 |
| ggplot2 | 4.0.2 |
| stringr | 1.6.0 |
| readr | 2.1.6 |
| egg | 0.4.5 |
| gridExtra | 2.3 |
| data.table | 1.18.6.1 |
| viridis | 0.6.5 |
| ggh4x | 0.3.1 |
| BSgenome.Mmusculus.UCSC.mm10 | 1.4.3 |
| Biostrings | 2.74.1 |
| GenomicRanges | 1.58.0 |
