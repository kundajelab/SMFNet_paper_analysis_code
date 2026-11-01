# SMFNet combined-motif benchmarking

## Overview

- **`plot_baselines.R`**: Plots frequency scatters and center-motif residual boxplot for Fig. 4b–c, as well as frequency, width and phasing correlations and mean errors for Fig. E7a–b, d–g.
- **`plot_SOX2_flanks.R`**: Plots SOX2-flank effects over motif spacing for Fig. 4e and E7c.
- **`plot_state_profiles.R`**: Plots SOX2 flanks at 20/60 bp and the three enhancer examples at 50 bp for Fig. 4d,i.
- **`plot_accessible_SMF_profiles.R`**: Plots accessible-conditioned SMF profiles with OCY4::SOX2 center and KLF or ESRRB flanks at 20/60 bp for Fig. 4g and E7i.
- **`plot_OCT4_footprints.R`**: Plots all-molecule and accessible-only footprint changes over spacing, as well as the BPNet comparison for Fig. 4h and E7h,j.
- **`plot_phasing.R`**: Plots boundary disperions versus frequency, residual boundary dispersion over distance and the spacing boxplot for Fig. 4j and E7k–m.

## Instructions

Run all commands from the repository root. Download the required archives from Zenodo and extract their contents into:

```text
input_files/Fig4_E7_combinedMotifs_insilico/
```

Outputs are written to `output/Fig4_E7_combinedMotifs_insilico/`. Run `compute_motif_metrics.py` once and reuse its results for the baseline, SOX2-flank and phasing analyses. Completed prerequisites do not need to be rerun. The existing input archives are unchanged.

### Combined prediction versus baselines

**Required archives:**

- `Three_motif_states.tar.gz`
- `Three_motif_flank_controls.tar.gz`
- `Combined_motif_analysis_inputs.tar.gz`

Run the following commands in order:

```bash
python figure_code/Fig4_E7_combinedMotifs_insilico/compute_motif_metrics.py
python figure_code/Fig4_E7_combinedMotifs_insilico/compute_baseline_statistics.py
Rscript figure_code/Fig4_E7_combinedMotifs_insilico/plot_baselines.R
```

### SOX2 flanks versus no flanks

**Required archives:**

- `Three_motif_states.tar.gz`
- `Three_motif_flank_controls.tar.gz`
- `Combined_motif_analysis_inputs.tar.gz`

Run the following commands in order:

```bash
python figure_code/Fig4_E7_combinedMotifs_insilico/compute_motif_metrics.py
Rscript figure_code/Fig4_E7_combinedMotifs_insilico/plot_SOX2_flanks.R
```

### Average and binned state profiles

**Required archives:**

- `Three_motif_states.tar.gz`
- `mESC_enhancer_example_states.tar.gz`

Run the following commands in order:

```bash
python figure_code/Fig4_E7_combinedMotifs_insilico/compute_state_profiles.py
Rscript figure_code/Fig4_E7_combinedMotifs_insilico/plot_state_profiles.R
```

### Accessible-molecule SMF profiles

**Required archives:**

- `Three_motif_states.tar.gz`
- `OCT4_footprint_sampled_reads.tar.gz`

Run the following commands in order:

```bash
python figure_code/Fig4_E7_combinedMotifs_insilico/compute_accessible_SMF_profiles.py
Rscript figure_code/Fig4_E7_combinedMotifs_insilico/plot_accessible_SMF_profiles.R
```

### OCT4 footprints and BPNet comparison

**Required archives:**

- `Three_motif_states.tar.gz`
- `OCT4_footprint_sampled_reads.tar.gz`
- `Combined_motif_analysis_inputs.tar.gz`

Run the following commands in order:

```bash
python figure_code/Fig4_E7_combinedMotifs_insilico/compute_OCT4_footprints.py
Rscript figure_code/Fig4_E7_combinedMotifs_insilico/plot_OCT4_footprints.R
```

### Nucleosome phasing

**Required archives:**

- `Three_motif_states.tar.gz`
- `Three_motif_flank_controls.tar.gz`
- `Combined_motif_analysis_inputs.tar.gz`

Run the following commands in order:

```bash
python figure_code/Fig4_E7_combinedMotifs_insilico/compute_motif_metrics.py
python figure_code/Fig4_E7_combinedMotifs_insilico/compute_phasing.py
Rscript figure_code/Fig4_E7_combinedMotifs_insilico/plot_phasing.R
```

## Dependencies

Versions used for validation:

### Python 3.10

| Package | Version |
|---------|---------|
| numpy | 1.26.4 |
| pandas | 2.3.3 |
| scipy | 1.15.3 |
| h5py | 3.15.1 |
| numba | 0.62.1 |
| tqdm | 4.67.1 |

### R 4.4

| Package | Version |
|---------|---------|
| data.table | 1.18.6.1 |
| ggplot2 | 4.0.2 |
| dplyr | 1.2.0 |
| tibble | 3.3.1 |
| egg | 0.4.5 |
| gridExtra | 2.3 |
| viridis | 0.6.5 |

