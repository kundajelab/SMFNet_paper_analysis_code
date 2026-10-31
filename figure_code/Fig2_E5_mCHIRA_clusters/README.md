# SMFNet mCHIRA cluster analysis

## Overview

- **`quantify_mCHIRA_clusters.py`**: Clusters wildtype CRE_6005 reads and assigns its mutants to the same clusters.
- **`plot_mCHIRA_clusters.R`**: Plots cluster fractions, single molecules and per-cluster SMF signal averages, and perturbation effects on cluster fractions for Fig. 2c-g + E5a-c.
- **`plot_mCHIRA_PRISM.R`**: Plots wildtype PRISM sequence logos for each cluster for Fig. 2e + E5c.
- **`plot_mCHIRA_PRISM_perturbations.R`**: Compares NRF1 and MAX::MYC PRISM scores with mutant-minus-wildtype changes in cluster fractions for Fig. E5f.
- **`plot_mCHIRA_bulk_ISM.R`**: Plots bulk ISM scores for CRE_6005 and its mutants for Fig. 2c + E5d-e.

## Instructions

Run all commands from the repository root. Download `mCHIRA_CRE6005_data.tar.gz` from Zenodo and extract its contents into:

```text
input_files/Fig2_E5_mCHIRA_clusters/
```

Run the following commands in order:

```bash
python figure_code/Fig2_E5_mCHIRA_clusters/quantify_mCHIRA_clusters.py
Rscript figure_code/Fig2_E5_mCHIRA_clusters/plot_mCHIRA_clusters.R
Rscript figure_code/Fig2_E5_mCHIRA_clusters/plot_mCHIRA_PRISM.R
Rscript figure_code/Fig2_E5_mCHIRA_clusters/plot_mCHIRA_PRISM_perturbations.R
Rscript figure_code/Fig2_E5_mCHIRA_clusters/plot_mCHIRA_bulk_ISM.R
```

## Dependencies

### Python 3.10

| Package | Version |
|---------|---------|
| numpy | 1.26.4 |
| pandas | 2.3.3 |
| h5py | 3.15.1 |
| scikit-learn | 1.7.2 |

### R 4.4

| Package | Version |
|---------|---------|
| dplyr | 1.2.1 |
| tidyr | 1.3.2 |
| ggplot2 | 4.0.3 |
| egg | 0.4.5 |
| gridExtra | 2.3.1 |
| hdf5r | 1.3.16 |
| ggseqlogo | 0.2.2 |
