# Overview

plot_Sox2_example.R: Plots Sox2 locus examples for Fig. 1c-d
plot_SMF-ONT_qc.R: Plots experimental QC for mESCH DNMT-TKO SMF-ONT data for Fig. E1b-d
plot_SMFNet_pos_comp.R: Plots comparison between SMFNet and Experiment occupancy at single positions for Fig. 1e
plot_SMFNet_metrics.R: Plots SMFNet metrics for mESC DNMT-TKO SMF-ONT model for Fig. E1e,i-j
compute_SMFNet_motif_avg.py, count_SMFNet_motif_states.py, plot_SMFNet_motif_states.R: Calculates state percentages at putative annotated TFBSs. Calculates average enrichment around annotated motifs. Plots composite plots and motif state correlations for Fig. E1f-h, k

# Instructions

plot_Sox2_example.R: 
- Retrieve Example_Sox2_locus_data.tar.gz from Zenodo
- Unpack data in /input_files/Fig1_E1_ONT_benchmarking/
- Run with Rscript ./figure_code/Fig1_E1_ONT_benchmarking/plot_sox2_example.R

plot_SMF-ONT_qc.R: 
- Retrieve WT_SMF-ONT_qc_data.tar.gz from Zenodo
- Unpack data in /input_files/Fig1_E1_ONT_benchmarking/
- Run with Rscript ./figure_code/Fig1_E1_ONT_benchmarking/plot_SMF-ONT_qc.R

plot_SMFNet_pos_comp.R:
- Retrieve SMFNet_mESC_basic_qc.tar.gz from Zenodo
- Unpack data in /input_files/Fig1_E1_ONT_benchmarking/
- Run with Rscript ./figure_code/Fig1_E1_ONT_benchmarking/plot_SMF-ONT_qc.R

plot_SMFNet_metrics.R:
- Retrieve SMFNet_mESC_basic_qc.tar.gz from Zenodo
- Unpack data in /input_files/Fig1_E1_ONT_benchmarking/
- Run with Rscript ./figure_code/Fig1_E1_ONT_benchmarking/plot_SMF-ONT_qc.R

plot_SMFNet_motif_states.R:
- Retrieve SMFNet_mESC_basic_qc.tar.gz, mESC_SMF-ONT_testRegions_reads.tar.gz, mESC_SMF-ONT_testRegions_states.tar.gz from Zenodo
- Unpack data in /input_files/Fig1_E1_ONT_benchmarking/
- Run with python figure_code/Fig1_E1_ONT_benchmarking/count_SMFNet_motif_states.py
- Run with python figure_code/Fig1_E1_ONT_benchmarking/compute_SMFNet_motif_avg.py
- Run with Rscript ./figure_code/Fig1_E1_ONT_benchmarking/plot_SMFNet_motif_states.R



grep "modisco" /g/krebs/schwammle/smf_model/ben_models/* 
how to check the files here recursive?
Use -r to search recursively through all files and subdirectories:

grep -rn "modisco" /g/krebs/schwammle/smf_model/ben_models/

-n includes line numbers. To show only filenames containing a match, use -rl instead.


Convert to PDF
library(ggplot2)
library(egg)
library(gridExtra)
library(hdf5r)

#Run from github repo: Rscript SMFNet_paper_analysis_code/figure_code/Fig1_E1_ONT_benchmarking/plot_sox2_example.R
input_dir <- file.path("input_files", "Fig1_E1_ONT_benchmarking")
output_dir <- file.path("output", "Fig1_E1_ONT_benchmarking")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
avg_dir <- file.path(output_dir, "avg_plots")
sm_dir <- file.path(output_dir, "sm_plots")
dir.create(avg_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(sm_dir, recursive = TRUE, showWarnings = FALSE)

where does the file path go?
These paths are relative to your current working directory, not the script’s location. Check it in R with:

getwd()

For example, if you run Rscript from /path/to/repo, your code uses:

/path/to/repo/input_files/Fig1_E1_ONT_benchmarking/
/path/to/repo/output/Fig1_E1_ONT_benchmarking/avg_plots/
/path/to/repo/output/Fig1_E1_ONT_benchmarking/sm_plots/
file.path() joins the folder names; dir.create() creates the output directories. This snippet doesn’t save any files yet.

To use a specific location regardless of where you run the script, set absolute paths:

input_dir <- "/your/path/input_files/Fig1_E1_ONT_benchmarking"
output_dir <- "/your/path/output/Fig1_E1_ONT_benchmarking"


Convert to PDF
# Dependencies
Python 3.12
R 4.4

plot_Sox2_example.R: 
- ggplot2 v4.0.3
- egg v0.4.5
- gridExtra v2.3.1 
- hdf5r v1.3.16

plot_SMF-ONT_qc.R:
- ggplot2 v4.0.3
- egg v0.4.5
- gridExtra v2.3.1
- data.table v1.18.6.1
- viridis v0.6.5
- ggrastr v1.0.2

plot_SMFNet_pos_comp.R:
- ggplot2 v4.0.3
- dplyr v1.2.1
- tidyr v1.3.2
- egg v0.4.5
- gridExtra v2.3.1

plot_SMFNet_metrics.R:
| Package | Version |
|---|---|
| numpy | 2.2.6 |
| pandas | 3.0.5 |
| h5py | 3.16.0 |
| tqdm | 4.70.1 |


| Package | Version |
|---|---|
| ggplot2 | 4.0.3 |
| viridis | 0.6.5 |
| egg | 0.4.5 |




this is for README.md Can you just change this to Python: show packages below and then R with version and packages bleow? independent of specific script. Also format nicely
…
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

Convert to PDF






