# Generative models of single-molecule footprints reveal principles of combinatorial transcription factor assembly on chromatin
Till Schwämmle#, Nathaniel Diamant#, Yuliya Badayeva, Rozemarijn Kleinendorst, Anshul Kundaje* and Arnaud R Krebs*

This repository contains code to reproduce computational analysis performed in Schwämmle, Diamant et al., bioRxiv 2026.

## Abstract
<p align="justify"> </p>

## Data availability
<p align="justify"> SMF-ONT data generated from F1 B6xSpret DNMT-TKO and 159 DNMT-TKO mESCs are available in GEO under GSE350193 and GSE350195.

Trained SMFNet models for DNMT-TKO mESC SMF-ONT and GM12878 Fiber-seq data are available on Hugging Face (https://huggingface.co/collections/kundajelab/smfnet). Additional models, processed training inputs, sampled reads and DNA sequences used for sampling are available on Zenodo (DOI: 10.5281/zenodo.23407172).
</p>

## Instructions
<p align="justify"> This repository contains the code needed to reproduce all manuscript figures from processed sequencing data and SMFNet-generated samples. The required input files and trained models are available on Zenodo and Hugging Face via the links in the manuscript’s Data Availability section. Before running an analysis, download the relevant files and place them in `/input_files/` and `/trained_models/`, respectively, preserving their original filenames.

The `/example_code/` directory provides examples for generating processed inputs. These examples require edits to file paths, dataset-specific parameters, and environment settings before they can be run.

The `figure_code/` directory contains analysis scripts for reproducing the manuscript figures. Follow the instructions provided in each analysis subdirectory.</p>
