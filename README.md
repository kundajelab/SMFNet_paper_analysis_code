# SMFNet reveals principles of combinatorial transcription factor assembly
Till Schwämmle#, Nathaniel Diamant#, Yuliya Badayeva, Rozemarijn Kleinendorst, Anshul Kundaje* and Arnaud R Krebs*

> [!WARNING]
> This repository is under construction. Code, documentation, and analysis instructions may be incomplete and are subject to change.

<p align="justify">This repository contains code to reproduce computational analysis performed in Schwämmle, Diamant et al. 2026, bioRxiv 2026 (https://doi.org/10.64898/2026.10.30.762954).</p>

## Abstract

<p align="justify">Combinatorial assembly of transcription factors (TFs) opens chromatin at cis-regulatory elements (cREs) to establish cell-type-specific gene expression. Dissecting the role of individual TFs and their collective function has been limited by the diversity of TF binding-site arrangements across cREs. A single regulatory sequence can support multiple configurations of TF and nucleosome occupancy, raising the possibility that distinct chromatin states reflect different TF assemblies. Predictive models linking regulatory sequence to chromatin accessibility collapse these distinct molecular states into population averages, obscuring how TF assemblies shape accessible chromatin. Here, we introduce SMFNet, a generative diffusion model that predicts distributions of single-molecule chromatin accessibility profiles from DNA sequence, and PRISM, which identifies sequence features associated with individual molecular states. Together, these approaches enable us to systematically  associate single molecule chromatin states at cREs with specific TF binding motifs. Experimental perturbation of these motifs redistributes the corresponding states as predicted by SMFNet, supporting their interpretation as TF-dependent configurations underlying chromatin opening. We then use controlled in silico motif screens to systematically dissect individual and combinatorial TF contributions to chromatin opening. Individual TFs show distinct chromatin-opening signatures, differing in accessible-region width and nucleosome-positioning precision even when associated with similar levels of chromatin accessibility. Combining motifs reveals partner- and spacing-dependent effects on chromatin accessibility, TF occupancy and nucleosome positioning. Genetic perturbation of TF binding sites using natural variation at endogenous cRE confirms that TF-specific properties shape chromatin accessibility in vivo. Together, our findings reveal sequence-encoded principles of combinatorial TF assembly during chromatin opening. </p>

## Data availability

<p align="justify">SMF-ONT data generated from F1 B6xSpret DNMT-TKO and 159 DNMT-TKO mESCs are available in GEO under GSE350193 and GSE350195.</p>

<p align="justify">Trained SMFNet models for DNMT-TKO mESC SMF-ONT and GM12878 Fiber-seq data are available on Hugging Face (<a href="https://huggingface.co/collections/kundajelab/smfnet">https://huggingface.co/collections/kundajelab/smfnet</a>). Additional models, processed training inputs, sampled reads and DNA sequences used for sampling are available on Zenodo (DOI: 10.5281/zenodo.23407172).</p>

## Instructions

<p align="justify">This repository contains the code needed to reproduce all manuscript figures from processed sequencing data and SMFNet-generated samples. The required input files and trained models are available on Zenodo and Hugging Face via the links in the manuscript’s Data Availability section. Before running an analysis, download the relevant files and place them in <code>/input_files/</code> and <code>/trained_models/</code>, respectively, preserving their original filenames.</p>

<p align="justify">The <code>/example_code/</code> directory provides examples for generating processed inputs. These examples require edits to file paths, dataset-specific parameters, and environment settings before they can be run.</p>

<p align="justify">The <code>/figure_code/</code> directory contains analysis scripts for reproducing the manuscript figures. Follow the instructions provided in each analysis subdirectory.</p>
