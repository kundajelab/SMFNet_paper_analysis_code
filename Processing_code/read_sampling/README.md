# SMFNet read sampling and in silico library generation

Run from the SMFNet_paper_analysis_code repository root. Place these scripts in `Processing_code/read_sampling/`. Edit paths and settings at the top of each script. Generated files go to `output/read_sampling/`, with directories created automatically.

## Inputs

The following files are available on zenodo.

| Archive | Relevant contents | Use |
|---------|-------------------|-----|
| `insilico_motif_CWMs.tar.gz` | `insertion_cwms.meme`, `selected_cwms.meme` | ChromBPNet CWMs used for background filtering and insertion. |
| `Combined_motif_analysis_inputs.tar.gz` | `background_noMotif.fasta` | Original 200 motif-depleted backgrounds. |
| `SMFNet-ONT_insilico_gc_inputs.tar.gz` | `gc_backgrounds/*.fasta` and GC-screen states | Original 200 backgrounds for each GC target from 30 to 70%. |
| `Combined_motif_original_FASTAs.tar.gz` | `fasta/*.fasta` | Original triplet and flank-control libraries. |
| `CTCF_PRISM_results.tar.gz` | `CTCF_example_sequence.fa` | The CTCF example sequence at chr1:191287950–191289998. |
| `Example_Sox2_locus_data.tar.gz` | `Fig1_E1_ONT_benchmarking/Example_SOX2_locus_regions.bed` | The paper's tiled Sox2 locus example: 21 windows of 2,048 bp. |
| `Endogenous_ISM_inputs.tar.gz` | `mESC_endogenous_regions.fa` | 11,682 active test-region sequences, each 2,048 bp, retained for the ISM scatter analysis versue ChromBPnet. |
| `F1_sampled_reads.tar.gz` | `F1_parental_sequences.fasta` | 37,590 BL6/SPRET sequences centered on ATAC peaks, each 2,048 bp. |
| `mCHIRA_CRE6005_data.tar.gz` | `mCHIRA_CRE6005_sequences.fasta` | 16 CRE6005 library sequences, each 2,048 bp. Used to sample the mCHIRA sequences. |
| `mESC_ac_testRegions_2048bp.bed` | 11,687 | Active test windows used for SMFNet benchmarking. |

`input_files/read_sampling/insertion_cwms.meme` is included locally in the GitHub repository.

To reuse the original backgrounds, put `background_noMotif.fasta` and the `gc_backgrounds/` directory under `output/read_sampling/`. Alternatively regenerate them with the two background scripts below. These are the paths the library generators read by default.

**Model:** download the DNMT-TKO mESC SMF-ONT model from the [SMFNet Hugging Face collection](https://huggingface.co/collections/kundajelab/smfnet). The model filename placeholder is `trained_models/DNMT-TKO_mESC_SMF-ONT_SMFNet_model.ckpt`. Set `CHECKPOINT` in `sample_reads.py` to the downloaded file under `trained_models/`.

### Sampled reads available on Zenodo

These archives are listed in the saved Zenodo upload inventory. The table identifies model-sampled reads, rather than experimental reads, state calls or summary metrics. Some archives contain both sampled and experimental data; the sampled-read members are named explicitly below.

| Archive | Sampled-read files inside | Dataset / extraction destination |
|---------|---------------------------|---------------------------------|
| `Example_Sox2_locus_data.tar.gz` | `Fig1_E1_ONT_benchmarking/Example_SOX2_locus_SMFNet_reads.h5` | Reads sampled for Sox2 locus example. |
| `mESC_SMF-ONT_testRegions_reads.tar.gz` | `mESC_SMF-ONT_testRegions_SMFNet_reads.h5` | mESC ONT test-region predictions. Used for benchmarking. |
| `mESC_Sonmezer_testRegions_data.tar.gz` | `mESC_Sonmezer_testRegions_SMFNet_reads.h5` | Short-read mESC model predictions. Used for benchmarking |
| `Extra_models_and_Fiber-seq_inputs.tar.gz` | `GM12878_Fiber-seq_testRegions_samples.h5` | GM12878 Fiber-seq examples. Used to plot Fiber-seq example regions. |
| `mCHIRA_CRE6005_data.tar.gz` | `mCHIRA_CRE6005_SMFNet_reads.h5` | CRE6005 WT and mutant reads restricted to the 252-bp variable region. |
| `mCHIRA_sampled_reads_and_clusters.tar.gz` | `mCHIRA_CRE6005_SMFNet_untrimmed_reads.h5`, `mCHIRA_CRE6005_SMFNet_reads.h5` | Full sampled output windows and the corresponding 252-bp reads: 1,024 reads for each of 16 constructs. Also includes cluster annotations. |
| `OCT4_footprint_sampled_reads.tar.gz` | `footprint_samples/three_motif_spacing_*bp/samples.h5` | Sampled-read subsets for the OCT4 footprint analysis at spacings 20–100 bp in 5-bp steps. |
| `F1_sampled_reads.tar.gz` | `F1_parental_samples.h5`, `F1_projected_common_samples.h5` | F1 parental predictions and their projection onto common coordinates.|


## Generate backgrounds

```bash
python Processing_code/read_sampling/make_background_fasta.py
python Processing_code/read_sampling/make_gc_background_fastas.py
```

The first script makes 200 sequences of length 2,048 with equal initial base probabilities. The GC screen makes 200 sequences for each target GC percentage. GC percentages describe generation probabilities, finite sequences need not have exactly that percentage.

Both scripts scan both strands for matches scoring at least 0.85 of each CWM's maximum score. They merge overlapping or touching 40-bp windows around matches, shuffle bases within those windows, and repeat until no qualifying matches remain. Shuffling preserves each sequence's base counts. While this removes strong matches to accessibility-related TF motifs, this approach does not guarantee removal of all biological TF binding sites. 

## Generate libraries

Run only the libraries needed:

```bash
python Processing_code/read_sampling/make_single_motif_fasta.py
python Processing_code/read_sampling/make_gc_motif_fasta.py
python Processing_code/read_sampling/make_homotypic_fastas.py
python Processing_code/read_sampling/make_triplet_fastas.py
python Processing_code/read_sampling/make_flank_control_fastas.py
python Processing_code/read_sampling/make_pairwise_fastas.py
```

| Script | Library |
|--------|---------|
| `make_single_motif_fasta.py` | One central motif or a no-motif control; 4,200 sequences. |
| `make_gc_motif_fasta.py` | Central motifs and controls across 41 GC targets; 172,200 sequences. |
| `make_homotypic_fastas.py` | Corrected symmetric placement of 1–7 copies, at 20–80-bp spacing in 10-bp steps; 28,200 sequences per spacing. |
| `make_triplet_fastas.py` | Symmetric flanks `(A, B, A)`, single-position controls and empty backgrounds; 20–100-bp spacing in 5-bp steps; 92,200 sequences per spacing. |
| `make_flank_control_fastas.py` | Flank-only `(A, noMotif, A)` controls for the triplet screen; 4,000 sequences per spacing. |
| `make_mesc_enhancer_fastas.py` | All 216 three-position combinations of noMotif, OCT4::SOX2, ZIC, SP-KLF, SOX2 and ESRRB; 20–80-bp spacing in 10-bp steps; 43,200 sequences per spacing. |

Counts assume the supplied 20 insertion motifs and 200 backgrounds. Spacing means **center-to-center distance**. Pair examples have one-sided boundary geometry; the triplet screen has symmetric flanks. Homotypic placement retains the corrected rounding for even copy numbers. Motifs are inserted as forward consensus sequences from CWM argmax values. Original sequence IDs and ordering are retained.

## Sample reads

```bash
python Processing_code/read_sampling/sample_reads.py
```

Set `CHECKPOINT` to a `Path` with a trained SMFNet model, `SOURCE` to `"fasta"`, and `INPUT_FILES` to any list of FASTAs for all generated libraries. To sample genomic regions, set `SOURCE = "bed"`, list BED files in the same setting, and set `REFERENCE_FASTA` to a `Path`.

Defaults match the in silico setup: 512 reads per region, 512 reads per chunk, 32 regions per batch, 32 steps and seed 0, using CUDA and bfloat16 computation with float32 parameters. Genomic regions were sampled with 1024 reads per region.

SMFNet centers/crops or pads FASTA records to the checkpoint's input length. 

## Dependencies

## Python 3.12

| Package | Version |
|---------|---------|
| numpy | 2.2.6 |
| pandas | 3.0.5 |
| h5py | 3.16.0 |
| pyfaidx | 0.9.0.4 |
| pytorch-lightning | 2.6.6 |
| hydra-core | 1.3.7 |
| hydra-zen | 0.16.0 |
| tqdm | 4.70.1 |
| tangermeme | 0.5.1 |
| torch | 2.5.1+cu121 |
