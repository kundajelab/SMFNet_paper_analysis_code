# SMFNet read sampling and in silico libraries

Run from the SMFNet_paper_analysis_code repository root. Place these scripts in `Processing_code/read_sampling/`. Edit paths and settings at the top of each script; there are no command-line arguments. Generated files go to `output/read_sampling/`, with directories created automatically.

## Inputs

The Zenodo deposition is https://zenodo.org/uploads/23407172. The following archive names were checked against the saved upload inventory; the draft was not accessible for a fresh online check.

| Archive | Relevant contents | Use |
|---------|-------------------|-----|
| `insilico_motif_CWMs.tar.gz` | `insertion_cwms.meme`, `selected_cwms.meme` | Use **insertion_cwms.meme** for background filtering and insertion. The selected CWMs are for logos. |
| `Combined_motif_analysis_inputs.tar.gz` | `background_noMotif.fasta` | Original 200 motif-depleted backgrounds; other archive members are not needed here. |
| `SMFNet-ONT_insilico_gc_inputs.tar.gz` | `gc_backgrounds/*.fasta` and GC-screen states | Original 200 backgrounds for each GC target from 30 to 70%. |
| `Combined_motif_original_FASTAs.tar.gz` | `fasta/*.fasta` | Original triplet, flank-control and mESC enhancer libraries; can be sampled directly. |

`input_files/read_sampling/insertion_cwms.meme` is included locally as a small GitHub-suitable input and also supplied in the CWM archive. Extract that archive into `input_files/read_sampling/` if needed.

To reuse the original backgrounds, put `background_noMotif.fasta` and the `gc_backgrounds/` directory under `output/read_sampling/`. Alternatively regenerate them with the two background scripts below. These are the paths the library generators read by default. An archived library can instead be listed directly in `INPUT_FILES` in the sampler.

**Model:** download the DNMT-TKO mESC SMF-ONT model from the [SMFNet Hugging Face collection](https://huggingface.co/collections/kundajelab/smfnet). The collection is currently private and the final model filename is not yet confirmed. The model filename placeholder is `trained_models/DNMT-TKO_mESC_SMF-ONT_SMFNet_model.ckpt`. Set `CHECKPOINT` in `sample_reads.py` to the downloaded file under `trained_models/`. The sampler expects an SMFNet-format checkpoint, removes compiled-model parameter prefixes and loads all weights strictly. Compatibility with the private Hugging Face artifact has not yet been checked.

BED sampling also needs the matching reference genome and a BED of the desired regions. `input_files/read_sampling/mESC_acRegions.bed` supplies all 71,324 original active bait regions in mm10 as three-column, 0-based half-open BED intervals. It is an unchanged copy of the active-region input (2,113 bp per interval), not restricted to the test chromosomes. The sampler recenters these intervals to the model input length. Its provenance and checksum are recorded in `input_files/read_sampling/source_inventory.tsv`. The synthetic FASTA examples need neither a genome nor training reads.

## Generate backgrounds

```bash
python Processing_code/read_sampling/make_background_fasta.py
python Processing_code/read_sampling/make_gc_background_fastas.py
```

The first script makes 200 sequences of length 2,048 with equal initial base probabilities and seed 33. The GC screen makes 200 sequences for each target GC percentage, using initial sequence seed `33 + GC` and one shared shuffle RNG seeded 33 across the GC screen. GC percentages describe generation probabilities; finite sequences need not have exactly that percentage.

Both scripts scan both strands for matches scoring at least 0.85 of each CWM's maximum score. They merge overlapping or touching 40-bp windows around matches, shuffle bases within those windows, and repeat until no qualifying matches remain. Shuffling preserves each sequence's base counts. This is CWM-score filtering, not a statistical motif-significance threshold, and it does not guarantee removal of all biological TF binding sites. It does not fix CG/GC dinucleotides. Filtering can take time; supplied backgrounds avoid rerunning it.

## Generate libraries

Run only the libraries needed:

```bash
python Processing_code/read_sampling/make_single_motif_fasta.py
python Processing_code/read_sampling/make_gc_motif_fasta.py
python Processing_code/read_sampling/make_homotypic_fastas.py
python Processing_code/read_sampling/make_triplet_fastas.py
python Processing_code/read_sampling/make_flank_control_fastas.py
python Processing_code/read_sampling/make_pairwise_fastas.py
python Processing_code/read_sampling/make_enhancer_fastas.py
```

| Script | Library |
|--------|---------|
| `make_single_motif_fasta.py` | One central motif or a no-motif control; 4,200 sequences. |
| `make_gc_motif_fasta.py` | Central motifs and controls across 41 GC targets; 172,200 sequences. |
| `make_homotypic_fastas.py` | Corrected symmetric placement of 1–7 copies, at 20–80-bp spacing in 10-bp steps; 28,200 sequences per spacing. |
| `make_triplet_fastas.py` | Symmetric flanks `(A, B, A)`, single-position controls and empty backgrounds; 20–100-bp spacing in 5-bp steps; 92,200 sequences per spacing. |
| `make_flank_control_fastas.py` | Flank-only `(A, noMotif, A)` controls for the triplet screen; 4,000 sequences per spacing. |
| `make_pairwise_fastas.py` | All ordered central/side pairs including no-motif controls, separately for left and right sides; 176,400 sequences per spacing. |
| `make_enhancer_fastas.py` | All 216 three-position combinations of noMotif, OCT4::SOX2, ZIC, SP-KLF, SOX2 and ESRRB; 20–80-bp spacing in 10-bp steps; 43,200 sequences per spacing. |

Counts assume the supplied 20 insertion motifs and 200 backgrounds. Spacing means **center-to-center distance**. Pair examples have one-sided boundary geometry; the triplet screen has symmetric flanks. Homotypic placement retains the corrected rounding for even copy numbers. Motifs are inserted as forward consensus sequences from CWM argmax values. Original sequence IDs and ordering are retained; scripts differ by library design, not by individual output file.

## Sample reads

Install SMFNet with the `smf_net.model.SMFNet` and `smf_net.evaluation` sampling APIs. This example was checked against the local September 14 SMFNet checkout (Git HEAD `6233d674f83e6f28e0a5a5bacb77a4f79e57896f`; local changes may also be present).

```bash
python Processing_code/read_sampling/sample_reads.py
```

The script has no preset model, input type, input files or reference genome. Set `CHECKPOINT` to a `Path`, `SOURCE` to `"fasta"` or `"bed"`, and `INPUT_FILES` to any list of FASTAs, or to `sorted(Path("output/read_sampling/fasta").glob("*.fasta"))` for all generated libraries. To sample genomic regions, set `SOURCE = "bed"`, list BED files in the same setting, and set `REFERENCE_FASTA` to a `Path`. For the supplied active regions, use `INPUT_FILES = [Path("input_files/read_sampling/mESC_acRegions.bed")]` and the mm10 reference FASTA. Unconfigured scripts stop with an explanatory error before loading a model or creating outputs. No separate sampling script is needed for each file.

Defaults match the in silico setup: 512 reads per region, 512 reads per chunk, 32 regions per batch, 32 steps and seed 0, using CUDA and bfloat16 computation with float32 parameters. Seed 0 restarts for each input file. Reduce chunk/batch sizes if memory is insufficient; changing batching, precision or library order can change stochastic samples. For CPU use, set `DEVICE = "cpu"` and `COMPUTE_DTYPE = torch.float32`.

SMFNet centers/crops or pads FASTA records to the checkpoint's input length. BED coordinates are 0-based and half-open; intervals of a different width are recentered to the model's input length, with the start clamped at zero. BED column 6 controls reverse complementation. BED-derived region IDs use the resulting genomic coordinates, not column 4 names. Prefer full-length windows away from chromosome edges and unique input regions.

Each input writes `output/read_sampling/samples/<input_stem>/samples.h5`, using the new `smfnet.samples` layout: `region_ids`, numbered `regions`, and `predicted/smf` plus `predicted/positions`. Calls are `-1` missing, `0` methylated/accessible, `1` unmethylated/protected. The model determines output length and assayable positions. H5 files are compressed and streamed in batches; only completed files receive the final name. Existing final files are not overwritten. Distinct input stems are required.

State calling is separate; use the general state caller on these H5 files before running figure analyses. Zenodo also lists `SMFNet-ONT_insilico_single_motif_states.tar.gz`, `SMFNet-ONT_insilico_homotypic_states.tar.gz`, `Three_motif_states.tar.gz`, `Three_motif_flank_controls.tar.gz` and `mESC_enhancer_example_states.tar.gz` for working from existing results.

## Dependencies

Sampling smoke-test environment: Python 3.12.14.

| Package | Version |
|---------|---------|
| torch | 2.14.0 |
| numpy | 2.2.6 |
| pandas | 3.0.5 |
| h5py | 3.16.0 |
| pyfaidx | 0.9.0.4 |
| pytorch-lightning | 2.6.6 |
| hydra-core | 1.3.7 |
| hydra-zen | 0.16.0 |
| tqdm | 4.70.1 |

Library generation additionally uses tangermeme. Original-generator comparisons used tangermeme 0.5.1, torch 2.5.1+cu121 and NumPy 1.26.4 in Python 3.10.19. The shared Python 3.12 environment with NumPy 2.2.6 is suitable for the sampling API; byte-for-byte background regeneration across dependency versions has not been established. Reuse the archived backgrounds when exact historical sequences are needed.

## Validation

Representative one-background outputs from all seven library generators were compared byte-for-byte against their source scripts, including one representative spacing and GC level. Both background generators were compared using one sequence and the original seeds. Real converted-checkpoint CPU smoke tests generated two reads each from FASTA and reverse-strand BED inputs in the new H5 layout, using two sampling steps. The CWM, analysis-input and original-FASTA archive checksums matched the saved Zenodo upload inventory. Full libraries and GPU sampling were not rerun. `source_inventory.tsv` records the source scripts used for this cleanup.
