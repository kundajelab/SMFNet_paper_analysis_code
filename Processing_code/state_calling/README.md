# Single-molecule chromatin state calling

`call_states.py` converts experimental or SMFNet-sampled methylation calls into chromatin-state annotations along each molecule. It supports both H5 layouts used in this repository and runs on CPUs without SMFNet, PyTorch or a GPU. The standard calling parameters are unchanged.

Place this folder at `Processing_code/state_calling/` and run commands from the SMFNet_paper_analysis_code repository root.

## How the calling works

This is **molecule-wide, annotation-independent footprint classification**, extending the logic of bin-based single-molecule calling. Krebs et al. (2017) classified individual molecules using binarized methylation in four predefined collection bins around transcription start sites. Here, classification follows the methylation pattern across each supplied molecule/window rather than requiring promoter coordinates, TF motifs or predefined collection bins. The relationship is conceptual: this is a stretch-length classifier, not a reimplementation of the original four-bin classifier. See [Krebs et al., Molecular Cell (2017)](https://doi.org/10.1016/j.molcel.2017.06.027), *Genome-wide Single-Molecule Footprinting Reveals High RNA Polymerase II Turnover at Paused Promoters*.

The experimental rationale is that SMF labeling is performed under saturating methyltransferase conditions: exposed, assayable DNA is efficiently methylated, while protein-protected DNA remains comparatively unmethylated. This makes stretches of methylated and unmethylated calls informative about accessibility and protection on individual molecules and supports a simple rule-based classification. Saturation concerns labeling of accessible substrate; it does not imply that every position on a molecule should be methylated. Reliable interpretation still depends on labeling efficiency, methylation-call quality and the spacing of assayable sites. The caller does not estimate or correct those properties itself.

The algorithm operates separately on each molecule:

1. Identify runs of protected calls (`1`) and estimate their lengths in base pairs from the measured-site positions. Boundaries are placed halfway between neighboring measured sites.
2. Label protected stretches shorter than 70 bp as short footprints (`TF`) and stretches of at least 70 bp as long protected footprints (`nucleosome`). Short, uncertain stretches touching a missing call or a read edge are handled using the 50-bp occupied-border threshold.
3. Identify accessible stretches, allowing short protected footprints to lie within a broader opening while retaining their footprint labels. Accessible stretches shorter than 80 bp receive the `linker` label; stretches of at least 80 bp receive the `accessible` label. The accessible-border threshold is 20 bp.

Calling uses valid measured sites and their genomic spacing; unassayed sites are not treated as protected calls. At least three valid sites are needed to estimate stretches. Missing calls and read boundaries limit the evidence available for calling short stretches. The output remains at the supplied measured positions, rather than creating measurements at unassayed bases.

`TF` and `nucleosome` are operational footprint-size classes. A short footprint does not identify a particular TF, and a long protected stretch is not proof of a specific nucleosome. Motif annotation and other comparisons are downstream analyses. The same calling rules are applied to experiment and model samples.

## Run with the standard settings

Edit the settings at the top of `call_states.py`. The input is deliberately unset:

```python
H5_IN = Path("input_files/Fig1_E1_ONT_benchmarking/mESC_SMF-ONT_testRegions_SMFNet_reads.h5")
OUT_DIR = Path("output/state_calling")
OUT_STR = "mESC_SMF-ONT_testRegions_SMFNet"
WORKERS = 8
```

Then run without command-line arguments:

```bash
python Processing_code/state_calling/call_states.py
```

The script creates the output directory. Replace `H5_IN` and `OUT_STR` for another dataset. For newly sampled reads, the input can instead be `output/read_sampling/samples/<input_stem>/samples.h5` from the [general sampling example](../read_sampling/README.md).

| Setting | Default | Meaning |
|--------|---------|---------|
| `OC_BORDER_VAL` | `50` | Border rule for short protected stretches, in bp. |
| `AC_BORDER_VAL` | `20` | Border rule for short accessible stretches, in bp. |
| `OC_THRESH` | `70` | Length separating short and long protected footprints, in bp. |
| `AC_THRESH` | `80` | Length separating linker and accessible stretches, in bp. |
| `SUBGROUP` | `predicted` | Read subgroup in the indexed H5 layout. |
| `MATRIX_ORDER` | `auto` | Infer matrix orientation from the positions dimension. |
| `CHUNK_READS` | `2048` | Maximum molecules per processing chunk. |
| `BLOCK_MB` | `256` | Additional approximate chunk-memory limit. |
| `WORKERS` | Up to 8 Slurm CPUs; otherwise 1 | Numba worker threads; set explicitly for interactive use. |
| `MAX_STRETCHES` | `1000000` | Maximum sampled stretch lengths for diagnostics or automatic thresholds. |
| `SEED` | `0` | Seed for sampling stretch lengths. |
| `OUT_STR` | `sm` | Output filename prefix. |
| `MAKE_PLOTS` | `True` | Write diagnostic density PDFs. |
| `OVERWRITE` | `False` | Allow replacing an existing output H5. |

The fixed 50/20/70/80-bp settings are the standard workflow. `OC_THRESH = "auto"` and/or `AC_THRESH = "auto"` optionally estimate thresholds from stretch-length distributions and require additional passes over the data. They are not needed to use the standard settings. `MAKE_PLOTS = False` skips diagnostic density plots without changing fixed-threshold state calls. `OVERWRITE = True` explicitly permits replacing an existing output H5.

## Input and output formats

Input accessibility calls must be hard calls: `-1` missing/unassayed, `0` methylated/accessible, and `1` unmethylated/protected. Probability matrices and endogenous methylation channels are not suitable substitutes.

Supported inputs:

- Original layout: `smf_mat/<region>` and `smf_pos/<region>`.
- Indexed samples layout: `region_ids`, with `regions/00000000/predicted/smf` and `regions/00000000/predicted/positions` for each numbered region. Use `SUBGROUP = "observed"` for experimental reads if an `observed` subgroup is present.

The matrix shape is molecules × measured positions. Non-square transposed matrices can be detected automatically. For square matrices stored as positions × molecules, specify `MATRIX_ORDER = "positions-by-reads"` explicitly. Positions must be finite and strictly increasing.

The output is `output/state_calling/<prefix>_smf_states.h5`, containing `smf_states/<region>` and `smf_pos/<region>`. State matrices preserve molecule order and use these codes:

| Code | State |
|------|-------|
| `-1` | Missing or uncertain/trimmed call |
| `0` | Accessible |
| `1` | Linker |
| `2` | Short protected footprint (`TF`) |
| `3` | Long protected footprint (`nucleosome`) |

Read IDs are not copied into the state H5; retain the input file to recover molecule identity by row order. Thresholds, border settings, matrix order and seed are recorded as H5 attributes. State datasets use gzip compression. The output H5 is written to a temporary file and renamed on completion.

By default the script also writes `<prefix>_states_oc_density.pdf` and `<prefix>_states_ac_density.pdf`. The first run compiles Numba kernels; later runs reuse the cache. Set `NUMBA_CACHE_DIR` to a writable location if the script directory is read-only. The memory limit controls processing chunks, not total process memory or all diagnostic buffers.

## State H5 files on Zenodo

Download archives from the [Zenodo deposition](https://zenodo.org/uploads/23407172). The table lists states called from **model-sampled reads**, using the saved upload inventory and local archive manifests. The private draft was not accessible for a fresh online check. Experimental states packaged alongside them are listed separately below.

| Archive | Model-sampled state H5 files | Extract into |
|---------|-----------------------------|--------------|
| `Example_Sox2_locus_data.tar.gz` | `Fig1_E1_ONT_benchmarking/Example_SOX2_locus_SMFNet_states.h5` | `input_files/` (archive includes the figure subdirectory) |
| `mESC_SMF-ONT_testRegions_states.tar.gz` | `mESC_SMF-ONT_testRegions_SMFNet_states.h5` | `input_files/Fig1_E1_ONT_benchmarking/` |
| `ChromBPNet_SMFNet_marginalizations.tar.gz` | `SMFNet_marginalization_states.h5` | `input_files/Fig2_E4_ChromBPNet_benchmark/` |
| `mCHIRA_CRE6005_data.tar.gz` | `mCHIRA_CRE6005_SMFNet_states.h5` | `input_files/Fig2_E5_mCHIRA_clusters/` |
| `SMFNet-ONT_insilico_gc_inputs.tar.gz` | `SMFNet-ONT_insilico_gc_test_states.h5` | `input_files/Fig3_E6_singleMotif_insilico/` |
| `SMFNet-ONT_insilico_homotypic_states.tar.gz` | `SMFNet-ONT_insilico_homotypic_motifs_spacing_{spacing}bp_states.h5`; spacing = 20, 30, 40, 50, 60, 70, 80 | `input_files/Fig3_E6_singleMotif_insilico/` |
| `Three_motif_states.tar.gz` | `three_motif_spacing_{spacing}bp_states.h5`; spacing = 20–100 in 5-bp steps (17 files) | `input_files/Fig4_E7_combinedMotifs_insilico/` |
| `Three_motif_flank_controls.tar.gz` | `flank_only_spacing_{spacing}bp_states.h5`; spacing = 20–100 in 5-bp steps (17 files) | `input_files/Fig4_E7_combinedMotifs_insilico/` |
| `mESC_enhancer_example_states.tar.gz` | `mesc_enhancer_spacing_{spacing}bp_states.h5`; spacing = 20, 30, 40, 50, 60, 70, 80 | `input_files/Fig4_E7_combinedMotifs_insilico/` |
| `F1_states.tar.gz` | `F1_predicted_common_states.h5` | `input_files/Fig5_E8_F1_analysis/` |

The braces above describe filename families, not literal filenames. The confirmed archive mappings cover 54 model-sampled state H5 files. These are already called states: use them directly for downstream analysis rather than passing them through `call_states.py` again.

The same archives also supply experimental state files:

- `Example_Sox2_locus_data.tar.gz`: `Fig1_E1_ONT_benchmarking/Example_SOX2_locus_Experiment_states.h5`.
- `mESC_SMF-ONT_testRegions_states.tar.gz`: `mESC_SMF-ONT_testRegions_Experiment_states.h5`.
- `mCHIRA_CRE6005_data.tar.gz`: `mCHIRA_CRE6005_Experiment_states.h5`.
- `F1_states.tar.gz`: `F1_B6_states.h5` and `F1_SPRET_states.h5`.

`SMFNet-ONT_insilico_single_motif_states.tar.gz`, containing `SMFNet-ONT_insilico_single_motif_example_states.h5`, exists locally but is absent from the saved upload inventory. Its upload status is therefore unconfirmed. Likewise, the contents of the broader `results_*.tar.gz` archives have not been inventoried here, so no additional state H5 members are inferred from those archive names. See the [sampled-read table](../read_sampling/README.md#sampled-reads-available-on-zenodo) for the upstream read H5s.

For example, download `mESC_SMF-ONT_testRegions_states.tar.gz` into Downloads and run:

```bash
mkdir -p input_files/Fig1_E1_ONT_benchmarking
tar -xzf "$HOME/Downloads/mESC_SMF-ONT_testRegions_states.tar.gz" \
  -C input_files/Fig1_E1_ONT_benchmarking
```

Use the corresponding destination in the table for other archives. Supplied state files keep their published filenames; newly computed files use the `OUT_STR` prefix described above. The local mCHIRA archive was updated to replace its sequence FASTA; its state H5s were unchanged.

## Dependencies

Shared reproduction targets: Python 3.12. These package versions follow the analysis environment targets; this documentation update does not install or change the environment.

| Package | Version |
|---------|---------|
| numpy | 2.2.6 |
| h5py | 3.15.1 |
| numba | 0.62.1 |
| scipy | 1.15.3 |
| matplotlib | 3.10.7 |
| seaborn | 0.13.2 |

Fixed-threshold calling with `MAKE_PLOTS = False` needs only NumPy, h5py and Numba. Diagnostics and automatic thresholds additionally use the plotting/scientific packages above.

## Validation

The existing validation compares intermediate stretch matrices, stretch lengths and final states against the original implementation, including real H5 regions and automatic thresholds. The recorded CPU benchmark measured a 61.2-fold improvement on 5,000 × 400 matrices with two threads after kernel warm-up; this excludes compilation and H5 I/O and is not an end-to-end runtime claim. The configuration cleanup preserves the calling kernels and default thresholds. A matrix-dimension typo was corrected (`positions.size`); both H5 layouts and a transposed input were checked against the saved reference state outputs.
