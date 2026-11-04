# PRISM calling for the Phc3 mCHIRA example

This example runs per-read in silico mutagenesis for wildtype `CRE_6005` from the Phc3 mCHIRA experiment. It follows the current [SMFNet PRISM notebook](https://github.com/kundajelab/SMFNet/blob/main/examples/prism.ipynb), using `load_model`, `densify_smf`, `prism` and `plot_prism`. It retains the calling settings and H5 result structure used for the paper's mCHIRA analysis. Clustering and cluster-average plots remain in the existing figure workflow.

Place this folder at `Processing_code/prism/`. Run from the SMFNet_paper_analysis_code repository root. Paths and settings are at the top of `call_PRISM.py`; there are no command-line arguments.

## Inputs

| Input | Where to obtain it |
|-------|--------------------|
| `input_files/prism/Phc3_CRE6005_PRISM_2048bp.fasta` | Small FASTA included with this example. Contains the exact 2,048-bp DNA window used by the original PRISM run. |
| `input_files/Fig2_E5_mCHIRA_clusters/mCHIRA_CRE6005_SMFNet_reads.h5` | Extract from `mCHIRA_CRE6005_data.tar.gz` or `mCHIRA_sampled_reads_and_clusters.tar.gz` on [Zenodo](https://zenodo.org/uploads/23407172). Both contain the same variable-region sampled reads. |
| `trained_models/mESC_Sonmezer_SMFNet_model` | Matching mESC Sonmezer short-read model, in SMFNet `save_pretrained` format. Set `MODEL` to this directory, a compatible checkpoint, or its eventual Hugging Face model ID. |

This run used the **Sonmezer short-read model**, not the ONT model in the upstream notebook. The final matching Hub ID is not yet confirmed; see the [SMFNet collection](https://huggingface.co/collections/kundajelab/smfnet). The historical local checkpoint was `smfnet_sonmezer_fold0_converted.ckpt`. A compatible `save_pretrained` export of that checkpoint has been prepared locally under `trained_models/mESC_Sonmezer_SMFNet_model/`, and reloading it through `load_model` preserves every model weight exactly. That export has not been uploaded by this step. Use the matching current-format Sonmezer model for external reproduction; the old Zenodo checkpoint may require conversion.

After downloading the archive:

```bash
mkdir -p input_files/Fig2_E5_mCHIRA_clusters
tar -xzf "$HOME/Downloads/mCHIRA_CRE6005_data.tar.gz" \
  -C input_files/Fig2_E5_mCHIRA_clusters
```

Use the dedicated PRISM FASTA supplied here. Its window is verified against the DNA stored in the original PRISM results. The separate read-sampling FASTA uses a different window and must not be substituted. The small FASTA's origin and checksum are recorded in `input_files/prism/source_inventory.tsv`.

## Run

Install SMFNet with its plotting dependencies following the [upstream instructions](https://github.com/kundajelab/SMFNet). This example was adapted against commit `6071b452f6bd5cced35341aab81eaa4b997ba1e5` and uses the public Python API directly:

```python
from smf_net.model import load_model
from smf_net.data import densify_smf
from smf_net.evaluation.prism import prism
from smf_net.evaluation.prism_plotting import plot_prism
```

Then run:

```bash
python Processing_code/prism/call_PRISM.py
```

The default scores all 1,024 wildtype reads and plots the first 10. Set `READ_INDICES = [0, 1]` for a small example; these are zero-based indices in the supplied read matrix. `PLOT_FIRST_N` controls plotting only. CUDA is the default; CPU execution is possible with `DEVICE = "cpu"` but is slower. No jobs are submitted automatically.

| Setting | Value |
|---------|-------|
| Region | `CRE_6005` |
| DNA input | 2,048 bp |
| Model output | Central 1,024 bp |
| Mutation interval, DNA coordinates | `[859, 1111)`; 252 bp |
| Same interval, output coordinates | `[347, 599)` |
| Mutation batch size | 64 |
| Unmasking orders | 8 |
| Order seed | 13 |
| Sampling mode | `one_pos_at_a_time` |
| Sample steps | 32 |
| Likelihood window | Full output (`None`) |
| Context mode | `full_read` |

The supplied read positions are local to the 252-bp variable region. The script adds 859 and calls `densify_smf` to place them in the DNA input coordinates. All other positions remain `-1`. Only these measured calls contribute to the score; no experimental or model calls are invented outside that window. With `one_pos_at_a_time`, trajectory length follows the valid calls rather than a fixed 32-step schedule; the original `sample_steps=32` argument is retained.

PRISM scores each read along seeded unmasking trajectories for the reference DNA and each single-base substitution. `deltas` are mutation minus reference log scores: positive values mean a substitution makes this read more likely, negative values mean it makes the read less likely. These are per-read trajectory scores, not changes in average accessibility. The upstream plotting helper displays the read and mutation effects.

## Outputs and existing results

`output/prism/CRE_6005_prism_results.h5` contains one `read_<index>` group per selected molecule. Each group preserves the original datasets:

- `scores`, `deltas`: mutation position × A/C/G/T.
- `reference_score`.
- `per_order_scores`, `per_order_reference_scores`.
- `dna`, `read`: the 252-bp mutation-window sequences and calls.

Mutation/output coordinates and the original read index are stored as group attributes. File attributes record input paths, model, seed and calling settings. The H5 is compressed and renamed from a temporary file on completion; existing final outputs are not overwritten. Plots are written as `CRE_6005_read_<index>.pdf` using the current upstream `plot_prism` style.

Existing complete results are available as `mCHIRA_CRE6005_PRISM.h5` in `mCHIRA_CRE6005_data.tar.gz` and `mCHIRA_WT_PRISM_complete.tar.gz`. They cover all 1,024 wildtype reads. They are not recomputed by downloading the inputs. The generated `read_<index>` groups are compatible with the existing clustering/plotting readers; configure their PRISM input path to use a newly generated result.

## Dependencies

Use Python 3.12 and the SMFNet repository's pixi environment/lockfile for the current API. Upstream pins PyTorch 2.5.1 and pytorch-lightning 2.6.0. The plotting extra includes matplotlib and tangermeme; this differs from some older figure environments. h5py and NumPy are used to read and write the supplied data. The upstream default and CPU pixi environments both include plotting. Open `pixi shell` in the SMFNet checkout, then change to the analysis repository root before running the example.

## Source versions

- Upstream API/example: `examples/prism.ipynb` at the commit above.
- Original calling script: `nate_model/scripts/20260918_run_prism_rmce.py`.
- Reference outputs: `nate_model/output/20260915_prism_cleanup/rmce_clusters/prism/CRE_6005_prism_results.h5`.

The example uses the current PRISM implementation, which has changed internally since the historical run. Matching inputs and settings alone does not establish bitwise equivalence across versions.

## Validation

The dedicated FASTA and prepared reads match the stored DNA/read windows for all 1,024 original wildtype PRISM results. The matching Sonmezer weights were exported and reloaded through the current `load_model` API with exact weight equality. A CPU check scored one read at one mutation position using all eight unmasking orders and the full original read likelihood: the reference score differed from the archive by 0.000210 and mutation deltas by at most 0.000253. This checks a bounded example, not all 252 positions or all molecules. The current upstream `plot_prism` also generated a PDF successfully from a complete archived read result. The full PRISM run has not been repeated.

The CPU check used Python 3.12.14, PyTorch 2.14.0, NumPy 2.2.6, h5py 3.16.0, pytorch-lightning 2.6.6, huggingface_hub 2.4.0 and safetensors 0.8.0. These are the validation environment versions, not replacements for the upstream lockfile. Temporary dependency installs were isolated from the existing environment.
