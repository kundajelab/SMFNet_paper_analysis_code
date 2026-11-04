# PRISM calling for the Phc3 mCHIRA example

This example runs per-read in silico mutagenesis for wildtype `CRE_6005` from the Phc3 mCHIRA experiment. It follows the current [SMFNet PRISM notebook](https://github.com/kundajelab/SMFNet/blob/main/examples/prism.ipynb), using `load_model`, `densify_smf`, `prism` and `plot_prism`. It retains the calling settings and H5 result structure used for the paper's mCHIRA/Phc3 promoter analysis.

Run from the SMFNet_paper_analysis_code repository root. Paths and settings are at the top of `call_PRISM.py`.

## Inputs

| Input | source |
|-------|--------------------|
| `input_files/prism/CTCF_example_sequence.fa` | CTCF example FASTA included with the GitHub repository. |
| `Phc3_CRE6005_PRISM_2048bp.fasta` | FASTA included with the GitHub repository. Contains the exact 2,048-bp DNA window used by the original PRISM run. |
| `mCHIRA_CRE6005_SMFNet_reads.h5` | Extract from `mCHIRA_CRE6005_data.tar.gz`. |
| `trained_models/mESC_Sonmezer_SMFNet_model` | Extract from `model_SMFNet_mESC_Sonmezer.tar.gz`. |

CTCF example used mESC SMF-ONT model supplied on huggingface. Phc3 promoter/mCHIR run used the **Sonmezer short-read model**, supplied via Zenodo, not the ONT model.

## Run

Install SMFNet with its plotting dependencies following the [upstream instructions](https://github.com/kundajelab/SMFNet). 

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

## mCHIRA/Phc3 promoter run

Set the input and output settings in `call_PRISM.py`:

```python
MODEL = "trained_models/model_SMFNet_mESC_Sonmezer.ckpt"
FASTA = Path("input_files/prism/Phc3_CRE6005_PRISM_2048bp.fasta")
READS_H5 = Path("input_files/prism/mCHIRA_CRE6005_SMFNet_reads.h5")
OUTPUT_DIR = Path("output/prism/Phc3_PRISM")
REGION = "CRE_6005"
READ_WINDOW_START = 859
READ_WINDOW_LENGTH = 252
MUTATION_START = 859
MUTATION_LENGTH = 252
BATCH_SIZE = 64
SAMPLE_STEPS = 32
N_UNMASKING_ORDERS = 8
SEED = 13
SAMPLING_MODE = "one_pos_at_a_time"
READ_INDICES = None
PLOT_FIRST_N = 10
```

After loading the model, DNA and sparse reads, the script performs this PRISM call for each selected molecule:

```python
result = prism(
    model, reads[index], dna,
    mutation_start=859,
    length=252,
    batch_size=64,
    sample_steps=32,
    n_unmasking_orders=8,
    seed=13,
    sampling_mode="one_pos_at_a_time",
    likelihood_window=None,
    context_mode="full_read",
    verbose=True,
)
```

Run:

```bash
python Processing_code/prism/call_PRISM.py
```

The mutation interval is `[859, 1111)` in DNA coordinates, corresponding to `[347, 599)` in the central 1,024-bp model output. `densify_smf` adds 859 to the supplied local read positions, all other calls remain `-1`. Only measured sites contribute to the score. With `one_pos_at_a_time`, trajectory length follows the valid calls rather than a fixed 32-step schedule.

Use the dedicated PRISM FASTA above, which matches the original result coordinates. Set `READ_INDICES = [0, 1]` for a smaller example. `PLOT_FIRST_N` controls plotting only. Output is `output/prism/Phc3_PRISM/CRE_6005_prism_results.h5`.



## CTCF example run

This example uses the mESC ONT model and the 2,048-bp sequence at `chr1:191287950-191289998`. Its FASTA identifier is `chr1_191287950_191289998`.

First sample reads with the current SMFNet sampling entry point. From the analysis repository root:

```bash
SMFNET_DIR=/path/to/SMFNet
python "$SMFNET_DIR/scripts/sample.py" \
  checkpoint=trained_models/DNMT-TKO_mESC_SMF-ONT_SMFNet_model.ckpt \
  regions.source=fasta \
  regions.fasta_path=input_files/prism/CTCF_example_sequence.fa \
  sampling.reads_per_region=1024 \
  sampling.reads_per_chunk=128 \
  sampling.regions_per_batch=1 \
  sampling.sample_steps=32 \
  sampling.seed=0 \
  hydra.run.dir=output/prism/CTCF_sampling
```

This generates a new set of molecules, so it will not replicate the molecules depicted in the SMFNet preprint exactly. Set the input and output settings in `call_PRISM.py`:

```python
MODEL = "trained_models/DNMT-TKO_mESC_SMF-ONT_SMFNet_model.ckpt"
FASTA = Path("input_files/prism/CTCF_example_sequence.fa")
READS_H5 = Path("output/prism/CTCF_sampling/samples.h5")
OUTPUT_DIR = Path("output/prism/CTCF_PRISM")
REGION = "chr1_191287950_191289998"
READ_WINDOW_START = 512
READ_WINDOW_LENGTH = 1024
MUTATION_START = 874
MUTATION_LENGTH = 300
BATCH_SIZE = 64
SAMPLE_STEPS = 300
N_UNMASKING_ORDERS = 1
SEED = 343
SAMPLING_MODE = "cosine"
READ_INDICES = [0, 1]
PLOT_FIRST_N = 10
```

After loading the model, DNA and sparse reads, the script performs this PRISM call for each selected molecule:

```python
result = prism(
    model, reads[index], dna,
    mutation_start=874,
    length=300,
    batch_size=64,
    sample_steps=300,
    n_unmasking_orders=1,
    seed=343,
    sampling_mode="cosine",
    likelihood_window=None,
    context_mode="full_read",
    verbose=True,
)
```

Run:

```bash
python Processing_code/prism/call_PRISM.py
```

The mutation interval is `[874, 1174)` in DNA coordinates, corresponding to `[362, 662)` in the central model output. `Set `READ_INDICES = None` to score all 1,024 molecules. Output is `output/prism/CTCF_PRISM/chr1_191287950_191289998_prism_results.h5`.

## Outputs and existing results

PRISM scores each molecule along seeded unmasking trajectories for the reference DNA and every single-base substitution. `deltas` are mutation-minus-reference log scores: positive values mean that a substitution makes this read more likely, and negative values mean it makes it less likely. They describe individual molecules, not changes in average accessibility.

Each output H5 contains one `read_<index>` group per selected molecule, with:

- `scores`, `deltas`: mutation position × A/C/G/T.
- `reference_score`.
- `per_order_scores`, `per_order_reference_scores`.
- `dna`, `read`: DNA and calls within the mutation interval.

Group attributes retain the mutation/output coordinates and original read index. File attributes record inputs, model and calling settings. The first `PLOT_FIRST_N` selected molecules also receive PDF plots through `plot_prism`.

Complete Phc3 results are supplied as `mCHIRA_CRE6005_PRISM.h5` in `mCHIRA_CRE6005_data.tar.gz` and `mCHIRA_WT_PRISM_complete.tar.gz`. Example CTCF scores are supplied in `CTCF_PRISM_results.tar.gz`. 


## Dependencies

## Python 3.12

| Package | Version |
|---------|---------|
| SMFNet | 0.1.0 |
| torch | 2.14.0 |
| pytorch-lightning | 2.6.6 |
| numpy | 2.2.6 |
| pandas | 3.0.5 |
| h5py | 3.16.0 |
| pyfaidx | 0.9.0.4 |
| scipy | 1.15.3 |
| matplotlib | 3.10.3 |
| tangermeme | 1.4.1 |
| huggingface_hub | 2.4.0 |
| safetensors | 0.8.0 |
| hydra-core | 1.3.7 |
| hydra-zen | 0.16.0 |
| tqdm | 4.70.1 |
