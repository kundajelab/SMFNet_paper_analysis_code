#Run from SMFNet_paper_analysis_code github repo: python ./Processing_code/prism/call_PRISM.py

from pathlib import Path

import h5py
import matplotlib
import numpy as np
import torch

from smf_net.data import densify_smf
from smf_net.evaluation.prism import prism
from smf_net.evaluation.prism_plotting import plot_prism
from smf_net.model import load_model

matplotlib.use("Agg")
import matplotlib.pyplot as plt

MODEL = "trained_models/mESC_Sonmezer_SMFNet_model"
FASTA = Path("input_files/prism/Phc3_CRE6005_PRISM_2048bp.fasta")
READS_H5 = Path("input_files/Fig2_E5_mCHIRA_clusters/mCHIRA_CRE6005_SMFNet_reads.h5")
OUTPUT_DIR = Path("output/prism")
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
DEVICE = "cuda"
READ_INDICES = None
PLOT_FIRST_N = 10


def load_inputs():
    sequences = {}
    name = None
    for line in FASTA.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            if name in sequences:
                raise ValueError(f"Duplicate FASTA ID: {name}")
            sequences[name] = ""
        elif line.strip():
            if name is None:
                raise ValueError("FASTA sequence is missing its identifier")
            sequences[name] += line.strip().upper()
    sequence = sequences[REGION]
    if len(sequence) != 2048 or set(sequence) - set("ACGT"):
        raise ValueError("Expected an unambiguous 2048-bp DNA sequence")
    dna = torch.tensor(["ACGT".index(base) for base in sequence], dtype=torch.long)

    with h5py.File(READS_H5, "r") as handle:
        if "region_ids" in handle:
            ids = list(handle["region_ids"].asstr()[:])
            if ids.count(REGION) != 1:
                raise ValueError(f"Expected one H5 region named {REGION}")
            group = handle[f"regions/{ids.index(REGION):08d}/predicted"]
            calls = group["smf"][:]
            positions = group["positions"][:]
        else:
            calls = handle["smf_mat"][REGION][:]
            positions = handle["smf_pos"][REGION][:]
    if calls.ndim != 2 or positions.ndim != 1:
        raise ValueError("Expected a read matrix and a site-position vector")
    if calls.shape[1] != len(positions) and calls.shape[0] == len(positions):
        calls = calls.T
    if calls.shape[1] != len(positions) or not len(positions):
        raise ValueError("Read matrix and positions do not match")
    if not np.isin(calls, [-1, 0, 1]).all():
        raise ValueError("Expected hard calls -1, 0 and 1")
    if (not np.isfinite(positions).all() or not np.equal(positions, positions.astype(int)).all()
            or np.any(np.diff(positions) <= 0) or positions.min() < 0 or positions.max() >= READ_WINDOW_LENGTH):
        raise ValueError("Expected local positions in the supplied 252-bp construct window")
    #Place the published variable-region calls in the original PRISM DNA coordinates.
    reads = densify_smf(
        torch.as_tensor(calls, dtype=torch.long),
        torch.as_tensor(positions, dtype=torch.long) + READ_WINDOW_START,
        region_len=len(dna),
    ).long()
    return dna, reads


def main():
    dna, reads = load_inputs()
    if not 512 <= MUTATION_START < MUTATION_START + MUTATION_LENGTH <= 1536:
        raise ValueError("Mutation interval must lie within the model output window")
    indices = list(range(len(reads))) if READ_INDICES is None else list(READ_INDICES)
    if not indices or len(set(indices)) != len(indices):
        raise ValueError("Select distinct read indices")
    if any(index < 0 or index >= len(reads) for index in indices):
        raise ValueError("A read index is outside the input matrix")
    output = OUTPUT_DIR / f"{REGION}_prism_results.h5"
    if output.exists():
        raise FileExistsError(output)
    model = load_model(MODEL).to(DEVICE).eval()
    if model.cfg.input_length != len(dna):
        raise ValueError("The selected model does not match the DNA input length")
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(".h5.tmp")

    with h5py.File(temporary, "w") as handle:
        handle.attrs.update(
            region_id=REGION, model=MODEL, fasta=str(FASTA), reads_h5=str(READS_H5),
            variable_start=MUTATION_START, variable_length=MUTATION_LENGTH,
            model_output_start=MUTATION_START - len(dna) // 4,
            sample_steps=SAMPLE_STEPS, n_unmasking_orders=N_UNMASKING_ORDERS,
            seed=SEED, sampling_mode=SAMPLING_MODE, context_mode="full_read",
        )
        for number, index in enumerate(indices):
            print(f"Read {index} ({number + 1}/{len(indices)})", flush=True)
            result = prism(
                model, reads[index], dna,
                mutation_start=MUTATION_START, length=MUTATION_LENGTH,
                batch_size=BATCH_SIZE, sample_steps=SAMPLE_STEPS,
                n_unmasking_orders=N_UNMASKING_ORDERS, seed=SEED,
                sampling_mode=SAMPLING_MODE, likelihood_window=None,
                context_mode="full_read", verbose=True,
            )
            group = handle.create_group(f"read_{index}")
            for name in ("scores", "deltas", "reference_score", "per_order_scores",
                         "per_order_reference_scores", "dna", "read"):
                values = getattr(result, name).detach().cpu().numpy()
                options = {"compression": "gzip", "compression_opts": 9, "shuffle": True} if values.ndim else {}
                group.create_dataset(name, data=values, **options)
            for name in ("mutation_start", "mutation_end", "output_start", "output_end",
                         "likelihood_start", "likelihood_end", "context_mode"):
                group.attrs[name] = getattr(result, name)
            group.attrs["read_index"] = index
            if number < PLOT_FIRST_N:
                figure = plot_prism(result, f"Phc3 {REGION} read {index}")
                figure.savefig(OUTPUT_DIR / f"{REGION}_read_{index}.pdf", bbox_inches="tight")
                plt.close(figure)
    temporary.replace(output)
    print(f"Wrote {output}")


if __name__ == "__main__":
    main()
