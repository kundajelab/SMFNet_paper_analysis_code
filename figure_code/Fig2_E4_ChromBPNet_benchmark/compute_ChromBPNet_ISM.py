#Run from the repository root: python figure_code/Fig2_E4_ChromBPNet_benchmark/compute_ChromBPNet_ISM.py
#Calculate count-head ISM for endogenous regions and plot the central 1024 bp.
from pathlib import Path

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import torch
from bpnetlite import BPNet
from bpnetlite.bpnet import CountWrapper
from pyfaidx import Fasta
from tangermeme.plot import plot_logo

INPUT_DIR = Path("input_files/Fig2_E4_ChromBPNet_benchmark")
OUTPUT_DIR = Path("output/Fig2_E4_ChromBPNet_benchmark")
FASTA_FILE = INPUT_DIR / "mESC_endogenous_regions.fa"
MODEL_FILE = INPUT_DIR / "model_ChromBPNet_mESC_nobias.h5"
OUTPUT_H5 = OUTPUT_DIR / "ChromBPNet_endogenous_ISM.h5"
BATCH_SIZE = 256
OUTPUT_LENGTH = 1024


def one_hot_encode(sequence, sequence_length=2048):
    """Encode A/C/G/T in the original 2048-bp input window."""
    bases = {base: i for i, base in enumerate("ACGT")}
    encoded = np.zeros((sequence_length, 4), dtype=np.float32)
    for i, base in enumerate(sequence.upper()[:sequence_length]):
        if base in bases:
            encoded[i, bases[base]] = 1
    return encoded


def read_sequences(path):
    #Preserve FASTA ordering and exclude sequences with ambiguous bases.
    with Fasta(str(path)) as fasta:
        for region in fasta.keys():
            sequence = str(fasta[region][:]).upper()
            if set(sequence).issubset(set("ACGT")):
                yield region, torch.tensor(one_hot_encode(sequence)).T


def chrombpnet_count_ism(model, sequence, batch_size=256):
    """Mean reference-minus-mutant count score over the three alternate bases."""
    model.eval()
    with torch.no_grad():
        reference = model(sequence).squeeze().item()

    length = sequence.shape[-1]
    reference_bases = sequence[0].argmax(dim=0).cpu().numpy()
    mutants = []
    for position in range(length):
        for base in range(4):
            if base == reference_bases[position]:
                continue
            mutant = sequence.clone()
            mutant[0, :, position] = 0
            mutant[0, base, position] = 1
            mutants.append(mutant)
    mutants = torch.cat(mutants, dim=0)

    with torch.no_grad():
        predictions = [model(mutants[i:i + batch_size]).cpu()
                       for i in range(0, len(mutants), batch_size)]
    delta = reference - torch.cat(predictions).squeeze().numpy()
    scores = np.zeros(length, dtype=np.float32)
    for position in range(length):
        scores[position] = np.mean(delta[position * 3:position * 3 + 3])
    return scores


def plot_region(region, scores, dna, plot_dir):
    center = len(scores) // 2
    start, end = center - OUTPUT_LENGTH // 2, center + OUTPUT_LENGTH // 2
    logo = (dna[start:end] * scores[start:end, None]).T
    fig, ax = plt.subplots(figsize=(12, 3))
    plot_logo(logo, ax=ax)
    ax.set_title(f"{region} — ChromBPNet ISM")
    fig.savefig(plot_dir / f"{region}_chrombpnet_ism.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    plot_dir = OUTPUT_DIR / "ChromBPNet_ISM_logos"
    plot_dir.mkdir(exist_ok=True)
    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"Loading ChromBPNet on {device}", flush=True)
    model = CountWrapper(BPNet.from_chrombpnet(filename=str(MODEL_FILE))).to(device)

    #Write each region immediately so the full collection need not stay in memory.
    with h5py.File(OUTPUT_H5, "w") as output:
        score_group = output.create_group("ism_sum")
        dna_group = output.create_group("one_hot_dna")
        for i, (region, sequence) in enumerate(read_sequences(FASTA_FILE), 1):
            print(f"[{i}] ISM {region}", flush=True)
            scores = chrombpnet_count_ism(model, sequence.unsqueeze(0).to(device), BATCH_SIZE)
            key = region.replace("/", "_").replace(":", "_")
            dna = sequence.T.numpy()
            score_group.create_dataset(key, data=scores)
            dna_group.create_dataset(key, data=dna)
            plot_region(region, scores, dna, plot_dir)
    print(f"Wrote {OUTPUT_H5}", flush=True)


if __name__ == "__main__":
    main()
