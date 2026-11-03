#Run from SMFNet_paper_analysis_code github repo: python ./Processing_code/motif_annotation/compute_ChromBPNet_DeepSHAP.py
from pathlib import Path
import numpy as np
import pandas as pd
import torch
from pyfaidx import Fasta
from bpnetlite.bpnet import BPNet, CountWrapper
from bpnetlite.attribute import deep_lift_shap
INPUT = Path("input_files/motif_annotation")
OUTPUT = Path("output/motif_annotation")

MODEL_FILE = Path("trained_models/model_ChromBPNet_mESC_nobias.h5")
FASTA_FILE = INPUT / "mm10.fa"
REGIONS_FILE = INPUT / "mESC_annotation_regions.bed"
OUTPUT_DIR = OUTPUT / "shap"
DEVICE = "cuda"
REGION_LIMIT = None


def read_regions(bed_file, fasta_file, length=2114, limit=None):
    sequences, regions, skipped = [], [], []
    lookup = np.full(256, -1, dtype=np.int16)
    lookup[list(b"ACGT")] = np.arange(4)
    with Fasta(str(fasta_file)) as fasta:
        bed = pd.read_csv(bed_file, sep="\t", header=None, usecols=[0, 1, 2])
        for chrom, start, end in bed.itertuples(index=False, name=None):
            center = (int(start) + int(end)) // 2
            start = max(0, center - length // 2)
            end = start + length
            seq = fasta[chrom][start:end].seq.upper()
            region = f"{chrom}:{start}-{end}"
            if len(seq) != length or any(base not in "ACGT" for base in seq):
                skipped.append(region)
                continue
            bases = lookup[np.frombuffer(seq.encode("ascii"), dtype=np.uint8)]
            sequences.append(np.eye(4, dtype=np.float32)[bases].T)
            regions.append(region)
            if limit is not None and len(regions) >= limit:
                break
    if not sequences:
        raise ValueError("No complete A/C/G/T windows found")
    return torch.from_numpy(np.stack(sequences)), np.asarray(regions), skipped


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    sequences, regions, skipped = read_regions(REGIONS_FILE, FASTA_FILE, limit=REGION_LIMIT)
    model = CountWrapper(BPNet.from_chrombpnet(filename=str(MODEL_FILE)))
    #Keep the original count target, seed, references and attribution convention.
    scores = deep_lift_shap(model, sequences, random_state=0, n_shuffles=20, device=DEVICE)
    hypothetical = deep_lift_shap(model, sequences, random_state=0, n_shuffles=20,
                                 hypothetical=True, device=DEVICE)
    for name, array in [("ohe", sequences.numpy()), ("attr", scores.cpu().numpy()),
                        ("attr_hypo", hypothetical.cpu().numpy()), ("region_ids", regions)]:
        np.savez_compressed(OUTPUT_DIR / f"{name}.npz", array)
    #Fi-NeMo peak IDs must match the array rows after sequence filtering.
    with (OUTPUT_DIR / "region_ids.narrowPeak").open("w") as handle:
        for i, region in enumerate(regions):
            chrom, coordinates = region.rsplit(":", 1)
            start, end = map(int, coordinates.split("-"))
            handle.write(f"{chrom}\t{start}\t{end}\tpeak_{i}\t1000\t.\t0\t0\t0\t1057\n")
    (OUTPUT_DIR / "skipped_regions.txt").write_text("".join(region + "\n" for region in skipped))
    print(f"Wrote {len(regions)} regions; skipped {len(skipped)}")


if __name__ == "__main__":
    main()
