#Run from the repository root: python figure_code/Fig1_E3_Codependency_benchmarking/compute_TFBS_coaccessibility.py
#Run prepare_TFBS_pairs.py first. Scores are calculated from existing state H5 files.
from pathlib import Path
import h5py
import numpy as np
import pandas as pd

input_dir = Path("input_files/Fig1_E3_Codependency_benchmarking")
output_dir = Path("output/Fig1_E3_Codependency_benchmarking/TFBS")
pair_file = output_dir / "TFBS_motif_pairs.tsv"
experiment_file = input_dir / "mESC_SMF-ONT_testRegions_Experiment_states.h5"
model_file = input_dir / "mESC_SMF-ONT_testRegions_SMFNet_states.h5"
MIN_OBSERVED_PAIR_READS = 100
MIN_FIRST_TFBS_BOUND_FRACTION = 0.1

def fast_corr(x, y):
    x = x - x.mean()
    y = y - y.mean()
    return np.sum(x * y) / np.sqrt(np.sum(x*x) * np.sum(y*y))

def fast_mi(x, y):
    xy = x * 2 + y
    counts = np.bincount(xy.astype(int), minlength=4).reshape(2,2)
    joint = counts / counts.sum()

    px = joint.sum(axis=1)
    py = joint.sum(axis=0)

    nz = joint > 0
    return np.sum(joint[nz] * np.log2(joint[nz] / (px[:, None] * py[None, :])[nz]))

def entropy(x):
    if len(x) == 0:
        return np.nan
    probs = np.bincount(x.astype(int)) / len(x)
    probs = probs[probs > 0]
    return -np.sum(probs * np.log2(probs))

def signed_mi_pair(x, y):
    valid = (x != -1) & (y != -1)
    xv = x[valid]
    yv = y[valid]

    if len(xv) == 0:
        return np.nan

    #binarize like before
    xv = np.where(xv <= 2, 1, 0)
    yv = np.where(yv <= 2, 1, 0)

    if np.std(xv) == 0 or np.std(yv) == 0:
        return np.nan

    Hx = entropy(xv)
    Hy = entropy(yv)

    mi = fast_mi(xv, yv)
    nmi = mi / np.sqrt(Hx * Hy) if Hx > 0 and Hy > 0 else 0

    corr = fast_corr(xv, yv)

    return nmi * np.sign(corr)


def main():
    pairs = pd.read_csv(pair_file, sep="\t")
    #Match the original pre-shift position filter; state coordinates can include -1.
    pairs = pairs.loc[(pairs["pos1"] >= -1) & (pairs["pos2"] >= -1)].copy()
    results = []
    with h5py.File(experiment_file, "r") as act_h5, h5py.File(model_file, "r") as pred_h5:
        for region, group in pairs.groupby("region", sort=True):
            if region not in act_h5["smf_states"] or region not in pred_h5["smf_states"]:
                continue
            matrices = {}
            for label, h5 in (("act", act_h5), ("pred", pred_h5)):
                positions = h5["smf_pos"][region][()].astype(int)
                matrix = h5["smf_states"][region][()]
                if matrix.shape[1] != len(positions):
                    raise ValueError(f"Unexpected state matrix orientation: {region}")
                matrices[label] = (matrix, {p: i for i, p in enumerate(positions)})
            for row in group.itertuples(index=False):
                #Reproduce the saved figure cohort using jointly observed experiment reads.
                #The first TFBS follows the original motif ordering, not genomic left/right order.
                if any(row.pos1 not in index or row.pos2 not in index for _, index in matrices.values()):
                    continue
                observed, index = matrices["act"]
                x = observed[:, index[row.pos1]]
                y = observed[:, index[row.pos2]]
                valid = (x != -1) & (y != -1)
                if valid.sum() < MIN_OBSERVED_PAIR_READS:
                    continue
                if np.mean(x[valid] == 2) < MIN_FIRST_TFBS_BOUND_FRACTION:
                    continue
                result = {key: getattr(row, key) for key in ("region", "motif1", "motif2", "pos1", "pos2", "distance")}
                #Entropy is calculated on jointly observed calls for each TFBS pair.
                for label, (matrix, positions) in matrices.items():
                    if row.pos1 not in positions or row.pos2 not in positions:
                        score = np.nan
                    else:
                        score = signed_mi_pair(matrix[:, positions[row.pos1]], matrix[:, positions[row.pos2]])
                    result[f"mi_{label}"] = score
                results.append(result)
    output_dir.mkdir(parents=True, exist_ok=True)
    output = output_dir / "TFBS_coaccessibility.tsv"
    pd.DataFrame(results, columns=["region", "motif1", "motif2", "pos1", "pos2", "distance", "mi_act", "mi_pred"]).to_csv(output, sep="\t", index=False)
    print(f"Wrote {len(results)} pairs to {output}")


if __name__ == "__main__":
    main()
