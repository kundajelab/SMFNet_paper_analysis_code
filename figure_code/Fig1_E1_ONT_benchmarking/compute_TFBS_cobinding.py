#Run from the repository root: python figure_code/Fig1_E3_Codependency_benchmarking/compute_TFBS_cobinding.py
#Run prepare_TFBS_pairs.py first. Co-binding uses its complete pair table, not the co-accessibility subset.
#Bound = state 2. Fractions and odds ratios use reads observed at both TFBS.
from pathlib import Path
import h5py
import numpy as np
import pandas as pd

input_dir = Path("input_files/Fig1_E3_Codependency_benchmarking")
analysis_dir = Path("output/Fig1_E3_Codependency_benchmarking/TFBS")
output_dir = analysis_dir / "cobinding"
pair_file = analysis_dir / "TFBS_motif_pairs.tsv"
experiment_file = input_dir / "mESC_SMF-ONT_testRegions_Experiment_states.h5"
model_file = input_dir / "mESC_SMF-ONT_testRegions_SMFNet_states.h5"
MIN_PAIR_READS = 100
MIN_BOUND_FRACTION = 0.05
PSEUDOCOUNT = 0.5


def log2_odds_ratio(x, y):
    counts = np.bincount((x == 2).astype(int) * 2 + (y == 2).astype(int), minlength=4).reshape(2, 2)
    counts = counts + PSEUDOCOUNT
    return np.log2((counts[1, 1] * counts[0, 0]) / (counts[1, 0] * counts[0, 1]))


def main():
    pairs = pd.read_csv(pair_file, sep="\t")
    #Match the original pre-shift site filter used to calculate state fractions.
    pairs = pairs.loc[(pairs.pos1 >= -1) & (pairs.pos2 >= -1)]
    selected, scores = [], []
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
                if any(row.pos1 not in index or row.pos2 not in index for _, index in matrices.values()):
                    continue
                stats = {key: getattr(row, key) for key in ("region", "motif1", "motif2", "pos1", "pos2", "distance")}
                vectors = {}
                for label, (matrix, index) in matrices.items():
                    x, y = matrix[:, index[row.pos1]], matrix[:, index[row.pos2]]
                    valid = (x != -1) & (y != -1)
                    x, y = x[valid], y[valid]
                    vectors[label] = (x, y)
                    stats[f"n_{label}"] = len(x)
                    stats[f"{label}_p2_pos1"] = np.mean(x == 2) if len(x) else np.nan
                    stats[f"{label}_p2_pos2"] = np.mean(y == 2) if len(y) else np.nan
                if stats["n_act"] < MIN_PAIR_READS or stats["n_pred"] < MIN_PAIR_READS:
                    continue
                fractions = [stats[f"{label}_p2_pos{site}"] for label in ("act", "pred") for site in (1, 2)]
                stats["min_bound"] = min(fractions)
                if stats["min_bound"] < MIN_BOUND_FRACTION:
                    continue
                selected.append(stats.copy())
                #The original odds-ratio calculation excludes negative state coordinates.
                if row.pos1 < 0 or row.pos2 < 0:
                    continue
                for label, (x, y) in vectors.items():
                    stats[f"log_odds_{label}"] = log2_odds_ratio(x, y)
                #Preserve the original complete-MI plotting cohort without calculating unused MI scores.
                stats["plot_eligible"] = all(0 < fraction < 1 for fraction in fractions)
                scores.append(stats)
    output_dir.mkdir(parents=True, exist_ok=True)
    if not scores:
        raise ValueError("No TFBS pairs pass the co-binding filters.")
    pd.DataFrame(selected).to_csv(output_dir / "TFBS_cobinding_pairs.tsv", sep="\t", index=False)
    pd.DataFrame(scores).to_csv(output_dir / "TFBS_cobinding.tsv", sep="\t", index=False)
    print(f"Selected {len(selected)} pairs; wrote {len(scores)} odds-ratio rows.")


if __name__ == "__main__":
    main()
