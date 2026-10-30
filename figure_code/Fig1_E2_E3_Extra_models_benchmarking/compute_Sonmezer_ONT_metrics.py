#!/usr/bin/env python3
#Run from the repository root: python figure_code/Fig1_E2_E3_Extra_models_benchmarking/compute_Sonmezer_ONT_metrics.py

from pathlib import Path
import h5py
import numpy as np
import pandas as pd
from tqdm import tqdm
from scipy.stats import pearsonr

INPUT_DIR = Path("input_files/Fig1_E2_E3_Extra_models_benchmarking")
OUTPUT_DIR = Path("output/Fig1_E2_E3_Extra_models_benchmarking")
#Use equal read counts for both comparisons and require coverage per site pair.
N_READS = 512
MIN_PAIR_READS = 50
SEED = 42


def region_index(h5):
    if h5.attrs["schema"] != "smfnet.samples":
        raise ValueError("Expected smfnet.samples")
    return {region: f"regions/{i:08d}" for i, region in enumerate(h5["region_ids"].asstr()[:])}


def read_region(h5, index, region):
    group = h5[index[region]]["predicted"]
    pos = group["positions"][:]
    mat = group["smf"][:]
    if mat.ndim != 2 or mat.shape[1] != len(pos):
        raise ValueError(f"{region}: incompatible matrix and positions")
    return mat, pos


#Match site coordinates before comparing the two models.
def align_regions(mat1, pos1, mat2, pos2):
    common, idx1, idx2 = np.intersect1d(pos1, pos2, return_indices=True)
    return mat1[:, idx1], mat2[:, idx2], common

def subsample(mat, n, rng):
    if mat.shape[0] < n:
        return None
    idx = rng.choice(mat.shape[0], size=n, replace=False)
    return mat[idx]

def pairwise_correlation(mat, min_reads=MIN_PAIR_READS):
    #Binary contingency counts exclude missing calls (-1) from each site pair.
    zero = (mat == 0).astype(np.float64)
    one = (mat == 1).astype(np.float64)

    n00 = zero.T @ zero
    n01 = zero.T @ one
    n10 = n01.T
    n11 = one.T @ one

    i, j = np.triu_indices(mat.shape[1], k=1)
    a, b = n00[i, j], n01[i, j]
    c, d = n10[i, j], n11[i, j]

    n = a + b + c + d
    #The phi coefficient equals Pearson correlation for binary calls.
    denominator = np.sqrt((a + b) * (c + d) * (a + c) * (b + d))

    corr = np.full(len(i), np.nan)
    valid = (n >= min_reads) & (denominator > 0)
    corr[valid] = (a[valid] * d[valid] - b[valid] * c[valid]) / denominator[valid]

    return i, j, corr

def correlation_correlation(mat1, mat2):
    i1, j1, corr1 = pairwise_correlation(mat1)
    i2, j2, corr2 = pairwise_correlation(mat2)

    if not (np.array_equal(i1, i2) and np.array_equal(j1, j2)):
        raise ValueError("Pairwise indices do not match")

    #Compare only site pairs with valid correlations in both read sets.
    valid = np.isfinite(corr1) & np.isfinite(corr2)
    n_pairs = int(valid.sum())

    if n_pairs < 2:
        return np.nan, n_pairs

    x, y = corr1[valid], corr2[valid]
    if np.std(x) == 0 or np.std(y) == 0:
        return np.nan, n_pairs

    return float(pearsonr(x, y).statistic), n_pairs

def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    results = []
    rng_model = np.random.default_rng(SEED)
    with h5py.File(INPUT_DIR / "mESC_Sonmezer_testRegions_SMFNet_reads.h5", "r") as hs, \
         h5py.File(INPUT_DIR / "mESC_SMF-ONT_testRegions_SMFNet_reads.h5", "r") as ho:
        son_index, ont_index = region_index(hs), region_index(ho)
        regions = sorted(son_index.keys() & ont_index.keys())
        for region in tqdm(regions, desc="Comparing regions"):
            son, son_pos = read_region(hs, son_index, region)
            ont, ont_pos = read_region(ho, ont_index, region)
            son, ont, pos = align_regions(son, son_pos, ont, ont_pos)
            if len(pos) < 2 or son.shape[0] < N_READS or ont.shape[0] < 2 * N_READS:
                continue
            son_sub = subsample(son, N_READS, rng_model)
            ont_sub = subsample(ont, N_READS, rng_model)
            #One disjoint ONT split, using the original first-split seed per region.
            split = np.random.default_rng(SEED).choice(ont.shape[0], 2 * N_READS, replace=False)
            #Compare models, then compare two non-overlapping ONT read subsets.
            comparisons = {
                "ont_sonmezer": (ont_sub, son_sub),
                "ont_ont": (ont[split[:N_READS]], ont[split[N_READS:]]),
            }
            for comparison, (first, second) in comparisons.items():
                correlation, n_pairs = correlation_correlation(first, second)
                results.append({
                    "region": region, "comparison": comparison,
                    "n_positions": len(pos), "n_sonmezer_reads": len(son),
                    "n_ont_reads": len(ont), "n_reads_per_comparison": N_READS,
                    "covariation_cor": correlation, "n_position_pairs": n_pairs,
                })
    #Write one row per region and comparison, including the number of valid pairs.
    columns = ["region", "comparison", "n_positions", "n_sonmezer_reads",
               "n_ont_reads", "n_reads_per_comparison", "covariation_cor", "n_position_pairs"]
    df = pd.DataFrame(results, columns=columns)
    path = OUTPUT_DIR / "Sonmezer_ONT_correlation_metrics.tsv"
    df.to_csv(path, sep="\t", index=False)
    print(df.groupby("comparison")["covariation_cor"].median())
    print("Wrote:", path)


if __name__ == "__main__":
    main()
