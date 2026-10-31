#Run from SMFNet_paper_analysis_code github repo: python ./figure_code/Fig3_E6_singleMotif_insilico/compute_single_motif_profiles.py

from pathlib import Path

INPUT_DIR = Path("input_files/Fig3_E6_singleMotif_insilico")
INTERMEDIATE_DIR = Path("output")
RESULTS_DIR = Path("output")
(INTERMEDIATE_DIR / "Fig3_E6_singleMotif_insilico").mkdir(parents=True, exist_ok=True)
import numpy as np
import h5py
import os
import pandas as pd
from collections import defaultdict
from numba import njit

H5_STATES = str(INPUT_DIR / "SMFNet-ONT_insilico_single_motif_example_states.h5")
out_dir = str(INTERMEDIATE_DIR / "Fig3_E6_singleMotif_insilico")
os.makedirs(out_dir, exist_ok=True)
N_BINS = 15

@njit(cache=True)
def extend_to_full(mat, smf_pos, window_center=512, full_len=1024):
    n = mat.shape[0]
    out = np.full((n, full_len), -1, dtype=mat.dtype)
    offset = window_center - full_len // 2
    for r in range(n):
        row = mat[r]
        out_row = out[r]
        n_valid = 0
        for i in range(row.shape[0]):
            if row[i] != -1:
                n_valid += 1
        if n_valid == 0:
            continue
        mapped_pos = np.empty(n_valid, dtype=np.int32)
        val = np.empty(n_valid, dtype=mat.dtype)
        j = 0
        for i in range(row.shape[0]):
            if row[i] != -1:
                mp = smf_pos[i] - offset
                if 0 <= mp < full_len:
                    mapped_pos[j] = mp
                    val[j] = row[i]
                    j += 1
        if j == 0:
            continue
        for i in range(j - 1):
            p1 = mapped_pos[i]
            p2 = mapped_pos[i + 1]
            if p2 <= p1:
                continue
            if val[i] == val[i + 1]:
                out_row[p1:p2] = val[i]
            else:
                mid = (p1 + p2) // 2
                out_row[p1:mid] = val[i]
                out_row[mid:p2] = val[i + 1]
        for i in range(j):
            out_row[mapped_pos[i]] = val[i]
    return out

def get_combo_from_region(region_name):
    parts = region_name.split("_", 1)
    return parts[1] if len(parts) == 2 else region_name

def mean_openness(X):
    valid = X != -1
    state3_fraction = np.divide(
        (X == 1).sum(axis=0),
        valid.sum(axis=0),
        out=np.full(X.shape[1], np.nan),
        where=valid.sum(axis=0) > 0,
    )
    return 1.0 - state3_fraction

with h5py.File(H5_STATES, "r") as f:
    all_regions = list(f["/smf_pos"].keys())
combo_to_regions = defaultdict(list)
for region in all_regions:
    combo_to_regions[get_combo_from_region(region)].append(region)
bin_rows = []
average_rows = []
for combo, regions in sorted(combo_to_regions.items()):
    print(f"Processing {combo}: {len(regions)} regions")
    X_list = []
    with h5py.File(H5_STATES, "r") as f:
        for region in regions:
            states = f["/smf_states/"][region][:]
            pos = f["/smf_pos/"][region][:]
            X = np.where(states == -1, -1, (states == 3).astype(np.int8))
            X = extend_to_full(X, pos)
            valid = X != -1
            row_mean = np.divide(
                (X == 1).sum(axis=1),
                valid.sum(axis=1),
                out=np.zeros(X.shape[0], dtype=float),
                where=valid.sum(axis=1) > 0,
            )
            keep_rows = (row_mean >= 0.05) & (row_mean <= 0.95)
            X_list.append(X[keep_rows])
    X_combo = np.vstack(X_list)
    n_reads, L = X_combo.shape
    print(f"  filtered reads: {n_reads}")
    avg_profile = mean_openness(X_combo)
    for pos_idx, value in enumerate(avg_profile):
        average_rows.append([combo, pos_idx + 1, value])
    valid = X_combo != -1
    state_score = np.divide(
        (X_combo == 1).sum(axis=1),
        valid.sum(axis=1),
        out=np.zeros(n_reads, dtype=float),
        where=valid.sum(axis=1) > 0,
    )
    X_sorted = X_combo[np.argsort(state_score)]
    bin_edges = np.linspace(0, n_reads, N_BINS + 1).astype(int)
    for b in range(N_BINS):
        i0, i1 = (bin_edges[b], bin_edges[b + 1])
        if i1 - i0 < 5:
            profile = np.full(L, np.nan)
        else:
            profile = mean_openness(X_sorted[i0:i1])
        pseudo_bin = N_BINS - 1 - b
        for pos_idx, value in enumerate(profile):
            bin_rows.append([combo, pseudo_bin, pos_idx + 1, value])
bin_df = pd.DataFrame(bin_rows, columns=["combo", "pseudo_bin", "position", "openness"])
average_df = pd.DataFrame(average_rows, columns=["combo", "position", "openness"])
bin_df.to_csv(
    os.path.join(out_dir, "example_three_motif_spacing_50bp_binned_profiles.tsv"),
    sep="\t",
    index=False,
)
average_df.to_csv(
    os.path.join(out_dir, "example_three_motif_spacing_50bp_average_profiles.tsv"),
    sep="\t",
    index=False,
)
print("DONE")
