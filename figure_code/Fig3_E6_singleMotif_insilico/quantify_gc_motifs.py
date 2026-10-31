#Run from SMFNet_paper_analysis_code github repo: python ./figure_code/Fig3_E6_singleMotif_insilico/quantify_gc_motifs.py

from pathlib import Path

INPUT_DIR = Path("input_files/Fig3_E6_singleMotif_insilico")
INTERMEDIATE_DIR = Path("output")
RESULTS_DIR = Path("output")
(INTERMEDIATE_DIR / "Fig3_E6_singleMotif_insilico").mkdir(parents=True, exist_ok=True)
import numpy as np
import pandas as pd
import h5py
from tqdm import tqdm
from numba import njit

H5_STATES = str(INPUT_DIR / "SMFNet-ONT_insilico_gc_test_states.h5")
OUT = str(INTERMEDIATE_DIR / "Fig3_E6_singleMotif_insilico/gc_single_motif_ac_summary.txt")

@njit(cache=True)
def extend_to_full(mat, smf_pos, window_center=512, full_len=1024):
    n = mat.shape[0]
    out = np.full((n, full_len), -1)
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
                if 1 <= mp <= full_len:
                    mapped_pos[j] = mp
                    val[j] = row[i]
                    j += 1
        if j == 0:
            continue
        out_row[: mapped_pos[0]] = val[0]
        for i in range(j - 1):
            p1, p2 = (mapped_pos[i], mapped_pos[i + 1])
            if p2 <= p1:
                continue
            if val[i] == val[i + 1]:
                out_row[p1:p2] = val[i]
            else:
                mid = (p1 + p2) // 2
                out_row[p1:mid] = val[i]
                out_row[mid:p2] = val[i + 1]
        out_row[mapped_pos[j - 1] :] = val[j - 1]
    return out

def stretch_around_x(vals, cols, x=512):
    vals = np.asarray(vals)
    cols = np.asarray(cols)

    bounds = np.r_[-np.inf, (cols[:-1] + cols[1:]) / 2, np.inf]
    idx = np.clip(np.searchsorted(bounds, x, side="right") - 1, 0, len(vals) - 1)

    if vals[idx] == 3:
        return np.array([np.nan, np.nan, 0])

    l = idx
    while l > 0 and vals[l - 1] != 3:
        l -= 1

    r = idx
    while r < len(vals) - 1 and vals[r + 1] != 3:
        r += 1

    left_edge = cols[0] - (cols[1] - cols[0]) / 2 if l == 0 else (cols[l - 1] + cols[l]) / 2
    right_edge = (
        cols[r] + (cols[r] - cols[r - 1]) / 2 if r == len(cols) - 1 else (cols[r] + cols[r + 1]) / 2
    )

    return np.array([left_edge, right_edge, right_edge - left_edge])

def boundary_euclidean_median(res):
    if len(res) < 2:
        return np.nan

    bounds = res[:, :2]
    diff = bounds[:, None, :] - bounds[None, :, :]
    dist = np.sqrt(np.sum(diff**2, axis=2))

    d = dist[np.triu_indices(len(bounds), k=1)]
    return np.median(d)

def extract_ac_info(h5, region, region_len=1024):
    smf_pos = h5["/smf_pos"][region][:]
    states_mat = h5["/smf_states"][region][:]
    center = region_len / 2
    binary_mat = (states_mat == 3).astype(int)
    extended_mat = extend_to_full(binary_mat, smf_pos, full_len=region_len)
    row_mean = np.mean(extended_mat, axis=1)
    keep_rows_nuc = (row_mean >= 0.05) & (row_mean <= 0.95)
    states_mat = states_mat[keep_rows_nuc]
    extended_mat = extended_mat[keep_rows_nuc]
    if states_mat.shape[0] == 0:
        return None
    center_idx = np.where(np.abs(np.arange(region_len) - center) <= 7.5)[0]
    keep_rows = np.sum(extended_mat[:, center_idx] == 1, axis=1) == 0
    perc_ac = np.mean(keep_rows)
    res = np.array(
        [stretch_around_x(states_mat[i], smf_pos, center) for i in np.where(keep_rows)[0]]
    )
    med_width = np.nanmedian(res[:, 2]) if len(res) else np.nan
    boundary_dispersion = boundary_euclidean_median(res)
    return {
        "region": region,
        "med_width": med_width,
        "perc_ac": perc_ac,
        "boundary_dispersion": boundary_dispersion,
    }

results = []
with h5py.File(H5_STATES, "r") as f:
    all_regions = list(f["/smf_pos"].keys())
    print(f"Processing {len(all_regions)} regions")
    for reg in tqdm(all_regions):
        try:
            res = extract_ac_info(f, reg)
            if res is not None:
                results.append(res)
        except Exception as e:
            print("Skipping", reg, ":", e)
pd.DataFrame(results).to_csv(OUT, sep="\t", index=False)
print("Wrote:", OUT)
