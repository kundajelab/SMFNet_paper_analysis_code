#Run from the repository root: python figure_code/Fig1_E3_Codependency_benchmarking/compute_CTCF_coaccessibility.py
#Step 2 of 3. Run count_CTCF_motif_states.py first
#Bound = state 2
from pathlib import Path
import os
import re
import time
import h5py
import numpy as np
import pandas as pd
from tqdm import tqdm

input_dir = Path("input_files/Fig1_E3_Codependency_benchmarking")
output_dir = Path("output/Fig1_E3_Codependency_benchmarking/CTCF")
MOTIF_FILE = output_dir / "CTCF_motif_state_counts.tsv"
ACT_STATE = input_dir / "mESC_SMF-ONT_testRegions_Experiment_states.h5"
PRED_STATE = input_dir / "mESC_SMF-ONT_testRegions_SMFNet_states.h5"
OUT_DIR = output_dir / "matrices"
OUT_SUMMARY = output_dir / "CTCF_coaccessibility_summary.tsv"


OUT_DIR.mkdir(parents=True, exist_ok=True)

MIN_ACT_BOUND_N = 50
MIN_COVERAGE_PER_POS = 50
MIN_MOLECULES_PER_SPLIT = 20

#Load and filter motifs

motif_df = pd.read_csv(MOTIF_FILE, sep="\t")

motif_df = motif_df[(motif_df["act_bound_n"] >= MIN_ACT_BOUND_N) & (motif_df["motif"] == "CTCF")].copy()

print(f"Motifs after act_bound_n >= {MIN_ACT_BOUND_N}: {len(motif_df)}")
print(f"Regions represented: {motif_df['region'].nunique()}")
if motif_df.empty:
    raise ValueError("No CTCF motifs pass the minimum observed bound-read count.")

#Helpers

def sanitize_name(x):
    x = str(x)
    x = x.replace(":", "_").replace("-", "_")
    x = re.sub(r"[^A-Za-z0-9_.]+", "_", x)
    return x

def entropy(x):
    """Shannon entropy of binary 1D array with NaNs."""
    x = x[~np.isnan(x)]
    if len(x) == 0:
        return np.nan

    x_int = x.astype(int)
    probs = np.bincount(x_int, minlength=2) / len(x_int)
    probs = probs[probs > 0]

    return -np.sum(probs * np.log2(probs))

def fast_corr(x, y):
    x = x - x.mean()
    y = y - y.mean()

    denom = np.sqrt(np.sum(x * x) * np.sum(y * y))

    if denom == 0:
        return np.nan

    return np.sum(x * y) / denom

def fast_mi(x, y):
    """
    Mutual information for binary arrays x/y.
    x and y must be valid non-NA binary arrays.
    """
    xy = x * 2 + y
    counts = np.bincount(xy.astype(int), minlength=4).reshape(2, 2)

    if counts.sum() == 0:
        return np.nan

    joint = counts / counts.sum()

    px = joint.sum(axis=1)
    py = joint.sum(axis=0)

    nz = joint > 0

    return np.sum(joint[nz] * np.log2(joint[nz] / (px[:, None] * py[None, :])[nz]))

def signed_mi(x, y, Hx, Hy):
    """Normalized mutual information multiplied by sign of correlation."""
    valid = ~np.isnan(x) & ~np.isnan(y)

    xv = x[valid]
    yv = y[valid]

    if len(xv) == 0:
        return np.nan

    if np.std(xv) == 0 or np.std(yv) == 0:
        return np.nan

    mi = fast_mi(xv, yv)

    if Hx > 0 and Hy > 0:
        nmi = mi / np.sqrt(Hx * Hy)
    else:
        nmi = 0

    corr = fast_corr(xv, yv)

    return nmi * np.sign(corr)

def compute_signed_mi_matrix(mat):
    """
    Compute signed MI between all columns.

    Original states:
      -1 = missing
       0/1/2 = accessible/bound/open-ish -> collapsed to 1
       3 = nucleosomal -> collapsed to 0
    """
    mat_collapsed = np.where(
        mat == -1,
        np.nan,
        np.where(mat <= 2, 1, 0)
    )

    n_pos = mat_collapsed.shape[1]

    mi_mat = np.full((n_pos, n_pos), np.nan, dtype=np.float32)

    entropies = np.array([
        entropy(mat_collapsed[:, i])
        for i in range(n_pos)
    ])

    for i in range(n_pos):
        xi = mat_collapsed[:, i]
        Hx = entropies[i]

        for j in range(i, n_pos):
            yj = mat_collapsed[:, j]
            Hy = entropies[j]

            val = signed_mi(xi, yj, Hx, Hy)

            mi_mat[i, j] = val
            mi_mat[j, i] = val

    return mi_mat

def save_mi_matrix(mi_mat, smf_pos, out_file):
    col_header = "\t".join(map(str, smf_pos))

    with open(out_file, "w") as f:
        f.write("\t" + col_header + "\n")

        for row_pos, row_vals in zip(smf_pos, mi_mat):
            row_str = "\t".join(
                "nan" if np.isnan(v) else f"{v:.4f}"
                for v in row_vals
            )
            f.write(f"{row_pos}\t{row_str}\n")

def orient_matrix(mat, pos, region, label):
    """Ensure matrix is molecules x positions."""
    if mat.ndim != 2:
        raise ValueError(f"{label} {region}: expected 2D matrix, got {mat.shape}")

    if mat.shape[1] == len(pos):
        return mat

    if mat.shape[0] == len(pos):
        return mat.T

    raise ValueError(
        f"{label} {region}: cannot match positions to matrix. "
        f"len(pos)={len(pos)}, mat.shape={mat.shape}"
    )

def subset_to_common_positions(mat_act, pos_act, mat_pred, pos_pred):
    """Restrict ACT and PRED to shared SMF positions."""
    common_pos = np.intersect1d(pos_act, pos_pred)

    if len(common_pos) == 0:
        return None, None, None

    act_idx = np.array([np.where(pos_act == p)[0][0] for p in common_pos])
    pred_idx = np.array([np.where(pos_pred == p)[0][0] for p in common_pos])

    mat_act = mat_act[:, act_idx]
    mat_pred = mat_pred[:, pred_idx]

    return mat_act, mat_pred, common_pos

def filter_positions_by_act_coverage(mat_act, mat_pred, pos):
    """Keep positions with sufficient ACT coverage."""
    coverage_per_pos = np.sum(mat_act != -1, axis=0)
    valid_pos_mask = coverage_per_pos >= MIN_COVERAGE_PER_POS

    mat_act = mat_act[:, valid_pos_mask]
    mat_pred = mat_pred[:, valid_pos_mask]
    pos = pos[valid_pos_mask]

    return mat_act, mat_pred, pos

def get_closest_position_index(pos, target_pos, max_dist=7.5):
    """Find closest SMF position to motif center."""
    closest_idx = int(np.argmin(np.abs(pos - target_pos)))
    closest_pos = float(pos[closest_idx])
    closest_dist = abs(closest_pos - target_pos)

    if closest_dist > max_dist:
        return None, closest_pos, closest_dist

    return closest_idx, closest_pos, closest_dist

def process_split(mat, pos, motif_idx, split_label, out_file):
    """Split molecules based on whether state at motif_idx is bound state 2."""
    motif_states = mat[:, motif_idx]

    if split_label == "bound":
        keep_rows = motif_states == 2
    elif split_label == "unbound":
        keep_rows = (motif_states != 2) & (motif_states != -1)
    else:
        raise ValueError("split_label must be 'bound' or 'unbound'")

    mat_split = mat[keep_rows, :]

    n_molecules = mat_split.shape[0]

    if n_molecules < MIN_MOLECULES_PER_SPLIT:
        return {
            "status": "too_few_molecules",
            "n_molecules": n_molecules,
            "n_positions": mat_split.shape[1]
        }

    mi_mat = compute_signed_mi_matrix(mat_split)
    save_mi_matrix(mi_mat, pos, out_file)

    return {
        "status": "ok",
        "n_molecules": n_molecules,
        "n_positions": mat_split.shape[1]
    }

#Main

summary_rows = []

with h5py.File(ACT_STATE, "r") as act_h5, h5py.File(PRED_STATE, "r") as pred_h5:

    act_states_group = act_h5["smf_states"]
    act_pos_group = act_h5["smf_pos"]

    pred_states_group = pred_h5["smf_states"]
    pred_pos_group = pred_h5["smf_pos"]

    for row_i, row in tqdm(
        motif_df.iterrows(),
        total=len(motif_df),
        desc="Processing motifs"
    ):
        region = row["region"]
        motif = row["motif"]
        rel_center = float(row["rel_center"])

        if region not in act_states_group:
            summary_rows.append({
                "region": region,
                "motif": motif,
                "rel_center": rel_center,
                "status": "missing_act_region"
            })
            continue

        if region not in pred_states_group:
            summary_rows.append({
                "region": region,
                "motif": motif,
                "rel_center": rel_center,
                "status": "missing_pred_region"
            })
            continue

        start_time = time.perf_counter()

        #Load ACT

        pos_act = np.asarray(act_pos_group[region][()], dtype=float)
        mat_act = np.asarray(act_states_group[region][()])

        mat_act = orient_matrix(mat_act, pos_act, region, "ACT")

        #Load PRED

        pos_pred = np.asarray(pred_pos_group[region][()], dtype=float)
        mat_pred = np.asarray(pred_states_group[region][()])

        mat_pred = orient_matrix(mat_pred, pos_pred, region, "PRED")

        #Match ACT/PRED positions

        mat_act, mat_pred, pos_common = subset_to_common_positions(
            mat_act,
            pos_act,
            mat_pred,
            pos_pred
        )

        if mat_act is None:
            summary_rows.append({
                "region": region,
                "motif": motif,
                "rel_center": rel_center,
                "status": "no_common_positions"
            })
            continue

        #Coverage filter based on ACT

        mat_act, mat_pred, pos_filtered = filter_positions_by_act_coverage(
            mat_act,
            mat_pred,
            pos_common
        )

        if mat_act.shape[1] < 2:
            summary_rows.append({
                "region": region,
                "motif": motif,
                "rel_center": rel_center,
                "status": "too_few_positions_after_coverage_filter",
                "n_positions": mat_act.shape[1]
            })
            continue

        #Find motif column after filtering

        motif_idx, motif_smf_pos, motif_smf_dist = get_closest_position_index(
            pos_filtered,
            rel_center,
            max_dist=7.5
        )

        if motif_idx is None:
            summary_rows.append({
                "region": region,
                "motif": motif,
                "rel_center": rel_center,
                "closest_smf_pos": motif_smf_pos,
                "closest_smf_dist": motif_smf_dist,
                "status": "motif_position_removed_or_too_far"
            })
            continue

        #Output filenames

        safe_region = sanitize_name(region)
        safe_motif = sanitize_name(motif)

        base = (
            f"{safe_region}"
            f"__motif_{safe_motif}"
            f"__center_{rel_center:.1f}"
            f"__smfpos_{motif_smf_pos:.1f}"
        )

        split_jobs = [
            ("ACT", "bound", mat_act),
            ("ACT", "unbound", mat_act),
            ("PRED", "bound", mat_pred),
            ("PRED", "unbound", mat_pred),
        ]

        for dataset, split_label, mat in split_jobs:
            out_file = os.path.join(
                OUT_DIR,
                f"{base}__{dataset}_{split_label}_mi_mat.txt"
            )

            res = process_split(
                mat=mat,
                pos=pos_filtered,
                motif_idx=motif_idx,
                split_label=split_label,
                out_file=out_file
            )

            summary_rows.append({
                "region": region,
                "motif": motif,
                "rel_center": rel_center,
                "motif_smf_pos": motif_smf_pos,
                "motif_smf_idx": motif_idx,
                "motif_smf_dist": motif_smf_dist,
                "dataset": dataset,
                "split": split_label,
                "status": res["status"],
                "n_molecules": res["n_molecules"],
                "n_positions": res["n_positions"],
                "out_file": out_file if res["status"] == "ok" else np.nan,
                "elapsed_sec": time.perf_counter() - start_time
            })

summary_df = pd.DataFrame(summary_rows)
summary_df.to_csv(OUT_SUMMARY, sep="\t", index=False)

print("Done")
print(f"Wrote summary: {OUT_SUMMARY}")
print(summary_df["status"].value_counts(dropna=False))
