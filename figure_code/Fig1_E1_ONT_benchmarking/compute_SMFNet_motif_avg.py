#Run from SMFNet_paper_analysis_code github repo: python ./figure_code/Fig1_E1_ONT_benchmarking/compute_SMFNet_motif_avg.py
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from tqdm import tqdm

INPUT_DIR = Path("input_files/Fig1_E1_ONT_benchmarking")
OUTPUT_DIR = Path("output/Fig1_E1_ONT_benchmarking")
WINDOW = 256
REL_GRID = np.arange(-WINDOW, WINDOW + 1)


def get_region_dense_profile(h5, region, full_region=True, region_index=None):
    if region_index is None:
        region_index = {rid: f"regions/{i:08d}" for i, rid in enumerate(h5["region_ids"].asstr()[:])}
    if region not in region_index:
        return None
    role = "observed" if full_region else "predicted"
    group = h5[region_index[region]][role]
    pos = np.asarray(group["positions"][:], dtype=int)
    mat = group["smf"][:]
    if len(pos) == 0:
        return None
    if mat.ndim != 2 or mat.shape[1] != len(pos):
        raise ValueError(f"{region}: expected molecules by positions, got {mat.shape}")

    #Preserve the original coordinates and experimental crop, without a one-base correction.
    if full_region:
        keep = (pos >= 512) & (pos <= 1535)
        if not keep.any():
            return None
        pos = pos[keep] - 512
        mat = mat[:, keep]
    avg = collapse_matrix_to_position_average(mat)
    return sparse_to_dense_step_function(pos, avg, full_start=0, full_end=1023)


def collapse_matrix_to_position_average(mat):
    mat = np.asarray(mat, dtype=float)
    mat[mat == -1] = np.nan

    with np.errstate(invalid="ignore"):
        avg = np.nanmean(mat, axis=0)

    return avg


def sparse_to_dense_step_function(pos, avg, full_start=0, full_end=2047):
    pos = np.asarray(pos, dtype=int)
    avg = np.asarray(avg, dtype=float)

    keep = ~np.isnan(pos)
    pos = pos[keep]
    avg = avg[keep]

    if len(pos) == 0:
        return np.full(full_end - full_start + 1, np.nan)

    order = np.argsort(pos)
    pos = pos[order]
    avg = avg[order]

    full_pos = np.arange(full_start, full_end + 1)
    dense = np.full(len(full_pos), np.nan)

    if len(pos) == 1:
        dense[:] = avg[0]
        return dense

    midpoints = np.floor((pos[:-1] + pos[1:]) / 2).astype(int)

    dense[full_pos <= midpoints[0]] = avg[0]

    for i in range(1, len(pos) - 1):
        dense[(full_pos > midpoints[i - 1]) & (full_pos <= midpoints[i])] = avg[i]

    dense[full_pos > midpoints[-1]] = avg[-1]

    return dense


def extract_motif_window(dense_profile, motif_center, window=256):
    motif_center = int(round(motif_center))

    out = np.full(window * 2 + 1, np.nan)

    start = motif_center - window
    end = motif_center + window

    src_start = max(start, 0)
    src_end = min(end, len(dense_profile) - 1)

    if src_start > src_end:
        return out

    dst_start = src_start - start
    dst_end = dst_start + (src_end - src_start)

    out[dst_start:dst_end + 1] = dense_profile[src_start:src_end + 1]

    return out


def build_composites_for_h5(h5_path, motif_df, dataset_label, full_region):
    rows = []

    with h5py.File(h5_path, "r") as h5:
        region_index = {rid: f"regions/{i:08d}" for i, rid in enumerate(h5["region_ids"].asstr()[:])}
        region_cache = {}

        for _, row in tqdm(
            motif_df.iterrows(),
            total=len(motif_df),
            desc=f"Extracting {dataset_label}"
        ):
            region = row["region"]
            motif_clean = row["motif_clean"]
            rel_center = row["rel_center"]

            if region not in region_cache:
                region_cache[region] = get_region_dense_profile(
                    h5,
                    region,
                    full_region=full_region,
                    region_index=region_index
                )

            dense = region_cache[region]

            if dense is None:
                continue

            window_values = extract_motif_window(
                dense_profile=dense,
                motif_center=rel_center,
                window=WINDOW
            )

            tmp = pd.DataFrame({
                "dataset": dataset_label,
                "region": region,
                "motif": motif_clean,
                "rel_pos": REL_GRID,
                "value": window_values
            })

            rows.append(tmp)

    if len(rows) == 0:
        return pd.DataFrame(columns=["dataset", "region", "motif", "rel_pos", "value"])

    return pd.concat(rows, ignore_index=True)


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    motifs = pd.read_csv(OUTPUT_DIR / "SMFNet_test_regions_motifs.tsv", sep="\t")
    motifs["motif"] = motifs["motif"].astype(str)
    motifs = motifs[~motifs["motif"].str.contains(r"-alt2$|-alt$|double", case=False, regex=True)].copy()
    motifs["motif_clean"] = motifs["motif"]
    motifs["rel_center"] = pd.to_numeric(motifs["rel_center"], errors="coerce")
    motifs = motifs.dropna(subset=["region", "motif_clean", "rel_center"]).copy()

    #Average reads excluding -1, extend to midpoint boundaries, then extract +/-256 bp.
    observed = build_composites_for_h5(
        INPUT_DIR / "mESC_SMF-ONT_testRegions_Experiment_reads.h5", motifs, "act", True
    )
    predicted = build_composites_for_h5(
        INPUT_DIR / "mESC_SMF-ONT_testRegions_SMFNet_reads.h5", motifs, "pred", False
    )
    profiles = pd.concat([observed, predicted], ignore_index=True)
    profiles.to_csv(OUTPUT_DIR / "SMFNet_motif_composites_long.tsv", sep="\t", index=False)

    #Weight motif instances equally; retain boundary padding and the original strand orientation.
    means = profiles.groupby(["dataset", "motif", "rel_pos"], as_index=False).agg(
        mean_value=("value", "mean"),
        n_regions=("region", "nunique"),
        n_values=("value", "count")
    )
    means.to_csv(OUTPUT_DIR / "SMFNet_motif_composites_mean.tsv", sep="\t", index=False)
    print(f"Motif instances: {len(motifs)}; profile rows: {len(profiles)}; mean rows: {len(means)}")


if __name__ == "__main__":
    main()
