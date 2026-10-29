#Run from SMFNet_paper_analysis_code github repo: python ./figure_code/Fig1_E1_ONT_benchmarking/count_SMFNet_motif_states.py
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from tqdm import tqdm

INPUT_DIR = Path("input_files/Fig1_E1_ONT_benchmarking")
OUTPUT_DIR = Path("output/Fig1_E1_ONT_benchmarking")


def remove_overlapping_motifs(df):
    df = df.sort_values("hit_coefficient_global", ascending=False).copy()

    kept_rows = []
    kept_intervals = []

    for _, row in df.iterrows():
        start = row["start"]
        end = row["end"]

        overlaps = False
        for kept_start, kept_end in kept_intervals:
            if start < kept_end and end > kept_start:
                overlaps = True
                break

        if not overlaps:
            kept_rows.append(row)
            kept_intervals.append((start, end))

    if len(kept_rows) == 0:
        return pd.DataFrame(columns=df.columns)

    return pd.DataFrame(kept_rows)


def get_state_vector_at_motif(h5, region, rel_center, max_dist=7.5):
    if region not in h5["smf_pos"]:
        return None

    if region not in h5["smf_states"]:
        return None

    smf_pos = np.asarray(h5["smf_pos"][region][:], dtype=float)
    states = np.asarray(h5["smf_states"][region][:])

    if len(smf_pos) == 0:
        return None

    closest_idx = int(np.argmin(np.abs(smf_pos - rel_center)))
    closest_pos = float(smf_pos[closest_idx])
    closest_dist = abs(closest_pos - rel_center)

    if closest_dist > max_dist:
        return None

    if states.ndim != 2:
        raise ValueError(f"Expected 2D smf_states for {region}, got shape {states.shape}")

    if states.shape[1] == len(smf_pos):
        state_vec = states[:, closest_idx]
    elif states.shape[0] == len(smf_pos):
        state_vec = states[closest_idx, :]
    else:
        raise ValueError(
            f"Cannot match smf_pos to smf_states for {region}: "
            f"len(smf_pos)={len(smf_pos)}, states.shape={states.shape}"
        )

    return {
        "closest_smf_pos": closest_pos,
        "closest_smf_idx": closest_idx,
        "closest_smf_dist": closest_dist,
        "state_vec": state_vec
    }


def summarize_states(state_vec, prefix):
    state_vec = np.asarray(state_vec)

    valid = state_vec[state_vec != -1]
    n_valid = len(valid)

    accessible_n = int(np.sum((valid == 0) | (valid == 1)))
    bound_n = int(np.sum(valid == 2))
    nucleosomal_n = int(np.sum(valid == 3))
    missing_n = int(np.sum(state_vec == -1))
    total_n = int(len(state_vec))

    if n_valid > 0:
        accessible_perc = accessible_n / n_valid
        bound_perc = bound_n / n_valid
        nucleosomal_perc = nucleosomal_n / n_valid
    else:
        accessible_perc = np.nan
        bound_perc = np.nan
        nucleosomal_perc = np.nan

    return {
        f"{prefix}_n_total": total_n,
        f"{prefix}_n_valid": n_valid,
        f"{prefix}_n_missing": missing_n,

        f"{prefix}_accessible_n": accessible_n,
        f"{prefix}_bound_n": bound_n,
        f"{prefix}_nucleosomal_n": nucleosomal_n,

        f"{prefix}_accessible_perc": accessible_perc,
        f"{prefix}_bound_perc": bound_perc,
        f"{prefix}_nucleosomal_perc": nucleosomal_perc
    }


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    regions = pd.read_csv(INPUT_DIR / "avg_TestRegions_SMFNet.tsv", sep="\t", usecols=["region"])
    regions = regions.drop_duplicates().copy()
    regions[["chr", "coords"]] = regions["region"].str.split(":", expand=True)
    regions[["region_start", "region_end"]] = regions["coords"].str.split("-", expand=True).astype(int)
    regions["central_start"] = regions["region_start"] + 512
    regions["central_end"] = regions["region_start"] + 1536
    regions = regions.drop(columns="coords")

    motifs = pd.read_csv(INPUT_DIR / "SMFNet_annotated_motif_hits.tsv", sep="\t")
    motifs = motifs[motifs["region"].isin(regions["region"])].merge(regions, on="region", how="inner")
    motifs["motif_center_abs"] = (motifs["start"] + motifs["end"]) / 2
    motifs = motifs[(motifs["motif_center_abs"] >= motifs["central_start"]) &
                    (motifs["motif_center_abs"] <= motifs["central_end"])].copy()
    motifs["rel_start"] = motifs["start"] - motifs["central_start"]
    motifs["rel_end"] = motifs["end"] - motifs["central_start"]
    motifs["rel_center"] = motifs["motif_center_abs"] - motifs["central_start"]

    #Keep the strongest non-overlapping hits within each central 1024-bp window.
    groups = [remove_overlapping_motifs(group) for _, group in motifs.groupby("region")]
    filtered = pd.concat(groups, ignore_index=True) if groups else motifs.iloc[:0].copy()
    columns = ["region", "chr", "region_start", "region_end", "central_start", "central_end",
               "start", "end", "rel_start", "rel_end", "rel_center", "motif", "motif_name",
               "pattern", "strand", "hit_coefficient_global", "hit_coefficient",
               "hit_similarity", "hit_importance"]
    filtered = filtered[[column for column in columns if column in filtered]].copy()
    filtered.to_csv(OUTPUT_DIR / "SMFNet_test_regions_motifs.tsv", sep="\t", index=False)

    #Use the closest assayed site within 7.5 bp in both datasets; exclude -1 from denominators.
    rows = []
    with h5py.File(INPUT_DIR / "mESC_SMF-ONT_testRegions_Experiment_states.h5", "r") as observed, \
         h5py.File(INPUT_DIR / "mESC_SMF-ONT_testRegions_SMFNet_states.h5", "r") as predicted:
        for _, row in tqdm(filtered.iterrows(), total=len(filtered), desc="Counting states at motifs"):
            act = get_state_vector_at_motif(observed, row["region"], row["rel_center"])
            pred = get_state_vector_at_motif(predicted, row["region"], row["rel_center"])
            if act is None or pred is None:
                continue
            result = row.to_dict()
            for prefix, hit in [("act", act), ("pred", pred)]:
                for key in ("closest_smf_pos", "closest_smf_idx", "closest_smf_dist"):
                    result[f"{prefix}_{key}"] = hit[key]
            #States 0 and 1 both count as accessible; 2 is bound and 3 nucleosomal.
            result.update(summarize_states(act["state_vec"], "act"))
            result.update(summarize_states(pred["state_vec"], "pred"))
            rows.append(result)
    pd.DataFrame(rows).to_csv(OUTPUT_DIR / "SMFNet_motif_state_counts.tsv", sep="\t", index=False)
    print(f"Test regions: {len(regions)}; filtered motifs: {len(filtered)}; matched motifs: {len(rows)}")


if __name__ == "__main__":
    main()
