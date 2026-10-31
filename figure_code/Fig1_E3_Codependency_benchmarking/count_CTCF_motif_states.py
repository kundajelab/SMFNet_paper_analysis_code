#Run from the repository root: python figure_code/Fig1_E3_Codependency_benchmarking/count_CTCF_motif_states.py
#Step 1 of 3.
from pathlib import Path
import pandas as pd
import numpy as np
import h5py
from tqdm import tqdm

input_dir = Path("input_files/Fig1_E3_Codependency_benchmarking")
output_dir = Path("output/Fig1_E3_Codependency_benchmarking/CTCF")
act_states = input_dir / "mESC_SMF-ONT_testRegions_Experiment_states.h5"
pred_states = input_dir / "mESC_SMF-ONT_testRegions_SMFNet_states.h5"
test_regions = input_dir / "avg_TestRegions_SMFNet.tsv"
motif_list = input_dir / "SMFNet_annotated_motif_hits.tsv"
output_dir.mkdir(parents=True, exist_ok=True)
out_motifs = output_dir / "CTCF_motifs.tsv"
out_motifs_states = output_dir / "CTCF_motif_state_counts.tsv"

#Load test regions
test_df = pd.read_csv(test_regions, sep="\t")

#Keep unique regions
regions_df = test_df[["region"]].drop_duplicates().copy()

#Parse region coordinates
regions_df[["chr", "coords"]] = regions_df["region"].str.split(":", expand=True)
regions_df[["region_start", "region_end"]] = regions_df["coords"].str.split("-", expand=True).astype(int)

#Regions are 2048 bp; central 1024 bp is +512 to +1536
regions_df["central_start"] = regions_df["region_start"] + 512
regions_df["central_end"] = regions_df["region_start"] + 1536

regions_df = regions_df.drop(columns=["coords"])

#Load motif list

motif_df = pd.read_csv(motif_list, sep="\t")

#Keep only motifs in test regions
motif_df = motif_df[motif_df["region"].isin(regions_df["region"])].copy()

#Add region coordinates
motif_df = motif_df.merge(
    regions_df,
    on="region",
    how="inner"
)

#Keep motifs whose center falls inside central 1024 bp

motif_df["motif_center_abs"] = (motif_df["start"] + motif_df["end"]) / 2

motif_df = motif_df[
    (motif_df["motif_center_abs"] >= motif_df["central_start"]) &
    (motif_df["motif_center_abs"] <= motif_df["central_end"])
].copy()

#Relative coordinates to central 1024 bp window
#central_start becomes 0
#central_end becomes 1024

motif_df["rel_start"] = motif_df["start"] - motif_df["central_start"]
motif_df["rel_end"] = motif_df["end"] - motif_df["central_start"]
motif_df["rel_center"] = motif_df["motif_center_abs"] - motif_df["central_start"]


#Remove overlapping motifs per region
#Keep motif with highest hit_coefficient_global
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

#Resolve overlaps across ALL motifs before selecting CTCF, preserving the original cohort.
motif_filtered = pd.concat(
    [remove_overlapping_motifs(group) for _, group in motif_df.groupby("region")],
    ignore_index=True
)
motif_filtered = motif_filtered.loc[motif_filtered["motif"] == "CTCF"].copy()

#Clean output columns

keep_cols = [
    "region",
    "chr",
    "region_start",
    "region_end",
    "central_start",
    "central_end",
    "start",
    "end",
    "rel_start",
    "rel_end",
    "rel_center",
    "motif",
    "motif_name",
    "pattern",
    "strand",
    "hit_coefficient_global",
    "hit_coefficient",
    "hit_similarity",
    "hit_importance"
]

keep_cols = [c for c in keep_cols if c in motif_filtered.columns]

motif_filtered = motif_filtered[keep_cols].copy()

motif_filtered.to_csv(out_motifs, sep="\t", index=False)

print(f"Input test regions: {regions_df['region'].nunique()}")
print(f"Motifs in test regions and central 1024 bp before overlap filtering: {len(motif_df)}")
print(f"CTCF motifs after overlap filtering: {len(motif_filtered)}")
print(f"Wrote: {out_motifs}")

#Add state counts / percentages at closest SMF position

def get_state_vector_at_motif(h5, region, rel_center, max_dist=7.5):
    """
    For one motif, find closest smf_pos in one H5 file and return the state vector
    at that position.

    Assumes:
      /smf_pos/{region} gives relative positions
      /smf_states/{region} gives state matrix
      smf_pos corresponds to one axis of smf_states

    Handles both possible matrix orientations:
      states shape = molecules x positions
      states shape = positions x molecules
    """

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

    #smf_pos indicates columns
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
    """
    Collapse states:
      -1 = missing, removed from denominator
       0/1 = Accessible
       2 = Bound
       3 = Nucleosomal
    """

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

state_rows = []

with h5py.File(act_states, "r") as act_h5, h5py.File(pred_states, "r") as pred_h5:

    for _, row in tqdm(
        motif_filtered.iterrows(),
        total=len(motif_filtered),
        desc="Counting states at motifs"
    ):
        region = row["region"]
        rel_center = row["rel_center"]

        act_hit = get_state_vector_at_motif(
            h5=act_h5,
            region=region,
            rel_center=rel_center,
            max_dist=7.5
        )

        pred_hit = get_state_vector_at_motif(
            h5=pred_h5,
            region=region,
            rel_center=rel_center,
            max_dist=7.5
        )

        #discard motif if no nearby SMF position in either act or pred
        if act_hit is None or pred_hit is None:
            continue

        out = row.to_dict()

        out["act_closest_smf_pos"] = act_hit["closest_smf_pos"]
        out["act_closest_smf_idx"] = act_hit["closest_smf_idx"]
        out["act_closest_smf_dist"] = act_hit["closest_smf_dist"]

        out["pred_closest_smf_pos"] = pred_hit["closest_smf_pos"]
        out["pred_closest_smf_idx"] = pred_hit["closest_smf_idx"]
        out["pred_closest_smf_dist"] = pred_hit["closest_smf_dist"]

        out.update(summarize_states(act_hit["state_vec"], prefix="act"))
        out.update(summarize_states(pred_hit["state_vec"], prefix="pred"))

        state_rows.append(out)

if not state_rows:
    raise ValueError("No CTCF motifs have a measured site within 7.5 bp in both datasets.")
motif_state_df = pd.DataFrame(state_rows)

motif_state_df.to_csv(out_motifs_states, sep="\t", index=False)

print(f"Motifs before SMF-position filtering: {len(motif_filtered)}")
print(f"Motifs after requiring closest SMF position within 7.5 bp in both act and pred: {len(motif_state_df)}")
print(f"Wrote: {out_motifs_states}")
