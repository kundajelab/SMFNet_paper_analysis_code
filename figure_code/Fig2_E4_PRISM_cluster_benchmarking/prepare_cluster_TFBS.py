#Run from the repository root: python figure_code/Fig2_E4_PRISM_cluster_benchmarking/prepare_cluster_TFBS.py
from pathlib import Path

Path("output/Fig2_E4_PRISM_cluster_benchmarking").mkdir(parents=True, exist_ok=True)

import pandas as pd
import numpy as np
import h5py
from intervaltree import IntervalTree

print("building interval trees...")

trees = {"chr1": IntervalTree(),
         "chr3": IntervalTree(),
         "chr6": IntervalTree()}

#Load file
df = pd.read_csv("input_files/Fig2_E4_PRISM_cluster_benchmarking/SMFNet_annotated_motif_hits_for_codependency_analysis.tsv", sep="\t")
H5_ACT = "input_files/Fig2_E4_PRISM_cluster_benchmarking/mESC_SMF-ONT_testRegions_Experiment_reads.h5"

#Filter chromosomes
print("loading motif file..")
df = df[df["chr"].isin(["chr1", "chr3", "chr6"])].copy()

#Compute motif centers and remove overlapping motifs
print("removing overlaps..")

df["center"] = ((df["start"] + df["end"]) / 2).astype(int)

def deduplicate(group):
    #Sort by score (descending)
    group = group.sort_values("hit_coefficient_global", ascending=False)

    tree = IntervalTree()
    kept_rows = []

    for row in group.itertuples(index=False):
        s1 = row.start
        e1 = row.end

        #Check overlap in O(log n)
        if tree.overlaps(s1, e1):
            continue

        #Accept motif
        kept_rows.append(row)
        tree[s1:e1] = True  #mark occupied region

    return pd.DataFrame(kept_rows)

df = pd.concat([deduplicate(group) for _, group in df.groupby("chr", sort=True)], ignore_index=True)

#Load H5 regions
print("loading h5 data..")
region_dict = {}

with h5py.File(H5_ACT, "r") as h5f:
    #The supplied sample H5 preserves the original 2048 bp measured positions.
    if h5f.attrs.get("schema") != "smfnet.samples":
        raise ValueError("Expected the published smfnet.samples H5 schema.")
    ids = h5f["region_ids"].asstr()[()]
    index = {region: f"regions/{i:08d}" for i, region in enumerate(ids)}
    for region_name in sorted(index):
        chrom, coords = region_name.split(":")
        r_start, r_end = map(int, coords.split("-"))
        smf_positions = h5f[index[region_name] + "/observed/positions"][()]

        region_dict[region_name] = {
            "chr": chrom,
            "start": r_start,
            "end": r_end,
            "smf_pos": smf_positions
        }

        if chrom in trees:
            trees[chrom][r_start:r_end] = region_name

#Map motifs to regions
def assign_region(df):
    regions = []

    for chrom, center in zip(df["chr"], df["center"]):
        hits = trees[chrom][center]

        if len(hits) == 0:
            regions.append(None)
        else:
            regions.append(list(hits)[0].data)

    return regions

print("Matching motifs to H5 regions")

df["region"] = assign_region(df)
df = df.dropna(subset=["region"]).copy()

#Compute relative positions
rel_2048 = []
rel_1024 = []
closest_smf_list = []
closest_smf_1024_list = []
smf_pos_dists = []

keep_mask = []

for idx, row in df.iterrows():
    region = region_dict[row["region"]]

    r_start = region["start"]
    r_end = region["end"]

    #Relative to 2048
    rel_pos = row["center"] - r_start

    #Centered 1024 window
    mid = (r_start + r_end) // 2
    win_start = mid - 512
    win_end = mid + 512

    rel_pos_1024 = row["center"] - win_start

    #Skip if outside 1024 window
    if rel_pos_1024 < 0 or rel_pos_1024 >= 1024:
        keep_mask.append(False)
        rel_2048.append(np.nan)
        rel_1024.append(np.nan)
        closest_smf_list.append(np.nan)
        closest_smf_1024_list.append(np.nan)
        smf_pos_dists.append(np.nan)
        continue

    #SMF matching
    smf_positions = region["smf_pos"]

    if len(smf_positions) == 0:
        keep_mask.append(False)
        rel_2048.append(np.nan)
        rel_1024.append(np.nan)
        closest_smf_list.append(np.nan)
        closest_smf_1024_list.append(np.nan)
        smf_pos_dists.append(np.nan)
        continue

    #Distance to nearest SMF site
    dists = np.abs(smf_positions - rel_pos)
    closest_smf = int(smf_positions[np.argmin(dists)])
    closest_smf_1024 = closest_smf - 512
    min_dist = dists.min()

    rel_2048.append(rel_pos)
    rel_1024.append(rel_pos_1024)
    closest_smf_list.append(closest_smf)
    closest_smf_1024_list.append(closest_smf_1024)
    smf_pos_dists.append(min_dist)

    keep_mask.append(True)

#Assign
df["rel_2048"] = rel_2048
df["rel_1024"] = rel_1024
df["smf_pos_2048"] = closest_smf_list
df["smf_pos_1024"] = closest_smf_1024_list
df["dist_to_smf_pos"] = smf_pos_dists

#Filter
df = df[keep_mask].copy()

print(f"Remaining motifs after SMF filter: {len(df)}")

df.to_csv("output/Fig2_E4_PRISM_cluster_benchmarking/cluster_TFBS.tsv", sep="\t", index=False)
