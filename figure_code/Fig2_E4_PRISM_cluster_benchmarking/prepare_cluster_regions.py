#Run from the repository root: python figure_code/Fig2_E4_PRISM_cluster_benchmarking/prepare_cluster_regions.py
from pathlib import Path

Path("output/Fig2_E4_PRISM_cluster_benchmarking").mkdir(parents=True, exist_ok=True)

import pandas as pd
import numpy as np

df = pd.read_csv("output/Fig2_E4_PRISM_cluster_benchmarking/cluster_TFBS.tsv", sep="\t")

WINDOW = 200
MIN_TFBS = 3
EXTEND = 200  #±200 → 400 bp window

cluster_rows = []

for region_id, group in df.groupby("region"):

    group = group.sort_values("center").reset_index(drop=True)
    centers = group["center"].values

    clusters = []
    current_cluster = [0]

    #Build clusters

    for i in range(1, len(centers)):
        if centers[i] - centers[i-1] <= WINDOW:
            current_cluster.append(i)
        else:
            clusters.append(current_cluster)
            current_cluster = [i]

    clusters.append(current_cluster)  #last cluster

    #Process clusters

    for clust in clusters:

        if len(clust) < MIN_TFBS:
            continue

        cluster_df = group.iloc[clust].copy()

        #Cluster center (median is more robust than mean)
        cluster_center = int(np.median(cluster_df["center"]))

        cluster_start = cluster_center - EXTEND
        cluster_end = cluster_center + EXTEND

        #Parse region start
        chrom, coords = region_id.split(":")
        region_start = int(coords.split("-")[0])

        #Relative position of cluster center within region
        cluster_rel = cluster_center - region_start - 512

        if cluster_rel - EXTEND < 0 or cluster_rel + EXTEND > 1024:
            continue

        #Keep motifs within 400 bp window
        mask = (
            (cluster_df["center"] >= cluster_start) &
            (cluster_df["center"] <= cluster_end)
        )

        cluster_final = cluster_df[mask].copy()

        #Optional: enforce again ≥3 motifs after trimming
        if len(cluster_final) < MIN_TFBS:
            continue

        #Annotate
        cluster_region = f"{chrom}:{cluster_start}-{cluster_end}"

        cluster_final["cluster_region"] = cluster_region
        cluster_final["cluster_center"] = cluster_center
        cluster_final["cluster_size"] = len(cluster_final)
        cluster_final["cluster_rel"] = cluster_rel

        cluster_rows.append(cluster_final)

#Combine + clean

out_df = pd.concat(cluster_rows, ignore_index=True)

#Remove duplicates (should be minimal now, but safe)
out_df = out_df.drop_duplicates(subset=["chr", "center", "cluster_region"])

#Save

out_df.to_csv("output/Fig2_E4_PRISM_cluster_benchmarking/enhancer_clusters_400bp.tsv", sep="\t", index=False)
