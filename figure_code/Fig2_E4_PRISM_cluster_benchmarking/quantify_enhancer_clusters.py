#Run from the repository root: python figure_code/Fig2_E4_PRISM_cluster_benchmarking/quantify_enhancer_clusters.py
from pathlib import Path

Path("output/Fig2_E4_PRISM_cluster_benchmarking").mkdir(parents=True, exist_ok=True)

import pandas as pd
import numpy as np
import h5py
from tqdm import tqdm
from sklearn.cluster import HDBSCAN, AgglomerativeClustering
from sklearn.metrics import pairwise_distances, silhouette_score
import matplotlib.pyplot as plt
import matplotlib as mpl
import os

EXTEND = 200
MIN_POS = 15
clust_df = pd.read_csv("output/Fig2_E4_PRISM_cluster_benchmarking/enhancer_clusters_400bp.tsv", sep="\t")
cluster_groups = clust_df.groupby("cluster_region")

states_act = "input_files/Fig2_E4_PRISM_cluster_benchmarking/mESC_SMF-ONT_testRegions_Experiment_states.h5"
states_pred = "input_files/Fig2_E4_PRISM_cluster_benchmarking/mESC_SMF-ONT_testRegions_SMFNet_states.h5"
out_h5 = "output/Fig2_E4_PRISM_cluster_benchmarking/enhancer_clusters_states.h5"
plot_dir = "output/Fig2_E4_PRISM_cluster_benchmarking/cluster_diagnostics/"

Path(plot_dir).mkdir(parents=True, exist_ok=True)

h5_act = h5py.File(states_act, "r")
h5_pred = h5py.File(states_pred, "r")

act_mat_grp = h5_act["smf_states"]
act_pos_grp = h5_act["smf_pos"]

pred_mat_grp = h5_pred["smf_states"]
pred_pos_grp = h5_pred["smf_pos"]

def filter_all_nuc_reads(mat, read_ids, max_frac=0.95):
    frac_1 = mat.mean(axis=1)

    keep = (
        (frac_1 <= max_frac)
    )

    return mat[keep], read_ids[keep]

print("Extracting states matrices from h5...")

cluster_results = {}
#Iterate over clusters in dataframe
for cluster_id, cluster_df in tqdm(cluster_groups, total=clust_df["cluster_region"].nunique()):

    region_id = cluster_df["region"].iloc[0]
    cluster_rel = int(cluster_df["cluster_rel"].iloc[0])

    #Window in 1024 coordinates
    win_start = cluster_rel - EXTEND
    win_end = cluster_rel + EXTEND

    states_mat_act = act_mat_grp[region_id][:]
    states_mat_pred = pred_mat_grp[region_id][:]

    read_ids_act = np.arange(states_mat_act.shape[0])
    read_ids_pred = np.arange(states_mat_pred.shape[0])

    #Find columns in cluster window
    smf_pos_act = act_pos_grp[region_id][:]
    smf_pos_pred = pred_pos_grp[region_id][:]

    #Window masks
    mask_act = (smf_pos_act >= win_start) & (smf_pos_act <= win_end)
    mask_pred = (smf_pos_pred >= win_start) & (smf_pos_pred <= win_end)

    pos_act = smf_pos_act[mask_act]
    pos_pred = smf_pos_pred[mask_pred]

    #Intersection of positions
    common_pos = np.intersect1d(pos_act, pos_pred)

    if len(common_pos) < MIN_POS:
        continue

    #Re-align both to identical positions
    mask_act_common = np.isin(smf_pos_act, common_pos)
    mask_pred_common = np.isin(smf_pos_pred, common_pos)

    sub_pos_act = smf_pos_act[mask_act_common]
    sub_mat_act = states_mat_act[:, mask_act_common]

    sub_pos_pred = smf_pos_pred[mask_pred_common]
    sub_mat_pred = states_mat_pred[:, mask_pred_common]

    cluster_results[cluster_id] = {
        "region": region_id,
        "cluster_rel": cluster_rel,

        "positions_act": sub_pos_act,
        "matrix_act": sub_mat_act,
        "read_ids_act": read_ids_act,

        "positions_pred": sub_pos_pred,
        "matrix_pred": sub_mat_pred,
        "read_ids_pred": read_ids_pred
    }

h5_act.close()
h5_pred.close()

def expand_reads_to_bp(
    mat,
    positions,
    read_ids,
    win_start,
    win_end,
    max_missing=0.15,
):
    """
    Convert sparse SMF states into complete per-base profiles.

    Filtering:
      - discard reads with -1 at the first or last SMF position
      - discard reads with > max_missing fraction of -1

    Expansion:
      - internal -1 positions are ignored
      - each valid state expands toward neighboring valid states
      - transitions occur at the midpoint between valid positions
      - output window is [win_start, win_end)
    """

    mat = mat.copy()
    positions = np.asarray(positions)
    read_ids = np.asarray(read_ids)

    if mat.shape[1] == 0:
        return (
            np.empty((0, win_end - win_start), dtype=np.int8),
            np.array([], dtype=read_ids.dtype),
            np.arange(win_start, win_end),
        )

    #Ensure columns are ordered by genomic position
    order = np.argsort(positions)
    positions = positions[order]
    mat = mat[:, order]

    #Recode states
    mat[np.isin(mat, [0, 1, 2])] = 0
    mat[mat == 3] = 1

    missing = mat == -1
    missing_fraction = missing.mean(axis=1)

    #Require known states at both sparse-window edges
    valid_reads = (
        (mat[:, 0] != -1) &
        (mat[:, -1] != -1) &
        (missing_fraction <= max_missing)
    )

    mat = mat[valid_reads]
    read_ids = read_ids[valid_reads]

    full_positions = np.arange(win_start, win_end)
    expanded = np.empty(
        (mat.shape[0], len(full_positions)),
        dtype=np.int8,
    )

    for i, row in enumerate(mat):
        called = row != -1
        called_positions = positions[called]
        called_states = row[called]

        #Insertion point of every base among the valid SMF positions
        right_idx = np.searchsorted(
            called_positions,
            full_positions,
            side="left",
        )

        left_idx = np.clip(right_idx - 1, 0, len(called_positions) - 1)
        right_idx = np.clip(right_idx, 0, len(called_positions) - 1)

        left_distance = np.abs(
            full_positions - called_positions[left_idx]
        )
        right_distance = np.abs(
            called_positions[right_idx] - full_positions
        )

        #Nearest valid call wins; an exact midpoint uses the left call
        use_right = right_distance < left_distance
        nearest_idx = np.where(use_right, right_idx, left_idx)

        expanded[i] = called_states[nearest_idx]

    return expanded, read_ids, full_positions

print("Expanding state matrices...")

clean_cluster_results = {}

for cluster_id, data in tqdm(
    cluster_results.items(),
    desc="Expanding clusters",
):
    cluster_rel = data["cluster_rel"]
    win_start = cluster_rel - EXTEND
    win_end = cluster_rel + EXTEND

    expanded_act, read_ids_act, bp_positions_act = expand_reads_to_bp(
        mat=data["matrix_act"],
        positions=data["positions_act"],
        read_ids=data["read_ids_act"],
        win_start=win_start,
        win_end=win_end,
        max_missing=0.15,
    )

    expanded_pred, read_ids_pred, bp_positions_pred = expand_reads_to_bp(
        mat=data["matrix_pred"],
        positions=data["positions_pred"],
        read_ids=data["read_ids_pred"],
        win_start=win_start,
        win_end=win_end,
        max_missing=0.15,
    )

    expanded_act, read_ids_act = filter_all_nuc_reads(
            expanded_act,
            read_ids_act,
        )

    expanded_pred, read_ids_pred = filter_all_nuc_reads(
            expanded_pred,
            read_ids_pred,
        )

    print(
        f"{cluster_id}: "
        f"ACT {len(read_ids_act)}/{len(data['read_ids_act'])}, "
        f"PRED {len(read_ids_pred)}/{len(data['read_ids_pred'])}"
    )

    clean_cluster_results[cluster_id] = {
        "region": data["region"],
        "cluster_rel": cluster_rel,

        #Both are the complete 400-bp coordinate axis
        "positions_act": bp_positions_act,
        "positions_pred": bp_positions_pred,

        "clean_matrix_act": expanded_act,
        "clean_read_ids_act": read_ids_act,

        "clean_matrix_pred": expanded_pred,
        "clean_read_ids_pred": read_ids_pred,
    }

MIN_READS = 100
MIN_CLUSTER_SIZE = 5
MAX_CLUSTERS = 10
MIN_CLUSTERS = 4
MIN_SAMPLES = 1

print("Clustering regions...")
cluster_model_results = {}

for cluster_id, data in tqdm(
    clean_cluster_results.items(),
    desc="Clustering regions",
):
    act_mat = data["clean_matrix_act"]
    pred_mat = data["clean_matrix_pred"]

    act_reads = data["clean_read_ids_act"]
    pred_reads = data["clean_read_ids_pred"]

    #Balanced training sample

    n = min(len(act_reads), len(pred_reads))

    if n < MIN_READS:
        print(f"Skipping {cluster_id}, n = {n}")
        continue

    rng = np.random.default_rng(0)

    act_idx = rng.choice(
        len(act_reads),
        size=n,
        replace=False,
    )
    pred_idx = rng.choice(
        len(pred_reads),
        size=n,
        replace=False,
    )

    act_sub = act_mat[act_idx]
    pred_sub = pred_mat[pred_idx]

    X_train = np.vstack([act_sub, pred_sub,]).astype(bool)

    #Fit HDBSCAN to the balanced training data

    clusterer = HDBSCAN(
        min_cluster_size=MIN_CLUSTER_SIZE,
        min_samples=MIN_SAMPLES,
        metric="hamming",
        algorithm="brute",
        cluster_selection_method="eom",
    )

    labels_train = clusterer.fit_predict(X_train)

    #HDBSCAN noise has label -1
    original_labels = sorted(
        label
        for label in np.unique(labels_train)
        if label != -1
    )

    if len(original_labels) == 0:
        print(f"Skipping {cluster_id}: HDBSCAN found no clusters")
        continue

    #Remap cluster labels to consecutive values 0, 1, ...
    label_map = {
        old_label: new_label
        for new_label, old_label in enumerate(original_labels)
    }

    labels_train = np.array([
        label_map[label] if label != -1 else -1
        for label in labels_train
    ])

    k = len(original_labels)

    if k < MIN_CLUSTERS:
        print(f"Skipping region: only {k} HDBSCAN clusters")
        continue

    #Restrict HDBSCAN result to at most MAX_CLUSTERS

    if k > MAX_CLUSTERS:

        #Calculate one preliminary medoid per HDBSCAN cluster
        preliminary_medoids = []

        for label in range(k):
            members = X_train[labels_train == label]

            distances = pairwise_distances(
                members,
                metric="hamming",
            )

            medoid_index = np.argmin(
                distances.sum(axis=1)
            )

            preliminary_medoids.append(
                members[medoid_index]
            )

        preliminary_medoids = np.vstack(
            preliminary_medoids
        )

        #Hamming distances between HDBSCAN cluster medoids
        medoid_distances = pairwise_distances(
            preliminary_medoids,
            metric="hamming",
        )

        #Merge the original HDBSCAN clusters down to MAX_CLUSTERS
        merger = AgglomerativeClustering(
            n_clusters=MAX_CLUSTERS,
            metric="precomputed",
            linkage="average",
        )

        merged_cluster_labels = merger.fit_predict(
            medoid_distances
        )

        #Map each original HDBSCAN label to its merged label
        merged_labels_train = np.full_like(
            labels_train,
            fill_value=-1,
        )

        for original_label in range(k):
            merged_label = merged_cluster_labels[
                original_label
            ]

            merged_labels_train[
                labels_train == original_label
            ] = merged_label

        labels_train = merged_labels_train

        #Ensure merged labels are consecutive 0, 1, ...
        merged_unique = sorted(
            label
            for label in np.unique(labels_train)
            if label != -1
        )

        merged_label_map = {
            old_label: new_label
            for new_label, old_label in enumerate(
                merged_unique
            )
        }

        labels_train = np.array([
            (
                merged_label_map[label]
                if label != -1
                else -1
            )
            for label in labels_train
        ])

        k = len(merged_unique)

        print(
            f"{cluster_id}: merged "
            f"{len(original_labels)} HDBSCAN clusters "
            f"into {k} final clusters"
        )

    #Calculate a Hamming medoid

    medoids = []

    for label in range(k):
        members = X_train[labels_train == label]

        #Pairwise Hamming distances among cluster members
        member_distances = pairwise_distances(
            members,
            metric="hamming",
        )

        #Member with the lowest total distance to other members
        medoid_index = np.argmin(
            member_distances.sum(axis=1)
        )
        medoid = members[medoid_index]

        #Distances of training members to the medoid
        distances_to_medoid = pairwise_distances(
            members,
            medoid.reshape(1, -1),
            metric="hamming",
        ).ravel()

        medoids.append(medoid)

    medoids = np.vstack(medoids)

    #Assign every experiment and SMFNet read to nearest medoid

    act_full = act_mat
    pred_full = pred_mat

    X_full = np.vstack([ act_full, pred_full,]).astype(bool)

    medoids = medoids.astype(bool)

    n_act = len(act_full)
    n_pred = len(pred_full)

    #Hamming distance from every read to every final medoid
    full_distances = pairwise_distances(
        X_full,
        medoids,
        metric="hamming",
    )

    #Every read is assigned to its nearest cluster
    labels_full = np.argmin(
        full_distances,
        axis=1,
    )

    labels_act = labels_full[:n_act]
    labels_pred = labels_full[n_act:]

    #Calculate final silhouette on balanced training reads

    X_train_distances = pairwise_distances(
        X_train.astype(bool),
        medoids,
        metric="hamming",
    )

    #Final labels for the balanced training data, including
    #Reads that HDBSCAN originally called noise
    labels_train_final = np.argmin(
        X_train_distances,
        axis=1,
    )

    represented_labels = np.unique(
        labels_train_final
    )

    if (
        len(represented_labels) >= 2
        and len(labels_train_final) > len(represented_labels)
    ):
        silhouette = silhouette_score(
            X_train.astype(bool),
            labels_train_final,
            metric="hamming",
        )
    else:
        silhouette = np.nan

    #Diagnostic only: fraction HDBSCAN originally called noise
    hdbscan_noise_fraction = np.mean(
        labels_train == -1
    )

    print(
        f"{cluster_id}: "
        f"k={k}, "
        f"silhouette={silhouette:.3f}, "
        f"HDBSCAN discovery noise="
        f"{hdbscan_noise_fraction:.1%}, "
        f"final unassigned=0.0%"
    )

    #Save region results

    cluster_model_results[cluster_id] = {
        "region": data["region"],
        "cluster_rel": data["cluster_rel"],

        "k": k,
        "silhouette": silhouette,

        "act_matrix": act_full,
        "pred_matrix": pred_full,

        "smf_pos_act": data["positions_act"],
        "smf_pos_pred": data["positions_pred"],

        "read_ids_act": act_reads,
        "read_ids_pred": pred_reads,

        #Balanced HDBSCAN training reads
        "act_reads_used": act_reads[act_idx],
        "pred_reads_used": pred_reads[pred_idx],

        #Labels for every read
        "labels_act_all": labels_act,
        "labels_pred_all": labels_pred,

        #Out-of-sample assignment information
        "cluster_medoids": medoids.astype(np.int8),
        "hdbscan_noise_fraction": hdbscan_noise_fraction,
    }

print("Writing cluster data to h5...")

with h5py.File(out_h5, "w") as h5f:
    for cluster_id, data in cluster_model_results.items():

        grp = h5f.create_group(str(cluster_id))

        grp.create_dataset(
            "act_matrix",
            data=data["act_matrix"],
            compression="gzip",
            compression_opts=9,
        )
        grp.create_dataset(
            "pred_matrix",
            data=data["pred_matrix"],
            compression="gzip",
            compression_opts=9,
        )

        grp.create_dataset(
            "smf_pos_act",
            data=data["smf_pos_act"],
        )
        grp.create_dataset(
            "smf_pos_pred",
            data=data["smf_pos_pred"],
        )

        grp.create_dataset(
            "read_ids_act",
            data=data["read_ids_act"],
        )
        grp.create_dataset(
            "read_ids_pred",
            data=data["read_ids_pred"],
        )

        grp.create_dataset(
            "act_reads_used",
            data=data["act_reads_used"],
        )
        grp.create_dataset(
            "pred_reads_used",
            data=data["pred_reads_used"],
        )

        grp.create_dataset(
            "labels_act",
            data=data["labels_act_all"],
        )
        grp.create_dataset(
            "labels_pred",
            data=data["labels_pred_all"],
        )

        grp.create_dataset(
            "cluster_medoids",
            data=data["cluster_medoids"],
            compression="gzip",
            compression_opts=9,
        )

        grp.create_dataset(
            "cluster_k",
            data=data["k"],
        )
        grp.create_dataset(
            "silhouette",
            data=data["silhouette"],
        )

        grp.create_dataset(
            "hdbscan_noise_fraction",
            data=data["hdbscan_noise_fraction"],
        )

        grp.create_dataset(
            "region_id",
            data=data["region"],
        )
        grp.create_dataset(
            "cluster_center_rel",
            data=data["cluster_rel"],
        )

#Calculate label proportions in long format

print("Calculating label proportions...")

summary_rows = []

for cluster_id, data in tqdm(
    cluster_model_results.items(),
    desc="Calculate fractions",
):
    labels_act = np.asarray(data["labels_act_all"])
    labels_pred = np.asarray(data["labels_pred_all"])

    k = data["k"]

    act_total = len(labels_act)
    pred_total = len(labels_pred)

    labels_to_report = list(range(k))

    for label in labels_to_report:
        act_count = np.sum(labels_act == label)
        pred_count = np.sum(labels_pred == label)

        #Fraction among all reads
        act_frac = (
            act_count / act_total
            if act_total > 0 else np.nan
        )
        pred_frac = (
            pred_count / pred_total
            if pred_total > 0 else np.nan
        )

        summary_rows.append({
            "cluster_id": cluster_id,
            "region": data["region"],
            "cluster_center_rel": data["cluster_rel"],
            "k": k,
            "silhouette": data["silhouette"],

            "label": label,

            "act_count": act_count,
            "pred_count": pred_count,

            "act_frac": act_frac,
            "pred_frac": pred_frac,

            "act_total_reads": act_total,
            "pred_total_reads": pred_total,
        })

summary_df = pd.DataFrame(summary_rows)

summary_df.to_csv(
    "output/Fig2_E4_PRISM_cluster_benchmarking/enhancer_cluster_fractions.tsv",
    sep="\t",
    index=False,
)

#Plot clusters

def plot_cluster(
    cluster_id,
    data,
    outdir,
    tfbs_positions,
):
    TFBS_HALF_WIDTH = 5
    BLOCK_HEIGHT = 20

    act_mat = data["act_matrix"].astype(float)
    pred_mat = data["pred_matrix"].astype(float)

    act_labels = np.asarray(data["labels_act_all"])
    pred_labels = np.asarray(data["labels_pred_all"])

    k = data["k"]

    #Regular clusters
    displayed_labels = list(range(k))
    n_panels = len(displayed_labels)
    act_pos = np.asarray(data["smf_pos_act"])
    pred_pos = np.asarray(data["smf_pos_pred"])

    if np.array_equal(act_pos, pred_pos):
        shared_pos = act_pos
        act_shared = act_mat
        pred_shared = pred_mat
    else:
        shared_pos = np.union1d(act_pos, pred_pos)

        act_shared = np.full(
            (len(act_mat), len(shared_pos)),
            np.nan,
        )
        pred_shared = np.full(
            (len(pred_mat), len(shared_pos)),
            np.nan,
        )

        act_columns = np.searchsorted(shared_pos, act_pos)
        pred_columns = np.searchsorted(shared_pos, pred_pos)

        act_shared[:, act_columns] = act_mat
        pred_shared[:, pred_columns] = pred_mat

    #Helpers

    def mean_for_label(mat, labels, label):
        selected = mat[labels == label]

        if len(selected) == 0:
            return np.full(mat.shape[1], np.nan)

        return np.nanmean(selected, axis=0)

    def heatmap_block(mat, labels, label, height, rng):
        selected = mat[labels == label]

        if len(selected) == 0:
            return np.full(
                (height, mat.shape[1]),
                np.nan,
            )

        #Randomly select reads when there are more than 20
        if len(selected) > height:
            selected_indices = rng.choice(
                len(selected),
                size=height,
                replace=False,
            )
            selected = selected[selected_indices]

        #Pad groups containing fewer than 20 reads
        if len(selected) < height:
            padding = np.full(
                (height - len(selected), mat.shape[1]),
                np.nan,
            )
            selected = np.vstack([selected, padding])

        return selected

    #Figure

    fig = plt.figure(
        figsize=(4 * n_panels, 10)
    )

    grid = fig.add_gridspec(
        3,
        n_panels,
        height_ratios=[1, 1, 3],
    )

    cmap = mpl.colors.ListedColormap([
        "lightgrey",
        "black",
    ])
    cmap.set_bad("white")

    x_start = shared_pos.min()
    x_end = shared_pos.max() + 1

    rng = np.random.default_rng(0)

    for panel, label in enumerate(displayed_labels):

        title = f"Label {label}"

        #Mean profiles

        ax_mean = fig.add_subplot(grid[0, panel])

        act_mean = mean_for_label(
            act_shared,
            act_labels,
            label,
        )
        pred_mean = mean_for_label(
            pred_shared,
            pred_labels,
            label,
        )

        ax_mean.plot(
            shared_pos,
            act_mean,
            color="blue",
            label="ACT",
        )
        ax_mean.plot(
            shared_pos,
            pred_mean,
            color="orange",
            label="PRED",
        )

        for tfbs_center in tfbs_positions:
            if x_start <= tfbs_center < x_end:
                ax_mean.axvspan(
                    tfbs_center - TFBS_HALF_WIDTH,
                    tfbs_center + TFBS_HALF_WIDTH,
                    color="grey",
                    alpha=0.25,
                    linewidth=0,
                )

        ax_mean.set_xlim(x_start, x_end)
        ax_mean.set_ylim(-0.05, 1.05)
        ax_mean.set_title(title)

        if panel == 0:
            ax_mean.set_ylabel(
                "Nucleosomal fraction (%)"
            )
            ax_mean.legend()

        #Label proportions

        ax_bar = fig.add_subplot(grid[1, panel])

        act_fraction = np.mean(act_labels == label)
        pred_fraction = np.mean(pred_labels == label)

        ax_bar.bar(
            [0, 1],
            [act_fraction, pred_fraction],
            color=["blue", "orange"],
        )

        ax_bar.set_xticks([0, 1])
        ax_bar.set_xticklabels(["ACT", "PRED"])
        ax_bar.set_ylim(-0.05, 1.05)
        ax_bar.set_title(
            f"{title} fraction"
        )

        if panel == 0:
            ax_bar.set_ylabel("Fraction of all reads")

        #Experiment and SMFNet heatmaps

        ax_heatmap = fig.add_subplot(grid[2, panel])

        act_block = heatmap_block(
            act_shared,
            act_labels,
            label,
            BLOCK_HEIGHT,
            rng,
        )
        pred_block = heatmap_block(
            pred_shared,
            pred_labels,
            label,
            BLOCK_HEIGHT,
            rng,
        )

        #Two blank rows separate experiment and SMFNet
        separator = np.full(
            (2, len(shared_pos)),
            np.nan,
        )

        combined = np.vstack([
            act_block,
            separator,
            pred_block,
        ])

        ax_heatmap.imshow(
            combined,
            aspect="auto",
            cmap=cmap,
            interpolation="nearest",
            vmin=0,
            vmax=1,
            extent=[
                x_start,
                x_end,
                combined.shape[0],
                0,
            ],
        )

        #TFBS locations on heatmaps
        for tfbs_center in tfbs_positions:
            if x_start <= tfbs_center < x_end:
                ax_heatmap.axvline(
                    tfbs_center,
                    color="grey",
                    alpha=0.5,
                    linewidth=0.7,
                )

        ax_heatmap.axhline(
            BLOCK_HEIGHT + 1,
            color="white",
            linewidth=2,
        )

        ax_heatmap.set_xlim(x_start, x_end)
        ax_heatmap.set_yticks([
            BLOCK_HEIGHT / 2,
            BLOCK_HEIGHT + 2 + BLOCK_HEIGHT / 2,
        ])
        ax_heatmap.set_yticklabels([
            "ACT",
            "PRED",
        ])
        ax_heatmap.set_xlabel(
            "Position in 1024-bp region"
        )
        ax_heatmap.set_title(title)

    fig.suptitle(
            (
                f"{cluster_id}\n"
                f"k={k}, "
                f"silhouette={data['silhouette']:.3f}"
            ),
            y=1.02,
        )
    os.makedirs(outdir, exist_ok=True)

    safe_cluster_id = str(cluster_id).replace("/", "_")
    outpath = os.path.join(
        outdir,
        f"{safe_cluster_id}.pdf",
    )

    fig.subplots_adjust(
        left=0.05,
        right=0.98,
        bottom=0.07,
        top=0.92,
        wspace=0.35,
        hspace=0.45,
    )

    fig.savefig(outpath)
    plt.close()

print("Plotting TFBS clusters...")

for cluster_id, data in tqdm(
    cluster_model_results.items(),
    desc="Plotting clusters",
):
    tfbs_df = clust_df[
        clust_df["cluster_region"] == cluster_id
    ]

    tfbs_positions = tfbs_df["rel_1024"].values

    plot_cluster(
        cluster_id,
        data,
        plot_dir,
        tfbs_positions,
    )
