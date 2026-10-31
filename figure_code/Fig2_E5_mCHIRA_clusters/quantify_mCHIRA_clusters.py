#Run from SMFNet_paper_analysis_code github repo: python ./figure_code/Fig2_E5_mCHIRA_clusters/quantify_mCHIRA_clusters.py

from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import sklearn
from sklearn.cluster import HDBSCAN, AgglomerativeClustering
from sklearn.metrics import pairwise_distances

ACT_H5 = "input_files/Fig2_E5_mCHIRA_clusters/mCHIRA_CRE6005_Experiment_states.h5"
PRED_H5 = "input_files/Fig2_E5_mCHIRA_clusters/mCHIRA_CRE6005_SMFNet_states.h5"

OUT_H5 = "output/Fig2_E5_mCHIRA_clusters/mCHIRA_CRE6005_clusters.h5"
SUMMARY_TSV = "output/Fig2_E5_mCHIRA_clusters/mCHIRA_CRE6005_cluster_fractions.tsv"

Path(OUT_H5).parent.mkdir(parents=True, exist_ok=True)

#Fit on balanced wildtype reads and assign mutants to the same medoids.
TRAIN_GROUP = "CRE_6005"

READ_LENGTH = 252
MAX_MISSING = 0.15
MAX_NUC_FRAC = 0.95

MIN_CLUSTER_SIZE = 5
MIN_SAMPLES = 1
MIN_CLUSTERS = 4
MAX_CLUSTERS = 8

RNG = np.random.default_rng(10)

def expand_reads(mat, pos):
    """
    Convert sparse SMF states to complete 252-bp binary profiles.
    0/1/2 -> 0
    3     -> 1
    """
    mat = mat.copy()
    pos = np.asarray(pos)

    keep_pos = (pos >= 0) & (pos < READ_LENGTH)
    pos = pos[keep_pos]
    mat = mat[:, keep_pos]

    order = np.argsort(pos)
    pos = pos[order]
    mat = mat[:, order]

    mat[np.isin(mat, [0, 1, 2])] = 0
    mat[mat == 3] = 1

    missing = mat == -1

    keep_reads = ((missing.mean(axis=1) <= MAX_MISSING))

    read_ids = np.arange(len(mat))[keep_reads]
    mat = mat[keep_reads]

    full_pos = np.arange(READ_LENGTH)
    expanded = np.empty((len(mat), READ_LENGTH), dtype=np.int8)

    for i, row in enumerate(mat):
        called = row != -1
        called_pos = pos[called]
        called_states = row[called]

        right = np.searchsorted(called_pos, full_pos, side="left")
        left = np.clip(right - 1, 0, len(called_pos) - 1)
        right = np.clip(right, 0, len(called_pos) - 1)

        left_dist = np.abs(full_pos - called_pos[left])
        right_dist = np.abs(called_pos[right] - full_pos)

        nearest = np.where(right_dist < left_dist, right, left)
        expanded[i] = called_states[nearest]

    keep = expanded.mean(axis=1) <= MAX_NUC_FRAC

    print("expanded last:", dict(zip(*np.unique(expanded[:, -1], return_counts=True))))

    return expanded[keep], read_ids[keep]

def preprocess_group(h5, name):
    return expand_reads(
        h5["smf_states"][name][:],
        h5["smf_pos"][name][:],
    )

def preprocess_training_group(h5_act, h5_pred):
    """
    Expand ACT and PRED wildtype reads at their respective assayed positions.
    """
    act_mat = h5_act["smf_states"][TRAIN_GROUP][:]
    pred_mat = h5_pred["smf_states"][TRAIN_GROUP][:]

    act_pos = h5_act["smf_pos"][TRAIN_GROUP][:]
    pred_pos = h5_pred["smf_pos"][TRAIN_GROUP][:]

    act_exp, act_ids = expand_reads(
        act_mat,
        act_pos,
    )
    pred_exp, pred_ids = expand_reads(
        pred_mat,
        pred_pos,
    )

    return act_exp, act_ids, pred_exp, pred_ids

def fit_clusters(act_mat, pred_mat):
    """
    Fit clusters using equal numbers of ACT and PRED CRE_6005 reads.
    Return final medoids only.
    """
    n = min(len(act_mat), len(pred_mat))

    act_sub = act_mat[RNG.choice(len(act_mat), n, replace=False)]
    pred_sub = pred_mat[RNG.choice(len(pred_mat), n, replace=False)]

    X = np.vstack([act_sub, pred_sub]).astype(bool)

    clusterer = HDBSCAN(
        min_cluster_size=MIN_CLUSTER_SIZE,
        min_samples=MIN_SAMPLES,
        metric="hamming",
        algorithm="brute",
        cluster_selection_method="eom",
    )

    labels = clusterer.fit_predict(X)

    valid_labels = sorted(x for x in np.unique(labels) if x != -1)
    label_map = {old: new for new, old in enumerate(valid_labels)}

    labels = np.array([
        label_map[x] if x != -1 else -1
        for x in labels
    ])

    k = len(valid_labels)

    if k > MAX_CLUSTERS:
        medoids = []

        for label in range(k):
            members = X[labels == label]
            d = pairwise_distances(members, metric="hamming")
            medoids.append(members[np.argmin(d.sum(axis=1))])

        medoids = np.vstack(medoids)

        d = pairwise_distances(medoids, metric="hamming")

        merger = AgglomerativeClustering(
            n_clusters=MAX_CLUSTERS,
            metric="precomputed",
            linkage="average",
        )

        merged = merger.fit_predict(d)

        new_labels = np.full_like(labels, -1)

        for old_label in range(k):
            new_labels[labels == old_label] = merged[old_label]

        labels = new_labels
        k = MAX_CLUSTERS

    if k < MIN_CLUSTERS:
        print(f"Warning: only {k} clusters found")

    medoids = []

    for label in sorted(x for x in np.unique(labels) if x != -1):
        members = X[labels == label]
        d = pairwise_distances(members, metric="jaccard")
        medoids.append(members[np.argmin(d.sum(axis=1))])

    return np.vstack(medoids).astype(bool)

def assign_clusters(mat, medoids):
    X = mat.astype(bool)

    d = pairwise_distances(
        X,
        medoids,
        metric="hamming",
    )

    return np.argmin(d, axis=1)

def write_group(out_parent, name, matrix, read_ids, labels):
    grp = out_parent.create_group(name)

    grp.create_dataset(
        "expanded_reads",
        data=matrix,
        compression="gzip",
    )

    grp.create_dataset(
        "read_ids",
        data=read_ids,
    )

    grp.create_dataset(
        "cluster_labels",
        data=labels,
    )

def add_summary(rows, sample, group, labels):
    n = len(labels)

    for cluster in range(labels.max() + 1):
        count = np.sum(labels == cluster)

        rows.append({
            "sample": sample,
            "group": group,
            "cluster": cluster,
            "count": count,
            "fraction": count / n,
        })

with h5py.File(ACT_H5, "r") as act, h5py.File(PRED_H5, "r") as pred:

    print(f"Training clusters on {TRAIN_GROUP}...")

    act_train, act_ids, pred_train, pred_ids = preprocess_training_group(
        act,
        pred,
    )

    medoids = fit_clusters(
        act_train,
        pred_train,
    )

    print(f"Found {len(medoids)} final clusters")

    act_groups = sorted(act["smf_states"].keys())
    pred_groups = sorted(pred["smf_states"].keys())

    summary = []

    with h5py.File(OUT_H5, "w") as out:

        out_act = out.create_group("act")
        out_pred = out.create_group("pred")

        for group in act_groups:

            if group == TRAIN_GROUP:
                matrix = act_train
                read_ids = act_ids
            else:
                matrix, read_ids = preprocess_group(act, group)

            labels = assign_clusters(matrix, medoids)

            write_group(
                out_act,
                group,
                matrix,
                read_ids,
                labels,
            )

            add_summary(
                summary,
                "act",
                group,
                labels,
            )

            print(f"ACT  {group}: {len(labels)} reads")

        for group in pred_groups:

            if group == TRAIN_GROUP:
                matrix = pred_train
                read_ids = pred_ids
            else:
                matrix, read_ids = preprocess_group(pred, group)

            labels = assign_clusters(matrix, medoids)

            write_group(
                out_pred,
                group,
                matrix,
                read_ids,
                labels,
            )

            add_summary(
                summary,
                "pred",
                group,
                labels,
            )

            print(f"PRED {group}: {len(labels)} reads")

summary_df = pd.DataFrame(summary)

summary_df.to_csv(
    SUMMARY_TSV,
    sep="\t",
    index=False,
)

print("\nCluster fractions:")
print(summary_df.to_string(index=False))

print(f"\nWrote: {OUT_H5}")
print(f"Wrote: {SUMMARY_TSV}")
