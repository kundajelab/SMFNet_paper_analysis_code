#Run from the repository root: python figure_code/Fig5_E8_F1_analysis/cluster_variant_windows.py
from pathlib import Path
Path("output/Fig5_E8_F1_analysis").mkdir(parents=True, exist_ok=True)

import os

import h5py
import numpy as np
import pandas as pd

from tqdm import tqdm

from sklearn.cluster import (
    HDBSCAN,
    AgglomerativeClustering
)

from sklearn.metrics import (
    pairwise_distances
)
#Input and output
BASE = "output/Fig5_E8_F1_analysis"

#Already expanded, aligned and filtered matrices.
INPUT_H5 = (
    f"{BASE}/"
    "F1_expanded_400bp_windows.h5"
)

#Used only to attach motif annotations.
WINDOWS = (
    f"{BASE}/"
    "F1_atacPeaks_trainRegions_variant_motif_windows.tsv"
)

OUT_H5 = (
    f"{BASE}/"
    "F1_all4_balanced_clusters.h5"
)

OUT_SUMMARY = (
    f"{BASE}/"
    "F1_all4_balanced_cluster_fractions.tsv"
)

OUT_DIAGNOSTICS = (
    f"{BASE}/"
    "F1_all4_clustering_diagnostics.tsv"
)
#Parameters
MIN_READS_PER_SOURCE = 50

MIN_CLUSTER_SIZE = 5
MIN_SAMPLES = 1

MAX_CLUSTERS = 8

SEED = 42

rng = np.random.default_rng(SEED)

CONDITIONS = (
    "actual_b6",
    "actual_spret",
    "predicted_b6",
    "predicted_spret"
)
#Load annotation table
annotations = pd.read_csv(
    WINDOWS,
    sep="\t"
)

annotations.columns = (
    annotations.columns.str.strip()
)

#Match directly to the identifiers stored
#In the expanded-window HDF5.
required_columns = (
    "region_id",
    "motif_start",
    "motif_end"
)

missing_columns = [
    col
    for col in required_columns
    if col not in annotations.columns
]

if missing_columns:
    raise ValueError(
        "Missing annotation columns: "
        + ", ".join(missing_columns)
    )

annotations["motif_start"] = pd.to_numeric(
    annotations["motif_start"]
)

annotations["motif_end"] = pd.to_numeric(
    annotations["motif_end"]
)

annotations = annotations.dropna(
    subset=[
        "region_id",
        "motif_start",
        "motif_end"
    ]
).copy()

annotations["motif_start"] = (
    annotations["motif_start"].astype(int)
)

annotations["motif_end"] = (
    annotations["motif_end"].astype(int)
)

#Retain all original annotation columns.
#If several rows share the same region and motif
#Coordinates, use the first and record how many
#Annotation rows matched.
annotation_lookup = {}

for key, group in annotations.groupby(
    [
        "region_id",
        "motif_start",
        "motif_end"
    ],
    sort=False
):
    annotation_lookup[key] = (
        group.iloc[0].to_dict(),
        len(group)
    )

print(
    f"Loaded {len(annotations):,} annotation rows"
)
#Helpers
def get_annotation(window_group):
    """
    Match HDF5 window to the annotation table.

    The input HDF5 identifies windows using:
        region_id
        motif_start
        motif_end
    """
    region = str(
        window_group.attrs["region_id"]
    )

    motif_start = int(
        window_group.attrs["motif_start"]
    )

    motif_end = int(
        window_group.attrs["motif_end"]
    )

    key = (
        region,
        motif_start,
        motif_end
    )

    if key not in annotation_lookup:
        return {}, 0

    annotation, n_matches = (
        annotation_lookup[key]
    )

    return annotation, n_matches

def get_medoids(X, labels):
    """
    Calculate the medoid of each non-noise cluster
    using Hamming distance.

    Cluster labels must be consecutive:
        0, 1, ..., k-1
    """
    unique_labels = sorted(
        int(x)
        for x in np.unique(labels)
        if x >= 0
    )

    result = []

    for label in unique_labels:
        members = X[labels == label]

        distances = pairwise_distances(
            members,
            metric="hamming"
        )

        medoid_index = np.argmin(
            distances.sum(axis=1)
        )

        result.append(
            members[medoid_index]
        )

    return np.vstack(result).astype(
        np.int8
    )

def assign_to_medoids(mat, medoids):
    """
    Assign every read to its nearest medoid.

    All reads receive a label, including:
      - reads not included in the balanced sample
      - reads labeled noise by HDBSCAN

    Returns labels 0 to k-1.
    """
    if len(mat) == 0:
        return np.empty(0, dtype=np.int32)

    distances = pairwise_distances(
        mat.astype(bool),
        medoids.astype(bool),
        metric="hamming"
    )

    return np.argmin(
        distances,
        axis=1
    ).astype(np.int32)

def cluster_balanced_reads(X):
    """
    Run HDBSCAN on equally sampled reads
    from all four conditions.

    Allow 1 to MAX_CLUSTERS clusters.

    If HDBSCAN finds more than MAX_CLUSTERS,
    merge its cluster medoids using hierarchical
    clustering.

    If HDBSCAN finds no clusters, use one medoid.
    """
    model = HDBSCAN(
        min_cluster_size=MIN_CLUSTER_SIZE,
        min_samples=MIN_SAMPLES,
        metric="hamming",
        algorithm="brute",
        cluster_selection_method="eom"
    )

    original_labels = model.fit_predict(
        X.astype(bool)
    )

    noise_fraction = float(
        np.mean(original_labels == -1)
    )

    non_noise_labels = sorted(
        int(label)
        for label in np.unique(original_labels)
        if label != -1
    )

    k_initial = len(non_noise_labels)
    #No clusters: use one medoid for every read.
    if k_initial == 0:
        d = pairwise_distances(
            X.astype(bool),
            metric="hamming"
        )

        medoid_idx = np.argmin(
            d.sum(axis=1)
        )

        medoids = X[
            medoid_idx:medoid_idx + 1
        ].astype(np.int8)

        final_labels = np.zeros(
            len(X),
            dtype=np.int32
        )

        return {
            "medoids": medoids,
            "balanced_hdbscan_labels": (
                original_labels.astype(np.int32)
            ),
            "balanced_final_labels": final_labels,
            "k_initial": 0,
            "k_final": 1,
            "noise_fraction": noise_fraction
        }

    #Remap original HDBSCAN labels to 0 ... k-1,
    #Preserving -1 for noise.
    label_map = {
        old: new
        for new, old in enumerate(
            non_noise_labels
        )
    }

    labels = np.array(
        [
            label_map[label]
            if label != -1
            else -1
            for label in original_labels
        ],
        dtype=np.int32
    )
    #Merge clusters if HDBSCAN found more than 8.
    if k_initial > MAX_CLUSTERS:
        initial_medoids = get_medoids(
            X,
            labels
        )

        medoid_distances = pairwise_distances(
            initial_medoids.astype(bool),
            metric="hamming"
        )

        merger = AgglomerativeClustering(
            n_clusters=MAX_CLUSTERS,
            metric="precomputed",
            linkage="average"
        )

        merged = merger.fit_predict(
            medoid_distances
        )

        merged_labels = np.full(
            len(labels),
            -1,
            dtype=np.int32
        )

        for old_label in range(k_initial):
            merged_labels[
                labels == old_label
            ] = merged[old_label]

        labels = merged_labels

    #Renumber final non-noise clusters
    #Consecutively, preserving noise.
    surviving_labels = sorted(
        int(label)
        for label in np.unique(labels)
        if label != -1
    )

    final_map = {
        old: new
        for new, old in enumerate(
            surviving_labels
        )
    }

    labels = np.array(
        [
            final_map[label]
            if label != -1
            else -1
            for label in labels
        ],
        dtype=np.int32
    )

    #Recalculate medoids using the balanced reads
    #After merging clusters.
    medoids = get_medoids(
        X,
        labels
    )

    #Assign all balanced reads to a final medoid,
    #Including HDBSCAN noise reads.
    final_labels = assign_to_medoids(
        X,
        medoids
    )

    return {
        "medoids": medoids,
        "balanced_hdbscan_labels": (
            original_labels.astype(np.int32)
        ),
        "balanced_final_labels": final_labels,
        "k_initial": k_initial,
        "k_final": len(medoids),
        "noise_fraction": noise_fraction
    }

def write_annotation_attrs(group, annotation):
    """
    Store annotation values as HDF5 attributes.

    Missing annotation values are skipped.
    """
    for column, value in annotation.items():
        if column in (
            "region_id",
            "motif_start",
            "motif_end"
        ):
            continue

        if pd.isna(value):
            continue

        name = f"annotation_{column}"

        if isinstance(
            value,
            (np.integer, int)
        ):
            group.attrs[name] = int(value)

        elif isinstance(
            value,
            (np.floating, float)
        ):
            group.attrs[name] = float(value)

        else:
            group.attrs[name] = str(value)
#Process input HDF5 windows
summary_rows = []
diagnostic_rows = []

n_total = 0
n_retained = 0
n_skipped_reads = 0

with (
    h5py.File(INPUT_H5, "r") as input_h5,
    h5py.File(OUT_H5, "w") as output_h5
):
    input_windows = input_h5["windows"]

    output_windows = output_h5.create_group(
        "windows"
    )

    output_h5.attrs["input_h5"] = INPUT_H5
    output_h5.attrs["min_reads_per_source"] = (
        MIN_READS_PER_SOURCE
    )
    output_h5.attrs["max_clusters"] = (
        MAX_CLUSTERS
    )
    output_h5.attrs["random_seed"] = SEED

    output_h5.attrs["clustering"] = (
        "Balanced four-condition HDBSCAN, "
        "maximum 8 medoid clusters; "
        "all reads assigned to nearest medoid"
    )

    #Preserve input preprocessing information
    #Where available.
    for attr in (
        "window_size",
        "max_missing",
        "max_nuc_frac",
        "actual_position_offset",
        "predicted_position_offset",
        "interpolation",
        "filter_positions"
    ):
        if attr in input_h5.attrs:
            output_h5.attrs[
                f"input_{attr}"
            ] = input_h5.attrs[attr]

    for window_id in tqdm(
        input_windows.keys(),
        total=len(input_windows),
        desc="Clustering windows"
    ):
        n_total += 1

        input_group = input_windows[
            window_id
        ]

        region = str(
            input_group.attrs["region_id"]
        )

        motif_start = int(
            input_group.attrs["motif_start"]
        )

        motif_end = int(
            input_group.attrs["motif_end"]
        )

        window_start = int(
            input_group.attrs["window_start"]
        )

        window_end = int(
            input_group.attrs["window_end"]
        )
        #Load all four already-filtered matrices.
        #No interpolation, missingness filtering,
        #Occupancy filtering or position matching.
        matrices = {
            condition: input_group[condition][:]
            for condition in CONDITIONS
        }

        read_counts = {
            condition: len(matrices[condition])
            for condition in CONDITIONS
        }

        #Validate that every input matrix really
        #Represents a 400-bp window.
        for condition, mat in matrices.items():
            if (
                mat.ndim != 2
                or mat.shape[1] != 400
            ):
                raise ValueError(
                    f"{window_id}/{condition}: "
                    f"expected reads x 400, "
                    f"got {mat.shape}"
                )

            if not np.all(
                np.isin(mat, [0, 1])
            ):
                raise ValueError(
                    f"{window_id}/{condition}: "
                    "expected binary expanded states "
                    "containing only 0 and 1"
                )
        #Only filter: >=50 reads per condition.
        n_balance = min(
            read_counts.values()
        )

        annotation, n_annotation_matches = (
            get_annotation(input_group)
        )

        diagnostic = {
            "window_id": window_id,
            "region_id": region,
            "motif_start": motif_start,
            "motif_end": motif_end,
            "window_start": window_start,
            "window_end": window_end,
            "n_actual_b6": (
                read_counts["actual_b6"]
            ),
            "n_actual_spret": (
                read_counts["actual_spret"]
            ),
            "n_predicted_b6": (
                read_counts["predicted_b6"]
            ),
            "n_predicted_spret": (
                read_counts["predicted_spret"]
            ),
            "n_balance": n_balance,
            "n_annotation_matches": (
                n_annotation_matches
            )
        }

        if n_balance < MIN_READS_PER_SOURCE:
            diagnostic["status"] = (
                "skipped_fewer_than_50_reads"
            )

            diagnostic_rows.append(
                diagnostic
            )

            n_skipped_reads += 1

            continue
        #Subsample every condition equally.
        #These indices refer to rows of the
        #Already-filtered input HDF5 matrices.
        balanced_indices = {}

        balanced_matrices = {}

        for condition in CONDITIONS:
            n_reads = read_counts[condition]

            if n_reads == n_balance:
                selected = np.arange(
                    n_reads
                )

            else:
                selected = rng.choice(
                    n_reads,
                    size=n_balance,
                    replace=False
                )

            balanced_indices[
                condition
            ] = selected.astype(
                np.int32
            )

            balanced_matrices[
                condition
            ] = matrices[condition][
                selected
            ]

        #Balanced clustering input:
        #Actual B6
        #Actual SPRET
        #Predicted B6
        #Predicted SPRET
        #Exactly n_balance reads per condition.
        X = np.vstack(
            [
                balanced_matrices[condition]
                for condition in CONDITIONS
            ]
        ).astype(np.int8)
        #Cluster the balanced sample.
        clustering = cluster_balanced_reads(
            X
        )

        medoids = clustering["medoids"]

        k = clustering["k_final"]
        #Reassign all reads from each condition
        #To the nearest final medoid.
        #This includes the sampled reads and
        #The unsampled reads.
        all_labels = {
            condition: assign_to_medoids(
                matrices[condition],
                medoids
            )
            for condition in CONDITIONS
        }

        #Write clustered window to HDF5.
        output_group = (
            output_windows.create_group(
                window_id
            )
        )

        #Copy original window metadata.
        for attr, value in (
            input_group.attrs.items()
        ):
            output_group.attrs[attr] = value

        #Attach new annotation fields.
        write_annotation_attrs(
            output_group,
            annotation
        )

        output_group.attrs[
            "n_annotation_matches"
        ] = n_annotation_matches

        output_group.attrs[
            "n_balance_per_source"
        ] = n_balance

        output_group.attrs[
            "n_clusters"
        ] = k

        output_group.attrs[
            "n_initial_hdbscan_clusters"
        ] = clustering["k_initial"]

        output_group.attrs[
            "hdbscan_noise_fraction"
        ] = clustering[
            "noise_fraction"
        ]

        #Preserve coordinates and the positions
        #Shared by all four input conditions.
        for dataset_name in (
            "relative_positions",
            "genomic_positions",
            "common_positions"
        ):
            if dataset_name in input_group:
                output_group.create_dataset(
                    dataset_name,
                    data=input_group[
                        dataset_name
                    ][:]
                )

        #Final cluster medoids:
        #n_clusters x 400
        output_group.create_dataset(
            "cluster_medoids",
            data=medoids,
            compression="gzip",
            compression_opts=4
        )

        #Original HDBSCAN labels of balanced
        #Samples, including -1 noise labels.
        output_group.create_dataset(
            "balanced_hdbscan_labels",
            data=clustering[
                "balanced_hdbscan_labels"
            ]
        )

        #Final medoid labels for the balanced
        #Clustering sample, with no -1 labels.
        output_group.create_dataset(
            "balanced_final_labels",
            data=clustering[
                "balanced_final_labels"
            ]
        )

        #Store the order in which the four source
        #Matrices were stacked for clustering.
        output_group.attrs[
            "balanced_condition_order"
        ] = ",".join(CONDITIONS)

        #Full matrices, full cluster labels, and
        #Indices of sampled reads for each source.
        for condition in CONDITIONS:
            output_group.create_dataset(
                condition,
                data=matrices[condition],
                compression="gzip",
                compression_opts=4
            )

            output_group.create_dataset(
                f"labels_{condition}",
                data=all_labels[condition]
            )

            output_group.create_dataset(
                f"balanced_indices_{condition}",
                data=balanced_indices[condition]
            )
        #Cluster-fraction summary.
        #Fractions are calculated over all reads
        #From each condition, not just the
        #Balanced clustering sample.
        for label in range(k):
            for condition in CONDITIONS:
                labels = all_labels[
                    condition
                ]

                count = int(
                    np.sum(labels == label)
                )

                total = len(labels)

                summary_row = {
                    "window_id": window_id,
                    "region_id": region,
                    "motif_start": motif_start,
                    "motif_end": motif_end,
                    "window_start": window_start,
                    "window_end": window_end,
                    "k": k,
                    "n_balance_per_source": (
                        n_balance
                    ),
                    "cluster": label,
                    "read_source": condition,
                    "count": count,
                    "fraction": (
                        count / total
                        if total > 0
                        else np.nan
                    ),
                    "total_reads": total
                }

                #Include annotation columns in
                #The summary TSV.
                for col, value in (
                    annotation.items()
                ):
                    if col not in summary_row:
                        summary_row[col] = value

                summary_rows.append(
                    summary_row
                )

        diagnostic.update({
            "status": "retained",
            "k_initial": clustering[
                "k_initial"
            ],
            "k_final": k,
            "hdbscan_noise_fraction": (
                clustering[
                    "noise_fraction"
                ]
            )
        })

        diagnostic_rows.append(
            diagnostic
        )

        n_retained += 1
#Save summaries
pd.DataFrame(
    summary_rows
).to_csv(
    OUT_SUMMARY,
    sep="\t",
    index=False,
    na_rep="NA"
)

pd.DataFrame(
    diagnostic_rows
).to_csv(
    OUT_DIAGNOSTICS,
    sep="\t",
    index=False,
    na_rep="NA"
)
#Final report
print(
    "\n========== CLUSTERING SUMMARY =========="
)

print(
    f"Input windows: {n_total:,}"
)

print(
    f"Retained windows: {n_retained:,}"
)

print(
    "Skipped (<50 reads in at least one condition): "
    f"{n_skipped_reads:,}"
)

print(
    f"\nClustered HDF5: {OUT_H5}"
)

print(
    f"Cluster fractions: {OUT_SUMMARY}"
)

print(
    f"Diagnostics: {OUT_DIAGNOSTICS}"
)
