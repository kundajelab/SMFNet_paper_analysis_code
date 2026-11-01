#Run from the repository root: python figure_code/Fig5_E8_F1_analysis/compute_variant_metrics.py
from pathlib import Path
Path("output/Fig5_E8_F1_analysis").mkdir(parents=True, exist_ok=True)

import numpy as np
import pandas as pd
import h5py
from tqdm import tqdm
from numba import njit
from pathlib import Path
from collections import Counter

BASE = Path("output/Fig5_E8_F1_analysis")

HITS = Path("input_files/Fig5_E8_F1_analysis/F1_annotated_motif_hits.tsv")
B6_H5 = Path("input_files/Fig5_E8_F1_analysis/F1_B6_states.h5")
SPRET_H5 = Path("input_files/Fig5_E8_F1_analysis/F1_SPRET_states.h5")

OUT = BASE / "F1_atacPeaks_variant_hits_ac_metrics_median_disp.tsv"

REGION_LEN = 2048
CENTER_START = 512
CENTER_END = 1536  #Exclusive

@njit(cache=False)
def extend_to_full(mat, smf_pos, full_len=2048):
    n = mat.shape[0]
    out = np.full((n, full_len), -1, dtype=np.int8)

    for r in range(n):
        row = mat[r]
        valid = row != -1
        pos = smf_pos[valid]
        val = row[valid]

        keep = (pos >= 0) & (pos < full_len)
        pos = pos[keep].astype(np.int32)
        val = val[keep]

        if len(pos) == 0:
            continue

        order = np.argsort(pos)
        pos = pos[order]
        val = val[order]

        #Preserve original behavior: interpolate only between
        #Measured positions, leaving the ends as -1.
        for i in range(len(pos) - 1):
            p1, p2 = pos[i], pos[i + 1]

            if p2 <= p1:
                continue

            if val[i] == val[i + 1]:
                out[r, p1:p2] = val[i]
            else:
                mid = (p1 + p2) // 2
                out[r, p1:mid] = val[i]
                out[r, mid:p2] = val[i + 1]

        out[r, pos] = val

    return out

def stretch_around_x(vals, cols, x):
    vals = np.asarray(vals)
    cols = np.asarray(cols)

    if len(cols) < 2:
        return np.array([np.nan, np.nan, np.nan])

    bounds = np.r_[
        -np.inf,
        (cols[:-1] + cols[1:]) / 2,
        np.inf
    ]

    idx = np.searchsorted(bounds, x, side="right") - 1
    idx = np.clip(idx, 0, len(vals) - 1)

    if vals[idx] == 3:
        return np.array([np.nan, np.nan, 0])

    l = idx
    while l > 0 and vals[l - 1] != 3:
        l -= 1

    r = idx
    while r < len(vals) - 1 and vals[r + 1] != 3:
        r += 1

    if l == 0:
        left_edge = cols[0] - (cols[1] - cols[0]) / 2
    else:
        left_edge = (cols[l - 1] + cols[l]) / 2

    if r == len(cols) - 1:
        right_edge = cols[r] + (cols[r] - cols[r - 1]) / 2
    else:
        right_edge = (cols[r] + cols[r + 1]) / 2

    return np.array([
        left_edge,
        right_edge,
        right_edge - left_edge
    ])

def boundary_euclidean(res):
    bounds = res[:, :2]
    n = len(bounds)

    if n < 2:
        return np.nan

    diff = bounds[:, None, :] - bounds[None, :, :]
    dist = np.sqrt(np.sum(diff**2, axis=2))
    
    d = dist[np.triu_indices(n, k=1)]
    return np.median(d)

def prepare_pair(b6_h5, spret_h5, region):
    """
    SMF positions in both files are already 0-based coordinates
    relative to the start of the 2048-bp region.
    """

    b6_pos = np.asarray(
        b6_h5["smf_pos"][region][:],
        dtype=np.int64
    )
    spret_pos = np.asarray(
        spret_h5["smf_pos"][region][:],
        dtype=np.int64
    )

    b6_states = b6_h5["smf_states"][region][:]
    spret_states = spret_h5["smf_states"][region][:]

    if b6_states.shape[1] != len(b6_pos):
        raise ValueError(f"B6 position mismatch: {region}")

    if spret_states.shape[1] != len(spret_pos):
        raise ValueError(f"SPRET position mismatch: {region}")

    #Preserve your original variability filter:
    #State 3 = nucleosome, all other observed states = accessible.
    #A position is variable only if both classes occur among
    #Its observed, nonmissing values.

    b6_bin = np.where(
        b6_states == -1,
        -1,
        (b6_states == 3).astype(np.int8)
    )

    spret_bin = np.where(
        spret_states == -1,
        -1,
        (spret_states == 3).astype(np.int8)
    )

    b6_var = (
        np.any(b6_bin == 0, axis=0) &
        np.any(b6_bin == 1, axis=0)
    )

    spret_var = (
        np.any(spret_bin == 0, axis=0) &
        np.any(spret_bin == 1, axis=0)
    )

    b6_pos = b6_pos[b6_var]
    spret_pos = spret_pos[spret_var]

    b6_states = b6_states[:, b6_var]
    spret_states = spret_states[:, spret_var]

    #Intersect positions AFTER independent variability filtering.
    common_pos = np.intersect1d(b6_pos, spret_pos)

    if len(common_pos) < 2:
        return None

    #Explicitly reorder both matrices into the identical,
    #Ascending common-position order.
    b6_order = np.argsort(b6_pos)
    spret_order = np.argsort(spret_pos)

    b6_pos = b6_pos[b6_order]
    spret_pos = spret_pos[spret_order]

    b6_states = b6_states[:, b6_order]
    spret_states = spret_states[:, spret_order]

    if (
        len(np.unique(b6_pos)) != len(b6_pos) or
        len(np.unique(spret_pos)) != len(spret_pos)
    ):
        raise ValueError(f"Duplicate SMF positions: {region}")

    b6_idx = np.searchsorted(b6_pos, common_pos)
    spret_idx = np.searchsorted(spret_pos, common_pos)

    b6_states = b6_states[:, b6_idx]
    spret_states = spret_states[:, spret_idx]

    return common_pos, b6_states, spret_states

def extract_ac_info(states_mat, smf_pos, motif_center, region_len=2048):

    binary_mat = np.where(
        states_mat == -1,
        -1,
        (states_mat == 3).astype(np.int8)
    )

    extended_mat = extend_to_full(
        binary_mat,
        smf_pos,
        full_len=region_len
    )

    #Original 5–95% nucleosome-occupancy filter.
    #Ignore unmeasured (-1) positions.
    valid = extended_mat != -1
    denom = valid.sum(axis=1)
    nuc = (extended_mat == 1).sum(axis=1)

    row_mean = np.full(states_mat.shape[0], np.nan)

    np.divide(
        nuc,
        denom,
        out=row_mean,
        where=denom > 0
    )

    keep_rows_nuc = (
        np.isfinite(row_mean) &
        (row_mean >= 0.05) &
        (row_mean <= 0.95)
    )

    states_mat = states_mat[keep_rows_nuc]
    extended_mat = extended_mat[keep_rows_nuc]

    if states_mat.shape[0] == 0:
        return None

    #Original ±7.5-bp accessibility window.
    center_idx = np.where(
        np.abs(np.arange(region_len) - motif_center) <= 7.5
    )[0]

    if len(center_idx) == 0:
        return None

    center_mat = extended_mat[:, center_idx]

    #The entire motif-center window must be measured
    #And accessible (binary state 0).
    keep_rows = (
        np.all(center_mat != -1, axis=1) &
        np.all(center_mat == 0, axis=1)
    )

    n_accessible = int(keep_rows.sum())
    n_reads = states_mat.shape[0]

    #Do not calculate the more expensive metrics for
    #Motifs that fail the requested read threshold.
    if n_accessible < 20:
        return {
            "n": n_reads,
            "n_accessible": n_accessible
        }

    res = np.array([
        stretch_around_x(
            states_mat[i],
            smf_pos,
            motif_center
        )
        for i in np.where(keep_rows)[0]
    ])

    widths = res[:, 2]
    med_width = (
        np.nanmedian(widths)
        if np.any(np.isfinite(widths))
        else np.nan
    )

    boundary_dispersion = (
        boundary_euclidean(res)
        if len(res) > 1
        else np.nan
    )

    return {
        "n": n_reads,
        "n_accessible": n_accessible,
        "freq": n_accessible / n_reads,
        "median_width": med_width,
        "boundary_dispersion": boundary_dispersion
    }
#Load and filter hits
df = pd.read_csv(HITS, sep="\t")

print(f"All hits: {len(df):,}")

#Retain only SNV and indel motifs, irrespective of case.
#This also matches compound labels containing these strings.
variant_mask = (
    df["motif_variant"]
    .fillna("")
    .astype(str)
    .str.contains(r"snv|indel", case=False, regex=True)
)

df = df.loc[variant_mask].copy()

print(f"SNV/indel hits: {len(df):,}")

#Extract the absolute mm10 start of each 2048-bp region.
region_start = (
    df["region_id"]
    .str.extract(r":(\d+)-\d+$")[0]
    .astype(np.int64)
)

#motif_start and motif_end are absolute mm10 positions,
#Assuming BED-style half-open [start, end) intervals.
df["motif_start"] = pd.to_numeric(
    df["motif_start"], errors="raise"
).astype(np.int64)

df["motif_end"] = pd.to_numeric(
    df["motif_end"], errors="raise"
).astype(np.int64)

relative_start = df["motif_start"] - region_start
relative_end = df["motif_end"] - region_start

#Require the entire motif to lie within the middle
#1024 bp, not merely its midpoint.
center_mask = (
    (relative_start >= CENTER_START) &
    (relative_end <= CENTER_END) &
    (relative_end > relative_start)
)

df = df.loc[center_mask].copy()

print(f"Variant hits entirely within central 1024 bp: {len(df):,}")
#Calculate metrics for every remaining hit
results = []
counts = Counter()

with h5py.File(B6_H5, "r") as b6_h5, \
     h5py.File(SPRET_H5, "r") as spret_h5:

    b6_regions = b6_h5["smf_states"]
    spret_regions = spret_h5["smf_states"]

    for region, region_hits in tqdm(
        df.groupby("region_id", sort=False),
        total=df["region_id"].nunique(),
        desc="Processing regions"
    ):

        if (
            region not in b6_regions or
            region not in spret_regions or
            region not in b6_h5["smf_pos"] or
            region not in spret_h5["smf_pos"]
        ):
            counts["missing_region"] += len(region_hits)
            continue

        pair = prepare_pair(
            b6_h5,
            spret_h5,
            region
        )

        if pair is None:
            counts["fewer_than_two_common_positions"] += len(region_hits)
            continue

        smf_pos, b6_states, spret_states = pair

        start = int(region.split(":")[1].split("-")[0])

        for _, row in region_hits.iterrows():

            #Motif center in the same 2048-relative coordinates
            #As smf_pos. Subtract 1 because motif_end is exclusive.
            motif_center = (
                (row["motif_start"] + row["motif_end"] - 1) / 2
                - start
            )

            b6 = extract_ac_info(
                b6_states,
                smf_pos,
                motif_center,
                region_len=REGION_LEN
            )

            spret = extract_ac_info(
                spret_states,
                smf_pos,
                motif_center,
                region_len=REGION_LEN
            )

            if b6 is None or spret is None:
                counts["no_reads_after_occupancy_filter"] += 1
                continue

            #Require >=20 ACCESSIBLE reads per allele,
            #After the original 5–95% read filter.
            if (
                b6["n_accessible"] < 20 or
                spret["n_accessible"] < 20
            ):
                counts["fewer_than_20_accessible"] += 1
                continue

            #Preserve EVERY original hit-table column.
            result = row.to_dict()

            #Append new metrics.
            result.update({
                "n_b6": b6["n"],
                "n_spret": spret["n"],

                "n_accessible_b6": b6["n_accessible"],
                "n_accessible_spret": spret["n_accessible"],

                "perc_ac_b6": b6["freq"],
                "perc_ac_spret": spret["freq"],
                "delta_perc_ac": spret["freq"] - b6["freq"],

                "median_width_b6": b6["median_width"],
                "median_width_spret": spret["median_width"],
                "delta_median_width": (
                    spret["median_width"] - b6["median_width"]
                ),

                "boundary_dispersion_b6": b6["boundary_dispersion"],
                "boundary_dispersion_spret": spret["boundary_dispersion"],
                "delta_boundary_dispersion": (
                    spret["boundary_dispersion"]
                    - b6["boundary_dispersion"]
                )
            })

            results.append(result)
            counts["retained"] += 1
#Write output
metric_cols = [
    "n_b6",
    "n_spret",
    "n_accessible_b6",
    "n_accessible_spret",
    "perc_ac_b6",
    "perc_ac_spret",
    "delta_perc_ac",
    "median_width_b6",
    "median_width_spret",
    "delta_median_width",
    "boundary_dispersion_b6",
    "boundary_dispersion_spret",
    "delta_boundary_dispersion"
]

out = pd.DataFrame(
    results,
    columns=list(df.columns) + metric_cols
)

out = out.replace([np.inf, -np.inf], np.nan)

#Do not drop otherwise eligible motifs just because a
#Width or dispersion metric is undefined.
out.to_csv(OUT, sep="\t", index=False)

print("\n=== SUMMARY ===")
print(f"Variant hits in central window: {len(df):,}")
print(f"Missing region: {counts['missing_region']:,}")
print(
    "Fewer than 2 shared variable SMF positions:",
    f"{counts['fewer_than_two_common_positions']:,}"
)
print(
    "No reads after occupancy filtering:",
    f"{counts['no_reads_after_occupancy_filter']:,}"
)
print(
    "Fewer than 20 accessible reads in either allele:",
    f"{counts['fewer_than_20_accessible']:,}"
)
print(f"Retained motifs: {counts['retained']:,}")
print(f"Retained unique regions: {out['region_id'].nunique():,}")
print(f"Output: {OUT}")
