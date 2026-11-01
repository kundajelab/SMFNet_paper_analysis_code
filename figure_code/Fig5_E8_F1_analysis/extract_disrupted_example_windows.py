#Run from the repository root: python figure_code/Fig5_E8_F1_analysis/extract_disrupted_example_windows.py
from pathlib import Path
Path("output/Fig5_E8_F1_analysis").mkdir(parents=True, exist_ok=True)

import os
import h5py
import numpy as np
import pandas as pd
from tqdm import tqdm
#Input / output
BASE = "output/Fig5_E8_F1_analysis"

METRICS = f"{BASE}/F1_atacPeaks_variant_hits_ac_metrics_median_disp.tsv"
B6_H5 = 'input_files/Fig5_E8_F1_analysis/F1_B6_states.h5'
SPRET_H5 = 'input_files/Fig5_E8_F1_analysis/F1_SPRET_states.h5'

OUT_H5 = f"{BASE}/F1_negative_delta_expanded_400bp_windows.h5"
OUT_TSV = f"{BASE}/F1_negative_delta_windows.tsv"
DIAGNOSTICS = f"{BASE}/F1_negative_delta_filter_diagnostics.tsv"

os.makedirs(os.path.dirname(OUT_H5), exist_ok=True)

WINDOW_SIZE = 400
MAX_DELTA = -0.2
MIN_POS = 15
MAX_MISSING = 0.35
MAX_NUC_FRAC = 0.95
MIN_READS = 50
#Select variant hits
df = pd.read_csv(METRICS, sep="\t")
df.columns = df.columns.str.strip()

for col in ("motif_start", "motif_end", "delta_perc_ac"):
    df[col] = pd.to_numeric(df[col], errors="coerce")

df = df.dropna(subset=["region_id", "motif_start", "motif_end", "delta_perc_ac"])
df = df[df["delta_perc_ac"] <= MAX_DELTA].copy()

#Region-relative motif midpoint determines the window center.
region_coords = df["region_id"].str.extract(r"^chr[^:]+:(\d+)-(\d+)$")
df = df.loc[region_coords[0].notna()].copy()

df["region_start"] = region_coords.loc[df.index, 0].astype(int)
df["region_end"] = region_coords.loc[df.index, 1].astype(int)
df["motif_start"] = df["motif_start"].astype(int)
df["motif_end"] = df["motif_end"].astype(int)

df["motif_center"] = (df["motif_start"] + df["motif_end"]) // 2
df["relative_motif_center"] = df["motif_center"] - df["region_start"]
df["relative_window_start"] = df["relative_motif_center"] - WINDOW_SIZE // 2
df["relative_window_end"] = df["relative_window_start"] + WINDOW_SIZE
df["window_start"] = df["region_start"] + df["relative_window_start"]
df["window_end"] = df["window_start"] + WINDOW_SIZE

#Require the complete 400-bp window to fit inside the HDF5 region.
inside = (
    (df["relative_window_start"] >= 0) &
    (df["window_end"] <= df["region_end"])
)
print(f"Hits with delta <= {MAX_DELTA}: {len(df):,}")
print(f"Hits outside the region boundary: {(~inside).sum():,}")
df = df[inside].copy().reset_index(drop=True)
#Only the KLF and OCT4::SOX2 examples shown in Fig. 5i-j and E8h-i.
selected = (
    ((df['region_id'] == 'chr1:46904135-46906183') & (df['window_start'] == 46905014)) |
    ((df['region_id'] == 'chr7:27621268-27623316') & (df['window_start'] == 27622166))
)
df = df.loc[selected].copy()
#HDF5 key handling
def get_region_keys(states, positions, suffixes):
    keys = set(states.keys()) & set(positions.keys())
    mapping = {}

    #Later suffixes are lower priority than earlier ones.
    for suffix in reversed(suffixes):
        for key in sorted(keys):
            if suffix:
                if key.endswith(suffix):
                    mapping[key[:-len(suffix)]] = key
            elif "|" not in key:
                mapping[key] = key

    return mapping
#Expand complete regions by nearest-neighbor calls
def expand_all_reads(mat, pos):
    if mat.ndim != 2 or mat.shape[1] != len(pos):
        raise ValueError("Incompatible SMF matrix and position dimensions")
    if len(pos) == 0:
        return None

    order = np.argsort(pos)
    pos, mat = pos[order], mat[:, order]

    pos, unique_idx = np.unique(pos, return_index=True)
    mat = mat[:, unique_idx]

    if not np.all(np.isin(mat, [-1, 0, 1, 2, 3])):
        raise ValueError("Unexpected SMF state")

    binary = np.full(mat.shape, -1, dtype=np.int8)
    binary[np.isin(mat, [0, 1, 2])] = 0
    binary[mat == 3] = 1

    grid_start = int(np.floor(pos.min()))
    full = np.arange(grid_start, int(np.ceil(pos.max())) + 1)

    expanded = np.full((len(mat), len(full)), -1, dtype=np.int8)

    for i, read in enumerate(binary):
        called = read != -1
        if not np.any(called):
            continue

        called_pos = pos[called]
        called_states = read[called]

        right = np.searchsorted(called_pos, full, side="left")
        left = np.clip(right - 1, 0, len(called_pos) - 1)
        right = np.clip(right, 0, len(called_pos) - 1)

        use_right = (
            np.abs(called_pos[right] - full) <
            np.abs(full - called_pos[left])
        )
        expanded[i] = called_states[np.where(use_right, right, left)]

    return {
        "expanded": expanded,
        "grid_start": grid_start,
        "positions": pos,
        "binary": binary
    }

def load_region(states, positions, key):
    mat = states[key][:]
    pos = positions[key][:].astype(float)
    return expand_all_reads(mat, pos)
#Filter one 400-bp window
def extract_window(data, start, end, shared):
    width = end - start
    empty = np.empty((0, width), dtype=np.int8)

    counts = {
        "raw_reads": 0,
        "shared_positions": len(shared),
        "failed_min_positions": len(shared) < MIN_POS,
        "failed_missing": 0,
        "failed_no_calls": 0,
        "failed_occupancy": 0,
        "retained_reads": 0
    }

    if data is None:
        counts["failed_min_positions"] = True
        return empty, counts

    binary = data["binary"]
    pos = data["positions"]
    counts["raw_reads"] = len(binary)

    if len(shared) < MIN_POS:
        return empty, counts

    indices = np.searchsorted(pos, shared)
    if np.any(indices >= len(pos)) or not np.array_equal(pos[indices], shared):
        raise ValueError("Shared positions cannot be aligned")

    missing_fraction = (binary[:, indices] == -1).mean(axis=1)
    keep_missing = missing_fraction <= MAX_MISSING
    counts["failed_missing"] = int((~keep_missing).sum())

    if not np.any(keep_missing):
        return empty, counts

    #Interpolation used the COMPLETE region.
    #Clip at its edges when extracting this window.
    target_pos = np.arange(start, end)
    source_idx = np.clip(
        target_pos - data["grid_start"],
        0,
        data["expanded"].shape[1] - 1
    )

    mat = data["expanded"][:, source_idx][keep_missing]

    has_calls = np.any(binary[keep_missing] != -1, axis=1)
    counts["failed_no_calls"] = int((~has_calls).sum())
    mat = mat[has_calls]

    if not len(mat):
        return empty, counts

    keep_occupancy = (mat == 1).mean(axis=1) <= MAX_NUC_FRAC
    counts["failed_occupancy"] = int((~keep_occupancy).sum())
    mat = mat[keep_occupancy]

    counts["retained_reads"] = len(mat)
    return mat, counts
#Process regions and write expanded HDF5
output_rows = []
diagnostics = []
window_index = 0

with (
    h5py.File(B6_H5, "r") as b6_file,
    h5py.File(SPRET_H5, "r") as spret_file,
    h5py.File(OUT_H5, "w") as output_h5
):
    b6_states, b6_pos = b6_file["smf_states"], b6_file["smf_pos"]
    spret_states, spret_pos = spret_file["smf_states"], spret_file["smf_pos"]

    b6_keys = get_region_keys(b6_states, b6_pos, ("", "|BL6", "|B6"))
    spret_keys = get_region_keys(spret_states, spret_pos, ("", "|SPRET"))

    common_regions = set(b6_keys) & set(spret_keys)
    df = df[df["region_id"].isin(common_regions)].copy()

    print(f"Motif windows with both actual sources: {len(df):,}")
    print(f"Regions: {df['region_id'].nunique():,}")

    windows_group = output_h5.create_group("windows")

    output_h5.attrs["window_size"] = WINDOW_SIZE
    output_h5.attrs["max_delta_perc_ac"] = MAX_DELTA
    output_h5.attrs["min_common_positions"] = MIN_POS
    output_h5.attrs["max_missing"] = MAX_MISSING
    output_h5.attrs["max_nuc_frac"] = MAX_NUC_FRAC
    output_h5.attrs["conditions"] = "actual_b6,actual_spret"
    output_h5.attrs["interpolation"] = (
        "Full-region nearest-neighbor interpolation with edge extrapolation"
    )

    for region, region_df in tqdm(
        df.groupby("region_id", sort=False),
        total=df["region_id"].nunique(),
        desc="Expanding regions"
    ):
        #Expand both complete regions only once.
        data_b6 = load_region(b6_states, b6_pos, b6_keys[region])
        data_spret = load_region(spret_states, spret_pos, spret_keys[region])

        common_positions = (
            np.intersect1d(data_b6["positions"], data_spret["positions"])
            if data_b6 is not None and data_spret is not None
            else np.array([], dtype=float)
        )

        for _, hit in region_df.iterrows():
            start = int(hit["relative_window_start"])
            end = int(hit["relative_window_end"])

            shared = common_positions[
                (common_positions >= start) &
                (common_positions < end)
            ]

            mat_b6, counts_b6 = extract_window(data_b6, start, end, shared)
            mat_spret, counts_spret = extract_window(data_spret, start, end, shared)

            diagnostic = {
                "region_id": region,
                "motif_start": int(hit["motif_start"]),
                "motif_end": int(hit["motif_end"]),
                "motif": hit["motif"],
                "delta_perc_ac": hit["delta_perc_ac"],
                "window_start": int(hit["window_start"]),
                "window_end": int(hit["window_end"]),
                "shared_positions": len(shared),
                **{f"b6_{k}": v for k, v in counts_b6.items()},
                **{f"spret_{k}": v for k, v in counts_spret.items()}
            }

            if len(shared) < MIN_POS:
                diagnostic["status"] = "insufficient_shared_positions"
                diagnostics.append(diagnostic)
                continue

            #Retain the original downstream minimum without computing unused clusters.
            if min(len(mat_b6), len(mat_spret)) < MIN_READS:
                diagnostic["status"] = "fewer_than_50_reads"
                diagnostics.append(diagnostic)
                continue

            #One HDF5 group per retained variant-centered window.
            window_id = f"window_{window_index:06d}"
            g = windows_group.create_group(window_id)

            for attr in (
                "region_id", "motif", "strand", "pattern_id",
                "motif_variant", "source_allele", "opposite_allele"
            ):
                if attr in hit.index and pd.notna(hit[attr]):
                    g.attrs[attr] = str(hit[attr])

            for attr in (
                "motif_start", "motif_end", "motif_center",
                "window_start", "window_end", "region_start",
                "relative_motif_center", "relative_window_start",
                "relative_window_end"
            ):
                g.attrs[attr] = int(hit[attr])

            for attr in ("delta_perc_ac", "perc_ac_b6", "perc_ac_spret"):
                if attr in hit.index and pd.notna(hit[attr]):
                    g.attrs[attr] = float(hit[attr])

            g.attrs["n_actual_b6"] = len(mat_b6)
            g.attrs["n_actual_spret"] = len(mat_spret)

            g.create_dataset(
                "relative_positions",
                data=np.arange(start, end, dtype=np.int32)
            )
            g.create_dataset(
                "genomic_positions",
                data=np.arange(
                    int(hit["window_start"]),
                    int(hit["window_end"]),
                    dtype=np.int64
                )
            )
            g.create_dataset(
                "common_positions",
                data=shared
            )

            g.create_dataset(
                "actual_b6",
                data=mat_b6,
                compression="gzip",
                compression_opts=4
            )
            g.create_dataset(
                "actual_spret",
                data=mat_spret,
                compression="gzip",
                compression_opts=4
            )

            result = hit.to_dict()
            result["h5_window_id"] = window_id
            result["common_positions_400bp"] = len(shared)
            result["actual_b6_n_reads_400bp"] = len(mat_b6)
            result["actual_spret_n_reads_400bp"] = len(mat_spret)
            output_rows.append(result)

            diagnostic["status"] = "retained"
            diagnostic["h5_window_id"] = window_id
            diagnostics.append(diagnostic)

            window_index += 1
#Save TSVs
pd.DataFrame(output_rows).to_csv(
    OUT_TSV, sep="\t", index=False, na_rep="NA"
)

pd.DataFrame(diagnostics).to_csv(
    DIAGNOSTICS, sep="\t", index=False, na_rep="NA"
)

print(f"\nRetained 400-bp windows: {window_index:,}")
print(f"Expanded matrices: {OUT_H5}")
print(f"Window annotations: {OUT_TSV}")
print(f"Filtering diagnostics: {DIAGNOSTICS}")
