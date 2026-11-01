#Run from the repository root: python figure_code/Fig5_E8_F1_analysis/compute_accessibility_windows.py
from pathlib import Path
Path("output/Fig5_E8_F1_analysis").mkdir(parents=True, exist_ok=True)

import h5py
import numpy as np
import pandas as pd
from tqdm import tqdm
#Input and output
WINDOWS = "output/Fig5_E8_F1_analysis/F1_atacPeaks_trainRegions_variant_motif_windows.tsv"

B6_H5 = "input_files/Fig5_E8_F1_analysis/F1_B6_states.h5"

SPRET_H5 = "input_files/Fig5_E8_F1_analysis/F1_SPRET_states.h5"

PRED_H5 = "input_files/Fig5_E8_F1_analysis/F1_predicted_common_states.h5"

OUT = "output/Fig5_E8_F1_analysis/F1_ac_pred_act.tsv"

OUT_H5 = "output/Fig5_E8_F1_analysis/F1_expanded_400bp_windows.h5"

DIAGNOSTICS_OUT = "output/Fig5_E8_F1_analysis/F1_ac_pred_act_filter_diagnostics.tsv"

SUMMARY_OUT = "output/Fig5_E8_F1_analysis/F1_ac_pred_act_filter_summary.tsv"

#Parameters
WINDOW_SIZE = 400
MIN_POS = 15
MAX_MISSING = 0.35
MAX_NUC_FRAC = 0.95
MOTIF_WIDTH = 15

#Predicted smf_pos coordinates are 1 bp behind
#The actual coordinates in the tested examples.
#Apply this offset before interpolation and
#Determining positions shared by all conditions.
PRED_POSITION_OFFSET = 1

CONDITIONS = (
    "actual_b6",
    "actual_spret",
    "predicted_b6",
    "predicted_spret"
)

METADATA_COLUMNS = (
    "strand",
    "motif",
    "pattern_id",
    "motif_variant",
    "source_allele",
    "source_seq",
    "source_cwm_score",
    "source_cwm_strand",
    "opposite_allele",
    "opposite_seq",
    "opposite_cwm_score",
    "opposite_cwm_strand",
    "ref_seq",
    "ref_cwm_score",
    "ref_cwm_strand",
    "cwm_delta"
)
#HDF5 key handling
def available_keys(states, positions):
    return set(states.keys()) & set(positions.keys())

def region_key_map(keys, suffixes):
    mapping = {}

    for key in sorted(keys):
        for suffix in suffixes:
            if suffix and key.endswith(suffix):
                mapping[key[:-len(suffix)]] = key

            elif not suffix and "|" not in key:
                mapping[key] = key

    return mapping
#Expand reads using all original positions
def expand_all_reads(mat, pos):
    """
    Expand full-region SMF reads onto a 1-bp grid.

    States:
        0, 1, 2 -> 0 (accessible)
        3       -> 1 (occupied)
        -1      -> missing

    Interpolate across missing measurements using
    all available full-region positions.

    Extrapolate each read's nearest valid measurement
    toward the edges of the full HDF5 position grid.

    Keep original binary calls for the subsequent
    shared-position missingness filter.
    """
    if mat.ndim != 2 or mat.shape[1] != len(pos):
        raise ValueError(
            "SMF matrix and position array have "
            "incompatible dimensions"
        )

    if len(pos) == 0:
        return None

    #Sort all original measured positions.
    order = np.argsort(pos)
    pos = pos[order]
    mat = mat[:, order].copy()

    #Remove duplicate positions.
    pos, unique_idx = np.unique(
        pos,
        return_index=True
    )
    mat = mat[:, unique_idx]

    if not np.all(np.isin(mat, [-1, 0, 1, 2, 3])):
        raise ValueError(
            "Unexpected SMF states detected"
        )

    #Binary conversion, preserving missing calls.
    binary = np.full(
        mat.shape,
        -1,
        dtype=np.int8
    )

    binary[np.isin(mat, [0, 1, 2])] = 0
    binary[mat == 3] = 1

    #Full region-relative 1-bp grid.
    grid_start = int(np.floor(pos.min()))
    grid_end = int(np.ceil(pos.max())) + 1

    full = np.arange(grid_start, grid_end)

    expanded = np.full(
        (len(binary), len(full)),
        -1,
        dtype=np.int8
    )

    #Each read is interpolated independently.
    for i, read in enumerate(binary):

        called = read != -1

        if not np.any(called):
            continue

        called_pos = pos[called]
        called_states = read[called]

        right = np.searchsorted(
            called_pos,
            full,
            side="left"
        )

        left = np.clip(
            right - 1,
            0,
            len(called_pos) - 1
        )

        right = np.clip(
            right,
            0,
            len(called_pos) - 1
        )

        use_right = (
            np.abs(called_pos[right] - full)
            < np.abs(full - called_pos[left])
        )

        expanded[i] = called_states[
            np.where(use_right, right, left)
        ]

    return {
        "expanded": expanded,
        "grid_start": grid_start,
        "positions": pos,
        "binary": binary
    }

def load_condition(states, positions, key, offset=0):
    """
    Load and expand one complete HDF5 region.

    Actual conditions:    offset=0
    Predicted conditions: offset=1

    The offset is applied before interpolation,
    shared-position matching or window extraction.
    """
    pos = (
        positions[key][:].astype(float)
        + offset
    )

    mat = states[key][:]

    return expand_all_reads(mat, pos)
#Find positions shared by all four conditions
def get_common_positions(region_data):
    """
    Intersect original measured positions across
    all four conditions AFTER the predicted +1 shift.

    These positions are used only for:
        - MIN_POS
        - MAX_MISSING

    Full-region interpolation still uses every
    available position from each condition.
    """
    if any(
        region_data[condition] is None
        for condition in CONDITIONS
    ):
        return np.array([], dtype=float)

    common = region_data[
        "actual_b6"
    ]["positions"].copy()

    for condition in CONDITIONS[1:]:
        common = np.intersect1d(
            common,
            region_data[condition]["positions"]
        )

        if len(common) == 0:
            break

    return common
#Trim to 400 bp and filter on common positions
def extract_window(data, start, end, common_positions):
    """
    Extract a 400-bp window from fully interpolated
    reads.

    Filters:
        1. At least 15 common measured positions
        2. At most 35% missing original calls at
           those common positions
        3. At least one original valid full-region call
        4. At most 95% occupancy in the expanded window
    """
    empty = np.empty(
        (0, end - start),
        dtype=np.int8
    )

    counts = {
        "raw_reads": 0,
        "common_positions_400bp": 0,
        "failed_min_positions": False,
        "failed_missing": 0,
        "after_missing": 0,
        "failed_no_calls": 0,
        "after_valid_calls": 0,
        "failed_occupancy": 0,
        "after_occupancy": 0
    }

    if data is None:
        counts["failed_min_positions"] = True
        return empty, counts

    pos = data["positions"]
    binary = data["binary"]
    expanded = data["expanded"]
    grid_start = data["grid_start"]

    counts["raw_reads"] = len(binary)

    #Only positions shared by all four conditions
    #Inside this motif-centered 400-bp window.
    shared = common_positions[
        (common_positions >= start)
        & (common_positions < end)
    ]

    n_shared = len(shared)

    counts["common_positions_400bp"] = n_shared

    if n_shared < MIN_POS:
        counts["failed_min_positions"] = True
        return empty, counts

    #Recover original columns corresponding to
    #The shared positions in this condition.
    shared_indices = np.searchsorted(pos, shared)

    #Check bounds before indexing.
    if np.any(shared_indices >= len(pos)):
        raise ValueError(
            "Common-position alignment failed: "
            "position outside original array"
        )

    if not np.array_equal(
        pos[shared_indices],
        shared
    ):
        raise ValueError(
            "Common-position alignment failed"
        )

    original_shared = binary[:, shared_indices]

    #Calculate missingness before interpolation,
    #Using the same positions for all conditions.
    missing_fraction = (
        original_shared == -1
    ).mean(axis=1)

    keep_missing = (
        missing_fraction <= MAX_MISSING
    )

    counts["failed_missing"] = int(
        (~keep_missing).sum()
    )

    counts["after_missing"] = int(
        keep_missing.sum()
    )

    if not np.any(keep_missing):
        return empty, counts

    #Trim the previously interpolated full-region
    #Reads to the required 400-bp window.
    #If the window extends beyond the complete
    #HDF5 grid, extrapolate its nearest boundary.
    target_positions = np.arange(start, end)

    source_indices = np.clip(
        target_positions - grid_start,
        0,
        expanded.shape[1] - 1
    )

    window = expanded[:, source_indices]

    #Apply missingness filter.
    window = window[keep_missing]

    #Remove reads that never had a valid original
    #Call anywhere in the complete HDF5 region.
    has_calls = np.any(
        binary[keep_missing] != -1,
        axis=1
    )

    counts["failed_no_calls"] = int(
        (~has_calls).sum()
    )

    window = window[has_calls]

    counts["after_valid_calls"] = len(window)

    if len(window) == 0:
        return empty, counts

    #Occupancy filter on the expanded 400-bp matrix.
    occupied_fraction = (
        window == 1
    ).mean(axis=1)

    keep_occupancy = (
        occupied_fraction <= MAX_NUC_FRAC
    )

    counts["failed_occupancy"] = int(
        (~keep_occupancy).sum()
    )

    window = window[keep_occupancy]

    counts["after_occupancy"] = len(window)

    return window, counts
#Central 15-bp accessibility
def accessibility(mat, start, center):
    """
    A read is accessible if all 15 central bases
    are zero after full-region interpolation.

    Returns:
        fraction accessible
        number of observed reads
        number of accessible reads
    """
    if len(mat) == 0:
        return np.nan, 0, 0

    motif_start = center - MOTIF_WIDTH // 2

    left = motif_start - start
    right = left + MOTIF_WIDTH

    if left < 0 or right > mat.shape[1]:
        return np.nan, 0, 0

    sub = mat[:, left:right]

    observed = np.all(
        sub != -1,
        axis=1
    )

    accessible = (
        observed
        & np.all(sub == 0, axis=1)
    )

    n_observed = int(observed.sum())
    n_accessible = int(accessible.sum())

    fraction = (
        n_accessible / n_observed
        if n_observed > 0
        else np.nan
    )

    return fraction, n_observed, n_accessible

def motif_coverage(mat, start, center):
    """
    Whole-region eligibility check.

    At least one retained read in each condition
    must cover the full interpolated central 15 bp.
    """
    if len(mat) == 0:
        return False

    motif_start = center - MOTIF_WIDTH // 2

    left = motif_start - start
    right = left + MOTIF_WIDTH

    if left < 0 or right > mat.shape[1]:
        return False

    sub = mat[:, left:right]

    return bool(
        np.any(
            np.all(sub != -1, axis=1)
        )
    )
#Read and filter motif table
df = pd.read_csv(
    WINDOWS,
    sep="\t"
)

df.columns = df.columns.str.strip()

df["motif_start"] = pd.to_numeric(
    df["motif_start"]
)

df["motif_end"] = pd.to_numeric(
    df["motif_end"]
)

#Keep chromosomes 1, 3 and 6.
region_info = df["region_id"].str.extract(
    r"^(chr(?:1|3|6)):(\d+)-(\d+)$"
)

df = df.loc[
    region_info[0].notna()
].copy()

df["region_start"] = (
    region_info.loc[df.index, 1].astype(int)
)

#Region-relative motif center.
df["relative_motif_center"] = (
    (
        df["motif_start"]
        + df["motif_end"]
    ) // 2
    - df["region_start"]
).astype(int)

#400-bp motif-centered windows.
df["relative_window_start"] = (
    df["relative_motif_center"] - 200
)

df["relative_window_end"] = (
    df["relative_window_start"] + WINDOW_SIZE
)

#Require the entire window to lie within 512–1536.
df = df[
    (df["relative_window_start"] >= 512)
    & (df["relative_window_end"] <= 1536)
].copy()

#Retain SNVs and indels.
if "motif_variant" in df.columns:
    df = df[
        df["motif_variant"]
        .str.lower()
        .isin(["snv", "indel"])
    ].copy()
#Initialize outputs
rows = []

filter_diagnostics = []

excluded_regions = 0
excluded_min_positions = 0
excluded_motif_coverage = 0
#Open HDF5 files
with (
    h5py.File(B6_H5, "r") as b6_file,
    h5py.File(SPRET_H5, "r") as spret_file,
    h5py.File(PRED_H5, "r") as predicted_file,
    h5py.File(OUT_H5, "w") as output_h5
):
    b6_states = b6_file["smf_states"]
    b6_positions = b6_file["smf_pos"]

    spret_states = spret_file["smf_states"]
    spret_positions = spret_file["smf_pos"]

    predicted_states = predicted_file["smf_states"]
    predicted_positions = predicted_file["smf_pos"]

    #Actual dataset keys.
    actual_b6_keys = region_key_map(
        available_keys(
            b6_states,
            b6_positions
        ),
        ("", "|BL6", "|B6")
    )

    actual_spret_keys = region_key_map(
        available_keys(
            spret_states,
            spret_positions
        ),
        ("", "|SPRET")
    )

    #Predicted dataset keys.
    prediction_keys = available_keys(
        predicted_states,
        predicted_positions
    )

    named_b6 = region_key_map(
        prediction_keys,
        ("|BL6", "|B6")
    )

    named_spret = region_key_map(
        prediction_keys,
        ("|SPRET",)
    )

    ref_keys = region_key_map(
        prediction_keys,
        ("|ref",)
    )

    alt_keys = region_key_map(
        prediction_keys,
        ("|alt",)
    )

    #Prefer named parental alleles.
    #Fallback: ref=B6 and alt=SPRET.
    predicted_b6_keys = dict(ref_keys)
    predicted_b6_keys.update(named_b6)

    predicted_spret_keys = dict(alt_keys)
    predicted_spret_keys.update(named_spret)

    if ref_keys or alt_keys:
        print(
            "WARNING: ref/alt fallback assumes "
            "ref=B6 and alt=SPRET."
        )

    #Require all four datasets.
    common_regions = (
        set(actual_b6_keys)
        & set(actual_spret_keys)
        & set(predicted_b6_keys)
        & set(predicted_spret_keys)
    )

    df = df[
        df["region_id"].isin(common_regions)
    ].copy()

    print(
        f"Found {len(df):,} motifs in "
        f"{df['region_id'].nunique():,} regions"
    )

    print(
        f"Applying +{PRED_POSITION_OFFSET} bp "
        "to both predicted position arrays"
    )

    windows_group = output_h5.create_group(
        "windows"
    )

    output_h5.attrs["window_size"] = WINDOW_SIZE
    output_h5.attrs["min_common_positions"] = MIN_POS
    output_h5.attrs["max_missing"] = MAX_MISSING
    output_h5.attrs["max_nuc_frac"] = MAX_NUC_FRAC

    #Document the coordinate correction.
    output_h5.attrs["actual_position_offset"] = 0

    output_h5.attrs["predicted_position_offset"] = (
        PRED_POSITION_OFFSET
    )

    output_h5.attrs["interpolation"] = (
        "Full-region nearest-neighbor interpolation "
        "with edge extrapolation"
    )

    output_h5.attrs["filter_positions"] = (
        "Intersection of original positions across "
        "all four conditions after predicted +1 bp shift"
    )

    window_index = 0
    #Process complete regions
    for region, region_df in tqdm(
        df.groupby("region_id", sort=False),
        total=df["region_id"].nunique(),
        desc="Processing regions"
    ):
        #Actual positions remain unchanged.
        #Predicted positions are shifted by +1
        #Before interpolation and shared-position
        #Matching.
        region_data = {
            "actual_b6": load_condition(
                b6_states,
                b6_positions,
                actual_b6_keys[region],
                offset=0
            ),

            "actual_spret": load_condition(
                spret_states,
                spret_positions,
                actual_spret_keys[region],
                offset=0
            ),

            "predicted_b6": load_condition(
                predicted_states,
                predicted_positions,
                predicted_b6_keys[region],
                offset=PRED_POSITION_OFFSET
            ),

            "predicted_spret": load_condition(
                predicted_states,
                predicted_positions,
                predicted_spret_keys[region],
                offset=PRED_POSITION_OFFSET
            )
        }

        #Shared positions in the corrected
        #Coordinate system.
        common_positions = get_common_positions(
            region_data
        )

        #Buffer all motifs until the entire region
        #Has passed the eligibility requirements.
        region_rows = []
        region_windows = []
        region_diagnostics = []

        region_valid = True
        region_failed_positions = False
        region_failed_coverage = False

        for row in region_df.itertuples(index=False):
            center = int(row.relative_motif_center)
            start = int(row.relative_window_start)
            end = int(row.relative_window_end)

            shared_window = common_positions[
                (common_positions >= start)
                & (common_positions < end)
            ]

            matrices = {}
            counts_by_condition = {}
            #Extract and filter all four conditions
            for condition in CONDITIONS:
                mat, counts = extract_window(
                    region_data[condition],
                    start,
                    end,
                    common_positions
                )

                matrices[condition] = mat
                counts_by_condition[condition] = counts

                fraction, n_obs, n_acc = accessibility(
                    mat,
                    start,
                    center
                )

                diagnostic = {
                    "region_id": region,
                    "motif_start": int(row.motif_start),
                    "motif_end": int(row.motif_end),
                    "condition": condition,
                    **counts,
                    "failed_missing_pct": (
                        100 * counts["failed_missing"]
                        / counts["raw_reads"]
                        if counts["raw_reads"] > 0
                        else np.nan
                    ),
                    "retained_pct": (
                        100 * counts["after_occupancy"]
                        / counts["raw_reads"]
                        if counts["raw_reads"] > 0
                        else np.nan
                    ),
                    "motif_observed": n_obs,
                    "motif_accessible": n_acc,
                    "motif_inaccessible": (
                        n_obs - n_acc
                    ),
                    "motif_accessible_pct": (
                        100 * fraction
                        if np.isfinite(fraction)
                        else np.nan
                    )
                }

                region_diagnostics.append(diagnostic)
            #Whole-region eligibility
            position_coverage_ok = (
                len(shared_window) >= MIN_POS
            )

            central_coverage_ok = all(
                motif_coverage(
                    matrices[condition],
                    start,
                    center
                )
                for condition in CONDITIONS
            )

            #Record failures but continue processing
            #Other motifs for complete diagnostics.
            if not position_coverage_ok:
                region_valid = False
                region_failed_positions = True

            if not central_coverage_ok:
                region_valid = False
                region_failed_coverage = True

            if not (
                position_coverage_ok
                and central_coverage_ok
            ):
                continue
            #Motif-level TSV results
            genomic_center = (
                center + int(row.region_start)
            )

            genomic_start = (
                start + int(row.region_start)
            )

            genomic_end = (
                end + int(row.region_start)
            )

            result = {
                "region_id": region,
                "motif_start": row.motif_start,
                "motif_end": row.motif_end,
                "motif_center": genomic_center,
                "window_start": genomic_start,
                "window_end": genomic_end,
                "relative_motif_center": center,
                "relative_window_start": start,
                "relative_window_end": end,
                "common_positions_400bp": len(
                    shared_window
                ),
                "h5_window_id": None
            }

            for col in METADATA_COLUMNS:
                if col in df.columns:
                    result[col] = getattr(row, col)

            for condition, mat in matrices.items():
                fraction, n_obs, n_acc = accessibility(
                    mat,
                    start,
                    center
                )

                result[
                    f"{condition}_perc_accessible"
                ] = fraction

                result[
                    f"{condition}_n"
                ] = n_obs

                result[
                    f"{condition}_n_accessible"
                ] = n_acc

                result[
                    f"{condition}_n_reads_400bp"
                ] = len(mat)

            result["actual_delta_accessible"] = (
                result["actual_spret_perc_accessible"]
                - result["actual_b6_perc_accessible"]
            )

            result["predicted_delta_accessible"] = (
                result["predicted_spret_perc_accessible"]
                - result["predicted_b6_perc_accessible"]
            )

            region_rows.append(result)

            region_windows.append({
                "matrices": matrices,
                "region_id": region,
                "motif_start": int(row.motif_start),
                "motif_end": int(row.motif_end),
                "motif_center": genomic_center,
                "window_start": genomic_start,
                "window_end": genomic_end,
                "relative_window_start": start,
                "relative_window_end": end,
                "relative_motif_center": center,
                "common_positions": shared_window.copy(),
                "motif": str(
                    getattr(row, "motif", "")
                )
            })
        #Save complete diagnostics for this region
        for diagnostic in region_diagnostics:
            diagnostic["region_retained"] = region_valid

            diagnostic["region_failed_positions"] = (
                region_failed_positions
            )

            diagnostic["region_failed_coverage"] = (
                region_failed_coverage
            )

        filter_diagnostics.extend(region_diagnostics)

        if not region_valid:
            excluded_regions += 1

            if region_failed_positions:
                excluded_min_positions += 1

            if region_failed_coverage:
                excluded_motif_coverage += 1

            continue
        #Write retained windows to HDF5
        for result, window_info in zip(
            region_rows,
            region_windows
        ):
            window_id = (
                f"window_{window_index:06d}"
            )

            window_group = windows_group.create_group(
                window_id
            )

            window_group.attrs["region_id"] = (
                window_info["region_id"]
            )

            window_group.attrs["motif"] = (
                window_info["motif"]
            )

            for attr in (
                "motif_start",
                "motif_end",
                "motif_center",
                "window_start",
                "window_end",
                "relative_window_start",
                "relative_window_end",
                "relative_motif_center"
            ):
                window_group.attrs[attr] = (
                    window_info[attr]
                )

            #Record the coordinate convention
            #For the stored matrices.
            window_group.attrs[
                "predicted_position_offset"
            ] = PRED_POSITION_OFFSET

            relative_positions = np.arange(
                window_info["relative_window_start"],
                window_info["relative_window_end"],
                dtype=np.int32
            )

            genomic_positions = np.arange(
                window_info["window_start"],
                window_info["window_end"],
                dtype=np.int64
            )

            window_group.create_dataset(
                "relative_positions",
                data=relative_positions
            )

            window_group.create_dataset(
                "genomic_positions",
                data=genomic_positions
            )

            #These positions are already aligned
            #After applying the predicted +1 shift.
            window_group.create_dataset(
                "common_positions",
                data=window_info["common_positions"]
            )

            #Four expanded read-by-400 matrices,
            #All aligned to the same coordinates.
            for condition, mat in (
                window_info["matrices"].items()
            ):
                window_group.create_dataset(
                    condition,
                    data=mat,
                    dtype=np.int8,
                    compression="gzip",
                    compression_opts=4
                )

            result["h5_window_id"] = window_id

            rows.append(result)

            window_index += 1

#Save accessibility TSV
results = pd.DataFrame(rows)

results.to_csv(
    OUT,
    sep="\t",
    index=False,
    na_rep="NA"
)

print(
    f"Wrote {len(results):,} windows to {OUT}"
)

print(
    f"Wrote expanded matrices to {OUT_H5}"
)
#Save complete filtering diagnostics
diagnostics = pd.DataFrame(filter_diagnostics)

diagnostics.to_csv(
    DIAGNOSTICS_OUT,
    sep="\t",
    index=False,
    na_rep="NA"
)

print(
    f"Saved filtering diagnostics to {DIAGNOSTICS_OUT}"
)
#Summarize read filtering
eligible = diagnostics[
    ~diagnostics["failed_min_positions"]
].copy()

summary_columns = [
    "raw_reads",
    "failed_missing",
    "after_missing",
    "failed_no_calls",
    "after_valid_calls",
    "failed_occupancy",
    "after_occupancy",
    "motif_observed",
    "motif_accessible",
    "motif_inaccessible"
]

if len(eligible) > 0:
    summary = eligible.groupby(
        "condition"
    )[summary_columns].sum()

    summary["missing_loss_pct"] = (
        100
        * summary["failed_missing"]
        / summary["raw_reads"].replace(0, np.nan)
    )

    summary["occupancy_loss_pct"] = (
        100
        * summary["failed_occupancy"]
        / summary["raw_reads"].replace(0, np.nan)
    )

    summary["retained_pct"] = (
        100
        * summary["after_occupancy"]
        / summary["raw_reads"].replace(0, np.nan)
    )

    summary["motif_accessible_pct"] = (
        100
        * summary["motif_accessible"]
        / summary["motif_observed"].replace(0, np.nan)
    )

else:
    summary = pd.DataFrame(
        columns=summary_columns
    )

summary.to_csv(
    SUMMARY_OUT,
    sep="\t",
    na_rep="NA"
)

print("\n========== FILTERING SUMMARY ==========")
print(summary.round(2).to_string())

print("\n========== REGION EXCLUSIONS ==========")
print(f"Excluded regions: {excluded_regions}")

print(
    "Failed shared-position requirement: "
    f"{excluded_min_positions}"
)

print(
    "Failed central motif coverage: "
    f"{excluded_motif_coverage}"
)

print(
    "\nNote: A region can fail both criteria. "
    "Insufficient shared positions can also cause "
    "the coverage check to fail."
)

print(f"\nSaved summary to {SUMMARY_OUT}")
