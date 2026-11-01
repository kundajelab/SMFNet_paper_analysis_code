#Run from the repository root: python figure_code/Fig4_E7_combinedMotifs_insilico/compute_state_profiles.py
from pathlib import Path
from collections import defaultdict
import os

#Keep the original non-AVX512 ordering of tied scores, also on newer cluster CPUs.
os.environ["NPY_DISABLE_CPU_FEATURES"] = ",".join(
    filter(
        None,
        [
            os.environ.get("NPY_DISABLE_CPU_FEATURES", ""),
            "AVX512F,AVX512_SKX,AVX512_CLX,AVX512_CNL,AVX512_ICL",
        ],
    )
)
import numpy as np
import pandas as pd
import h5py
from numba import njit

#NumPy 1.26.4 preserves the original ordering of tied read scores between bins.
if np.__version__ != "1.26.4":
    raise RuntimeError("Use NumPy 1.26.4 to reproduce the original binned profiles.")
INPUT = Path("input_files/Fig4_E7_combinedMotifs_insilico")
OUT = Path("output/Fig4_E7_combinedMotifs_insilico/state_profiles")
OUT.mkdir(parents=True, exist_ok=True)
N_BINS = 15


#Interpolate states between measured sites while leaving unmeasured ends as missing.
@njit(cache=False)
def extend_to_full(mat, smf_pos, window_center=512, full_len=1024):
    n = mat.shape[0]
    out = np.full((n, full_len), -1, dtype=mat.dtype)
    offset = window_center - full_len // 2

    for r in range(n):
        row = mat[r]
        out_row = out[r]

        n_valid = 0
        for i in range(row.shape[0]):
            if row[i] != -1:
                n_valid += 1

        if n_valid == 0:
            continue

        mapped_pos = np.empty(n_valid, dtype=np.int32)
        val = np.empty(n_valid, dtype=mat.dtype)

        j = 0
        for i in range(row.shape[0]):
            if row[i] != -1:
                mp = smf_pos[i] - offset
                if 0 <= mp < full_len:
                    mapped_pos[j] = mp
                    val[j] = row[i]
                    j += 1

        if j == 0:
            continue

        for i in range(j - 1):
            p1, p2 = mapped_pos[i], mapped_pos[i + 1]

            if p2 <= p1:
                continue

            if val[i] == val[i + 1]:
                out_row[p1:p2] = val[i]
            else:
                mid = (p1 + p2) // 2
                out_row[p1:mid] = val[i]
                out_row[mid:p2] = val[i + 1]

        for i in range(j):
            out_row[mapped_pos[i]] = val[i]

    return out


def get_combo_from_region(region_name):
    parts = region_name.split("_", 1)
    return parts[1] if len(parts) == 2 else region_name


#Calculate one minus the state-3 fraction at each position, excluding missing calls.
def mean_openness(X):
    valid = X != -1

    state3_fraction = np.divide(
        (X == 1).sum(axis=0),
        valid.sum(axis=0),
        out=np.full(X.shape[1], np.nan),
        where=valid.sum(axis=0) > 0,
    )

    return 1.0 - state3_fraction


def calculate_profiles(filename, analysis, spacing, selected, bin_profiles=True):
    averages, bins = [], []
    with h5py.File(INPUT / filename, "r") as f:
        groups = defaultdict(list)
        for region in f["smf_pos"]:
            if region.split("_", 2)[2] in selected:
                groups[get_combo_from_region(region)].append(region)
        assert len(groups) == len(selected)
        for combo, regions in sorted(groups.items()):
            reads = []
            for region in regions:
                states = f["smf_states"][region][:]
                pos = f["smf_pos"][region][:]
                X = np.where(states == -1, -1, (states == 3).astype(np.int8))
                X = extend_to_full(X, pos)
                valid = X != -1
                score = np.divide(
                    (X == 1).sum(axis=1),
                    valid.sum(axis=1),
                    out=np.zeros(X.shape[0], dtype=float),
                    where=valid.sum(axis=1) > 0,
                )
                #Retain reads with 5-95% state-3 coverage before pooling backgrounds.
                keep = (score >= 0.05) & (score <= 0.95)
                if np.any(keep):
                    reads.append(X[keep])
            X = np.vstack(reads)
            n_reads, length = X.shape
            for pos, value in enumerate(mean_openness(X), 1):
                averages.append([analysis, spacing, combo, pos, value])
            if not bin_profiles:
                continue
            valid = X != -1
            score = np.divide(
                (X == 1).sum(axis=1),
                valid.sum(axis=1),
                out=np.zeros(n_reads, dtype=float),
                where=valid.sum(axis=1) > 0,
            )
            #Sort reads by state-3 coverage and divide them into 15 approximately equal bins.
            X = X[np.argsort(score)]
            edges = np.linspace(0, n_reads, N_BINS + 1).astype(int)
            for b in range(N_BINS):
                lo, hi = edges[b], edges[b + 1]
                profile = mean_openness(X[lo:hi]) if hi - lo >= 5 else np.full(length, np.nan)
                for pos, value in enumerate(profile, 1):
                    bins.append([analysis, spacing, combo, N_BINS - 1 - b, pos, value])
    return averages, bins


def main():
    averages, bins = [], []
    #Prepare the SOX2-flank spacing panels and the selected enhancer examples.
    tasks = [
        (f"three_motif_spacing_{s}bp_states.h5", "SOX2_flanks", s, {"SOX2_OCT4::SOX2_SOX2"}, True)
        for s in [20, 60]
    ]
    #The grey overlay comes from the original center-only enhancer profiles.
    tasks.extend(
        (
            f"mesc_enhancer_spacing_{s}bp_states.h5",
            "SOX2_flanks",
            s,
            {"noMotif_OCT4::SOX2_noMotif"},
            False,
        )
        for s in [20, 60]
    )
    tasks.append(
        (
            "mesc_enhancer_spacing_50bp_states.h5",
            "enhancer",
            50,
            {
                "noMotif_OCT4::SOX2_noMotif",
                "SP-KLF_OCT4::SOX2_noMotif",
                "OCT4::SOX2_OCT4::SOX2_noMotif",
            },
            True,
        )
    )
    for filename, analysis, spacing, selected, bin_profiles in tasks:
        a, b = calculate_profiles(filename, analysis, spacing, selected, bin_profiles)
        averages.extend(a)
        bins.extend(b)
        print("Profiles finished:", analysis, spacing, flush=True)
    pd.DataFrame(averages, columns=["analysis", "spacing", "combo", "position", "openness"]).to_csv(
        OUT / "average_profiles.tsv", sep="\t", index=False
    )
    pd.DataFrame(
        bins, columns=["analysis", "spacing", "combo", "pseudo_bin", "position", "openness"]
    ).to_csv(OUT / "binned_profiles.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
