#Run from the repository root: python figure_code/Fig4_E7_combinedMotifs_insilico/compute_motif_metrics.py
from pathlib import Path
import h5py
import numpy as np
import pandas as pd
from numba import njit
from scipy.spatial.distance import pdist
from tqdm import tqdm

ROOT = Path("output/Fig4_E7_combinedMotifs_insilico/selected_pairs_baseline")
ROOT.mkdir(parents=True, exist_ok=True)
STATES = Path("input_files/Fig4_E7_combinedMotifs_insilico")
PAIRS = Path(
    "input_files/Fig4_E7_combinedMotifs_insilico/motif_pair_counts_training_within_100bp.txt"
)


@njit(cache=False)
def _add_segment(total, central_nuc, left, right, value):
    if right <= left:
        return total, central_nuc
    total += (right - left) * value
    if value == 1:
        central_nuc += max(0, min(right, 520) - max(left, 505))
    return total, central_nuc


@njit(cache=False)

#Measure state-3 coverage and the accessible stretch around the central motif.
def quant_kernel(states, pos):
    n, k = states.shape
    occupancy = np.empty(n, dtype=np.float64)
    open_mask = np.zeros(n, dtype=np.bool_)
    bounds = np.full((n, 2), np.nan)
    widths = np.full(n, np.nan)
    mapidx = np.empty(k, dtype=np.int64)
    mapped = np.empty(k, dtype=np.int64)
    m = 0
    for j in range(k):
        if 1 <= pos[j] <= 1024:
            mapidx[m] = j
            mapped[m] = int(pos[j])
            m += 1
    center_idx = 0
    for j in range(k - 1):
        if (pos[j] + pos[j + 1]) / 2 <= 512:
            center_idx = j + 1
    for i in range(n):
        if m == 0:
            occupancy[i] = -1.0
            continue
        total = 0
        central_nuc = 0
        value = int(states[i, mapidx[0]] == 3)
        total, central_nuc = _add_segment(total, central_nuc, 0, mapped[0], value)
        for j in range(m - 1):
            p1 = mapped[j]
            p2 = mapped[j + 1]
            if p2 <= p1:
                continue
            v1 = int(states[i, mapidx[j]] == 3)
            v2 = int(states[i, mapidx[j + 1]] == 3)
            if v1 == v2:
                total, central_nuc = _add_segment(total, central_nuc, p1, p2, v1)
            else:
                mid = (p1 + p2) // 2
                total, central_nuc = _add_segment(total, central_nuc, p1, mid, v1)
                total, central_nuc = _add_segment(total, central_nuc, mid, p2, v2)
        value = int(states[i, mapidx[m - 1]] == 3)
        total, central_nuc = _add_segment(total, central_nuc, mapped[m - 1], 1024, value)
        occupancy[i] = total / 1024.0
        #Exclude nearly uniform reads before testing accessibility at the motif center.
        if occupancy[i] < 0.05 or occupancy[i] > 0.95:
            continue
        if central_nuc != 0:
            continue
        open_mask[i] = True
        idx = center_idx
        if states[i, idx] == 3:
            widths[i] = 0
            continue
        l = idx
        r = idx
        while l > 0 and states[i, l - 1] != 3:
            l -= 1
        while r < k - 1 and states[i, r + 1] != 3:
            r += 1
        left = pos[0] - (pos[1] - pos[0]) / 2 if l == 0 else (pos[l - 1] + pos[l]) / 2
        right = pos[r] + (pos[r] - pos[r - 1]) / 2 if r == k - 1 else (pos[r] + pos[r + 1]) / 2
        bounds[i, 0] = left
        bounds[i, 1] = right
        widths[i] = right - left
    return occupancy, open_mask, bounds, widths


#Keep selected flanked constructs and all controls needed for their baselines.
def wanted(region, allowed, motifs):
    _, _, left, center, right = region.split("_", 4)
    is_bab = left == right != "noMotif" and center != "noMotif"
    selected = is_bab and tuple(sorted((left, center))) in allowed
    is_a = left == right == "noMotif" and center in motifs
    is_f = left == right and left in motifs and center == "noMotif"
    is_bg = left == center == right == "noMotif"
    return selected or is_a or is_f or is_bg


def main(spacing):
    pairs = pd.read_csv(PAIRS, sep="\t")
    pairs = pairs[pairs.n_pairs > 1000].copy()
    pairs[["motif_1", "motif_2"]] = pairs[["motif_1", "motif_2"]].replace({"KLF": "SP-KLF"})
    allowed = {tuple(sorted((r.motif_1, r.motif_2))) for r in pairs.itertuples()}
    motifs = set(pairs.motif_1) | set(pairs.motif_2)
    folder = ROOT / "per_spacing"
    folder.mkdir(exist_ok=True)
    dest = folder / f"spacing_{spacing}bp_quantification.txt"
    #Resume from existing per-spacing tables when a previous run has already finished regions.
    rows = pd.read_csv(dest, sep="\t").to_dict("records") if dest.exists() else []
    rows = [row for row in rows if wanted(row["region"], allowed, motifs)]
    done = {row["region"] for row in rows}
    sources = [
        ("three_motif", f"three_motif_spacing_{spacing}bp_states.h5"),
        ("flank_only", f"flank_only_spacing_{spacing}bp_states.h5"),
    ]
    for source, filename in sources:
        with h5py.File(STATES / filename, "r") as h5:
            for region in tqdm(h5["smf_pos"], desc=f"{spacing} bp {source}", mininterval=30):
                if region in done:
                    continue
                bg, _, left, center, right = region.split("_", 4)
                if not wanted(region, allowed, motifs):
                    continue
                states = h5["smf_states"][region][:]
                pos = np.asarray(h5["smf_pos"][region][:], dtype=np.float64)
                occ, opened, bounds, widths = quant_kernel(states, pos)
                #Summarize accessible fraction, median stretch width and pairwise boundary dispersion.
                retained = (occ >= 0.05) & (occ <= 0.95)
                nr, no = int(retained.sum()), int(opened.sum())
                rows.append(
                    dict(
                        spacing=spacing,
                        region=region,
                        bg_region=bg,
                        left=left,
                        center=center,
                        right=right,
                        source=source,
                        med_width=float(np.nanmedian(widths[opened])) if no else np.nan,
                        perc_ac=no / nr if nr else np.nan,
                        boundary_dispersion=float(np.median(pdist(bounds[opened])))
                        if no >= 2
                        else np.nan,
                        n_total=len(occ),
                        n_lost_low=int((occ < 0.05).sum()),
                        n_lost_high=int((occ > 0.95).sum()),
                        n_retained=nr,
                        n_open=no,
                    )
                )
    pd.DataFrame(rows).to_csv(str(dest) + ".tmp", sep="\t", index=False)
    Path(str(dest) + ".tmp").replace(dest)
    print("Wrote", dest, "rows:", len(rows), flush=True)


if __name__ == "__main__":
    for spacing in range(20, 101, 5):
        main(spacing)

    #Combine spacing tables for the baseline, SOX2 and phasing analyses.
    with (ROOT / "combined_quantification.txt").open("w") as output:
        for i, spacing in enumerate(range(20, 101, 5)):
            table = pd.read_csv(
                ROOT / "per_spacing" / f"spacing_{spacing}bp_quantification.txt", sep="\t"
            )
            table.to_csv(output, sep="\t", index=False, header=i == 0)
