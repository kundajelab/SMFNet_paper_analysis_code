#Run from the repository root: python figure_code/Fig1_E3_Codependency_benchmarking/compute_coaccessibility.py
#Calculating co-accessibility from states h5 files. Step 1 out of 3.
#Region coverage thresholds reproduce the selection used for the supplied matrices.
from pathlib import Path
import h5py
import numpy as np
import pandas as pd

#Input and output paths; edit these settings to use recalculated inputs.
COVERAGE_FILE = Path("input_files/Fig1_E3_Codependency_benchmarking/coaccessibility_region_coverage.tsv")
MIN_REGION_READS = 200
MIN_REGION_POSITIONS = 150
OUTPUT_DIR = Path("output/Fig1_E3_Codependency_benchmarking/coaccessibility_matrices")

input_dir = Path("input_files/Fig1_E3_Codependency_benchmarking")
experiment = input_dir / "mESC_SMF-ONT_testRegions_Experiment_states.h5"
model = input_dir / "mESC_SMF-ONT_testRegions_SMFNet_states.h5"


def entropy(x):
    """Shannon entropy of 1D array x"""
    x = x[~np.isnan(x)]
    if len(x) == 0:
        return np.nan
    x_int = x.astype(int)
    probs = np.bincount(x_int) / len(x_int)
    probs = probs[probs > 0]
    return -np.sum(probs * np.log2(probs))

def fast_corr(x, y):
    x = x - x.mean()
    y = y - y.mean()
    return np.sum(x * y) / np.sqrt(np.sum(x*x) * np.sum(y*y))

def fast_mi(x, y):
    xy = x * 2 + y  #maps (0,0)->0, (0,1)->1, (1,0)->2, (1,1)->3
    counts = np.bincount(xy.astype(int), minlength=4).reshape(2,2)
    joint = counts / counts.sum()

    px = joint.sum(axis=1)
    py = joint.sum(axis=0)

    nz = joint > 0
    return np.sum(joint[nz] * np.log2(joint[nz] / (px[:, None] * py[None, :])[nz]))

def signed_mi(x, y, Hx, Hy):
    """Compute mutual information between two 1D arrays, times sign of correlation.
    Returns NaN if one column has no variance."""
    valid = ~np.isnan(x) & ~np.isnan(y)
    xv, yv = x[valid], y[valid]
    if len(xv) == 0:
        return np.nan
    if np.std(xv) == 0 or np.std(yv) == 0:
        return np.nan
    mi = fast_mi(xv, yv)
    nmi = mi / np.sqrt(Hx * Hy) if Hx > 0 and Hy > 0 else 0
    corr = fast_corr(xv, yv)
    return nmi * np.sign(corr)

def compute_signed_mi_matrix(mat):
    """Compute signed MI between all positions (columns) in mat."""
    mat_collapsed = np.where(mat == -1, np.nan, np.where(mat <= 2, 1, 0))
    n_pos = mat_collapsed.shape[1]
    mi_mat = np.zeros((n_pos, n_pos), dtype=np.float32)
    entropies = np.array([entropy(mat_collapsed[:, i]) for i in range(n_pos)])
    for i in range(n_pos):
        xi = mat_collapsed[:, i]
        Hx = entropies[i]
        for j in range(i, n_pos):
            Hy = entropies[j]
            xy = mat_collapsed[:, j]
            mi_mat[i,j] = signed_mi(xi, xy, Hx, Hy)
            mi_mat[j,i] = mi_mat[i,j]
    return mi_mat

def main():
    for path in (COVERAGE_FILE, experiment, model):
        if not path.is_file():
            raise FileNotFoundError(f"Missing input: {path}. Run from the repository root.")
    #Counts describe the original 2048 bp observed-read windows.
    coverage = pd.read_csv(COVERAGE_FILE, sep="\t")
    selected = coverage.loc[
        (coverage["cov"] >= MIN_REGION_READS) &
        (coverage["n_pos"] >= MIN_REGION_POSITIONS)
    ]
    regions = sorted(selected["region"].drop_duplicates().tolist())
    print(f"Selected {len(regions)} regions with at least {MIN_REGION_READS} reads "
          f"and {MIN_REGION_POSITIONS} measured positions.")
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    with h5py.File(experiment, "r") as act_h5, h5py.File(model, "r") as pred_h5:
        for number, region in enumerate(regions, 1):
            pos_act = act_h5["smf_pos"][region][()]
            pos_pred = pred_h5["smf_pos"][region][()]
            mat_act = act_h5["smf_states"][region][()]
            mat_pred = pred_h5["smf_states"][region][()]
            common = np.intersect1d(pos_act, pos_pred)
            act_mask = np.isin(pos_act, common)
            pred_mask = np.isin(pos_pred, common)
            if not np.array_equal(pos_act[act_mask], pos_pred[pred_mask]):
                raise ValueError(f"Experiment/model position order differs: {region}")
            mat_act = mat_act[:, act_mask]
            mat_pred = mat_pred[:, pred_mask]
            #Keep shared sites with at least 50 observed state calls.
            keep = np.sum(mat_act != -1, axis=0) >= 50
            positions = pos_act[act_mask][keep]
            for label, mat in (("ACT", mat_act[:, keep]), ("PRED", mat_pred[:, keep])):
                matrix = compute_signed_mi_matrix(mat)
                path = OUTPUT_DIR / f"{region}_{label}_mi_mat.txt"
                with path.open("w") as handle:
                    handle.write("\t" + "\t".join(map(str, positions)) + "\n")
                    for position, row in zip(positions, matrix):
                        handle.write(str(position) + "\t" + "\t".join(f"{x:.4f}" for x in row) + "\n")
            print(f"{number}/{len(regions)}: {region}", flush=True)


if __name__ == "__main__":
    main()
