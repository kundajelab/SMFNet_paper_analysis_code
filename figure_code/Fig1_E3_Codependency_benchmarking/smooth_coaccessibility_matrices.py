#Run from the repository root: python figure_code/Fig1_E3_Codependency_benchmarking/smooth_coaccessibility_matrices.py
#Run compute_coaccessibility.py first; reads its unsmoothed output matrices.
#Outputs go to output/Fig1_E3_Codependency_benchmarking/coaccessibility_smoothed_matrices.
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter


#Read the previous calculation and write the next pipeline output.
INPUT_DIR = Path("output/Fig1_E3_Codependency_benchmarking/coaccessibility_matrices")
OUTPUT_DIR = Path("output/Fig1_E3_Codependency_benchmarking/coaccessibility_smoothed_matrices")


def sparse_to_dense(matrix, positions, full_len=1024, radius=20):
    #Average overlapping 41-by-41 bp neighbourhoods; retain unobserved pixels as NaN.
    total = np.zeros((full_len, full_len), dtype=np.float32)
    count = np.zeros_like(total)
    for i, pi in enumerate(positions):
        pi = int(pi)
        if pi < 0 or pi >= full_len:
            continue
        for j, pj in enumerate(positions):
            pj = int(pj)
            if pj < 0 or pj >= full_len or np.isnan(matrix[i, j]):
                continue
            x = slice(max(0, pi - radius), min(full_len, pi + radius + 1))
            y = slice(max(0, pj - radius), min(full_len, pj + radius + 1))
            total[x, y] += matrix[i, j]
            count[x, y] += 1
    dense = np.divide(total, count, where=count > 0)
    dense[count == 0] = np.nan
    return dense


def main():
    files = sorted(INPUT_DIR.glob("*_mi_mat.txt"))
    if not files:
        raise FileNotFoundError(f"No raw matrices in {INPUT_DIR}. Run compute_coaccessibility.py first from the repository root.")
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    #Process one matrix at a time instead of keeping the entire dataset in memory.
    for number, path in enumerate(files, 1):
        raw = pd.read_csv(path, sep="\t", header=None)
        positions = raw.iloc[0, 1:].to_numpy(dtype=int)
        matrix = raw.iloc[1:, 1:].to_numpy(dtype=float)
        dense = sparse_to_dense(matrix, positions)
        #Preserve the original Gaussian filter, including NaN propagation at gaps.
        dense = gaussian_filter(dense, sigma=0.8, mode="nearest", truncate=3.0)
        output = OUTPUT_DIR / path.name.replace("_mi_mat.txt", "_mi_1024.txt")
        with output.open("w") as handle:
            handle.write("\t" + "\t".join(map(str, range(1024))) + "\n")
            for i, row in enumerate(dense):
                values = "\t".join("nan" if np.isnan(x) else f"{x:.6f}" for x in row)
                handle.write(f"{i}\t{values}\n")
        print(f"{number}/{len(files)}: {path.name}", flush=True)


if __name__ == "__main__":
    main()
