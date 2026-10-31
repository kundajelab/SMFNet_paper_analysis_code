#Run from the repository root: python figure_code/Fig1_E3_Codependency_benchmarking/extract_coaccessibility_distances.py
#Run compute_coaccessibility.py first; reads its unsmoothed output matrices.
#Writes output/Fig1_E3_Codependency_benchmarking/coaccessibility_distances.tsv.
import csv
from pathlib import Path
import numpy as np
import pandas as pd


#Read the previous calculation and write the next pipeline output.
INPUT_DIR = Path("output/Fig1_E3_Codependency_benchmarking/coaccessibility_matrices")
OUTPUT = Path("output/Fig1_E3_Codependency_benchmarking/coaccessibility_distances.tsv")


def main():
    files = sorted(INPUT_DIR.glob("*_mi_mat.txt"))
    if not files:
        raise FileNotFoundError(f"No raw matrices in {INPUT_DIR}. Run compute_coaccessibility.py first from the repository root.")
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    #Write each region immediately rather than accumulating every pair in memory.
    with OUTPUT.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["region", "label", "pos1", "pos2", "distance", "mi"])
        for number, path in enumerate(files, 1):
            region, label = path.name.removesuffix("_mi_mat.txt").rsplit("_", 1)
            raw = pd.read_csv(path, sep="\t", header=None)
            row_pos = raw.iloc[1:, 0].to_numpy(dtype=int)
            col_pos = raw.iloc[0, 1:].to_numpy(dtype=int)
            matrix = raw.iloc[1:, 1:].to_numpy(dtype=float)
            #The strict upper triangle excludes self-pairs and duplicate pairs.
            for i in range(matrix.shape[0]):
                for j in range(i + 1, matrix.shape[1]):
                    value = matrix[i, j]
                    if not np.isnan(value):
                        pos1, pos2 = int(row_pos[i]), int(col_pos[j])
                        writer.writerow([region, label, pos1, pos2, abs(pos2 - pos1), value])
            print(f"{number}/{len(files)}: {path.name}", flush=True)


if __name__ == "__main__":
    main()
