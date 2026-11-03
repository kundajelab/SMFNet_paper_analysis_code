#Run from SMFNet_paper_analysis_code github repo: python ./Processing_code/motif_annotation/trim_motif_CWMs.py
import re
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import logomaker
OUTPUT = Path("output/motif_annotation")

CWM_FILE = OUTPUT / "modisco/chrombpnet_motifs_cwm.meme"
OUTPUT_FILE = OUTPUT / "modisco/chrombpnet_motifs_cwm_trim.meme"
LOGO_DIR = OUTPUT / "modisco/logos_trimmed"


def trim_matrix(matrix, window=5, threshold=0.3):
    importance = np.abs(matrix).sum(axis=1)
    smooth = np.convolve(importance, np.ones(window) / window, mode="same")
    indices = np.flatnonzero(smooth > threshold * smooth.max())
    if not len(indices):
        return matrix, 0, len(matrix)
    start, end = indices[0], indices[-1] + 1
    return matrix[start:end], start, len(matrix) - end


def parse_meme(path):
    lines = Path(path).read_text().splitlines()
    motifs, i = [], 0
    while i < len(lines):
        if not lines[i].startswith("MOTIF "):
            i += 1
            continue
        name = lines[i].split()[1]
        i += 1
        while i < len(lines) and not lines[i].lstrip().startswith("letter-probability"):
            i += 1
        if i == len(lines):
            raise ValueError(f"No matrix for {name}")
        width = int(re.search(r"w=\s*(\d+)", lines[i]).group(1))
        matrix = np.asarray([list(map(float, row.split())) for row in lines[i + 1:i + 1 + width]])
        if matrix.shape != (width, 4):
            raise ValueError(f"Invalid matrix for {name}")
        motifs.append((name, matrix))
        i += width + 1
    return motifs


def write_meme(motifs, path):
    with Path(path).open("w") as handle:
        handle.write("MEME version 4\n\n")
        for name, matrix in motifs:
            handle.write(f"MOTIF {name}\nletter-probability matrix: alength= 4 w= {len(matrix)} nsites= 1\n")
            for row in matrix:
                handle.write(" ".join(f"{value:.6f}" for value in row) + "\n")
            handle.write("\n")


def main():
    OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)
    LOGO_DIR.mkdir(parents=True, exist_ok=True)
    trimmed, report = [], []
    for name, matrix in parse_meme(CWM_FILE):
        values, left, right = trim_matrix(matrix)
        trimmed.append((name, values))
        report.append(dict(motif_name=name, original_width=len(matrix), trimmed_width=len(values), left_trim=left, right_trim=right))
        fig, ax = plt.subplots(figsize=(8, 2))
        logomaker.Logo(logomaker.transform_matrix(pd.DataFrame(values, columns=list("ACGT"))), ax=ax)
        ax.set_title(name)
        ax.set_xticks([])
        ax.set_yticks([])
        fig.tight_layout()
        fig.savefig(LOGO_DIR / (re.sub(r"[^A-Za-z0-9_.-]", "_", name) + ".png"))
        plt.close(fig)
    write_meme(trimmed, OUTPUT_FILE)
    pd.DataFrame(report).to_csv(OUTPUT_FILE.with_suffix(".trimming.tsv"), sep="\t", index=False)


if __name__ == "__main__":
    main()
