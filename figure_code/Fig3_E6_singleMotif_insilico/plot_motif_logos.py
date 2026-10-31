#Run from SMFNet_paper_analysis_code github repo: python ./figure_code/Fig3_E6_singleMotif_insilico/plot_motif_logos.py

from pathlib import Path

INPUT_DIR = Path("input_files/Fig3_E6_singleMotif_insilico")
INTERMEDIATE_DIR = Path("output")
RESULTS_DIR = Path("output")
(INTERMEDIATE_DIR / "Fig3_E6_singleMotif_insilico").mkdir(parents=True, exist_ok=True)
import os, re
import numpy as np
import matplotlib.pyplot as plt
from tangermeme.plot import plot_logo

MEME_IN = str(INPUT_DIR / "selected_cwms.meme")
MEME_OUT = str(INTERMEDIATE_DIR / "Fig3_E6_singleMotif_insilico/cwms_trimmed.meme")
PLOT_DIR = str(RESULTS_DIR / "Fig3_E6_singleMotif_insilico/motif_logos")
MAX_LEN = 16
os.makedirs(PLOT_DIR, exist_ok=True)
with open(MEME_IN) as f:
    lines = f.readlines()
header = []
motifs = []
i = 0
while i < len(lines):
    if not lines[i].startswith("MOTIF"):
        if not motifs:
            header.append(lines[i])
        i += 1
        continue
    motif_line = lines[i].rstrip()
    matrix_line = lines[i + 1].rstrip()
    name = motif_line.split()[1]
    w = int(re.search("w=\\s*(\\d+)", matrix_line).group(1))
    mat = np.array([list(map(float, lines[i + 2 + j].split())) for j in range(w)])
    if w > MAX_LEN:
        scores = [mat[s : s + MAX_LEN].max(axis=1).sum() for s in range(w - MAX_LEN + 1)]
        start = int(np.argmax(scores))
        mat = mat[start : start + MAX_LEN]
        print(f"{name}: {w} -> 16 bp, kept {start + 1}-{start + MAX_LEN}")
    else:
        print(f"{name}: {w} bp unchanged")
    motifs.append((motif_line, matrix_line, mat))
    i += 2 + w
with open(MEME_OUT, "w") as f:
    f.writelines(header)
    for motif_line, matrix_line, mat in motifs:
        matrix_line = re.sub("w=\\s*\\d+", f"w= {len(mat)}", matrix_line)
        f.write(motif_line + "\n")
        f.write(matrix_line + "\n")
        for row in mat:
            f.write(" ".join((f"{x:.6f}" for x in row)) + "\n")
        f.write("\n")
        safe_name = motif_line.split()[1].replace("::", "_").replace("-", "_")
        fig, ax = plt.subplots(figsize=(5, 2))
        plot_logo(mat.T, ax=ax)
        ax.set_title(motif_line.split()[1])
        ax.set_xlabel("Position")
        ax.set_ylabel("CWM score")
        plt.tight_layout()
        fig.savefig(os.path.join(PLOT_DIR, safe_name + ".pdf"), dpi=200)
        plt.close(fig)
print("Wrote:", MEME_OUT)
print("Plots:", PLOT_DIR)
