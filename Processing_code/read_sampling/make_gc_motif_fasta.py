#Run from SMFNet_paper_analysis_code github repo: python ./Processing_code/read_sampling/make_gc_motif_fasta.py

import os
import glob
import numpy as np
from tangermeme.io import read_meme

BACKGROUND_DIR = "output/read_sampling/gc_backgrounds"
MOTIF_CWM = "input_files/read_sampling/insertion_cwms.meme"
OUT_FASTA = "output/read_sampling/fasta/gc_single_motifs.fasta"

alphabet = np.array(["A", "C", "G", "T"])

def read_fasta(path):
    seqs, name, buf = {}, None, []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            elif line:
                buf.append(line)
    if name is not None:
        seqs[name] = "".join(buf)
    return seqs

meme = read_meme(MOTIF_CWM)

motifs = {}
for motif_id, cwm in meme.items():
    tf = "_".join(motif_id.split("_")[2:])
    cwm = np.asarray(cwm)
    motifs[tf] = "".join(alphabet[np.argmax(cwm, axis=0)])

gc_fastas = sorted(glob.glob(os.path.join(BACKGROUND_DIR, "*.fasta")))

print("GC FASTAs:", len(gc_fastas))
print("Motifs:", len(motifs))

os.makedirs(os.path.dirname(OUT_FASTA), exist_ok=True)

with open(OUT_FASTA, "w") as out:
    for fasta in gc_fastas:
        gc_label = os.path.basename(fasta).replace(".fasta", "").split("_")[-1]
        backgrounds = read_fasta(fasta)

        print(f"{gc_label}: {len(backgrounds)} backgrounds")

        for bg_name, bg in backgrounds.items():
            center = len(bg) // 2

            out.write(f">{bg_name}_{gc_label}_noMotif\n{bg}\n")

            for motif_name, motif in motifs.items():
                seq = list(bg)
                start = center - len(motif) // 2
                seq[start:start + len(motif)] = motif

                out.write(
                    f">{bg_name}_{gc_label}_{motif_name}\n"
                    f"{''.join(seq)}\n"
                )

print("Wrote:", OUT_FASTA)
