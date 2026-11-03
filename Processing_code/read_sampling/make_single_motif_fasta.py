#Run from SMFNet_paper_analysis_code github repo: python ./Processing_code/read_sampling/make_single_motif_fasta.py

import os
import numpy as np
from tangermeme.io import read_meme

BACKGROUND_FASTA = "output/read_sampling/background_noMotif.fasta"
MOTIF_CWM = "input_files/read_sampling/insertion_cwms.meme"
OUT_DIR = "output/read_sampling/fasta/"

SPACINGS = [50]
alphabet = np.array(["A","C","G","T"])
os.makedirs(OUT_DIR, exist_ok=True)

def read_fasta(path):
    seqs, name, buf = {}, None, []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if name is not None: seqs[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            elif line:
                buf.append(line)
    if name is not None: seqs[name] = "".join(buf)
    return seqs

backgrounds = read_fasta(BACKGROUND_FASTA)
meme = read_meme(MOTIF_CWM)

motifs = {}
for motif_id, cwm in meme.items():
    tf = "_".join(motif_id.split("_")[2:])
    cwm = np.asarray(cwm)
    motifs[tf] = "".join(alphabet[np.argmax(cwm, axis=0)])

motif_ids = list(motifs)

combinations = [("noMotif", "noMotif", "noMotif")]
combinations += [("noMotif", m, "noMotif") for m in motif_ids]

print("Backgrounds:", len(backgrounds))
print("Motifs:", len(motif_ids))
print("Combinations per spacing:", len(combinations))

for spacing in SPACINGS:
    out_fasta = os.path.join(OUT_DIR, "single_motifs.fasta")

    with open(out_fasta, "w") as out:
        for bg_name, bg in backgrounds.items():
            center = len(bg)//2
            centers = [center-spacing, center, center+spacing]

            for comb in combinations:
                seq = list(bg)

                for pos, motif_name in zip(centers, comb):
                    if motif_name == "noMotif": continue
                    motif = motifs[motif_name]
                    start = pos - len(motif)//2
                    seq[start:start+len(motif)] = motif

                out.write(
                    f">{bg_name}_{spacing}_{comb[0]}_{comb[1]}_{comb[2]}\n"
                    f"{''.join(seq)}\n"
                )

    print("Wrote:", out_fasta)
