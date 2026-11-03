#Run from SMFNet_paper_analysis_code github repo: python ./Processing_code/read_sampling/make_homotypic_fastas.py

import os
import numpy as np
from tangermeme.io import read_meme

BACKGROUND_FASTA = "output/read_sampling/background_noMotif.fasta"
MOTIF_CWM = "input_files/read_sampling/insertion_cwms.meme"
OUT_DIR = "output/read_sampling/fasta/"

SPACINGS = [20,30,40,50,60,70,80]
COPY_NUMBERS = range(1,8)
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

print("Backgrounds:", len(backgrounds))
print("Motifs:", len(motifs))

for spacing in SPACINGS:
    out_fasta = os.path.join(OUT_DIR, f"homotypic_spacing_{spacing}bp.fasta")

    with open(out_fasta, "w") as out:
        for bg_name, bg in backgrounds.items():
            center = len(bg)//2

            out.write(f">{bg_name}_{spacing}_0copy_noMotif\n{bg}\n")

            for motif_name, motif in motifs.items():
                for n_copies in COPY_NUMBERS:
                    offsets = (np.arange(n_copies) - (n_copies-1)/2) * spacing
                    centers = np.rint(center + offsets).astype(int)

                    seq = list(bg)
                    for pos in centers:
                        start = pos - len(motif)//2
                        seq[start:start+len(motif)] = motif

                    out.write(
                        f">{bg_name}_{spacing}_{n_copies}copy_{motif_name}\n"
                        f"{''.join(seq)}\n"
                    )

    print("Wrote:", out_fasta)
