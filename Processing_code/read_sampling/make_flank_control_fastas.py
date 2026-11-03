#Run from SMFNet_paper_analysis_code github repo: python ./Processing_code/read_sampling/make_flank_control_fastas.py

import os
import numpy as np
from tangermeme.io import read_meme

BACKGROUND_FASTA = "output/read_sampling/background_noMotif.fasta"
MOTIF_CWM = "input_files/read_sampling/insertion_cwms.meme"
OUT_DIR = "output/read_sampling/fasta"

SPACINGS = range(20, 101, 5)
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

                if name in seqs:
                    raise ValueError(f"Duplicate background ID: {name}")

            elif line:
                buf.append(line)

    if name is not None:
        seqs[name] = "".join(buf)

    return seqs


def main():
    backgrounds = read_fasta(BACKGROUND_FASTA)
    meme = read_meme(MOTIF_CWM)

    #Same motif names and consensus extraction as your original script.
    motifs = {}

    for motif_id, cwm in meme.items():
        tf = "_".join(motif_id.split("_")[2:])
        cwm = np.asarray(cwm)

        if cwm.ndim != 2 or cwm.shape[0] != 4 or cwm.shape[1] == 0:
            raise ValueError(f"Expected a nonempty (4, L) CWM: {motif_id}")

        if not tf or tf == "noMotif" or tf in motifs:
            raise ValueError(f"Invalid or duplicate parsed motif name: {tf!r}")

        motifs[tf] = "".join(alphabet[np.argmax(cwm, axis=0)])

    if not backgrounds or not motifs:
        raise ValueError("Backgrounds and motifs must both be nonempty.")

    #Only the missing B–B controls: identical flanks, no central insertion.
    combinations = [(b, "noMotif", b) for b in motifs]

    os.makedirs(OUT_DIR, exist_ok=True)

    print("Backgrounds:", len(backgrounds))
    print("Motifs:", len(motifs))
    print("Combinations per spacing:", len(combinations))
    print("Sequences per FASTA:", len(backgrounds) * len(combinations))

    for spacing in SPACINGS:
        #Separate filenames preserve your existing triplet FASTAs.
        out_fasta = os.path.join(
            OUT_DIR,
            f"flank_only_spacing_{spacing}bp.fasta",
        )

        with open(out_fasta, "w") as out:
            for bg_name, bg in backgrounds.items():
                center = len(bg) // 2
                centers = [center - spacing, center, center + spacing]

                for comb in combinations:
                    seq = list(bg)

                    for pos, motif_name in zip(centers, comb):
                        if motif_name == "noMotif":
                            continue

                        motif = motifs[motif_name]
                        start = pos - len(motif) // 2
                        end = start + len(motif)

                        if start < 0 or end > len(bg):
                            raise ValueError(
                                f"Motif outside background: {bg_name}, "
                                f"{spacing} bp, {motif_name}"
                            )

                        seq[start:end] = motif

                    out.write(
                        f">{bg_name}_{spacing}_{comb[0]}_{comb[1]}_{comb[2]}\n"
                        f"{''.join(seq)}\n"
                    )

        print("Wrote:", out_fasta)


if __name__ == "__main__":
    main()
