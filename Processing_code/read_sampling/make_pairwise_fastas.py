#Run from SMFNet_paper_analysis_code github repo: python ./Processing_code/read_sampling/make_pairwise_fastas.py

import os
from itertools import product

import numpy as np
from tangermeme.io import read_meme

BACKGROUND_FASTA = "output/read_sampling/background_noMotif.fasta"
MOTIF_CWM = "input_files/read_sampling/insertion_cwms.meme"
OUT_DIR = "output/read_sampling/fasta//"

#Center-to-center distance between the central and side motifs.
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

    #All ordered pairs, including identical motifs and noMotif.
    #The first motif stays central; the second moves left or right.
    motif_names = [*motifs, "noMotif"]
    combinations = list(product(motif_names, repeat=2))
    sides = [("left", -1), ("right", 1)]

    os.makedirs(OUT_DIR, exist_ok=True)

    records_per_background = 2 * len(combinations)

    print("Backgrounds:", len(backgrounds))
    print("Motifs:", len(motifs))
    print("Combinations per background per spacing:", records_per_background)
    print("Sequences per FASTA:", len(backgrounds) * records_per_background)
    print(
        "Total sequences:",
        len(backgrounds) * records_per_background * len(SPACINGS),
    )

    for spacing in SPACINGS:
        out_fasta = os.path.join(
            OUT_DIR,
            f"pairwise_spacing_{spacing}bp.fasta",
        )

        with open(out_fasta, "w") as out:
            for bg_name, bg in backgrounds.items():
                center = len(bg) // 2

                for central_motif, side_motif in combinations:
                    for side, direction in sides:
                        seq = list(bg)

                        placements = [
                            (center, central_motif),
                            (center + direction * spacing, side_motif),
                        ]
                        occupied = []

                        for pos, motif_name in placements:
                            if motif_name == "noMotif":
                                continue

                            motif = motifs[motif_name]
                            start = pos - len(motif) // 2
                            end = start + len(motif)

                            if start < 0 or end > len(bg):
                                raise ValueError(
                                    f"Motif outside background: {bg_name}, "
                                    f"{spacing} bp, {side}, {motif_name}"
                                )

                            for previous_start, previous_end in occupied:
                                if start < previous_end and end > previous_start:
                                    raise ValueError(
                                        f"Overlapping motifs: {bg_name}, "
                                        f"{spacing} bp, {side}, "
                                        f"{central_motif}, {side_motif}"
                                    )

                            occupied.append((start, end))
                            seq[start:end] = motif

                        out.write(
                            f">{bg_name}_{spacing}"
                            f"_center-{central_motif}"
                            f"_{side}-{side_motif}\n"
                            f"{''.join(seq)}\n"
                        )

        print("Wrote:", out_fasta)


if __name__ == "__main__":
    main()
