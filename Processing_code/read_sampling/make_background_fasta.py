#Run from SMFNet_paper_analysis_code github repo: python ./Processing_code/read_sampling/make_background_fasta.py

from pathlib import Path

import numpy as np
from tangermeme.utils import random_one_hot
from tangermeme.io import read_meme
from tqdm import tqdm

OUT_FASTA = "output/read_sampling/background_noMotif.fasta"
MOTIF_CWM = "input_files/read_sampling/insertion_cwms.meme"

N_SEQS = 200
L = 2048
SEED = 33
THRESHOLD = 0.85
SHUFFLE_LEN = 40

rng = np.random.default_rng(SEED)

alphabet = np.array(["A","C","G","T"])
base_idx = {b:i for i,b in enumerate(alphabet)}

meme = read_meme(MOTIF_CWM)

motifs = {}
for name, cwm in meme.items():
    cwm = np.asarray(cwm)
    motifs[name] = {"cwm": cwm, "max_score": cwm.max(axis=0).sum()}

def revcomp_idx(idx):
    return 3 - idx[::-1]

def score_window(idx, cwm):
    return cwm[idx, np.arange(len(idx))].sum()

def find_hits(seq):
    idx = np.array([base_idx[b] for b in seq])
    hits = []

    for name, m in motifs.items():
        cwm = m["cwm"]
        max_score = m["max_score"]
        w = cwm.shape[1]

        for start in range(len(idx) - w + 1):
            x = idx[start:start+w]
            sf = score_window(x, cwm) / max_score
            sr = score_window(revcomp_idx(x), cwm) / max_score

            if sf >= THRESHOLD:
                hits.append({"name":name, "start":start, "length":w, "score":sf, "strand":"+"})
            if sr >= THRESHOLD:
                hits.append({"name":name, "start":start, "length":w, "score":sr, "strand":"-"})

    return hits

def remove_motif_like_hits(seq):
    seq = list(seq)
    n_shuffle = 0
    round_i = 0

    while True:
        hits = find_hits("".join(seq))
        if not hits:
            return "".join(seq), n_shuffle

        round_i += 1
        windows = []

        for hit in hits:
            center = hit["start"] + hit["length"] // 2
            left = min(max(center - SHUFFLE_LEN // 2, 0), len(seq) - SHUFFLE_LEN)
            right = left + SHUFFLE_LEN
            windows.append((left, right))

        #merge overlapping shuffle windows
        windows = sorted(set(windows))
        merged = []

        for left, right in windows:
            if not merged or left > merged[-1][1]:
                merged.append([left, right])
            else:
                merged[-1][1] = max(merged[-1][1], right)

        print(f"round {round_i}: {len(hits)} motif hits -> {len(merged)} shuffle regions")

        for left, right in merged:
            window = seq[left:right]
            rng.shuffle(window)
            seq[left:right] = window
            n_shuffle += 1

probs = np.array([[0.25,0.25,0.25,0.25]])
X = random_one_hot((N_SEQS,4,L), probs=probs, random_state=SEED).float()

total_shuffles = 0

Path(OUT_FASTA).parent.mkdir(parents=True, exist_ok=True)
with open(OUT_FASTA,"w") as f:
    for i, onehot in tqdm(enumerate(X), total=len(X)):
        idx = onehot.argmax(dim=0).cpu().numpy()
        seq = "".join(alphabet[idx])

        seq, n_shuffle = remove_motif_like_hits(seq)
        total_shuffles += n_shuffle

        f.write(f">region{i+1}\n{seq}\n")

print("Wrote:", OUT_FASTA)
print("Total motif-like regions reshuffled:", total_shuffles)
