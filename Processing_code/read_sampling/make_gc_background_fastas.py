#Run from SMFNet_paper_analysis_code github repo: python ./Processing_code/read_sampling/make_gc_background_fastas.py

import os
import numpy as np
from tangermeme.utils import random_one_hot
from tangermeme.io import read_meme
from tqdm import tqdm

OUT_DIR = "output/read_sampling/gc_backgrounds//"
MOTIF_CWM = "input_files/read_sampling/insertion_cwms.meme"

GC_CONTENTS = range(30,71,1)
N_SEQS = 200
L = 2048
SEED = 33
THRESHOLD = 0.85
SHUFFLE_LEN = 40

os.makedirs(OUT_DIR,exist_ok=True)
rng = np.random.default_rng(SEED)

alphabet = np.array(["A","C","G","T"])
base_idx = {b:i for i,b in enumerate(alphabet)}

meme = read_meme(MOTIF_CWM)
motifs = {}

for name,cwm in meme.items():
    cwm = np.asarray(cwm)
    motifs[name] = {"cwm":cwm,"max_score":cwm.max(axis=0).sum()}

def revcomp_idx(idx):
    return 3-idx[::-1]

def score_window(idx,cwm):
    return cwm[idx,np.arange(len(idx))].sum()

def find_hits(seq):
    idx = np.array([base_idx[b] for b in seq])
    hits = []

    for name,m in motifs.items():
        cwm,max_score = m["cwm"],m["max_score"]
        w = cwm.shape[1]

        for start in range(len(idx)-w+1):
            x = idx[start:start+w]
            sf = score_window(x,cwm)/max_score
            sr = score_window(revcomp_idx(x),cwm)/max_score

            if sf >= THRESHOLD:
                hits.append({"name":name,"start":start,"length":w,"score":sf,"strand":"+"})
            if sr >= THRESHOLD:
                hits.append({"name":name,"start":start,"length":w,"score":sr,"strand":"-"})

    return hits

def remove_motif_like_hits(seq):
    seq = list(seq)
    n_shuffle = 0

    while True:
        hits = find_hits("".join(seq))
        if not hits:
            return "".join(seq),n_shuffle

        windows = []
        for hit in hits:
            center = hit["start"] + hit["length"]//2
            left = min(max(center-SHUFFLE_LEN//2,0),len(seq)-SHUFFLE_LEN)
            windows.append((left,left+SHUFFLE_LEN))

        windows = sorted(set(windows))
        merged = []

        for left,right in windows:
            if not merged or left > merged[-1][1]:
                merged.append([left,right])
            else:
                merged[-1][1] = max(merged[-1][1],right)

        for left,right in merged:
            window = seq[left:right]
            rng.shuffle(window)
            seq[left:right] = window
            n_shuffle += 1

for gc in GC_CONTENTS:
    print(f"Generating GC {gc}%")

    gc_frac = gc/100
    probs = np.array([[(1-gc_frac)/2,gc_frac/2,gc_frac/2,(1-gc_frac)/2]])

    X = random_one_hot(
        (N_SEQS,4,L),
        probs=probs,
        random_state=SEED+gc
    ).float()

    out_fasta = os.path.join(OUT_DIR,f"background_200_GC{gc}.fasta")
    total_shuffles = 0

    with open(out_fasta,"w") as f:
        for i,onehot in tqdm(enumerate(X),total=len(X),desc=f"GC {gc}%"):
            idx = onehot.argmax(dim=0).cpu().numpy()
            seq = "".join(alphabet[idx])

            seq,n_shuffle = remove_motif_like_hits(seq)
            total_shuffles += n_shuffle

            f.write(f">region{i+1}_GC{gc}\n{seq}\n")

    print(f"GC {gc}% | reshuffled {total_shuffles} regions | {out_fasta}")
