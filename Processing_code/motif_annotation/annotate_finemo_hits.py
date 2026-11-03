#Run from SMFNet_paper_analysis_code github repo: python ./Processing_code/motif_annotation/annotate_finemo_hits.py
import re
from pathlib import Path
import numpy as np
import pandas as pd
INPUT = Path("input_files/motif_annotation")
OUTPUT = Path("output/motif_annotation")

HITS_FILE = OUTPUT / "finemo/hits.tsv"
REGION_IDS_FILE = OUTPUT / "shap/region_ids.npz"
REGIONS_FILE = INPUT / "mESC_annotation_regions_2048bp.bed"
NAMES_FILE = OUTPUT / "modisco/motif_names.tsv"
NAMED_MEME_FILE = None
OUTPUT_FILE = OUTPUT / "finemo/annotated_hits.tsv"


def read_region_ids(path):
    with np.load(path) as data:
        key = "arr_0" if "arr_0" in data else "region_ids"
        return data[key].astype(str)


def match_regions(region_ids, bed_path):
    bed = pd.read_csv(bed_path, sep="\t", header=None, usecols=[0, 1, 2], names=["chr", "start", "end"])
    if not ((bed.end - bed.start) == 2048).all():
        raise ValueError("Expected 2048 bp annotation regions")
    bed["order"] = np.arange(len(bed))
    groups = {chrom: group.sort_values("start", kind="stable") for chrom, group in bed.groupby("chr", sort=False)}
    matched = []
    #Only potentially overlapping intervals are examined; ties retain BED order.
    for region in region_ids:
        chrom, coordinates = region.rsplit(":", 1)
        start, end = map(int, coordinates.split("-"))
        if chrom not in groups:
            matched.append(np.nan)
            continue
        group = groups[chrom]
        starts = group.start.to_numpy()
        lo, hi = np.searchsorted(starts, [start - 2048, end], side="left")
        candidates = group.iloc[lo:hi]
        overlap = np.minimum(end, candidates.end.to_numpy()) - np.maximum(start, candidates.start.to_numpy())
        if not len(overlap) or overlap.max() <= 0:
            matched.append(np.nan)
            continue
        winners = candidates.iloc[np.flatnonzero(overlap == overlap.max())]
        best = winners.loc[winners.order.idxmin()]
        matched.append(f"{chrom}:{best.start}-{best.end}")
    return matched


def labels_from_meme(path):
    labels = {}
    for line in Path(path).read_text().splitlines():
        if not line.startswith("MOTIF "):
            continue
        name = line.split()[1]
        match = re.fullmatch(r"(?:(pos_patterns|neg_patterns)[._])?(pattern_\d+)_(.+)", name)
        if match:
            group, pattern, label = match.groups()
            labels[f"{group or 'pos_patterns'}.{pattern}"] = label
    if not labels:
        raise ValueError("No named motifs found in MEME file")
    return labels


def main():
    if NAMED_MEME_FILE:
        labels = labels_from_meme(NAMED_MEME_FILE)
    else:
        names = pd.read_csv(NAMES_FILE, sep="\t")
        if names.motif_name.duplicated().any():
            raise ValueError("Motif IDs must be unique")
        labels = dict(zip(names.motif_name, names.motif))
    region_ids = read_region_ids(REGION_IDS_FILE)
    hits = pd.read_csv(HITS_FILE, sep="\t")
    if not hits.peak_id.between(0, len(region_ids) - 1).all() or not (hits.peak_id == hits.peak_id.astype(int)).all():
        raise ValueError("Hit peak IDs do not index the attribution arrays")
    hits["pattern"] = hits.motif_name.str.extract(r"(pattern_\d+)", expand=False)
    hits["motif"] = hits.motif_name.map(labels)
    mapping = pd.DataFrame(dict(peak_id=np.arange(len(region_ids)), region=match_regions(region_ids, REGIONS_FILE)))
    hits = hits.merge(mapping, on="peak_id", how="left", validate="many_to_one", sort=False)
    OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)
    hits.to_csv(OUTPUT_FILE, sep="\t", index=False)
    print(f"Wrote {len(hits)} hits; {hits.motif.isna().sum()} unnamed; {hits.region.isna().sum()} without a matching region")


if __name__ == "__main__":
    main()
