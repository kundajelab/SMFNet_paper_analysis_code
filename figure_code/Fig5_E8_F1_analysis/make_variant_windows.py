#Run from the repository root: python figure_code/Fig5_E8_F1_analysis/make_variant_windows.py
from pathlib import Path
Path("output/Fig5_E8_F1_analysis").mkdir(parents=True, exist_ok=True)
import pandas as pd

#Load your complete input file
df = pd.read_csv("output/Fig5_E8_F1_analysis/F1_atacPeaks_variant_hits_ac_metrics_filter.tsv", sep="\t")
df.columns = df.columns.str.strip()

#Parse region coordinates
region = df["region_id"].str.extract(
    r"^(chr(?:1|3|6)):(\d+)-(\d+)$"
)
df = df.loc[region[0].notna()].copy()
df["region_start"] = region.loc[df.index, 1].astype(int)
df["region_end"] = region.loc[df.index, 2].astype(int)

#Retain SNVs and indels
df = df[
    df["motif_variant"].str.lower().isin(["snv", "indel"])
].copy()

#Create 400-bp windows around motif centers
df["motif_center"] = (
    (df["motif_start"] + df["motif_end"]) // 2
)
df["window_start"] = df["motif_center"] - 200
df["window_end"] = df["window_start"] + 400

#Calculate window positions relative to region_id start
df["relative_window_start"] = (
    df["window_start"] - df["region_start"]
)
df["relative_window_end"] = (
    df["window_end"] - df["region_start"]
)

#Require both window boundaries to be within 512–1536 bp
df = df[
    (df["relative_window_start"] >= 512) &
    (df["relative_window_end"] <= 1536)
].copy()

#Select output columns
columns = [
    "region_id",
    "motif_start",
    "motif_end",
    "strand",
    "motif",
    "pattern_id",
    "motif_variant",
    "source_allele",
    "source_seq",
    "source_cwm_score",
    "source_cwm_strand",
    "opposite_allele",
    "opposite_seq",
    "opposite_cwm_score",
    "opposite_cwm_strand",
    "ref_seq",
    "ref_cwm_score",
    "ref_cwm_strand",
    "cwm_delta",
    "window_start",
    "window_end",
    "relative_window_start",
    "relative_window_end",
]

df[columns].to_csv(
    "output/Fig5_E8_F1_analysis/F1_atacPeaks_trainRegions_variant_motif_windows.tsv",
    sep="\t",
    index=False
)

print(f"Saved {len(df)} motif windows.")
