#Run from SMFNet_paper_analysis_code github repo: python ./Processing_code/motif_annotation/name_motifs.py
from pathlib import Path
import pandas as pd
import h5py
OUTPUT = Path("output/motif_annotation")

MODISCO_FILE = OUTPUT / "modisco/modisco_results.h5"
REPORT_DIR = OUTPUT / "modisco/report"
OUTPUT_FILE = OUTPUT / "modisco/motif_names.tsv"


def main():
    rows = []
    with h5py.File(MODISCO_FILE) as handle:
        for group in ("pos_patterns", "neg_patterns"):
            for pattern in handle.get(group, {}):
                key = f"{group}.{pattern}"
                path = REPORT_DIR / "tomtom" / f"{key}.tomtom.tsv"
                row = dict(motif_name=key, motif="unassigned", jaspar_id="", q_value=float("nan"))
                if not path.exists():
                    raise FileNotFoundError(f"Tomtom results missing: {path}")
                matches = pd.read_csv(path, sep="\t")
                matches["q-value"] = pd.to_numeric(matches["q-value"], errors="coerce")
                matches = matches.dropna(subset=["q-value"]).sort_values("q-value", kind="stable")
                if len(matches):
                    best = matches.iloc[0]
                    target = str(best["Target ID"])
                    row.update(motif=target.rsplit("_MA", 1)[0], jaspar_id=target, q_value=float(best["q-value"]))
                rows.append(row)
    OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(OUTPUT_FILE, sep="\t", index=False)
    print(f"Wrote {len(rows)} proposed labels to {OUTPUT_FILE}; review these before annotating hits.")


if __name__ == "__main__":
    main()
