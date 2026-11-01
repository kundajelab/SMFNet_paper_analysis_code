#Run from the repository root: python figure_code/Fig4_E7_combinedMotifs_insilico/compute_phasing.py
from pathlib import Path
import json
import numpy as np
import pandas as pd

OUT = Path("output/Fig4_E7_combinedMotifs_insilico/median_global_residuals")
OUT.mkdir(parents=True, exist_ok=True)
q = pd.read_csv(
    Path(
        "output/Fig4_E7_combinedMotifs_insilico/selected_pairs_baseline/combined_quantification.txt"
    ),
    sep="\t",
    usecols=["spacing", "bg_region", "left", "center", "right", "perc_ac", "boundary_dispersion"],
)

#Compare each flanked construct with its center-only control on the same background.
a = q[(q.left == "noMotif") & (q.right == "noMotif") & (q.center != "noMotif")]
b = q[(q.left == q.right) & (q.left != "noMotif") & (q.center != "noMotif")]
keys = ["spacing", "flank_motif", "center_motif"]
d = (
    b[["spacing", "left", "center"]]
    .drop_duplicates()
    .rename(columns={"left": "flank_motif", "center": "center_motif"})
    .sort_values(keys)
    .reset_index(drop=True)
)
d["point"] = np.arange(len(d))
assert not a.duplicated(["spacing", "bg_region", "center"]).any()
p = b.merge(
    a[["spacing", "bg_region", "center", "perc_ac", "boundary_dispersion"]],
    on=["spacing", "bg_region", "center"],
    suffixes=("", "_A"),
    validate="many_to_one",
)
p = p.rename(columns={"left": "flank_motif", "center": "center_motif"}).merge(
    d[keys + ["point"]], on=keys, validate="many_to_one"
)
p["freq_delta"] = p.perc_ac - p.perc_ac_A
p["dispersion_delta"] = p.boundary_dispersion - p.boundary_dispersion_A
assert not p.duplicated(["bg_region", "point"]).any()

#Arrange matched-background changes as matrices before taking medians for each combination.
X = (
    p.pivot(index="bg_region", columns="point", values="freq_delta")
    .sort_index()
    .reindex(columns=d.point)
    .to_numpy()
)
Y = (
    p.pivot(index="bg_region", columns="point", values="dispersion_delta")
    .sort_index()
    .reindex(columns=d.point)
    .to_numpy()
)
assert X.shape == Y.shape and np.isfinite(X).all() and np.isfinite(Y).all()
d["freq_delta"] = np.median(X, axis=0)
d["dispersion_delta"] = np.median(Y, axis=0)


def fit(x, y):
    slope = np.dot(x - x.mean(), y - y.mean()) / np.dot(x - x.mean(), x - x.mean())
    return y.mean() - slope * x.mean(), slope


#Fit one global relationship between frequency and dispersion changes across all spacings.
intercept, slope = fit(d.freq_delta.to_numpy(), d.dispersion_delta.to_numpy())
d["residual"] = d.dispersion_delta - intercept - slope * d.freq_delta
d["n_backgrounds"] = X.shape[0]
d.to_csv(OUT / "all_spacings_scatter_data.txt", sep="\t", index=False)
focus = ["SP-KLF", "ESRRB", "OCT4::SOX2"]
selected = (d.center_motif == "OCT4::SOX2") & d.flank_motif.isin(focus)
idx = np.flatnonzero(selected)
nboot = 2000
seed = 20261021
rng = np.random.default_rng(seed)
res = np.empty((nboot, len(idx)))
coefs = np.empty((nboot, 2))

#Resample whole backgrounds and refit the global model to estimate residual uncertainty.
for i in range(nboot):
    draw = rng.integers(X.shape[0], size=X.shape[0])
    x = np.median(X[draw], axis=0)
    y = np.median(Y[draw], axis=0)
    a, b = fit(x, y)
    coefs[i] = a, b
    res[i] = y[idx] - a - b * x[idx]
    if (i + 1) % 200 == 0:
        print("Bootstrap", i + 1, flush=True)

#Export spacing-dependent residuals for the selected OCT4::SOX2 combinations.
plot = d.loc[selected].copy()
plot["ci_low"], plot["ci_high"] = np.quantile(res, [0.025, 0.975], axis=0)
plot.to_csv(OUT / "OCT4_SOX2_global_residual_spacing_CIs.txt", sep="\t", index=False)
pd.DataFrame({"term": ["intercept", "slope"], "estimate": [intercept, slope]}).to_csv(
    OUT / "global_all_spacings_linear_model.txt", sep="\t", index=False
)
info = dict(
    n_points=len(d),
    n_combinations=len(d[keys[1:]].drop_duplicates()),
    n_backgrounds=X.shape[0],
    n_bootstrap=nboot,
    seed=seed,
    intercept=intercept,
    slope=slope,
    r_squared=float(
        1 - np.sum(d.residual**2) / np.sum((d.dispersion_delta - d.dispersion_delta.mean()) ** 2)
    ),
    method="Per-background BAB minus center-only deltas, then separate medians for frequency and dispersion at each spacing/combination. OLS fit across all retained median points. Shared background bootstrap resampling; medians and global model refitted in each draw. Residual is median dispersion delta minus global fitted value at median frequency delta.",
)
(OUT / "statistics.json").write_text(json.dumps(info, indent=2))
print(info)
