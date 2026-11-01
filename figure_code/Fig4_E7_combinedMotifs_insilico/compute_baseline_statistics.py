#Run from the repository root: python figure_code/Fig4_E7_combinedMotifs_insilico/compute_baseline_statistics.py
from pathlib import Path
import numpy as np
import pandas as pd

INPUT = Path("input_files/Fig4_E7_combinedMotifs_insilico")
METRICS_DIR = Path("output/Fig4_E7_combinedMotifs_insilico/selected_pairs_baseline")
OUT = Path("output/Fig4_E7_combinedMotifs_insilico/baselines")
OUT.mkdir(parents=True, exist_ok=True)
METRICS = {"freq": "perc_ac", "width": "med_width", "dispersion": "boundary_dispersion"}
N_BOOT = 2000


#Original Markov frequency and log-additive width/dispersion baselines.
def action(p0, pi):
    u = np.divide(np.maximum(pi - p0, 0), 1 - p0, out=np.zeros_like(p0), where=p0 < 1)
    v = np.divide(np.maximum(p0 - pi, 0), p0, out=np.zeros_like(p0), where=p0 > 0)
    return u, v


#Match each flanked construct to center-only, flank-only and no-motif controls.
def calculate(df, pairs):
    keys = ["spacing", "bg_region"]
    cols = list(METRICS.values())
    bab = df[(df.left == df.right) & (df.left != "noMotif") & (df.center != "noMotif")].copy()
    a = df[(df.left == "noMotif") & (df.right == "noMotif") & (df.center != "noMotif")]
    f = df[(df.left == df.right) & (df.left != "noMotif") & (df.center == "noMotif")]
    bg = df[(df.left == "noMotif") & (df.center == "noMotif") & (df.right == "noMotif")]
    for control, join, suffix in [
        (a, keys + ["center"], "a"),
        (f, keys + ["left"], "f"),
        (bg, keys, "bg"),
    ]:
        table = control[join + cols].rename(columns={c: c + "_" + suffix for c in cols})
        bab = bab.merge(table, on=join, how="left", validate="many_to_one")
    bab["combo"] = bab.left + "_" + bab.center + "_" + bab.right
    bab["pair_training_count"] = [
        pairs[tuple(sorted((b, a)))] for b, a in zip(bab.left, bab.center)
    ]
    outputs = []
    for measure, col in METRICS.items():
        out = bab[
            ["spacing", "bg_region", "left", "center", "right", "combo", "pair_training_count"]
        ].copy()
        combined, p0, pa, pf = [
            bab[c].to_numpy(float) for c in [col, col + "_bg", col + "_a", col + "_f"]
        ]
        inputs = np.column_stack([p0, pa, pf])
        finite = np.isfinite(inputs).all(axis=1)
        valid = (
            finite & ((inputs >= 0) & (inputs <= 1)).all(axis=1)
            if measure == "freq"
            else finite & (inputs > 0).all(axis=1)
        )
        expected = np.full(len(bab), np.nan)
        af, fa = expected.copy(), expected.copy()
        #For frequency, average the two orders of the Markov-style motif effects.
        if measure == "freq":
            ua, va = action(p0[valid], pa[valid])
            uf, vf = action(p0[valid], pf[valid])
            af[valid] = pa[valid] * (1 - vf) + (1 - pa[valid]) * uf
            fa[valid] = pf[valid] * (1 - va) + (1 - pf[valid]) * ua
            expected[valid] = (af[valid] + fa[valid]) / 2
            obs_valid = np.isfinite(combined) & (combined >= 0) & (combined <= 1)
        else:
            #For width and dispersion, add control effects on the log scale relative to background.
            expected[valid] = np.exp(np.log(pa[valid]) + np.log(pf[valid]) - np.log(p0[valid]))
            obs_valid = np.isfinite(combined) & (combined >= 0)
        eligible = valid & np.isfinite(expected) & obs_valid
        out["measure"] = measure
        out["baseline"] = "symmetric_markov" if measure == "freq" else "log_additive_two_factor"
        out["background"], out["center_only"], out["flanks_only"] = p0, pa, pf
        out["combined"], out["expected"] = combined, expected
        out["delta"] = np.where(eligible, combined - expected, np.nan)
        out["expected_eligible"], out["eligible"] = valid & np.isfinite(expected), eligible
        out["ineligible_reason"] = np.select(
            [~finite, finite & ~valid, ~np.isfinite(expected), ~obs_valid],
            [
                "missing_or_nonfinite_control",
                "invalid_control_range" if measure == "freq" else "nonpositive_control",
                "nonfinite_expectation",
                "missing_or_invalid_observed",
            ],
            default="",
        )
        out["p_A_then_F"], out["p_F_then_A"] = af, fa
        out["order_lower"], out["order_upper"] = np.minimum(af, fa), np.maximum(af, fa)
        out["order_gap"] = np.abs(af - fa)
        outputs.append(out)
    return pd.concat(outputs, ignore_index=True)


def pearson_rows(x, y):
    xc = x - x.mean(axis=-1, keepdims=True)
    yc = y - y.mean(axis=-1, keepdims=True)
    return (xc * yc).sum(axis=-1) / np.sqrt((xc * xc).sum(axis=-1) * (yc * yc).sum(axis=-1))


def bootstrap_medians(values, counts):
    """Exact medians of bootstrap samples, using their integer multiplicities."""
    n, n_columns = values.shape
    assert np.all(counts.sum(axis=1) == n)
    order = np.argsort(values, axis=0)
    sorted_values = np.take_along_axis(values, order, axis=0)
    result = np.empty((len(counts), n_columns))
    for start in range(0, len(counts), 64):
        stop = min(start + 64, len(counts))
        cumulative = np.cumsum(counts[start:stop, order], axis=1)
        lower = (cumulative >= (n + 1) // 2).argmax(axis=1)
        upper = (cumulative >= (n + 2) // 2).argmax(axis=1)
        result[start:stop] = (
            np.take_along_axis(sorted_values, lower, axis=0)
            + np.take_along_axis(sorted_values, upper, axis=0)
        ) / 2
    #Independent check against explicitly expanded bootstrap samples.
    for i in range(min(3, len(counts))):
        sample = np.repeat(np.arange(n), counts[i])
        assert np.allclose(result[i], np.median(values[sample], axis=0), rtol=0, atol=1e-14)
    return result


def main():
    counts = pd.read_csv(INPUT / "motif_pair_counts_training_within_100bp.txt", sep="\t")
    #Keep motif pairs observed more than 1,000 times in the training regions.
    counts = counts[counts.n_pairs > 1000].copy()
    counts[["motif_1", "motif_2"]] = counts[["motif_1", "motif_2"]].replace({"KLF": "SP-KLF"})
    pairs = {tuple(sorted((r.motif_1, r.motif_2))): r.n_pairs for r in counts.itertuples()}
    quant = pd.read_csv(METRICS_DIR / "combined_quantification.txt", sep="\t")
    summaries, statistics = [], []
    for spacing in range(20, 101, 5):
        raw = calculate(quant[quant.spacing == spacing], pairs)
        valid = raw[raw.eligible].copy()
        keys = [
            "measure",
            "spacing",
            "combo",
            "left",
            "center",
            "right",
            "pair_training_count",
            "baseline",
        ]
        summary = (
            valid.groupby(keys, sort=True)
            .agg(
                combined_median=("combined", "median"),
                expected_median=("expected", "median"),
                n_eligible=("bg_region", "size"),
            )
            .reset_index()
        )
        summary["delta"] = summary.combined_median - summary.expected_median
        summaries.append(summary)
        #Resample whole matched backgrounds, keeping all 184 combinations together.
        for measure, t in valid.groupby("measure", sort=True):
            xt = (
                t.pivot(index="bg_region", columns="combo", values="expected")
                .sort_index()
                .sort_index(axis=1)
            )
            yt = t.pivot(index="bg_region", columns="combo", values="combined").reindex_like(xt)
            x, y = xt.to_numpy(), yt.to_numpy()
            assert x.shape == y.shape == (200, 184)
            assert np.isfinite(x).all() and np.isfinite(y).all()
            seed = 20261021 + spacing
            counts = np.random.default_rng(seed).multinomial(
                200, np.full(200, 1 / 200), size=N_BOOT
            )
            bx, by = bootstrap_medians(x, counts), bootstrap_medians(y, counts)
            mx, my = np.median(x, axis=0), np.median(y, axis=0)
            #Report signed mean error in percentage points for frequency and bp for other measures.
            scale = 100 if measure == "freq" else 1
            r = float(pearson_rows(mx, my))
            rlo, rhi = np.quantile(pearson_rows(bx, by), [0.025, 0.975])
            elo, ehi = np.quantile(scale * (by - bx).mean(axis=1), [0.025, 0.975])
            statistics.append(
                dict(
                    measure=measure,
                    spacing=spacing,
                    pearson_r=r,
                    pearson_ci_low=rlo,
                    pearson_ci_high=rhi,
                    mean_error=scale * np.mean(my - mx),
                    mean_error_ci_low=elo,
                    mean_error_ci_high=ehi,
                    n_bootstrap=N_BOOT,
                    seed=seed,
                )
            )
        print("Finished baselines:", spacing, "bp", flush=True)
    pd.concat(summaries, ignore_index=True).sort_values(["measure", "spacing", "combo"]).to_csv(
        OUT / "combination_medians.tsv", sep="\t", index=False
    )
    pd.DataFrame(statistics).to_csv(OUT / "baseline_statistics.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
