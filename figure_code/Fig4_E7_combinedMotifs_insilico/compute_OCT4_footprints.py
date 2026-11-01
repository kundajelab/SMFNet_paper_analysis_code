#Run from the repository root: python figure_code/Fig4_E7_combinedMotifs_insilico/compute_OCT4_footprints.py
from pathlib import Path
import h5py
import numpy as np
import pandas as pd
from numba import njit
from scipy.stats import pearsonr

INPUT = Path("input_files/Fig4_E7_combinedMotifs_insilico")
SAMPLES = INPUT / "footprint_samples"
STATES = INPUT
CALLS = Path("output/Fig4_E7_combinedMotifs_insilico/oct4_footprints/raw_calls")
ROOT = Path("output/Fig4_E7_combinedMotifs_insilico/oct4_footprints")
CALLS.mkdir(parents=True, exist_ok=True)
CACHE = CALLS
WIDTHS = {"SP-KLF": 11, "ESRRB": 12}
CALL_LABELS = ["accessible", "bound", "long_protected", "mixed_center", "unresolved"]
THRESHOLDS = [70]
CENTER = np.arange(505, 520)


#Original raw protection caller and shared-site masking outside inserted flanks.
@njit(cache=False)
def call(mat, pos, cutoff=70.0):
    n = len(mat)
    state = np.full(n, 4, np.int8)
    length = np.full(n, np.nan)
    outer = np.full(n, np.nan)
    left = np.full(n, np.nan)
    right = np.full(n, np.nan)
    center = np.where((pos >= 504.5) & (pos < 519.5))[0]
    if not len(center):
        return state, length, outer, left, right
    for i in range(n):
        #Missing measurements are never skipped across to infer a boundary.
        vals = mat[i, center]
        if np.any(vals < 0):
            continue
        if np.all(vals == 0):
            state[i] = 0
            continue
        if not np.all(vals == 1):
            state[i] = 3
            continue
        a = center[0]
        b = center[-1]
        if np.any(mat[i, a : b + 1] != 1):
            continue
        while a > 0 and mat[i, a - 1] == 1:
            a -= 1
        while b < len(pos) - 1 and mat[i, b + 1] == 1:
            b += 1
        if a == 0 or b == len(pos) - 1:
            continue
        if mat[i, a - 1] != 0 or mat[i, b + 1] != 0:
            continue
        lo = (pos[a - 1] + pos[a]) / 2
        hi = (pos[b] + pos[b + 1]) / 2
        #Require midpoint expansion from the masked sites to cover 10-bp core.
        if lo > 507 or hi < 517:
            state[i] = 3
            continue
        length[i] = hi - lo + 1
        outer[i] = pos[b + 1] - pos[a - 1] - 1
        left[i] = pos[a - 1]
        right[i] = pos[b + 1]
        state[i] = 1 if length[i] < cutoff else 2
    return state, length, outer, left, right


#Compare center-only (A) and flanked (BAB) samples on matched sites for each background.
def call_samples(spacing):
    dest = CALLS / f"calls_{spacing}bp.h5"
    if dest.exists():
        return
    rows = []
    with h5py.File(SAMPLES / f"three_motif_spacing_{spacing}bp/samples.h5", "r") as h, h5py.File(
        dest.with_suffix(".writing.h5"), "w"
    ) as out:
        ids = {v: i for i, v in enumerate(h["region_ids"].asstr()[:])}
        for bi in range(1, 201):
            bg = f"region{bi}"
            an = f"{bg}_{spacing}_noMotif_OCT4::SOX2_noMotif"
            ag = h["regions"][f"{ids[an]:08d}"]
            ap = ag["predicted"]["positions"][:]
            am = ag["predicted"]["smf"][:]
            ad = ag["dna"][:]
            for partner, w in WIDTHS.items():
                bn = f"{bg}_{spacing}_{partner}_OCT4::SOX2_{partner}"
                g = h["regions"][f"{ids[bn]:08d}"]
                bp = g["predicted"]["positions"][:]
                bm = g["predicted"]["smf"][:]
                bd = g["dna"][:]
                assert np.isin(am, [-1, 0, 1]).all() and np.isin(bm, [-1, 0, 1]).all()
                assert int(ag.attrs["model_output_offset"]) == 512
                intervals = [
                    (512 - spacing - w // 2, 512 - spacing - w // 2 + w),
                    (512 + spacing - w // 2, 512 + spacing - w // 2 + w),
                ]
                allowed = np.ones(1024, bool)
                for lo, hi in intervals:
                    allowed[lo:hi] = False
                assert np.array_equal(ad[allowed], bd[allowed]), (bg, partner, spacing)
                #One-base padding excludes dinucleotides spanning insertion junctions.
                for lo, hi in intervals:
                    allowed[max(0, lo - 1) : min(1024, hi + 1)] = False
                shared = np.intersect1d(ap, bp)
                shared = shared[allowed[shared]]
                shared = np.array(
                    [
                        p
                        for p in shared
                        if np.array_equal(ad[max(0, p - 1) : p + 2], bd[max(0, p - 1) : p + 2])
                    ],
                    dtype=np.int64,
                )
                assert np.all(np.diff(shared) > 0)
                meta = dict(
                    bg_region=bg,
                    partner=partner,
                    spacing=spacing,
                    n_shared_sites=len(shared),
                    n_center_sites=int(((shared >= 504.5) & (shared < 519.5)).sum()),
                    n_original_A=len(ap),
                    n_original_BAB=len(bp),
                )
                for cond, p, m in [("A", ap, am), ("BAB", bp, bm)]:
                    idx = np.searchsorted(p, shared)
                    assert np.array_equal(p[idx], shared)
                    mat = m[:, idx]
                    st, ln, outer, left, right = call(mat, shared)
                    z = out.create_group(f"{partner}/{bg}/{cond}")
                    for k, v in [
                        ("positions", shared),
                        ("state", st),
                        ("length", ln),
                        ("outer_gap", outer),
                        ("left_accessible", left),
                        ("right_accessible", right),
                    ]:
                        z.create_dataset(k, data=v, compression="gzip")
                    z.attrs.update(meta)
                    r = dict(**meta, condition=cond, n_reads=len(m))
                    for j, name in enumerate(CALL_LABELS):
                        r[name] = int((st == j).sum())
                    r["bound_outer70"] = int(((st == 1) & (outer < 70)).sum())
                    #More conservative: accessible witnesses between center and insertion edges.
                    r["bound_inner_witnesses"] = int(
                        (
                            (st == 1)
                            & (left >= intervals[0][1] + 1)
                            & (right < intervals[1][0] - 1)
                        ).sum()
                    )
                    assert sum(r[k] for k in CALL_LABELS) == r["n_reads"]
                    rows.append(r)
            if bi % 50 == 0:
                print(spacing, bi, flush=True)
        out.attrs.update(complete=True, spacing=spacing)
    dest.with_suffix(".writing.h5").replace(dest)


#Require a covered center without nucleosome state 3 to classify a read as accessible.
def central_calls(states, pos):
    indices = np.searchsorted((pos[:-1] + pos[1:]) // 2, CENTER, side="right")
    measured = (pos >= 505) & (pos < 520)
    indices[pos[measured] - 505] = np.flatnonzero(measured)
    central = states[:, indices]
    covered = (pos[0] <= CENTER[0]) and (pos[-1] >= CENTER[-1])
    known = (central >= 0).all(axis=1) & covered
    nucleosome = (central == 3).any(axis=1)
    return known & ~nucleosome, known, nucleosome


#Count binding using the strict outer span below 70 bp, for all and accessible reads.
def count_binding(s):
    rows = []
    audit = []
    with h5py.File(CACHE / f"calls_{s}bp.h5") as h, h5py.File(
        STATES / f"three_motif_spacing_{s}bp_states.h5"
    ) as statefile:
        assert h.attrs["complete"]
        for partner in ["SP-KLF", "ESRRB"]:
            for bi in range(1, 201):
                bg = f"region{bi}"
                for cond in ["A", "BAB"]:
                    g = h[f"{partner}/{bg}/{cond}"]
                    pos = g["positions"][:]
                    st = g["state"][:]
                    ln = g["length"][:]
                    outer = g["outer_gap"][:]
                    left = g["left_accessible"][:]
                    right = g["right_accessible"][:]
                    n = len(st)
                    #Conservative span between the observed accessible boundary coordinates.
                    outer = right - left
                    name = (
                        f"{bg}_{s}_noMotif_OCT4::SOX2_noMotif"
                        if cond == "A"
                        else f"{bg}_{s}_{partner}_OCT4::SOX2_{partner}"
                    )
                    pp = statefile["smf_pos"][name][:]
                    ss = statefile["smf_states"][name][:]
                    ids = np.searchsorted(pp, pos)
                    assert np.array_equal(pp[ids], pos) and len(ss) == n
                    keep, known, nucleosome = central_calls(ss[:, ids], pos)
                    #Both length estimates refer to the same raw bounded central protected run.
                    finite = np.isfinite(ln)
                    assert np.array_equal(finite, (st == 1) | (st == 2))
                    assert (left[finite] < 504.5).all() and (right[finite] >= 519.5).all()
                    assert (outer[finite] >= ln[finite]).all()
                    assert int(g.attrs["n_center_sites"]) > 0
                    #Raw unresolved/ambiguous are never silently counted as unbound in resolved-only summaries.
                    raw_unresolved = (st == 3) | (st == 4)
                    audit.append(
                        dict(
                            partner=partner,
                            bg_region=bg,
                            spacing=s,
                            condition=cond,
                            n_reads=n,
                            n_open=int(keep.sum()),
                            n_missing_center=int((~known).sum()),
                            n_raw_unresolved=int(raw_unresolved.sum()),
                            n_raw_unresolved_open=int((raw_unresolved & keep).sum()),
                        )
                    )
                    for estimator, length in [("outer_bound", outer)]:
                        last = 0
                        for t in THRESHOLDS:
                            bound = finite & (length < t)
                            assert int(bound.sum()) >= last
                            last = int(bound.sum())
                            r = dict(
                                partner=partner,
                                bg_region=bg,
                                spacing=s,
                                condition=cond,
                                estimator=estimator,
                                threshold=t,
                                n_reads=n,
                                n_open=int(keep.sum()),
                                n_bound=int(bound.sum()),
                                n_bound_open=int((bound & keep).sum()),
                                n_resolved=int((~raw_unresolved).sum()),
                                n_resolved_open=int((~raw_unresolved & keep).sum()),
                            )
                            r["all_binding"] = r["n_bound"] / n
                            r["open_binding"] = (
                                r["n_bound_open"] / r["n_open"] if r["n_open"] else np.nan
                            )
                            r["resolved_binding"] = (
                                r["n_bound"] / r["n_resolved"] if r["n_resolved"] else np.nan
                            )
                            r["resolved_open_binding"] = (
                                r["n_bound_open"] / r["n_resolved_open"]
                                if r["n_resolved_open"]
                                else np.nan
                            )
                            rows.append(r)
    pd.DataFrame(rows).to_csv(ROOT / f"counts_{s}bp.tsv.gz", sep="\t", index=False)
    pd.DataFrame(audit).to_csv(ROOT / f"audit_{s}bp.tsv", sep="\t", index=False)
    print("Finished", s, flush=True)


#Bootstrap matched backgrounds to estimate uncertainty in the median BAB-minus-A change.
def summarize_binding():
    d = pd.concat(
        [pd.read_csv(ROOT / f"counts_{s}bp.tsv.gz", sep="\t") for s in range(20, 101, 5)],
        ignore_index=True,
    )
    a = pd.concat(
        [pd.read_csv(ROOT / f"audit_{s}bp.tsv", sep="\t") for s in range(20, 101, 5)],
        ignore_index=True,
    )
    d.to_csv(ROOT / "per_background_threshold_counts.tsv.gz", sep="\t", index=False)
    a.to_csv(ROOT / "read_accounting.tsv", sep="\t", index=False)
    rng = np.random.default_rng(20261023)
    draw = rng.integers(0, 200, (2000, 200))
    rows = []
    for key, g in d.groupby(["partner", "spacing", "estimator", "threshold"]):
        aa = g[g.condition == "A"].set_index("bg_region").sort_index()
        bb = g[g.condition == "BAB"].set_index("bg_region").sort_index()
        assert list(aa.index) == list(bb.index) and len(aa) == 200
        for metric in ["all_binding", "open_binding"]:
            x = 100 * aa[metric].to_numpy()
            y = 100 * bb[metric].to_numpy()
            ok = np.isfinite(x) & np.isfinite(y)
            for cond, values in [("BAB_minus_A", y - x)]:
                v = np.where(ok, values, np.nan)
                lo, hi = np.nanquantile(np.nanmedian(v[draw], axis=1), [0.025, 0.975])
                rows.append(
                    dict(
                        zip(["partner", "spacing", "estimator", "threshold"], key),
                        metric=metric,
                        condition=cond,
                        n_backgrounds=int(ok.sum()),
                        value=np.nanmedian(v),
                        low=lo,
                        high=hi,
                    )
                )
    pd.DataFrame(rows).to_csv(ROOT / "threshold_summary.tsv", sep="\t", index=False)


#Match BPNet count changes to SMFNet binding changes and compare medians across spacing.
def compare_bpnet():
    predictions = pd.read_csv(INPUT / "oct4_footprints/BPNet_matched_predictions.tsv.gz", sep="\t")
    v = predictions[
        (predictions.condition == "center_with_flanks") & predictions.partner.isin(WIDTHS)
    ].copy()
    a = predictions[predictions.condition == "center_only"][
        ["background", "spacing_bp", "predicted_counts"]
    ]
    v = v.merge(a, on=["background", "spacing_bp"], suffixes=("", "_A"), validate="many_to_one")
    v["bpnet_delta"] = v.predicted_counts - v.predicted_counts_A
    v = v.rename(columns={"background": "bg_region", "spacing_bp": "spacing"})
    d = pd.read_csv(ROOT / "per_background_threshold_counts.tsv.gz", sep="\t")
    keys = ["partner", "bg_region", "spacing", "estimator", "threshold"]
    paired = d[d.condition == "BAB"].merge(
        d[d.condition == "A"][keys + ["all_binding", "open_binding"]],
        on=keys,
        suffixes=("", "_A"),
        validate="one_to_one",
    )
    paired = paired.merge(
        v[["partner", "bg_region", "spacing", "bpnet_delta"]],
        on=["partner", "bg_region", "spacing"],
        validate="many_to_one",
    )
    rows = []
    for (partner, spacing), g in paired.groupby(["partner", "spacing"]):
        assert len(g) == 200
        rows.append(
            dict(
                partner=partner,
                spacing=spacing,
                bpnet_median=g.bpnet_delta.median(),
                smfnet_median=100 * (g.all_binding - g.all_binding_A).median(),
            )
        )
    b = pd.DataFrame(rows)
    b.to_csv(ROOT / "bpnet_delta_medians.tsv", sep="\t", index=False)
    r = [
        dict(partner=partner, pearson_r=float(pearsonr(g.bpnet_median, g.smfnet_median).statistic))
        for partner, g in b.groupby("partner")
    ]
    pd.DataFrame(r).to_csv(ROOT / "bpnet_delta_correlations.tsv", sep="\t", index=False)


if __name__ == "__main__":
    for spacing in range(20, 101, 5):
        call_samples(spacing)
        count_binding(spacing)
    summarize_binding()
    compare_bpnet()
