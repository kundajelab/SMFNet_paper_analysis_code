#Run from the repository root: python figure_code/Fig4_E7_combinedMotifs_insilico/compute_accessible_SMF_profiles.py
import os

os.environ["OPENBLAS_NUM_THREADS"] = "1"
from pathlib import Path
import json
import numpy as np, pandas as pd, h5py
from scipy.ndimage import gaussian_filter1d

ROOT = Path("output/Fig4_E7_combinedMotifs_insilico/oct4_accessible_profiles")
ROOT.mkdir(parents=True, exist_ok=True)
SAMPLES = Path("input_files/Fig4_E7_combinedMotifs_insilico/footprint_samples")
STATES = Path("input_files/Fig4_E7_combinedMotifs_insilico")
CENTER = np.arange(505, 520)


#Identify accessible reads from the central state calls, keeping unknown coverage separate.
def central_calls(states, pos):
    indices = np.searchsorted((pos[:-1] + pos[1:]) // 2, CENTER, side="right")
    measured = (pos >= 505) & (pos < 520)
    indices[pos[measured] - 505] = np.flatnonzero(measured)
    central = states[:, indices]
    covered = (pos[0] <= CENTER[0]) and (pos[-1] >= CENTER[-1])
    known = (central >= 0).all(axis=1) & covered
    nucleosome = (central == 3).any(axis=1)
    return known & ~nucleosome, known, nucleosome


#Pool measured protection and coverage across matched backgrounds for A and BAB constructs.
def extract(s):
    totals = {
        (p, c, k): [np.zeros(1024, np.int64), np.zeros(1024, np.int64), np.zeros(1024, np.int64)]
        for p in ["KLF", "ESRRB"]
        for c in ["A", "BAB"]
        for k in ["accessible"]
    }
    audit = []
    with h5py.File(SAMPLES / f"three_motif_spacing_{s}bp/samples.h5") as raw, h5py.File(
        STATES / f"three_motif_spacing_{s}bp_states.h5"
    ) as state:
        ids = {name: i for i, name in enumerate(raw["region_ids"].asstr()[:])}
        for bi, partner in [(bi, partner) for bi in range(1, 201) for partner in ["KLF", "ESRRB"]]:
            bg = f"region{bi}"
            motif = "SP-KLF" if partner == "KLF" else partner
            ag = raw["regions"][f"{ids[f'{bg}_{s}_noMotif_OCT4::SOX2_noMotif']:08d}"]
            bg_group = raw["regions"][f"{ids[f'{bg}_{s}_{motif}_OCT4::SOX2_{motif}']:08d}"]
            ap = ag["predicted"]["positions"][:]
            ad = ag["dna"][:]
            bp = bg_group["predicted"]["positions"][:]
            bd = bg_group["dna"][:]
            shared = np.intersect1d(ap, bp)

            #Same coordinate and same informative CG/GC context, including junctions.
            #Context represented as CG/GC status, not unrelated neighbouring bases.
            def context(d, p):
                return (
                    bool(p > 0 and d[p - 1] == 2 and d[p] == 1),
                    bool(p + 1 < len(d) and d[p] == 1 and d[p + 1] == 2),
                )

            shared = np.array(
                [p for p in shared if context(ad, p) == context(bd, p)], dtype=np.int64
            )
            assert len(shared) > 1 and np.all(np.diff(shared) > 0)
            assert all(
                any(context(ad, p)) for p in shared if 0 < p < len(ad) - 1
            ), "Expected cytosine CG/GC coordinates"
            for cond in ["A", "BAB"]:
                name = (
                    f"{bg}_{s}_noMotif_OCT4::SOX2_noMotif"
                    if cond == "A"
                    else f"{bg}_{s}_{motif}_OCT4::SOX2_{motif}"
                )
                z = raw["regions"][f"{ids[name]:08d}"]["predicted"]
                pos = z["positions"][:]
                mat = z["smf"][:]
                ss = state["smf_states"][name][:]
                sp = state["smf_pos"][name][:]
                assert np.array_equal(pos, sp) and mat.shape == ss.shape
                assert ((ss < 0) | ((ss >= 2) == mat)).all()
                ix = np.searchsorted(pos, shared)
                assert np.array_equal(pos[ix], shared)
                mat = mat[:, ix]
                ss = ss[:, ix]
                accessible, known, nuc = central_calls(ss, shared)
                nonaccessible = known & nuc
                assert not np.any(accessible & nonaccessible)
                assert int(accessible.sum() + nonaccessible.sum() + (~known).sum()) == len(mat)
                audit.append(
                    dict(
                        partner=partner,
                        bg_region=bg,
                        condition=cond,
                        spacing=s,
                        n_total=len(mat),
                        n_accessible=int(accessible.sum()),
                        n_nonaccessible=int(nonaccessible.sum()),
                        n_unknown=int((~known).sum()),
                        n_shared_sites=len(shared),
                    )
                )
                #Use only center-accessible molecules in the protection profiles.
                for cohort, keep in [("accessible", accessible)]:
                    selected = mat[keep]
                    occ, cov, bgs = totals[(partner, cond, cohort)]
                    observed = (selected >= 0).sum(axis=0)
                    occ[shared] += (selected == 1).sum(axis=0)
                    cov[shared] += observed
                    bgs[shared] += observed > 0
    rows = []
    for (partner, cond, cohort), (occ, cov, bgs) in totals.items():
        assert (occ <= cov).all()
        signal = np.divide(occ, cov, out=np.full(1024, np.nan), where=cov > 0)
        rows.append(
            pd.DataFrame(
                dict(
                    partner=partner,
                    condition=cond,
                    cohort=cohort,
                    spacing=s,
                    position=np.arange(1024),
                    pos=np.arange(1024) - 512,
                    n_occupied=occ,
                    n_observed=cov,
                    n_backgrounds_observed=bgs,
                    signal=signal,
                )
            )
        )
    pd.concat(rows).to_csv(ROOT / f"direct_{s}bp.tsv.gz", sep="\t", index=False)
    pd.DataFrame(audit).to_csv(ROOT / f"audit_{s}bp.tsv", sep="\t", index=False)
    print("Done", s, flush=True)


#Combine the 20-bp and 60-bp results and retain read-accounting summaries.
def finish():
    d = pd.concat(
        [pd.read_csv(ROOT / f"direct_{s}bp.tsv.gz", sep="\t") for s in [20, 60]], ignore_index=True
    )
    a = pd.concat(
        [pd.read_csv(ROOT / f"audit_{s}bp.tsv", sep="\t") for s in [20, 60]], ignore_index=True
    )
    assert len(d) == 2 * 2 * 2 * 1024 and len(a) == 2 * 2 * 2 * 200
    d.to_csv(ROOT / "direct_site_profiles.tsv.gz", sep="\t", index=False)
    a.to_csv(ROOT / "molecule_accounting.tsv", sep="\t", index=False)
    rows = []
    for key, g in d.groupby(["partner", "condition", "cohort", "spacing"]):
        g = g.sort_values("position")
        for sigma in [4.0]:
            #Smooth protection counts and coverage separately, then divide to account for missing sites.
            numerator = gaussian_filter1d(
                g.n_occupied.to_numpy(float), sigma=sigma, truncate=4, mode="constant", cval=0
            )
            denominator = gaussian_filter1d(
                g.n_observed.to_numpy(float), sigma=sigma, truncate=4, mode="constant", cval=0
            )
            y = np.divide(numerator, denominator, out=np.full(1024, np.nan), where=denominator > 0)
            assert np.all((y[np.isfinite(y)] >= 0) & (y[np.isfinite(y)] <= 1 + 1e-12))
            z = g[["position", "pos"]].copy()
            for k, v in zip(["partner", "condition", "cohort", "spacing"], key):
                z[k] = v
            z["sigma_bp"] = sigma
            z["signal"] = y
            z["smoothed_support"] = denominator
            rows.append(z)
    pd.concat(rows).to_csv(ROOT / "smoothed_profiles.tsv.gz", sep="\t", index=False)
    f = a.groupby(["partner", "condition", "spacing"], as_index=False)[
        ["n_total", "n_accessible", "n_nonaccessible", "n_unknown"]
    ].sum()
    for kind in ["accessible", "nonaccessible", "unknown"]:
        f[kind + "_percent"] = 100 * f["n_" + kind] / f.n_total
    f.to_csv(ROOT / "cohort_fractions.tsv", sep="\t", index=False)


if __name__ == "__main__":
    for spacing in [20, 60]:
        extract(spacing)
    finish()
