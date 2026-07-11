"""Sensitivity of the composite ranking to the choice of weights.

The composite score weights structural similarity, model confidence, and pocket
quality (default 0.5 / 0.2 / 0.3). A reviewer will reasonably ask whether the
ranking depends on those particular numbers. This script answers empirically:
it scores a pool of proteins (the benchmark set, reusing its cached structures
and TM matrix), then re-ranks them under *every* reasonable weight vector on the
simplex and measures how much the ranking changes relative to the default.

Stability is quantified by Kendall's tau (rank correlation to the default
ranking) and by top-10 overlap (Jaccard). High values across the simplex show
the ranking is driven by the data, not by the exact weights.

Run after the benchmark:  python benchmark/weight_sensitivity.py
"""

from __future__ import annotations

import json
import sys
from itertools import product
from pathlib import Path

import numpy as np
from scipy.stats import kendalltau

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from pipeline.parse import parse_structure          # noqa: E402
from pipeline.pocket_geom import detect_pockets      # noqa: E402
from pipeline.structio import load_atoms             # noqa: E402

RESULTS = ROOT / "benchmark" / "results"
FIGS = RESULTS / "figures"
DEFAULT_W = {"structural_similarity": 0.5, "model_confidence": 0.2, "pocket_quality": 0.3}
TOPK = 10


def load_components() -> tuple[list[dict], str]:
    """Return (component pool, description).

    Prefers a realistic candidate pool (benchmark/results/candidate_pool.json,
    built by build_candidate_pool.py) — actual rescue candidates whose component
    scores have a genuine spread. Falls back to the benchmark enzymes, whose
    scores are near-maximal and therefore a poor pool for ranking sensitivity.
    """
    pool_path = RESULTS / "candidate_pool.json"
    if pool_path.exists():
        pool = json.loads(pool_path.read_text(encoding="utf-8"))
        if len(pool) >= 10:
            return pool, "realistic candidate pool"

    cache = RESULTS / "components.json"
    if cache.exists():
        return json.loads(cache.read_text(encoding="utf-8")), "benchmark enzymes (fallback)"

    dataset = json.loads((RESULTS / "dataset.json").read_text(encoding="utf-8"))
    tm = np.load(RESULTS / "tm_matrix.npy")
    n = len(dataset)
    comps = []
    for i, d in enumerate(dataset):
        row = tm[i].copy()
        row[i] = -1
        struct = float(row.max()) if n > 1 else 0.0
        rec = parse_structure(Path(d["path"]))
        plddt = rec["mean_plddt"] / 100.0 if rec else 0.0
        coords, radii = load_atoms(d["path"])
        pk = detect_pockets(coords, radii, 1.2, 5.0, 5, 30)
        cav = pk.get("top_druggability")
        comps.append({
            "accession": d["accession"],
            "structural_similarity": round(max(struct, 0.0), 4),
            "model_confidence": round(plddt, 4),
            "pocket_quality": (round(cav, 4) if cav is not None else None),
        })
        print(f"  [{i+1}/{n}] {d['accession']}")
    cache.write_text(json.dumps(comps, indent=2), encoding="utf-8")
    return comps, "benchmark enzymes (computed)"


def composite(item: dict, w: dict) -> float:
    num = den = 0.0
    for k, wk in w.items():
        v = item.get(k)
        if v is None:
            continue
        num += wk * v
        den += wk
    return num / den if den else 0.0


def rank_order(items, w):
    scores = [composite(it, w) for it in items]
    # argsort descending -> ranking positions
    return np.argsort(np.argsort(-np.asarray(scores)))


def simplex_weights(step: float = 0.05, wmin: float = 0.05):
    """All (w1,w2,w3) on the simplex with each >= wmin, in `step` increments."""
    grid = np.round(np.arange(wmin, 1 - 2 * wmin + 1e-9, step), 3)
    out = []
    for w1, w2 in product(grid, grid):
        w3 = round(1 - w1 - w2, 3)
        if w3 >= wmin - 1e-9:
            out.append((float(w1), float(w2), float(w3)))
    return out


def reasonable_weights(step: float = 0.025):
    """Weight vectors in a *reasonable* neighbourhood of the default.

    Structure stays the dominant term, as it must for a structure-based method:
    structural_similarity in [0.40, 0.60], model_confidence in [0.10, 0.30],
    pocket_quality in [0.20, 0.40]. This is the region a user would plausibly pick.
    """
    g1 = np.round(np.arange(0.40, 0.60 + 1e-9, step), 3)
    g2 = np.round(np.arange(0.10, 0.30 + 1e-9, step), 3)
    out = []
    for w1 in g1:
        for w2 in g2:
            w3 = round(1 - w1 - w2, 3)
            if 0.20 - 1e-9 <= w3 <= 0.40 + 1e-9:
                out.append((float(w1), float(w2), float(w3)))
    return out


def stability(items, weights, default_rank, default_top):
    taus, jacc, records = [], [], []
    for (w1, w2, w3) in weights:
        w = {"structural_similarity": w1, "model_confidence": w2, "pocket_quality": w3}
        r = rank_order(items, w)
        tau = kendalltau(default_rank, r).statistic
        top = set(np.argsort([-composite(it, w) for it in items])[:TOPK])
        j = len(default_top & top) / len(default_top | top)
        taus.append(float(tau))
        jacc.append(float(j))
        records.append({"w": [w1, w2, w3], "kendall_tau": round(float(tau), 4),
                        "top10_jaccard": round(float(j), 4)})
    return np.asarray(taus), np.asarray(jacc), records


def _stats(taus, jacc):
    return {
        "n_weight_vectors": int(len(taus)),
        "kendall_tau": {
            "median": round(float(np.median(taus)), 4),
            "min": round(float(taus.min()), 4),
            "p05": round(float(np.percentile(taus, 5)), 4),
            "frac_above_0.9": round(float((taus > 0.9).mean()), 4),
        },
        "top10_jaccard": {
            "median": round(float(np.median(jacc)), 4),
            "min": round(float(jacc.min()), 4),
            "p05": round(float(np.percentile(jacc, 5)), 4),
        },
    }


def main():
    for s in (sys.stdout, sys.stderr):
        try:
            s.reconfigure(encoding="utf-8", errors="replace")
        except Exception:
            pass

    items, pool_desc = load_components()
    n = len(items)
    print(f"[weights] pool of {n} items ({pool_desc})")

    default_rank = rank_order(items, DEFAULT_W)
    default_top = set(np.argsort([-composite(it, DEFAULT_W) for it in items])[:TOPK])

    reasonable = reasonable_weights()
    simplex = simplex_weights()
    r_taus, r_jacc, r_rec = stability(items, reasonable, default_rank, default_top)
    s_taus, s_jacc, s_rec = stability(items, simplex, default_rank, default_top)

    summary = {
        "n_items": n,
        "pool": pool_desc,
        "default_weights": DEFAULT_W,
        "reasonable_neighbourhood": _stats(r_taus, r_jacc),
        "full_simplex_stress_test": _stats(s_taus, s_jacc),
    }
    (RESULTS / "weight_sensitivity.json").write_text(
        json.dumps({"summary": summary, "reasonable_grid": r_rec}, indent=2),
        encoding="utf-8")

    _make_figure(reasonable, r_taus, simplex, s_taus, summary, pool_desc)

    print("\n==== WEIGHT SENSITIVITY ====")
    print(f"pool: {n} items ({pool_desc})")
    for name, st in (("reasonable neighbourhood", summary["reasonable_neighbourhood"]),
                     ("full simplex (stress test)", summary["full_simplex_stress_test"])):
        kt = st["kendall_tau"]
        jc = st["top10_jaccard"]
        print(f"[{name}] {st['n_weight_vectors']} weightings | "
              f"Kendall tau median {kt['median']} (min {kt['min']}, "
              f"{kt['frac_above_0.9']*100:.0f}% > 0.9) | "
              f"top-{TOPK} Jaccard median {jc['median']}")
    print(f"written to {RESULTS/'weight_sensitivity.json'}")


def _make_figure(reasonable, r_taus, simplex, s_taus, summary, pool_desc):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.4))

    # (a) Kendall tau: reasonable neighbourhood vs full simplex.
    ax1.hist(s_taus, bins=20, color="#c6dbef", edgecolor="white",
             label="full simplex (stress test)")
    ax1.hist(r_taus, bins=20, color="#2b8cbe", edgecolor="white",
             label="reasonable neighbourhood")
    ax1.axvline(np.median(r_taus), color="red", ls="--",
                label=f"reasonable median = {np.median(r_taus):.3f}")
    ax1.set_xlabel("Kendall's tau vs. default-weight ranking")
    ax1.set_ylabel("number of weight vectors")
    ax1.set_title("Ranking stability to weight choice")
    ax1.legend(fontsize=8)

    # (b) tau over the simplex, with the reasonable band outlined.
    w1 = [w[0] for w in simplex]
    w3 = [w[2] for w in simplex]
    sc = ax2.scatter(w1, w3, c=s_taus, cmap="viridis", s=42,
                     vmin=float(np.min(s_taus)), vmax=1.0)
    # outline reasonable region
    ax2.add_patch(plt.Rectangle((0.40, 0.20), 0.20, 0.20, fill=False,
                                edgecolor="red", lw=1.8, ls="--",
                                label="reasonable region"))
    ax2.scatter([0.5], [0.3], marker="*", s=240, c="red", edgecolor="white",
                zorder=5, label="default (0.5, ·, 0.3)")
    ax2.set_xlabel("weight on structural similarity")
    ax2.set_ylabel("weight on pocket quality")
    ax2.set_title("Rank correlation over the weight simplex")
    ax2.legend(fontsize=8, loc="upper right")
    fig.colorbar(sc, ax=ax2, label="Kendall's tau")

    fig.suptitle(f"Composite-weight sensitivity ({pool_desc})", fontsize=10, y=1.02)
    fig.tight_layout()
    fig.savefig(FIGS / "fig_weight_sensitivity.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
