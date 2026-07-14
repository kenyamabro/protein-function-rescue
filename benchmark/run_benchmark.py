"""Validation benchmark for the Protein Function Rescue pipeline.

Question: when we *hide* a characterised enzyme's annotation and predict its
function from its nearest structural neighbour, how often is the prediction
correct — and does structure recover functions that sequence search misses?

Design
------
1. Build a labelled gold-standard set: several EC families (full 4-level EC
   numbers) across all enzyme classes, with multiple diverse members each, all
   with AlphaFold models.
2. Leave-one-out: for every protein, find its best hit among all the others by
   (a) TM-align structural similarity and (b) local sequence identity.
3. Score EC agreement at levels 1-4 for both, giving structure-vs-sequence
   top-1 accuracy.
4. "Twilight zone" analysis: among correct structural recoveries, how many occur
   below 30 % sequence identity (where homology is undetectable from sequence).
5. Calibration: precision as a function of the TM-score acceptance threshold.

Outputs: benchmark/results/metrics.json, per_query.csv, and figures/*.png.
Everything is cached so the benchmark can resume.

Run:  python benchmark/run_benchmark.py            (uses defaults below)
"""

from __future__ import annotations

import json
import sys
from itertools import combinations
from pathlib import Path

import numpy as np
import requests
from tqdm import tqdm

# Reuse the pipeline's tested components.
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from pipeline import download as download_mod          # noqa: E402
from pipeline.structio import load_ca                  # noqa: E402
from tmtools import tm_align                            # noqa: E402
from Bio.Align import PairwiseAligner, substitution_matrices  # noqa: E402

# --------------------------------------------------------------------------- #
# Configuration
# --------------------------------------------------------------------------- #
# Full EC numbers spanning enzyme classes 1-6, each a distinct, well-populated
# structural family, chosen so leave-one-out recovery is possible. Families that
# return fewer than MIN_MEMBERS usable entries are dropped automatically.
EC_FAMILIES = [
    # Oxidoreductases (class 1)
    "1.1.1.1",    # alcohol dehydrogenase
    "1.1.1.27",   # L-lactate dehydrogenase
    "1.1.1.37",   # malate dehydrogenase
    "1.1.1.49",   # glucose-6-phosphate dehydrogenase
    "1.11.1.6",   # catalase
    "1.11.1.7",   # peroxidase
    "1.15.1.1",   # superoxide dismutase
    "1.8.1.7",    # glutathione reductase
    # Transferases (class 2)
    "2.6.1.1",    # aspartate aminotransferase
    "2.1.1.45",   # thymidylate synthase
    "2.7.1.1",    # hexokinase
    "2.7.1.40",   # pyruvate kinase
    "2.4.2.1",    # purine nucleoside phosphorylase
    "2.5.1.18",   # glutathione transferase
    # Hydrolases (class 3)
    "3.2.1.1",    # alpha-amylase
    "3.5.2.6",    # beta-lactamase
    "3.4.21.4",   # trypsin
    "3.1.1.7",    # acetylcholinesterase
    "3.1.3.1",    # alkaline phosphatase
    "3.2.1.17",   # lysozyme
    "3.1.1.3",    # triacylglycerol lipase
    # Lyases (class 4)
    "4.2.1.1",    # carbonic anhydrase
    "4.1.2.13",   # fructose-bisphosphate aldolase
    "4.2.1.11",   # enolase
    "4.3.1.3",    # histidine ammonia-lyase
    # Isomerases (class 5)
    "5.3.1.1",    # triosephosphate isomerase
    "5.3.1.9",    # glucose-6-phosphate isomerase
    "5.4.2.11",   # phosphoglycerate mutase
    "5.1.1.1",    # alanine racemase
    # Ligases (class 6)
    "6.3.1.2",    # glutamine synthetase
    "6.2.1.1",    # acetate--CoA ligase
    "6.3.4.2",    # CTP synthase
]
MEMBERS_PER_FAMILY = 8
MIN_MEMBERS = 3                      # drop families with fewer usable members
LEN_MIN, LEN_MAX = 80, 520          # raised cap -> larger enzymes, ~a few hundred proteins
TWILIGHT_IDENTITY = 30.0            # % identity below which sequence search fails

BDIR = ROOT / "benchmark"
RESULTS = BDIR / "results"
STRUCT_DIR = BDIR / "structures"
FIGS = RESULTS / "figures"
for d in (RESULTS, STRUCT_DIR, FIGS):
    d.mkdir(parents=True, exist_ok=True)

UNIPROT = "https://rest.uniprot.org/uniprotkb/search"


# --------------------------------------------------------------------------- #
# 1. Build the labelled dataset
# --------------------------------------------------------------------------- #
def build_dataset() -> list[dict]:
    manifest_path = RESULTS / "dataset.json"
    if manifest_path.exists():
        data = json.loads(manifest_path.read_text(encoding="utf-8"))
        print(f"[bench] using cached dataset ({len(data)} proteins)")
        return data

    seen: set[str] = set()
    dataset: list[dict] = []
    for ec in EC_FAMILIES:
        params = {
            "query": (
                f"(ec:{ec}) AND (reviewed:true) AND (existence:1) "
                f"AND (fragment:false) AND (length:[{LEN_MIN} TO {LEN_MAX}])"
            ),
            "fields": "accession,ec,protein_name,organism_name,sequence,length",
            "format": "tsv",
            "size": MEMBERS_PER_FAMILY * 3,
            "sort": "annotation_score desc",
        }
        r = requests.get(UNIPROT, params=params, timeout=90)
        r.raise_for_status()
        kept = 0
        for row in r.text.strip().splitlines()[1:]:
            cols = row.split("\t")
            if len(cols) < 6:
                continue
            acc, ecs, name, org, seq, length = cols[:6]
            if acc in seen or not seq:
                continue
            try:
                path = download_mod.download_single(acc, STRUCT_DIR)
            except Exception:
                continue  # no AlphaFold model
            seen.add(acc)
            dataset.append({
                "accession": acc, "ec": ec, "name": name, "organism": org,
                "sequence": seq, "length": int(length or 0), "path": str(path),
            })
            kept += 1
            if kept >= MEMBERS_PER_FAMILY:
                break
        print(f"[bench] EC {ec}: {kept} members")

    # Drop families with too few members for a valid leave-one-out test.
    counts: dict[str, int] = {}
    for d in dataset:
        counts[d["ec"]] = counts.get(d["ec"], 0) + 1
    usable = {ec for ec, c in counts.items() if c >= MIN_MEMBERS}
    dropped = sorted(set(counts) - usable)
    if dropped:
        print(f"[bench] dropping families with <{MIN_MEMBERS} members: {dropped}")
    dataset = [d for d in dataset if d["ec"] in usable]

    manifest_path.write_text(json.dumps(dataset, indent=2), encoding="utf-8")
    print(f"[bench] dataset built: {len(dataset)} proteins across "
          f"{len(usable)} families")
    return dataset


# --------------------------------------------------------------------------- #
# 2. Pairwise structural and sequence similarity
# --------------------------------------------------------------------------- #
def structural_matrix(dataset: list[dict]) -> np.ndarray:
    cache = RESULTS / "tm_matrix.npy"
    n = len(dataset)
    if cache.exists():
        m = np.load(cache)
        if m.shape == (n, n):
            print("[bench] using cached TM matrix")
            return m

    # Resumable checkpointing: this is an O(n^2) TM-align sweep that can take
    # hours for a few-hundred-protein set, so we save progress periodically and
    # resume where we left off if interrupted.
    ckpt = RESULTS / "tm_matrix.partial.npy"
    prog = RESULTS / "tm_progress.txt"
    pairs = list(combinations(range(n), 2))
    if ckpt.exists() and np.load(ckpt).shape == (n, n) and prog.exists():
        tm = np.load(ckpt)
        start = int(prog.read_text().strip() or 0)
        print(f"[bench] resuming TM matrix from pair {start}/{len(pairs)}")
    else:
        tm = np.zeros((n, n))
        start = 0

    coords = [load_ca(d["path"]) for d in tqdm(dataset, desc="load-ca")]
    CKPT_EVERY = 500
    for k in tqdm(range(start, len(pairs)), desc="tm-align", initial=start,
                  total=len(pairs)):
        i, j = pairs[k]
        (ci, si), (cj, sj) = coords[i], coords[j]
        if len(si) >= 10 and len(sj) >= 10:
            try:
                res = tm_align(ci, cj, si, sj)
                score = max(res.tm_norm_chain1, res.tm_norm_chain2)
            except Exception:
                score = 0.0
            tm[i, j] = tm[j, i] = score
        if (k + 1) % CKPT_EVERY == 0:
            np.save(ckpt, tm)
            prog.write_text(str(k + 1))

    np.save(cache, tm)
    ckpt.unlink(missing_ok=True)
    prog.unlink(missing_ok=True)
    return tm


def _directed_identity(aligner: PairwiseAligner, a: str, b: str) -> float:
    """One directed local-alignment identity calculation."""
    try:
        aln = aligner.align(a, b)[0]
    except Exception:
        return 0.0
    idn = 0
    for (s1, e1), (s2, e2) in zip(*aln.aligned):
        for x, y in zip(a[s1:e1], b[s2:e2]):
            if x == y:
                idn += 1
    return 100.0 * idn / max(1, min(len(a), len(b)))


def _identity(aligner: PairwiseAligner, a: str, b: str) -> float:
    """Symmetric local % identity, normalised by the shorter sequence.

    Local alignment can select a different tied optimum when the arguments are
    reversed. Averaging both directions makes the pairwise matrix independent
    of dataset ordering while retaining the original accessible baseline.
    """
    return 0.5 * (
        _directed_identity(aligner, a, b)
        + _directed_identity(aligner, b, a)
    )


def sequence_matrix(dataset: list[dict]) -> np.ndarray:
    cache = RESULTS / "seq_matrix.npy"
    n = len(dataset)
    if cache.exists():
        m = np.load(cache)
        if m.shape == (n, n):
            print("[bench] using cached sequence-identity matrix")
            return m

    aligner = PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1
    aligner.mode = "local"

    seqs = [d["sequence"] for d in dataset]
    ident = np.zeros((n, n))
    for i, j in tqdm(list(combinations(range(n), 2)), desc="seq-align"):
        v = _identity(aligner, seqs[i], seqs[j])
        ident[i, j] = ident[j, i] = v
    np.save(cache, ident)
    return ident


# --------------------------------------------------------------------------- #
# 3. Evaluation
# --------------------------------------------------------------------------- #
def ec_agreement(ec_a: str, ec_b: str) -> int:
    """Number of leading EC levels that match (0-4)."""
    pa, pb = ec_a.split("."), ec_b.split(".")
    lvl = 0
    for x, y in zip(pa, pb):
        if x == y and x != "-":
            lvl += 1
        else:
            break
    return lvl


def evaluate(dataset, tm, seq):
    n = len(dataset)
    rows = []
    for i in range(n):
        tm_row = tm[i].copy(); tm_row[i] = -1
        sq_row = seq[i].copy(); sq_row[i] = -1
        j_tm = int(np.argmax(tm_row))
        j_sq = int(np.argmax(sq_row))
        rows.append({
            "accession": dataset[i]["accession"],
            "true_ec": dataset[i]["ec"],
            "name": dataset[i]["name"],
            "length": dataset[i]["length"],
            "struct_hit": dataset[j_tm]["accession"],
            "struct_hit_ec": dataset[j_tm]["ec"],
            "tm_score": round(float(tm_row[j_tm]), 4),
            "struct_hit_identity": round(float(seq[i, j_tm]), 2),
            "struct_agree": ec_agreement(dataset[i]["ec"], dataset[j_tm]["ec"]),
            "seq_hit": dataset[j_sq]["accession"],
            "seq_hit_ec": dataset[j_sq]["ec"],
            "seq_identity": round(float(sq_row[j_sq]), 2),
            "seq_agree": ec_agreement(dataset[i]["ec"], dataset[j_sq]["ec"]),
        })
    return rows


def summarise(rows):
    n = len(rows)
    acc = {"structure": {}, "sequence": {}}
    for lvl in range(1, 5):
        acc["structure"][lvl] = round(sum(r["struct_agree"] >= lvl for r in rows) / n, 3)
        acc["sequence"][lvl] = round(sum(r["seq_agree"] >= lvl for r in rows) / n, 3)

    # Paired exact-reaction comparison. McNemar's exact test uses only
    # discordant queries and is appropriate because both methods are evaluated
    # on the same 227 proteins.
    from scipy.stats import binomtest
    struct4 = [r["struct_agree"] >= 4 for r in rows]
    seq4 = [r["seq_agree"] >= 4 for r in rows]
    both_correct = sum(s and q for s, q in zip(struct4, seq4))
    struct_only = sum(s and not q for s, q in zip(struct4, seq4))
    seq_only = sum(not s and q for s, q in zip(struct4, seq4))
    both_wrong = sum(not s and not q for s, q in zip(struct4, seq4))
    discordant = struct_only + seq_only
    mcnemar_p = (binomtest(struct_only, discordant, 0.5, alternative="two-sided").pvalue
                 if discordant else 1.0)

    # Twilight zone: correct (>=3-level) structural recoveries at low identity.
    correct3 = [r for r in rows if r["struct_agree"] >= 3]
    twilight = [r for r in correct3 if r["struct_hit_identity"] < TWILIGHT_IDENTITY]
    twilight_frac = round(len(twilight) / len(correct3), 3) if correct3 else 0.0

    # Precision vs TM threshold (target: >=3-level EC agreement).
    calib = []
    for thr in np.round(np.arange(0.3, 0.95, 0.05), 2):
        acc_thr = [r for r in rows if r["tm_score"] >= thr]
        if not acc_thr:
            continue
        prec = sum(r["struct_agree"] >= 3 for r in acc_thr) / len(acc_thr)
        calib.append({
            "threshold": float(thr),
            "precision": round(prec, 3),
            "coverage": round(len(acc_thr) / n, 3),
        })

    return {
        "n_proteins": n,
        "n_families": len({r["true_ec"] for r in rows}),
        "topN_accuracy_by_ec_level": acc,
        "paired_exact_ec_comparison": {
            "test": "exact McNemar (two-sided)",
            "both_correct": both_correct,
            "structure_only_correct": struct_only,
            "sequence_only_correct": seq_only,
            "both_wrong": both_wrong,
            "accuracy_difference_percentage_points": round(
                100 * (sum(struct4) - sum(seq4)) / n, 3),
            "p_value": mcnemar_p,
        },
        "twilight_zone": {
            "identity_cutoff_pct": TWILIGHT_IDENTITY,
            "n_correct_level3": len(correct3),
            "n_correct_below_cutoff": len(twilight),
            "fraction_below_cutoff": twilight_frac,
        },
        "precision_vs_tm_threshold": calib,
    }


# --------------------------------------------------------------------------- #
# 4. Figures
# --------------------------------------------------------------------------- #
def make_figures(rows, summary):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # (a) EC-level accuracy: structure vs sequence.
    levels = [1, 2, 3, 4]
    s = [summary["topN_accuracy_by_ec_level"]["structure"][l] for l in levels]
    q = [summary["topN_accuracy_by_ec_level"]["sequence"][l] for l in levels]
    x = np.arange(len(levels)); w = 0.38
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.bar(x - w / 2, s, w, label="Structure (TM-align)", color="#2b8cbe")
    ax.bar(x + w / 2, q, w, label="Sequence (local identity)", color="#fdae61")
    ax.set_xticks(x); ax.set_xticklabels([f"Level {l}" for l in levels])
    ax.set_ylabel("Top-1 EC agreement (fraction)"); ax.set_ylim(0, 1)
    ax.set_title("Leave-one-out function recovery"); ax.legend()
    for xi, v in zip(x - w / 2, s): ax.text(xi, v + 0.02, f"{v:.2f}", ha="center", fontsize=8)
    for xi, v in zip(x + w / 2, q): ax.text(xi, v + 0.02, f"{v:.2f}", ha="center", fontsize=8)
    fig.tight_layout(); fig.savefig(FIGS / "fig_ec_accuracy.png", dpi=150); plt.close(fig)

    # (b) TM vs sequence identity, coloured by correctness (the twilight-zone plot).
    fig, ax = plt.subplots(figsize=(6, 4.2))
    for r in rows:
        ok = r["struct_agree"] >= 3
        ax.scatter(r["struct_hit_identity"], r["tm_score"],
                   c=("#1a9850" if ok else "#d73027"), s=26,
                   edgecolor="k", linewidth=0.3, alpha=0.85)
    ax.axvline(TWILIGHT_IDENTITY, ls="--", c="gray", lw=1)
    ax.axhline(0.5, ls=":", c="gray", lw=1)
    ax.set_xlabel("Sequence identity to structural hit (%)")
    ax.set_ylabel("TM-score to structural hit")
    ax.set_title("Correct structural calls can have low pairwise identity")
    ax.text(2, 0.96, "low identity to\nstructural hit (<30%)", fontsize=8, color="gray")
    from matplotlib.lines import Line2D
    ax.legend(handles=[
        Line2D([0], [0], marker="o", ls="", mfc="#1a9850", mec="k", label="correct (≥3 EC levels)"),
        Line2D([0], [0], marker="o", ls="", mfc="#d73027", mec="k", label="incorrect"),
    ], loc="lower right", fontsize=8)
    fig.tight_layout(); fig.savefig(FIGS / "fig_tm_vs_identity.png", dpi=150); plt.close(fig)

    # (c) Precision vs TM threshold.
    calib = summary["precision_vs_tm_threshold"]
    if calib:
        thr = [c["threshold"] for c in calib]
        prec = [c["precision"] for c in calib]
        cov = [c["coverage"] for c in calib]
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.plot(thr, prec, "-o", color="#2b8cbe", label="Precision (≥3-level EC)")
        ax.plot(thr, cov, "-s", color="#999999", label="Coverage (fraction accepted)")
        ax.axvline(0.5, ls="--", c="red", lw=1, label="default threshold 0.5")
        ax.set_xlabel("TM-score acceptance threshold"); ax.set_ylim(0, 1.02)
        ax.set_ylabel("Fraction"); ax.set_title("Calibration of the TM-score cutoff")
        ax.legend(fontsize=8)
        fig.tight_layout(); fig.savefig(FIGS / "fig_precision_threshold.png", dpi=150); plt.close(fig)


# --------------------------------------------------------------------------- #
def main():
    dataset = build_dataset()
    if len(dataset) < 10:
        sys.exit("dataset too small — check network/UniProt access")
    tm = structural_matrix(dataset)
    seq = sequence_matrix(dataset)
    rows = evaluate(dataset, tm, seq)
    summary = summarise(rows)

    import csv
    with open(RESULTS / "per_query.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    (RESULTS / "metrics.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    make_figures(rows, summary)

    a = summary["topN_accuracy_by_ec_level"]
    tz = summary["twilight_zone"]
    print("\n==================== BENCHMARK SUMMARY ====================")
    print(f"proteins: {summary['n_proteins']}  families: {summary['n_families']}")
    print(f"{'EC level':<10}{'structure':>12}{'sequence':>12}")
    for lvl in range(1, 5):
        print(f"{lvl:<10}{a['structure'][lvl]:>12.3f}{a['sequence'][lvl]:>12.3f}")
    print(f"\ntwilight zone: {tz['n_correct_below_cutoff']}/{tz['n_correct_level3']} "
          f"correct structural recoveries are below {tz['identity_cutoff_pct']}% identity "
          f"({tz['fraction_below_cutoff']*100:.0f}%)")
    print("figures + metrics written to benchmark/results/")
    print("==========================================================")


if __name__ == "__main__":
    for s in (sys.stdout, sys.stderr):
        try:
            s.reconfigure(encoding="utf-8", errors="replace")
        except Exception:
            pass
    main()
