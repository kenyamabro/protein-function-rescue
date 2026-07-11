"""Build a realistic candidate pool for the weight-sensitivity analysis.

The composite score is used to rank *rescue candidates* — unannotated proteins
that passed the TM-score threshold. Those have a genuine spread of component
values (structural similarity 0.5-0.9, variable confidence, and pockets that are
present or absent). This script assembles such a pool by running the real search
+ scoring on a batch of M. tuberculosis hypothetical proteins, so the weight
sensitivity is measured on the distribution the composite actually sees — not on
a set of well-folded enzymes whose scores are all near-maximal.

Output: benchmark/results/candidate_pool.json
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import requests
from tmtools import tm_align
from tqdm import tqdm

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from pipeline import download as download_mod            # noqa: E402
from pipeline import tmalign as tmalign_mod              # noqa: E402
from pipeline.config import load_config                  # noqa: E402
from pipeline.parse import parse_structure               # noqa: E402
from pipeline.pocket_geom import detect_pockets          # noqa: E402
from pipeline.structio import load_ca, load_atoms        # noqa: E402

N_CANDIDATES = 40         # target pool size (stop once reached)
MAX_PROCESSED = 90        # hard cap on accessions examined (bounds runtime)
REF_SAMPLE = 45           # subsample of the reference library (speed; spread is preserved)
LEN_MAX = 460
MIN_TM = 0.45             # slightly below the pipeline cutoff to populate the low end
POOL_DIR = ROOT / "benchmark" / "candidate_structures"
OUT = ROOT / "benchmark" / "results" / "candidate_pool.json"
UNIPROT = "https://rest.uniprot.org/uniprotkb/search"


def fetch_hypothetical_accessions(n: int) -> list[str]:
    params = {
        "query": (
            "(proteome:UP000001584) AND (protein_name:uncharacterized) "
            f"AND (length:[80 TO {LEN_MAX}])"
        ),
        "fields": "accession",
        "format": "tsv",
        "size": n * 2,
    }
    r = requests.get(UNIPROT, params=params, timeout=90)
    r.raise_for_status()
    return [ln.split("\t")[0] for ln in r.text.strip().splitlines()[1:]]


def main():
    for s in (sys.stdout, sys.stderr):
        try:
            s.reconfigure(encoding="utf-8", errors="replace")
        except Exception:
            pass

    config = load_config(ROOT / "config.yaml")
    POOL_DIR.mkdir(parents=True, exist_ok=True)

    # Reference library of known enzymes (reuse the pipeline's cached one),
    # subsampled for speed — a realistic score spread does not need all 119.
    import random
    random.seed(0)
    manifest = tmalign_mod.ensure_reference_library(config)
    random.shuffle(manifest)
    refs = []
    for entry in manifest[:REF_SAMPLE]:
        coords, seq = load_ca(entry["path"])
        if len(seq) >= 10:
            refs.append((coords, seq))
    print(f"[pool] {len(refs)} reference structures (subsampled)")

    accs = fetch_hypothetical_accessions(MAX_PROCESSED)
    pool = []
    for processed, acc in enumerate(tqdm(accs, desc="candidates")):
        if len(pool) >= N_CANDIDATES or processed >= MAX_PROCESSED:
            break
        try:
            path = download_mod.download_single(acc, POOL_DIR)
        except Exception:
            continue
        qc, qs = load_ca(path)
        if len(qs) < 20:
            continue
        # Best structural similarity to the reference library.
        best = 0.0
        for (rc, rs) in refs:
            try:
                res = tm_align(qc, rc, qs, rs)
                best = max(best, res.tm_norm_chain1, res.tm_norm_chain2)
            except Exception:
                continue
        if best < MIN_TM:
            continue  # not a candidate (no confident structural hit)
        rec = parse_structure(path)
        plddt = rec["mean_plddt"] / 100.0 if rec else 0.0
        coords, radii = load_atoms(path)
        pk = detect_pockets(coords, radii, 1.2, 5.0, 5, 30)
        pool.append({
            "accession": acc,
            "structural_similarity": round(float(best), 4),
            "model_confidence": round(float(plddt), 4),
            "pocket_quality": (round(pk["top_druggability"], 4)
                               if pk.get("top_druggability") is not None else None),
        })

    OUT.write_text(json.dumps(pool, indent=2), encoding="utf-8")
    ss = np.array([p["structural_similarity"] for p in pool])
    npocket = sum(1 for p in pool if p["pocket_quality"] is None)
    print(f"\n[pool] {len(pool)} candidates written to {OUT}")
    print(f"[pool] structural_similarity: min {ss.min():.3f} med "
          f"{np.median(ss):.3f} max {ss.max():.3f} std {ss.std():.3f}")
    print(f"[pool] candidates without a detected pocket: {npocket}")


if __name__ == "__main__":
    main()
