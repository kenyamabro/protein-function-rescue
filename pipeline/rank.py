"""Stage 6 — score, rank and export the candidates.

Combines three normalised signals into one composite score:
  * structural_similarity — Foldseek TM-score to a known fold (0-1)
  * model_confidence      — mean pLDDT / 100 (0-1)
  * pocket_quality        — native cavity score, or fpocket druggability (0-1)

Missing components (e.g. no fpocket) are dropped and the remaining weights are
renormalised. Writes results/results.json, results/candidates.csv and web/data.js.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd

from .config import afdb_entry_url, afdb_file_version, afdb_structure_url


def _composite(components: dict[str, float | None], weights: dict[str, float]) -> float:
    num = 0.0
    den = 0.0
    for key, weight in weights.items():
        val = components.get(key)
        if val is None:
            continue
        num += weight * val
        den += weight
    return round(num / den, 4) if den else 0.0


def build_candidates(
    proteins: list[dict[str, Any]],
    hits: dict[str, dict[str, Any]],
    pockets: dict[str, dict[str, Any]],
    config: dict[str, Any],
) -> list[dict[str, Any]]:
    """Assemble ranked candidate records.

    A candidate = an *originally unannotated* protein that has a confident
    structural match to a *known* fold. Those are the rescued functions.
    """
    weights = config["ranking"]["weights"]
    version = afdb_file_version(config)
    candidates: list[dict[str, Any]] = []

    for p in proteins:
        acc = p["accession"]
        hit = hits.get(acc)
        if hit is None:
            continue  # no structural match to a known fold -> not a candidate
        if not p.get("originally_unannotated", False):
            continue  # already characterised -> nothing to rescue

        pocket = pockets.get(acc, {})
        components = {
            "structural_similarity": hit["tm_score"],
            "model_confidence": p["mean_plddt"] / 100.0,
            "pocket_quality": pocket.get("top_druggability"),
        }
        composite = _composite(components, weights)

        candidates.append(
            {
                "accession": acc,
                "name": p.get("name", "Uncharacterized protein"),
                "length": p["length"],
                "mean_plddt": p["mean_plddt"],
                "originally_unannotated": True,
                "best_structural_match": {
                    "target": hit["target"],
                    "tm_score": hit["tm_score"],
                    "tm_score_query": hit.get("tm_score_query"),
                    "tm_score_reference": hit.get("tm_score_reference"),
                    "query_length": hit.get("query_length", p.get("length")),
                    "reference_length": hit.get("reference_length"),
                    "rmsd": hit.get("rmsd"),
                    "evalue": hit["evalue"],
                    "description": hit["description"],
                },
                "pocket": {
                    "n_pockets": pocket.get("n_pockets", 0),
                    "top_druggability": pocket.get("top_druggability"),
                    "method": pocket.get("method", "fpocket"),
                },
                "scores": {
                    "structural_similarity": round(components["structural_similarity"], 4),
                    "model_confidence": round(components["model_confidence"], 4),
                    "pocket_quality": components["pocket_quality"],
                    "composite": composite,
                },
                "afdb_entry_url": afdb_entry_url(acc),
                "structure_cif_url": afdb_structure_url(acc, "cif", version),
                "structure_pdb_url": afdb_structure_url(acc, "pdb", version),
            }
        )

    candidates.sort(key=lambda c: c["scores"]["composite"], reverse=True)
    return candidates


def export(candidates: list[dict[str, Any]], config: dict[str, Any]) -> None:
    """Write JSON, CSV and the web viewer's data.js."""
    results_dir = Path(config["paths"]["results_dir"])
    results_dir.mkdir(parents=True, exist_ok=True)

    payload = {
        "organism": config["organism"]["name"],
        "proteome_id": config["organism"]["proteome_id"],
        "n_candidates": len(candidates),
        "candidates": candidates,
    }

    json_path = results_dir / "results.json"
    json_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    # Flat CSV for spreadsheets / supplementary data.
    rows = [
        {
            "accession": c["accession"],
            "name": c["name"],
            "length": c["length"],
            "mean_plddt": c["mean_plddt"],
            "match_target": c["best_structural_match"]["target"],
            "match_description": c["best_structural_match"]["description"],
            "tm_score": c["best_structural_match"]["tm_score"],
            "evalue": c["best_structural_match"]["evalue"],
            "n_pockets": c["pocket"]["n_pockets"],
            "top_druggability": c["pocket"]["top_druggability"],
            "composite_score": c["scores"]["composite"],
            "afdb_entry_url": c["afdb_entry_url"],
        }
        for c in candidates
    ]
    pd.DataFrame(rows).to_csv(results_dir / "candidates.csv", index=False)

    # data.js so web/index.html works by double-click (no server / CORS needed).
    web_data = Path(config["paths"]["web_data"])
    web_data.parent.mkdir(parents=True, exist_ok=True)
    web_data.write_text(
        "window.RESCUE_DATA = " + json.dumps(payload) + ";\n", encoding="utf-8"
    )

    print(
        f"[rank] wrote {len(candidates)} candidates -> {json_path}, "
        f"{results_dir / 'candidates.csv'}, {web_data}"
    )
