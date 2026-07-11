"""Functional-residue verification for a structural match.

Global fold similarity is necessary but not sufficient evidence for shared
function. This module strengthens a candidate->reference assignment by checking
whether the reference protein's *annotated functional residues* — catalytic
(active site), metal-binding, and ligand-binding residues from UniProt — are
structurally conserved in the candidate.

Method: superpose the candidate onto the reference using the TM-align rotation/
translation, map each annotated reference residue to its nearest candidate Ca,
and report whether a candidate residue occupies that position (within a distance
cutoff) and whether its identity is conserved. A high fraction of conserved
active-site/binding residues is much stronger evidence than fold alone.

Assumes AlphaFold models (single chain, residues numbered 1..L contiguously),
which holds for AlphaFold DB entries.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import requests
from tmtools import tm_align

from . import download as download_mod
from .structio import load_ca

UNIPROT_ENTRY = "https://rest.uniprot.org/uniprotkb/{acc}.json"

# UniProt feature types that mark a functional residue.
_FUNCTIONAL_TYPES = {"Active site", "Binding site", "Metal binding", "Site"}

# Conservative substitution groups (for "chemically similar" scoring).
_GROUPS = [
    set("AVLIM"), set("FYW"), set("ST"), set("NQ"),
    set("DE"), set("KRH"), set("GP"), set("C"),
]


def _same_group(a: str, b: str) -> bool:
    if a == b:
        return True
    return any(a in g and b in g for g in _GROUPS)


def fetch_functional_sites(accession: str) -> list[dict[str, Any]]:
    """Return annotated functional residues for a UniProt accession.

    Each item: {num (1-based), type, description, ligand}.
    """
    resp = requests.get(UNIPROT_ENTRY.format(acc=accession), timeout=60)
    resp.raise_for_status()
    data = resp.json()
    sites: list[dict[str, Any]] = []
    for feat in data.get("features", []):
        ftype = feat.get("type", "")
        if ftype not in _FUNCTIONAL_TYPES:
            continue
        loc = feat.get("location", {})
        start = (loc.get("start") or {}).get("value")
        end = (loc.get("end") or {}).get("value")
        if start is None:
            continue
        ligand = ""
        if isinstance(feat.get("ligand"), dict):
            ligand = feat["ligand"].get("name", "")
        desc = feat.get("description", "")
        for num in range(int(start), int(end) + 1):
            sites.append(
                {"num": num, "type": ftype, "description": desc, "ligand": ligand}
            )
    return sites


def _best_superposition(qc: np.ndarray, rc: np.ndarray, qs: str, rs: str):
    """TM-align the candidate (q) onto the reference (r); return transformed qc.

    Tries both transform conventions and keeps the one giving the smaller mean
    nearest-neighbour distance, so the result is robust to the u/t convention.
    """
    res = tm_align(qc, rc, qs, rs)
    u, t = np.asarray(res.u), np.asarray(res.t)
    cand_a = (qc @ u.T) + t
    cand_b = (qc @ u) + t

    def mean_nn(x):
        # mean distance from each reference Ca to the nearest transformed cand Ca
        d = np.linalg.norm(rc[:, None, :] - x[None, :, :], axis=2)
        return float(d.min(axis=1).mean())

    return (cand_a if mean_nn(cand_a) <= mean_nn(cand_b) else cand_b), res


def verify_match(
    candidate_path: str | Path,
    reference_acc: str,
    reference_path: str | Path,
    cutoff: float = 4.0,
) -> dict[str, Any]:
    """Check conservation of the reference's functional residues in the candidate."""
    qc, qs = load_ca(candidate_path)          # candidate (query)
    rc, rs = load_ca(reference_path)           # reference (known function)
    if len(qs) < 10 or len(rs) < 10:
        return {"error": "structure too short", "reference": reference_acc}

    sites = fetch_functional_sites(reference_acc)
    if not sites:
        return {
            "reference": reference_acc,
            "n_functional_residues": 0,
            "note": "no annotated functional residues in the reference",
        }

    cand_t, aln = _best_superposition(qc, rc, qs, rs)

    per_site = []
    n_aligned = n_identical = n_similar = 0
    for s in sites:
        idx_r = s["num"] - 1
        if idx_r < 0 or idx_r >= len(rc):
            continue
        ref_xyz = rc[idx_r]
        d = np.linalg.norm(cand_t - ref_xyz, axis=1)
        j = int(d.argmin())
        dist = float(d[j])
        aligned = dist <= cutoff
        ref_res = rs[idx_r]
        cand_res = qs[j] if aligned else "-"
        identical = aligned and cand_res == ref_res
        similar = aligned and _same_group(cand_res, ref_res)
        n_aligned += aligned
        n_identical += identical
        n_similar += similar
        per_site.append({
            "ref_residue": f"{ref_res}{s['num']}",
            "type": s["type"],
            "ligand": s["ligand"],
            "description": s["description"],
            "candidate_residue": (f"{cand_res}{j + 1}" if aligned else None),
            "distance_A": round(dist, 2),
            "position_conserved": aligned,
            "identity_conserved": identical,
            "chemistry_conserved": similar,
        })

    n = len(per_site)
    return {
        "reference": reference_acc,
        "n_functional_residues": n,
        "n_position_conserved": n_aligned,
        "n_identity_conserved": n_identical,
        "n_chemistry_conserved": n_similar,
        "frac_position_conserved": round(n_aligned / n, 3) if n else 0.0,
        "frac_identity_conserved": round(n_identical / n, 3) if n else 0.0,
        "tm_score": round(max(aln.tm_norm_chain1, aln.tm_norm_chain2), 4),
        "cutoff_A": cutoff,
        "sites": per_site,
    }


def verify_results(
    results_json: str | Path,
    ref_dir: str | Path = "references",
    cutoff: float = 4.0,
) -> list[dict[str, Any]]:
    """Verify functional residues for every candidate in a results.json file."""
    results = json.loads(Path(results_json).read_text(encoding="utf-8"))
    ref_dir = Path(ref_dir)
    out = []
    for c in results.get("candidates", []):
        acc = c["accession"]
        ref_acc = c["best_structural_match"]["target"]
        cand_path = _candidate_structure(c, ref_dir.parent)
        try:
            ref_path = download_mod.download_single(ref_acc, ref_dir)
        except Exception as exc:
            out.append({"accession": acc, "reference": ref_acc, "error": str(exc)})
            continue
        summary = verify_match(cand_path, ref_acc, ref_path, cutoff)
        summary["accession"] = acc
        summary["candidate_name"] = c.get("name", "")
        summary["match_description"] = c["best_structural_match"]["description"]
        out.append(summary)
    return out


def _candidate_structure(candidate: dict[str, Any], project_root: Path) -> Path:
    """Locate (or fetch) the candidate structure file."""
    acc = candidate["accession"]
    local = project_root / "data" / "structures" / f"AF-{acc}.cif"
    if local.exists():
        return local
    return download_mod.download_single(acc, project_root / "data" / "structures")


if __name__ == "__main__":
    import sys
    for s in (sys.stdout, sys.stderr):
        try:
            s.reconfigure(encoding="utf-8", errors="replace")
        except Exception:
            pass

    root = Path(__file__).resolve().parent.parent
    summaries = verify_results(root / "results" / "results.json")
    out_path = root / "results" / "functional_residues.json"
    out_path.write_text(json.dumps(summaries, indent=2), encoding="utf-8")

    print("\n==== FUNCTIONAL-RESIDUE VERIFICATION ====")
    for s in summaries:
        if s.get("error"):
            print(f"{s['accession']}: error ({s['error']})")
            continue
        if s.get("n_functional_residues", 0) == 0:
            print(f"{s['accession']} -> {s['reference']}: no annotated sites")
            continue
        print(f"{s['accession']} -> {s['reference']}  "
              f"({s['match_description'][:50]})")
        print(f"   functional residues: {s['n_functional_residues']}  "
              f"position-conserved: {s['n_position_conserved']} "
              f"({s['frac_position_conserved']*100:.0f}%)  "
              f"identity-conserved: {s['n_identity_conserved']} "
              f"({s['frac_identity_conserved']*100:.0f}%)")
    print(f"\nwritten to {out_path}")
