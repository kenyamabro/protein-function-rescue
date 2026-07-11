"""Structural search via TM-align (native Windows alternative to Foldseek).

TM-align is the gold-standard structural aligner; the TM-score it returns is the
same metric Foldseek reports. We align each *unannotated* query protein against a
reference library of *functionally characterised* proteins and keep the best hit.

Foldseek is faster and scales to millions of targets; TM-align is exhaustive and
per-pair exact. For a laptop project against a curated reference set they give the
same kind of result. (Set structural_search.backend: foldseek in config to use
Foldseek instead once you have a Linux environment.)
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import requests
from tmtools import tm_align
from tqdm import tqdm

from . import download as download_mod
from .config import accession_from_filename
from .structio import load_ca

UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search"


# ---------------------------------------------------------------------------
# Reference library of known-function proteins
# ---------------------------------------------------------------------------
def _fetch_reference_accessions(count: int) -> list[dict[str, str]]:
    """Curated set of reviewed, experimentally-evidenced enzymes (diverse folds)."""
    params = {
        "query": "(reviewed:true) AND (ec:*) AND (existence:1) AND (fragment:false)",
        "fields": "accession,protein_name,ec",
        "format": "tsv",
        "size": min(count, 500),
        "sort": "annotation_score desc",
    }
    resp = requests.get(UNIPROT_SEARCH, params=params, timeout=90)
    resp.raise_for_status()
    refs: list[dict[str, str]] = []
    for row in resp.text.strip().splitlines()[1:]:
        cols = (row.split("\t") + ["", "", ""])[:3]
        acc, name, ec = cols[0], cols[1], cols[2]
        refs.append({"accession": acc, "name": name or acc, "ec": ec})
    return refs[:count]


def ensure_reference_library(config: dict[str, Any]) -> list[dict[str, Any]]:
    """Download the reference structures (once) and return their manifest.

    Manifest entries: {accession, name, ec, path}. Cached in references/manifest.json.
    """
    cfg = config["structural_search"]
    ref_dir = Path(cfg["reference_dir"])
    ref_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = ref_dir / "manifest.json"

    if manifest_path.exists():
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        if len(manifest) >= 1:
            print(f"[search] using cached reference library ({len(manifest)} proteins)")
            return manifest

    print("[search] building reference library of known-function proteins…")
    wanted = _fetch_reference_accessions(int(cfg["reference_count"]))
    manifest: list[dict[str, Any]] = []
    for ref in tqdm(wanted, desc="ref-download", unit="struct"):
        try:
            path = download_mod.download_single(ref["accession"], ref_dir)
        except Exception:
            continue  # no AlphaFold model for this accession — skip
        manifest.append({**ref, "path": str(path)})

    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"[search] reference library ready: {len(manifest)} proteins")
    return manifest


# ---------------------------------------------------------------------------
# Search
# ---------------------------------------------------------------------------
def _load_reference_coords(manifest: list[dict[str, Any]]) -> list[dict[str, Any]]:
    refs = []
    for entry in tqdm(manifest, desc="ref-load", unit="struct"):
        coords, seq = load_ca(entry["path"])
        if len(seq) < 10:
            continue
        refs.append({**entry, "coords": coords, "seq": seq})
    return refs


def _query_accessions(config: dict[str, Any]) -> set[str] | None:
    """Restrict queries to originally-unannotated proteins, if annotation exists."""
    proteins_path = Path(config["paths"]["results_dir"]) / "proteins.json"
    if not proteins_path.exists():
        return None
    proteins = json.loads(proteins_path.read_text(encoding="utf-8"))
    unannotated = {
        p["accession"] for p in proteins if p.get("originally_unannotated")
    }
    return unannotated or None


def run_tmalign(config: dict[str, Any]) -> dict[str, dict[str, Any]]:
    """Return {query_accession: {target, evalue, tm_score, description}}."""
    cfg = config["structural_search"]
    min_tm = float(cfg["min_tmscore"])

    manifest = ensure_reference_library(config)
    refs = _load_reference_coords(manifest)
    if not refs:
        raise RuntimeError("reference library is empty — cannot search")

    structures_dir = Path(config["paths"]["structures_dir"])
    query_paths = sorted(structures_dir.glob("AF-*.cif"))
    only = _query_accessions(config)
    if only is not None:
        query_paths = [p for p in query_paths if accession_from_filename(p.name) in only]

    ref_accs = {r["accession"] for r in refs}
    best: dict[str, dict[str, Any]] = {}

    for qpath in tqdm(query_paths, desc="tm-align", unit="query"):
        acc = accession_from_filename(qpath.name)
        if acc is None or acc in ref_accs:
            continue
        qc, qs = load_ca(qpath)
        if len(qs) < 10:
            continue

        top = None
        for ref in refs:
            try:
                res = tm_align(qc, ref["coords"], qs, ref["seq"])
            except Exception:
                continue
            # Higher of the two length-normalisations = best mutual fold match.
            tm = max(res.tm_norm_chain1, res.tm_norm_chain2)
            if top is None or tm > top["tm_score"]:
                top = {
                    "target": ref["accession"],
                    "tm_score": round(float(tm), 4),
                    "rmsd": round(float(res.rmsd), 3),
                    "evalue": None,  # TM-align has no e-value; kept for schema parity
                    "description": ref["name"]
                    + (f" (EC {ref['ec']})" if ref.get("ec") else ""),
                }
        if top is not None and top["tm_score"] >= min_tm:
            best[acc] = top

    print(f"[search] {len(best)} queries with a structural hit (TM >= {min_tm})")
    return best
