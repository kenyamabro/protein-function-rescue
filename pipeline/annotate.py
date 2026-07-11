"""Stage 3 — annotate proteins with their current UniProt description.

For each accession we ask UniProt for the protein name, evidence level and
keywords, then flag the ones that are still *unannotated* (uncharacterized /
hypothetical / DUF). Those are the targets worth rescuing structurally.
"""

from __future__ import annotations

from typing import Any, Iterable

import requests
from tqdm import tqdm

UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search"
_BATCH = 100


def _chunks(items: list[str], size: int) -> Iterable[list[str]]:
    for i in range(0, len(items), size):
        yield items[i : i + size]


def _is_unannotated(name: str, existence: str, markers: list[str]) -> bool:
    low = name.lower()
    if any(marker in low for marker in markers):
        return True
    # UniProt "Predicted" / "Uncertain" existence = no experimental support.
    return existence.lower() in {"predicted", "uncertain"}


def fetch_annotations(
    accessions: list[str], markers: list[str]
) -> dict[str, dict[str, Any]]:
    """Return {accession: {name, existence, keywords, originally_unannotated}}."""
    out: dict[str, dict[str, Any]] = {}
    fields = "accession,protein_name,protein_existence,keyword"

    for batch in tqdm(
        list(_chunks(accessions, _BATCH)), desc="annotate", unit="batch"
    ):
        query = " OR ".join(f"accession:{acc}" for acc in batch)
        params = {
            "query": query,
            "fields": fields,
            "format": "tsv",
            "size": _BATCH,
        }
        resp = requests.get(UNIPROT_SEARCH, params=params, timeout=60)
        resp.raise_for_status()

        lines = resp.text.strip().splitlines()
        for row in lines[1:]:  # skip header
            cols = row.split("\t")
            # accession, protein names, existence, keywords
            while len(cols) < 4:
                cols.append("")
            acc, name, existence, keywords = cols[0], cols[1], cols[2], cols[3]
            out[acc] = {
                "name": name or "Uncharacterized protein",
                "existence": existence,
                "keywords": keywords,
                "originally_unannotated": _is_unannotated(name, existence, markers),
            }

    # Any accession UniProt didn't return: treat as unknown/unannotated.
    for acc in accessions:
        out.setdefault(
            acc,
            {
                "name": "Uncharacterized protein",
                "existence": "",
                "keywords": "",
                "originally_unannotated": True,
            },
        )
    return out


def annotate_proteins(
    proteins: list[dict[str, Any]], config: dict[str, Any]
) -> list[dict[str, Any]]:
    """Merge UniProt annotation into the parsed protein records (in place)."""
    markers = [m.lower() for m in config["filters"]["unannotated_markers"]]
    accessions = [p["accession"] for p in proteins]
    ann = fetch_annotations(accessions, markers)

    n_unannotated = 0
    for p in proteins:
        info = ann.get(p["accession"], {})
        p["name"] = info.get("name", "Uncharacterized protein")
        p["existence"] = info.get("existence", "")
        p["keywords"] = info.get("keywords", "")
        p["originally_unannotated"] = bool(info.get("originally_unannotated", True))
        n_unannotated += p["originally_unannotated"]

    print(
        f"[annotate] {n_unannotated}/{len(proteins)} proteins are currently "
        f"unannotated (the rescue candidates)"
    )
    return proteins
