"""Configuration loading and small shared helpers."""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any

import yaml

# AlphaFold DB file / URL conventions (v4).
AFDB_FILES_BASE = "https://alphafold.ebi.ac.uk/files"
# Matches both the AFDB tar names (AF-P9WGR1-F1-model_v6.cif.gz) and our
# version-agnostic local names (AF-P9WGR1.cif). UniProt accessions are
# alphanumeric with no internal '-', so the capture stops at the first '-'/'.'.
_ACCESSION_RE = re.compile(r"AF-([A-Z0-9]+)")


def load_config(path: str | Path = "config.yaml") -> dict[str, Any]:
    """Read the YAML config into a plain dict."""
    with open(path, "r", encoding="utf-8") as fh:
        return yaml.safe_load(fh)


def ensure_dirs(config: dict[str, Any]) -> None:
    """Create the data/results directories referenced by the config."""
    for key in ("data_dir", "structures_dir", "results_dir"):
        Path(config["paths"][key]).mkdir(parents=True, exist_ok=True)
    Path(config["paths"]["web_data"]).parent.mkdir(parents=True, exist_ok=True)


def accession_from_filename(name: str) -> str | None:
    """Extract the UniProt accession from an AlphaFold filename.

    e.g. 'AF-P9WGE7-F1-model_v4.cif' -> 'P9WGE7'
    """
    match = _ACCESSION_RE.search(name)
    return match.group(1) if match else None


def afdb_structure_url(accession: str, fmt: str = "cif", version: int = 6) -> str:
    """Direct URL to a single AlphaFold model file (cif or pdb).

    ``version`` is the AlphaFold DB release (currently 6). It is configurable via
    ``afdb.file_version`` in config.yaml so a future release bump is a one-liner.
    """
    return f"{AFDB_FILES_BASE}/AF-{accession}-F1-model_v{version}.{fmt}"


def afdb_file_version(config: dict[str, Any]) -> int:
    """Read the configured AlphaFold file version (defaults to 6)."""
    return int(config.get("afdb", {}).get("file_version", 6))


def afdb_entry_url(accession: str) -> str:
    """Human-facing AlphaFold DB entry page."""
    return f"https://alphafold.ebi.ac.uk/entry/{accession}"
