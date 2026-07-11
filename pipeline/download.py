"""Stage 1 — download an organism's AlphaFold DB structures.

Fetches the proteome `.tar` from the AlphaFold DB FTP and extracts the
per-protein `.cif` models into ``paths.structures_dir``. Use ``--limit`` to grab
only the first N structures for a quick test run.
"""

from __future__ import annotations

import gzip
import shutil
import tarfile
from pathlib import Path
from typing import Any

import requests
from tqdm import tqdm

from .config import accession_from_filename, afdb_structure_url

_CHUNK = 1 << 20  # 1 MiB
AFDB_API = "https://alphafold.ebi.ac.uk/api/prediction"


def _resolve_cif_url(accession: str) -> str:
    """Ask the AFDB API for a structure's current cif URL (version-proof)."""
    resp = requests.get(f"{AFDB_API}/{accession}", timeout=60)
    resp.raise_for_status()
    payload = resp.json()
    if payload and payload[0].get("cifUrl"):
        return payload[0]["cifUrl"]
    return afdb_structure_url(accession, "cif")


def _stream_download(url: str, dest: Path) -> None:
    """Download ``url`` to ``dest`` with a progress bar."""
    with requests.get(url, stream=True, timeout=120) as resp:
        resp.raise_for_status()
        total = int(resp.headers.get("content-length", 0))
        with open(dest, "wb") as fh, tqdm(
            total=total, unit="B", unit_scale=True, desc=dest.name, leave=False
        ) as bar:
            for chunk in resp.iter_content(chunk_size=_CHUNK):
                fh.write(chunk)
                bar.update(len(chunk))


def download_single(accession: str, structures_dir: Path) -> Path:
    """Fetch one AlphaFold `.cif` model by UniProt accession (helper / fallback)."""
    structures_dir.mkdir(parents=True, exist_ok=True)
    dest = structures_dir / f"AF-{accession}.cif"
    if dest.exists():
        return dest
    resp = requests.get(_resolve_cif_url(accession), timeout=60)
    resp.raise_for_status()
    dest.write_bytes(resp.content)
    return dest


def download_proteome(config: dict[str, Any], limit: int | None = None) -> list[Path]:
    """Download and extract the configured proteome. Returns extracted .cif paths."""
    data_dir = Path(config["paths"]["data_dir"])
    structures_dir = Path(config["paths"]["structures_dir"])
    data_dir.mkdir(parents=True, exist_ok=True)
    structures_dir.mkdir(parents=True, exist_ok=True)

    url = config["organism"]["proteome_tar_url"]
    tar_path = data_dir / Path(url).name

    if not tar_path.exists():
        print(f"[download] fetching proteome tar: {url}")
        _stream_download(url, tar_path)
    else:
        print(f"[download] using cached tar: {tar_path}")

    extracted: list[Path] = []
    with tarfile.open(tar_path, "r") as tar:
        # Keep only the structure models (skip PAE json, etc.).
        members = [m for m in tar.getmembers() if m.name.endswith(".cif.gz")]
        members.sort(key=lambda m: m.name)
        if limit is not None:
            members = members[:limit]

        for member in tqdm(members, desc="extract", unit="struct"):
            acc = accession_from_filename(member.name)
            if acc is None:
                continue
            out_path = structures_dir / f"AF-{acc}.cif"
            if out_path.exists():
                extracted.append(out_path)
                continue
            src = tar.extractfile(member)
            if src is None:
                continue
            with gzip.open(src, "rb") as gz, open(out_path, "wb") as fh:
                shutil.copyfileobj(gz, fh)
            extracted.append(out_path)

    print(f"[download] {len(extracted)} structures available in {structures_dir}")
    return extracted
