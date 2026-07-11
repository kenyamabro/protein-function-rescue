"""Stage 5 (optional) — binding-pocket detection with fpocket.

For each candidate we convert its `.cif` to PDB, run ``fpocket``, and read back
the number of detected pockets and the best druggability score. A prominent,
druggable cavity is extra evidence that an unannotated protein is an enzyme.

fpocket is optional: if the binary is missing (or --skip-pocket is passed), the
pocket component is simply left out of the ranking.
"""

from __future__ import annotations

import re
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Any

import gemmi

_DRUG_RE = re.compile(r"Druggability Score\s*:\s*([0-9.]+)")
_POCKET_RE = re.compile(r"^Pocket\s+\d+\s*:", re.MULTILINE)


class PocketUnavailable(RuntimeError):
    pass


def _resolve_binary(binary: str) -> str:
    found = shutil.which(binary)
    if found is None and not Path(binary).exists():
        raise PocketUnavailable(
            f"fpocket binary '{binary}' not found. Install it or run with "
            f"--skip-pocket. See docs/SETUP_TOOLS.md."
        )
    return found or binary


def _cif_to_pdb(cif_path: Path, pdb_path: Path) -> None:
    st = gemmi.read_structure(str(cif_path))
    st.setup_entities()
    pdb_path.write_text(st.make_pdb_string())


def _analyse_one(binary: str, cif_path: Path) -> dict[str, Any]:
    """Run fpocket on a single structure, returning pocket stats."""
    with tempfile.TemporaryDirectory() as tmp:
        tmp_dir = Path(tmp)
        pdb_path = tmp_dir / (cif_path.stem + ".pdb")
        _cif_to_pdb(cif_path, pdb_path)

        proc = subprocess.run(
            [binary, "-f", str(pdb_path)], capture_output=True, text=True
        )
        if proc.returncode != 0:
            return {"n_pockets": 0, "top_druggability": None}

        info = tmp_dir / (pdb_path.stem + "_out") / (pdb_path.stem + "_info.txt")
        if not info.exists():
            return {"n_pockets": 0, "top_druggability": None}

        text = info.read_text(encoding="utf-8", errors="ignore")
        drugs = [float(m) for m in _DRUG_RE.findall(text)]
        n_pockets = len(_POCKET_RE.findall(text))
        return {
            "n_pockets": n_pockets,
            "top_druggability": round(max(drugs), 3) if drugs else None,
        }


def run_pockets(
    config: dict[str, Any], accessions: list[str]
) -> dict[str, dict[str, Any]]:
    """Run fpocket for the given candidate accessions. {accession: pocket_stats}."""
    binary = _resolve_binary(config["pocket"]["binary"])
    structures_dir = Path(config["paths"]["structures_dir"])

    out: dict[str, dict[str, Any]] = {}
    for acc in accessions:
        cif = structures_dir / f"AF-{acc}.cif"
        if not cif.exists():
            continue
        try:
            out[acc] = _analyse_one(binary, cif)
        except Exception as exc:  # pragma: no cover - fpocket edge cases
            print(f"[pocket] {acc}: skipped ({exc})")
            out[acc] = {"n_pockets": 0, "top_druggability": None}

    print(f"[pocket] analysed {len(out)} candidates")
    return out
