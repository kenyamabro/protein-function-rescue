"""Stage 2 — parse structures for confidence, sequence and length.

AlphaFold stores the per-residue pLDDT confidence in the B-factor column, so the
mean pLDDT is just the average CA B-factor. Structures below the configured
``filters.min_mean_plddt`` are dropped here.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import gemmi
from tqdm import tqdm

from .config import accession_from_filename


def _one_letter(res_name: str) -> str:
    info = gemmi.find_tabulated_residue(res_name)
    if info is not None and info.is_amino_acid():
        code = info.one_letter_code.upper()
        return code if code.isalpha() else "X"
    return "X"


def parse_structure(path: Path) -> dict[str, Any] | None:
    """Return {accession, length, mean_plddt, sequence} for one model, or None."""
    acc = accession_from_filename(path.name)
    if acc is None:
        return None

    st = gemmi.read_structure(str(path))
    if len(st) == 0:
        return None
    model = st[0]

    plddts: list[float] = []
    sequence: list[str] = []
    for chain in model:
        for res in chain:
            ca = None
            for atom in res:
                if atom.name == "CA":
                    ca = atom
                    break
            if ca is None:
                continue
            plddts.append(ca.b_iso)
            sequence.append(_one_letter(res.name))

    if not plddts:
        return None

    return {
        "accession": acc,
        "length": len(sequence),
        "mean_plddt": round(sum(plddts) / len(plddts), 2),
        "sequence": "".join(sequence),
    }


def parse_all(config: dict[str, Any]) -> list[dict[str, Any]]:
    """Parse every structure in ``structures_dir``, keeping only confident models."""
    structures_dir = Path(config["paths"]["structures_dir"])
    min_plddt = float(config["filters"]["min_mean_plddt"])

    paths = sorted(structures_dir.glob("AF-*.cif"))
    proteins: list[dict[str, Any]] = []
    dropped = 0
    for path in tqdm(paths, desc="parse", unit="struct"):
        rec = parse_structure(path)
        if rec is None:
            continue
        if rec["mean_plddt"] < min_plddt:
            dropped += 1
            continue
        proteins.append(rec)

    print(
        f"[parse] kept {len(proteins)} structures "
        f"(dropped {dropped} below pLDDT {min_plddt})"
    )
    return proteins
