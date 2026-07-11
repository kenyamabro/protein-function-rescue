"""Shared structure I/O built on gemmi (robust for AlphaFold mmCIF files).

Used by the native Windows stages (TM-align search, geometric pockets). gemmi
reads AlphaFold `.cif` reliably here, where Biopython's CIF parser did not.
"""

from __future__ import annotations

from pathlib import Path

import gemmi
import numpy as np

# Van der Waals radii (Angstrom) for the elements seen in protein models.
_VDW = {"C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80, "P": 1.80, "H": 1.20}
_DEFAULT_VDW = 1.70


def one_letter(res_name: str) -> str:
    info = gemmi.find_tabulated_residue(res_name)
    if info is not None and info.is_amino_acid():
        code = info.one_letter_code.upper()
        return code if code.isalpha() else "X"
    return "X"


def load_ca(path: str | Path) -> tuple[np.ndarray, str]:
    """Return (CA coordinates [N,3], one-letter sequence) for a structure."""
    st = gemmi.read_structure(str(path))
    coords: list[list[float]] = []
    seq: list[str] = []
    if len(st) == 0:
        return np.zeros((0, 3)), ""
    for chain in st[0]:
        for res in chain:
            for atom in res:
                if atom.name == "CA":
                    p = atom.pos
                    coords.append([p.x, p.y, p.z])
                    seq.append(one_letter(res.name))
                    break
    return np.asarray(coords, dtype=float), "".join(seq)


def load_atoms(path: str | Path) -> tuple[np.ndarray, np.ndarray]:
    """Return (all heavy-atom coordinates [M,3], their vdW radii [M])."""
    st = gemmi.read_structure(str(path))
    coords: list[list[float]] = []
    radii: list[float] = []
    if len(st) == 0:
        return np.zeros((0, 3)), np.zeros((0,))
    for chain in st[0]:
        for res in chain:
            for atom in res:
                elem = atom.element.name if atom.element else atom.name[:1]
                if elem == "H":
                    continue  # heavy atoms only
                p = atom.pos
                coords.append([p.x, p.y, p.z])
                radii.append(_VDW.get(elem, _DEFAULT_VDW))
    return np.asarray(coords, dtype=float), np.asarray(radii, dtype=float)
