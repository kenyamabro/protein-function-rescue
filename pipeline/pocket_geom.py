"""Geometric binding-pocket detection (native Windows alternative to fpocket).

A compact LIGSITE-style algorithm:
  1. Rasterise the protein onto a 3D grid (voxels within vdW radius = "protein").
  2. For each empty (solvent) voxel, count how many of 7 scan directions have
     protein on *both* sides — a "protein-solvent-protein" (PSP) sandwich.
  3. Voxels sandwiched in >= `min_directions` directions are pocket voxels.
  4. Cluster adjacent pocket voxels; each cluster is a candidate pocket.

The reported "cavity_score" (0-1) is a geometric proxy for druggability derived
from the largest pocket's volume and buriedness. It is NOT fpocket's trained
druggability score — it is an unsupervised geometric estimate, and labelled as
such in the output.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
from scipy import ndimage
from tqdm import tqdm

from .structio import load_atoms

# 7 scan directions: 3 axes + 4 body diagonals (classic LIGSITE set).
_DIRECTIONS = [
    (1, 0, 0), (0, 1, 0), (0, 0, 1),
    (1, 1, 1), (1, 1, -1), (1, -1, 1), (-1, 1, 1),
]


def _rasterise(coords: np.ndarray, radii: np.ndarray, spacing: float, margin: float):
    """Boolean protein-occupancy grid + the grid origin."""
    origin = coords.min(axis=0) - margin
    extent = coords.max(axis=0) + margin
    dims = np.ceil((extent - origin) / spacing).astype(int) + 1
    grid = np.zeros(dims, dtype=bool)

    # Mark voxels within each atom's vdW radius.
    idx = np.floor((coords - origin) / spacing).astype(int)
    max_r = int(np.ceil(radii.max() / spacing))
    offs = range(-max_r, max_r + 1)
    ball = [
        (dx, dy, dz)
        for dx in offs for dy in offs for dz in offs
    ]
    ball = np.array(ball)
    ball_dist = np.linalg.norm(ball * spacing, axis=1)

    for (i, j, k), r in zip(idx, radii):
        near = ball[ball_dist <= r]
        vox = near + np.array([i, j, k])
        m = np.all((vox >= 0) & (vox < dims), axis=1)
        vox = vox[m]
        grid[vox[:, 0], vox[:, 1], vox[:, 2]] = True
    return grid, origin


def _protein_ahead(protein: np.ndarray, direction: tuple[int, int, int]) -> np.ndarray:
    """Boolean grid: is there any protein voxel strictly ahead along `direction`?

    Reachability via the doubling technique (O(log diameter) grid shifts).
    """
    d = np.array(direction)
    ahead = np.roll(protein, tuple(d), axis=(0, 1, 2))
    ahead = _zero_wrap(ahead, d)
    step = d.copy()
    max_dim = max(protein.shape)
    reach = 1
    while reach < max_dim:
        shifted = np.roll(ahead, tuple(step), axis=(0, 1, 2))
        shifted = _zero_wrap(shifted, step)
        ahead = ahead | shifted
        step = step * 2
        reach *= 2
    return ahead


def _zero_wrap(arr: np.ndarray, shift: np.ndarray) -> np.ndarray:
    """Blank out the wrap-around region introduced by np.roll for a given shift."""
    out = arr
    for axis, s in enumerate(shift):
        if s == 0:
            continue
        sl = [slice(None)] * 3
        if s > 0:
            sl[axis] = slice(0, min(s, arr.shape[axis]))
        else:
            sl[axis] = slice(max(arr.shape[axis] + s, 0), arr.shape[axis])
        idx = tuple(sl)
        out = out.copy()
        out[idx] = False
    return out


def detect_pockets(
    coords: np.ndarray,
    radii: np.ndarray,
    spacing: float,
    margin: float,
    min_directions: int,
    min_voxels: int,
) -> dict[str, Any]:
    """Run the LIGSITE-style scan and summarise the pockets found."""
    if len(coords) < 20:
        return {"n_pockets": 0, "top_druggability": None, "method": "geometry"}

    protein, _ = _rasterise(coords, radii, spacing, margin)
    solvent = ~protein

    psp = np.zeros(protein.shape, dtype=np.uint8)
    for d in _DIRECTIONS:
        fwd = _protein_ahead(protein, d)
        bwd = _protein_ahead(protein, tuple(-np.array(d)))
        psp += (fwd & bwd).astype(np.uint8)

    pocket_vox = solvent & (psp >= min_directions)
    labels, n = ndimage.label(pocket_vox)
    if n == 0:
        return {"n_pockets": 0, "top_druggability": None, "method": "geometry"}

    sizes = ndimage.sum(np.ones_like(labels), labels, index=range(1, n + 1))
    sizes = np.asarray(sizes)
    keep = sizes >= min_voxels
    n_pockets = int(keep.sum())
    if n_pockets == 0:
        return {"n_pockets": 0, "top_druggability": None, "method": "geometry"}

    top_label = int(np.argmax(sizes)) + 1
    top_size = float(sizes[top_label - 1])
    top_buried = float(psp[labels == top_label].mean()) / len(_DIRECTIONS)  # 0-1

    volume = top_size * (spacing ** 3)  # Angstrom^3
    # Squash volume (~50-800 A^3 typical) and combine with buriedness -> 0-1.
    vol_term = 1.0 / (1.0 + np.exp(-(volume - 250.0) / 120.0))
    cavity_score = round(float(0.6 * vol_term + 0.4 * top_buried), 3)

    return {
        "n_pockets": n_pockets,
        "top_druggability": cavity_score,  # geometric proxy (see module docstring)
        "top_volume_A3": round(volume, 1),
        "method": "geometry",
    }


def run_pockets_geometry(
    config: dict[str, Any], accessions: list[str]
) -> dict[str, dict[str, Any]]:
    cfg = config["pocket"]
    spacing = float(cfg.get("grid_spacing", 1.2))
    margin = float(cfg.get("grid_margin", 5.0))
    min_dirs = int(cfg.get("min_directions", 5))
    min_vox = int(cfg.get("min_voxels", 30))

    structures_dir = Path(config["paths"]["structures_dir"])
    out: dict[str, dict[str, Any]] = {}
    for acc in tqdm(accessions, desc="pockets", unit="struct"):
        cif = structures_dir / f"AF-{acc}.cif"
        if not cif.exists():
            continue
        coords, radii = load_atoms(cif)
        try:
            out[acc] = detect_pockets(
                coords, radii, spacing, margin, min_dirs, min_vox
            )
        except Exception as exc:  # pragma: no cover
            print(f"[pocket] {acc}: skipped ({exc})")
            out[acc] = {"n_pockets": 0, "top_druggability": None, "method": "geometry"}

    print(f"[pocket] analysed {len(out)} candidates (geometry backend)")
    return out
