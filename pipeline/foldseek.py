"""Stage 4 — structural search with Foldseek.

Runs ``foldseek easy-search`` over the downloaded structures against a target
database of *functionally known* proteins (PDB or Swiss-Prot), and keeps the
best hit per query above the configured TM-score / e-value thresholds.

Build the target database once, e.g.:
    foldseek databases PDB foldseek_db/pdb tmp
"""

from __future__ import annotations

import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Any

from .config import accession_from_filename

# Foldseek tabular columns we request, in order.
_FORMAT = "query,target,evalue,alntmscore,theader"


class FoldseekError(RuntimeError):
    pass


def _resolve_binary(binary: str) -> str:
    found = shutil.which(binary)
    if found is None and not Path(binary).exists():
        raise FoldseekError(
            f"Foldseek binary '{binary}' not found. Install it and/or set "
            f"foldseek.binary in config.yaml. See docs/SETUP_TOOLS.md."
        )
    return found or binary


def run_foldseek(config: dict[str, Any]) -> dict[str, dict[str, Any]]:
    """Return {query_accession: {target, evalue, tm_score, description}}."""
    fs = config["foldseek"]
    binary = _resolve_binary(fs["binary"])
    target_db = Path(fs["target_db"])
    if not target_db.parent.exists():
        raise FoldseekError(
            f"Foldseek target DB '{target_db}' not found. Build it once with:\n"
            f"    {fs['binary']} databases PDB {target_db} tmp"
        )

    structures_dir = Path(config["paths"]["structures_dir"])
    results_dir = Path(config["paths"]["results_dir"])
    out_m8 = results_dir / "foldseek_raw.m8"

    with tempfile.TemporaryDirectory() as tmp:
        cmd = [
            binary,
            "easy-search",
            str(structures_dir),
            str(target_db),
            str(out_m8),
            tmp,
            "--format-output",
            _FORMAT,
        ]
        print(f"[foldseek] {' '.join(cmd)}")
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode != 0:
            raise FoldseekError(
                f"foldseek failed (exit {proc.returncode}):\n{proc.stderr[-2000:]}"
            )

    return _best_hits(out_m8, config)


def _best_hits(m8_path: Path, config: dict[str, Any]) -> dict[str, dict[str, Any]]:
    """Parse the m8 file and keep the best-scoring accepted hit per query."""
    fs = config["foldseek"]
    min_tm = float(config["structural_search"]["min_tmscore"])
    max_e = float(fs["max_evalue"])

    best: dict[str, dict[str, Any]] = {}
    if not m8_path.exists():
        return best

    with open(m8_path, "r", encoding="utf-8") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            query, target, evalue, tm = parts[0], parts[1], parts[2], parts[3]
            description = parts[4] if len(parts) > 4 else target

            acc = accession_from_filename(query) or query
            try:
                evalue_f = float(evalue)
                tm_f = float(tm)
            except ValueError:
                continue
            if tm_f < min_tm or evalue_f > max_e:
                continue

            prev = best.get(acc)
            if prev is None or tm_f > prev["tm_score"]:
                best[acc] = {
                    "target": target,
                    "evalue": evalue_f,
                    "tm_score": round(tm_f, 4),
                    "description": description.strip(),
                }

    print(f"[foldseek] {len(best)} queries with an accepted structural hit")
    return best
