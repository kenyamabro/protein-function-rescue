"""Pipeline orchestrator / CLI.

Each stage caches its output to the results directory, so stages can be run
independently and re-used:

    python -m pipeline.run all --limit 50      # quick end-to-end test
    python -m pipeline.run all                 # full run
    python -m pipeline.run foldseek            # just re-run one stage
    python -m pipeline.run all --skip-pocket   # omit the fpocket stage

Cache files (in results/):
    proteins.json        parsed + annotated proteins
    foldseek_hits.json   best structural hit per protein
    pockets.json         fpocket stats per candidate
    results.json         final ranked candidates  (also web/data.js, candidates.csv)
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

if __package__ in {None, ""}:
    import os

    package_root = Path(__file__).resolve().parent.parent
    if str(package_root) not in sys.path:
        sys.path.insert(0, str(package_root))
    from pipeline import annotate as annotate_mod
    from pipeline import download as download_mod
    from pipeline import foldseek as foldseek_mod
    from pipeline import parse as parse_mod
    from pipeline import pocket as pocket_mod
    from pipeline import pocket_geom as pocket_geom_mod
    from pipeline import rank as rank_mod
    from pipeline import tmalign as tmalign_mod
    from pipeline.config import ensure_dirs, load_config
else:
    from . import annotate as annotate_mod
    from . import download as download_mod
    from . import foldseek as foldseek_mod
    from . import parse as parse_mod
    from . import pocket as pocket_mod
    from . import pocket_geom as pocket_geom_mod
    from . import rank as rank_mod
    from . import tmalign as tmalign_mod
    from .config import ensure_dirs, load_config

# "search" is the structural-search stage (tmalign or foldseek backend).
# "foldseek" is kept as an alias for backward compatibility.
STAGES = ["download", "parse", "annotate", "search", "foldseek", "pocket", "rank", "all"]


def _cache(config: dict[str, Any], name: str) -> Path:
    return Path(config["paths"]["results_dir"]) / name


def _load_cache(config: dict[str, Any], name: str) -> Any:
    path = _cache(config, name)
    if not path.exists():
        sys.exit(
            f"error: {path} not found — run the earlier stage(s) first "
            f"(e.g. `python -m pipeline.run all`)."
        )
    return json.loads(path.read_text(encoding="utf-8"))


def _save_cache(config: dict[str, Any], name: str, data: Any) -> None:
    _cache(config, name).write_text(json.dumps(data, indent=2), encoding="utf-8")


def do_download(config: dict[str, Any], args: argparse.Namespace) -> None:
    download_mod.download_proteome(config, limit=args.limit)


def do_parse(config: dict[str, Any], args: argparse.Namespace) -> None:
    proteins = parse_mod.parse_all(config)
    _save_cache(config, "proteins.json", proteins)


def do_annotate(config: dict[str, Any], args: argparse.Namespace) -> None:
    proteins = _load_cache(config, "proteins.json")
    proteins = annotate_mod.annotate_proteins(proteins, config)
    _save_cache(config, "proteins.json", proteins)


def do_search(config: dict[str, Any], args: argparse.Namespace) -> None:
    backend = config["structural_search"].get("backend", "tmalign").lower()
    if backend == "foldseek":
        print("[search] backend: foldseek")
        hits = foldseek_mod.run_foldseek(config)
    else:
        print("[search] backend: tmalign (native)")
        hits = tmalign_mod.run_tmalign(config)
    _save_cache(config, "foldseek_hits.json", hits)


def do_pocket(config: dict[str, Any], args: argparse.Namespace) -> None:
    if args.skip_pocket or not config["pocket"].get("enabled", True):
        print("[pocket] skipped")
        _save_cache(config, "pockets.json", {})
        return
    proteins = _load_cache(config, "proteins.json")
    hits = _load_cache(config, "foldseek_hits.json")
    # Only bother running fpocket on actual candidates (unannotated + has a hit).
    unannotated = {p["accession"] for p in proteins if p.get("originally_unannotated")}
    candidates = [acc for acc in hits.keys() if acc in unannotated]

    backend = config["pocket"].get("backend", "geometry").lower()
    if backend == "fpocket":
        print("[pocket] backend: fpocket")
        try:
            pockets = pocket_mod.run_pockets(config, candidates)
        except pocket_mod.PocketUnavailable as exc:
            print(f"[pocket] {exc}\n[pocket] continuing without pocket scores")
            pockets = {}
    else:
        print("[pocket] backend: geometry (native)")
        pockets = pocket_geom_mod.run_pockets_geometry(config, candidates)
    _save_cache(config, "pockets.json", pockets)


def do_rank(config: dict[str, Any], args: argparse.Namespace) -> None:
    proteins = _load_cache(config, "proteins.json")
    hits = _load_cache(config, "foldseek_hits.json")
    pockets_path = _cache(config, "pockets.json")
    pockets = (
        json.loads(pockets_path.read_text(encoding="utf-8"))
        if pockets_path.exists()
        else {}
    )
    candidates = rank_mod.build_candidates(proteins, hits, pockets, config)
    rank_mod.export(candidates, config)


def do_all(config: dict[str, Any], args: argparse.Namespace) -> None:
    do_download(config, args)
    do_parse(config, args)
    do_annotate(config, args)
    try:
        do_search(config, args)
    except foldseek_mod.FoldseekError as exc:
        sys.exit(
            f"\nerror: {exc}\n\nThe structural-search stage is required for results. "
            f"Either install Foldseek (docs/SETUP_TOOLS.md) or set "
            f"structural_search.backend: tmalign in config.yaml, then re-run "
            f"`python -m pipeline.run search` followed by `pocket` and `rank`."
        )
    do_pocket(config, args)
    do_rank(config, args)


DISPATCH = {
    "download": do_download,
    "parse": do_parse,
    "annotate": do_annotate,
    "search": do_search,
    "foldseek": do_search,  # backward-compatible alias
    "pocket": do_pocket,
    "rank": do_rank,
    "all": do_all,
}


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="Protein Function Rescue pipeline")
    parser.add_argument("stage", choices=STAGES, help="pipeline stage to run")
    parser.add_argument("--config", default="config.yaml", help="path to config.yaml")
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="only process the first N structures (for quick test runs)",
    )
    parser.add_argument(
        "--skip-pocket", action="store_true", help="omit the fpocket stage"
    )
    args = parser.parse_args(argv)

    # Make console output robust to non-ASCII regardless of the OS locale
    # (Windows consoles default to cp1252 and choke on some Unicode).
    for stream in (sys.stdout, sys.stderr):
        try:
            stream.reconfigure(encoding="utf-8", errors="replace")
        except (AttributeError, ValueError):
            pass

    config = load_config(args.config)
    ensure_dirs(config)
    try:
        DISPATCH[args.stage](config, args)
    except foldseek_mod.FoldseekError as exc:
        sys.exit(f"\nerror (foldseek): {exc}")
    except pocket_mod.PocketUnavailable as exc:
        sys.exit(f"\nerror (fpocket): {exc}")


if __name__ == "__main__":
    main()
