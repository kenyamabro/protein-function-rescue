# Protein Function Rescue

**Structure-based rediscovery of "hypothetical proteins" using the AlphaFold Protein Structure Database.**

Every sequenced genome contains proteins annotated only as *"uncharacterized"*, *"hypothetical"*,
or *"domain of unknown function (DUF)"*. Their amino-acid sequence matched nothing known, so no
biological function was ever assigned. But **structure is conserved long after sequence similarity
fades** — a protein can be a near-perfect 3D match to a well-studied enzyme while looking like noise
at the sequence level.

This project mines a whole organism's [AlphaFold DB](https://alphafold.ebi.ac.uk/) proteome, finds
unannotated proteins whose *shape* strongly matches a protein of **known** function, characterises
their binding pockets, and produces a ranked list of **candidate functional assignments** — plus a
browsable web viewer to explore them in 3D.

> The output is genuinely new data: previously-unannotated proteins with newly-proposed functions,
> ranked by confidence. That is a real, citable contribution to structural biology.

---

## Why this is scientifically sound

- **Consumes** pre-computed AlphaFold predictions — no GPU and no structure-prediction model needed.
  Runs on an ordinary laptop.
- Every prediction carries a **pLDDT confidence** (per-residue) that we filter on, so we never trust
  a low-confidence model.
- Function is transferred only from **experimentally/functionally characterised** targets, and always
  reported as a *hypothesis* with its supporting evidence.

## The pipeline

```
download  →  parse       →  annotate     →  search          →  pocket          →  rank
(AFDB)       (pLDDT, seq)    (UniProt name)   (structural hit)   (binding cavity)   (composite score)
```

1. **download** – fetch the target proteome's AlphaFold structures (one `.tar` from the AFDB FTP).
2. **parse** – read each model, compute mean pLDDT, extract sequence/length.
3. **annotate** – pull each protein's UniProt name/keywords; flag the *unannotated* ones.
4. **search** – structural search of the unannotated proteins against a library of **known** folds.
5. **pocket** – detect and score each hit's binding cavity.
6. **rank** – combine structural match + model confidence + pocket quality into one score, and
   emit `results/results.json` and `web/data.js`.

### Two interchangeable backends (set in `config.yaml`)

The **search** and **pocket** stages each have two backends, so the whole thing runs **natively on
Windows with no external binaries** — while keeping the faster Linux tools available:

| Stage  | `native` backend (default)                          | Linux/macOS backend |
|--------|-----------------------------------------------------|---------------------|
| search | **TM-align** (`tmtools`) vs. an auto-built reference library of known enzymes | **Foldseek** vs. PDB/Swiss-Prot |
| pocket | **LIGSITE-style** geometric detector (NumPy/SciPy)  | **fpocket**         |

TM-align is the gold-standard structural aligner that *defines* the TM-score Foldseek reports, so the
native path is scientifically equivalent for a curated reference set; Foldseek's advantage is speed and
searching millions of targets. Switch with `structural_search.backend: foldseek` / `pocket.backend: fpocket`.

## Quick start (native backends — no external tools)

```bash
# 1. Install Python deps (includes the native TM-align + pocket backends)
python -m pip install -r requirements.txt

# 2. Run the whole pipeline on a small test slice first
python -m pipeline.run all --limit 50

# 3. Then run it for real
python -m pipeline.run all
```

That's it — the default config uses the TM-align and geometric-pocket backends, which need no external
binaries. The first `search` run auto-downloads a small reference library of known-function enzymes into
`references/` (cached afterwards). The heavy analysis runs **once, offline**; the result is a static site.

**Performance note.** TM-align is *exhaustive* (every query vs. every reference), so it's thorough but
not fast: budget roughly a minute per candidate against a 120-protein library. Keep runs quick by using
`--limit`, lowering `structural_search.reference_count`, or — for a whole proteome — switching to the
Foldseek backend. (This is exactly the speed/scale trade-off Foldseek was built to solve.)

### Optional: the faster Linux tools

To use Foldseek + fpocket instead (faster, and Foldseek scales to millions of targets), install them
(see [docs/SETUP_TOOLS.md](docs/SETUP_TOOLS.md)), then set `structural_search.backend: foldseek` and
`pocket.backend: fpocket` in `config.yaml`.

## Viewing the results

Open `web/index.html` in a browser (double-click works). It reads `web/data.js`, shows a sortable,
searchable table of candidates, and renders each structure in 3D on click. Deploy it for free on
GitHub Pages by serving the `web/` folder.

A **sample `web/data.js`** ships with the repo so the interface works before you run anything.

## Running only part of the pipeline

Each stage caches to `results/`, so you can re-run stages independently:

```bash
python -m pipeline.run download --limit 200
python -m pipeline.run parse
python -m pipeline.run annotate
python -m pipeline.run search        # structural search (tmalign or foldseek backend)
python -m pipeline.run pocket        # optional; --skip-pocket to omit
python -m pipeline.run rank
```

## Choosing a different organism

Edit `config.yaml` → `organism`. Point `proteome_tar_url` at any proteome listed on the
[AlphaFold DB FTP](https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/). Understudied pathogens and
extremophiles have the most unannotated proteins — i.e. the most potential findings.

## Limitations (read before claiming anything)

- Predicted structures are **models**, not experiments. Low-pLDDT regions are unreliable.
- A structural match suggests a *possible* function; it is a hypothesis for wet-lab follow-up, not
  proof.
- Foldseek TM-scores and fpocket druggability are heuristics — treat the ranking as a triage tool.

## License

MIT — see [LICENSE](LICENSE). AlphaFold data is under the
[EMBL-EBI terms](https://alphafold.ebi.ac.uk/assets/License-Disclaimer.pdf) and is for theoretical
modelling only (no clinical use).
