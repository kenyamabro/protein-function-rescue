# Protein Function Rescue

**Structure-based rediscovery of "hypothetical proteins" using the AlphaFold Protein Structure Database.**

Every sequenced genome contains proteins annotated only as *"uncharacterized"*, *"hypothetical"*,
or *"domain of unknown function (DUF)"*. Their amino-acid sequence matched nothing known, so no
biological function was ever assigned. But **structure is often conserved after sequence similarity
fades**: a protein can share a fold or domain with a well-studied enzyme while showing little pairwise
sequence identity.

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

TM-align provides exact global pairwise alignments for the curated reference set. Foldseek instead uses
its 3Di representation, prefilters, and local alignment machinery to search very large databases; the
backends are complementary but not methodologically equivalent. Switch with
`structural_search.backend: foldseek` / `pocket.backend: fpocket`.

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
not fast, and runtime grows sharply with protein length. In the pinned pilot, nine queries against 119
references took about 19 minutes on one tested laptop; a 1,327-residue query dominated the run. Keep runs quick by using
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

## Reproducing the paper's pilot

The manuscript's pilot (Section 3) uses a fixed set of 10 *M. tuberculosis* uncharacterized proteins,
pinned in [`results/pilot_accessions.txt`](results/pilot_accessions.txt) (9 pass the pLDDT filter →
4 candidates). To reproduce it exactly, keep only those structures in `data/structures/` and run
`parse → annotate → search → pocket → rank` — the full `pipeline.run download` fetches the entire
proteome instead. The canonical 4-candidate output is tracked in `results/sample_results.json`, and
`results/results.json` is kept in sync with it.

## Manuscript, figures, and slides

The paper lives in [`paper/`](paper/): `protein_function_rescue.md`, the rendered `.docx`/`.pdf`, and a
`.pptx` deck. Regenerate everything with:

```bash
python -m pip install -r paper/requirements-paper.txt   # pandoc, docx2pdf, playwright
python -m playwright install chromium                   # one-time, for the Figure 2 screenshot

python paper/make_figure1.py     # pipeline schematic (Figure 1)
python paper/make_figure2.py     # screenshot of the live 3D web viewer (Figure 2)
python paper/render.py           # Markdown -> DOCX (+ PDF if MS Word / a LaTeX engine is present)
```

**Figure 2 — the viewer screenshot.** `paper/make_figure2.py` generates a clean `web/data.js` from
`results/sample_results.json`, serves `web/` on `127.0.0.1:8137`, opens it in headless Chromium with
WebGL enabled (so 3Dmol.js renders the structure), waits for the AlphaFold model to load, and writes
`paper/figures/fig2_viewer.png`. If the 3D panel comes out blank, increase the `time.sleep(...)` in that
script or check your network (the viewer fetches the structure from AlphaFold DB). Benchmark figures
(3–6) are written by `benchmark/run_benchmark.py` and `benchmark/weight_sensitivity.py`. The slide deck
is built with `NODE_PATH=$(npm root -g) node paper/make_slides.js` (after `npm install -g pptxgenjs`).

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
- TM-scores, the native geometric cavity score, and fpocket druggability (when that backend is used)
  are different heuristics; treat the ranking as a triage tool.

## Provenance and AI assistance

This project was built with assistance from Anthropic Claude (initial design, implementation,
benchmark, figures and drafting) and OpenAI Codex (independent verification, corrections, testing,
paired statistical analysis and manuscript revision). The human author reviews, verifies, and takes
full responsibility for the code, analyses and claims; neither AI system is an author. See the
generative-AI disclosure in the manuscript Methods and cite via [CITATION.cff](CITATION.cff).

## Citing

If you use this software or its results, please cite it via [CITATION.cff](CITATION.cff) (GitHub shows a
"Cite this repository" button). Connecting the repository to [Zenodo](https://zenodo.org) and cutting a
release mints a citable DOI automatically from that file.

## License

MIT — see [LICENSE](LICENSE). AlphaFold data is under the
[EMBL-EBI terms](https://alphafold.ebi.ac.uk/assets/License-Disclaimer.pdf) and is for theoretical
modelling only (no clinical use).
