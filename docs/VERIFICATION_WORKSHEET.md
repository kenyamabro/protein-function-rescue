# Verification worksheet: manuscript Tiers 1–2

Run every command from the repository root in PowerShell. Record the date, commit, Python version, and whether the working tree was clean before interpreting any result.

```powershell
Set-Location 'C:\Users\DELL\Documents\Projects\Programming\Apps\protein-function-rescue'
git rev-parse HEAD
git status --short
python --version
python -m pip freeze | Set-Content verification-pip-freeze.txt
```

Do not overwrite the tracked/canonical outputs without preserving a snapshot:

```powershell
$stamp = Get-Date -Format 'yyyyMMdd-HHmmss'
$snapshot = "verification-snapshots\$stamp"
New-Item -ItemType Directory -Force $snapshot | Out-Null
Copy-Item results\*.json,results\*.csv,benchmark\results\*.json,benchmark\results\*.csv $snapshot -ErrorAction SilentlyContinue
Get-FileHash results\sample_results.json,results\results.json,benchmark\results\metrics.json,benchmark\results\per_query.csv |
  Format-Table -AutoSize
```

## 1. References

For each reference, open the DOI, PubMed record, or publisher page and compare title, authors, journal, year, volume/issue, and pages/article number with the manuscript. Never use a search-result snippet as the final authority.

| # | Manuscript citation | Verification target | Status / correction |
|---|---|---|---|
| 1 | Perdigão et al., PNAS (2015) 112(52):15898–15903 | DOI and PubMed | [ ] |
| 2 | Illergård et al., Proteins (2009) 77(3):499–508 | DOI and PubMed | [ ] |
| 3 | Jumper et al., Nature (2021) 596:583–589 | Publisher | [ ] |
| 4 | Varadi et al., NAR (2022) 50(D1):D439–D444 | Publisher; distinguish original paper from later update | [ ] |
| 5 | van Kempen et al., Nature Biotechnology (2024) 42:243–246 | Publisher | [ ] |
| 6 | Zhang & Skolnick, NAR (2005) 33(7):2302–2309 | Publisher/PubMed | [ ] |
| 7 | Zhang & Skolnick, Proteins (2004) 57(4):702–710 | DOI/PubMed | [ ] |
| 8 | Wojdyr, JOSS (2022) 7(73):4200 | JOSS article page | [ ] |
| 9 | UniProt Consortium, NAR (2023) 51(D1):D523–D531 | Publisher | [ ] |
| 10 | Hendlich et al., J Mol Graph Model (1997) 15(6):359–363 | DOI/PubMed | [ ] |
| 11 | Le Guilloux et al., BMC Bioinformatics (2009) 10:168 | Publisher/PubMed | [ ] |
| 12 | Gligorijević et al., Nature Communications (2021) 12:3168 | Publisher | [ ] |
| 13 | Thornton et al., Nature Medicine (2021) 27:1666–1669 | Publisher | [ ] |
| 14 | Akdel et al., NSMB (2022) 29:1056–1067 | Publisher | [ ] |

Acceptance rule: all six bibliographic fields agree with a primary record. If a detail cannot be confirmed, correct it or remove the citation before signing the disclosure.

## 2. Reproduce the pinned pilot

The repository currently contains exactly 10 files in `data/structures`, matching the 10 accessions in `results/pilot_accessions.txt`. Confirm this before running:

```powershell
$wanted = Get-Content results\pilot_accessions.txt
$present = Get-ChildItem data\structures -Filter 'AF-*.cif' |
  ForEach-Object { $_.BaseName -replace '^AF-','' }
Compare-Object ($wanted | Sort-Object) ($present | Sort-Object)
```

Expected: no output. Any output means the input set is not the pinned pilot.

Re-run the stages:

```powershell
python -m pipeline.run parse
python -m pipeline.run annotate
python -m pipeline.run search
python -m pipeline.run pocket
python -m pipeline.run rank
```

Compare the fresh result with the canonical record:

```powershell
python -c "import json; p=json.load(open('results/results.json')); print(p['n_candidates']); [(print(c['accession'], c['best_structural_match']['target'], c['best_structural_match']['tm_score'], c['scores']['composite'])) for c in p['candidates']]"
git diff -- results/results.json results/proteins.json results/foldseek_hits.json results/pockets.json web/data.js
```

Expected canonical values:

| Rank | Candidate → reference | TM-score | Composite |
|---:|---|---:|---:|
| 1 | P9WIT1 → A0A0H3KZS3 | 0.9042 | 0.9329 |
| 2 | P9WJG7 → A0A009IHW8 | 0.8618 | 0.8440 |
| 3 | O08343 → A0A075TJ05 | 0.7154 | 0.8423 |
| 4 | P9WQ67 → A0A0D3LQY6 | 0.6663 | 0.8164 |

Also expect nine proteins to pass the mean-pLDDT filter and four hits at TM-score ≥ 0.5. Record any drift caused by current UniProt/AlphaFold data rather than silently replacing the canonical output.

## 3. Reproduce and audit the benchmark

The full benchmark is documented as a 2–3 hour resumable run:

```powershell
python benchmark/run_benchmark.py
```

Then independently compute the headline fields from the generated files:

```powershell
python -c "import json; m=json.load(open('benchmark/results/metrics.json')); print('n=',m['n_proteins'],'families=',m['n_families']); print('structure L4=',m['topN_accuracy_by_ec_level']['structure']['4']); print('sequence L4=',m['topN_accuracy_by_ec_level']['sequence']['4']); print('twilight=',m['twilight_zone']['n_correct_below_cutoff'],'/',m['twilight_zone']['n_correct_level3']); print('TM>=0.5=',next(x for x in m['precision_vs_tm_threshold'] if x['threshold']==0.5))"
python -c "import csv; r=list(csv.DictReader(open('benchmark/results/per_query.csv',newline=''))); print('rows=',len(r)); print('columns=',sorted(r[0]))"
git diff -- benchmark/results/metrics.json benchmark/results/per_query.csv benchmark/results/figures
```

Expected: 227 proteins, 29 families, exact four-level EC accuracy 0.965 (structure) and 0.938 (sequence), 24 of 219 correct level-3 calls below 30% identity, and at TM ≥ 0.5 precision 0.982 with coverage 0.982.

Audit at least five randomly selected rows and every apparent error in `per_query.csv` back to `dataset.json`. For multi-EC proteins, record both the benchmark family EC and every EC currently attached to the protein; classify a mismatch as “benchmark-wrong” only after ruling out a legitimate secondary function.

## 4. O08343 end-to-end residue check

Regenerate the residue report:

```powershell
python -m pipeline.functional_residues
python -c "import json; x=next(v for v in json.load(open('results/functional_residues.json')) if v['accession']=='O08343'); print(x['tm_score'],x['n_position_conserved'],x['n_identity_conserved']); [(print(s['ref_residue'],'->',s['candidate_residue'],s['distance_A'])) for s in x['sites']]"
```

Expected tracked output: TM-score 0.7154; 8/8 annotated reference positions within 4 Å; 4/8 residue identities conserved. Specifically, reference H111 → candidate H5 (1.45 Å), H113 → H7 (1.41 Å), and D378 → D208 (1.05 Å).

Critical wording check: do **not** write that “His111/His113/Asp378 are conserved in O08343.” Those are reference-protein numbers. State that the reference residues map structurally to O08343 H5/H7/D208, after confirming the numbering directly in both mmCIF files and in a visual superposition.

Manual inspection checklist:

- [ ] Load `data/structures/AF-O08343.cif` and the reference structure for A0A075TJ05 in PyMOL/ChimeraX.
- [ ] Reproduce the TM-align transform or load an exported superposition.
- [ ] Display reference H111, H113, D378 and candidate H5, H7, D208 as sticks.
- [ ] Confirm the reported Cα distances and inspect side-chain orientation/metal geometry.
- [ ] Check whether candidate H134 (reference H287) also supports a metal-binding site.
- [ ] Describe the check as a coarse nearest-Cα proxy unless side-chain geometry and coordination are separately demonstrated.

## 5. Other biological claims

### P9WIT1 / Rv2280

- [ ] Confirm the current UniProt annotation and evidence for “FAD-linked oxidoreductase” and EC 1.-.-.-.
- [ ] Confirm the reference A0A0H3KZS3 annotation and EC 1.1.99.39.
- [ ] Decide whether the match refines function or merely preserves a broad oxidoreductase fold.
- [ ] Note that the tracked functional-site report has 5/6 positional matches but **0/6 identity or chemistry matches**; this weakens a reaction-specific claim.

### P9WQ67 / Rv3778c

- [ ] Confirm the candidate is currently uncharacterized.
- [ ] Confirm the reference is a PLP-dependent aromatic amino-acid aminotransferase with multiple EC numbers.
- [ ] Inspect the PLP lysine and other catalytic residues, not only nearest-Cα occupancy.
- [ ] Note that the tracked report has 11/11 positional matches but **0/11 identity or chemistry matches**; fold evidence alone does not establish the specific reaction.

### P9WJG7 / ArfB

- [ ] Confirm length = 50 aa and lack of detected pocket.
- [ ] Compare both TM-score normalizations, aligned length, RMSD, and coverage.
- [ ] Treat the max-normalized TM-score of 0.8618 as potentially inflated by the short query.
- [ ] Note that the tracked residue report finds 0/6 reference functional sites within 4 Å (distances 36–55 Å). Label this hit spurious/self-flagged unless stronger independent evidence emerges.

## 6. Code assumptions to defend

### Residue numbering

`pipeline/functional_residues.py` uses `reference residue number - 1` as the coordinate/sequence-array index and reports candidate index `+ 1`. It does not read author residue IDs, insertion codes, chain breaks, or fragments when mapping functional sites.

- [ ] For all four candidate/reference pairs, inspect mmCIF residue IDs and confirm one chain numbered contiguously 1..L with no gaps or insertion codes.
- [ ] If that condition fails, fix the mapper before using Table 2.

### TM-score normalization

`pipeline/tmalign.py` ranks with `max(tm_norm_chain1, tm_norm_chain2)`. This asks whether either length normalization supports a fold match and favors short-to-long/domain matches. It is intentionally generous, but it can make a short fragment look strong.

- [ ] Record both normalizations for all four pilot hits.
- [ ] Report query/reference lengths and aligned coverage beside P9WJG7.
- [ ] Repeat the four-hit decision using `min(...)` or query-normalized TM-score as a sensitivity check.

## 7. One-sentence oral defenses

- **pLDDT:** “pLDDT is AlphaFold’s per-residue confidence score; we filtered models below mean 70 to avoid comparing globally unreliable predictions, while retaining local scores for inspection.”
- **TM-score:** “TM-score is a length-normalized structural-similarity measure; values above roughly 0.5 usually indicate the same overall fold, but normalization and coverage matter for short or partial matches.”
- **Pocket detector:** “The native detector marks solvent grid points enclosed by protein–solvent–protein transitions in several directions, clusters them into cavities, and uses a bounded heuristic score whose saturation limits discrimination.”
- **Composite score:** “The rank combines TM-score, scaled mean pLDDT, and pocket quality with weights 0.5/0.2/0.3, renormalizing when a component is missing; nearby weights barely change the ordering (median Kendall τ about 0.98).”
- **Residue check:** “After a TM-align superposition, each annotated reference site is matched to the nearest candidate Cα within 4 Å; this checks positional compatibility, not side-chain geometry or catalytic equivalence.”

## 8. Sign-off

- [ ] All 14 citations checked against primary records.
- [ ] Pilot re-run from exactly the pinned 10 structures; deviations explained.
- [ ] Benchmark re-run; `metrics.json` independently reconciled with `per_query.csv`.
- [ ] O08343 and one benchmark family inspected end to end.
- [ ] Four pilot biological interpretations checked against current literature and structures.
- [ ] Residue-numbering and TM-normalization assumptions tested.
- [ ] Multi-EC benchmark caveat documented.
- [ ] Novelty/positioning claims narrowed to what the cited literature supports.
- [ ] AI disclosure is literally accurate.

Only sign the disclosure after every applicable box is checked and the manuscript has been updated to reflect failures or caveats.
