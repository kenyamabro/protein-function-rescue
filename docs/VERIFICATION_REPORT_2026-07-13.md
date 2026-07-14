# Independent verification report — 2026-07-13

> **Remediation status (2026-07-13):** The code, results, web viewer, README, and manuscript were revised in response to this audit. The report below records the pre-remediation findings; it is retained as an audit trail rather than a description of the current manuscript.

Repository commit: `bf68b9ff45d5b07df87cdd62e7ee5c503b0af81d`

Environment used: Python 3.12.13 with the packages declared in `requirements.txt`. Original result files were preserved under `verification-snapshots/20260713-182212/` before rerunning the pipeline.

## Executive verdict

The software pipeline and its reported pilot/benchmark summary numbers are reproducible from the checked-in inputs and caches. The manuscript's biological interpretation is substantially stronger than its evidence supports, however. O08343 supports a broad TatD/metallo-hydrolase-fold interpretation, but the pipeline does not establish an ochratoxinase/amidohydrolase-2 function. P9WIT1 and P9WQ67 lack identity or even conservative chemistry at every checked reference functional residue, so reaction-specific function transfer is not defensible. P9WJG7 is a short-fragment normalization artefact and should not be presented as a candidate functional assignment.

The manuscript should not be submitted or signed in its present form. The numerical results may remain after the caveats below are incorporated, but the biological and novelty claims require revision.

## 1. Citation audit

All 14 cited works exist. Author, journal, year, volume, and the manuscript's main page/article-number details were confirmed against publisher, PubMed, or journal records, with these qualifications:

1. Perdigão et al. is correct: *PNAS* 112(52), 15898–15903 (2015), DOI `10.1073/pnas.1508380112`.
2. Illergård et al. is correct: *Proteins* 77(3), 499–508 (2009), DOI `10.1002/prot.22458`.
3. Jumper et al. is correct: *Nature* 596, 583–589 (2021), DOI `10.1038/s41586-021-03819-2`.
4. Varadi et al. is correct as the 2022 issue citation: *NAR* 50(D1), D439–D444, DOI `10.1093/nar/gkab1061`. It was published online in November 2021; 2022 is the issue year.
5. van Kempen et al. is correct as the 2024 issue citation: *Nature Biotechnology* 42, 243–246, DOI `10.1038/s41587-023-01773-0`. It was published online in 2023.
6. Zhang & Skolnick TM-align is correct: *NAR* 33(7), 2302–2309 (2005), DOI `10.1093/nar/gki524`.
7. Zhang & Skolnick TM-score is correct: *Proteins* 57(4), 702–710 (2004), DOI `10.1002/prot.20264`. A 2007 erratum exists (DOI `10.1002/prot.21643`).
8. Wojdyr is correct: *JOSS* 7(73), article 4200 (2022), DOI `10.21105/joss.04200`.
9. UniProt Consortium is correct: *NAR* 51(D1), D523–D531 (2023), DOI `10.1093/nar/gkac1052`.
10. Hendlich et al. needs a page correction: PubMed indexes *J Mol Graph Model* 15(6), **359–363, 389** (1997), DOI `10.1016/S1093-3263(98)00002-3`. The manuscript omits page 389.
11. Le Guilloux et al. is correct: *BMC Bioinformatics* 10, article 168 (2009), DOI `10.1186/1471-2105-10-168`.
12. Gligorijević et al. is correct: *Nature Communications* 12, article 3168 (2021), DOI `10.1038/s41467-021-23303-9`.
13. Thornton, Laskowski & Borkakoti is correct: *Nature Medicine* 27, 1666–1669 (2021), DOI `10.1038/s41591-021-01533-0`.
14. Akdel et al. is correct: *Nature Structural & Molecular Biology* 29, 1056–1067 (2022), DOI `10.1038/s41594-022-00849-w`.

Conclusion: no reference needs to be dropped, but LIGSITE's pagination should be corrected and the online-publication/issue-year distinction for Foldseek and AlphaFold DB should be understood.

## 2. Fresh pilot reproduction

The 10 files in `data/structures/` matched `results/pilot_accessions.txt` exactly.

Fresh commands executed:

```powershell
python -m pipeline.run parse
python -m pipeline.run annotate
python -m pipeline.run search
python -m pipeline.run pocket
python -m pipeline.run rank
```

Observed preprocessing:

- 10 pinned structures parsed.
- 9 passed mean pLDDT ≥ 70; 1 was dropped.
- Current UniProt annotation logic flagged all 9 retained proteins as unannotated.
- Exhaustive search against the cached 119-reference library produced 4 hits at TM-score ≥ 0.5.
- Search runtime was 18 min 37 s on this environment. P9WN15 (1,327 residues) dominated runtime. The README's rough one-minute-per-candidate estimate does not hold for long proteins.

Fresh output exactly reproduced the canonical ranking:

| Rank | Candidate | Reference | TM-score | Composite |
|---:|---|---|---:|---:|
| 1 | P9WIT1 | A0A0H3KZS3 | 0.9042 | 0.9329 |
| 2 | P9WJG7 | A0A009IHW8 | 0.8618 | 0.8440 |
| 3 | O08343 | A0A075TJ05 | 0.7154 | 0.8423 |
| 4 | P9WQ67 | A0A0D3LQY6 | 0.6663 | 0.8164 |

No tracked result diff was produced by the fresh run.

## 3. Benchmark verification

The benchmark driver was executed. It loaded the checked-in/cached 227-protein dataset, TM matrix, and sequence-identity matrix and regenerated `metrics.json`, `per_query.csv`, and the figures.

The summary was then recomputed independently from all 227 CSV rows:

- 227 proteins, 29 single-labelled EC families.
- Structure accuracy at EC levels 1–4: 0.965 at every level.
- Sequence accuracy: 0.965, 0.943, 0.943, 0.938 at levels 1–4.
- Correct structure calls at EC level ≥3: 219.
- Correct structure calls below 30% identity: 24/219 = 0.110.
- At TM-score ≥0.5: 223/227 accepted; precision 0.982; coverage 0.982.

These values exactly match `metrics.json`.

### Matrix-integrity spot check

Twenty-three pairs were recomputed directly, including every structural-error pair plus 15 deterministic random pairs. TM-align values matched the cached matrix to floating-point precision.

The full sequence matrix was also recomputed in reverse argument order. This exposed an implementation defect:

- `_identity(aligner, a, b)` is not always equal to `_identity(aligner, b, a)` for local alignments.
- The code calculates only the `i,j` direction and mirrors it as though it were symmetric.
- 1,400 of 25,651 pairs differed when argument order was reversed.
- Maximum absolute difference was 14.894 percentage points; mean absolute difference was 0.091 points.
- One nearest-sequence hit changed.
- The reported EC-level sequence accuracies and 24/219 twilight count did **not** change.

Example: P0C037→P84140 gives 17.02% in one order, while P84140→P0C037 gives 12.77% in the other.

Required fix: define a genuinely symmetric sequence score (for example a documented mean/max/min of both directions, or identity over an explicitly reconstructed alignment with a symmetric normalization), regenerate the matrix, and add a symmetry test. The current headlines happen to be stable, but the implementation is not methodologically clean.

### Benchmark interpretation caveats

- “Structure matched sequence overall” is imprecise: exact-EC accuracy is 96.5% versus 93.8%. Say that both were high and the observed 2.7-point difference was not subjected to a significance/paired uncertainty analysis.
- The “twilight-zone” result counts low identity to the **structural hit**, not cases in which the sequence baseline necessarily failed. It should not automatically be described as 24 unique wins over sequence.
- The dataset assigns each entry the single EC family under which it was fetched. Multi-functional/multi-EC enzymes can therefore be labelled incorrectly for valid secondary functions.
- The current run reused cached all-pairs matrices. Their integrity was spot-checked, not recomputed in full from zero during this audit.
- The baseline is local pairwise BLOSUM62 alignment, not BLAST, profile HMM, or modern protein embeddings. Claims should say “simple pairwise-sequence baseline,” never “sequence methods generally.”

## 4. Functional-residue and numbering audit

All eight AlphaFold models in the four candidate/reference pairs were checked. Each contains one chain A whose Cα-bearing residues are numbered uniquely and contiguously from 1 through L. The current `residue number - 1` assumption therefore holds for these four pairs.

It remains unsafe as a general implementation assumption because `load_ca()` discards author residue IDs, chain IDs, gaps, and insertion codes.

Fresh residue results:

| Candidate | Reference | Positional ≤4 Å | Identity | Chemistry |
|---|---|---:|---:|---:|
| P9WIT1 | A0A0H3KZS3 | 5/6 | 0/6 | 0/6 |
| P9WJG7 | A0A009IHW8 | 0/6 | 0/6 | 0/6 |
| O08343 | A0A075TJ05 | 8/8 | 4/8 | 4/8 |
| P9WQ67 | A0A0D3LQY6 | 11/11 | 0/11 | 0/11 |

The nearest-Cα method measures spatial occupancy, not an alignment correspondence. A random or chemically unrelated candidate residue can count as “position conserved.” Accordingly, the 5/6 and 11/11 positional figures are not positive catalytic-residue evidence when identity and conservative chemistry are both 0%.

## 5. TM-score normalization audit

The search ranks by the larger of the two length-normalized TM-scores:

| Candidate → reference | Lengths | Query-normalized | Reference-normalized | RMSD Å |
|---|---:|---:|---:|---:|
| P9WIT1 → A0A0H3KZS3 | 459 → 1022 | 0.9042 | 0.4213 | 2.647 |
| P9WJG7 → A0A009IHW8 | 50 → 269 | 0.8618 | 0.1801 | 1.138 |
| O08343 → A0A075TJ05 | 264 → 480 | 0.7154 | 0.4228 | 3.461 |
| P9WQ67 → A0A0D3LQY6 | 398 → 423 | 0.6663 | 0.6325 | 4.014 |

Consequences:

- Only P9WQ67 exceeds 0.5 under both normalizations.
- P9WJG7 is plainly a short-fragment/domain match, not a same-full-fold match.
- P9WIT1 aligns strongly as a 459-residue query to part of a 1,022-residue reference; describing it as a “near-perfect 3D match” without coverage/normalization is misleading.
- O08343 likewise does not cover the full 480-residue reference.

At minimum, both TM-scores, both lengths, aligned coverage, and RMSD should be reported. A hit threshold based only on `max()` is too permissive for automatic full-function transfer.

## 6. Biological claim audit

### O08343 (TatD)

Current UniProt independently annotates O08343 as a reviewed, 264-aa “Uncharacterized metal-dependent hydrolase TatD,” EC 3.1.-.-, belonging to the TatD-type metallo-dependent hydrolase family. UniProt already annotates candidate metal-binding residues 5, 7, 93, 134, 158, and 208 by similarity.

The matched reference A0A075TJ05 is a 480-aa ochratoxinase, EC 3.4.17.-. Its reference H111/H113/D378 map geometrically to O08343 H5/H7/D208. Thus:

- The original residue wording is wrong if it implies O08343 contains H111/H113/D378.
- H5/H7/D208 are already present in the candidate's UniProt annotation; the pipeline has not newly discovered them.
- The result supports a shared metallo-hydrolase architecture and metal-site compatibility.
- It does **not** establish ochratoxinase, carboxypeptidase, or “Zn-amidohydrolase” substrate specificity.

Safe wording: “O08343 recovered a metallo-hydrolase-family structural match consistent with its existing TatD annotation; four reference-site residue identities, including candidate H5/H7/D208, support metal-site compatibility but not transfer of the reference's specific reaction.”

### P9WIT1 (Rv2280)

Current UniProt calls P9WIT1 a reviewed 459-aa “Uncharacterized FAD-linked oxidoreductase Rv2280,” EC 1.-.-.-. The 1,022-aa reference is D-2-hydroxyglutarate dehydrogenase, EC 1.1.99.39.

Broad oxidoreductase-class compatibility is biologically plausible. A D2HGDH-specific assignment is not supported because:

- Reference-normalized TM-score is only 0.4213.
- All six checked functional sites have 0% identity and 0% conservative chemistry.
- Several reported nearest-Cα matches are chemically incompatible (for example Fe–S cysteines mapping to Gly/Glu/Val).

Safe wording: “The match is consistent with an oxidoreductase-related fold but does not refine P9WIT1 to D-2-hydroxyglutarate dehydrogenase.”

### P9WQ67 (Rv3778c)

The 398-aa candidate is currently uncharacterized. The 423-aa reference is a PLP-dependent aromatic amino-acid aminotransferase with EC 2.6.1.57, 2.6.1.58, and 2.6.1.70. Both TM normalizations exceed 0.5, making this the cleanest whole-chain fold match of the four.

Nevertheless, all 11 reference feature positions map to non-identical and non-conservative residues under the current proxy. The output does not verify the PLP catalytic machinery. The manuscript may call this an aminotransferase-fold hypothesis, but not a demonstrated PLP aminotransferase function or any specific EC reaction without a correct sequence/structure alignment of catalytic residues.

### P9WJG7 (ArfB)

The candidate is only 50 aa; the reference is 269 aa. TM-score is 0.8618 when normalized by the 50-aa query but 0.1801 when normalized by the reference. No pocket is detected and all six reference functional sites are 36–55 Å from the transformed candidate.

Verdict: spurious short-fragment match. Remove it from the biological candidate table or mark it explicitly as a negative-control/known failure of max-normalization. It must not be counted as a credible rescued function.

## 7. Pocket and ranking audit

The native pocket score is an unsupervised geometric proxy, not fpocket druggability. The logistic volume term and buriedness make the pilot scores saturate at 0.966–0.968 for the three proteins with pockets, so pocket quality barely distinguishes them.

The weight-sensitivity JSON is internally consistent: within the declared reasonable neighborhood, median Kendall τ is 0.9769 (minimum 0.9487); over the full simplex, median is 0.9128. This demonstrates stability for the constructed 40-item candidate pool, not general robustness across organisms or reference libraries.

The ranking labels and manuscript should consistently say “geometric cavity score,” not “druggability,” when the geometry backend is used.

## 8. Framing and novelty

The statement that published workflows “typically depend on HPC/GPUs/large pre-indexed databases” is too broad and is not established by references 5, 12–14 alone. Foldseek is specifically a fast CPU-oriented search tool and exposes a webserver; DeepFRI also exposes a webserver. The defensible novelty is narrower:

“This implementation combines precomputed AlphaFold models, exhaustive TM-align against a small curated reference library, a transparent geometric cavity heuristic, and a static viewer in a Windows-compatible teaching workflow without requiring a local GPU or compiled command-line bioinformatics executable.”

Avoid claims that TM-align is scientifically equivalent to Foldseek. Foldseek uses a 3Di structural alphabet, local alignments, sequence information, and prefilters; its paper reports different sensitivity/precision behavior from TM-align.

## 9. Required manuscript changes before sign-off

1. Correct the LIGSITE page range.
2. Replace the O08343 residue numbering and narrow the claim to broad metallo-hydrolase compatibility.
3. Remove the D2HGDH-specific conclusion for P9WIT1.
4. Recast P9WQ67 as an aminotransferase-fold hypothesis pending catalytic-site verification.
5. Remove or explicitly flag P9WJG7 as a normalization-driven false positive.
6. Report both TM normalizations, lengths, and coverage for every pilot hit.
7. Fix and regenerate the symmetric sequence baseline; add a regression test.
8. Describe 96.5% versus 93.8% without claiming statistical superiority or equivalence unless a paired uncertainty/significance analysis is added.
9. Clarify that 24/219 are correct structural calls below 30% identity to their structural hits, not necessarily 24 cases where structure beat sequence.
10. Replace “druggability” with “geometric cavity score” for native-backend results.
11. Narrow the novelty claim and remove “scientifically equivalent to Foldseek.”
12. Update the runtime estimate to warn about strong protein-length dependence.

Once those changes and the sequence-score fix are made and rechecked, the AI-disclosure sentence can truthfully state that citations, analyses, and claims were reviewed and verified.
