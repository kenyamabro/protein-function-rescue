/* Deck: Protein Function Rescue — paper + program overview.
 * Build:  NODE_PATH=$(npm root -g) node paper/make_slides.js
 */
const pptxgen = require("pptxgenjs");
const P = new pptxgen();
P.layout = "LAYOUT_WIDE"; // 13.33 x 7.5 in
P.author = "Protein Function Rescue";

const ROOT = "C:/Users/DELL/Documents/Projects/Programming/Apps/protein-function-rescue";
const FIG = {
  pipeline: ROOT + "/paper/figures/fig1_pipeline.png",
  viewer: ROOT + "/paper/figures/fig2_viewer.png",
  acc: ROOT + "/benchmark/results/figures/fig_ec_accuracy.png",
  twilight: ROOT + "/benchmark/results/figures/fig_tm_vs_identity.png",
  weights: ROOT + "/benchmark/results/figures/fig_weight_sensitivity.png",
};

// Palette (matches the app: deep navy + teal)
const BG = "0F1420", PANEL = "1B2437", PANEL2 = "232E47";
const TEAL = "5EEAD4", BLUE = "7DD3FC", TEXT = "E6EBF5", MUTED = "9AA8C7", WHITE = "FFFFFF";
const HEAD = "Cambria", BODY = "Calibri", MONO = "Courier New";

const dark = () => { const s = P.addSlide(); s.background = { color: BG }; return s; };

function title(s, t, sub) {
  s.addText(t, { x: 0.6, y: 0.42, w: 12.1, h: 0.7, fontFace: HEAD, fontSize: 32,
    bold: true, color: WHITE, align: "left", margin: 0 });
  if (sub) s.addText(sub, { x: 0.62, y: 1.16, w: 12.1, h: 0.4, fontFace: BODY,
    fontSize: 15, color: TEAL, align: "left", margin: 0 });
}
// White-bg figure framed as an intentional card.
function figCard(s, path, x, y, w, h) {
  s.addShape("roundRect", { x: x - 0.12, y: y - 0.12, w: w + 0.24, h: h + 0.24,
    rectRadius: 0.08, fill: { color: WHITE }, line: { color: PANEL2, width: 1 },
    shadow: { type: "outer", color: "000000", opacity: 0.35, blur: 8, offset: 3, angle: 90 } });
  s.addImage({ path, x, y, w, h });
}
function chip(s, label, x, y, w, color) {
  s.addShape("roundRect", { x, y, w, h: 0.44, rectRadius: 0.09, fill: { color: PANEL2 },
    line: { color, width: 1 } });
  s.addText(label, { x, y, w, h: 0.44, fontFace: BODY, fontSize: 12.5, bold: true,
    color, align: "center", valign: "middle", margin: 0 });
}

/* 1 — Title */
(() => {
  const s = dark();
  s.addShape("roundRect", { x: 0.6, y: 2.0, w: 0.9, h: 0.9, rectRadius: 0.12,
    fill: { color: PANEL }, line: { color: TEAL, width: 1.5 } });
  s.addText("🧬", { x: 0.6, y: 2.0, w: 0.9, h: 0.9, fontSize: 34, align: "center", valign: "middle", margin: 0 });
  s.addText("Protein Function Rescue", { x: 1.7, y: 2.05, w: 11.0, h: 1.0, fontFace: HEAD,
    fontSize: 46, bold: true, color: WHITE, margin: 0 });
  s.addText("Structure-based rediscovery of hypothetical proteins from the AlphaFold Database",
    { x: 1.72, y: 3.15, w: 10.8, h: 0.6, fontFace: BODY, fontSize: 19, color: TEAL, margin: 0 });
  s.addText("A transparent, laptop-scale pipeline — no GPU, no HPC, fully reproducible",
    { x: 1.72, y: 3.75, w: 10.8, h: 0.5, fontFace: BODY, fontSize: 14, color: MUTED, margin: 0 });
  s.addText("[Student Name]  ·  Methods / proof-of-concept", { x: 1.72, y: 5.9, w: 10, h: 0.4,
    fontFace: BODY, fontSize: 13, color: MUTED, margin: 0 });
  s.addNotes("Protein Function Rescue: a laptop-scale, transparent pipeline that turns AlphaFold structures into ranked, evidence-linked functional hypotheses for proteins that sequence methods leave unannotated.");
})();

/* 2 — Problem */
(() => {
  const s = dark();
  title(s, "The dark proteome", "A quarter to a half of every proteome has no assigned function");
  const rows = [
    ["●", "\"Hypothetical\", \"uncharacterized\", DUF", "Sequence matched nothing known — no function assigned."],
    ["●", "Sequence homology fails first", "Below ~20–25% identity, function transfer by sequence breaks down."],
    ["●", "These proteins matter", "Enzymes, transporters, and untapped drug/vaccine targets in pathogens."],
  ];
  let y = 1.9;
  rows.forEach(r => {
    s.addText(r[0], { x: 0.7, y, w: 0.4, h: 0.5, color: TEAL, fontSize: 16, align: "center", margin: 0 });
    s.addText([{ text: r[1] + "  ", options: { bold: true, color: TEXT } },
               { text: r[2], options: { color: MUTED } }],
      { x: 1.15, y, w: 7.2, h: 0.7, fontFace: BODY, fontSize: 15, margin: 0, valign: "top" });
    y += 0.95;
  });
  // Stat cards
  const stat = (x, big, lab, col) => {
    s.addShape("roundRect", { x, y: 2.0, w: 2.05, h: 2.6, rectRadius: 0.1, fill: { color: PANEL }, line: { color: PANEL2, width: 1 } });
    s.addText(big, { x, y: 2.45, w: 2.05, h: 0.9, fontFace: HEAD, fontSize: 34, bold: true, color: col, align: "center", margin: 0 });
    s.addText(lab, { x: x + 0.12, y: 3.5, w: 1.81, h: 0.9, fontFace: BODY, fontSize: 12.5, color: MUTED, align: "center", margin: 0 });
  };
  stat(8.8, "200M+", "predicted structures in the AlphaFold DB", TEAL);
  stat(11.0, "25–50%", "of a typical proteome is unannotated", BLUE);
  s.addNotes("The gap between known sequences and known functions keeps widening. AlphaFold gives us 200M+ structures to attack it with.");
})();

/* 3 — Insight */
(() => {
  const s = dark();
  title(s, "Structure outlives sequence", "The biological basis for the whole approach");
  s.addText([
    { text: "Structure is empirically ", options: { color: TEXT } },
    { text: "3–10× more conserved", options: { color: TEAL, bold: true } },
    { text: " than sequence.\n\nSo a protein can be a near-perfect ", options: { color: TEXT } },
    { text: "3D match", options: { color: BLUE, bold: true } },
    { text: " to a well-studied enzyme while looking like ", options: { color: TEXT } },
    { text: "noise", options: { color: TEXT, italic: true } },
    { text: " at the sequence level. Comparing shape, not sequence, recovers those hidden relationships.", options: { color: TEXT } },
  ], { x: 0.7, y: 2.0, w: 7.0, h: 3.0, fontFace: BODY, fontSize: 18, lineSpacingMultiple: 1.15, margin: 0, valign: "top" });
  s.addShape("roundRect", { x: 8.3, y: 2.0, w: 4.4, h: 3.2, rectRadius: 0.12, fill: { color: PANEL }, line: { color: TEAL, width: 1.5 } });
  s.addText("3–10×", { x: 8.3, y: 2.5, w: 4.4, h: 1.2, fontFace: HEAD, fontSize: 60, bold: true, color: TEAL, align: "center", margin: 0 });
  s.addText("structure conserved longer than sequence\n(Illergård et al., 2009)", { x: 8.5, y: 3.9, w: 4.0, h: 1.0, fontFace: BODY, fontSize: 14, color: MUTED, align: "center", margin: 0 });
  s.addNotes("This is the entire biological motivation: fold outlasts sequence, so structural search reaches into the sequence 'twilight zone'.");
})();

/* 4 — Novelty */
(() => {
  const s = dark();
  title(s, "Why this pipeline exists", "The novelty is accessibility, transparency, and reproducibility");
  s.addText("Published structure-based annotation workflows typically depend on:", { x: 0.7, y: 1.85, w: 12, h: 0.4, fontFace: BODY, fontSize: 15, color: MUTED, margin: 0 });
  const cons = ["High-performance computing", "GPUs", "Large pre-indexed databases", "Complex software environments"];
  cons.forEach((c, i) => chip(s, c, 0.7 + i * 3.05, 2.35, 2.85, MUTED));
  s.addShape("roundRect", { x: 0.7, y: 3.15, w: 12.0, h: 1.5, rectRadius: 0.12, fill: { color: PANEL }, line: { color: TEAL, width: 2 } });
  s.addText([
    { text: "Our question:  ", options: { color: TEAL, bold: true } },
    { text: "Can a transparent, dependency-light workflow achieve meaningful annotation performance on commodity hardware while remaining fully reproducible?", options: { color: WHITE } },
  ], { x: 1.0, y: 3.25, w: 11.4, h: 1.3, fontFace: BODY, fontSize: 18, italic: true, valign: "middle", margin: 0 });
  chip(s, "Interpretable — unlike learned predictors (DeepFRI)", 0.7, 5.2, 5.9, BLUE);
  chip(s, "Accessible — unlike speed-first search (Foldseek)", 6.8, 5.2, 5.9, BLUE);
  s.addNotes("Existing work proves structure-based annotation works, but on heavy infrastructure. We ask whether it can be done transparently on a laptop — complementary to DeepFRI (interpretable) and Foldseek (accessible).");
})();

/* 5 — Pipeline */
(() => {
  const s = dark();
  title(s, "The pipeline", "Six cached stages — consumes AlphaFold DB, no GPU, no compiled binaries");
  figCard(s, FIG.pipeline, 1.15, 2.15, 11.0, 3.13);
  s.addText("download → parse → annotate → search → pocket → rank.  Each stage caches to disk and is independently re-runnable.",
    { x: 0.7, y: 5.7, w: 12, h: 0.5, fontFace: BODY, fontSize: 13.5, color: MUTED, align: "center", margin: 0 });
  s.addNotes("The pipeline consumes pre-computed AlphaFold models. The two compute-heavy stages — search and pocket — have native, dependency-free implementations.");
})();

/* 6 — Backends */
(() => {
  const s = dark();
  title(s, "Two interchangeable backends", "The same workflow runs on a laptop — or scales on a Linux server");
  const rows = [
    ["Stage", "Native (default, any OS)", "High-performance (Linux/macOS)"],
    ["search", "TM-align  (tmtools, pure Python)", "Foldseek"],
    ["pocket", "LIGSITE-style geometry  (NumPy/SciPy)", "fpocket"],
  ];
  const colX = [0.7, 3.0, 8.0], colW = [2.2, 4.9, 4.7];
  let y = 2.1;
  rows.forEach((r, ri) => {
    const head = ri === 0;
    r.forEach((cell, ci) => {
      s.addShape("roundRect", { x: colX[ci], y, w: colW[ci], h: 0.85, rectRadius: 0.06,
        fill: { color: head ? PANEL2 : PANEL }, line: { color: PANEL2, width: 1 } });
      s.addText(cell, { x: colX[ci] + 0.15, y, w: colW[ci] - 0.3, h: 0.85, fontFace: ci === 0 ? MONO : BODY,
        fontSize: head ? 14 : 15, bold: head || ci === 0, color: head ? TEAL : (ci === 1 ? WHITE : MUTED),
        valign: "middle", margin: 0 });
    });
    y += 0.98;
  });
  s.addText("Switch with one line in config.yaml. TM-align is the gold-standard aligner that defines the TM-score Foldseek reports.",
    { x: 0.7, y: 5.4, w: 12, h: 0.6, fontFace: BODY, fontSize: 14, color: MUTED, margin: 0 });
  s.addNotes("Native backends need no external binaries. Foldseek/fpocket are optional for speed and scale.");
})();

/* 7 — Pilot / viewer */
(() => {
  const s = dark();
  title(s, "Pilot: Mycobacterium tuberculosis", "Ranked candidates, each explorable in 3D");
  figCard(s, FIG.viewer, 0.7, 1.85, 8.2, 5.25); // 2560x1640 -> aspect 1.561
  s.addText([
    { text: "Top hit — P9WIT1\n", options: { bold: true, color: TEAL, fontSize: 17 } },
    { text: "\"Uncharacterized FAD-linked oxidoreductase\"\n\n", options: { color: TEXT, fontSize: 14 } },
    { text: "matches ", options: { color: MUTED, fontSize: 14 } },
    { text: "D-2-hydroxyglutarate dehydrogenase\n", options: { color: WHITE, bold: true, fontSize: 14 } },
    { text: "at TM-score ", options: { color: MUTED, fontSize: 14 } },
    { text: "0.90", options: { color: BLUE, bold: true, fontSize: 14 } },
    { text: "  —  both are oxidoreductases (EC 1.x).", options: { color: MUTED, fontSize: 14 } },
  ], { x: 9.4, y: 2.2, w: 3.3, h: 4.4, fontFace: BODY, valign: "top", margin: 0, lineSpacingMultiple: 1.1 });
  s.addNotes("The viewer colours each structure by pLDDT confidence and links every candidate to its structural match, pocket, and per-component scores.");
})();

/* 8 — Validation */
(() => {
  const s = dark();
  title(s, "Does it get the right answer?", "Leave-one-out benchmark: 227 enzymes, 29 families, 6 classes");
  const stat = (x, big, lab, col) => {
    s.addShape("roundRect", { x, y: 1.95, w: 3.0, h: 2.0, rectRadius: 0.1, fill: { color: PANEL }, line: { color: col, width: 1.5 } });
    s.addText(big, { x, y: 2.2, w: 3.0, h: 1.0, fontFace: HEAD, fontSize: 40, bold: true, color: col, align: "center", margin: 0 });
    s.addText(lab, { x: x + 0.15, y: 3.15, w: 2.7, h: 0.7, fontFace: BODY, fontSize: 13, color: MUTED, align: "center", margin: 0 });
  };
  stat(0.7, "96.5%", "exact-EC recovery (structure)", TEAL);
  stat(3.85, "93.8%", "sequence-identity baseline", BLUE);
  s.addText([
    { text: "Structure ties at the class level and is ", options: { color: TEXT } },
    { text: "better at pinpointing the exact reaction", options: { color: TEAL, bold: true } },
    { text: ".  We claim complementarity, not dominance.", options: { color: TEXT } },
  ], { x: 0.7, y: 4.2, w: 6.2, h: 1.5, fontFace: BODY, fontSize: 15, valign: "top", margin: 0, lineSpacingMultiple: 1.1 });
  figCard(s, FIG.acc, 7.4, 2.0, 5.2, 3.47);
  s.addNotes("Structure recovers the exact 4-level EC number for 96.5% vs a pairwise sequence baseline at 93.8% (tied at class level).");
})();

/* 9 — Twilight zone */
(() => {
  const s = dark();
  title(s, "Where structure wins", "Recovering function in the sequence \"twilight zone\"");
  figCard(s, FIG.twilight, 0.7, 1.95, 7.0, 4.9);
  const pts = [
    ["24 / 219", "correct recoveries below 30% identity"],
    ["2.8%", "identity — a lipase, still recovered correctly"],
    ["5 vs 0", "cases where structure beat sequence (never the reverse)"],
  ];
  let y = 2.2;
  pts.forEach(p => {
    s.addText(p[0], { x: 8.1, y, w: 4.6, h: 0.5, fontFace: HEAD, fontSize: 24, bold: true, color: TEAL, margin: 0 });
    s.addText(p[1], { x: 8.12, y: y + 0.55, w: 4.5, h: 0.7, fontFace: BODY, fontSize: 14, color: MUTED, margin: 0 });
    y += 1.5;
  });
  s.addNotes("Below 30% identity sequence homology is unreliable; structure still recovers correct function — including a lipase matched at under 3% identity.");
})();

/* 10 — Functional residues */
(() => {
  const s = dark();
  title(s, "Beyond fold: catalytic residues", "Verifying the active site, not just the shape");
  s.addText("For each match we superpose the candidate onto its reference and check whether the reference's annotated catalytic, metal- and ligand-binding residues are structurally conserved.",
    { x: 0.7, y: 1.9, w: 12, h: 0.8, fontFace: BODY, fontSize: 15, color: MUTED, margin: 0 });
  // Two contrast cards
  s.addShape("roundRect", { x: 0.7, y: 2.9, w: 5.9, h: 3.4, rectRadius: 0.12, fill: { color: PANEL }, line: { color: TEAL, width: 1.5 } });
  s.addText([{ text: "O08343 — strong support\n", options: { bold: true, color: TEAL, fontSize: 18 } },
    { text: "matched to a Zn amidohydrolase\n\n", options: { color: MUTED, fontSize: 14 } },
    { text: "8 / 8 functional residues position-conserved\n", options: { color: TEXT, fontSize: 15 } },
    { text: "Zn-His111, Zn-His113, catalytic Asp378\nconserved in identity", options: { color: WHITE, bold: true, fontSize: 15 } }],
    { x: 1.0, y: 3.15, w: 5.3, h: 3.0, fontFace: BODY, valign: "top", margin: 0, lineSpacingMultiple: 1.1 });
  s.addShape("roundRect", { x: 6.8, y: 2.9, w: 5.9, h: 3.4, rectRadius: 0.12, fill: { color: PANEL }, line: { color: "F0A868", width: 1.5 } });
  s.addText([{ text: "P9WJG7 — self-flagged\n", options: { bold: true, color: "F0A868", fontSize: 18 } },
    { text: "a 50-residue protein\n\n", options: { color: MUTED, fontSize: 14 } },
    { text: "0 / 6 functional residues conserved\n", options: { color: TEXT, fontSize: 15 } },
    { text: "reference active-site residues 36–55 Å away\nafter superposition", options: { color: WHITE, fontSize: 15 } },
    { text: "\n\nThe method down-weights its own weak call.", options: { color: MUTED, italic: true, fontSize: 13 } }],
    { x: 7.1, y: 3.15, w: 5.3, h: 3.0, fontFace: BODY, valign: "top", margin: 0, lineSpacingMultiple: 1.1 });
  s.addNotes("Residue-level verification strengthens confident assignments (conserved zinc-binding histidines) and automatically flags spurious ones.");
})();

/* 11 — Robustness */
(() => {
  const s = dark();
  title(s, "Robust and calibrated", "The ranking does not hinge on the score weights");
  figCard(s, FIG.weights, 0.7, 1.95, 8.0, 3.2);
  s.addText([
    { text: "Kendall τ ≈ 0.98\n", options: { color: TEAL, bold: true, fontSize: 22, fontFace: HEAD } },
    { text: "across all reasonable weightings\n(top-10 unchanged)\n\n", options: { color: MUTED, fontSize: 13 } },
    { text: "Precision 0.982\n", options: { color: BLUE, bold: true, fontSize: 22, fontFace: HEAD } },
    { text: "at the default TM ≥ 0.5 cutoff\n(1.00 at TM ≥ 0.75)", options: { color: MUTED, fontSize: 13 } },
  ], { x: 9.0, y: 2.1, w: 3.7, h: 3.2, fontFace: BODY, valign: "top", margin: 0, lineSpacingMultiple: 1.05 });
  s.addText("Weights tune emphasis, not outcome — and the confidence threshold keeps precision high while abstaining on ambiguous folds.",
    { x: 0.7, y: 5.6, w: 12, h: 0.6, fontFace: BODY, fontSize: 13.5, color: MUTED, margin: 0 });
  s.addNotes("Across the reasonable weight neighbourhood the ranking is essentially unchanged (median Kendall tau 0.98); precision at the default threshold is 0.982.");
})();

/* 12 — How to use */
(() => {
  const s = dark();
  title(s, "Using it yourself", "Three commands on a laptop — no GPU, no admin rights");
  const steps = [
    ["1", "Install", "python -m pip install -r requirements.txt"],
    ["2", "Run the pipeline", "python -m pipeline.run all"],
    ["3", "Explore results", "open web/index.html   (double-click)"],
  ];
  let y = 1.95;
  steps.forEach(st => {
    s.addShape("roundRect", { x: 0.7, y, w: 0.7, h: 0.9, rectRadius: 0.1, fill: { color: PANEL2 }, line: { color: TEAL, width: 1.5 } });
    s.addText(st[0], { x: 0.7, y, w: 0.7, h: 0.9, fontFace: HEAD, fontSize: 26, bold: true, color: TEAL, align: "center", valign: "middle", margin: 0 });
    s.addText(st[1], { x: 1.6, y: y + 0.02, w: 4.0, h: 0.9, fontFace: BODY, fontSize: 17, bold: true, color: WHITE, valign: "middle", margin: 0 });
    s.addShape("roundRect", { x: 5.7, y: y + 0.1, w: 7.0, h: 0.7, rectRadius: 0.06, fill: { color: "0B0F18" }, line: { color: PANEL2, width: 1 } });
    s.addText(st[2], { x: 5.9, y: y + 0.1, w: 6.7, h: 0.7, fontFace: MONO, fontSize: 14, color: TEAL, valign: "middle", margin: 0 });
    y += 1.15;
  });
  s.addText([
    { text: "Swap organism:  ", options: { color: MUTED } },
    { text: "edit the organism block in config.yaml", options: { color: TEXT } },
    { text: "     ·     Faster/at scale:  ", options: { color: MUTED } },
    { text: "set backend: foldseek / fpocket", options: { color: TEXT } },
  ], { x: 0.7, y: 5.7, w: 12, h: 0.6, fontFace: BODY, fontSize: 14, margin: 0 });
  s.addNotes("Install deps, run one command, open the viewer. Organism and backends are single config edits.");
})();

/* 13 — Limitations */
(() => {
  const s = dark();
  title(s, "Honest limitations", "What the method does — and does not — show");
  const items = [
    ["Predictions are hypotheses", "Structural matches suggest function; they are not experimental proof. None validated in the wet lab."],
    ["Curated benchmark", "227 well-behaved enzymes — an optimistic estimate vs. the messier real dark proteome."],
    ["EC-recovery is a proxy", "Captures enzyme function, not non-enzymatic roles, multi-domain or moonlighting proteins."],
    ["Simple sequence baseline", "Compared to pairwise identity, not profile/HMM search — a stronger baseline would narrow the gap."],
  ];
  let y = 1.95;
  items.forEach(it => {
    s.addText("▸", { x: 0.7, y, w: 0.35, h: 0.5, color: TEAL, fontSize: 15, margin: 0 });
    s.addText([{ text: it[0] + "  —  ", options: { bold: true, color: WHITE } },
               { text: it[1], options: { color: MUTED } }],
      { x: 1.1, y, w: 11.5, h: 0.9, fontFace: BODY, fontSize: 15.5, valign: "top", margin: 0, lineSpacingMultiple: 1.05 });
    y += 1.05;
  });
  s.addNotes("The framing is deliberately careful: hypotheses not proof, curated benchmark, EC proxy, simple baseline.");
})();

/* 14 — Conclusion */
(() => {
  const s = dark();
  s.addText("A credible, accessible tool", { x: 0.7, y: 1.1, w: 12, h: 0.9, fontFace: HEAD, fontSize: 40, bold: true, color: WHITE, margin: 0 });
  s.addText("Structure-based function annotation — transparent, validated, and reproducible on a laptop.",
    { x: 0.72, y: 2.05, w: 11.5, h: 0.6, fontFace: BODY, fontSize: 18, color: TEAL, margin: 0 });
  const cards = [
    ["96.5%", "exact-EC recovery"],
    ["227", "enzymes benchmarked"],
    ["0.982", "precision at TM ≥ 0.5"],
    ["0", "GPUs / HPC required"],
  ];
  cards.forEach((c, i) => {
    const x = 0.7 + i * 3.05;
    s.addShape("roundRect", { x, y: 3.1, w: 2.85, h: 2.2, rectRadius: 0.12, fill: { color: PANEL }, line: { color: i === 3 ? TEAL : PANEL2, width: i === 3 ? 2 : 1 } });
    s.addText(c[0], { x, y: 3.5, w: 2.85, h: 1.0, fontFace: HEAD, fontSize: 38, bold: true, color: i % 2 ? BLUE : TEAL, align: "center", margin: 0 });
    s.addText(c[1], { x: x + 0.1, y: 4.5, w: 2.65, h: 0.7, fontFace: BODY, fontSize: 13, color: MUTED, align: "center", margin: 0 });
  });
  s.addText("Interpretable, complementary to DeepFRI and Foldseek — every assignment traceable to a structural match and conserved catalytic residues.",
    { x: 0.7, y: 5.7, w: 12, h: 0.7, fontFace: BODY, fontSize: 14, color: MUTED, margin: 0 });
  s.addNotes("Takeaway: meaningful, interpretable, validated structure-based annotation on commodity hardware.");
})();

P.writeFile({ fileName: ROOT + "/paper/protein_function_rescue_slides.pptx" })
  .then(f => console.log("wrote", f));
