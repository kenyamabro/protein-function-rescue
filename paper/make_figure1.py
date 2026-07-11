"""Generate the pipeline schematic (Figure 1) for the manuscript."""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

OUT = Path(__file__).resolve().parent / "figures" / "fig1_pipeline.png"
OUT.parent.mkdir(parents=True, exist_ok=True)

stages = [
    ("download", "AlphaFold DB\nstructures"),
    ("parse", "mean pLDDT,\nsequence"),
    ("annotate", "UniProt name;\nflag unannotated"),
    ("search", "structural hit\n(known fold)"),
    ("pocket", "binding cavity"),
    ("rank", "composite\nscore"),
]

fig, ax = plt.subplots(figsize=(11, 3.1))
ax.set_xlim(0, len(stages) * 2.0)
ax.set_ylim(0, 3)
ax.axis("off")

box_w, box_h, y = 1.7, 1.15, 1.2
centers = []
for i, (name, sub) in enumerate(stages):
    x = i * 2.0 + 0.15
    cx = x + box_w / 2
    centers.append(cx)
    highlight = name in ("search", "pocket")
    face = "#dceefb" if not highlight else "#cfe8d6"
    edge = "#2b6cb0" if not highlight else "#2f855a"
    ax.add_patch(FancyBboxPatch((x, y), box_w, box_h,
                                boxstyle="round,pad=0.03,rounding_size=0.12",
                                linewidth=1.6, edgecolor=edge, facecolor=face))
    ax.text(cx, y + box_h - 0.32, name, ha="center", va="center",
            fontsize=12, fontweight="bold", color=edge)
    ax.text(cx, y + 0.32, sub, ha="center", va="center", fontsize=8.3, color="#333")
    if highlight:
        ax.text(cx, y - 0.28, "2 backends", ha="center", va="center",
                fontsize=7.5, style="italic", color="#2f855a")

for a, b in zip(centers[:-1], centers[1:]):
    ax.add_patch(FancyArrowPatch((a + box_w / 2, y + box_h / 2),
                                 (b - box_w / 2, y + box_h / 2),
                                 arrowstyle="-|>", mutation_scale=16,
                                 linewidth=1.4, color="#555"))

ax.text(len(stages), 2.75, "Protein Function Rescue pipeline",
        ha="center", fontsize=13, fontweight="bold")
ax.text(len(stages), 0.18,
        "search: TM-align (native) or Foldseek     pocket: LIGSITE-style geometry (native) or fpocket",
        ha="center", fontsize=8.2, color="#2f855a")

fig.tight_layout()
fig.savefig(OUT, dpi=150, bbox_inches="tight")
print("wrote", OUT)
