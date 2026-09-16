"""Figure 18: why is a transcriptome-wide association found with one ancestry's models and not the other?

Usage: python3 src/fig18_specific_hits.py

Every association significant with one model set and not the other is placed in one of two
situations. Either the other ancestry never produced a usable model for that gene, in which
case no comparison of weights, frequencies or LD is even defined, or both ancestries modelled
the gene and the disagreement can be attributed to the largest component of the decomposition.

(a) The composition for each disorder, split by which model set found the association.
(b) The same totals as a single ranking, which is the summary the text quotes.
Reads results/figdata/fig18_specific_hits.tsv.
"""
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import viz_style as vs

ROOT = Path(__file__).resolve().parents[1]
r = pd.read_csv(ROOT / "results/figdata/fig18_specific_hits.tsv", sep="\t")

ORDER = [("no usable model in the other ancestry", vs.COVERAGE, "no usable model\nin the other ancestry"),
         ("weights", vs.PHI_W, "weights"),
         ("allele frequency", vs.PHI_D, "allele frequency"),
         ("linkage disequilibrium", vs.PHI_R, "linkage disequilibrium"),
         ("decomposition undefined", vs.MUTED, "decomposition undefined")]
TRAITS = ["Schizophrenia", "Bipolar disorder", "Major depression", "PTSD"]
SIDES = [("European models", vs.EUR), ("Yoruba models", vs.YRI)]

# Both panels iterate over ORDER only, and panel b takes its percentage denominator from the
# classes it finds there, so a class missing from ORDER would vanish from the figure AND
# quietly change every percentage. This figure has already been redesigned once for silently
# dropping rows; refuse to draw rather than do it again.
unknown = sorted(set(r.klass) - {k for k, _, _ in ORDER})
if unknown:
    raise SystemExit(f"classes absent from ORDER would be dropped silently: {unknown}")

vs.apply()
fig = plt.figure(figsize=(vs.DOUBLE, 3.7))
ax = fig.add_axes([0.275, 0.175, 0.40, 0.55])
bx = fig.add_axes([0.775, 0.175, 0.195, 0.55])

# Labels use get_yaxis_transform: x in axes fraction, y in data. Placing them in data
# coordinates put the trait names in the same region as the row labels and the bar values.
T = ax.get_yaxis_transform()
y = 0.0
for trait in TRAITS:
    y -= 0.55
    first_y = y
    for side, scol in SIDES:
        sub = r[(r.trait == trait) & (r.found_with == side)]
        x0 = 0.0
        for klass, color, _ in ORDER:
            n = int((sub.klass == klass).sum())
            if n:
                ax.barh([y], [n], left=x0, height=0.72, color=color, linewidth=0, zorder=3)
                if n >= 4:
                    ax.text(x0 + n / 2, y, str(n), ha="center", va="center", fontsize=7,
                            color=vs.SURFACE, zorder=4)
                x0 += n
        ax.text(-0.015, y, side.replace(" models", ""), transform=T, ha="right", va="center",
                fontsize=7.5, color=vs.INK_2)
        y -= 1.0
    # Trait name sits on its own line above its two rows. Beside them it cannot fit: the name
    # needs about 1.05 in and the row labels another 0.55 in, in a 0.86 in margin.
    ax.text(-0.30, first_y + 0.92, trait, transform=T, ha="left", va="center",
            fontsize=8, color=vs.INK)
    y -= 0.75

ax.set_ylim(y + 0.4, 0.8)
ax.set_yticks([])
ax.spines["left"].set_visible(False)
ax.set_xlabel("associations found with one model set and not the other")
ax.grid(True, axis="x"); ax.set_axisbelow(True)
# Shortened: the long form reached across the gap and collided with panel b's letter.
ax.set_title("Why an association is ancestry-specific", fontsize=8.5,
             color=vs.INK, loc="left", pad=64)
# Two rows: the first label alone is about 1.9 in against a 2.86 in panel, so three
# entries cannot share a line.
LEG = [(0.0, 1.235, vs.COVERAGE, "no usable model in the other ancestry"),
       (0.0, 1.150, vs.PHI_W, "weights"),
       (0.34, 1.150, vs.PHI_D, "allele frequency"),
       (0.0, 1.065, vs.MUTED, "decomposition undefined")]
for x, yy, color, label in LEG:
    ax.add_patch(plt.Rectangle((x, yy - 0.024), 0.024, 0.048, transform=ax.transAxes, color=color,
                               lw=0, clip_on=False))
    ax.text(x + 0.034, yy, label, transform=ax.transAxes, fontsize=7.5, color=vs.INK_2, va="center")

counts = [(lab, int((r.klass == k).sum()), c) for k, c, lab in ORDER]
counts = [c for c in counts if c[1] > 0]
tot = sum(c[1] for c in counts)
yy = np.arange(len(counts))[::-1]
for (label, n, color), y_ in zip(counts, yy):
    bx.barh([y_], [n / tot * 100], height=0.6, color=color, linewidth=0, zorder=3)
    bx.text(n / tot * 100 + 1.5, y_, f"{n}  ({n/tot:.0%})", va="center", fontsize=7.5, color=vs.INK)
bx.set_yticks([])
bx.spines["left"].set_visible(False)
bx.set_xlim(0, 100)
bx.set_xticks([0, 50, 100]); bx.set_xticklabels(["0%", "50%", "100%"])
bx.xaxis.ibeji_value_of = lambda v: v
bx.set_xlabel(f"share of all {tot}")
bx.grid(True, axis="x"); bx.set_axisbelow(True)
bx.set_title("All disorders", fontsize=8.5, color=vs.INK, loc="left", pad=64)
bx.text(0.0, 1.055, "colours as in panel a", transform=bx.transAxes, fontsize=7.5, color=vs.INK_2, va="bottom")

for a_, letter in ((ax, "a"), (bx, "b")):
    fig.text(max(a_.get_position().x0 - 0.052, 0.010), 0.975, letter, fontsize=9, fontweight="bold",
             color=vs.INK, va="top")
vs.save(fig, "fig18_specific_hits")
print(pd.DataFrame(counts, columns=["class", "n", "color"]).drop(columns=["color"]).to_string(index=False))
print(f"total {tot}")
