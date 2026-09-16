"""Figure 9: how often do the two models actually choose the same predictors?

Usage: python3 src/fig09_shared_predictors.py

(a) One column per gene with usable models in both sets, sorted by the number of SNPs the
    two models share. The shared SNPs form the band straddling the axis; the SNPs only the
    European model chose rise above it and the SNPs only the Yoruba model chose fall below.
    The axis is scaled to the largest model so that no column is cut off.
(b) The same overlap read as a cumulative curve: for any level on the horizontal axis, the
    height is the share of genes at or below it.
Reads results/figdata/fig09_shared_predictors.tsv.
"""
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import viz_style as vs

ROOT = Path(__file__).resolve().parents[1]
d = pd.read_csv(ROOT / "results/figdata/fig09_shared_predictors.tsv", sep="\t")
# Sorted by the COUNT of shared SNPs, with a stable sort and no secondary key: an earlier
# version broke ties by model size, which produced a repeating sawtooth that read as an artifact.
d = d.sort_values("shared", ascending=False, kind="mergesort").reset_index(drop=True)
x = np.arange(len(d))
half = d.shared / 2.0

vs.apply()
fig = plt.figure(figsize=(vs.DOUBLE, 3.6))
ax = fig.add_axes([0.075, 0.185, 0.55, 0.6])
bx = fig.add_axes([0.73, 0.185, 0.24, 0.6])

ax.bar(x, d.eur_only, bottom=half, width=1.0, color=vs.EUR, linewidth=0, zorder=3)
ax.bar(x, -d.yri_only, bottom=-half, width=1.0, color=vs.YRI, linewidth=0, zorder=3)
ax.bar(x, d.shared, bottom=-half, width=1.0, color=vs.COVERAGE, linewidth=0, zorder=4)
ax.axhline(0, color=vs.INK_2, lw=0.6, zorder=5)
n_zero = int((d.shared == 0).sum())
first_zero = int((d.shared > 0).sum())
top = float((d.eur_only + half).max()); bot = float((d.yri_only + half).max())
span = max(top, bot) * 1.06
ax.axvline(first_zero, color=vs.INK_2, lw=0.7, ls=(0, (3, 2)), zorder=6)
ax.text(first_zero + 8, -span * 0.86, f"from here rightward the two models share\nno SNP at all: "
        f"{n_zero} of {len(d)} genes", fontsize=7.5, color=vs.INK_2, ha="left", va="bottom")
ax.set_xlim(-2, len(d) + 2)
ax.set_ylim(-span, span)
step = 50 if span > 120 else 25
ticks = [t for t in range(0, int(span), step)]
ax.set_yticks([-t for t in ticks[::-1]] + ticks[1:])
ax.set_yticklabels([f"{t}" for t in ticks[::-1]] + [f"{t}" for t in ticks[1:]])
ax.yaxis.ibeji_value_of = np.abs
ax.set_ylabel("SNPs in the model")
ax.set_xlabel("the 562 genes with a usable model in both sets, sorted by shared SNPs")
ax.set_title("Two models for the same gene, built from different SNPs", fontsize=8.5, color=vs.INK, loc="left", pad=30)
for i, (label, color) in enumerate((("European model only", vs.EUR),
                                    ("chosen by both", vs.COVERAGE),
                                    ("Yoruba model only", vs.YRI))):
    xx = 0.0 + 0.335 * i
    ax.add_patch(plt.Rectangle((xx, 1.045), 0.022, 0.045, transform=ax.transAxes, color=color, lw=0, clip_on=False))
    ax.text(xx + 0.032, 1.068, label, transform=ax.transAxes, fontsize=7.5, color=vs.INK_2, va="center")

xs = np.sort(d.jaccard.values)
ys = np.arange(1, len(xs) + 1) / len(xs)
bx.step(np.concatenate([[0], xs]), np.concatenate([[0], ys]), where="post", color=vs.COVERAGE, lw=2.0, zorder=3)
med = float(d.jaccard.median())
frac0 = (d.jaccard == 0).mean()
bx.plot([0, 0], [0, frac0], color=vs.COVERAGE, lw=2.0, zorder=3)
bx.scatter([0], [frac0], s=18, color=vs.COVERAGE, edgecolors=vs.SURFACE, linewidths=0.7, zorder=5)
bx.annotate(f"{frac0:.0%} share nothing", xy=(0, frac0), xytext=(0.09, frac0 - 0.16), fontsize=7.5,
            color=vs.INK_2, arrowprops=dict(arrowstyle="-", color=vs.INK_2, lw=0.6, shrinkA=1, shrinkB=3))
bx.plot([med, med], [0, 0.5], color=vs.INK, lw=0.7, ls=(0, (3, 2)), zorder=4)
bx.text(med + 0.02, 0.5, f"half the genes\nshare {med:.1%} or less", fontsize=7.5, color=vs.INK, va="center")
bx.set_xlim(-0.01, 0.55)
bx.set_ylim(0, 1.02)
bx.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
bx.set_yticklabels(["0%", "25%", "50%", "75%", "100%"])
bx.yaxis.ibeji_value_of = lambda v: v * 100
bx.set_xlabel("share of a gene's SNPs\nchosen by both models")
bx.set_ylabel("genes at or below")
bx.grid(True, axis="y"); bx.set_axisbelow(True)
bx.set_title("Overlap is small", fontsize=8.5, color=vs.INK, loc="left", pad=30)
for a_, letter in ((ax, "a"), (bx, "b")):
    fig.text(a_.get_position().x0 - 0.058, 0.985, letter, fontsize=9, fontweight="bold", color=vs.INK, va="top")
vs.save(fig, "fig09_shared_predictors")
print(f"genes {len(d)}, zero-overlap {n_zero} ({n_zero/len(d):.1%}), median jaccard {med:.4f}, "
      f"median shared {d.shared.median():.0f}, max shared {d.shared.max()}")
