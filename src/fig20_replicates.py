"""Figure 20: does the answer depend on which 87 Europeans you happened to draw?

Usage: python3 src/fig20_replicates.py

The same decomposition is run twice, against the same Yoruba sample, using two independent
random draws of 87 European individuals. Anything that changes between them is a property of
the draw, not of ancestry.

(a) Each gene's weights share in one draw against the other.
(b) Which component is largest for a gene, draw against draw.
(c) The median share of each component, per draw.

Any pair of completed EUR87 draws is used, so this extends as further draws finish.
"""
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats

import viz_style as vs

ROOT = Path(__file__).resolve().parents[1]
D = ROOT / "results" / "decomposition"
COMP = [("weights", vs.PHI_W), ("allele frequency", vs.PHI_D), ("LD", vs.PHI_R)]

draws = sorted(f for f in D.glob("EUR87_r*_vs_YRI87.tsv"))
assert len(draws) >= 2, "need at least two European draws decomposed against YRI87"


def shares(path):
    d = pd.read_csv(path, sep="\t").dropna(subset=["phi_w", "phi_D", "phi_R"])
    t = d[["phi_w", "phi_D", "phi_R"]].abs()
    out = t.div(t.sum(axis=1), axis=0)
    out.columns = ["s_w", "s_D", "s_R"]
    out["gene"] = d.gene.values
    out["dom"] = t.values.argmax(axis=1)
    return out


sa, sb = shares(draws[0]), shares(draws[1])
m = sa.merge(sb, on="gene", suffixes=("_1", "_2"))
agree = float((m.dom_1 == m.dom_2).mean())

vs.apply()
fig = plt.figure(figsize=(vs.DOUBLE, 3.25))
ax = fig.add_axes([0.075, 0.20, 0.235, 0.52])
bx = fig.add_axes([0.435, 0.20, 0.175, 0.52])
cx = fig.add_axes([0.755, 0.20, 0.205, 0.52])

ax.plot([0, 1], [0, 1], color=vs.MUTED, lw=0.7, zorder=2)
ax.scatter(m.s_w_1, m.s_w_2, s=7, color=vs.PHI_W, alpha=0.45, linewidths=0, zorder=3)
rho = stats.spearmanr(m.s_w_1, m.s_w_2)[0]
ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_aspect("equal")
ax.set_xlabel("weights share, draw 1")
ax.set_ylabel("weights share, draw 2")
ax.grid(True); ax.set_axisbelow(True)
ax.set_title("The same gene, two draws", fontsize=8.5, color=vs.INK, loc="left", pad=40)
ax.text(0.0, 1.03, f"{len(m)} genes;  Spearman {rho:.2f}\nmedian change {float((m.s_w_1-m.s_w_2).abs().median()):.2f}",
        transform=ax.transAxes, fontsize=7.5, color=vs.INK_2, va="bottom")

M = np.zeros((3, 3), int)
for i in range(3):
    for j in range(3):
        M[i, j] = int(((m.dom_1 == i) & (m.dom_2 == j)).sum())
ramp = plt.matplotlib.colors.LinearSegmentedColormap.from_list("c", ["#eef2fb", "#9ec5f4", "#2a78d6", "#123f77"])
bx.imshow(M, cmap=ramp, vmin=0, vmax=M.max(), origin="upper")
for i in range(3):
    for j in range(3):
        bx.text(j, i, str(M[i, j]), ha="center", va="center", fontsize=7.5,
                color=vs.SURFACE if M[i, j] > M.max() * 0.45 else vs.INK)
# Symbols, horizontal. Three rotated words cannot fit across a 1.25 in panel at any
# rotation, the same trap as Figure 19; the symbols are already the paper's notation.
short = [r"$\phi_w$", r"$\phi_D$", r"$\phi_R$"]
bx.set_xticks(range(3)); bx.set_yticks(range(3))
bx.set_xticklabels(short, fontsize=8)
bx.set_yticklabels(short, fontsize=8)
bx.tick_params(length=0)
for sp in bx.spines.values():
    sp.set_visible(False)
bx.set_xlabel("largest in draw 2")
bx.set_ylabel("largest in draw 1")
bx.set_title("Which component is largest", fontsize=8.5, color=vs.INK, loc="left", pad=40)
bx.text(0.0, 1.03, f"agrees for {agree:.0%} of genes", transform=bx.transAxes, fontsize=7.5,
        color=vs.INK_2, va="bottom")

for i, (name, color) in enumerate(COMP):
    col = ["s_w", "s_D", "s_R"][i]
    v1, v2 = float(sa[col].median()), float(sb[col].median())
    y = 2 - i
    cx.plot([v1, v2], [y, y], color=color, lw=2.0, solid_capstyle="round", zorder=3)
    cx.scatter([v1, v2], [y, y], s=26, color=color, edgecolors=vs.SURFACE, linewidths=0.8, zorder=4)
    cx.text(max(v1, v2) + 0.035, y, f"{v1:.2f} / {v2:.2f}", va="center", fontsize=7.5, color=vs.INK)
    cx.text(-0.02, y + 0.30, name, va="center", fontsize=7.5, color=vs.INK_2,
            transform=cx.get_yaxis_transform())
cx.set_xlim(0, 0.95); cx.set_ylim(-0.6, 2.75)
cx.set_yticks([])
cx.spines["left"].set_visible(False)
cx.set_xlabel("median share of the total")
cx.grid(True, axis="x"); cx.set_axisbelow(True)
cx.set_title("The population answer", fontsize=8.5, color=vs.INK, loc="left", pad=40)
# Said in words because the two dots nearly coincide: that near-coincidence is the finding,
# but without saying so it reads as a single dot or a rendering fault.
cx.text(0.0, 1.03, "one dot per draw\nthey nearly coincide",
        transform=cx.transAxes, fontsize=7.5, color=vs.INK_2, va="bottom")

for a_, letter in ((ax, "a"), (bx, "b"), (cx, "c")):
    fig.text(max(a_.get_position().x0 - 0.058, 0.008), 0.975, letter, fontsize=9, fontweight="bold",
             color=vs.INK, va="top")
vs.save(fig, "fig20_replicates")
print(f"genes {len(m)}, dominant agrees {agree:.1%}, weights-share spearman {rho:.3f}")
print("median shares draw1:", {c: round(float(sa[c].median()), 3) for c in ("s_w", "s_D", "s_R")})
print("median shares draw2:", {c: round(float(sb[c].median()), 3) for c in ("s_w", "s_D", "s_R")})
