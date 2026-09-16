"""Figure 20: does the answer depend on which 87 Europeans you happened to draw?

Usage: python3 src/fig20_replicates.py

The same decomposition is run against the same Yoruba sample using independent random draws
of 87 European individuals. Anything that changes between them is a property of the draw,
not of ancestry.

(a) Each gene's weights share in one draw against another, pooled over every pair of draws.
(b) Which component is largest for a gene, draw against draw, pooled over every pair.
(c) The median share of each component, one dot per draw.

Every completed EUR87 draw is used. Panels a and b pool all C(n, 2) pairs and panel c plots
one dot per draw, so with two draws this is a single pair and with five it is ten pairs.
(An earlier version globbed all the draws but then used only the first two, in every panel,
and both this docstring and the figure registry claimed it extended itself. It did not.)
"""
import itertools
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


per_draw = [shares(p) for p in draws]
n_draws = len(per_draw)

# Pool every unordered pair. With two draws this is the single pair the figure used before,
# so the rendering is unchanged; with five draws it is ten pairs and the agreement rate and
# the scatter both rest on far more comparisons.
merged = [a.merge(b, on="gene", suffixes=("_1", "_2")) for a, b in itertools.combinations(per_draw, 2)]
m = pd.concat(merged, ignore_index=True)
n_pairs = len(merged)
n_genes = m.gene.nunique()
agree = float((m.dom_1 == m.dom_2).mean())
rho = stats.spearmanr(m.s_w_1, m.s_w_2)[0]
med_change = float((m.s_w_1 - m.s_w_2).abs().median())

vs.apply()
fig = plt.figure(figsize=(vs.DOUBLE, 3.25))
ax = fig.add_axes([0.075, 0.20, 0.235, 0.52])
bx = fig.add_axes([0.435, 0.20, 0.175, 0.52])
cx = fig.add_axes([0.755, 0.20, 0.205, 0.52])

# Pooling multiplies the point count by the number of pairs, so thin the ink to match or the
# cloud turns into a solid block and the density near the diagonal stops reading.
alpha = 0.45 if n_pairs == 1 else max(0.10, 0.45 / np.sqrt(n_pairs))
msize = 7 if n_pairs == 1 else max(3.0, 7 / np.sqrt(n_pairs))

ax.plot([0, 1], [0, 1], color=vs.MUTED, lw=0.7, zorder=2)
ax.scatter(m.s_w_1, m.s_w_2, s=msize, color=vs.PHI_W, alpha=alpha, linewidths=0, zorder=3)
ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_aspect("equal")
# With a single pair the axes really are draw 1 and draw 2, which is what the caption says.
# Pooling over more pairs makes those names meaningless, so they lose their numbers. Panel b
# below plots the same two draws on the opposite axes, so it reuses these in reverse.
lab_1, lab_2 = ("draw 1", "draw 2") if n_pairs == 1 else ("one draw", "the other")
ax.set_xlabel(f"weights share, {lab_1}")
ax.set_ylabel(f"weights share, {lab_2}")
ax.grid(True); ax.set_axisbelow(True)
title_a = "The same gene, two draws" if n_pairs == 1 else f"The same gene, every pair of {n_draws} draws"
ax.set_title(title_a, fontsize=8.5, color=vs.INK, loc="left", pad=40)
count_a = f"{len(m)} genes" if n_pairs == 1 else f"{n_genes} genes, {n_pairs} pairs"
ax.text(0.0, 1.03, f"{count_a};  Spearman {rho:.2f}\nmedian change {med_change:.2f}",
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
bx.set_xlabel(f"largest in {lab_2}")
bx.set_ylabel(f"largest in {lab_1}")
bx.set_title("Which component is largest", fontsize=8.5, color=vs.INK, loc="left", pad=40)
bx.text(0.0, 1.03, f"agrees for {agree:.0%} of genes", transform=bx.transAxes, fontsize=7.5,
        color=vs.INK_2, va="bottom")

widest = 0.0
for i, (name, color) in enumerate(COMP):
    col = ["s_w", "s_D", "s_R"][i]
    meds = [float(d[col].median()) for d in per_draw]
    widest = max(widest, max(meds) - min(meds))
    y = 2 - i
    cx.plot([min(meds), max(meds)], [y, y], color=color, lw=2.0, solid_capstyle="round", zorder=3)
    cx.scatter(meds, [y] * len(meds), s=26, color=color, edgecolors=vs.SURFACE, linewidths=0.8, zorder=4)
    label = " / ".join(f"{v:.2f}" for v in meds) if n_draws == 2 else f"{min(meds):.2f} to {max(meds):.2f}"
    cx.text(max(meds) + 0.035, y, label, va="center", fontsize=7.5, color=vs.INK)
    cx.text(-0.02, y + 0.30, name, va="center", fontsize=7.5, color=vs.INK_2,
            transform=cx.get_yaxis_transform())
cx.set_xlim(0, 0.95); cx.set_ylim(-0.6, 2.75)
cx.set_yticks([])
cx.spines["left"].set_visible(False)
cx.set_xlabel("median share of the total")
cx.grid(True, axis="x"); cx.set_axisbelow(True)
cx.set_title("The population answer", fontsize=8.5, color=vs.INK, loc="left", pad=40)
# Said in words because the dots nearly coincide: that near-coincidence is the finding,
# but without saying so it reads as a single dot or a rendering fault. The magnitude is
# already on the panel, in the per-component labels beside the dots. Spelling it out a
# second time here ran the line past the figure edge (check 2 caught it) and widened the
# whole image, so the measured spread goes to stdout instead.
cx.text(0.0, 1.03, "one dot per draw\nthey nearly coincide",
        transform=cx.transAxes, fontsize=7.5, color=vs.INK_2, va="bottom")

for a_, letter in ((ax, "a"), (bx, "b"), (cx, "c")):
    fig.text(max(a_.get_position().x0 - 0.058, 0.008), 0.975, letter, fontsize=9, fontweight="bold",
             color=vs.INK, va="top")
vs.save(fig, "fig20_replicates")
print(f"draws {n_draws} ({', '.join(p.name.split('_vs_')[0] for p in draws)}), pairs {n_pairs}")
print(f"genes {n_genes}, pooled comparisons {len(m)}, dominant agrees {agree:.1%}, "
      f"weights-share spearman {rho:.3f}, median change {med_change:.3f}")
print(f"widest spread of a median share across draws: {widest:.3f}")
for p, d in zip(draws, per_draw):
    print(f"  {p.name.split('_vs_')[0]:10s}", {c: round(float(d[c].median()), 3) for c in ("s_w", "s_D", "s_R")})
