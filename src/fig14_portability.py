"""Figure 14: what decides whether a model still works in the other population?

Usage: python3 src/fig14_portability.py

(a) Cross-population accuracy against the number of SNPs the two models share, in both
    directions, as a median with the middle half shaded.
(b) Partial rank correlation of each candidate predictor with cross-population accuracy,
    holding the model's accuracy at home fixed, with bootstrap intervals.

The outcome is cross-population accuracy itself, not the drop from home accuracy. Home
accuracy is a nested cross-validation estimate while the cross-population number comes from
the final model evaluated in an independent sample, so their difference is not a clean loss:
it is negative for 28% of genes one way and 48% the other, tracking shared-SNP count.
Reads results/figdata/fig14_portability.tsv.
"""
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats

import viz_style as vs

ROOT = Path(__file__).resolve().parents[1]
m = pd.read_csv(ROOT / "results/figdata/fig14_portability.tsv", sep="\t")
RNG = np.random.default_rng(1)

DIRS = [("r2_eur_model_in_yri", "home_E", vs.EUR, "European model in Yoruba"),
        ("r2_yri_model_in_eur", "home_Y", vs.YRI, "Yoruba model in European")]
# Kept short: a y tick label extends leftward out of its own panel and lands on whatever
# sits beside it, which is how the first version wrote across panel a.
PREDS = [("n_shared", "shared SNPs"),
         ("fst_weighted", "frequency divergence"),
         ("ld_divergence", "LD divergence")]


def partial_spearman(df, x, y, controls):
    s = df.dropna(subset=[x, y] + controls)
    R = {c: stats.rankdata(s[c]) for c in [x, y] + controls}
    C = np.column_stack([R[c] for c in controls] + [np.ones(len(s))])
    res = lambda v: v - C @ np.linalg.lstsq(C, v, rcond=None)[0]
    return float(stats.pearsonr(res(R[x]), res(R[y]))[0])


def boot_ci(df, x, y, controls, n=800):
    vals = np.empty(n)
    idx = np.arange(len(df))
    for i in range(n):
        vals[i] = partial_spearman(df.iloc[RNG.choice(idx, len(idx), replace=True)], x, y, controls)
    return float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))


vs.apply()
fig = plt.figure(figsize=(vs.DOUBLE, 3.6))
ax = fig.add_axes([0.06, 0.185, 0.32, 0.60])
bx = fig.add_axes([0.615, 0.185, 0.335, 0.60])

EDGES = [-0.5, 0.5, 1.5, 2.5, 4.5, 8.5, 16.5, 10000]
LABELS = ["0", "1", "2", "3-4", "5-8", "9-16", "17+"]
xs = np.arange(len(LABELS))
for acc, home, color, label in DIRS:
    s = m.dropna(subset=[acc]).copy()
    s["b"] = pd.cut(s.n_shared, EDGES, labels=False)
    g = s.groupby("b").agg(med=(acc, "median"), lo=(acc, lambda v: v.quantile(0.25)),
                           hi=(acc, lambda v: v.quantile(0.75)), n=(acc, "size"))
    gx = xs[g.index.values]
    ax.fill_between(gx, g.lo, g.hi, color=color, alpha=0.16, lw=0, zorder=2)
    ax.plot(gx, g.med, color=color, lw=2.0, zorder=3, solid_capstyle="round")
    ax.scatter(gx, g.med, s=16, color=color, edgecolors=vs.SURFACE, linewidths=0.7, zorder=4)
ax.axhline(0, color=vs.INK_2, lw=0.6, zorder=1)
ax.set_xticks(xs); ax.set_xticklabels(LABELS)
ax.set_xlabel("SNPs chosen by both models")
ax.set_ylabel("accuracy in the other population (signed r²)")
ax.grid(True, axis="y"); ax.set_axisbelow(True)
ax.set_title("Sharing predictors is what carries a model across", fontsize=8.5, color=vs.INK, loc="left", pad=30)
ax.text(0.0, 1.03, "median, with the middle half shaded", transform=ax.transAxes, fontsize=7.5,
        color=vs.INK_2, va="bottom")
# Placed in the open lower-right of the panel: both curves rise away from it, so the key
# sits on the surface rather than on top of the shaded bands as it did in the first version.
for i, (acc, home, color, label) in enumerate(DIRS):
    y = 0.215 - 0.085 * i
    ax.add_patch(plt.Rectangle((0.40, y - 0.016), 0.035, 0.030, transform=ax.transAxes, color=color, lw=0,
                               zorder=6))
    ax.text(0.45, y, label.replace("\n", " "), transform=ax.transAxes, fontsize=7.5, color=vs.INK_2,
            va="center", zorder=6)

rows = []
for acc, home, color, label in DIRS:
    s = m.dropna(subset=[acc, home]).copy()
    for pred, pname in PREDS:
        r = partial_spearman(s, pred, acc, [home])
        lo, hi = boot_ci(s, pred, acc, [home])
        rows.append(dict(direction=label.replace("\n", " "), predictor=pname, rho=r, lo=lo, hi=hi, color=color, n=len(s)))
res = pd.DataFrame(rows)

bx.axvline(0, color=vs.INK_2, lw=0.7, zorder=2)
ypos, yticks, ylabels = [], [], []
row_y = 0.0
# Reversed so the strongest predictor reads at the top rather than the bottom.
for pi, (pred, pname) in enumerate(reversed(PREDS)):
    block = res[res.predictor == pname]
    centre = row_y + 0.5
    for k, (_, r) in enumerate(block.iterrows()):
        yy = row_y + (1 - k) * 0.42
        bx.plot([r.lo, r.hi], [yy, yy], color=r.color, lw=2.0, solid_capstyle="round", zorder=3)
        bx.scatter([r.rho], [yy], s=22, color=r.color, edgecolors=vs.SURFACE, linewidths=0.8, zorder=4)
    yticks.append(centre - 0.29); ylabels.append(pname)
    row_y += 1.35
bx.set_yticks(yticks); bx.set_yticklabels(ylabels)
bx.tick_params(axis="y", length=0)
bx.set_ylim(-0.45, row_y - 0.45)
bx.set_xlim(-0.35, 0.8)
bx.set_xlabel("partial rank correlation with\naccuracy in the other population")
bx.grid(True, axis="x"); bx.set_axisbelow(True)
bx.set_title("Holding accuracy at home fixed", fontsize=8.5, color=vs.INK, loc="left", pad=30)
bx.text(0.0, 1.03, "bar: 95% bootstrap interval;  colours as in panel a", transform=bx.transAxes,
        fontsize=7.5, color=vs.INK_2, va="bottom")
for a_, letter in ((ax, "a"), (bx, "b")):
    # Clamped: a fixed offset put the letter off the canvas once panel a moved left.
    fig.text(max(a_.get_position().x0 - 0.052, 0.010), 0.985, letter, fontsize=9, fontweight="bold",
             color=vs.INK, va="top")
vs.save(fig, "fig14_portability")

res.drop(columns=["color"]).to_csv(ROOT / "results/figdata/fig14_summary.tsv", sep="\t", index=False)
print(res.drop(columns=["color"]).round(3).to_string(index=False))
