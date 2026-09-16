"""Figure 13: crossing the ocean. How much of each model's accuracy survives in the other population.

Usage: python3 src/fig13_transfer.py

(a) Every usable European model (EUR87_r1) applied to the Yoruba individuals.
(b) Every usable Yoruba model (YRI87) applied to the European individuals (EUR358).
x: the model's cross-validated R^2 in its own population (square-root axis, declared transform).
y: signed r^2 between its predictions and observed expression in the other population
   (negative when the prediction runs the wrong way).
Hexagons count genes. The diagonal is perfect transfer, the horizontal line is no signal, and
the thick line is the median transfer among models of similar home accuracy.
Reads results/transfer/*.tsv written by src/10_transfer_all.R.
"""
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, LogNorm

import viz_style as vs

ROOT = Path(__file__).resolve().parents[1]
T = ROOT / "results" / "transfer"
OUT = ROOT / "results" / "figdata"
XMAX = np.sqrt(0.8)
YLIM = (-0.25, 0.8)


def panel(ax, path, line_color, title):
    d = pd.read_csv(path, sep="\t")
    d = d[d["tested"]].copy()
    x = np.sqrt(d["source_cv_r2"].clip(lower=0))
    y = d["r2_signed"]
    cmap = LinearSegmentedColormap.from_list("dens", ["#eef0f3", "#b9bdc6", "#6f7581", "#2f3440"])
    hb = ax.hexbin(x, y, gridsize=(30, 22), extent=(0, XMAX, *YLIM), cmap=cmap, norm=LogNorm(vmin=1), mincnt=1,
                   linewidths=0, rasterized=True)
    # running median of transfer by home accuracy, in deciles of the source R^2
    d["bin"] = pd.qcut(d["source_cv_r2"], 10, labels=False, duplicates="drop")
    g = d.groupby("bin").agg(x=("source_cv_r2", "median"), med=("r2_signed", "median"),
                             q25=("r2_signed", lambda v: v.quantile(0.25)), q75=("r2_signed", lambda v: v.quantile(0.75)))
    gx = np.sqrt(g["x"])
    ax.fill_between(gx, g["q25"], g["q75"], color=line_color, alpha=0.18, lw=0, zorder=3)
    ax.plot(gx, g["med"], color=line_color, lw=2.0, zorder=4, solid_capstyle="round")
    ax.scatter(gx, g["med"], s=12, color=line_color, edgecolors=vs.SURFACE, linewidths=0.6, zorder=5)
    # Full transfer means transfer r^2 == home R^2. With x at sqrt(R^2) that is the curve y = x^2,
    # so it is drawn as a curve through many points, never as a straight segment between the ends.
    xs = np.linspace(0, XMAX, 200)
    ax.plot(xs, xs ** 2, color=vs.MUTED, lw=0.7, zorder=2)
    ax.axhline(0, color=vs.INK_2, lw=0.6, zorder=2)
    ticks = [0, 0.05, 0.1, 0.2, 0.4, 0.6]
    ax.set_xticks(np.sqrt(ticks))
    ax.set_xticklabels([f"{t:g}" for t in ticks])
    ax.xaxis.ibeji_value_of = np.square
    ax.set_xlim(0, XMAX)
    ax.set_ylim(*YLIM)
    ax.grid(True, axis="y")
    ax.set_axisbelow(True)
    ax.set_title(title, fontsize=8, color=vs.INK, loc="left", pad=30)
    ax.text(0.0, 1.03, f"{len(d):,} models;  transfer r² median {y.median():.3f}, mean {y.mean():.3f}\n"
            f"{(d['p'] < 0.05).mean():.0%} reach nominal p < 0.05 in the other population", transform=ax.transAxes,
            fontsize=7, color=vs.INK_2, va="bottom")
    return hb, {"models": len(d), "median": y.median(), "mean": y.mean(), "frac_p05": (d["p"] < 0.05).mean(),
                "frac_negative": (y < 0).mean(), "top_decile_median_transfer": g["med"].iloc[-1],
                "top_decile_median_home_r2": g["x"].iloc[-1]}


vs.apply()
fig = plt.figure(figsize=(vs.DOUBLE, 4.1))
ax_a = fig.add_axes([0.08, 0.235, 0.37, 0.55])
ax_b = fig.add_axes([0.53, 0.235, 0.36, 0.55])
cax = fig.add_axes([0.915, 0.3, 0.012, 0.36])

hb, sa = panel(ax_a, T / "EUR87_r1_into_YRI87.tsv", vs.EUR, "European models, predicting Yoruba expression")
_, sb = panel(ax_b, T / "YRI87_into_EUR358.tsv", vs.YRI, "Yoruba models, predicting European expression")
ax_a.set_xlabel("cross-validated R² at home (European draw)")
ax_b.set_xlabel("cross-validated R² at home (Yoruba sample)")
ax_a.set_ylabel("signed r² in the other population")

cb = fig.colorbar(hb, cax=cax)
cb.set_label("models per hexagon", fontsize=7.5, color=vs.INK_2)
cb.ax.tick_params(labelsize=7)
cb.ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v:g}"))
cb.outline.set_visible(False)

fig.text(0.08, 0.035, "thick line: median transfer among models of similar home accuracy (shaded: middle half)\n"
         "thin curve: accuracy fully kept;  horizontal line: no signal", fontsize=7, color=vs.INK_2)
for a_, letter in ((ax_a, "a"), (ax_b, "b")):
    fig.text(a_.get_position().x0 - 0.055, 0.975, letter, fontsize=10, fontweight="bold", color=vs.INK, va="top")
vs.save(fig, "fig13_transfer")

pd.DataFrame([{"direction": "EUR87_r1 into YRI87", **sa}, {"direction": "YRI87 into EUR358", **sb}]).to_csv(
    OUT / "fig13_summary.tsv", sep="\t", index=False)
print(pd.DataFrame([{"direction": "EUR->YRI", **sa}, {"direction": "YRI->EUR", **sb}]).round(3).to_string(index=False))
