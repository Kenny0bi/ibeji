"""Figure 7: what sample size alone buys, measured against what a different draw already buys.

Usage: python3 src/fig07_sample_size.py [large_set] [small_set] [null_a] [null_b]
Default: EUR358 EUR87_r1 EUR87_r1 EUR87_r2

Training on more individuals of the SAME ancestry should improve prediction. The question is
how much of any improvement is real. Two independent draws of 87 Europeans already differ from
each other, so the honest reference for the large-versus-small comparison is the difference
between two equal-sized draws, not zero.

(a) Cross-validated accuracy at the small size against the large size, per gene.
(b) The paired per-gene difference for the real comparison, drawn against the equal-size null.
    The separation between the two curves is what sample size actually buys.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, LogNorm
from scipy import stats

import viz_style as vs

ROOT = Path(__file__).resolve().parents[1]
R = ROOT / "results" / "models"
LARGE, SMALL, NA, NB = (sys.argv[1:5] + ["EUR358", "EUR87_r1", "EUR87_r1", "EUR87_r2"][len(sys.argv) - 1:])[:4]
NICE = {"EUR358": "358 Europeans", "EUR87_r1": "87 Europeans (draw 1)", "EUR87_r2": "87 Europeans (draw 2)"}


def load(s):
    files = sorted((R / s).glob("chr*.summary.tsv"))
    if len(files) < 22:
        raise SystemExit(f"{s} has only {len(files)}/22 chromosomes; not ready")
    d = pd.concat([pd.read_csv(f, sep="\t") for f in files], ignore_index=True)
    d["usable"] = (d.cv_r2 > 0.01) & (d.cv_pval < 0.05) & (d.n_model > 0)
    return d[["gene", "cv_r2", "usable"]].rename(columns={"cv_r2": s, "usable": f"u_{s}"})


def paired(big, small):
    m = load(big).merge(load(small), on="gene")
    keep = m[f"u_{big}"] & m[f"u_{small}"]
    return m[keep], m[keep][big] - m[keep][small]


real, d_real = paired(LARGE, SMALL)
null, d_null = paired(NB, NA)
MACHINERY_TEST = LARGE != "EUR358"
if MACHINERY_TEST:
    print(f"*** MACHINERY TEST: large set is {LARGE}, not EUR358. Output is NOT a result. ***")

vs.apply()
fig = plt.figure(figsize=(vs.DOUBLE, 3.3))
ax = fig.add_axes([0.085, 0.215, 0.30, 0.54])
bx = fig.add_axes([0.545, 0.215, 0.41, 0.54])

MAXV = 0.8
cmap = LinearSegmentedColormap.from_list("d", ["#eef0f3", "#b9bdc6", "#6f7581", "#2f3440"])
hb = ax.hexbin(np.sqrt(real[SMALL].clip(0, MAXV)), np.sqrt(real[LARGE].clip(0, MAXV)),
               gridsize=26, extent=(0, np.sqrt(MAXV), 0, np.sqrt(MAXV)), cmap=cmap,
               norm=LogNorm(vmin=1), mincnt=1, linewidths=0, rasterized=True)
xs = np.linspace(0, np.sqrt(MAXV), 200)
ax.plot(xs, xs, color=vs.MUTED, lw=0.7, zorder=3)
# Thinned: on a square-root axis 0.05 and 0.1 sit at 0.224 and 0.316, and their labels collide.
ticks = [0, 0.1, 0.2, 0.4, 0.6]
ax.set_xticks(np.sqrt(ticks)); ax.set_xticklabels([f"{t:g}" for t in ticks])
ax.set_yticks(np.sqrt(ticks)); ax.set_yticklabels([f"{t:g}" for t in ticks])
ax.xaxis.ibeji_value_of = np.square
ax.yaxis.ibeji_value_of = np.square
ax.set_xlim(0, np.sqrt(MAXV)); ax.set_ylim(0, np.sqrt(MAXV)); ax.set_aspect("equal")
ax.set_xlabel(f"cross-validated R², {NICE.get(SMALL, SMALL)}")
ax.set_ylabel(f"cross-validated R², {NICE.get(LARGE, LARGE)}")
ax.set_title("Per gene, small against large", fontsize=8.5, color=vs.INK, loc="left", pad=30)
ax.text(0.0, 1.03, f"{len(real):,} genes usable in both\nline: equal accuracy", transform=ax.transAxes,
        fontsize=7.5, color=vs.INK_2, va="bottom")

# Density key. Without it the grey ramp is an unexplained colour, which the manual checklist
# forbids, and Figure 10's equivalent hexbin panel carries the same key. The label sits beside
# the bar rather than under it, because a stacked label runs off the bottom of a 3.3 in figure.
cax = fig.add_axes([0.105, 0.052, 0.115, 0.018])
cb = fig.colorbar(hb, cax=cax, orientation="horizontal")
cb.outline.set_visible(False)
cb.ax.tick_params(labelsize=7, length=2, colors=vs.INK_2, pad=1.5)
# Plain "1" and "10", not 10^0 and 10^1. Figure 10's colorbar shows the same quantity and
# uses this formatter; two different labellings of "genes per hexagon" in one paper is a
# reader's problem, not a style preference.
cb.ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v:g}"))
fig.text(0.232, 0.061, "genes per hexagon", fontsize=7.5, color=vs.INK_2, va="center")

bins = np.linspace(-0.35, 0.35, 121)
centres = 0.5 * (bins[:-1] + bins[1:])
def smooth(y, k=7):
    return np.convolve(np.pad(y.astype(float), k // 2, mode="edge"), np.ones(k) / k, mode="valid")
for vals, color, label in ((d_null, vs.MUTED, f"two draws of 87\n(what noise alone gives)"),
                           (d_real, vs.EUR, f"87 against {NICE.get(LARGE, LARGE).split()[0]}\n(sample size)")):
    h, _ = np.histogram(np.clip(vals, bins[0], bins[-1]), bins=bins)
    y = smooth(h) / smooth(h).max()
    bx.fill_between(centres, 0, y, color=color, alpha=0.30, lw=0, zorder=2)
    bx.plot(centres, y, color=color, lw=1.6, zorder=3)
    bx.plot([vals.median()] * 2, [0, 1.06], color=color, lw=0.9, ls=(0, (3, 2)), zorder=4)
bx.axvline(0, color=vs.INK_2, lw=0.7, zorder=1)
# Headroom for the legend. The curves are normalised to peak at 1.0 and the median rules are
# drawn to 1.06, so the legend rows have to clear both. At the old limit of 1.24 the second
# row sat at data y 0.968 to 1.038 and the blue curve ran straight through it. Measured, not
# eyeballed: the curve reaches 1.000 inside that row's x range.
bx.set_xlim(bins[0], bins[-1]); bx.set_ylim(0, 1.35)
bx.set_yticks([])
bx.spines["left"].set_visible(False)
bx.set_xlabel("change in cross-validated R² per gene")
bx.set_title("Against what a different draw already gives", fontsize=8.5, color=vs.INK, loc="left", pad=30)
bx.text(0.0, 1.03, f"dashed: medians, {d_null.median():+.3f} null and {d_real.median():+.3f} real",
        transform=bx.transAxes, fontsize=7.5, color=vs.INK_2, va="bottom")
for i, (color, label) in enumerate(((vs.MUTED, "two draws of 87, noise alone"),
                                    (vs.EUR, f"87 against {NICE.get(LARGE, LARGE)}"))):
    y = 0.93 - 0.075 * i
    bx.add_patch(plt.Rectangle((0.035, y - 0.035), 0.030, 0.055, transform=bx.transAxes, color=color,
                               alpha=0.65, lw=0, zorder=5))
    bx.text(0.082, y - 0.008, label, transform=bx.transAxes, fontsize=7.5, color=vs.INK_2, va="center", zorder=5)

for a_, letter in ((ax, "a"), (bx, "b")):
    fig.text(max(a_.get_position().x0 - 0.062, 0.008), 0.975, letter, fontsize=9, fontweight="bold",
             color=vs.INK, va="top")
vs.save(fig, "fig07_sample_size" + ("_MACHINERYTEST" if MACHINERY_TEST else ""))

w_real = stats.wilcoxon(real[LARGE], real[SMALL])
print(f"real  {SMALL} -> {LARGE}: n={len(real):,} median {d_real.median():+.4f} improving {(d_real>0).mean():.1%} p={w_real.pvalue:.2g}")
print(f"null  {NA} -> {NB}: n={len(null):,} median {d_null.median():+.4f} improving {(d_null>0).mean():.1%}")
