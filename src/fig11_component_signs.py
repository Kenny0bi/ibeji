"""Figure 11: how large is each ingredient, and which way does it push?

Usage: python3 src/fig11_component_signs.py

(a) The signed distribution of each component of the gap between the two models, drawn as a
    density ridge on a common axis. Each curve is scaled to its own height so the narrow
    components stay visible; width on the horizontal axis is the quantity being compared.
    Values outside the axis are counted at the edge rather than drawn, so no artificial
    mode appears there.
(b) The share of genes for which each component is positive, with a 95% interval. A component
    that only reflects training noise should sit at half; one that reflects a systematic
    difference between the populations should not.
Reads results/decomposition/EUR87_r1_vs_YRI87.tsv.
"""
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import viz_style as vs

ROOT = Path(__file__).resolve().parents[1]
d = pd.read_csv(ROOT / "results/decomposition/EUR87_r1_vs_YRI87.tsv", sep="\t")
d = d.dropna(subset=["phi_w", "phi_D", "phi_R"])

COMPS = [("phi_w", vs.PHI_W, "weights", r"$\phi_w$"),
         ("phi_D", vs.PHI_D, "allele frequency", r"$\phi_D$"),
         ("phi_R", vs.PHI_R, "linkage disequilibrium", r"$\phi_R$")]
XLIM = 4.0
bins = np.linspace(-XLIM, XLIM, 161)
centers = 0.5 * (bins[:-1] + bins[1:])


def smooth(y, k=7):
    w = np.ones(k) / k
    return np.convolve(np.pad(y, k // 2, mode="edge"), w, mode="valid")


vs.apply()
fig = plt.figure(figsize=(vs.DOUBLE, 3.5))
ax = fig.add_axes([0.20, 0.175, 0.48, 0.63])
bx = fig.add_axes([0.775, 0.175, 0.195, 0.63])

rows, gap = [], 1.0
ax.axvline(0, color=vs.INK_2, lw=0.7, zorder=2)
for i, (col, color, name, sym) in enumerate(COMPS):
    v = d[col].values
    base = (len(COMPS) - 1 - i) * gap
    # Values outside the axis are counted and annotated, never piled into the end bins:
    # an earlier version did that and produced two walls that read as modes.
    inside = v[np.abs(v) <= XLIM]
    h, _ = np.histogram(inside, bins=bins)
    y = smooth(h.astype(float))
    y = y / y.max() * (gap * 0.86)
    ax.fill_between(centers, base, base + y, color=color, alpha=0.5, lw=0, zorder=3)
    ax.plot(centers, base + y, color=color, lw=1.2, zorder=4)
    med = float(np.median(v))
    ax.plot([med, med], [base, base + gap * 0.5], color=vs.INK, lw=0.9, zorder=5)
    n_lo, n_hi = int((v < -XLIM).sum()), int((v > XLIM).sum())
    # Raised clear of the left-hand label block, which reaches the axis edge at this scale.
    for n_out, xpos, ha in ((n_lo, -XLIM, "left"), (n_hi, XLIM, "right")):
        if n_out:
            ax.text(xpos, base + gap * 0.55, f"{n_out} beyond", fontsize=7, color=vs.INK_2, ha=ha, va="bottom")
    # identity comes from the filled ridge immediately to the right of this label
    ax.text(-XLIM * 1.30, base + gap * 0.30, f"{sym}  {name}", fontsize=7.5, color=vs.INK_2, ha="left", va="center")
    ax.text(-XLIM * 1.30, base + gap * 0.11, f"median {med:+.2f}", fontsize=7, color=vs.INK, ha="left", va="center")
    out = int((np.abs(v) > XLIM).sum())
    p = float((v > 0).mean()); se = float(np.sqrt(p * (1 - p) / len(v)))
    rows.append(dict(component=name, sym=sym, color=color, n=len(v), median=med, p_pos=p, se=se,
                     lo=p - 1.96 * se, hi=p + 1.96 * se, beyond=out,
                     iqr=float(np.percentile(v, 75) - np.percentile(v, 25))))

ax.set_xlim(-XLIM, XLIM)
ax.set_ylim(-0.08, len(COMPS) * gap)
ax.set_yticks([])
ax.spines["left"].set_visible(False)
ax.set_xlabel("component of $\\log V_{\\mathrm{YRI}} - \\log V_{\\mathrm{EUR}}$")
ax.set_title("Size and direction of each ingredient", fontsize=8.5, color=vs.INK, loc="left", pad=30)
beyond = sum(r["beyond"] for r in rows)
ax.text(0.0, 1.03, f"{len(d)} genes;  each curve scaled to its own height\n"
        f"black tick: the median;  {beyond} values fall outside the axis and are counted, not drawn",
        transform=ax.transAxes, fontsize=7.5, color=vs.INK_2, va="bottom")

for i, r in enumerate(rows):
    y = (len(rows) - 1 - i)
    bx.plot([r["lo"], r["hi"]], [y, y], color=r["color"], lw=2.0, solid_capstyle="round", zorder=3)
    bx.scatter([r["p_pos"]], [y], s=26, color=r["color"], edgecolors=vs.SURFACE, linewidths=0.8, zorder=4)
    bx.text(r["p_pos"], y + 0.22, f"{r['p_pos']:.0%}", fontsize=7.5, color=vs.INK, ha="center", va="bottom")
bx.axvline(0.5, color=vs.INK_2, lw=0.7, ls=(0, (3, 2)), zorder=2)
bx.text(0.507, -0.62, "half", fontsize=7.5, color=vs.INK_2, ha="left", va="bottom")
bx.set_xlim(0.38, 0.65)
bx.set_ylim(-0.75, len(rows) - 0.25)
bx.set_yticks([])
bx.spines["left"].set_visible(False)
bx.set_xticks([0.4, 0.5, 0.6])
bx.set_xticklabels(["40%", "50%", "60%"])
bx.xaxis.ibeji_value_of = lambda v: v * 100
bx.set_xlabel("genes where the\ncomponent is positive")
bx.set_title("Which way?", fontsize=8.5, color=vs.INK, loc="left", pad=30)
bx.text(0.0, 1.03, "bar: 95% interval", transform=bx.transAxes, fontsize=7.5, color=vs.INK_2, va="bottom")

for a_, letter in ((ax, "a"), (bx, "b")):
    fig.text(a_.get_position().x0 - 0.075, 0.985, letter, fontsize=9, fontweight="bold", color=vs.INK, va="top")
vs.save(fig, "fig11_component_signs")

out = pd.DataFrame(rows).drop(columns=["color"])
out.to_csv(ROOT / "results/figdata/fig11_summary.tsv", sep="\t", index=False)
print(out.round(4).to_string(index=False))
