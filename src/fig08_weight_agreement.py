"""Figure 8: when both models pick the same SNP, do they agree about it?

Usage: python3 src/fig08_weight_agreement.py

(a) Every SNP selected by BOTH the European and the Yoruba model for the same gene, one dot,
    placed by the weight each model gave it. The few SNPs the two models push in opposite
    directions are ringed. Whether one model systematically gave larger weights is reported
    as a number rather than drawn, since the split is close to even.
(b) The same agreement measured within bins of how large the weight is, since a weight near
    zero carries little information about either model.
Reads results/figdata/fig08_shared_weights.tsv.
"""
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import viz_style as vs

ROOT = Path(__file__).resolve().parents[1]
d = pd.read_csv(ROOT / "results/figdata/fig08_shared_weights.tsv", sep="\t")
d["mag"] = np.sqrt(d.weight_E.abs() * d.weight_Y.abs())
d["concordant"] = np.sign(d.weight_E) == np.sign(d.weight_Y)
d["eur_heavier"] = d.weight_E.abs() > d.weight_Y.abs()
r_all = np.corrcoef(d.weight_E, d.weight_Y)[0, 1]

vs.apply()
fig = plt.figure(figsize=(vs.DOUBLE, 3.5))
ax = fig.add_axes([0.075, 0.155, 0.42, 0.66])
bx = fig.add_axes([0.62, 0.155, 0.35, 0.66])

lim = float(np.ceil(max(d.weight_E.abs().max(), d.weight_Y.abs().max()) * 10) / 10)
ax.axhline(0, color=vs.GRID, lw=0.6, zorder=1)
ax.axvline(0, color=vs.GRID, lw=0.6, zorder=1)
ax.plot([-lim, lim], [-lim, lim], color=vs.MUTED, lw=0.7, zorder=2)
ax.scatter(d.weight_E, d.weight_Y, s=5, color=vs.PHI_W, alpha=0.55, linewidths=0, zorder=3)
disc = d[~d.concordant]
ax.scatter(disc.weight_E, disc.weight_Y, s=22, facecolors="none", edgecolors=vs.INK, linewidths=0.7, zorder=4)
ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
ax.set_xlabel("weight in the European model")
ax.set_ylabel("weight in the Yoruba model")
ax.set_aspect("equal")

# quadrant counts: the two agreeing quadrants carry almost everything
q = {"++": ((d.weight_E > 0) & (d.weight_Y > 0)).sum(), "--": ((d.weight_E < 0) & (d.weight_Y < 0)).sum(),
     "+-": ((d.weight_E > 0) & (d.weight_Y < 0)).sum(), "-+": ((d.weight_E < 0) & (d.weight_Y > 0)).sum()}
# The two agreeing quadrants are labelled off the diagonal so the reference line never
# crosses their text; the two disagreeing quadrants are labelled in their far corners.
pad = lim * 0.94
off = lim * 0.42
for (key, xy, ha, va) in (("++", (pad, off), "right", "center"), ("--", (-pad, -off), "left", "center"),
                          ("+-", (pad, -pad), "right", "bottom"), ("-+", (-pad, pad), "left", "top")):
    ax.text(*xy, f"{q[key]:,}", ha=ha, va=va, fontsize=8, color=vs.INK_2)

ax.set_title("Weights on the SNPs both models selected", fontsize=8.5, color=vs.INK, loc="left", pad=26)
ax.text(0.0, 1.04, f"{len(d):,} SNP-gene pairs;  r = {r_all:.2f};  {d.concordant.mean():.1%} agree in sign\n"
        f"ringed: the {len(disc)} pushed in opposite directions", transform=ax.transAxes, fontsize=7.5,
        color=vs.INK_2, va="bottom")
print(f"European weight was the larger one for {d.eur_heavier.mean():.1%} of shared SNPs")
# (b) agreement within bins of weight magnitude
d["bin"] = pd.qcut(d["mag"], 6, labels=False, duplicates="drop")
g = d.groupby("bin").apply(lambda s: pd.Series({
    "mag": s["mag"].median(), "n": len(s), "sign": s["concordant"].mean(),
    "r": np.corrcoef(s.weight_E, s.weight_Y)[0, 1]}))
bx.plot(g["mag"], g["r"], color=vs.PHI_W, lw=2.0, zorder=3, solid_capstyle="round")
bx.scatter(g["mag"], g["r"], s=16, color=vs.PHI_W, edgecolors=vs.SURFACE, linewidths=0.7, zorder=4)
bx.plot(g["mag"], g["sign"], color=vs.COVERAGE, lw=2.0, zorder=3, solid_capstyle="round")
bx.scatter(g["mag"], g["sign"], s=16, color=vs.COVERAGE, edgecolors=vs.SURFACE, linewidths=0.7, zorder=4)
bx.text(g["mag"].iloc[-1], g["sign"].iloc[-1] + 0.035, "agree in sign", ha="right", fontsize=7.5, color=vs.COVERAGE)
# Placed in the open band between the two curves: below the sign-agreement line, well above
# the correlation line at this x. Sitting it just under its own curve put the curve through it.
bx.text(g["mag"].iloc[1], 0.73, "correlation of the\ntwo weights", ha="left", va="center",
        fontsize=7.5, color=vs.PHI_W)
bx.set_xscale("log")
bx.set_ylim(0, 1.05)
bx.set_xlabel("size of the weight (geometric mean of |w|)")
bx.set_ylabel("agreement")
bx.grid(True, axis="y"); bx.set_axisbelow(True)
bx.set_title("Larger weights agree more", fontsize=8.5, color=vs.INK, loc="left", pad=26)
bx.text(0.0, 1.04, f"six equal groups of {int(g['n'].mean()):,} SNP-gene pairs", transform=bx.transAxes,
        fontsize=7.5, color=vs.INK_2, va="bottom")
for a_, letter in ((ax, "a"), (bx, "b")):
    fig.text(a_.get_position().x0 - 0.055, 0.985, letter, fontsize=9, fontweight="bold", color=vs.INK, va="top")
vs.save(fig, "fig08_weight_agreement")
g.round(3).to_csv(ROOT / "results/figdata/fig08_summary.tsv", sep="\t")
print(g.round(3).to_string())
