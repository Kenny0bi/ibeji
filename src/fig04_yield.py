"""Figure 4: how many genes get a usable model, and how much does that depend on the draw?

Usage: python3 src/fig04_yield.py

(a) The funnel from genes attempted to usable models, per training set. Most genes that clear
    cross-validation still end with an empty elastic-net fit at this sample size, and that is
    where the yield is actually lost.
(b) Usable models per set. The European draws of 87 are shown as a group so the spread between
    two random draws of one population can be read against the gap to the Yoruba sample.

Any training set with all 22 chromosomes is included automatically, so this figure extends
itself as the remaining sets finish.
"""
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import viz_style as vs

ROOT = Path(__file__).resolve().parents[1]
R = ROOT / "results" / "models"

ORDER = ["YRI87", "EUR87_r1", "EUR87_r2", "EUR87_r3", "EUR87_r4", "EUR87_r5", "EUR358"]
LABEL = {"YRI87": "Yoruba, 87", "EUR358": "European, 358"}
for i in range(1, 6):
    LABEL[f"EUR87_r{i}"] = f"European, 87 (draw {i})"


def load(set_name):
    files = sorted((R / set_name).glob("chr*.summary.tsv"))
    if len(files) < 22:
        return None
    d = pd.concat([pd.read_csv(f, sep="\t") for f in files], ignore_index=True)
    passed = d[(d.cv_r2 > 0.01) & (d.cv_pval < 0.05)]
    usable = passed[passed.n_model > 0]
    return dict(set=set_name, attempted=len(d), passed_cv=len(passed), usable=len(usable),
                empty_final=len(passed) - len(usable), median_r2=float(usable.cv_r2.median()),
                median_snps=float(usable.n_model.median()))


rows = [r for r in (load(s) for s in ORDER) if r is not None]
d = pd.DataFrame(rows)
print(d.to_string(index=False))

vs.apply()
fig = plt.figure(figsize=(vs.DOUBLE, 3.7))
ax = fig.add_axes([0.20, 0.215, 0.40, 0.46])
bx = fig.add_axes([0.715, 0.215, 0.255, 0.46])

STAGES = [("attempted", vs.MUTED, "genes attempted"),
          ("passed_cv", vs.SEQ[7], "cleared cross-validation"),
          ("usable", vs.SEQ[11], "usable model")]
T = ax.get_yaxis_transform()
y = 0.0
for _, r in d.iterrows():
    for k, (col, color, _) in enumerate(STAGES):
        ax.barh([y - k * 0.26], [r[col]], height=0.24, color=color, linewidth=0, zorder=3)
        ax.text(r[col] + 300, y - k * 0.26, f"{int(r[col]):,}", va="center", fontsize=7, color=vs.INK_2)
    ax.text(-0.015, y - 0.26, LABEL[r["set"]], transform=T, ha="right", va="center",
            fontsize=7.5, color=vs.INK)
    y -= 1.15
ax.set_ylim(y + 0.55, 0.42)
ax.set_yticks([])
ax.spines["left"].set_visible(False)
ax.set_xlim(0, 24500)
ax.set_xticks([0, 5000, 10000, 15000, 20000])
ax.set_xticklabels(["0", "5k", "10k", "15k", "20k"])
ax.xaxis.ibeji_value_of = lambda v: v
ax.set_xlabel("genes")
ax.grid(True, axis="x"); ax.set_axisbelow(True)
ax.set_title("Where the genes are lost", fontsize=8.5, color=vs.INK, loc="left", pad=62)
LEG = [(0.0, 1.150), (0.0, 1.060), (0.45, 1.060)]
for (x, yy), (_, color, label) in zip(LEG, STAGES):
    ax.add_patch(plt.Rectangle((x, yy - 0.025), 0.024, 0.05, transform=ax.transAxes, color=color,
                               lw=0, clip_on=False))
    ax.text(x + 0.034, yy, label, transform=ax.transAxes, fontsize=7.5, color=vs.INK_2, va="center")

eur87 = d[d.set.str.startswith("EUR87")]
others = d[~d.set.str.startswith("EUR87")]
xs = []
ticklabels = []
if len(eur87):
    lo, hi = eur87.usable.min(), eur87.usable.max()
    # The axis is deliberately zoomed, because the spread between draws is the subject. The
    # shaded band carries that spread across the panel so the Yoruba point can be read against
    # it directly, rather than leaving a truncated axis to exaggerate a small difference.
    bx.axhspan(lo, hi, color=vs.EUR, alpha=0.11, lw=0, zorder=1)
    bx.plot([0, 0], [lo, hi], color=vs.EUR, lw=2.2, solid_capstyle="round", zorder=3)
    bx.scatter([0] * len(eur87), eur87.usable, s=30, color=vs.EUR, edgecolors=vs.SURFACE,
               linewidths=0.8, zorder=4)
    xs.append(0); ticklabels = [f"European, 87\n{len(eur87)} draws"]
for j, (_, r) in enumerate(others.iterrows(), start=1):
    color = vs.YRI if r["set"] == "YRI87" else vs.PHI_W
    bx.scatter([j], [r.usable], s=30, color=color, edgecolors=vs.SURFACE, linewidths=0.8, zorder=4)
    xs.append(j); ticklabels.append(LABEL[r["set"]].replace(", ", "\n"))
bx.set_xlim(-0.6, max(xs) + 0.6)
bx.set_xticks(xs); bx.set_xticklabels(ticklabels, fontsize=7.5)
bx.tick_params(axis="x", length=0)
bx.set_ylabel("usable models")
bx.grid(True, axis="y"); bx.set_axisbelow(True)
bx.set_title("Yield at equal sample size", fontsize=8.5, color=vs.INK, loc="left", pad=62)
if len(eur87) > 1:
    spread = int(eur87.usable.max() - eur87.usable.min())
    yri = d.loc[d.set == "YRI87", "usable"]
    inside = bool(len(yri) and eur87.usable.min() <= yri.iloc[0] <= eur87.usable.max())
    # Kept to short lines: panel b is 1.83 in wide and a single long line overran the canvas.
    note = f"shaded: spread of the\n{len(eur87)} European draws, {spread} genes"
    if inside:
        note += "\nthe Yoruba set is inside it"
    bx.text(0.0, 1.055, note, transform=bx.transAxes, fontsize=7.5, color=vs.INK_2, va="bottom")

for a_, letter in ((ax, "a"), (bx, "b")):
    fig.text(max(a_.get_position().x0 - 0.062, 0.008), 0.975, letter, fontsize=9, fontweight="bold",
             color=vs.INK, va="top")
vs.save(fig, "fig04_yield")
d.to_csv(ROOT / "results" / "figdata" / "fig04_yield.tsv", sep="\t", index=False)
