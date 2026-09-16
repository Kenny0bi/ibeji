"""Figure 1: the twin design, drawn one person at a time.

Every dot is one of the 445 individuals. The European arm sits above the mirror
axis and the Yoruba arm below it, drawn as a matched pair.
Filled dots are the 87 Europeans in draw 1 of the downsampled sets (EUR87_r1); hollow dots are the rest.
The weight stems in the middle are the real DNAJB7 models (EUR87_r1 above, YRI87 below).
"""
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

import viz_style as vs

ROOT = Path(__file__).resolve().parents[1]
S = ROOT / "data" / "processed" / "samples"
RAW = ROOT / "data" / "raw"

SP = 0.95          # dot spacing in axis units
NCOL = 10
DOT = 13           # marker area (pt^2)


def ids(name):
    return pd.read_csv(S / f"{name}.txt", sep="\t")["#IID"].tolist()


def arrow(ax, a, b, rad=0.0, color=None):
    ax.add_patch(FancyArrowPatch(a, b, arrowstyle="-|>", mutation_scale=7, color=color or vs.MUTED,
                                 linewidth=0.8, connectionstyle=f"arc3,rad={rad}", shrinkA=0, shrinkB=0))


def people(ax, members, x0, y0, color, highlighted=None, down=False):
    """One dot per person. With `highlighted`, members of that set are filled, the rest hollow."""
    k = np.arange(len(members))
    xs = x0 + (k % NCOL) * SP
    ys = y0 - (k // NCOL) * SP if down else y0 + (k // NCOL) * SP
    if highlighted is None:
        ax.scatter(xs, ys, s=DOT, color=color, linewidths=0, zorder=3)
    else:
        on = np.array([m in highlighted for m in members])
        ax.scatter(xs[on], ys[on], s=DOT, color=color, linewidths=0, zorder=3)
        ax.scatter(xs[~on], ys[~on], s=DOT * 0.8, facecolors="none", edgecolors=color, linewidths=0.6, zorder=3)
    return xs.min(), xs.max(), ys.min(), ys.max()


def weights_glyph(ax, x0, base, color, weights, flip=False, width=14.0, height=3.4, scale=None):
    """Real elastic-net weights for one gene (DNAJB7): one thin stem per SNP the model uses.

    Heights use a square-root scale so small weights stay visible at this size; stems point
    away from the mirror axis for both models. No heads: at 76 SNPs they would pile up.
    """
    n = len(weights)
    xs = x0 + np.arange(n) * (width / max(n - 1, 1))
    scale = scale or np.abs(weights).max()
    ax.plot([xs[0] - 0.4, xs[-1] + 0.4], [base, base], color=vs.GRID, lw=0.8, zorder=1)
    direction = -1 if flip else 1
    nz = weights != 0
    h = direction * np.sqrt(np.abs(weights[nz]) / scale) * height
    ax.vlines(xs[nz], base, base + h, color=color, lw=0.9, zorder=2)


def panel(ax, x0, y0, w, h):
    ax.add_patch(FancyBboxPatch((x0, y0), w, h, boxstyle="round,pad=0,rounding_size=1.0",
                                linewidth=0.8, edgecolor=vs.GRID, facecolor="#faf9f6", zorder=0))


def main():
    vs.apply()
    fig, ax = plt.subplots(figsize=(vs.DOUBLE, 3.05))
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 42)
    ax.set_aspect("equal")
    ax.axis("off")

    psam = pd.read_csv(RAW / "phase3_corrected.psam", sep="\t").set_index("#IID")
    eur, yri, draw = ids("EUR358"), ids("YRI87"), set(ids("EUR87_r1"))
    axis_y = 21

    # mirror axis
    ax.plot([1, 68], [axis_y, axis_y], color=vs.MUTED, lw=0.7, zorder=1)
    ax.text(20.0, axis_y, " 17.3M SNVs common in either sample ", ha="center", va="center",
            fontsize=7.5, color=vs.INK_2, bbox=dict(boxstyle="round,pad=0.25", fc=vs.SURFACE, ec="none"), zorder=5)

    # European arm, above the axis
    x = 2.0
    label_y = axis_y + 2.2 + (NCOL - 1) * SP + 1.4
    for pop in ("FIN", "TSI", "CEU", "GBR"):
        members = [i for i in eur if psam.loc[i, "Population"] == pop]
        _, xmax, _, _ = people(ax, members, x, axis_y + 2.2, vs.EUR, highlighted=draw)
        ax.text((x + xmax) / 2, label_y, f"{pop}  {len(members)}", ha="center", va="bottom", fontsize=7.5, color=vs.INK_2)
        x = xmax + 2.1
    eur_right = x - 2.1
    ax.text(2.0, 40.6, "European ancestry", ha="left", va="top", fontsize=8.5, fontweight="bold", color=vs.INK)
    ax.text(2.0, 37.6, "EUR358 in total", ha="left", va="top", fontsize=7.5, color=vs.INK_2)
    ax.scatter([24.0], [39.9], s=DOT, color=vs.EUR, linewidths=0)
    ax.text(24.9, 39.9, "in random draw 1 of 5 (EUR87)", ha="left", va="center", fontsize=7.5, color=vs.INK_2)
    ax.scatter([24.0], [37.9], s=DOT * 0.8, facecolors="none", edgecolors=vs.EUR, linewidths=0.6)
    ax.text(24.9, 37.9, "not in that draw", ha="left", va="center", fontsize=7.5, color=vs.INK_2)

    # Yoruba arm, mirrored below
    yx0 = 2.0 + 10.5
    xmin, xmax, ymin, _ = people(ax, yri, yx0, axis_y - 2.2, vs.YRI, down=True)
    ax.text(xmax + 1.8, axis_y - 3.4, "all 87 Yoruba with genotypes,", ha="left", va="center", fontsize=7.5, color=vs.INK_2)
    ax.text(xmax + 1.8, axis_y - 5.6, "the same n as each draw", ha="left", va="center", fontsize=7.5, color=vs.INK_2)
    ax.text(2.0, 4.6, "Yoruba, Ibadan", ha="left", va="bottom", fontsize=8.5, fontweight="bold", color=vs.INK)
    ax.text(2.0, 2.6, "YRI87", ha="left", va="bottom", fontsize=7.5, color=vs.INK_2)

    # twin models
    arrow(ax, (eur_right + 0.8, 28), (48.5, 28.5))
    arrow(ax, (xmax + 1.5, 11.0), (48.5, 12.8))
    detail = pd.read_csv(ROOT / "results" / "figdata" / "locus" / "DNAJB7_snps_detail.tsv", sep="\t")
    scale = max(detail["w_EUR"].abs().max(), detail["w_YRI"].abs().max())
    weights_glyph(ax, 50.0, 28.5, vs.EUR, detail["w_EUR"].to_numpy(), scale=scale)
    weights_glyph(ax, 50.0, 13.5, vs.YRI, detail["w_YRI"].to_numpy(), flip=True, scale=scale)
    ax.text(57, 25.2, f"DNAJB7, {len(detail)} SNPs", ha="center", fontsize=7, color=vs.INK_2)
    ax.text(57, 36.6, "Elastic net per gene", ha="center", fontsize=8, fontweight="bold", color=vs.INK)
    ax.text(57, 33.9, r"weights $w_{\mathrm{EUR}}$ on cis-SNPs", ha="center", fontsize=7.5, color=vs.INK_2)
    ax.text(57, 7.6, r"weights $w_{\mathrm{YRI}}$ on cis-SNPs", ha="center", fontsize=7.5, color=vs.INK_2)
    ax.text(57, 4.6, "α = 0.5, nested cross-validation", ha="center", fontsize=7.5, color=vs.INK_2)
    ax.text(58.5, axis_y, " ~20,800 genes per set ", ha="center", va="center", fontsize=7.5, color=vs.INK_2,
            bbox=dict(boxstyle="round,pad=0.25", fc=vs.SURFACE, ec="none"), zorder=5)

    # analyses
    panel(ax, 71, 22.8, 28, 18.6)
    ax.text(85, 40.5, "Variance decomposition", ha="center", va="top", fontsize=8, fontweight="bold", color=vs.INK)
    ax.text(85, 36.9, r"$V = w^{\top} D^{1/2}\, R\, D^{1/2}\, w$", ha="center", va="center", fontsize=8.5, color=vs.INK)
    ax.text(85, 34.2, r"$\Delta = \log V_{\mathrm{YRI}} - \log V_{\mathrm{EUR}}$", ha="center", va="center", fontsize=8, color=vs.INK)
    ax.text(85, 31.6, r"$\Delta = \phi_w + \phi_D + \phi_R$", ha="center", va="center", fontsize=8, color=vs.INK)
    segs = [("weights", r"$w$", vs.PHI_W), ("allele freq.", r"$D$", vs.PHI_D), ("LD", r"$R$", vs.PHI_R)]
    bx = 73.0
    for label, sym, color in segs:
        ax.add_patch(plt.Rectangle((bx, 28.2), 7.6, 1.4, color=color, lw=0, zorder=2))
        ax.text(bx + 3.8, 27.7, sym, ha="center", va="top", fontsize=8, color=vs.INK_2)
        bx += 8.0
    ax.text(85, 24.3, "w weights, D frequency, R LD", ha="center", va="center", fontsize=7, color=vs.INK_2)

    panel(ax, 71, 0.4, 28, 18.8)
    ax.text(85, 17.4, "Summary-statistic TWAS", ha="center", va="top", fontsize=8, fontweight="bold", color=vs.INK)
    ax.text(85, 14.4, "European LD reference", ha="center", va="center", fontsize=7.5, color=vs.INK_2)
    traits = ["ASD", "SCZ", "BIP", "MDD", "PTSD"]
    cx = 71.8
    for t in traits:
        ax.add_patch(FancyBboxPatch((cx, 9.0), 5.0, 3.4, boxstyle="round,pad=0,rounding_size=0.8",
                                    lw=0.8, edgecolor=vs.INK_2, facecolor=vs.SURFACE, zorder=2))
        ax.text(cx + 2.5, 10.7, t, ha="center", va="center", fontsize=7, color=vs.INK, zorder=3)
        cx += 5.3
    ax.text(85, 7.2, "+ GWAS coverage of each model", ha="center", va="center", fontsize=7, color=vs.INK_2)
    ax.text(85, 4.6, "autism + 4 psychiatric disorders", ha="center", va="center", fontsize=7, color=vs.INK_2)
    ax.text(85, 2.0, "with large public GWAS", ha="center", va="center", fontsize=7, color=vs.INK_2)

    arrow(ax, (65.5, 30.5), (70.6, 32))
    arrow(ax, (65.5, 11.5), (70.6, 10))
    arrow(ax, (65.5, 26.0), (70.6, 14.5), rad=-0.12)
    arrow(ax, (65.5, 16.0), (70.6, 27.5), rad=0.12)

    vs.save(fig, "fig01_design")


if __name__ == "__main__":
    main()
