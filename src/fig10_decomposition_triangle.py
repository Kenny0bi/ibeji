"""Figure 10: where each gene's disagreement comes from, against the noise floor.

Usage: python3 src/fig10_decomposition_triangle.py [cross_file] [floor_file]

Every gene is placed in a triangle by the shares of its disagreement due to weights (φ_w),
allele frequencies (φ_D) and LD (φ_R), using absolute values normalized to sum to 1:
  corner "weights" = all of the gap from weights, and so on.
(a) European model vs Yoruba model (cross-ancestry).
(b) European draw 1 vs European draw 2 (same ancestry). Here φ_D = φ_R = 0 by construction,
    so every gene sits on the weights corner: the noise floor for the weight share.
Panel b therefore does not show shape; its size is summarized as a bar under each triangle,
the median |φ_w|, so the cross-ancestry weight disagreement can be read against what two
European draws already produce. Degenerate genes (V = 0 in some mix) are excluded.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, LogNorm
from matplotlib.patches import Polygon

import viz_style as vs

ROOT = Path(__file__).resolve().parents[1]
D = ROOT / "results" / "decomposition"
OUT = ROOT / "results" / "figdata"
cross_file = Path(sys.argv[1]) if len(sys.argv) > 1 else D / "EUR87_r1_vs_YRI87.tsv"
floor_file = Path(sys.argv[2]) if len(sys.argv) > 2 else D / "EUR87_r1_vs_EUR87_r2.tsv"
H = np.sqrt(3) / 2


def shares(d):
    a = d[["phi_w", "phi_D", "phi_R"]].abs().to_numpy()
    s = a.sum(axis=1, keepdims=True)
    return a / np.where(s == 0, 1, s)


def to_xy(sh):
    # corners: weights (0, 0), allele frequency (1, 0), LD (0.5, H)
    return sh[:, 1] + sh[:, 2] / 2, sh[:, 2] * H


def triangle(ax, sh, title):
    tri = Polygon([(0, 0), (1, 0), (0.5, H)], closed=True, facecolor="#faf9f6", edgecolor=vs.MUTED, lw=0.7, zorder=0)
    ax.add_patch(tri)
    for t in (0.2, 0.4, 0.6, 0.8):
        for p, q in (((t, 0), (t / 2, t * H)), ((t, 0), (0.5 + t / 2, (1 - t) * H)), ((t / 2, t * H), (1 - t / 2, t * H))):
            ax.plot([p[0], q[0]], [p[1], q[1]], color=vs.GRID, lw=0.5, zorder=1)
    x, y = to_xy(sh)
    cmap = LinearSegmentedColormap.from_list("dens", ["#e8f0fb"] + vs.SEQ[3:])
    hb = ax.hexbin(x, y, gridsize=24, extent=(0, 1, 0, H), cmap=cmap, norm=LogNorm(vmin=1), mincnt=1,
                   linewidths=0.2, edgecolors="#faf9f6", zorder=2)
    hb.set_clip_path(tri)  # hexagons never spill past the triangle's edges
    med = np.median(sh, axis=0)
    mx, my = to_xy(med[None, :])
    ax.scatter(mx, my, s=34, facecolors=vs.SURFACE, edgecolors=vs.INK, linewidths=1.1, zorder=4)
    ax.annotate("median gene", (mx[0], my[0]), xytext=(8, 6), textcoords="offset points", fontsize=7, color=vs.INK,
                bbox=dict(boxstyle="round,pad=0.2", fc=vs.SURFACE, ec="none", alpha=0.9), zorder=5)
    for (cx, cy), color, label, ha, dy in (((0, 0), vs.PHI_W, "weights", "right", -0.045),
                                          ((1, 0), vs.PHI_D, "allele frequency", "left", -0.045),
                                          ((0.5, H), vs.PHI_R, "LD", "center", 0.035)):
        ax.scatter([cx], [cy], s=26, color=color, edgecolors=vs.SURFACE, linewidths=0.8, zorder=6, clip_on=False)
        ax.text(cx, cy + dy, label, ha=ha, va="top" if dy < 0 else "bottom", fontsize=7.5, color=vs.INK)
    ax.set_xlim(-0.12, 1.12)
    ax.set_ylim(-0.12, H + 0.1)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(title, fontsize=8, color=vs.INK, loc="left", pad=12)
    ax.text(0.0, 1.0, "each gene placed by its shares of $|\\phi_w|$, $|\\phi_D|$ and $|\\phi_R|$", transform=ax.transAxes,
            fontsize=7, color=vs.INK_2, va="bottom")
    return hb, med


cross = pd.read_csv(cross_file, sep="\t")
cross = cross[~cross["degenerate"]]
floor = pd.read_csv(floor_file, sep="\t") if floor_file.exists() else None
if floor is not None:
    floor = floor[~floor["degenerate"]]
floor_is_draft = floor is not None and len(floor) < 200

vs.apply()
fig = plt.figure(figsize=(vs.DOUBLE, 3.7))
ax_a = fig.add_axes([0.02, 0.22, 0.46, 0.66])
ax_bar = fig.add_axes([0.66, 0.2, 0.31, 0.52])
cax = fig.add_axes([0.1, 0.13, 0.3, 0.016])

sh_cross = shares(cross)
hb, med = triangle(ax_a, sh_cross, f"European vs Yoruba model: {len(cross):,} genes")
cb = fig.colorbar(hb, cax=cax, orientation="horizontal")
cb.set_label("genes per hexagon", fontsize=7.5, color=vs.INK_2)
cb.ax.tick_params(labelsize=7)
cb.ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v:g}"))
cb.outline.set_visible(False)

# (b) sizes: cross-ancestry component medians next to the within-ancestry floor
rows = [(r"$\phi_w$  weights", cross["phi_w"].abs().median(), floor["phi_w"].abs().median() if floor is not None else np.nan, vs.PHI_W),
        (r"$\phi_D$  allele frequency", cross["phi_D"].abs().median(), 0.0, vs.PHI_D),
        (r"$\phi_R$  LD", cross["phi_R"].abs().median(), 0.0, vs.PHI_R)]
ypos = np.arange(len(rows))[::-1]
for yy, (label, c_val, f_val, color) in zip(ypos, rows):
    ax_bar.barh(yy + 0.17, c_val, height=0.3, color=color, lw=0)
    ax_bar.barh(yy - 0.17, f_val if np.isfinite(f_val) else 0, height=0.3, color=color, alpha=0.35, lw=0)
    ax_bar.text(c_val + 0.01, yy + 0.17, f"{c_val:.2f}", va="center", fontsize=7.5, color=vs.INK)
    floor_label = "pending" if not np.isfinite(f_val) else ("0 by construction" if f_val == 0 else f"{f_val:.2f}")
    ax_bar.text((f_val if np.isfinite(f_val) else 0) + 0.01, yy - 0.17, floor_label, va="center", fontsize=7.5,
                color=vs.INK_2)
ax_bar.set_yticks(ypos)
ax_bar.set_yticklabels([r[0] for r in rows], fontsize=8, color=vs.INK)
ax_bar.tick_params(axis="y", length=0)
ax_bar.set_xlabel("median |component| (log variance)")
ax_bar.grid(True, axis="x")
ax_bar.set_axisbelow(True)
ax_bar.set_xlim(0, max(r[1] for r in rows) * 1.35)
ax_bar.set_title("Size of each component", fontsize=8, color=vs.INK, loc="left", pad=28)
ax_bar.text(0.0, 1.02, "solid: European vs Yoruba\nlight: two European draws (noise floor)"
            + ("  DRAFT" if floor_is_draft else ""), transform=ax_bar.transAxes, fontsize=7, color=vs.INK_2,
            va="bottom")

for a_, letter in ((ax_a, "a"), (ax_bar, "b")):
    bb = a_.get_position()
    fig.text(0.005 if a_ is ax_a else 0.50, 0.985, letter, fontsize=10, fontweight="bold", color=vs.INK, va="top")
vs.save(fig, "fig10_decomposition_triangle")

summary = {"genes": len(cross), "median_share_w": med[0], "median_share_D": med[1], "median_share_R": med[2],
           "median_abs_phi_w": rows[0][1], "median_abs_phi_D": rows[1][1], "median_abs_phi_R": rows[2][1],
           "floor_genes": len(floor) if floor is not None else 0, "floor_median_abs_phi_w": rows[0][2]}
pd.DataFrame([summary]).to_csv(OUT / "fig10_summary.tsv", sep="\t", index=False)
print(pd.DataFrame([summary]).T.to_string(header=False))
