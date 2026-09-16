"""Figure 6: do the two models agree on which genes are predictable?

Usage: python3 src/fig06_r2_concordance.py [eur_set] [yri_set]

Each hexagon counts genes by their cross-validated R^2 in the European set (x) and the
Yoruba set (y), on square-root axes so the many small values are not crushed into a corner.
Dashed lines mark the usability threshold (R^2 = 0.01); the diagonal is equal accuracy.
Known strong LCL eQTL genes are ringed and named. The counts of genes usable in only one set
or in both are printed on the figure. Numbers quoted in the text go to results/figdata/.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, LogNorm
from scipy.stats import spearmanr

import viz_style as vs

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "figdata"
eur_set = sys.argv[1] if len(sys.argv) > 1 else "EUR87_r1"
yri_set = sys.argv[2] if len(sys.argv) > 2 else "YRI87"

CONTROLS = {"ENSG00000164308": "ERAP2", "ENSG00000196735": "HLA-DQA1", "ENSG00000124587": "PEX6",
            "ENSG00000134202": "GSTM3", "ENSG00000064886": "CHI3L2"}


def load(set_name):
    d = ROOT / "results" / "models" / set_name
    files = sorted(d.glob("chr*.summary.tsv"))
    s = pd.concat([pd.read_csv(f, sep="\t") for f in files])
    s["usable"] = s["cv_r2"].notna() & (s["cv_r2"] > 0.01) & (s["cv_pval"] < 0.05) & (s["n_model"] > 0)
    return s.set_index("gene"), len(files)


E, nE = load(eur_set)
Y, nY = load(yri_set)
draft = min(nE, nY) < 22
m = E[["cv_r2", "usable"]].join(Y[["cv_r2", "usable"]], lsuffix="_E", rsuffix="_Y", how="inner").dropna()
m["cv_r2_E"] = m["cv_r2_E"].clip(lower=0)
m["cv_r2_Y"] = m["cv_r2_Y"].clip(lower=0)
rho, rho_p = spearmanr(m["cv_r2_E"], m["cv_r2_Y"])
both = int((m["usable_E"] & m["usable_Y"]).sum())
only_E = int((m["usable_E"] & ~m["usable_Y"]).sum())
only_Y = int((~m["usable_E"] & m["usable_Y"]).sum())
both_u = m[m["usable_E"] & m["usable_Y"]]
rho_both, _ = spearmanr(both_u["cv_r2_E"], both_u["cv_r2_Y"]) if len(both_u) > 2 else (np.nan, np.nan)

vs.apply()
fig, ax = plt.subplots(figsize=(vs.SINGLE, 3.55))
fig.subplots_adjust(left=0.16, right=0.8, top=0.78, bottom=0.14)
cmap = LinearSegmentedColormap.from_list("dens", ["#e8f0fb"] + vs.SEQ[3:])
hb = ax.hexbin(np.sqrt(m["cv_r2_E"]), np.sqrt(m["cv_r2_Y"]), gridsize=38, extent=(0, 0.92, 0, 0.92),
               cmap=cmap, norm=LogNorm(), mincnt=1, linewidths=0, rasterized=True)
# Axis values are drawn on a square-root scale: positions are sqrt(R^2), ticks show R^2.
ticks = [0, 0.01, 0.05, 0.1, 0.2, 0.4, 0.6]
ax.set_xticks(np.sqrt(ticks))
ax.set_yticks(np.sqrt(ticks))
ax.set_xticklabels([f"{t:g}" for t in ticks])
ax.set_yticklabels([f"{t:g}" for t in ticks])
# Declare the transform so the tick check verifies labels against R^2 = position^2.
ax.xaxis.ibeji_value_of = np.square
ax.yaxis.ibeji_value_of = np.square
ax.set_xlim(0, 0.92)
ax.set_ylim(0, 0.92)
ax.set_aspect("equal")
ax.plot([0, 0.92], [0, 0.92], color=vs.MUTED, lw=0.6)
thr = np.sqrt(0.01)
ax.axvline(thr, color=vs.MUTED, lw=0.6, ls=(0, (3, 2)))
ax.axhline(thr, color=vs.MUTED, lw=0.6, ls=(0, (3, 2)))
ax.set_xlabel(f"cross-validated R², European model ({eur_set})")
ax.set_ylabel(f"cross-validated R², Yoruba model ({yri_set})")

m_base = m.copy()
m_base.index = [g.split(".")[0] for g in m_base.index]
for ens, name in CONTROLS.items():
    if ens not in m_base.index:
        continue
    r = m_base.loc[ens]
    x, y = np.sqrt(r["cv_r2_E"]), np.sqrt(r["cv_r2_Y"])
    ax.scatter([x], [y], s=18, facecolors="none", edgecolors=vs.INK, linewidths=0.8, zorder=4)
    right = x > 0.65
    ax.annotate(name, (x, y), xytext=(-5 if right else 5, 3), textcoords="offset points", fontsize=7,
                fontstyle="italic", color=vs.INK, ha="right" if right else "left")

cax = fig.add_axes([0.83, 0.3, 0.022, 0.4])
cb = fig.colorbar(hb, cax=cax)
cb.set_label("genes per hexagon", fontsize=7.5, color=vs.INK_2)
cb.ax.tick_params(labelsize=7)
cb.ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v:g}"))
cb.outline.set_visible(False)

fig.text(0.16, 0.965, f"{len(m):,} genes fit in both sets; Spearman ρ = {rho:.2f}", fontsize=8, color=vs.INK, va="top")
fig.text(0.16, 0.925, f"usable in both {both:,}   only European {only_E:,}   only Yoruba {only_Y:,}",
         fontsize=7.5, color=vs.INK_2, va="top")
fig.text(0.16, 0.885, "dashed: usability threshold R² = 0.01;  line: equal accuracy",
         fontsize=7, color=vs.INK_2, va="top")
fig.text(0.16, 0.855, "ringed: known strong eQTL genes;  axes on a square-root scale",
         fontsize=7, color=vs.INK_2, va="top")
if draft:
    fig.text(0.97, 0.965, "DRAFT", fontsize=7.5, color=vs.INK_2, va="top", ha="right")
vs.save(fig, "fig06_r2_concordance")

pd.DataFrame([{"eur_set": eur_set, "yri_set": yri_set, "genes_both_fit": len(m), "spearman_all": rho,
               "spearman_all_p": rho_p, "usable_both": both, "usable_only_eur": only_E, "usable_only_yri": only_Y,
               "spearman_usable_both": rho_both}]).to_csv(OUT / "fig06_summary.tsv", sep="\t", index=False)
print(f"genes {len(m):,}; rho {rho:.3f}; usable both {both}, only EUR {only_E}, only YRI {only_Y}; rho among usable-both {rho_both:.3f}")
