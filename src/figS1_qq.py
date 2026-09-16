"""Figure S1: QQ plots and genomic inflation for every trait and both model sets.

Usage: python3 src/figS1_qq.py

One panel per trait. Observed against expected -log10 p under the null, with the European
model set and the Yoruba model set drawn together so their calibration can be compared
directly. Genomic inflation is the median chi-square statistic divided by its null median.
Reads results/twas/<set>/twas_all_traits.tsv.
"""
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats

import viz_style as vs

ROOT = Path(__file__).resolve().parents[1]
TRAITS = [("ASD", "Autism"), ("SCZ", "Schizophrenia"), ("BIP", "Bipolar disorder"),
          ("MDD", "Major depression"), ("PTSD", "PTSD")]
SETS = [("EUR87_r1", vs.EUR, "European models"), ("YRI87", vs.YRI, "Yoruba models")]

data = {s: pd.read_csv(ROOT / f"results/twas/{s}/twas_all_traits.tsv", sep="\t") for s, _, _ in SETS}

vs.apply()
fig = plt.figure(figsize=(vs.DOUBLE, 2.15))
n = len(TRAITS)
left, right, gap = 0.055, 0.995, 0.019
w = (right - left - gap * (n - 1)) / n
axes, rows = [], []
for i, (code, name) in enumerate(TRAITS):
    ax = fig.add_axes([left + i * (w + gap), 0.215, w, 0.475])
    axes.append(ax)
    hi = 0.0
    for s, color, slabel in SETS:
        d = data[s]
        p = d.loc[d.trait == code, "p"].dropna().values
        p = np.clip(p, 1e-300, 1.0)
        obs = -np.log10(np.sort(p))
        exp = -np.log10((np.arange(1, len(p) + 1) - 0.5) / len(p))
        lam = float(np.median(stats.chi2.isf(p, df=1)) / stats.chi2.ppf(0.5, df=1))
        ax.plot(exp, obs, color=color, lw=1.4, solid_capstyle="round", zorder=3)
        hi = max(hi, obs.max(), exp.max())
        rows.append(dict(trait=name, model_set=slabel, n_genes=len(p), lambda_gc=lam, max_obs=float(obs.max())))
    lim = hi * 1.06
    ax.plot([0, lim], [0, lim], color=vs.MUTED, lw=0.6, zorder=2)
    ax.set_xlim(0, lim); ax.set_ylim(0, lim)
    ax.set_aspect("equal")
    ax.grid(True); ax.set_axisbelow(True)
    ax.set_title(name, fontsize=8, color=vs.INK, loc="left", pad=16)
    lam_txt = "  ".join(f"{lam_:.2f}" for lam_ in
                        [r["lambda_gc"] for r in rows if r["trait"] == name])
    ax.text(0.0, 1.02, f"$\\lambda$ {lam_txt}", transform=ax.transAxes, fontsize=7.5, color=vs.INK_2, va="bottom")
    if i == 0:
        ax.set_ylabel("observed $-\\log_{10} p$")
axes[2].set_xlabel("expected $-\\log_{10} p$ under the null")

# The key sits on its own line above every panel title; a centred note on the same line
# collided with the middle panel titles in the first version.
for i, (s, color, slabel) in enumerate(SETS):
    x = 0.055 + 0.175 * i
    fig.patches.append(plt.Rectangle((x, 0.935), 0.015, 0.040, transform=fig.transFigure, color=color, lw=0))
    fig.text(x + 0.021, 0.955, slabel, fontsize=7.5, color=vs.INK_2, va="center")
fig.text(0.44, 0.955, "$\\lambda$ is listed European then Yoruba", fontsize=7.5, color=vs.INK_2, va="center")
vs.save(fig, "figS1_qq")

out = pd.DataFrame(rows)
out.to_csv(ROOT / "results/figdata/figS1_summary.tsv", sep="\t", index=False)
print(out.round(3).to_string(index=False))
