"""Figure 5: the genome read twice, once by each model.

Usage: python3 src/fig05_accuracy_mirror.py [eur_set] [yri_set]

Every usable model is one thin stem at its gene's genome position. Stems above the axis are
the European model's cross-validated R^2, stems below are the Yoruba model's. Chromosomes
alternate a faint band so position stays readable without tick clutter. The best-predicted
genes in each set are named.
Numbers quoted in the text go to results/figdata/fig05_summary.tsv.
"""
import json
import sys
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import viz_style as vs

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "figdata"
eur_set = sys.argv[1] if len(sys.argv) > 1 else "EUR87_r1"
yri_set = sys.argv[2] if len(sys.argv) > 2 else "YRI87"
GAP = 12_000_000


def usable(set_name):
    files = sorted((ROOT / "results" / "models" / set_name).glob("chr*.summary.tsv"))
    s = pd.concat([pd.read_csv(f, sep="\t") for f in files])
    chroms = sorted(int(f.name.split(".")[0][3:]) for f in files)
    u = s[s["cv_r2"].notna() & (s["cv_r2"] > 0.01) & (s["cv_pval"] < 0.05) & (s["n_model"] > 0)]
    return u, chroms


def symbols(ids):
    cache_path = OUT / "gene_symbols.json"
    cache = json.loads(cache_path.read_text()) if cache_path.exists() else {}
    for g in ids:
        base = g.split(".")[0]
        if base in cache:
            continue
        try:
            with urllib.request.urlopen(f"https://grch37.rest.ensembl.org/lookup/id/{base}?content-type=application/json", timeout=15) as r:
                cache[base] = json.load(r).get("display_name") or base
        except Exception:
            cache[base] = base
    cache_path.write_text(json.dumps(cache, indent=1))
    return {g: cache[g.split(".")[0]] for g in ids}


ann = pd.read_csv(ROOT / "data" / "processed" / "expr" / "gene_annotation.tsv", sep="\t").set_index("gene")
chr_len = ann.groupby("chr")["coord"].max()
offset = (chr_len + GAP).cumsum().shift(fill_value=0)

uE, chrE = usable(eur_set)
uY, chrY = usable(yri_set)
shared_chroms = sorted(set(chrE) & set(chrY))
partial = len(shared_chroms) < 22
for u in (uE, uY):
    u["x"] = [offset[ann.loc[g, "chr"]] + ann.loc[g, "coord"] for g in u["gene"]]

def spaced_top(u, k=3, min_sep=200_000_000):
    """The best-predicted genes, keeping only ones at least min_sep apart so labels cannot pile up."""
    picked = []
    for _, r in u.sort_values("cv_r2", ascending=False).iterrows():
        if all(abs(r["x"] - q["x"]) >= min_sep for q in picked):
            picked.append(r)
        if len(picked) == k:
            break
    return pd.DataFrame(picked)


top = pd.concat([spaced_top(uE).assign(side="EUR"), spaced_top(uY).assign(side="YRI")])
names = symbols(top["gene"].tolist())

vs.apply()
fig, ax = plt.subplots(figsize=(vs.DOUBLE, 2.9))
fig.subplots_adjust(left=0.07, right=0.99, top=0.86, bottom=0.18)
lim = max(uE["cv_r2"].max(), uY["cv_r2"].max()) * 1.3
for c in range(1, 23):
    if c % 2 == 0:
        ax.axvspan(offset[c], offset[c] + chr_len[c], color="#f4f3ef", lw=0, zorder=0)
ax.vlines(uE["x"], 0, uE["cv_r2"], color=vs.EUR, lw=0.45, zorder=2)
ax.vlines(uY["x"], 0, -uY["cv_r2"], color=vs.YRI, lw=0.45, zorder=2)
ax.axhline(0, color=vs.INK_2, lw=0.6, zorder=3)
ax.set_xlim(-GAP / 2, offset[22] + chr_len[22] + GAP / 2)
ax.set_ylim(-lim, lim)
ax.set_xticks([])
ax.spines["bottom"].set_visible(False)
# Ticks every 0.2 so the one-decimal |value| labels are always exact.
ax.yaxis.set_major_locator(plt.MultipleLocator(0.2))
ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{abs(v):.1f}"))
ax.set_ylabel("cross-validated R²")

# Chromosome numbers sit just below the plotting area (data x, axes y), so no stem or gene label can reach them.
below = ax.get_xaxis_transform()
for c in range(1, 23):
    ax.text(offset[c] + chr_len[c] / 2, -0.015, str(c) if c < 19 or c % 2 == 1 else "", transform=below,
            ha="center", va="top", fontsize=7, color=vs.MUTED)
ax.text(-GAP / 2, -0.015, "chr", transform=below, ha="right", va="top", fontsize=7, color=vs.MUTED)

x_right = offset[22] + chr_len[22]
for _, r in top.iterrows():
    y = r["cv_r2"] if r["side"] == "EUR" else -r["cv_r2"]
    near_right = r["x"] > 0.85 * x_right
    # A ringed dot on the stem tip ties each name to exactly one gene.
    ax.scatter([r["x"]], [y], s=16, facecolors=vs.SURFACE, edgecolors=vs.INK, linewidths=0.8, zorder=5)
    ax.annotate(names[r["gene"]], (r["x"], y), xytext=(-5 if near_right else 5, 4 if y > 0 else -4),
                textcoords="offset points", fontsize=7, color=vs.INK, ha="right" if near_right else "left",
                va="bottom" if y > 0 else "top", fontstyle="italic")

def label(y, color, text, va):
    ax.scatter([0.012], [y], s=14, color=color, transform=ax.transAxes, clip_on=False, zorder=4)
    ax.text(0.022, y, text, transform=ax.transAxes, fontsize=8, color=vs.INK, va="center")

label(1.07, vs.EUR, f"above: European model ({eur_set}), {len(uE):,} usable genes, median R² {uE['cv_r2'].median():.2f}", "top")
label(-0.12, vs.YRI, f"below: Yoruba model ({yri_set}), {len(uY):,} usable genes, median R² {uY['cv_r2'].median():.2f}", "bottom")
if partial:
    ax.text(1.0, 1.07, f"DRAFT: training incomplete (EUR chromosomes {len(chrE)}, YRI chromosomes {len(chrY)})",
            transform=ax.transAxes, ha="right", va="center", fontsize=7.5, color=vs.INK_2)

vs.save(fig, "fig05_accuracy_mirror")
pd.DataFrame([{"eur_set": eur_set, "yri_set": yri_set, "eur_chromosomes": len(chrE), "yri_chromosomes": len(chrY),
               "eur_usable": len(uE), "yri_usable": len(uY), "eur_median_r2": uE["cv_r2"].median(),
               "yri_median_r2": uY["cv_r2"].median()}]).to_csv(OUT / "fig05_summary.tsv", sep="\t", index=False)
print(top[["side", "gene", "cv_r2"]].assign(name=top["gene"].map(names)).to_string(index=False))
