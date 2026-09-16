"""Figure 15: how common each model's SNPs are in the other population.

Usage: python3 src/fig15_model_snps.py [eur_set] [yri_set]

Above the axis: every SNP used by a usable European model, binned by its minor allele
frequency in the Yoruba sample. Below: every SNP used by a usable Yoruba model, binned by its
frequency in the European sample. The first bin (MAF below 1%, rare or absent there) is
called out, because those SNPs are the ones a GWAS in the other population barely sees.
Frequencies are streamed from the union .afreq files, keeping only model SNPs.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import viz_style as vs

ROOT = Path(__file__).resolve().parents[1]
G = ROOT / "data" / "processed" / "geno"
OUT = ROOT / "results" / "figdata"
eur_set = sys.argv[1] if len(sys.argv) > 1 else "EUR87_r1"
yri_set = sys.argv[2] if len(sys.argv) > 2 else "YRI87"


def model_snps(set_name):
    d = ROOT / "results" / "models" / set_name
    chroms = sorted(int(f.name.split(".")[0][3:]) for f in d.glob("chr*.weights.tsv"))
    s = pd.concat([pd.read_csv(f, sep="\t") for f in d.glob("chr*.summary.tsv")])
    usable = set(s.loc[s["cv_r2"].notna() & (s["cv_r2"] > 0.01) & (s["cv_pval"] < 0.05) & (s["n_model"] > 0), "gene"])
    w = pd.concat([pd.read_csv(f, sep="\t") for f in d.glob("chr*.weights.tsv")])
    w = w[w["gene"].isin(usable)]
    return set(w["varID"]), chroms, len(usable)


def freqs_for(ids, pop):
    parts = []
    for chunk in pd.read_csv(G / f"union_freq_{pop}.afreq", sep="\t", usecols=["ID", "ALT_FREQS"], chunksize=2_000_000):
        parts.append(chunk[chunk["ID"].isin(ids)])
    f = pd.concat(parts)
    return np.minimum(f["ALT_FREQS"].to_numpy(), 1 - f["ALT_FREQS"].to_numpy())


snpsE, chrE, nE = model_snps(eur_set)
snpsY, chrY, nY = model_snps(yri_set)
shared_chr = set(chrE) & set(chrY)
draft = len(shared_chr) < 22
if draft:
    keep = lambda ids: {i for i in ids if int(i.split(":")[0]) in shared_chr}
    snpsE, snpsY = keep(snpsE), keep(snpsY)

mafE_in_Y = freqs_for(snpsE, "YRI87")
mafY_in_E = freqs_for(snpsY, "EUR358")

edges = np.concatenate([[0, 0.01], np.arange(0.05, 0.5001, 0.05)])
hE, _ = np.histogram(mafE_in_Y, bins=edges)
hY, _ = np.histogram(mafY_in_E, bins=edges)
fracE, fracY = hE / hE.sum(), hY / hY.sum()
rareE, rareY = fracE[0], fracY[0]

vs.apply()
fig, ax = plt.subplots(figsize=(vs.SINGLE, 2.9))
fig.subplots_adjust(left=0.17, right=0.97, top=0.83, bottom=0.24)
lefts, widths = edges[:-1], np.diff(edges)
ax.bar(lefts, fracE, width=widths * 0.9, align="edge", color=vs.EUR, lw=0)
ax.bar(lefts, -fracY, width=widths * 0.9, align="edge", color=vs.YRI, lw=0)
ax.axhline(0, color=vs.INK_2, lw=0.6)
lim = max(fracE.max(), fracY.max()) * 1.3
ax.set_ylim(-lim, lim)
ax.set_xlim(0, 0.5)
ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{abs(v) * 100:.0f}%"))
ax.set_ylabel("share of model SNPs")
ax.set_xlabel("minor allele frequency in the other population")
ax.text(0.5, -0.27, "bins: below 1% (rare or absent), 1 to 5%, then 5% steps\nbar height: share of that model set's SNPs in the bin",
        transform=ax.transAxes, ha="center", va="top", fontsize=7, color=vs.INK_2)
ax.grid(True, axis="y")
ax.set_axisbelow(True)
ax.annotate(f"{rareE:.0%} rare or absent\nin the Yoruba sample", (0.005, fracE[0]), xytext=(0.07, lim * 0.78),
            fontsize=7.5, color=vs.INK, arrowprops=dict(arrowstyle="-", color=vs.MUTED, lw=0.6))
ax.annotate(f"{rareY:.0%} rare or absent\nin the European sample", (0.005, -fracY[0]), xytext=(0.07, -lim * 0.78),
            fontsize=7.5, color=vs.INK, va="top", arrowprops=dict(arrowstyle="-", color=vs.MUTED, lw=0.6))
for yy, color, text in ((1.17, vs.EUR, f"above: SNPs in European models ({len(snpsE):,})"),
                        (1.06, vs.YRI, f"below: SNPs in Yoruba models ({len(snpsY):,})")):
    ax.scatter([0.0], [yy], s=12, color=color, transform=ax.transAxes, clip_on=False)
    ax.text(0.03, yy, text, transform=ax.transAxes, fontsize=7.5, color=vs.INK, va="center")
if draft:
    ax.text(1.0, 1.17, f"DRAFT: {len(shared_chr)} chr", transform=ax.transAxes, ha="right", va="center", fontsize=7,
            color=vs.INK_2)
vs.save(fig, "fig15_model_snps")

pd.DataFrame({"bin_left": lefts, "bin_right": edges[1:], "frac_eur_model_snps_by_yri_maf": fracE,
              "frac_yri_model_snps_by_eur_maf": fracY}).to_csv(OUT / "fig15_model_snps.tsv", sep="\t", index=False)
print(f"EUR model SNPs {len(snpsE):,} (usable genes {nE}); YRI model SNPs {len(snpsY):,} (usable genes {nY}); "
      f"chromosomes compared {len(shared_chr)}; rare in other population: EUR models {rareE:.3f}, YRI models {rareY:.3f}")
