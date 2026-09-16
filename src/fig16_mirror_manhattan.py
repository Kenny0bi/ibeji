"""Figures 16 and 17: TWAS in the mirror.

Usage: python3 src/fig16_mirror_manhattan.py [eur_set] [yri_set]

Figure 16: autism. Every gene tested with the European model is a dot above the axis at
           -log10 p; every gene tested with the Yoruba model is a dot below. A gene
           significant with both models is joined across the axis by a thin line.
Figure 17: the same mirror for schizophrenia, bipolar disorder, major depression and PTSD.

Significance is Bonferroni within each trait and model set, drawn as a threshold line on
each side. If a set's genome-wide TWAS file is missing, the script uses whatever
per-chromosome TWAS files exist and marks the figure as a draft.
"""
import json
import sys
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

import viz_style as vs

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "figdata"
eur_set = sys.argv[1] if len(sys.argv) > 1 else "EUR87_r1"
yri_set = sys.argv[2] if len(sys.argv) > 2 else "YRI87"
GAP = 12_000_000
TRAIT_NAMES = {"ASD": "Autism", "SCZ": "Schizophrenia", "BIP": "Bipolar disorder",
               "MDD": "Major depression", "PTSD": "PTSD"}

ann = pd.read_csv(ROOT / "data" / "processed" / "expr" / "gene_annotation.tsv", sep="\t").set_index("gene")
chr_len = ann.groupby("chr")["coord"].max()
offset = (chr_len + GAP).cumsum().shift(fill_value=0)
X_END = offset[22] + chr_len[22]


def load_twas(set_name):
    d = ROOT / "results" / "twas" / set_name
    full = d / "twas_all_traits.tsv"
    if full.exists():
        return pd.read_csv(full, sep="\t"), False
    parts = sorted(d.glob("twas_chr*.tsv"))
    if not parts:
        raise SystemExit(f"no TWAS results for {set_name}")
    t = pd.concat([pd.read_csv(p, sep="\t") for p in parts])
    return t, True


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


def prepare(t, trait):
    d = t[(t["trait"] == trait) & t["p"].notna()].copy()
    d["x"] = [offset[ann.loc[g, "chr"]] + ann.loc[g, "coord"] for g in d["gene"]]
    d["chr"] = [int(ann.loc[g, "chr"]) for g in d["gene"]]
    d["nlp"] = -np.log10(d["p"].clip(lower=1e-300))
    n = len(d)
    d["sig"] = d["p"] < 0.05 / n
    return d, -np.log10(0.05 / n)


def mirror(ax, E, Y, thrE, thrY, label_k=3, fontsize=7, show_chr=True):
    for c in range(1, 23):
        if c % 2 == 0:
            ax.axvspan(offset[c], offset[c] + chr_len[c], color="#f4f3ef", lw=0, zorder=0)
    for d, color, sign in ((E, vs.EUR, 1), (Y, vs.YRI, -1)):
        # One opacity for every dot: the gray chromosome bands mark boundaries, so no second, unexplained cue.
        ax.scatter(d["x"], sign * d["nlp"], s=3, color=color, linewidths=0, zorder=2, rasterized=True)
        sig = d[d["sig"]]
        ax.scatter(sig["x"], sign * sig["nlp"], s=9, color=color, edgecolors=vs.SURFACE, linewidths=0.4, zorder=3)
    both = E[E["sig"]].merge(Y[Y["sig"]], on="gene", suffixes=("_E", "_Y"))
    ax.vlines(both["x_E"], -both["nlp_Y"], both["nlp_E"], color=vs.INK_2, lw=0.6, zorder=1)
    ax.axhline(0, color=vs.INK_2, lw=0.6, zorder=4)
    ax.axhline(thrE, color=vs.MUTED, lw=0.6, ls=(0, (3, 2)), zorder=1)
    ax.axhline(-thrY, color=vs.MUTED, lw=0.6, ls=(0, (3, 2)), zorder=1)
    top = max(E["nlp"].max(), thrE) * 1.25
    bot = max(Y["nlp"].max(), thrY) * 1.25
    ax.set_ylim(-bot, top)
    ax.set_xlim(-GAP / 2, X_END + GAP / 2)
    ax.set_xticks([])
    ax.spines["bottom"].set_visible(False)
    # Whole-number tick positions, so the |value| formatter can never round a tick like 2.5 to "2".
    ax.yaxis.set_major_locator(MaxNLocator(integer=True, nbins=6))
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{abs(v):.0f}"))
    if show_chr:
        below = ax.get_xaxis_transform()
        for c in range(1, 23):
            ax.text(offset[c] + chr_len[c] / 2, -0.015, str(c) if c < 19 or c % 2 == 1 else "", transform=below,
                    ha="center", va="top", fontsize=7, color=vs.MUTED)
    lab = pd.concat([E[E["sig"]].nsmallest(label_k, "p").assign(side=1), Y[Y["sig"]].nsmallest(label_k, "p").assign(side=-1)])
    names = symbols(lab["gene"].tolist()) if len(lab) else {}
    placed = []
    for _, r in lab.iterrows():
        if any(abs(r["x"] - px) < 150_000_000 and ps == r["side"] for px, ps in placed):
            continue
        y = r["side"] * r["nlp"]
        near_right = r["x"] > 0.85 * X_END
        # There is room for one label per side at the right edge. A second one would sit on the first
        # label or on the data, so only the stronger hit is named (lab is ordered by p within each side).
        if near_right and any(px > 0.85 * X_END and ps == r["side"] for px, ps in placed):
            continue
        placed.append((r["x"], r["side"]))
        ax.scatter([r["x"]], [y], s=16, facecolors="none", edgecolors=vs.INK, linewidths=0.8, zorder=5)
        ax.annotate(names[r["gene"]], (r["x"], y), xytext=(-5 if near_right else 5, 3 * r["side"]),
                    textcoords="offset points", fontsize=fontsize, fontstyle="italic", color=vs.INK,
                    ha="right" if near_right else "left", va="bottom" if y > 0 else "top")
    return both


def summary_row(trait, E, Y, both):
    return {"trait": trait, "genes_eur_model": len(E), "genes_yri_model": len(Y),
            "sig_eur_model": int(E["sig"].sum()), "sig_yri_model": int(Y["sig"].sum()),
            "sig_both": len(both), "sig_yri_only": int(Y["sig"].sum()) - len(both),
            "sig_eur_only": int(E["sig"].sum()) - len(both)}


tE, draftE = load_twas(eur_set)
tY, draftY = load_twas(yri_set)
if draftE or draftY:
    # Draft mode: compare only chromosomes present for both sets.
    chroms = set(ann.loc[tE["gene"].unique(), "chr"]) & set(ann.loc[tY["gene"].unique(), "chr"])
    tE = tE[tE["gene"].map(ann["chr"]).isin(chroms)]
    tY = tY[tY["gene"].map(ann["chr"]).isin(chroms)]
draft = draftE or draftY
rows = []

# ---- Figure 16: autism
vs.apply()
fig, ax = plt.subplots(figsize=(vs.DOUBLE, 3.1))
fig.subplots_adjust(left=0.07, right=0.99, top=0.85, bottom=0.17)
E, thrE = prepare(tE, "ASD")
Y, thrY = prepare(tY, "ASD")
both = mirror(ax, E, Y, thrE, thrY)
rows.append(summary_row("ASD", E, Y, both))
ax.set_ylabel("−log10 p")
for yy, color, text in ((1.08, vs.EUR, f"above: European model ({eur_set}), {len(E):,} genes, {int(E['sig'].sum())} significant"),
                        (-0.13, vs.YRI, f"below: Yoruba model ({yri_set}), {len(Y):,} genes, {int(Y['sig'].sum())} significant")):
    ax.scatter([0.012], [yy], s=14, color=color, transform=ax.transAxes, clip_on=False)
    ax.text(0.022, yy, text, transform=ax.transAxes, fontsize=8, color=vs.INK, va="center")
ax.text(1.0, 1.08, "dashed: Bonferroni threshold; gray line: significant with both models"
        + ("   DRAFT: partial chromosomes" if draft else ""), transform=ax.transAxes, ha="right", va="center",
        fontsize=7.5, color=vs.INK_2)
vs.save(fig, "fig16_mirror_manhattan_ASD")

# ---- Figure 17: four disorders
fig, axes = plt.subplots(2, 2, figsize=(vs.DOUBLE, 4.4))
fig.subplots_adjust(left=0.07, right=0.99, top=0.855, bottom=0.05, hspace=0.3, wspace=0.14)
for a, trait in zip(axes.ravel(), ["SCZ", "BIP", "MDD", "PTSD"]):
    E, thrE = prepare(tE, trait)
    Y, thrY = prepare(tY, trait)
    both = mirror(a, E, Y, thrE, thrY, label_k=2, fontsize=7, show_chr=False)
    rows.append(summary_row(trait, E, Y, both))
    a.set_title(f"{TRAIT_NAMES[trait]}: {int(E['sig'].sum())} above, {int(Y['sig'].sum())} below, {len(both)} both",
                fontsize=8, color=vs.INK, loc="left", pad=3)
    a.tick_params(axis="y", labelsize=7)
# Legend on its own small axes so the markers stay round (figure coordinates are not square).
lax = fig.add_axes([0.07, 0.915, 0.42, 0.075])
lax.set_xlim(0, 1); lax.set_ylim(0, 1); lax.axis("off")
for yy, color, text in ((0.75, vs.EUR, f"above each axis: European model ({eur_set})"),
                        (0.25, vs.YRI, f"below each axis: Yoruba model ({yri_set})")):
    lax.scatter([0.01], [yy], s=14, color=color)
    lax.text(0.03, yy, text, fontsize=7.5, color=vs.INK, va="center")
fig.text(0.99, 0.943, "dashed: Bonferroni threshold; gray line: significant with both models",
         ha="right", va="center", fontsize=7.5, color=vs.INK_2)
if draft:
    fig.text(0.99, 0.975, "DRAFT: partial chromosomes", ha="right", va="center", fontsize=7.5, color=vs.INK_2)
vs.save(fig, "fig17_mirror_manhattan_four")

pd.DataFrame(rows).to_csv(OUT / "fig16_17_summary.tsv", sep="\t", index=False)
print(pd.DataFrame(rows).to_string(index=False))
