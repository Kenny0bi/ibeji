"""Figure 19: do the two model sets find the same cross-disorder genes?

Usage: python3 src/fig19_cross_disorder.py

(a, b) Bonferroni-significant genes shared between each pair of disorders, for the European
       and the Yoruba model set. Counts are small, so they are printed rather than left to a
       colour scale alone.
(c) Every gene significant in three or more disorders, and which model set found it.
Reads results/twas/<set>/twas_all_traits.tsv.
"""
import json
import urllib.request
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

import viz_style as vs

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "figdata"
TRAITS = [("SCZ", "SCZ"), ("BIP", "BIP"), ("MDD", "MDD"), ("PTSD", "PTSD")]
SETS = [("EUR87_r1", vs.EUR, "European models"), ("YRI87", vs.YRI, "Yoruba models")]


def symbols(ids):
    """Same resolver the other figures use: Ensembl GRCh37 REST, cached, ID as fallback."""
    cache_path = OUT / "gene_symbols.json"
    cache = json.loads(cache_path.read_text()) if cache_path.exists() else {}
    for g in ids:
        base = g.split(".")[0]
        if base in cache:
            continue
        try:
            with urllib.request.urlopen(
                    f"https://grch37.rest.ensembl.org/lookup/id/{base}?content-type=application/json",
                    timeout=15) as r:
                cache[base] = json.load(r).get("display_name") or base
        except Exception:
            cache[base] = base
    cache_path.write_text(json.dumps(cache, indent=1))
    return {g: cache[g.split(".")[0]] for g in ids}


sig = {}
for s, _, _ in SETS:
    d = pd.read_csv(ROOT / f"results/twas/{s}/twas_all_traits.tsv", sep="\t")
    for code, _ in TRAITS:
        sub = d[d.trait == code]
        sig[(s, code)] = set(sub.loc[sub.p < 0.05 / len(sub), "gene"])

vs.apply()
fig = plt.figure(figsize=(vs.DOUBLE, 3.1))
ax_a = fig.add_axes([0.075, 0.30, 0.20, 0.45])
ax_b = fig.add_axes([0.345, 0.30, 0.20, 0.45])
ax_c = fig.add_axes([0.765, 0.30, 0.175, 0.45])

ramp = LinearSegmentedColormap.from_list("cnt", ["#eef2fb", "#9ec5f4", "#3987e5", "#184f95"])
n = len(TRAITS)
vmax = 8
for ax, (s, color, label) in zip((ax_a, ax_b), SETS):
    M = np.full((n, n), np.nan)
    for i, j in combinations(range(n), 2):
        M[j, i] = len(sig[(s, TRAITS[i][0])] & sig[(s, TRAITS[j][0])])
    ax.imshow(np.ma.masked_invalid(M), cmap=ramp, vmin=0, vmax=vmax, origin="upper")
    for i, j in combinations(range(n), 2):
        v = int(M[j, i])
        ax.text(i, j, str(v), ha="center", va="center", fontsize=7.5,
                color=vs.SURFACE if v >= 5 else vs.INK)
    # Only the triangle's visible rows and columns get ticks. Defining all four left stray
    # label text that collided even though the limits hid the empty row and column.
    # Short codes, horizontal. Rotated disorder names cannot fit three across a 1.4 in axis:
    # they collided at every rotation and panel size tried.
    ax.set_xticks(range(n - 1)); ax.set_yticks(range(1, n))
    ax.set_xticklabels([t[1] for t in TRAITS[:-1]], fontsize=7.5)
    ax.set_yticklabels([t[1] for t in TRAITS[1:]], fontsize=7.5)
    ax.tick_params(length=0)
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.set_xlim(-0.5, n - 1.5); ax.set_ylim(n - 0.5, 0.5)
    ax.add_patch(plt.Rectangle((0.0, 1.10), 0.10, 0.075, transform=ax.transAxes, color=color, lw=0, clip_on=False))
    ax.text(0.14, 1.138, label, transform=ax.transAxes, fontsize=8, color=vs.INK, va="center")

ax_a.text(0.0, 1.30, "Genes significant in both disorders of a pair", transform=ax_a.transAxes,
          fontsize=8.5, color=vs.INK, va="bottom")
fig.text(0.075, 0.055, "SCZ schizophrenia,  BIP bipolar disorder,  MDD major depression,  "
         "PTSD post-traumatic stress disorder", fontsize=7.5, color=vs.INK_2)

multi = {}
for s, _, _ in SETS:
    for g in set().union(*[sig[(s, c)] for c, _ in TRAITS]):
        k = [c for c, _ in TRAITS if g in sig[(s, c)]]
        if len(k) >= 3:
            multi.setdefault(g, {})[s] = k
names = symbols(sorted(multi))
genes = sorted(multi, key=lambda g: (-max(len(v) for v in multi[g].values()), names[g]))
for gi, g in enumerate(genes):
    y = len(genes) - 1 - gi
    for si, (s, color, _) in enumerate(SETS):
        ks = multi[g].get(s, [])
        for ci, (code, _) in enumerate(TRAITS):
            if code in ks:
                ax_c.scatter([ci + (si - 0.5) * 0.26], [y], s=26, color=color, linewidths=0, zorder=3)
    ax_c.text(-0.85, y, names[g], ha="right", va="center", fontsize=7.5, fontstyle="italic", color=vs.INK)
ax_c.set_xticks(range(n))
ax_c.set_xticklabels([t[1] for t in TRAITS], fontsize=7.5)
ax_c.set_yticks([])
ax_c.set_xlim(-0.6, n - 0.4); ax_c.set_ylim(-0.7, len(genes) - 0.3)
ax_c.tick_params(length=0)
for sp in ax_c.spines.values():
    sp.set_visible(False)
ax_c.grid(True, axis="x"); ax_c.set_axisbelow(True)
# Titles start left of the narrow panel so they do not run off the figure edge.
ax_c.text(-0.62, 1.30, "Genes in three or more disorders", transform=ax_c.transAxes,
          fontsize=8.5, color=vs.INK, va="bottom")
ax_c.text(-0.62, 1.10, "left dot European, right dot Yoruba", transform=ax_c.transAxes,
          fontsize=7.5, color=vs.INK_2, va="bottom")

for a_, letter in ((ax_a, "a"), (ax_b, "b"), (ax_c, "c")):
    fig.text(max(a_.get_position().x0 - 0.062, 0.008), 0.975, letter, fontsize=9, fontweight="bold",
             color=vs.INK, va="top")
vs.save(fig, "fig19_cross_disorder")

rows = []
for s, _, label in SETS:
    for i, j in combinations(range(n), 2):
        rows.append(dict(model_set=label, pair=f"{TRAITS[i][1]} x {TRAITS[j][1]}",
                         n=len(sig[(s, TRAITS[i][0])] & sig[(s, TRAITS[j][0])])))
pd.DataFrame(rows).to_csv(OUT / "fig19_summary.tsv", sep="\t", index=False)
print(pd.DataFrame(rows).to_string(index=False))
print("\ngenes in 3+ disorders:", {names[g]: {k.split('87')[0]: v for k, v in multi[g].items()} for g in genes})
