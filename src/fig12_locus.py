"""Figure 12: anatomy of one disagreement.

Usage: python3 src/fig12_locus.py <ensembl_gene_id> <display_name> [eur_set] [yri_set]

One gene, read top to bottom, following V = w' D^{1/2} R D^{1/2} w:
  (a) weights: European model as stems above the axis, Yoruba model below, on genome position
  (b) allele frequency: for every SNP either model uses, a vertical link from its European
      to its Yoruba minor allele frequency
  (c) LD: r^2 among those SNPs, European triangle above, Yoruba below
  (d) the decomposition: delta split into phi_w, phi_D, phi_R as a signed waterfall
All four panels share one ordering of SNPs, so a SNP can be followed straight down.
"""
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import Rectangle

import viz_style as vs

ROOT = Path(__file__).resolve().parents[1]
P = ROOT / "tools" / "plink2"
G = ROOT / "data" / "processed" / "geno"
S = ROOT / "data" / "processed" / "samples"
TMP = ROOT / "results" / "figdata" / "locus"
TMP.mkdir(parents=True, exist_ok=True)

gene, name = sys.argv[1], sys.argv[2]
eur_set = sys.argv[3] if len(sys.argv) > 3 else "EUR87_r1"
yri_set = sys.argv[4] if len(sys.argv) > 4 else "YRI87"


def weights(set_name):
    chrom = int(ann.loc[gene, "chr"])
    w = pd.read_csv(ROOT / "results" / "models" / set_name / f"chr{chrom}.weights.tsv", sep="\t")
    return w[w["gene"] == gene].set_index("varID")["weight"]


def dosages(pop, snps):
    snp_file = TMP / f"{name}_snps.txt"
    pd.Series(snps).to_csv(snp_file, index=False, header=False)
    out = TMP / f"{name}_{pop}"
    subprocess.run([str(P), "--pfile", str(G / "geuv445"), "--keep", str(S / f"{pop}.txt"),
                    "--extract", str(snp_file), "--export", "A", "--export-allele", str(alt_file),
                    "--threads", "1", "--out", str(out)], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    raw = pd.read_csv(f"{out}.raw", sep="\t")
    x = raw.iloc[:, list(raw.columns).index("PHENOTYPE") + 1:]
    x.columns = [c.rsplit("_", 1)[0] for c in x.columns]
    return x[snps].to_numpy(dtype=float)


ann = pd.read_csv(ROOT / "data" / "processed" / "expr" / "gene_annotation.tsv", sep="\t").set_index("gene")
wE, wY = weights(eur_set), weights(yri_set)
snps = sorted(set(wE.index) | set(wY.index), key=lambda s: int(s.split(":")[1]))
alt_file = TMP / f"{name}_alt.txt"
pd.DataFrame({"id": snps, "alt": [s.split(":")[3] for s in snps]}).to_csv(alt_file, sep="\t", index=False, header=False)

xE, xY = dosages("EUR358", snps), dosages("YRI87", snps)
pE, pY = xE.mean(0) / 2, xY.mean(0) / 2
mafE, mafY = np.minimum(pE, 1 - pE), np.minimum(pY, 1 - pY)
rE = np.nan_to_num(np.corrcoef(xE, rowvar=False) ** 2)
rY = np.nan_to_num(np.corrcoef(xY, rowvar=False) ** 2)
we = wE.reindex(snps).fillna(0).to_numpy()
wy = wY.reindex(snps).fillna(0).to_numpy()
pos = np.array([int(s.split(":")[1]) for s in snps])
n = len(snps)

dec = pd.read_csv(ROOT / "results" / "decomposition" / f"{eur_set}_vs_{yri_set}.tsv", sep="\t").set_index("gene").loc[gene]

vs.apply()
MINUS = "−"


def signed(v):
    return f"{v:+.2f}".replace("-", MINUS)


fig = plt.figure(figsize=(vs.DOUBLE, 5.7))
ax = fig.add_axes([0.08, 0.60, 0.60, 0.25])
axf = fig.add_axes([0.08, 0.40, 0.60, 0.14], sharex=ax)
axl = fig.add_axes([0.08, 0.10, 0.60, 0.24], sharex=ax)
cax = fig.add_axes([0.08, 0.055, 0.18, 0.012])
axd = fig.add_axes([0.80, 0.40, 0.185, 0.45])
xs = np.arange(n)

# (a) weights on an evenly spaced SNP index
lim = max(np.abs(we).max(), np.abs(wy).max()) * 1.25
shared = (we != 0) & (wy != 0)
for k in xs[shared]:
    ax.axvspan(k - 0.45, k + 0.45, color=vs.GRID, lw=0, zorder=0)
for w, color, sign in ((we, vs.EUR, 1), (wy, vs.YRI, -1)):
    nz = w != 0
    mag = np.abs(w[nz]) * sign
    ax.vlines(xs[nz], 0, mag, color=color, lw=1.1, zorder=2)
    pos_w = w[nz] > 0
    ax.scatter(xs[nz][pos_w], mag[pos_w], s=14, color=color, edgecolors=vs.SURFACE, linewidths=0.6, zorder=3)
    ax.scatter(xs[nz][~pos_w], mag[~pos_w], s=14, color=vs.SURFACE, edgecolors=color, linewidths=0.9, zorder=3)
ax.axhline(0, color=vs.INK_2, lw=0.6, zorder=1)
ax.set_xlim(-1, n)
ax.set_ylim(-lim, lim)
ax.set_xticks([])
ax.spines["bottom"].set_visible(False)
ax.set_ylabel("|weight|")
ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{abs(v):.1f}"))
ax.scatter([0.018], [0.925], s=16, color=vs.EUR, transform=ax.transAxes, clip_on=False, zorder=4)
ax.text(0.032, 0.925, f"European model ({eur_set}), {int((we != 0).sum())} SNPs", transform=ax.transAxes,
        fontsize=8, color=vs.INK, va="center")
ax.scatter([0.018], [0.075], s=16, color=vs.YRI, transform=ax.transAxes, clip_on=False, zorder=4)
ax.text(0.032, 0.075, f"Yoruba model ({yri_set}), {int((wy != 0).sum())} SNPs", transform=ax.transAxes,
        fontsize=8, color=vs.INK, va="center")
ax.text(1.0, 1.03, f"filled = positive weight, hollow = negative,  gray band = SNP in both models ({int(shared.sum())})",
        transform=ax.transAxes, fontsize=7.5, color=vs.INK_2, va="bottom", ha="right")

# (b) minor allele frequency in each sample, joined per SNP
axf.vlines(xs, mafE, mafY, color=vs.PHI_D, lw=1.0, alpha=0.9, zorder=1)
axf.scatter(xs, mafE, s=11, color=vs.EUR, zorder=3, linewidths=0)
axf.scatter(xs, mafY, s=11, color=vs.YRI, zorder=3, linewidths=0)
axf.set_ylim(-0.03, 0.53)
axf.set_yticks([0, 0.25, 0.5])
axf.set_ylabel("MAF")
axf.set_xticks([])
axf.spines["bottom"].set_visible(False)
axf.grid(True, axis="y")
axf.set_axisbelow(True)
from matplotlib.lines import Line2D
axf.legend([Line2D([], [], marker="o", ls="", color=vs.EUR, markersize=3.5),
            Line2D([], [], marker="o", ls="", color=vs.YRI, markersize=3.5),
            Line2D([], [], color=vs.PHI_D, lw=1.0, alpha=0.9)],
           ["European sample", "Yoruba sample", "frequency difference"],
           loc="lower right", bbox_to_anchor=(1.0, 1.0), ncol=3, fontsize=7.5, handlelength=1.4,
           columnspacing=1.2, borderaxespad=0.1)

# (c) mirrored LD triangles on the same SNP index
i, j = np.triu_indices(n, k=1)
cx, cy = (i + j) / 2.0, (j - i) / 2.0
diamond = np.stack([np.stack([cx, cy - 0.5], 1), np.stack([cx + 0.5, cy], 1),
                    np.stack([cx, cy + 0.5], 1), np.stack([cx - 0.5, cy], 1)], axis=1)
cmap = LinearSegmentedColormap.from_list("ld", ["#fbfcff"] + vs.SEQ[2:])
for r2, sign in ((rE, 1), (rY, -1)):
    v = diamond.copy(); v[:, :, 1] *= sign
    pc = PolyCollection(v, array=r2[i, j], cmap=cmap, norm=Normalize(0, 1), edgecolors="none", rasterized=True)
    axl.add_collection(pc)
axl.axhline(0, color=vs.SURFACE, lw=1.0)
half = (n - 1) / 2
axl.set_ylim(-half - 0.5, half + 0.5)
axl.axis("off")
axl.text(0.0, 0.97, "LD in the European sample", transform=axl.transAxes, fontsize=7.5, color=vs.INK, va="top")
axl.text(0.0, 0.03, "LD in the Yoruba sample", transform=axl.transAxes, fontsize=7.5, color=vs.INK, va="bottom")
cb = fig.colorbar(pc, cax=cax, orientation="horizontal")
cb.set_ticks([0, 0.5, 1])
cb.ax.tick_params(labelsize=7.5, length=2)
cb.outline.set_visible(False)
cax.text(1.06, 0.5, "r²", transform=cax.transAxes, fontsize=7.5, color=vs.INK_2, va="center")

# (d) decomposition waterfall
steps = [("w", dec["phi_w"], vs.PHI_W), ("D", dec["phi_D"], vs.PHI_D), ("R", dec["phi_R"], vs.PHI_R)]
level = 0.0
tops = [0.0]
for k, (label, val, color) in enumerate(steps):
    axd.bar(k, val, bottom=level, width=0.66, color=color, lw=0, zorder=2)
    level += val
    tops.append(level)
axd.bar(3, dec["delta"], width=0.66, color=vs.INK_2, lw=0, zorder=2)
axd.axhline(0, color=vs.INK_2, lw=0.6, zorder=1)
axd.set_xticks(range(4))
axd.set_xticklabels([f"$\\phi_w$\n{signed(dec['phi_w'])}", f"$\\phi_D$\n{signed(dec['phi_D'])}",
                     f"$\\phi_R$\n{signed(dec['phi_R'])}", f"$\\Delta$\n{signed(dec['delta'])}"], fontsize=7.5)
axd.set_xlim(-0.6, 3.6)
lo, hi = min(tops + [dec["delta"]]), max(tops + [dec["delta"]])
span = hi - lo
axd.set_ylim(lo - 0.22 * span - 0.05, hi + 0.22 * span + 0.05)
axd.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v:.2f}".replace("-", MINUS)))
axd.set_ylabel(r"$\log V_{\mathrm{YRI}} - \log V_{\mathrm{EUR}}$", fontsize=8)
axd.grid(True, axis="y")
axd.set_axisbelow(True)
axd.text(0.5, 1.03, "Where the gap comes from", transform=axd.transAxes, ha="center", va="bottom",
         fontsize=8, color=vs.INK)
keyax = fig.add_axes([0.79, 0.08, 0.19, 0.22])
keyax.set_xlim(0, 10); keyax.set_ylim(0, 10); keyax.axis("off")
keyax.text(0, 9.6, "Key for panel d", fontsize=8, color=vs.INK, va="top")
entries = [(vs.PHI_W, r"$\phi_w$", "weights differ"),
           (vs.PHI_D, r"$\phi_D$", "allele frequencies differ"),
           (vs.PHI_R, r"$\phi_R$", "LD differs"),
           (vs.INK_2, r"$\Delta$", "total gap (their sum)")]
for r, (color, sym, text) in enumerate(entries):
    y = 7.2 - r * 2.0
    keyax.add_patch(Rectangle((0, y - 0.55), 1.1, 1.1, color=color, lw=0))
    keyax.text(1.6, y, sym, fontsize=8, color=vs.INK, va="center")
    keyax.text(3.0, y, text, fontsize=7.5, color=vs.INK_2, va="center")

chrom = int(ann.loc[gene, "chr"])
fig.text(0.08, 0.985, name, fontsize=10, fontweight="bold", color=vs.INK, va="top")
fig.text(0.08, 0.952,
         f"chr{chrom}:{pos.min() / 1e6:.2f} to {pos.max() / 1e6:.2f} Mb,  {n} SNPs used by either model",
         fontsize=7.5, color=vs.INK_2, va="top")
fig.text(0.08, 0.927,
         f"European model predicts Yoruba expression at r² {dec['r2_eur_model_in_yri']:.2f};  "
         f"Yoruba model predicts European expression at r² {dec['r2_yri_model_in_eur']:.2f}",
         fontsize=7.5, color=vs.INK_2, va="top")


def panel_label(a, letter, dx=-0.06):
    b = a.get_position()
    fig.text(b.x0 + dx, b.y1 + 0.012, letter, fontsize=10, fontweight="bold", color=vs.INK, va="bottom")


panel_label(ax, "a"); panel_label(axf, "b"); panel_label(axl, "c"); panel_label(axd, "d", dx=-0.09)

vs.save(fig, f"fig12_locus_{name}")
pd.DataFrame({"varID": snps, "pos": pos, "w_EUR": we, "w_YRI": wy, "maf_EUR": mafE, "maf_YRI": mafY}).to_csv(
    TMP / f"{name}_snps_detail.tsv", sep="\t", index=False)
print(dec[["delta", "phi_w", "phi_D", "phi_R"]].to_string())
