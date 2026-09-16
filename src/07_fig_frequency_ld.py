"""Figures 2 and 3: the genetic ground under each twin model.

Figure 2 (a) Genome map: every 1 Mb tile of chr1 to chr22, colored by how far its share of
             Yoruba-only common SNVs sits from the genome-wide share.
         (b) Mirrored allele-frequency spectrum: European above the axis, Yoruba
             below, binned on Yoruba allele counts (k/174) so the discrete
             frequencies of an 87-person sample cannot create stripes.
Figure 3 (a) Mirrored LD triangles around DNAJB7 (chr22), European above, Yoruba below.
         (b) Mean r^2 by distance on chr22.
         Both panels use EUR87_r1 against YRI87, so sample size cannot bias r^2.

Numbers quoted in the text are written to results/figdata/.
"""
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.colors import LinearSegmentedColormap, PowerNorm
from matplotlib.patches import Polygon, Rectangle
from matplotlib.ticker import FuncFormatter

import viz_style as vs

ROOT = Path(__file__).resolve().parents[1]
P = ROOT / "tools" / "plink2"
G = ROOT / "data" / "processed" / "geno"
S = ROOT / "data" / "processed" / "samples"
OUT = ROOT / "results" / "figdata"
OUT.mkdir(parents=True, exist_ok=True)
PFILE = str(G / "geuv445")
YRI_CHR = 174


def plink(*args):
    # 2 threads: this runs beside GReX training, which owns the 4 physical cores.
    subprocess.run([str(P), "--pfile", PFILE, "--threads", "2", "--memory", "3000", *args], check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def union_freqs():
    for pop in ("EUR358", "YRI87"):
        out = G / f"union_freq_{pop}"
        if not Path(f"{out}.afreq").exists():
            plink("--keep", str(S / f"{pop}.txt"), "--freq", "--out", str(out))


# ----------------------------------------------------------------------------- Figure 2

def spectrum_counts():
    cache = OUT / "fig02_spectrum.tsv"
    classes = OUT / "fig02_frequency_classes.tsv"
    if cache.exists() and classes.exists():
        return pd.read_csv(cache, sep="\t"), pd.read_csv(classes, sep="\t").iloc[0]
    edges = (np.arange(0, 88) + 0.5) / YRI_CHR  # bin k (1..87) holds MAF in [(k-0.5)/174, (k+0.5)/174)
    nb = len(edges) - 1
    h = {k: np.zeros(nb, dtype=np.int64) for k in ("E_shared", "E_only", "Y_shared", "Y_only")}
    counts = {"both": 0, "EUR_only": 0, "YRI_only": 0}
    kw = dict(sep="\t", usecols=["ID", "ALT_FREQS"], chunksize=2_000_000)
    for ce, cy in zip(pd.read_csv(G / "union_freq_EUR358.afreq", **kw), pd.read_csv(G / "union_freq_YRI87.afreq", **kw)):
        assert (ce["ID"].values == cy["ID"].values).all()
        fe, fy = ce["ALT_FREQS"].to_numpy(), cy["ALT_FREQS"].to_numpy()
        me, my = np.minimum(fe, 1 - fe), np.minimum(fy, 1 - fy)
        oe, oy = me >= 0.01, my >= 0.01
        ie = np.clip(np.searchsorted(edges, me, side="right") - 1, 0, nb - 1)
        iy = np.clip(np.searchsorted(edges, my, side="right") - 1, 0, nb - 1)
        h["E_shared"] += np.bincount(ie[oe & oy], minlength=nb)[:nb]
        h["E_only"] += np.bincount(ie[oe & ~oy], minlength=nb)[:nb]
        h["Y_shared"] += np.bincount(iy[oe & oy], minlength=nb)[:nb]
        h["Y_only"] += np.bincount(iy[~oe & oy], minlength=nb)[:nb]
        counts["both"] += int((oe & oy).sum())
        counts["EUR_only"] += int((oe & ~oy).sum())
        counts["YRI_only"] += int((~oe & oy).sum())
    spec = pd.DataFrame({"maf_center": np.arange(1, 88) / YRI_CHR, **h})
    spec.to_csv(cache, sep="\t", index=False)
    total = sum(counts.values())
    row = {**counts, "total": total,
           "frac_both": counts["both"] / total, "frac_EUR_only": counts["EUR_only"] / total,
           "frac_YRI_only": counts["YRI_only"] / total,
           "yri_common": counts["both"] + counts["YRI_only"], "eur_common": counts["both"] + counts["EUR_only"],
           "frac_yri_common_not_eur": counts["YRI_only"] / (counts["both"] + counts["YRI_only"]),
           "frac_eur_common_not_yri": counts["EUR_only"] / (counts["both"] + counts["EUR_only"])}
    pd.DataFrame([row]).to_csv(classes, sep="\t", index=False)
    return spec, pd.Series(row)


def mirrored_spectrum(ax, spec, cls):
    x = spec["maf_center"].to_numpy()
    w = 0.78 / YRI_CHR
    faint = 0.38
    ax.bar(x, spec["E_shared"], width=w, color=vs.EUR, lw=0)
    ax.bar(x, spec["E_only"], width=w, bottom=spec["E_shared"], color=vs.EUR, alpha=faint, lw=0)
    ax.bar(x, -spec["Y_shared"], width=w, color=vs.YRI, lw=0)
    ax.bar(x, -spec["Y_only"], width=w, bottom=-spec["Y_shared"], color=vs.YRI, alpha=faint, lw=0)
    ax.axhline(0, color=vs.INK_2, lw=0.7)
    ax.set_xlim(0, 0.505)
    ax.set_xlabel("Minor allele frequency in that sample")
    ax.set_ylabel("SNVs per frequency bin")
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{abs(v) / 1e6:.1f}M" if v else "0"))
    top = (spec["E_shared"] + spec["E_only"]).max()
    bot = (spec["Y_shared"] + spec["Y_only"]).max()
    ax.set_ylim(-bot * 1.12, top * 1.55)
    ax.grid(True, axis="y")
    ax.set_axisbelow(True)

    tr = ax.transAxes
    ax.text(1.0, 0.99, "European sample (EUR358)", transform=tr, ha="right", va="top", fontsize=8,
            fontweight="bold", color=vs.INK)
    ax.text(1.0, 0.905, f"{cls['frac_eur_common_not_yri']:.0%} of its common SNVs are not common in Ibadan",
            transform=tr, ha="right", va="top", fontsize=7.5, color=vs.INK_2)
    ax.text(1.0, 0.01, "Yoruba sample, Ibadan (YRI87)", transform=tr, ha="right", va="bottom", fontsize=8,
            fontweight="bold", color=vs.INK)
    ax.text(1.0, 0.095, f"{cls['frac_yri_common_not_eur']:.0%} of its common SNVs are not common in the European sample",
            transform=tr, ha="right", va="bottom", fontsize=7.5, color=vs.INK_2)

    # Key drawn with the exact bar colors: a solid indigo+camwood pair and a light pair,
    # so the swatches are literally the same fills as the bars they describe.
    key_ax = ax.inset_axes([0.70, 0.40, 0.30, 0.20])
    key_ax.set_xlim(0, 10); key_ax.set_ylim(0, 4); key_ax.axis("off")
    for row, (alpha, label) in enumerate(((1.0, "common in both samples"), (faint, "common in this sample only"))):
        y = 2.6 - row * 1.8
        key_ax.add_patch(Rectangle((0.0, y), 0.55, 0.9, color=vs.EUR, alpha=alpha, lw=0))
        key_ax.add_patch(Rectangle((0.65, y), 0.55, 0.9, color=vs.YRI, alpha=alpha, lw=0))
        key_ax.text(1.55, y + 0.45, label, fontsize=7.5, color=vs.INK_2, va="center")


def genome_map(fig, ax, cax, min_snvs=500, lim=0.15):
    """Every 1 Mb tile of the genome, colored by how far its Yoruba-only share sits from the genome-wide share.

    Rows are chromosomes 1 to 22; columns are megabases. Tiles with fewer than `min_snvs`
    common SNVs (centromeres, gaps, chromosome ends) are left blank rather than shown noisy.
    """
    t = pd.read_csv(OUT / "fig02_tiles.tsv", sep="\t")
    genome = t["yri_only"].sum() / t["total"].sum()
    keep = t[t["total"] >= min_snvs]
    M = np.full((22, int(t["tile_mb"].max()) + 1), np.nan)
    M[keep["chrom"].to_numpy() - 1, keep["tile_mb"].to_numpy()] = keep["frac_yri_only"].to_numpy() - genome
    cmap = LinearSegmentedColormap.from_list("twin_div", [vs.EUR, "#f0efec", vs.YRI])
    cmap.set_bad(vs.SURFACE)
    im = ax.imshow(M, aspect="auto", cmap=cmap, vmin=-lim, vmax=lim, interpolation="nearest", rasterized=True)
    ax.set_yticks(range(22))
    ax.set_yticklabels([f"chr{c}" for c in range(1, 23)], fontsize=7)
    ax.tick_params(axis="y", length=0)
    ax.set_xticks(np.arange(0, M.shape[1], 50))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{int(v)} Mb"))
    for sp in ax.spines.values():
        sp.set_visible(False)
    cb = fig.colorbar(im, cax=cax, orientation="horizontal")
    cb.set_ticks([-lim, 0, lim])
    cb.set_ticklabels([f"{lim * 100:.0f}+ points fewer", f"genome-wide {genome:.0%}", f"{lim * 100:.0f}+ points more"])
    cb.ax.tick_params(labelsize=7, length=2)
    cb.outline.set_visible(False)
    cax.text(-0.04, 0.5, f"blank tiles: centromeres and gaps (fewer than {min_snvs} common SNVs)",
             transform=cax.transAxes, ha="right", va="center", fontsize=7, color=vs.INK_2)
    out = keep.assign(dev=keep["frac_yri_only"] - genome)
    summary = {"genome_share_yri_only": genome, "tiles_shown": len(keep), "tiles_total": len(t),
               "sd_dev": out["dev"].std(), "tiles_above_plus10": int((out["dev"] > 0.10).sum()),
               "tiles_below_minus10": int((out["dev"] < -0.10).sum())}
    pd.DataFrame([summary]).to_csv(OUT / "fig02_map_summary.tsv", sep="\t", index=False)
    return genome, summary


def figure2():
    spec, cls = spectrum_counts()
    vs.apply()
    fig = plt.figure(figsize=(vs.DOUBLE, 6.8))
    ax_m = fig.add_axes([0.08, 0.525, 0.9, 0.405])
    cax = fig.add_axes([0.60, 0.455, 0.33, 0.01])
    ax_s = fig.add_axes([0.11, 0.06, 0.87, 0.3])
    genome, summary = genome_map(fig, ax_m, cax)
    fig.text(0.075, 0.985, "Where in the genome the Yoruba sample carries common variants the European sample does not",
             fontsize=8, color=vs.INK, va="top")
    fig.text(0.075, 0.96, "Each tile is 1 Mb: its share of Yoruba-only common SNVs, relative to the genome-wide share",
             fontsize=7.5, color=vs.INK_2, va="top")
    mirrored_spectrum(ax_s, spec, cls)
    for a_, letter in ((ax_m, "a"), (ax_s, "b")):
        bb = a_.get_position()
        fig.text(0.008, bb.y1 + 0.01, letter, fontsize=10, fontweight="bold", color=vs.INK, va="bottom")
    vs.save(fig, "fig02_frequency_landscape")
    return cls


# ----------------------------------------------------------------------------- Figure 3

def chr22_common_both():
    for pop in ("EUR87_r1", "YRI87"):
        out = G / f"chr22_freq_{pop}"
        if not Path(f"{out}.afreq").exists():
            plink("--keep", str(S / f"{pop}.txt"), "--chr", "22", "--freq", "--out", str(out))
    e = pd.read_csv(G / "chr22_freq_EUR87_r1.afreq", sep="\t", usecols=["ID", "ALT_FREQS"])
    y = pd.read_csv(G / "chr22_freq_YRI87.afreq", sep="\t", usecols=["ID", "ALT_FREQS"]).set_index("ID").loc[e["ID"]]
    me = np.minimum(e["ALT_FREQS"].to_numpy(), 1 - e["ALT_FREQS"].to_numpy())
    my = np.minimum(y["ALT_FREQS"].to_numpy(), 1 - y["ALT_FREQS"].to_numpy())
    ids = e["ID"].to_numpy()[(me >= 0.05) & (my >= 0.05)]
    pos = np.array([int(i.split(":")[1]) for i in ids])
    order = np.argsort(pos)
    return ids[order], pos[order]


# Centered on DNAJB7 (GRCh37 chr22:41,255,553 to 41,258,130), the worked example gene in Figure 12,
# so the LD contrast shown here is the one that shapes that gene's models.
DNAJB7_CENTER = (41_255_553 + 41_258_130) // 2


def region_r2(ids, pos, span=250_000, n_max=150, center=DNAJB7_CENTER):
    s = center - span // 2
    sel = np.where((pos >= s) & (pos < s + span))[0]
    sel = sel[np.linspace(0, len(sel) - 1, min(n_max, len(sel))).round().astype(int)]
    snp_file = G / "ld_region_chr22_snps.txt"
    pd.Series(ids[sel]).to_csv(snp_file, index=False, header=False)
    r2 = {}
    for pop in ("EUR87_r1", "YRI87"):
        out = G / f"ld_region_{pop}"
        plink("--keep", str(S / f"{pop}.txt"), "--extract", str(snp_file), "--export", "A", "--out", str(out))
        raw = pd.read_csv(f"{out}.raw", sep="\t")
        start = list(raw.columns).index("PHENOTYPE") + 1
        x = raw.iloc[:, start:].to_numpy(dtype=float)
        r = np.corrcoef(x, rowvar=False)
        r2[pop] = np.nan_to_num(r ** 2)
    return r2, pos[sel]


def ld_triangles(ax, r2, positions):
    n = r2["EUR87_r1"].shape[0]
    i, j = np.triu_indices(n, k=1)
    cx, cy = (i + j) / 2.0, (j - i) / 2.0
    diamond = np.stack([np.stack([cx, cy - 0.5], 1), np.stack([cx + 0.5, cy], 1),
                        np.stack([cx, cy + 0.5], 1), np.stack([cx - 0.5, cy], 1)], axis=1)
    # Ramp starts at near-white so pairs in weak LD recede and the haplotype blocks stand out.
    cmap = plt.matplotlib.colors.LinearSegmentedColormap.from_list("ld", ["#fbfcff"] + vs.SEQ[2:])
    # Linear scale: at n = 87 pairs with no real linkage still show r^2 near 1/n, and a
    # compressive scale turns that sampling noise into texture. Linear keeps only real blocks.
    norm = plt.matplotlib.colors.Normalize(vmin=0, vmax=1)
    for pop, sign in (("EUR87_r1", 1), ("YRI87", -1)):
        verts = diamond.copy()
        verts[:, :, 1] *= sign
        pc = PolyCollection(verts, array=r2[pop][i, j], cmap=cmap, norm=norm, edgecolors="none", rasterized=True)
        ax.add_collection(pc)
    ax.axhline(0, color=vs.SURFACE, lw=1.2)
    half = (n - 1) / 2
    ax.set_xlim(-1, n)
    ax.set_ylim(-half - 3, half + 3)
    ax.set_aspect(0.5)
    ax.axis("off")
    ax.text(-1, half * 0.95, "European sample\n(EUR87, draw 1)", fontsize=8, color=vs.INK, va="top", ha="left", fontweight="bold")
    ax.text(-1, -half * 0.95, "Yoruba sample\n(YRI87)", fontsize=8, color=vs.INK, va="bottom", ha="left", fontweight="bold")
    mb0, mb1 = positions.min() / 1e6, positions.max() / 1e6
    ax.text(n, -half * 0.95, f"around DNAJB7, chr22:{mb0:.2f} to {mb1:.2f} Mb\n{n} SNVs common in both",
            fontsize=7.5, color=vs.INK_2, va="bottom", ha="right")
    return pc


def decay_curves(ax):
    rng = np.random.default_rng(2026)
    ids, _ = chr22_common_both()
    pick = np.sort(rng.choice(ids, size=min(6000, len(ids)), replace=False))
    snp_file = G / "ld_chr22_snps.txt"
    pd.Series(pick).to_csv(snp_file, index=False, header=False)
    edges = np.unique(np.round(np.logspace(np.log10(0.5), np.log10(500), 26), 3))
    curves = {}
    for pop, color, label in (("EUR87_r1", vs.EUR, "European"), ("YRI87", vs.YRI, "Yoruba")):
        out = G / f"ld_chr22_{pop}"
        plink("--keep", str(S / f"{pop}.txt"), "--extract", str(snp_file), "--r2-unphased",
              "--ld-window-kb", "500", "--ld-window", "999999", "--ld-window-r2", "0", "--out", str(out))
        d = pd.read_csv(f"{out}.vcor", sep="\t")
        r2col = [c for c in d.columns if "R2" in c.upper()][0]
        dist = (d["POS_B"] - d["POS_A"]).abs() / 1000
        idx = np.digitize(dist, edges)
        g = pd.DataFrame({"bin": idx, "r2": d[r2col]}).groupby("bin")["r2"].agg(["mean", "count"])
        g = g[(g.index > 0) & (g.index < len(edges))]
        mid = np.sqrt(edges[g.index - 1] * edges[g.index])
        c = pd.DataFrame({"distance_kb": mid, "mean_r2": g["mean"].to_numpy(), "pairs": g["count"].to_numpy()})
        c.to_csv(OUT / f"fig03_ld_decay_{pop}.tsv", sep="\t", index=False)
        curves[pop] = c
        ax.plot(c["distance_kb"], c["mean_r2"], color=color, lw=1.6)
        ax.scatter(c["distance_kb"].iloc[[0]], c["mean_r2"].iloc[[0]], s=14, color=color, zorder=3,
                   edgecolors=vs.SURFACE, linewidths=0.8)
    ax.set_xscale("log")
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    ax.set_xlabel("Distance between SNVs (kb)")
    ax.set_ylabel("Mean r²")
    ax.grid(True, axis="y")
    ax.set_axisbelow(True)
    e1 = curves["EUR87_r1"]; y1 = curves["YRI87"]
    k = int(np.argmin(np.abs(e1["distance_kb"] - 5)))
    ax.annotate("European", (e1["distance_kb"].iloc[k], e1["mean_r2"].iloc[k]), xytext=(5, 4),
                textcoords="offset points", fontsize=7.5, color=vs.INK)
    ax.annotate("Yoruba", (y1["distance_kb"].iloc[k], y1["mean_r2"].iloc[k]), xytext=(-5, -9),
                textcoords="offset points", fontsize=7.5, color=vs.INK, ha="right")
    return curves


def figure3():
    ids, pos = chr22_common_both()
    r2, positions = region_r2(ids, pos)
    vs.apply()
    fig = plt.figure(figsize=(vs.DOUBLE, 2.7))
    ax_t = fig.add_axes([0.02, 0.03, 0.58, 0.93])
    cax = fig.add_axes([0.615, 0.22, 0.011, 0.56])
    ax_d = fig.add_axes([0.775, 0.19, 0.21, 0.72])
    pc = ld_triangles(ax_t, r2, positions)
    cb = fig.colorbar(pc, cax=cax, orientation="vertical")
    cb.set_label("r² between SNV pairs", fontsize=7.5, color=vs.INK_2)
    cb.set_ticks([0, 0.25, 0.5, 0.75, 1.0])
    cb.ax.tick_params(labelsize=7.5)
    cb.outline.set_visible(False)
    curves = decay_curves(ax_d)
    fig.text(0.005, 0.975, "a", fontsize=10, fontweight="bold", color=vs.INK, va="top")
    fig.text(0.695, 0.975, "b", fontsize=10, fontweight="bold", color=vs.INK, va="top")
    vs.save(fig, "fig03_ld_decay")
    summary = pd.DataFrame([{"region_start": int(positions.min()), "region_end": int(positions.max()), "n_snvs": len(positions),
                             "mean_r2_region_EUR87_r1": float(r2["EUR87_r1"][np.triu_indices(len(positions), 1)].mean()),
                             "mean_r2_region_YRI87": float(r2["YRI87"][np.triu_indices(len(positions), 1)].mean())}])
    summary.to_csv(OUT / "fig03_region_summary.tsv", sep="\t", index=False)
    return summary, curves


if __name__ == "__main__":
    union_freqs()
    cls = figure2()
    print(cls.to_string())
    summary, curves = figure3()
    print(summary.to_string(index=False))
    for k, v in curves.items():
        print(k, v.iloc[[0, len(v) // 2, -1]].to_string(index=False))
