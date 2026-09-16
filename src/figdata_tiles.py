"""Per-megabase counts of common SNVs by sharing class, for the Figure 2 genome map.

For every autosome and every 1 Mb tile: SNVs with MAF >= 0.01 in both samples, only in
EUR358, and only in YRI87. Streams the two union frequency files in chunks (both list
the same variants in the same order), so memory stays flat.
Output: results/figdata/fig02_tiles.tsv
"""
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
G = ROOT / "data" / "processed" / "geno"
OUT = ROOT / "results" / "figdata"
OUT.mkdir(parents=True, exist_ok=True)

TILE = 1_000_000
N_TILES = 250  # chr1 is 249 Mb
counts = np.zeros((23, N_TILES, 3), dtype=np.int64)  # [chr][tile][both, eur_only, yri_only]

kw = dict(sep="\t", usecols=["#CHROM", "ID", "ALT_FREQS"], chunksize=2_000_000)
ky = dict(sep="\t", usecols=["ID", "ALT_FREQS"], chunksize=2_000_000)
for ce, cy in zip(pd.read_csv(G / "union_freq_EUR358.afreq", **kw), pd.read_csv(G / "union_freq_YRI87.afreq", **ky)):
    assert (ce["ID"].values == cy["ID"].values).all()
    chrom = ce["#CHROM"].to_numpy().astype(int)
    pos = ce["ID"].str.split(":", n=2, expand=True)[1].astype(np.int64).to_numpy()
    fe, fy = ce["ALT_FREQS"].to_numpy(), cy["ALT_FREQS"].to_numpy()
    oe, oy = np.minimum(fe, 1 - fe) >= 0.01, np.minimum(fy, 1 - fy) >= 0.01
    cls = np.select([oe & oy, oe & ~oy, ~oe & oy], [0, 1, 2], default=-1)
    ok = cls >= 0
    np.add.at(counts, (chrom[ok], pos[ok] // TILE, cls[ok]), 1)

rows = []
for c in range(1, 23):
    for t in range(N_TILES):
        b, e, y = counts[c, t]
        if b + e + y:
            rows.append((c, t, b, e, y))
df = pd.DataFrame(rows, columns=["chrom", "tile_mb", "both", "eur_only", "yri_only"])
df["total"] = df[["both", "eur_only", "yri_only"]].sum(axis=1)
df["frac_yri_only"] = df["yri_only"] / df["total"]
df["frac_eur_only"] = df["eur_only"] / df["total"]
df.to_csv(OUT / "fig02_tiles.tsv", sep="\t", index=False)
print(f"tiles: {len(df):,}; SNVs counted: {df['total'].sum():,}")
print(df["frac_yri_only"].describe().round(3).to_string())
