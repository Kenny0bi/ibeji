"""Harmonize the five GWAS to the ibeji genotype variant IDs.

Every GWAS is GRCh37, the same build as 1000 Genomes phase 3. For each trait:
  - compute a z-score for the reported effect allele
      ASD: log(OR) / SE      SCZ, BIP, MDD: BETA / SE      PTSD: Z as reported
  - match on chromosome, position and alleles to the geuv445 variants (chr:pos:ref:alt)
  - flip the sign when the effect allele is REF, so every z is for the ALT allele
  - allow strand-complement matches for non-palindromic SNPs
  - drop palindromic SNPs (A/T, C/G), where strand cannot be resolved from alleles alone

Output: data/processed/gwas/<trait>.tsv.gz with columns varID, z
and a harmonization report.
"""
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
OUT = ROOT / "data" / "processed" / "gwas"
OUT.mkdir(parents=True, exist_ok=True)

COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}

TRAITS = {
    "ASD": dict(file="iPSYCH-PGC_ASD_Nov2017.gz", chrom="CHR", pos="BP", ea="A1", oa="A2", effect=("OR", "SE")),
    "SCZ": dict(file="PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv.gz", chrom="CHROM", pos="POS", ea="A1", oa="A2", effect=("BETA", "SE")),
    "BIP": dict(file="pgc-bip2021-all.vcf.tsv.gz", chrom="#CHROM", pos="POS", ea="A1", oa="A2", effect=("BETA", "SE")),
    "MDD": dict(file="pgc-mdd2025_no23andMe_eur_v3-49-24-11.tsv.gz", chrom="#CHROM", pos="POS", ea="EA", oa="NEA", effect=("BETA", "SE")),
    "PTSD": dict(file="eur_ptsd_pcs_v4_aug3_2021.vcf.gz", chrom="#CHROM", pos="POS", ea="A1", oa="A2", effect=("Z",)),
}


def load_reference() -> pd.DataFrame:
    """Read geuv445.pvar compactly: skip the '##' metadata block, keep chrom/pos/ref/alt.

    Variant IDs are chr:pos:ref:alt by construction, so they are rebuilt only for matched rows.
    """
    path = ROOT / "data" / "processed" / "geno" / "geuv445.pvar"
    with open(path) as fh:
        n_meta = 0
        for line in fh:
            if line.startswith("##"):
                n_meta += 1
            else:
                break
    pvar = pd.read_csv(path, sep="\t", skiprows=n_meta, usecols=["#CHROM", "POS", "REF", "ALT"],
                       dtype={"#CHROM": "int8", "POS": "int32", "REF": "category", "ALT": "category"})
    return pvar.rename(columns={"#CHROM": "chrom", "POS": "pos", "REF": "ref", "ALT": "alt"})


def count_meta_lines(path: Path) -> int:
    """PGC VCF-style tsv files open with '##' metadata lines before the header."""
    import gzip
    n = 0
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("##"):
                n += 1
            else:
                break
    return n


def harmonize(trait: str, spec: dict, ref: pd.DataFrame) -> dict:
    cols = [spec["chrom"], spec["pos"], spec["ea"], spec["oa"], *spec["effect"]]
    path = RAW / spec["file"]
    g = pd.read_csv(path, sep="\t", comment=None, skiprows=count_meta_lines(path), usecols=cols,
                    dtype={spec["chrom"]: str}, low_memory=False)
    g = g.rename(columns={spec["chrom"]: "chrom", spec["pos"]: "pos", spec["ea"]: "ea", spec["oa"]: "oa"})
    n_raw = len(g)
    g["chrom"] = pd.to_numeric(g["chrom"].str.replace("chr", "", regex=False), errors="coerce")
    g = g[g["chrom"].between(1, 22)]
    g["chrom"] = g["chrom"].astype("int8")
    g["pos"] = g["pos"].astype("int32")
    g["ea"] = g["ea"].str.upper()
    g["oa"] = g["oa"].str.upper()
    g = g[g["ea"].isin(COMP) & g["oa"].isin(COMP)]

    if spec["effect"] == ("OR", "SE"):
        g["z"] = np.log(g["OR"]) / g["SE"]
    elif spec["effect"] == ("BETA", "SE"):
        g["z"] = g["BETA"] / g["SE"]
    else:
        g["z"] = g["Z"]
    g = g[np.isfinite(g["z"])]

    palindromic = g["ea"].map(COMP) == g["oa"]
    n_pal = int(palindromic.sum())
    g = g[~palindromic]

    m = g.merge(ref, on=["chrom", "pos"], how="inner")
    m["ref"] = m["ref"].astype(str)
    m["alt"] = m["alt"].astype(str)
    m["varID"] = m["chrom"].astype(str) + ":" + m["pos"].astype(str) + ":" + m["ref"] + ":" + m["alt"]
    same = (m["ea"] == m["alt"]) & (m["oa"] == m["ref"])
    swap = (m["ea"] == m["ref"]) & (m["oa"] == m["alt"])
    same_c = (m["ea"].map(COMP) == m["alt"]) & (m["oa"].map(COMP) == m["ref"])
    swap_c = (m["ea"].map(COMP) == m["ref"]) & (m["oa"].map(COMP) == m["alt"])
    m["z_alt"] = np.select([same | same_c, swap | swap_c], [m["z"], -m["z"]], default=np.nan)
    m = m[np.isfinite(m["z_alt"])]
    m = m.drop_duplicates("varID", keep=False)

    m[["varID", "z_alt"]].rename(columns={"z_alt": "z"}).to_csv(OUT / f"{trait}.tsv.gz", sep="\t", index=False)
    row = dict(trait=trait, rows_raw=n_raw, palindromic_dropped=n_pal, matched=len(m),
                flipped=int((swap | swap_c).loc[m.index].sum()),
                strand_complement=int((same_c | swap_c).loc[m.index].sum()),
                reference_variants=len(ref), coverage=len(m) / len(ref))
    pd.DataFrame([row]).to_csv(OUT / f"{trait}.report.tsv", sep="\t", index=False)
    print(row, flush=True)
    return row


if __name__ == "__main__":
    reference = load_reference()
    for t, s in TRAITS.items():
        if (OUT / f"{t}.report.tsv").exists():
            print(f"{t}: already harmonized, skipping", flush=True)
            continue
        harmonize(t, s, reference)
    reports = [pd.read_csv(OUT / f"{t}.report.tsv", sep="\t") for t in TRAITS if (OUT / f"{t}.report.tsv").exists()]
    report = pd.concat(reports)
    report.to_csv(OUT / "harmonization_report.tsv", sep="\t", index=False)
    print(report.to_string(index=False))
