"""Build the ancestry-specific association table that Figure 18 draws.

Usage: python3 src/figdata_specific_hits.py [eur_set] [yri_set]
Default: EUR87_r1 YRI87

Every gene and trait association that is Bonferroni-significant with one model set and not
the other is classified into exactly one of:

  no usable model in the other ancestry
      The other ancestry never produced a usable model for that gene, so weights,
      frequencies and LD cannot be compared at all.
  weights / allele frequency / linkage disequilibrium
      Both ancestries modelled the gene, and this is the largest absolute component of the
      Shapley decomposition of the gap.
  decomposition undefined
      Both models exist, but some mixed combination has V = 0 because all weight sits on
      SNVs monomorphic in the other sample, so the gene cannot be placed on the log scale.
      This is the degenerate case the Methods say is reported separately.

Bonferroni is taken within each trait and model set over the genes actually TESTED, that is
those with a defined p, which is the same rule fig16_mirror_manhattan.py uses. Genes with no
usable SNV in a GWAS carry p = NA from 06_twas.R and must not inflate the denominator.

Writes results/figdata/fig18_specific_hits.tsv.

This script did not exist until 2026-09-16: the table was previously built by hand and was
wrong in two ways. It held 129 of the 134 associations, and it classified one degenerate
gene (ENSG00000132394.6, schizophrenia, European models, all three components NaN) as
weights-dominated, which is not a judgement call because there is no component to be
largest. Both errors reached the manuscript.
"""
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
TWAS = ROOT / "results" / "twas"
MODELS = ROOT / "results" / "models"
DEC = ROOT / "results" / "decomposition"
OUT = ROOT / "results" / "figdata"

EUR_SET = sys.argv[1] if len(sys.argv) > 1 else "EUR87_r1"
YRI_SET = sys.argv[2] if len(sys.argv) > 2 else "YRI87"

TRAIT_NAME = {"ASD": "Autism", "SCZ": "Schizophrenia", "BIP": "Bipolar disorder",
              "MDD": "Major depression", "PTSD": "PTSD"}
COMPONENT = {"phi_w": "weights", "phi_D": "allele frequency", "phi_R": "linkage disequilibrium"}
PHI = list(COMPONENT)


def significant(set_name):
    """Bonferroni-significant genes per trait, with the denominator over tested genes."""
    d = pd.read_csv(TWAS / set_name / "twas_all_traits.tsv", sep="\t")
    out = {}
    for trait, g in d.groupby("trait"):
        tested = g[g["p"].notna()]
        out[trait] = set(tested.loc[tested["p"] < 0.05 / len(tested), "gene"])
    return out


def usable(set_name):
    files = sorted((MODELS / set_name).glob("chr*.summary.tsv"))
    if len(files) < 22:
        raise SystemExit(f"{set_name} has only {len(files)}/22 chromosomes; not ready")
    d = pd.concat([pd.read_csv(f, sep="\t") for f in files], ignore_index=True)
    return set(d[(d.cv_r2 > 0.01) & (d.cv_pval < 0.05) & (d.n_model > 0)].gene)


def require_ready(*sets):
    """Fail with a clear message, not a traceback, when an input is not ready yet.

    The check has to come before anything is read. The 22-chromosome guard inside usable()
    never fired when this was first written, because significant() runs first and died on a
    missing TWAS file instead.
    """
    for s in sets:
        n = len(list((MODELS / s).glob("chr*.summary.tsv")))
        if n < 22:
            raise SystemExit(f"{s} has only {n}/22 chromosomes; not ready")
        if not (TWAS / s / "twas_all_traits.tsv").exists():
            raise SystemExit(f"{s} has no TWAS output yet; run src/06_twas.R {s} first")
    path = DEC / f"{EUR_SET}_vs_{YRI_SET}.tsv"
    if not path.exists():
        raise SystemExit(f"{path.name} not found; run src/05_decompose.R {EUR_SET} {YRI_SET} first")


require_ready(EUR_SET, YRI_SET)
sig_e, sig_y = significant(EUR_SET), significant(YRI_SET)
use_e, use_y = usable(EUR_SET), usable(YRI_SET)
dec = pd.read_csv(DEC / f"{EUR_SET}_vs_{YRI_SET}.tsv", sep="\t").set_index("gene")

records = []
for code in ["SCZ", "BIP", "MDD", "PTSD", "ASD"]:
    for gene, found_with, usable_other in (
            [(g, "European models", use_y) for g in sorted(sig_e[code] - sig_y[code])]
            + [(g, "Yoruba models", use_e) for g in sorted(sig_y[code] - sig_e[code])]):
        if gene not in usable_other:
            klass = "no usable model in the other ancestry"
        elif gene in dec.index and dec.loc[gene, PHI].notna().all():
            klass = COMPONENT[dec.loc[gene, PHI].abs().idxmax()]
        else:
            klass = "decomposition undefined"
        records.append({"trait": TRAIT_NAME[code], "found_with": found_with,
                        "gene": gene, "klass": klass})

out = pd.DataFrame(records)
OUT.mkdir(parents=True, exist_ok=True)
out.to_csv(OUT / "fig18_specific_hits.tsv", sep="\t", index=False)

print(f"{len(out)} ancestry-specific associations ({EUR_SET} vs {YRI_SET})")
print(out.klass.value_counts().to_string())
both = out[out.klass != "no usable model in the other ancestry"]
print(f"\nboth ancestries modelled the gene: {len(both)}")
print(f"\n{pd.crosstab(out.trait, out.found_with).to_string()}")
