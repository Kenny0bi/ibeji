"""Build the model-overlap tables that Figures 8 and 9 draw.

Usage: python3 src/figdata_model_overlap.py [eur_set] [yri_set]
Default: EUR87_r1 YRI87

Both figures rest on the same computation, the intersection of the two models' predictors
for every gene usable in both sets, so one script writes both tables.

  results/figdata/fig08_shared_weights.tsv
      One row per SNV selected by BOTH models for the same gene, with the weight each model
      gave it: gene, varID, weight_E, weight_Y.
  results/figdata/fig09_shared_predictors.tsv
      One row per gene: how many SNVs each model used, how many they share, how many are
      private to each, and the Jaccard index of the two predictor sets.

A gene is usable when its cross-validated R2 exceeds 0.01 with p < 0.05 and its final fit
retained at least one non-zero weight, which is the same rule used everywhere else.

These tables existed before this script did: they were built by hand on 2026-09-15 and had
no generator, so the claim that every figure comes from a named script was false for them.
They were checked before this script was written and both reproduce exactly, so this is a
reproducibility fix and not a correction. Compare with figdata_specific_hits.py, where the
hand-built table was also wrong.
"""
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
MODELS = ROOT / "results" / "models"
OUT = ROOT / "results" / "figdata"

EUR_SET = sys.argv[1] if len(sys.argv) > 1 else "EUR87_r1"
YRI_SET = sys.argv[2] if len(sys.argv) > 2 else "YRI87"


def require_ready(*sets):
    for s in sets:
        n = len(list((MODELS / s).glob("chr*.summary.tsv")))
        if n < 22:
            raise SystemExit(f"{s} has only {n}/22 chromosomes; not ready")


def usable(set_name):
    d = pd.concat([pd.read_csv(f, sep="\t") for f in sorted((MODELS / set_name).glob("chr*.summary.tsv"))],
                  ignore_index=True)
    return set(d[(d.cv_r2 > 0.01) & (d.cv_pval < 0.05) & (d.n_model > 0)].gene)


def weights(set_name):
    return pd.concat([pd.read_csv(f, sep="\t") for f in sorted((MODELS / set_name).glob("chr*.weights.tsv"))],
                     ignore_index=True)


require_ready(EUR_SET, YRI_SET)
both = usable(EUR_SET) & usable(YRI_SET)
w_e = weights(EUR_SET)
w_y = weights(YRI_SET)
w_e = w_e[w_e.gene.isin(both)]
w_y = w_y[w_y.gene.isin(both)]

# Figure 8: the SNVs both models chose, with each model's weight on them.
shared = w_e.merge(w_y, on=["gene", "varID"], suffixes=("_E", "_Y"))
shared = shared.rename(columns={"weight_E": "weight_E", "weight_Y": "weight_Y"})
shared = shared[["gene", "varID", "weight_E", "weight_Y"]]

# Figure 9: per-gene predictor counts. A gene with a usable model in both sets always has at
# least one SNV in each, so the union is never empty, but guard the division anyway.
by_e = w_e.groupby("gene").varID.apply(set)
by_y = w_y.groupby("gene").varID.apply(set)
rows = []
# sorted(), not the set itself. Iterating a set of strings follows hash order, which Python
# randomises per process, so this wrote a different byte stream on every run. The table's
# values were right either way, but a generator whose output is not reproducible does not
# make the figure reproducible, which is the whole point of writing it.
for gene in sorted(both):
    e, y = by_e.get(gene, set()), by_y.get(gene, set())
    n_shared, n_union = len(e & y), len(e | y)
    rows.append({"gene": gene, "n_eur": len(e), "n_yri": len(y), "shared": n_shared,
                 "eur_only": len(e) - n_shared, "yri_only": len(y) - n_shared,
                 "jaccard": n_shared / n_union if n_union else 0.0})
predictors = pd.DataFrame(rows)

OUT.mkdir(parents=True, exist_ok=True)
shared.to_csv(OUT / "fig08_shared_weights.tsv", sep="\t", index=False)
predictors.to_csv(OUT / "fig09_shared_predictors.tsv", sep="\t", index=False)

print(f"{len(both)} genes usable in both {EUR_SET} and {YRI_SET}")
print(f"fig08_shared_weights.tsv: {len(shared):,} SNV and gene pairs chosen by both models")
print(f"fig09_shared_predictors.tsv: {len(predictors)} genes, "
      f"{int(predictors.n_eur.sum()):,} European and {int(predictors.n_yri.sum()):,} Yoruba SNVs, "
      f"{int((predictors.shared == 0).sum())} genes sharing none "
      f"({100 * (predictors.shared == 0).mean():.1f}%)")
