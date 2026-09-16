"""Build the portability table that Figure 16 draws (script name fig14_portability.py).

Usage: python3 src/figdata_portability.py [eur_set] [yri_set]
Default: EUR87_r1 YRI87

The table is the per-gene decomposition joined to each model's accuracy at home, plus the
quantities derived from the two:

  home_E, home_Y   cross-validated R2 of the European and Yoruba model for that gene
  loss_E, loss_Y   home accuracy minus accuracy in the other population
  keep_E, keep_Y   the share of home accuracy that survives the crossing

loss_E and loss_Y are written because the table is the record of what was computed, but they
are NOT used as an outcome anywhere in the study. Home accuracy is a nested cross-validation
estimate while the cross-population figure comes from the final model fitted on everyone, so
their difference is not a like-for-like measurement; it is negative for 28.5% of genes in one
direction and 47.6% in the other. The Methods say so explicitly.

Writes results/figdata/fig14_portability.tsv.

This table existed before this script did: it was built by hand on 2026-09-15 with no
generator. It was checked before this script was written and reproduces exactly, so this is
a reproducibility fix and not a correction.
"""
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
MODELS = ROOT / "results" / "models"
DEC = ROOT / "results" / "decomposition"
OUT = ROOT / "results" / "figdata"

EUR_SET = sys.argv[1] if len(sys.argv) > 1 else "EUR87_r1"
YRI_SET = sys.argv[2] if len(sys.argv) > 2 else "YRI87"


def home_accuracy(set_name):
    files = sorted((MODELS / set_name).glob("chr*.summary.tsv"))
    if len(files) < 22:
        raise SystemExit(f"{set_name} has only {len(files)}/22 chromosomes; not ready")
    d = pd.concat([pd.read_csv(f, sep="\t") for f in files], ignore_index=True)
    return d.set_index("gene").cv_r2


path = DEC / f"{EUR_SET}_vs_{YRI_SET}.tsv"
if not path.exists():
    raise SystemExit(f"{path.name} not found; run src/05_decompose.R {EUR_SET} {YRI_SET} first")

dec = pd.read_csv(path, sep="\t")
home_e, home_y = home_accuracy(EUR_SET), home_accuracy(YRI_SET)

out = dec.copy()
out["home_E"] = out.gene.map(home_e)
out["home_Y"] = out.gene.map(home_y)
out["loss_E"] = out.home_E - out.r2_eur_model_in_yri
out["loss_Y"] = out.home_Y - out.r2_yri_model_in_eur
out["keep_E"] = out.r2_eur_model_in_yri / out.home_E
out["keep_Y"] = out.r2_yri_model_in_eur / out.home_Y

OUT.mkdir(parents=True, exist_ok=True)
out.to_csv(OUT / "fig14_portability.tsv", sep="\t", index=False)

print(f"{len(out)} genes decomposed in {EUR_SET} against {YRI_SET}")
print(f"median home R2: European {out.home_E.median():.4f}, Yoruba {out.home_Y.median():.4f}")
top_e = out.nlargest(max(1, len(out) // 10), "home_E")
top_y = out.nlargest(max(1, len(out) // 10), "home_Y")
print(f"top tenth by home accuracy, median kept across the crossing: "
      f"European {top_e.keep_E.median():.3f}, Yoruba {top_y.keep_Y.median():.3f}")
