"""Validate the ibeji TWAS implementation against the official S-PrediXcan.

Usage: python3 src/09_validate_twas.py <set> <chromosomes> <trait>
  e.g. python3 src/09_validate_twas.py YRI87 22 SCZ

Needs: src/09_export_predictdb.R output for the set/chromosomes, and
src/06_twas.R output (results/twas/<set>/twas_chr<..>.tsv).
Builds a PredictDB sqlite file, runs MetaXcan's SPrediXcan.py on the harmonized
GWAS, and compares gene z-scores with our implementation.
"""
import sqlite3
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
set_name, chroms, trait = sys.argv[1], sys.argv[2], sys.argv[3]
tag = "_".join(chroms.split(","))
vdir = ROOT / "results" / "validation" / f"{set_name}_chr{tag}"

db = vdir / "model.db"
if db.exists():
    db.unlink()
con = sqlite3.connect(db)
pd.read_csv(vdir / "weights.tsv", sep="\t").to_sql("weights", con, index=False)
pd.read_csv(vdir / "extra.tsv", sep="\t").to_sql("extra", con, index=False)
con.execute("CREATE INDEX weights_rsid ON weights (rsid)")
con.execute("CREATE INDEX weights_gene ON weights (gene)")
con.commit()
con.close()

model_snps = set(pd.read_csv(vdir / "weights.tsv", sep="\t", usecols=["rsid"])["rsid"])
g = pd.read_csv(ROOT / "data" / "processed" / "gwas" / f"{trait}.tsv.gz", sep="\t")
g = g[g["varID"].isin(model_snps)].copy()
parts = g["varID"].str.split(":", expand=True)
g["A2"] = parts[2]
g["A1"] = parts[3]
gwas_path = vdir / f"gwas_{trait}.tsv.gz"
g.rename(columns={"varID": "SNP", "z": "Z"})[["SNP", "A1", "A2", "Z"]].to_csv(gwas_path, sep="\t", index=False)

out = vdir / f"spredixcan_{trait}.csv"
cmd = [sys.executable, str(ROOT / "tools" / "MetaXcan" / "software" / "SPrediXcan.py"),
       "--model_db_path", str(db), "--covariance", str(vdir / "covariance.txt.gz"),
       "--gwas_file", str(gwas_path), "--snp_column", "SNP", "--effect_allele_column", "A1",
       "--non_effect_allele_column", "A2", "--zscore_column", "Z", "--keep_non_rsid",
       "--output_file", str(out), "--throw"]
res = subprocess.run(cmd, capture_output=True, text=True)
if res.returncode != 0:
    print(res.stdout[-2000:])
    print(res.stderr[-4000:])
    sys.exit(res.returncode)

theirs = pd.read_csv(out)
ours = pd.read_csv(ROOT / "results" / "twas" / set_name / f"twas_chr{tag}.tsv", sep="\t")
ours = ours[ours["trait"] == trait]
m = theirs.merge(ours, on="gene", suffixes=("_spredixcan", "_ours")).dropna(subset=["zscore", "z"])
r = np.corrcoef(m["zscore"], m["z"])[0, 1]
diff = (m["zscore"] - m["z"]).abs()
report = pd.DataFrame([{
    "set": set_name, "chromosomes": chroms, "trait": trait, "genes_compared": len(m),
    "pearson_r": r, "max_abs_diff": diff.max(), "median_abs_diff": diff.median(),
    "n_snps_used_match": float(np.mean(m["n_snps_used"] == m["n_used"])),
}])
report.to_csv(vdir / f"validation_{trait}.tsv", sep="\t", index=False)
print(report.to_string(index=False))
worst = m.assign(diff=diff).sort_values("diff", ascending=False).head(5)
print(worst[["gene", "zscore", "z", "n_snps_used", "n_used", "diff"]].to_string(index=False))
