"""Build the training sample lists for ibeji.

Matches GEUVADIS expression samples to 1000 Genomes phase 3 genotypes, then
writes PLINK 2 --keep files for the European and Yoruba sets plus five random
European downsamples matched to the Yoruba sample size.
"""
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
OUT = ROOT / "data" / "processed" / "samples"
OUT.mkdir(parents=True, exist_ok=True)

EUR_POPS = {"CEU", "FIN", "GBR", "TSI"}
N_DOWNSAMPLE_REPS = 5

expr_header = pd.read_csv(
    RAW / "GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz", sep="\t", nrows=0
)
expr_samples = list(expr_header.columns[4:])

psam = pd.read_csv(RAW / "phase3_corrected.psam", sep="\t").set_index("#IID")
matched = [s for s in expr_samples if s in psam.index]
missing = [s for s in expr_samples if s not in psam.index]

pops = psam.loc[matched, "Population"]
eur = sorted(pops[pops.isin(EUR_POPS)].index)
yri = sorted(pops[pops == "YRI"].index)


def write_keep(name: str, ids: list[str]) -> None:
    pd.DataFrame({"#IID": ids}).to_csv(OUT / f"{name}.txt", sep="\t", index=False)


write_keep("ALL", sorted(eur + yri))
write_keep("EUR358", eur)
write_keep("YRI87", yri)

rows = []
for rep in range(1, N_DOWNSAMPLE_REPS + 1):
    draw = sorted(pd.Series(eur).sample(n=len(yri), random_state=rep))
    write_keep(f"EUR87_r{rep}", draw)
    comp = pops.loc[draw].value_counts().to_dict()
    rows.append({"set": f"EUR87_r{rep}", "n": len(draw), **comp})

summary = pd.DataFrame(
    [{"set": "EUR358", "n": len(eur), **pops.loc[eur].value_counts().to_dict()},
     {"set": "YRI87", "n": len(yri), "YRI": len(yri)}] + rows
).fillna(0)
summary.to_csv(OUT / "set_composition.tsv", sep="\t", index=False)

print(f"expression samples: {len(expr_samples)}; matched to 1KG phase 3: {len(matched)}")
print(f"missing from phase 3 ({len(missing)}): {', '.join(missing)}")
print(summary.to_string(index=False))
