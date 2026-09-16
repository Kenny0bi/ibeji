"""Export the real numbers the Manim animation shows, and prove they match the analysis.

For DNAJB7 (EUR87_r1 vs YRI87) this computes log V for all 8 mixes of weights, allele
frequencies and LD, recomputes the three Shapley components from those 8 values, and
checks them against results/decomposition/EUR87_r1_vs_YRI87.tsv. The animation reads only
the JSON written here, so it cannot show a number the analysis did not produce.
Uses the dosage exports written by src/fig12_locus.py.
"""
import itertools
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
LOC = ROOT / "results" / "figdata" / "locus"
GENE, NAME = "ENSG00000172404.4", "DNAJB7"

detail = pd.read_csv(LOC / f"{NAME}_snps_detail.tsv", sep="\t")
snps = detail["varID"].tolist()


def dosages(pop):
    raw = pd.read_csv(LOC / f"{NAME}_{pop}.raw", sep="\t")
    x = raw.iloc[:, list(raw.columns).index("PHENOTYPE") + 1:]
    x.columns = [c.rsplit("_", 1)[0] for c in x.columns]
    return x[snps].to_numpy(dtype=float)


def ingredients(x):
    p = x.mean(axis=0) / 2
    d = 2 * p * (1 - p)
    with np.errstate(invalid="ignore", divide="ignore"):
        r = np.corrcoef(x, rowvar=False)
    r = np.nan_to_num(r)
    np.fill_diagonal(r, 1.0)
    return d, r


dE, RE = ingredients(dosages("EUR358"))
dY, RY = ingredients(dosages("YRI87"))
wE = detail["w_EUR"].to_numpy()
wY = detail["w_YRI"].to_numpy()


def logV(w, d, R):
    s = np.sqrt(d) * w
    return math.log(float(s @ R @ s))


players = ["w", "D", "R"]
vertex = {}
for bits in itertools.product([0, 1], repeat=3):
    w = wY if bits[0] else wE
    d = dY if bits[1] else dE
    R = RY if bits[2] else RE
    vertex["".join(map(str, bits))] = logV(w, d, R)

# Shapley value for each player: average marginal change over the 6 orders.
orders = list(itertools.permutations(range(3)))
paths = []
phi = np.zeros(3)
for order in orders:
    state = [0, 0, 0]
    steps = []
    for k in order:
        before = vertex["".join(map(str, state))]
        state[k] = 1
        after = vertex["".join(map(str, state))]
        steps.append({"player": players[k], "change": after - before})
        phi[k] += (after - before) / len(orders)
    paths.append({"order": [players[k] for k in order], "steps": steps})

dec = pd.read_csv(ROOT / "results" / "decomposition" / "EUR87_r1_vs_YRI87.tsv", sep="\t").set_index("gene").loc[GENE]
check = {"phi_w": abs(phi[0] - dec["phi_w"]), "phi_D": abs(phi[1] - dec["phi_D"]),
         "phi_R": abs(phi[2] - dec["phi_R"]), "delta": abs((vertex["111"] - vertex["000"]) - dec["delta"])}
worst = max(check.values())
if worst > 1e-6:
    raise SystemExit(f"animation numbers disagree with the decomposition output: {check}")

out = {
    "gene": NAME, "eur_set": "EUR87_r1", "yri_set": "YRI87", "n_snps": len(snps),
    "n_nonzero_eur": int((wE != 0).sum()), "n_nonzero_yri": int((wY != 0).sum()),
    "vertex_logV": vertex, "delta": vertex["111"] - vertex["000"],
    "phi": {"w": phi[0], "D": phi[1], "R": phi[2]}, "paths": paths,
    "weights_eur": wE.tolist(), "weights_yri": wY.tolist(),
    "max_abs_difference_vs_decomposition": worst,
}
(ROOT / "animation" / "dnajb7_numbers.json").write_text(json.dumps(out, indent=1))
print(f"8 vertices exported; delta {out['delta']:+.3f}; phi w {phi[0]:+.3f} D {phi[1]:+.3f} R {phi[2]:+.3f}; "
      f"max difference vs decomposition {worst:.2e}")
for p in paths:
    print(" -> ".join(f"{s['player']} {s['change']:+.3f}" for s in p["steps"]))
