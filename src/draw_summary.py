"""Every number in the paper that depends on how many European draws have been trained.

Usage: python3 src/draw_summary.py

The study trains EUR87 repeatedly with different seeds, so a set of statements in the paper
move as draws accumulate: the yield span, the within-ancestry noise floor and its range, the
cross-ancestry weights term, and the whole of the draw-dependence subsection. This script
computes all of them in one pass from the finished sets, so the manuscript can be updated
from measured values rather than from memory, and rechecked at any time.

It also makes one check explicit rather than leaving it to be noticed. The paper claims the
noise floor of the weights term exceeds the frequency and LD components combined, and says
this does not depend on which control pair or which draw is used. That holds only while the
SMALLEST floor still exceeds the LARGEST frequency plus LD sum. The margin has been narrowing
as draws are added (0.545 against 0.395 at three draws, 0.510 against 0.405 at four), so the
claim is re-tested here and reported as PASS or FAIL.

Conventions match the figures: degenerate genes, where a mixed combination has V = 0 and the
components are undefined, are excluded throughout, exactly as fig10 does.
"""
import itertools
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
MODELS = ROOT / "results" / "models"
DEC = ROOT / "results" / "decomposition"
YRI = "YRI87"


def complete_draws():
    out = []
    for d in sorted(MODELS.glob("EUR87_r*")):
        if len(list(d.glob("chr*.summary.tsv"))) == 22:
            out.append(d.name)
    return out


def yield_row(set_name):
    d = pd.concat([pd.read_csv(f, sep="\t") for f in sorted((MODELS / set_name).glob("chr*.summary.tsv"))],
                  ignore_index=True)
    cv = d[(d.cv_r2 > 0.01) & (d.cv_pval < 0.05)]
    usable = cv[cv.n_model > 0]
    return {"set": set_name, "attempted": len(d), "cleared_cv": len(cv), "usable": len(usable),
            "empty_final": len(cv) - len(usable), "median_r2": float(usable.cv_r2.median())}


def load(path):
    d = pd.read_csv(path, sep="\t")
    return d[~d["degenerate"]]


def shares(path):
    d = load(path)
    t = d[["phi_w", "phi_D", "phi_R"]].abs()
    o = t.div(t.sum(axis=1), axis=0)
    o.columns = ["s_w", "s_D", "s_R"]
    o["gene"] = d.gene.values
    o["dom"] = t.values.argmax(axis=1)
    return o


draws = complete_draws()
print(f"{len(draws)} complete European draws: {', '.join(draws)}\n")

print("=" * 78)
print("YIELD  (paper: the draws span X to Y usable, and A to B clear cross-validation)")
print("=" * 78)
rows = [yield_row(s) for s in draws + [YRI]]
y = pd.DataFrame(rows)
print(y.to_string(index=False))
eur = y[y.set.str.startswith("EUR87")]
yri = y[y.set == YRI].iloc[0]
print(f"\n  usable span across draws: {eur.usable.min():,} to {eur.usable.max():,} "
      f"(width {eur.usable.max() - eur.usable.min()})")
print(f"  Yoruba usable: {yri.usable:,}  -> inside the span: "
      f"{eur.usable.min() <= yri.usable <= eur.usable.max()}")
print(f"  cleared cross-validation, all sets: {min(y.cleared_cv):,} to {max(y.cleared_cv):,}")
frac = [(r.empty_final / r.cleared_cv) for r in y.itertuples()]
print(f"  empty final fit as a share of cleared: {min(frac):.0%} to {max(frac):.0%}")

print("\n" + "=" * 78)
print("NOISE FLOOR  (paper: the floor lies between X and Y)")
print("=" * 78)
floors = []
for a, b in itertools.combinations(draws, 2):
    p = DEC / f"{a}_vs_{b}.tsv"
    if not p.exists():
        print(f"  {a} vs {b}: MISSING, not yet computed")
        continue
    d = load(p)
    floors.append((f"{a}_vs_{b}", len(d), float(d.phi_w.abs().median())))
for name, n, v in floors:
    print(f"  {name:28s} n={n:4d}  median |phi_w| = {v:.3f}")
fv = [v for _, _, v in floors]
print(f"\n  {len(floors)} of {len(draws) * (len(draws) - 1) // 2} possible pairs present")
print(f"  floor range {min(fv):.3f} to {max(fv):.3f}, median {float(np.median(fv)):.3f}")

print("\n" + "=" * 78)
print("CROSS-ANCESTRY COMPONENTS  (paper: weights term lies between X and Y)")
print("=" * 78)
cross = []
for s in draws:
    p = DEC / f"{s}_vs_{YRI}.tsv"
    if not p.exists():
        print(f"  {s}: MISSING")
        continue
    d = load(p)
    w, dd, r = (float(d[c].abs().median()) for c in ("phi_w", "phi_D", "phi_R"))
    cross.append((s, len(d), w, dd, r))
    print(f"  {s:12s} n={len(d):4d}  |phi_w|={w:.3f}  |phi_D|={dd:.3f}  |phi_R|={r:.3f}  D+R={dd + r:.3f}")
cw = [w for _, _, w, _, _ in cross]
sums = [dd + r for _, _, _, dd, r in cross]
print(f"\n  weights term across draws: {min(cw):.3f} to {max(cw):.3f}")
print(f"  frequency plus LD sum:     {min(sums):.3f} to {max(sums):.3f}")

print("\n" + "=" * 78)
print("THE MARGIN CLAIM  (paper: the floor exceeds frequency and LD combined,")
print("                   whichever control pair and whichever draw is used)")
print("=" * 78)
lo_floor, hi_sum = min(fv), max(sums)
print(f"  smallest floor      {lo_floor:.3f}")
print(f"  largest D+R sum     {hi_sum:.3f}")
print(f"  margin              {lo_floor - hi_sum:+.3f}")
print(f"  VERDICT: {'PASS, the claim holds for every pairing' if lo_floor > hi_sum else 'FAIL, the sentence must be rewritten'}")

print("\n" + "=" * 78)
print("DRAW DEPENDENCE  (paper: pooling all pairs gives N comparisons over G genes,")
print("                  the largest component changes for X percent)")
print("=" * 78)
per = {s: shares(DEC / f"{s}_vs_{YRI}.tsv") for s in draws if (DEC / f"{s}_vs_{YRI}.tsv").exists()}
merged = [per[a].merge(per[b], on="gene", suffixes=("_1", "_2")) for a, b in itertools.combinations(per, 2)]
m = pd.concat(merged, ignore_index=True)
M = np.zeros((3, 3), int)
for i in range(3):
    for j in range(3):
        M[i, j] = int(((m.dom_1 == i) & (m.dom_2 == j)).sum())
agree = float((m.dom_1 == m.dom_2).mean())
print(f"  {len(merged)} pairs, {len(m):,} pooled comparisons, {m.gene.nunique()} distinct genes")
print(f"  agreement {agree:.1%}, so the largest component changes for {1 - agree:.1%}")
print(f"  matrix diagonal: weights {M[0, 0]}, frequency {M[1, 1]}, LD {M[2, 2]}; switched {int(M.sum() - np.trace(M))}")
print(f"  weights-share Spearman {stats.spearmanr(m.s_w_1, m.s_w_2)[0]:.2f}, "
      f"median absolute change {float((m.s_w_1 - m.s_w_2).abs().median()):.2f}")
print("  per-pair agreement:")
for (a, b), mm in zip(itertools.combinations(per, 2), merged):
    print(f"    {a} vs {b}: {len(mm)} genes, {(mm.dom_1 == mm.dom_2).mean():.1%}")

print("\n  median share of the total, per draw (Figure 13 panel c):")
for s in per:
    d = per[s]
    print(f"    {s:12s} weights {d.s_w.median():.3f}  frequency {d.s_D.median():.3f}  LD {d.s_R.median():.3f}")
sw = [per[s].s_w.median() for s in per]
sd = [per[s].s_D.median() for s in per]
sr = [per[s].s_R.median() for s in per]
print(f"    ranges: weights {min(sw):.3f} to {max(sw):.3f}, "
      f"frequency {min(sd):.3f} to {max(sd):.3f}, LD {min(sr):.3f} to {max(sr):.3f}")
