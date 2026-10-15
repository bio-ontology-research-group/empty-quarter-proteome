#!/usr/bin/env python3
"""Step 4 -- environmental-gradient test for the metal-resistance layer.

The PI's exploratory XRF scan (06_resistome_safety/amr_metal/metal_xrf_corr.tsv, old
15,469 basis) found arsenic-gene load tracking the salinity proxy (r=+0.26) -- a
non-cognate association, since the XRF panel measures neither As nor Hg. Here the
same question is asked on the CURRENT basis with the section-2.3 method: per-sample
community carriage (occurrence and whole-community abundance) of STRESS:ARSENIC /
STRESS:MERCURY (+ any_AMR as a control trait) against the salt and climate axes,
partial Spearman controlling ranked lat+lon (phase_e_partial_osmolyte.py precedent),
BH-FDR within block. Ported from osmolyte_clean/scripts/osmo_5_gradients.py.

In:  results/community_by_sample.tsv (amr_2_community.py), repo env_matrix.tsv
Out: results/gradient_tests.tsv

Run: python3 scripts/amr_4_gradients.py
"""
import csv, collections
import numpy as np
from scipy.stats import rankdata, spearmanr, pearsonr

M    = "/ibex/scratch/projects/c2014/rob/metagenomics"
B    = f"{M}/06_metabolism/resistome_clean"
RES  = f"{B}/results"
REPO = "/ibex/project/c2014/alelopezv/repos/empty-quarter-mags"
ENV  = f"{REPO}/results/darkmatter_ecology/env_matrix.tsv"

SALT    = ["salinity_proxy", "S", "Cl", "Na", "Ca"]
CLIMATE = ["aridity_index", "precip", "temp"]
TRAITS  = ["STRESS:ARSENIC", "STRESS:MERCURY", "any_AMR"]
METRICS = ["occurrence_pct", "abund_whole_pct"]

env = {r["sample"]: r for r in csv.DictReader(open(ENV), delimiter="\t")}
print(f"[grad] env matrix: {len(env)} samples")

def envval(s, ax):
    v = env.get(s, {}).get(ax, "")
    if v in ("", "NA"):
        return None
    if ax in SALT and env[s].get("has_xrf") != "1":
        return None
    try:
        return float(v)
    except ValueError:
        return None

d = collections.defaultdict(dict)
for r in csv.DictReader(open(f"{RES}/community_by_sample.tsv"), delimiter="\t"):
    for m in METRICS:
        d[(r["sample"], m)][r["trait"]] = float(r[m])
samples = sorted({s for s, m in d})
print(f"[grad] samples with traits: {len(samples)}")

def partial_spearman(x, y, lat, lon):
    rx, ry = rankdata(x), rankdata(y)
    C = np.column_stack([np.ones(len(x)), rankdata(lat), rankdata(lon)])
    bx, *_ = np.linalg.lstsq(C, rx, rcond=None)
    by, *_ = np.linalg.lstsq(C, ry, rcond=None)
    return pearsonr(rx - C @ bx, ry - C @ by)

def bh(pvals):
    n = len(pvals)
    order = sorted(range(n), key=lambda i: pvals[i])
    q, prev = [0.0] * n, 1.0
    for rank in range(n - 1, -1, -1):
        i = order[rank]
        prev = min(prev, pvals[i] * n / (rank + 1))
        q[i] = prev
    return q

rows = []
for block, axes in (("salt", SALT), ("climate", CLIMATE)):
    stash = []
    for ax in axes:
        sub = [s for s in samples if envval(s, ax) is not None
               and env.get(s, {}).get("lat") not in ("", None)
               and env.get(s, {}).get("lon") not in ("", None)]
        if len(sub) < 30:
            continue
        y   = np.array([envval(s, ax) for s in sub])
        lat = np.array([float(env[s]["lat"]) for s in sub])
        lon = np.array([float(env[s]["lon"]) for s in sub])
        for t in TRAITS:
            for m in METRICS:
                x = np.array([d[(s, m)][t] for s in sub])
                if np.all(x == x[0]):
                    continue
                rho_m, p_m = spearmanr(x, y)
                rho_p, p_p = partial_spearman(x, y, lat, lon)
                stash.append(dict(block=block, axis=ax, n=len(sub), trait=t, metric=m,
                                  rho_marginal=rho_m, p_marginal=p_m,
                                  rho_partial=rho_p, p_partial=p_p))
    for r, q in zip(stash, bh([r["p_partial"] for r in stash])):
        r["q_partial"] = q
    rows += stash

cols = ["block", "axis", "n", "trait", "metric",
        "rho_marginal", "p_marginal", "rho_partial", "p_partial", "q_partial"]
with open(f"{RES}/gradient_tests.tsv", "w") as o:
    o.write("\t".join(cols) + "\n")
    for r in rows:
        o.write("\t".join(f"{r[c]:.4g}" if isinstance(r[c], float) else str(r[c])
                          for c in cols) + "\n")
print(f"[grad] wrote {RES}/gradient_tests.tsv ({len(rows)} tests)")
print("\nFDR-significant (q<0.05), partial:")
for r in rows:
    if r["q_partial"] < 0.05:
        print(f"  {r['block']:8s} {r['axis']:15s} {r['trait']:16s} {r['metric']:16s} "
              f"rho_p={r['rho_partial']:+.3f} p={r['p_partial']:.2g} q={r['q_partial']:.2g} n={r['n']}")
