#!/usr/bin/env python3
"""Step 5 -- environmental-gradient tests for BOTH §2.3 trait sets (repair + osmolytes).

Compartment (Surface/Deep) tested sun exposure. It never tested salt or aridity, which vary
BETWEEN sites. Here: per-sample community trait values vs continuous environmental axes.

  axes (from the paper repo env_matrix, XRF gate has_xrf):
    salt block   : salinity_proxy, S, Cl, Na, Ca          (XRF samples, n~251)
    climate block: aridity_index, precip, temp            (n~456-547)
  method: partial Spearman controlling ranked lat+lon (the PI's precedent,
    code/phase_e_partial_osmolyte.py / PARTIAL.md), BH-FDR within each block.
  traits: per-sample occurrence_pct and abund_whole_pct for every pathway/group of
    osmolyte_clean + uv_survival_clean, PLUS two repertoire-depth summaries:
      osmo_breadth   = community mean # parallel osmolyte routes (0-4)
                       (= sum of the 4 route prevalences; identity: mean of a sum of
                        indicators = sum of the means)
      repair_breadth = community mean # repair systems (0-8), analogously.

Out: results/gradient_tests.tsv (+ console summary of FDR-significant rows).
"""
import csv, os, math, collections
import numpy as np
from scipy.stats import rankdata, spearmanr, pearsonr, t as tdist

M    = "/ibex/scratch/projects/c2014/rob/metagenomics"
REPO = "/ibex/project/c2014/alelopezv/repos/empty-quarter-mags"
OSMO = f"{M}/06_metabolism/osmolyte_clean/results"
UV   = f"{M}/06_metabolism/uv_survival_clean/results"
ENV  = f"{REPO}/results/darkmatter_ecology/env_matrix.tsv"

SALT    = ["salinity_proxy", "S", "Cl", "Na", "Ca"]
CLIMATE = ["aridity_index", "precip", "temp"]
OSMO_ROUTES  = ["trehalose_any", "ectoine", "glucosylglycerol", "mannosylglycerate"]
REPAIR_GROUPS = ["base_excision_repair", "nucleotide_excision_repair", "SOS_translesion",
                 "dna_protection_oxidative", "dsb_nhej", "mismatch_repair",
                 "photoprotective_pigment", "photoreactivation"]
METRICS = ["occurrence_pct", "abund_whole_pct"]

# ---- env ---------------------------------------------------------------------
env = {}
for r in csv.DictReader(open(ENV), delimiter="\t"):
    env[r["sample"]] = r
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

# ---- per-sample trait values -------------------------------------------------
def load(res):
    d = collections.defaultdict(dict)      # (sample, metric) -> trait -> value
    for r in csv.DictReader(open(f"{res}/community_by_sample.tsv"), delimiter="\t"):
        for m in METRICS:
            d[(r["sample"], m)][r["trait"]] = float(r[m])
    return d

osmo, uv = load(OSMO), load(UV)
samples = sorted({s for s, m in osmo} & {s for s, m in uv})
print(f"[grad] samples with both trait sets: {len(samples)}")

traits = {}                                # (label, metric) -> {sample: value}
for m in METRICS:
    for t in OSMO_ROUTES:
        traits[(f"osmo:{t}", m)] = {s: osmo[(s, m)][t] for s in samples}
    traits[("osmo:breadth_0to4", m)] = {s: sum(osmo[(s, m)][t] for t in OSMO_ROUTES) / 100 for s in samples}
    for t in REPAIR_GROUPS:
        traits[(f"uv:{t}", m)] = {s: uv[(s, m)][t] for s in samples}
    traits[("uv:breadth_0to8", m)] = {s: sum(uv[(s, m)][t] for t in REPAIR_GROUPS) / 100 for s in samples}

# ---- partial Spearman (control ranked lat+lon) -------------------------------
def partial_spearman(x, y, lat, lon):
    rx, ry = rankdata(x), rankdata(y)
    C = np.column_stack([np.ones(len(x)), rankdata(lat), rankdata(lon)])
    bx, *_ = np.linalg.lstsq(C, rx, rcond=None)
    by, *_ = np.linalg.lstsq(C, ry, rcond=None)
    rho = pearsonr(rx - C @ bx, ry - C @ by)[0]
    df = len(x) - np.linalg.matrix_rank(C) - 1
    p = 2 * tdist.sf(abs(rho) * np.sqrt(df / max(1 - rho*rho, np.finfo(float).tiny)), df)
    return rho, p

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
        for (label, metric), tv in traits.items():
            x = np.array([tv[s] for s in sub])
            if np.all(x == x[0]):
                continue
            rho_m, p_m = spearmanr(x, y)
            rho_p, p_p = partial_spearman(x, y, lat, lon)
            stash.append(dict(block=block, axis=ax, n=len(sub), trait=label, metric=metric,
                              rho_marginal=rho_m, p_marginal=p_m,
                              rho_partial=rho_p, p_partial=p_p))
    for r, q in zip(stash, bh([r["p_partial"] for r in stash])):
        r["q_partial"] = q
    rows += stash

cols = ["block", "axis", "n", "trait", "metric",
        "rho_marginal", "p_marginal", "rho_partial", "p_partial", "q_partial"]
with open('/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p4_review_20261013/results/round4_corrections/gradient_tests.tsv', "w") as o:
    o.write("\t".join(cols) + "\n")
    for r in rows:
        o.write("\t".join(f"{r[c]:.4g}" if isinstance(r[c], float) else str(r[c]) for c in cols) + "\n")

sig = [r for r in rows if r["q_partial"] < 0.05]
print(f"\n[grad] {len(rows)} tests; FDR q<0.05 (partial, lat/lon-controlled): {len(sig)}")
print(f"{'axis':16s} {'trait':32s} {'metric':16s} {'n':>4s} {'rho_p':>7s} {'q':>9s}")
for r in sorted(sig, key=lambda r: r["q_partial"]):
    print(f"{r['axis']:16s} {r['trait']:32s} {r['metric']:16s} {r['n']:4d} "
          f"{r['rho_partial']:+7.3f} {r['q_partial']:9.2e}")
print(f"\n[grad] wrote {OSMO}/gradient_tests.tsv")
