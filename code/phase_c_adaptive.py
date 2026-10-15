#!/usr/bin/env python3
"""Historical Phase C environment-association versus pairwise-NG86 diagnostic.

Per dark family (those with dN/dS), test abundance association with aridity (Spearman), salinity
(Spearman, XRF subset) and compartment (Kruskal-Wallis); BH-FDR per test; structured = q<0.05 in >=1
axis. It then records unclustered Mann--Whitney and prevalence correlations. These dependent pairwise summaries
do not establish evolutionary constraint, adaptation, or selection. Writes phase_b_family_labels.tsv and the
historically named SELECTION_ADAPTIVE.md.
"""
import os, csv, math
import numpy as np
from scipy.stats import spearmanr, kruskal, mannwhitneyu

HERE = os.path.dirname(os.path.abspath(__file__))
ECO = os.path.join(HERE, "..", "results", "darkmatter_ecology")
ENVF = os.path.join(ECO, "env_matrix.tsv")
ABUND = os.path.join(ECO, "darkfam_abundance.tsv")
DNDS = os.path.join(ECO, "dnds.tsv")

def fnum(x):
    try: return float(x)
    except (TypeError, ValueError): return None

env = {r["sample"]: r for r in csv.DictReader(open(ENVF), delimiter="\t")}
dnds = {r["family_id"]: float(r["median_dNdS"]) for r in csv.DictReader(open(DNDS), delimiter="\t")
        if r["set"] == "dark" and r["median_dNdS"]}

# load abundance header + iterate rows (only families with dN/dS)
with open(ABUND) as fh:
    samples = fh.readline().rstrip("\n").split("\t")[1:]
arid = np.array([fnum(env.get(s, {}).get("aridity_index")) if s in env else None for s in samples], dtype=object)
sal  = np.array([fnum(env.get(s, {}).get("salinity_proxy")) if s in env else None for s in samples], dtype=object)
comp = np.array([env.get(s, {}).get("compartment", "") if s in env else "" for s in samples])
def mask_num(arr):
    m = np.array([v is not None for v in arr]); return m, np.array([v if v is not None else np.nan for v in arr], float)
m_ar, v_ar = mask_num(arid); m_sal, v_sal = mask_num(sal)
comps = sorted(set(c for c in comp if c))

labels = {}   # fid -> dict
def bh(pdict):
    items = [(f, p) for f, p in pdict.items() if p is not None and not math.isnan(p)]
    items.sort(key=lambda x: x[1]); n = len(items); q = {}
    prev = 1.0
    for i in range(n-1, -1, -1):
        f, p = items[i]; val = min(prev, p*n/(i+1)); q[f] = val; prev = val
    return q

p_ar, p_sal, p_comp, prev_n = {}, {}, {}, {}
with open(ABUND) as fh:
    next(fh)
    for line in fh:
        c = line.rstrip("\n").split("\t")
        fid = c[0]
        if fid not in dnds: continue
        tpm = np.array([float(x) if x else 0.0 for x in c[1:]], float)
        prev_n[fid] = int((tpm > 0).sum())
        # aridity
        idx = m_ar & (~np.isnan(v_ar))
        if idx.sum() >= 20 and (tpm[idx] > 0).sum() >= 5 and np.ptp(tpm[idx]) > 0:
            r, p = spearmanr(tpm[idx], v_ar[idx]); p_ar[fid] = p
        # salinity
        idx = m_sal & (~np.isnan(v_sal))
        if idx.sum() >= 20 and (tpm[idx] > 0).sum() >= 5 and np.ptp(tpm[idx]) > 0:
            r, p = spearmanr(tpm[idx], v_sal[idx]); p_sal[fid] = p
        # compartment
        groups = [tpm[comp == g] for g in comps]
        groups = [g for g in groups if len(g) >= 5 and np.ptp(g) > 0]
        if len(groups) >= 2 and sum((g > 0).sum() for g in groups) >= 5:
            try:
                H, p = kruskal(*groups); p_comp[fid] = p
            except ValueError:
                pass

q_ar, q_sal, q_comp = bh(p_ar), bh(p_sal), bh(p_comp)
out = []
def w(x): out.append(x); print(x)
lab = open(os.path.join(ECO, "selection", "phase_b_family_labels.tsv"), "w") \
    if os.path.isdir(os.path.join(ECO, "selection")) else open(os.path.join(ECO, "phase_b_family_labels.tsv"), "w")
lab.write("family_id\tstructured\tn_present\tq_aridity\tq_salinity\tq_compartment\n")
struct_set = set()
for fid in dnds:
    qa, qs, qc = q_ar.get(fid), q_sal.get(fid), q_comp.get(fid)
    s = 1 if any(q is not None and q < 0.05 for q in (qa, qs, qc)) else 0
    if s: struct_set.add(fid)
    lab.write(f"{fid}\t{s}\t{prev_n.get(fid,0)}\t{qa if qa else 'NA'}\t{qs if qs else 'NA'}\t{qc if qc else 'NA'}\n")
lab.close()

vs = [dnds[f] for f in dnds if f in struct_set]
vc = [dnds[f] for f in dnds if f not in struct_set]
import statistics as st
w("# Historical Phase C environment-association versus pairwise-NG86 diagnostic")
w("\n**Status:** exploratory unclustered output; not promoted as an adaptation or selection test.\n")
w(f"dark families with dN/dS: {len(dnds)}")
w(f"environmentally STRUCTURED (FDR<0.05 in aridity/salinity/compartment): {len(vs)}")
w(f"not passing that association screen: {len(vc)}")
if vs and vc:
    U, p = mannwhitneyu(vs, vc, alternative="two-sided")
    w(f"\nHistorical medians: association-screen-positive {st.median(vs):.3f}; screen-negative {st.median(vc):.3f}.")
    w(f"Unclustered Mann--Whitney p={p:.2e}. Orthogroup/pair dependence is not modelled, so this does not support "
      "a difference in evolutionary constraint or adaptation.")
# descriptive occupancy correlation
allf = list(dnds); pv = np.array([prev_n.get(f, 0) for f in allf]); dv = np.array([dnds[f] for f in allf])
r, p = spearmanr(pv, dv)
w(f"\nHistorical occupancy diagnostic: Spearman(prevalence, pairwise NG86 dN/dS) r={r:+.3f}, nominal p={p:.2e}. "
  "This is descriptive and not an endemism or selection result.")
open(os.path.join(ECO, "SELECTION_ADAPTIVE.md"), "w").write("\n".join(out) + "\n")
print("\nwritten -> results/darkmatter_ecology/SELECTION_ADAPTIVE.md + phase_b_family_labels.tsv")
