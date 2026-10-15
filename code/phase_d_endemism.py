#!/usr/bin/env python3
"""Quantify within-EQ occupancy of the frozen sequence-dark orthogroups.

The sequence-dark label is reference-relative to the frozen UniRef50 screen.
Current global sequence/structure reference refreshes are separate analyses.
This script measures site/trip occupancy and reports the historical,
unclustered occupancy-versus-selection comparison descriptively.
"""
import os, csv
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ECO = os.path.join(HERE, "..", "results", "darkmatter_ecology")
env = {r["sample"]: r for r in csv.DictReader(open(os.path.join(ECO, "env_matrix.tsv")), delimiter="\t")}
# sample -> site, trip from env_matrix (has site, trip cols)
site = {s: env[s].get("site", "") for s in env}
trip = {s: env[s].get("trip", "") for s in env}

with open(os.path.join(ECO, "darkfam_abundance.tsv")) as fh:
    samples = fh.readline().rstrip("\n").split("\t")[1:]
sidx = {s: i for i, s in enumerate(samples)}

# pairwise NG86 summaries for a descriptive occupancy diagnostic
dnds = {}
dn = os.path.join(ECO, "dnds.tsv")
if os.path.exists(dn):
    for r in csv.DictReader(open(dn), delimiter="\t"):
        if r["set"] == "dark" and r["median_dNdS"]: dnds[r["family_id"]] = float(r["median_dNdS"])

n_samp, n_site, n_trip = {}, {}, {}
with open(os.path.join(ECO, "darkfam_abundance.tsv")) as fh:
    next(fh)
    for line in fh:
        c = line.rstrip("\n").split("\t")
        fid = c[0]; vals = c[1:]
        present = [samples[i] for i, v in enumerate(vals) if v and float(v) > 0]
        n_samp[fid] = len(present)
        n_site[fid] = len(set(site.get(s, "") for s in present if site.get(s)))
        n_trip[fid] = len(set(trip.get(s, "") for s in present if trip.get(s)))

fams = list(n_samp)
ns = np.array([n_site[f] for f in fams]); nt = np.array([n_trip[f] for f in fams])
nsamp = np.array([n_samp[f] for f in fams])
out = []
def w(x): out.append(x); print(x)
w("# Phase D — within-EQ occupancy of the dark proteome")
w(f"dark families: {len(fams)} (sequence-dark relative to the frozen UniRef50 screen)")
w("\n## site occupancy (60 sites total)")
for thr, lab in [(1, "single-site"), (3, "<=3 sites"), (10, "<=10 sites"), (30, ">=30 sites (EQ-wide)")]:
    if lab.startswith(">="):
        n = int((ns >= 30).sum())
    else:
        n = int((ns <= thr).sum())
    w(f"  {lab:28s}: {n:6d} ({100*n/len(fams):.1f}%)")
w(f"  median sites/family = {int(np.median(ns))} ; median trips/family = {int(np.median(nt))} (of 5)")
w(f"\n  single-trip families (temporally restricted): {int((nt==1).sum())} ({100*(nt==1).mean():.1f}%)")

# occupancy versus pairwise-NG86 diagnostic; no clustered inference
if dnds:
    from scipy.stats import spearmanr
    keys = [f for f in fams if f in dnds]
    r, p = spearmanr([n_site[f] for f in keys], [dnds[f] for f in keys])
    w(f"\n## occupancy vs pairwise NG86 (historical unclustered diagnostic): Spearman(n_sites, dN/dS) r={r:+.3f} p={p:.2e}")
    w("  Report effect sizes descriptively; orthogroup dependence is not modelled.")
    site_restricted = [dnds[f] for f in keys if n_site[f] <= 3]
    wide = [dnds[f] for f in keys if n_site[f] >= 30]
    if site_restricted and wide:
        import statistics as st
        from scipy.stats import mannwhitneyu
        U, pp = mannwhitneyu(site_restricted, wide, alternative="two-sided")
        w(f"  <=3-site median dN/dS={st.median(site_restricted):.3f} vs >=30-site {st.median(wide):.3f}  unclustered MWU p={pp:.2e}")

open(os.path.join(ECO, "ENDEMISM.md"), "w").write("\n".join(out) + "\n")
# occupancy distribution figure
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(6, 4))
ax.hist(ns, bins=range(0, 62, 2), color="#117733", alpha=0.85)
ax.set_xlabel("number of sites a dark family occupies (of 60)"); ax.set_ylabel("number of families")
ax.set_title("Site occupancy of dark orthogroups (60-site survey)")
single_site_percent = 100 * float((ns == 1).sum()) / len(ns)
eq_wide_percent = 100 * float((ns >= 30).sum()) / len(ns)
ax.text(
    0.98,
    0.95,
    f"{single_site_percent:.1f}% single-site; {eq_wide_percent:.1f}% in >=30 sites",
    transform=ax.transAxes,
    ha="right",
    va="top",
    fontsize=9,
)
for ext in ("pdf", "png"):
    fig.savefig(os.path.join(HERE, "..", "figures", f"fig_darkmatter_endemism_sites.{ext}"), bbox_inches="tight", dpi=200)
print("\nwritten -> results/darkmatter_ecology/ENDEMISM.md + figures/fig_darkmatter_endemism_sites")
