#!/usr/bin/env python3
"""Historical Phase D MGnify90-relative no-hit analysis.

Joins the archived DIAMOND search of 23,382 UniRef50-dark representatives against MGnify90. A hit or no-hit is
only a result under that named database/search rule. Pairwise NG86 and occupancy comparisons are unclustered
descriptive diagnostics; they do not establish geographic absence, endemism, or selection.
"""
import os, csv
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ECO = os.path.join(HERE, "..", "results", "darkmatter_ecology")
HITS = os.path.join(ECO, "endemism", "mgnify_hits.m8")
META = os.path.join(ECO, "family_ec.tsv")          # has all 23,382 family ids in col1

allfam = set()
for i, l in enumerate(open(META)):
    if i == 0: continue
    allfam.add(l.split("\t")[0])

best_pid = {}
if os.path.exists(HITS):
    for line in open(HITS):
        c = line.rstrip("\n").split("\t")
        q, pid = c[0], float(c[2])
        if q not in best_pid or pid > best_pid[q]: best_pid[q] = pid
shared = set(best_pid)                              # >=1 MGnify hit
endemic = allfam - shared

# selection
dnds = {}
dn = os.path.join(ECO, "dnds.tsv")
if os.path.exists(dn):
    for r in csv.DictReader(open(dn), delimiter="\t"):
        if r["set"] == "dark" and r["median_dNdS"]: dnds[r["family_id"]] = float(r["median_dNdS"])
# within-EQ occupancy from labels
prevn = {}
lab = os.path.join(ECO, "selection", "phase_b_family_labels.tsv")
if not os.path.exists(lab): lab = os.path.join(ECO, "phase_b_family_labels.tsv")
if os.path.exists(lab):
    for r in csv.DictReader(open(lab), delimiter="\t"): prevn[r["family_id"]] = int(r["n_present"])

import statistics as st
out = []
def w(x): out.append(x); print(x)
w("# Historical Phase D comparison against MGnify90")
w("\n**Status:** historical descriptive screen; not promoted as an endemism or selection result.\n")
w(f"- Sequence-dark families searched: {len(allfam):,} (all were already no-hits in the frozen UniRef50 screen).")
w(f"- Families with a detected MGnify90 homolog under the historical rule: {len(shared):,} ({100*len(shared)/len(allfam):.1f}%).")
w(f"- Families with no detected MGnify90 homolog under that rule: {len(endemic):,} ({100*len(endemic)/len(allfam):.1f}%).")
if best_pid:
    pv = sorted(best_pid.values())
    w(f"- Median best-hit amino-acid identity among detected homologs: {st.median(pv):.1f}%.")
w("\nThe no-hit class is **MGnify90-relative**. It does not establish geographic absence, desert specificity, or endemism.")

if dnds:
    from scipy.stats import mannwhitneyu
    ds = [dnds[f] for f in shared if f in dnds]; de = [dnds[f] for f in endemic if f in dnds]
    if ds and de:
        U, p = mannwhitneyu(de, ds, alternative="two-sided")
        w(f"\nHistorical unclustered diagnostic: median pairwise NG86 dN/dS {st.median(de):.3f} in the no-hit class "
          f"versus {st.median(ds):.3f} in the hit class (Mann--Whitney p={p:.2e}). This is descriptive and does not "
          "establish a difference in evolutionary constraint or selection.")
if prevn:
    pe = [prevn[f] for f in endemic if f in prevn]; ps = [prevn[f] for f in shared if f in prevn]
    if pe and ps:
        w(f"Historical within-EQ occupancy diagnostic: median samples/sites present {st.median(pe):.0f} in the "
          f"MGnify90 no-hit class versus {st.median(ps):.0f} in the hit class.")

open(os.path.join(ECO, "GLOBAL_ENDEMISM.md"), "w").write("\n".join(out) + "\n")
# figure: best-hit identity distribution + named-database hit/no-hit split
if best_pid:
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    ax[0].hist(list(best_pid.values()), bins=30, color="#4477aa")
    ax[0].set_xlabel("best MGnify hit %identity"); ax[0].set_ylabel("families")
    ax[0].set_title("Identity of detected MGnify90 homologs")
    ax[1].bar(["MGnify90\nhit", "MGnify90\nno hit"], [len(shared), len(endemic)], color=["#4477aa", "#cc6677"])
    ax[1].set_ylabel("families"); ax[1].set_title("Named-database detection result")
    fig.text(0.5, 0.01, "Historical MGnify90-relative screen; no-hit is not geographic absence", ha="center", fontsize=8)
    fig.tight_layout(rect=[0, 0.04, 1, 1])
    for ext in ("pdf", "png"):
        kwargs = {"bbox_inches": "tight", "dpi": 200}
        if ext == "pdf": kwargs["metadata"] = {"CreationDate": None, "ModDate": None}
        fig.savefig(os.path.join(HERE, "..", "figures", f"fig_darkmatter_global_endemism.{ext}"), **kwargs)
    print("wrote fig_darkmatter_global_endemism")
print("\nwritten -> results/darkmatter_ecology/GLOBAL_ENDEMISM.md")
