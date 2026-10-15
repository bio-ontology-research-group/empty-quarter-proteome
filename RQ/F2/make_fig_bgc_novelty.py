#!/usr/bin/env python3
"""RQ-F2 historical figure: % of EQ GCFs with no MIBiG member.

This is an operational, fragmentation-sensitive upper bound, not a biological novelty estimate. Run after pulling
the result tables
from IBEX into RQ/F2/results/.

  python3 RQ/F2/make_fig_bgc_novelty.py

Inputs : RQ/F2/results/result_summary.json
         RQ/F2/results/novelty_by_chemclass.tsv
Outputs: RQ/F2/fig_bgc_novelty.{pdf,png}
"""
import os
import json
import csv

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "results")

with open(os.path.join(RES, "result_summary.json")) as fh:
    S = json.load(fh)

rows = []
with open(os.path.join(RES, "novelty_by_chemclass.tsv")) as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        r["n_eq_bgc"] = int(r["n_eq_bgc"])
        r["pct_gcf_novel"] = float(r["pct_gcf_novel"])
        r["n_gcf"] = int(r["n_gcf"])
        rows.append(r)
rows = [r for r in rows if r["n_eq_bgc"] > 0][:10]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.6),
                               gridspec_kw={"width_ratios": [1, 1.5]})

# Panel A: overall no-MIBiG-member vs shared GCFs
novel = S["n_gcf_novel"]
shared = S["n_gcf_shared_with_mibig"]
ax1.bar(["no MIBiG\nmember", "shared\nwith MIBiG"], [novel, shared],
        color=["#c0392b", "#7f8c8d"])
ax1.set_ylabel("Empty Quarter GCFs")
ax1.set_title("A. {:.1f}% of {:,} EQ GCFs lack a MIBiG member".format(
    S["pct_gcf_novel"], S["n_gcf_eq"]), fontsize=11)
for i, v in enumerate([novel, shared]):
    ax1.text(i, v, "{:,}".format(v), ha="center", va="bottom", fontsize=9)

# Panel B: no-MIBiG-member GCFs by biosynthetic class
labels = [r["class"] for r in rows]
vals = [r["pct_gcf_novel"] for r in rows]
ncl = [r["n_gcf"] for r in rows]
y = range(len(labels))
ax2.barh(list(y), vals, color="#c0392b")
ax2.set_yticks(list(y))
ax2.set_yticklabels(labels, fontsize=9)
ax2.invert_yaxis()
ax2.set_xlabel("% of GCFs with no MIBiG member")
ax2.set_xlim(0, 100)
ax2.set_title("B. Named-reference no-hit fraction by class", fontsize=11)
for i, (v, n) in enumerate(zip(vals, ncl)):
    ax2.text(min(v + 1, 99), i, "{:.1f}%  (n={:,} GCFs)".format(v, n),
             va="center", fontsize=8)

fig.suptitle(
    "Historical EQ/MIBiG 3.1 BiG-SLICE clustering (T={}); 90.1% contig-edge regions"
    .format(S["threshold"]), fontsize=12)
fig.text(0.5, 0.01, "Fragmentation-sensitive upper bound; not active manuscript evidence", ha="center", fontsize=8)
fig.tight_layout(rect=[0, 0.04, 1, 0.95])
for ext in ("pdf", "png"):
    kwargs = {"dpi": 200}
    if ext == "pdf":
        kwargs["metadata"] = {"CreationDate": None, "ModDate": None}
    fig.savefig(os.path.join(HERE, "fig_bgc_novelty." + ext), **kwargs)
print("wrote fig_bgc_novelty.{pdf,png}")
