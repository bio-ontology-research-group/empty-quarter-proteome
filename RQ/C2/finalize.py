#!/usr/bin/env python3
"""Finalize the historical exploratory RQ-C2 output.

This computes the archived alpha/accessory association and draws a visibly labelled provenance figure. It does
not promote the raw-union/dRep result to the current manuscript catalogue. Heavy clustering ran on IBEX (job
47843640, scripts at workbench 140_pangenome_openness.*); this only processes the small pulled summaries.
Run: python3 RQ/C2/finalize.py  (from repo root)
"""
import json, csv, os
import statistics as st
from scipy.stats import fisher_exact
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
d = json.load(open(os.path.join(HERE, "enrichment_summary.json")))
t = d["enrichment_2x2"]
OR, p = fisher_exact([[t["accessory_dark"], t["accessory_known"]], [t["core_dark"], t["core_known"]]])
alpha = [float(r["heaps_alpha"]) for r in csv.DictReader(open(os.path.join(HERE, "species_openness.tsv")), delimiter="\t")
         if r.get("heaps_alpha") not in (None, "", "NA")]
pct_open = 100 * sum(a < 1 for a in alpha) / len(alpha)
acc_dark = t["accessory_dark"] / (t["accessory_dark"] + t["accessory_known"])
core_dark = t["core_dark"] / (t["core_dark"] + t["core_known"])
print(f"species={d['n_species']} genomes={d['n_genomes_total']} open={pct_open:.1f}% medianAlpha={st.median(alpha):.3f} "
      f"core/pan={d['median_core_pan_ratio']:.3f} OR={OR:.2f} p={p:.1e}")

fig, (a, b) = plt.subplots(1, 2, figsize=(10, 4))
a.hist(alpha, bins=30, color="#4878a8", edgecolor="k", lw=0.3)
a.axvline(1.0, ls="--", color="red", label="α=1 (open↔closed)")
a.axvline(st.median(alpha), ls=":", color="k", label=f"median {st.median(alpha):.2f}")
a.set_xlabel("Heaps' law α (per species)"); a.set_ylabel("# species")
a.set_title(f"(a) {pct_open:.0f}% with fitted α<1", loc="left", weight="bold", fontsize=10)
a.legend(fontsize=8, frameon=False)
b.bar(["accessory", "core"], [100 * acc_dark, 100 * core_dark],
      color=["#c97b3c", "#4878a8"], edgecolor="k")
b.set_ylabel("% dark (UniRef50 no-hit) genes")
b.set_title(f"(b) historical no-hit association\nOR={OR:.2f}, nominal Fisher p<1e-300", loc="left", weight="bold", fontsize=10)
fig.text(0.5, 0.01, "Exploratory historical catalogue; not active manuscript evidence", ha="center", fontsize=8)
fig.tight_layout(rect=[0, 0.04, 1, 1])
for ext in ("pdf", "png"):
    kwargs = {"bbox_inches": "tight", "dpi": 200}
    if ext == "pdf":
        kwargs["metadata"] = {"CreationDate": None, "ModDate": None}
    fig.savefig(os.path.join(HERE, f"fig_pangenome_openness.{ext}"), **kwargs)
print("wrote fig_pangenome_openness.{pdf,png}")
