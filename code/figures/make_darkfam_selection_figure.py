#!/usr/bin/env python3
"""Section 1.4 figure — coding constraint in the dark families, against a matched control.

One panel, by design: the Results text already states the medians and the percentages, so the
figure earns its place only by showing what a sentence cannot, namely that the two whole
distributions coincide across their range. Dark families are drawn above the axis and the length-
and size-matched annotated control below it, on one shared scale (the "mirror" grammar), so the
fairness of the comparison is visible rather than asserted. Occupancy, persistence and phylum
coherence stay in prose.

Numbers are asserted against the registered Section 1.4 values before drawing.
Reads M/06_darkmatter/darkfam_basis_clean/results/selection/dnds.tsv.
"""
import os, csv, statistics
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "..", "..", "figures")
WB = "/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/darkfam_basis_clean/results"
from palette import FOCAL, CONTROL, GREY, INK, MUTED
T1, T2, T3 = 8.0, 7.5, 6.8
plt.rcParams.update({"font.size": T1, "axes.spines.top": False, "axes.spines.right": False})

d = {"dark": [], "ctrl": []}
for r in csv.DictReader(open(os.path.join(WB, "selection", "dnds.tsv")), delimiter="\t"):
    d[r["set"]].append(float(r["median_dNdS"]))
med_d, med_c = statistics.median(d["dark"]), statistics.median(d["ctrl"])
lt1_d = 100 * sum(1 for x in d["dark"] if x < 1) / len(d["dark"])
lt1_c = 100 * sum(1 for x in d["ctrl"] if x < 1) / len(d["ctrl"])
assert (len(d["dark"]), len(d["ctrl"])) == (97200, 95647)
assert round(med_d, 3) == 0.425 and round(med_c, 3) == 0.388
assert round(lt1_d, 1) == 96.7 and round(lt1_c, 1) == 99.6
print(f"asserted: medians {med_d:.3f}/{med_c:.3f}; below 1 {lt1_d:.1f}%/{lt1_c:.1f}%")

fig, ax = plt.subplots(figsize=(5.6, 3.3))
bins = np.arange(0, 1.62, 0.03)
ctr = (bins[:-1] + bins[1:]) / 2
hd = np.histogram(d["dark"], bins=bins, density=True)[0]
hc = np.histogram(d["ctrl"], bins=bins, density=True)[0]
ax.fill_between(ctr, hd, color=FOCAL, alpha=0.55, lw=0)
ax.fill_between(ctr, -hc, color=CONTROL, alpha=0.55, lw=0)
ax.plot(ctr, hd, color=FOCAL, lw=1.0); ax.plot(ctr, -hc, color=CONTROL, lw=1.0)
ax.axhline(0, color=INK, lw=0.8)
ax.axvline(1.0, ls="--", lw=0.9, color=GREY)
ax.text(1.03, 0.93, "neutral expectation:\nno protein product", transform=ax.get_xaxis_transform(),
        fontsize=T3, color=MUTED, va="top")
top, bot = hd.max(), -hc.max()
ax.plot([med_d, med_d], [0, top * 0.72], color=FOCAL, lw=1.1)
ax.plot([med_c, med_c], [0, bot * 0.72], color=CONTROL, lw=1.1)
ax.text(med_d + 0.02, top * 0.74, f"{med_d:.3f}", fontsize=T2, color=FOCAL, fontweight="bold")
ax.text(med_c - 0.02, bot * 0.74, f"{med_c:.3f}", fontsize=T2, color=CONTROL, fontweight="bold",
        ha="right", va="top")
ax.text(0.02, 0.97, f"dark families\nn = {len(d['dark']):,}   {lt1_d:.1f}% below 1",
        transform=ax.transAxes, fontsize=T2, color=FOCAL, va="top", fontweight="bold")
ax.text(0.98, 0.05, f"matched annotated control\nn = {len(d['ctrl']):,}   {lt1_c:.1f}% below 1",
        transform=ax.transAxes, fontsize=T2, color=CONTROL, va="bottom", ha="right",
        fontweight="bold")
ax.text(0.98, 0.72, "3.3% of dark and 0.4% of control\nfamilies lie above 1",
        transform=ax.transAxes, fontsize=T3, color=MUTED, ha="right", va="top")
ax.set_yticks([]); ax.set_xlim(0, 1.6)
ax.set_xlabel("per-family median dN/dS"); ax.set_ylabel("share of families")
ax.spines["left"].set_visible(False)
fig.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(os.path.join(FIG, f"fig_darkfam_selection.{ext}"), bbox_inches="tight", dpi=200)
print("wrote figures/fig_darkfam_selection.{pdf,png}")
