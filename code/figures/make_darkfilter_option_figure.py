#!/usr/bin/env python3
"""Section 2.5 option panel — the gradient acts as a filter (per-axis association census).

DESIGN DECISION (2026-08-23, session 09-display). Candidate third panel for fig:environment: dark
families associated with each geochemical/climate axis at q<0.05, ordered; the % negative annotation on
each bar is the "filter, not recruitment" reading. The six all-zero trace-metal axes are drawn as
explicit zero rows because the nulls are load-bearing (vanadium legend). Bars #4477aa; zero rows grey text.

Reads  M/06_darkmatter/darkfam_ecology_clean/results/indicators_summary.tsv.
Writes figures/options/fig_darkfam_filter_opt.pdf. Asserted: S 54,227 / Ca 47,174 / Cl 46,972 /
salinity 37,280 / aridity 23,031; V/Mo/Ni/Co/Mn/Zn = 0.
"""
import os, csv
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.join(HERE, "..", "..")
OPT = os.path.join(REPO, "figures", "options")
WB = "/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/darkfam_ecology_clean/results"
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
from palette import C, INK, MUTED, GRID

PRETTY = {"S": "sulfur", "Cl": "chloride", "Ca": "calcium", "Na": "sodium", "K": "potassium",
          "Fe": "iron", "Cu": "copper", "salinity_proxy": "salinity proxy",
          "aridity_index": "aridity index", "V": "vanadium", "Mo": "molybdenum", "Ni": "nickel",
          "Co": "cobalt", "Mn": "manganese", "Zn": "zinc"}
rows = [r for r in csv.DictReader(open(os.path.join(WB, "indicators_summary.tsv")), delimiter="\t")
        if r["axis"] != "plantroot_vs_other"]
d = {r["axis"]: (int(r["n_q05"]), r["pct_negative_q05"], int(r["n_tested"])) for r in rows}
assert d["S"][0] == 54227 and d["Ca"][0] == 47174 and d["Cl"][0] == 46972 \
    and d["salinity_proxy"][0] == 37280 and d["aridity_index"][0] == 23031, d
assert all(d[k][0] == 0 for k in ("V", "Mo", "Ni", "Co", "Mn", "Zn"))
print("asserted: per-axis q<0.05 counts and the six trace-metal zeros")

order = sorted(d, key=lambda k: -d[k][0])
fig, a = plt.subplots(figsize=(7.2, 5.2))
ys = range(len(order) - 1, -1, -1)
for y, k in zip(ys, order):
    n, pctneg, ntest = d[k]
    if n:
        a.barh([y], [n], height=0.62, color=C, lw=0)
        a.text(n + 700, y, f"{n:,}  ({float(pctneg):.0f}% negative)", va="center", fontsize=8, color=INK)
    else:
        a.text(700, y, "0", va="center", fontsize=8, color=MUTED)
a.set_yticks(list(ys))
a.set_yticklabels([PRETTY.get(k, k) for k in order], fontsize=8.5)
a.set_xlabel("dark families associated at q < 0.05 (of 94,568–100,293 tested)")
a.set_xlim(0, 72000)
a.grid(axis="x", color=GRID, lw=0.6); a.set_axisbelow(True)
a.set_title("The gradient excludes: salt and aridity filter the dark repertoire", loc="left", fontsize=10)
fig.tight_layout()
os.makedirs(OPT, exist_ok=True)
fig.savefig(os.path.join(OPT, "fig_darkfam_filter_opt.pdf"), bbox_inches="tight")
print("wrote figures/options/fig_darkfam_filter_opt.pdf")
