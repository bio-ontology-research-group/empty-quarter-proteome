#!/usr/bin/env python3
"""Best-hit ANI per cultured RAK-PR genome against the 5,438 soil MAG species representatives.

Form: stacked histogram (magnitude by identity), one stack per batch. Genomes with no hit at the
75 % screen cannot be placed on an ANI axis, so they get their own detached bar rather than being
dropped or parked at 75.
Colour: categorical slots 1-3 of the reference palette (blue / orange / aqua), the documented
three-slot set that clears the all-pairs CVD and normal-vision floors in light mode.
"""
import csv, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

P = "/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/cultured_comparison_clean"
BATCHES = ["enrichment", "isolate", "soil"]
LABELS = {"enrichment": "enrichment", "isolate": "isolate", "soil": "thermophilic (soil batch)"}
COLOR = {"enrichment": "#2a78d6", "isolate": "#eb6834", "soil": "#1baf7a"}
INK, INK2, GRID = "#1f2328", "#57606a", "#d8dee4"

rows = list(csv.DictReader(open(os.path.join(P, "results/cultured_best_hit.tsv")), delimiter="\t"))
hits = {b: [float(r["ani"]) for r in rows if r["batch"] == b and r["ani"]] for b in BATCHES}
nohit = {b: sum(1 for r in rows if r["batch"] == b and not r["ani"]) for b in BATCHES}
n_shared = sum(1 for r in rows if r["bin"] == "shared")
n_related = sum(1 for r in rows if r["bin"] == "related")
n_none = sum(1 for r in rows if r["bin"] == "none")

edges = np.arange(75, 100.5, 1.0)
fig, (axL, axR) = plt.subplots(
    1, 2, figsize=(9.6, 4.5), gridspec_kw={"width_ratios": [1, 9], "wspace": 0.06})

# detached "no hit" stack
bottom = 0.0
for b in BATCHES:
    axL.bar(0, nohit[b], bottom=bottom, width=0.62, color=COLOR[b], edgecolor="white", linewidth=2)
    bottom += nohit[b]
axL.set_xlim(-0.6, 0.6); axL.set_xticks([0]); axL.set_xticklabels(["no hit"], fontsize=9, color=INK2)
axL.set_ylabel("cultured genomes", fontsize=10, color=INK2)

# ANI histogram
counts = [np.histogram(hits[b], bins=edges)[0] for b in BATCHES]
bottom = np.zeros(len(edges) - 1)
for b, c in zip(BATCHES, counts):
    axR.bar(edges[:-1], c, bottom=bottom, width=0.92, align="edge",
            color=COLOR[b], edgecolor="white", linewidth=2, label=LABELS[b])
    bottom += c
ymax = max(bottom.max(), sum(nohit.values())) * 1.16
for ax in (axL, axR):
    ax.set_ylim(0, ymax)
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_color(GRID)
    ax.tick_params(colors=INK2, labelsize=9, length=3)
    ax.set_axisbelow(True)
    ax.yaxis.grid(True, color=GRID, linewidth=0.6)
axR.set_yticklabels([])
axR.set_xlim(74.6, 100.4)
axR.set_xticks(np.arange(75, 101, 5))
axR.set_xlabel("best-hit ANI to a soil MAG species representative (%)", fontsize=10, color=INK2)

for x in (80, 95):
    axR.axvline(x, color=INK2, linewidth=1.0, linestyle=(0, (4, 3)), zorder=3)
axR.text(87.5, ymax * 0.955, f"related  {n_related}", ha="center", fontsize=9.5, color=INK)
axR.text(97.6, ymax * 0.955, f"shared species  {n_shared}", ha="center", fontsize=9.5, color=INK)
axL.text(0, ymax * 0.955, f"none  {n_none}", ha="center", fontsize=9.5, color=INK)

axR.legend(frameon=False, fontsize=9.5, labelcolor=INK, loc="upper left",
           bbox_to_anchor=(0.30, 0.86), handlelength=1.1, handleheight=1.1)
fig.suptitle(f"{n_shared} of {len(rows)} cultured RAK-PR genomes share a species with a soil MAG",
             x=0.055, y=0.965, ha="left", fontsize=12.5, color=INK, weight="bold")
fig.text(0.055, 0.905,
         "skani 0.3.1, best hit per cultured genome against the 5,438 dereplicated soil MAG",
         ha="left", fontsize=8.8, color=INK2)
fig.text(0.055, 0.855,
         "species representatives. Shared species: ANI >= 95 and query aligned fraction >= 30 %.",
         ha="left", fontsize=8.8, color=INK2)
fig.subplots_adjust(top=0.79, bottom=0.13, left=0.075, right=0.985)
out = os.path.join(P, "figures/ani_histogram.png")
fig.savefig(out, dpi=200, facecolor="white")
print("wrote", out)
print(f"shared {n_shared}  related {n_related}  none {n_none} (of which no hit at 75 %: "
      f"{sum(nohit.values())})")
