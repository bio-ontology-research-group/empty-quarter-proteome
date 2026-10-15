#!/usr/bin/env python3
"""§2.5 figures on the rebuilt 106,042-family universe. Writes pdf+png to ../figures/.

Replaces figs 1-2 of make_darkmatter_figures.py (June, 23,382 object): reads the precomputed
ordination and PERMANOVA tables of the workbench rerun instead of recomputing from the matrix.
  1. fig_darkmatter_pcoa      — PCoA of dark-family composition, by aridity + compartment
  2. fig_darkmatter_permanova — PERMANOVA R2 per env block (marginal vs |space)
Inputs: M/06_darkmatter/darkfam_ecology_clean/results/{pcoa_coords,pcoa_explained,permanova_blocks}.tsv
        + results/darkmatter_ecology/env_matrix.tsv (aridity join).
Run: python3 code/make_darkfam_env_figures.py
"""
import os, csv
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ECO = os.path.join(HERE, "..", "..", "results", "darkmatter_ecology")
WB = "/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/darkfam_ecology_clean/results"
FIG = os.path.join(HERE, "..", "..", "figures")
os.makedirs(FIG, exist_ok=True)
plt.rcParams.update({"font.size": 11, "axes.spines.top": False, "axes.spines.right": False})

def save(fig, name):
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(FIG, f"{name}.{ext}"), bbox_inches="tight", dpi=200)
    plt.close(fig); print("wrote", name)

def fnum(x):
    try: return float(x)
    except (TypeError, ValueError): return None

env = {r["sample"]: r for r in csv.DictReader(open(os.path.join(ECO, "env_matrix.tsv")), delimiter="\t")}

# ---------- Fig 1: PCoA (precomputed coordinates) ----------
rows = list(csv.DictReader(open(os.path.join(WB, "pcoa_coords.tsv")), delimiter="\t"))
expl = {r["axis"]: float(r["explained"]) for r in
        csv.DictReader(open(os.path.join(WB, "pcoa_explained.tsv")), delimiter="\t")}
coords = np.array([[float(r["PC1"]), float(r["PC2"])] for r in rows])
arid = np.array([fnum(env.get(r["sample"], {}).get("aridity_index")) for r in rows], float)
comp = [r["compartment"] for r in rows]
fig, ax = plt.subplots(1, 2, figsize=(12, 4.6), gridspec_kw={"wspace": 0.38})
m = ~np.isnan(arid)
sc = ax[0].scatter(coords[m, 0], coords[m, 1], c=arid[m], cmap="viridis_r", s=14, alpha=0.8)
ax[0].set_title("Dark-family composition by aridity"); fig.colorbar(sc, ax=ax[0], label="aridity index (low=arid)")
from palette import COMPARTMENT, FOCAL, FOCAL_LIGHT
for cmp, cc in COMPARTMENT.items():
    idx = [i for i, c in enumerate(comp) if c == cmp]
    if idx: ax[1].scatter(coords[idx, 0], coords[idx, 1], c=cc, s=14, alpha=0.8, label=cmp)
ax[1].legend(frameon=False, fontsize=8); ax[1].set_title("by compartment")
for a in ax:
    a.set_xlabel(f"PCoA1 ({expl['PC1']*100:.1f}%)"); a.set_ylabel(f"PCoA2 ({expl['PC2']*100:.1f}%)")
save(fig, "fig_darkmatter_pcoa")

# ---------- Fig 2: PERMANOVA R2 ----------
want = ("compartment", "space(latlon)", "climate", "climate|space",
        "geochem", "geochem|space", "vanadium")
blocks, r2, kind = [], [], []
for r in csv.DictReader(open(os.path.join(WB, "permanova_blocks.tsv")), delimiter="\t"):
    if r["block"] in want and r["block"] not in blocks:
        blocks.append(r["block"]); r2.append(float(r["R2_marginal"]) * 100)
        kind.append("|space" if "|space" in r["block"] else "marginal")
fig, ax = plt.subplots(figsize=(7, 4))
ax.bar(range(len(blocks)), r2, color=[FOCAL if k == "marginal" else FOCAL_LIGHT for k in kind])
ax.set_xticks(range(len(blocks))); ax.set_xticklabels(blocks, rotation=40, ha="right", fontsize=8)
ax.set_ylabel("variance explained R² (%)")
ax.set_title("Environment structures the dark families\n(PERMANOVA; all p=0.001 except vanadium)")
save(fig, "fig_darkmatter_permanova")
print("figures done")
