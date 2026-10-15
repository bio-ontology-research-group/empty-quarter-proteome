#!/usr/bin/env python3
"""Dark-matter ecology paper figures (matplotlib). Writes pdf+png to ../figures/. Separate from
make_figures.py (other pillars). Run: python3 code/make_darkmatter_figures.py
  1. fig_darkmatter_pcoa        — PCoA of dark-proteome composition, by aridity + compartment
  2. fig_darkmatter_permanova   — PERMANOVA R2 per env block (marginal vs |space)
  3. fig_darkmatter_dnds        — bounded dark-family dN/dS distribution and size sensitivity
  4. fig_darkmatter_endemism    — descriptive pairwise NG86 vs family prevalence
"""
import os, csv
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import eco_stats as es
from make_darkmatter_dnds_figure import build_figure as build_dnds_figure

HERE = os.path.dirname(os.path.abspath(__file__))
ECO = os.path.join(HERE, "..", "results", "darkmatter_ecology")
FIG = os.path.join(HERE, "..", "figures")
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

# ---------- load community matrix ----------
with open(os.path.join(ECO, "darkfam_abundance.tsv")) as fh:
    samples = fh.readline().rstrip("\n").split("\t")[1:]
    data = [[float(x) if x else 0.0 for x in line.rstrip("\n").split("\t")[1:]] for line in fh]
M = np.array(data, float).T
M = M[:, (M > 0).sum(0) >= 5]

# ---------- Fig 1: PCoA ----------
keep = [i for i, s in enumerate(samples) if s in env]
Mk = M[keep]; sk = [samples[i] for i in keep]
coords, expl = es.pcoa(es.bray_curtis(es.hellinger(Mk)), k=2)
arid = np.array([fnum(env[s].get("aridity_index")) for s in sk], float)
comp = [env[s].get("compartment", "") for s in sk]
fig, ax = plt.subplots(1, 2, figsize=(11, 4.6))
m = ~np.isnan(arid)
sc = ax[0].scatter(coords[m, 0], coords[m, 1], c=arid[m], cmap="viridis_r", s=14, alpha=0.8)
ax[0].set_title("Dark-proteome composition by aridity"); fig.colorbar(sc, ax=ax[0], label="aridity index (low=arid)")
for cmp, cc in {"Surface": "#d95f02", "Deep": "#7570b3", "PlantRoot": "#1b9e77"}.items():
    idx = [i for i, c in enumerate(comp) if c == cmp]
    if idx: ax[1].scatter(coords[idx, 0], coords[idx, 1], c=cc, s=14, alpha=0.8, label=cmp)
ax[1].legend(frameon=False, fontsize=8); ax[1].set_title("by compartment")
for a in ax:
    a.set_xlabel(f"PCoA1 ({expl[0]*100:.1f}%)"); a.set_ylabel(f"PCoA2 ({expl[1]*100:.1f}%)")
save(fig, "fig_darkmatter_pcoa")

# ---------- Fig 2: PERMANOVA R2 ----------
blocks, r2, kind = [], [], []
sm = os.path.join(ECO, "phase_b_summary.txt")
if os.path.exists(sm):
    want = ("compartment", "space(latlon)", "climate", "climate|space",
            "geochem", "geochem|space", "vanadium")
    for ln in open(sm):
        p = ln.split()
        if len(p) >= 5 and p[0] in want and p[0] not in blocks:
            try:
                v = float(p[2]) * 100        # only commit after the numeric parse succeeds
            except ValueError:
                continue
            blocks.append(p[0]); r2.append(v)
            kind.append("|space" if "|space" in p[0] else "marginal")
if blocks:
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.bar(range(len(blocks)), r2, color=["#4477aa" if k == "marginal" else "#aa3377" for k in kind])
    ax.set_xticks(range(len(blocks))); ax.set_xticklabels(blocks, rotation=40, ha="right", fontsize=8)
    ax.set_ylabel("variance explained R² (%)")
    ax.set_title("Environment structures the dark proteome\n(PERMANOVA; all p=0.001 except vanadium)")
    save(fig, "fig_darkmatter_permanova")

# ---------- Fig 3 & 4: dN/dS ----------
dn = os.path.join(ECO, "dnds.tsv")
if os.path.exists(dn):
    rows = list(csv.DictReader(open(dn), delimiter="\t"))
    build_dnds_figure(
        Path(dn),
        Path(ECO) / "darkfam_meta.tsv",
        Path(FIG),
    )

    lab = os.path.join(ECO, "selection", "phase_b_family_labels.tsv")
    if not os.path.exists(lab): lab = os.path.join(ECO, "phase_b_family_labels.tsv")
    if os.path.exists(lab):
        prev_l = {r["family_id"]: int(r["n_present"]) for r in csv.DictReader(open(lab), delimiter="\t")}
        dd = {r["family_id"]: float(r["median_dNdS"]) for r in rows if r["set"] == "dark" and r["median_dNdS"]}
        keys = [f for f in dd if f in prev_l]
        xs = np.array([prev_l[f] for f in keys]); ys = np.array([dd[f] for f in keys])
        order = np.argsort(xs); xs, ys = xs[order], ys[order]
        nb = 12; edges = np.linspace(0, len(xs), nb + 1).astype(int)
        bx = [xs[edges[i]:edges[i+1]].mean() for i in range(nb) if edges[i+1] > edges[i]]
        by = [np.median(ys[edges[i]:edges[i+1]]) for i in range(nb) if edges[i+1] > edges[i]]
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.scatter(xs, ys, s=3, alpha=0.07, color="#332288")
        ax.plot(bx, by, "o-", color="#cc6677", label="binned median")
        ax.set_xlabel("family prevalence (n samples present)"); ax.set_ylabel("median dN/dS")
        ax.set_title("Prevalence versus pairwise NG86 estimate\n(descriptive; unclustered)"); ax.legend(frameon=False)
        save(fig, "fig_darkmatter_endemism")
print("figures done")
