#!/usr/bin/env python3
"""Two figures for the differential-function scan.

figures/passing_heatmap.png  passing features x contrast, signed effect. Two panels
    because the two effect units cannot share a scale: percentage points of carriage
    (compartment, trip) on the left, partial Spearman rho (chemistry) on the right.
    Diverging blue/red with a neutral gray midpoint; no rainbow, no dual axis.

figures/top6_paired.png      the six strongest paired results as slopegraphs, one line
    per site (Surface vs Deep or PlantRoot vs bulk) or per site x compartment unit
    (trip pairs), so the reader sees the pairing rather than a summary bar.
"""
import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import RES, FIG, load_design

MAX_KO_ROWS = 30

# dataviz reference palette: diverging = blue <-> red with a gray neutral midpoint
BLUE, RED, GRAY = "#2a78d6", "#e34948", "#f0efec"
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#b9b8b3"
DIVERGING = LinearSegmentedColormap.from_list("bl_gr_rd", [BLUE, GRAY, RED])
plt.rcParams.update({
    "font.size": 8, "axes.edgecolor": MUTED, "axes.labelcolor": INK2,
    "xtick.color": INK2, "ytick.color": INK2, "text.color": INK,
    "figure.facecolor": "#fcfcfb", "axes.facecolor": "#fcfcfb",
    "axes.spines.top": False, "axes.spines.right": False,
})


def short(name, feat, n=44):
    lab = (name or feat).split(" [EC:")[0]
    return feat + "  " + (lab[:n] + "..." if len(lab) > n else lab)


def heatmap():
    keep = pd.read_csv(os.path.join(RES, "passing_modules.tsv"), sep="\t")
    if keep.empty:
        print("nothing passed; no heatmap")
        return
    allt = pd.concat([pd.read_csv(os.path.join(RES, f"{a}_tests.tsv"), sep="\t")
                      for a in ("compartment", "trip", "chemistry")], ignore_index=True)
    ids = keep[["basis", "feature"]].drop_duplicates()
    sub = allt.merge(ids, on=["basis", "feature"])
    sub["row"] = [short(n, f) for n, f in zip(sub["name"].fillna(""), sub["feature"])]

    # every passing module, but only the strongest KOs, or the figure becomes a poster
    strength = (keep.assign(a=keep["effect"].abs())
                    .groupby(["basis", "feature"])["a"].max()
                    .sort_values(ascending=False))
    mods = [k for k in strength.index if k[0] == "module"]
    kos = [k for k in strength.index if k[0] == "ko"][:MAX_KO_ROWS]
    order = [k for k in strength.index if k in set(mods) | set(kos)]
    n_ko_all = len([k for k in strength.index if k[0] == "ko"])
    rowlab = {(b, f): short(
        sub[(sub.basis == b) & (sub.feature == f)]["name"].fillna("").iloc[0], f)
        for b, f in order}
    rows = [rowlab[k] for k in order]

    panels = [("percentage points of carriage", ["compartment", "trip"]),
              ("partial Spearman rho", ["chemistry"])]
    mats = []
    for _, axes in panels:
        s = sub[sub["axis"].isin(axes)]
        m = s.pivot_table(index="row", columns="contrast", values="effect", aggfunc="first")
        mats.append(m.reindex(rows))

    widths = [max(1, m.shape[1]) for m in mats]
    fig, axs = plt.subplots(1, 2, figsize=(1.0 + 0.42 * sum(widths) + 5.2,
                                           1.4 + 0.20 * len(rows)),
                            gridspec_kw={"width_ratios": widths}, sharey=True)
    for ax, m, (unit, _) in zip(axs, mats, panels):
        v = np.nanmax(np.abs(m.values)) if m.size and np.isfinite(m.values).any() else 1.0
        im = ax.imshow(m.values, cmap=DIVERGING, aspect="auto",
                       norm=TwoSlopeNorm(vmin=-v, vcenter=0.0, vmax=v))
        ax.set_xticks(range(m.shape[1]))
        ax.set_xticklabels(m.columns, rotation=45, ha="right")
        ax.set_title(unit, fontsize=8, color=INK2, pad=6)
        ax.set_xticks(np.arange(-0.5, m.shape[1], 1), minor=True)
        ax.set_yticks(np.arange(-0.5, len(rows), 1), minor=True)
        ax.grid(which="minor", color="#fcfcfb", linewidth=1.4)
        ax.tick_params(which="minor", length=0)
        fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
    axs[0].set_yticks(range(len(rows)))
    axs[0].set_yticklabels(rows, fontsize=6.5)
    extra = (f"; {len(mods)} modules and the {len(kos)} strongest of {n_ko_all} KOs"
             if n_ko_all > len(kos) else f"; {len(mods)} modules, {len(kos)} KOs")
    fig.suptitle("Features passing an effect floor and BH q < 0.05 on any axis" + extra,
                 fontsize=9.5, y=0.995)
    fig.tight_layout(rect=[0, 0, 1, 0.98])
    out = os.path.join(FIG, "passing_heatmap.png")
    fig.savefig(out, dpi=200)
    print("wrote", out, "rows", len(rows))


def top6():
    keep = pd.read_csv(os.path.join(RES, "passing_modules.tsv"), sep="\t")
    car = pd.read_csv(os.path.join(RES, "carriage_passing.tsv"), sep="\t", index_col=0)
    design = load_design().reindex(car.columns)
    paired = keep[keep["axis"].isin(["compartment", "trip"])].copy()
    if paired.empty:
        print("no compartment or trip feature passed; the paired figure has nothing to "
              "draw (chemistry results are correlations, not pairs)")
        return
    paired["a"] = paired["effect"].abs()
    top = paired.sort_values("a", ascending=False).drop_duplicates("feature").head(6)
    if top.empty:
        print("nothing to plot for top6")
        return

    agg = car.T.groupby([design["site"], design["compartment"], design["trip"]]).median()
    fig, axs = plt.subplots(2, 3, figsize=(10.5, 6.4))
    for ax, (_, r) in zip(axs.ravel(), top.iterrows()):
        f, contrast = r["feature"], r["contrast"]
        s = agg[f]
        if r["axis"] == "compartment":
            w = s.unstack("compartment")
            if contrast == "surface_deep":
                left, right, ln, rn = w.get("Surface"), w.get("Deep"), "Surface", "Deep"
            else:
                left = w.get("PlantRoot")
                right = w[[c for c in ("Surface", "Deep") if c in w]].median(axis=1)
                ln, rn = "PlantRoot", "bulk (S/D)"
            pair = pd.concat([left, right], axis=1).dropna()
            pair.columns = [ln, rn]
            pair = pair.groupby(level="site").median()
        else:
            ta, tb = contrast.split("_vs_")
            w = s.unstack("trip")
            pair = w[[ta, tb]].dropna()
            ln, rn = ta, tb
        for _, row in pair.iterrows():
            up = row.iloc[0] >= row.iloc[1]
            ax.plot([0, 1], row.values, color=RED if up else BLUE,
                    lw=1.0, alpha=0.45, solid_capstyle="round")
        ax.plot([0, 1], [pair.iloc[:, 0].median(), pair.iloc[:, 1].median()],
                color=INK, lw=2.4, marker="o", ms=5, zorder=5)
        ax.set_xticks([0, 1])
        ax.set_xticklabels([ln, rn])
        ax.set_xlim(-0.25, 1.25)
        ax.set_ylabel("carriage, % of mapped community", fontsize=7)
        ax.set_title(short(r.get("name", ""), f, 34) +
                     f"\n{r['axis']}, median diff {r['effect']:.1f} pp, q = {r['q']:.1e}",
                     fontsize=7.5)
        ax.grid(axis="y", color="#eceae5", lw=0.7)
        ax.set_axisbelow(True)
    for ax in axs.ravel()[len(top):]:
        ax.axis("off")
    fig.suptitle("Six strongest paired contrasts; one line per pair, black line is the median",
                 fontsize=9.5)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    out = os.path.join(FIG, "top6_paired.png")
    fig.savefig(out, dpi=200)
    print("wrote", out)


if __name__ == "__main__":
    os.makedirs(FIG, exist_ok=True)
    heatmap()
    top6()
