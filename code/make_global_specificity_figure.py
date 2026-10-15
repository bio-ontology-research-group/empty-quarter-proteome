#!/usr/bin/env python3
"""Figure — database-relative EQ specificity: fraction of 23,382 UniRef50-dark family reps with no detected
homolog in each named global/desert catalogue under one criterion (DIAMOND --sensitive e<1e-3).
Two panels: (a) no-hit % per catalogue grouped by type; (b) no-hit % vs catalogue depth. Persistence at deep
references is descriptive and cannot distinguish geographic restriction from incomplete sampling.
Reads results/global_specificity/{global,desert}_specificity_summary.tsv. Missing rows
(SPIRE / combined still running) are skipped gracefully. Writes figures/fig_darkmatter_global_specificity.*
and results/global_specificity/MULTICATALOG.md.
"""
import os, csv
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.join(HERE, "..")
GS = os.path.join(REPO, "results", "global_specificity")
FIG = os.path.join(REPO, "figures")
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})

NQ = 23382
# label -> (pretty, group, n_target_proteins).
META = {
    "mgnify":           ("MGnify90 (all biomes)",        "global", 718_000_000),
    "gem":              ("GEM (52k MAGs)",               "global", 122_978_191),
    "spire":            ("SPIRE (rep. genomes)",         "global", None),
    "arid":             ("arid soils",                   "desert", 30_569_562),
    "mojave":           ("Mojave",                       "desert", 54_361_006),
    "antarctic_colddesert": ("Antarctic cold-desert",    "desert", 49_206_550),
    "desert_autotroph": ("trace-gas-autotrophy desert",  "desert", 7_200_701),
    "atacama":          ("Atacama",                      "desert", 6_193_359),
    "deserts_combined": ("ALL deserts (union)",          "desert", 147_531_178),
}
# Preserve the whole-million label used by the manuscript and companion result
# record while retaining the exact target count above for plotting/provenance.
DEPTH_LABEL = {"deserts_combined": "147M"}
ROWS = {}

for fn in ("global_specificity_summary.tsv", "desert_specificity_summary.tsv"):
    p = os.path.join(GS, fn)
    if not os.path.exists(p):
        continue
    for r in csv.DictReader(open(p), delimiter="\t"):
        ROWS[r["catalog"]] = dict(n_with_hit=int(r["n_with_hit"]), n_endemic=int(r["n_endemic"]),
                                  pct=float(r["pct_endemic"]))

# order: global first (by depth desc), then deserts (by depth desc; combined emphasised)
present = [k for k in META if k in ROWS]
def depth(k): return META[k][2] or 0
glob = sorted([k for k in present if META[k][1] == "global"], key=depth, reverse=True)
des = sorted([k for k in present if META[k][1] == "desert"], key=depth, reverse=True)
order = glob + des

# ---- table ----
ml = ["# Multi-catalogue EQ-specificity of the dark proteome\n",
      "Fraction of the 23,382 UniRef50-dark family representatives with **no detected homolog** in each named "
      "catalogue under the same criterion (DIAMOND blastp `--sensitive`, e-value <1e-3). Higher means more "
      "database-relative EQ specificity under this search rule. Catalogue depth is shown because shallow "
      "catalogues inflate no-hit fractions.\n",
      "| Catalogue | Type | # target proteins | EQ dark-family no-hits |",
      "|-----------|------|-------------------:|------------------------:|"]
for k in order:
    pretty, grp, n = META[k]
    nd = DEPTH_LABEL.get(k, f"{n/1e6:.0f}M" if n else "—")
    ml.append(f"| {pretty} | {grp} | {nd} | **{ROWS[k]['pct']:.1f}%** ({ROWS[k]['n_endemic']:,}) |")
ml.append("\nThe no-hit fraction remains substantial in the deepest named comparisons — MGnify90 (718M protein "
          "clusters, all biomes), GEM (123M proteins from 52k MAGs)"
          + (", and the union of five desert catalogues (147M proteins)" if "deserts_combined" in ROWS else "")
          + ". This is a database-relative result, not evidence that homologs are geographically absent, and it "
          "cannot distinguish true restriction from incomplete reference sampling. The historical no-hit/hit "
          "dN/dS contrast is descriptive and cannot validate biogeographic restriction.\n")
open(os.path.join(GS, "MULTICATALOG.md"), "w").write("\n".join(ml) + "\n")

# ---- figure ----
COL = {"global": "#4878a8", "desert": "#c97b3c"}
fig, (a, b) = plt.subplots(1, 2, figsize=(12, 5.2), gridspec_kw={"width_ratios": [1.25, 1]})

# (a) bars
yv = range(len(order))
a.barh(list(yv), [ROWS[k]["pct"] for k in order],
       color=[COL[META[k][1]] for k in order], edgecolor="k", lw=0.3)
a.set_yticks(list(yv)); a.set_yticklabels([META[k][0] for k in order], fontsize=9)
a.invert_yaxis()
for i, k in enumerate(order):
    a.text(ROWS[k]["pct"] + 1, i, f"{ROWS[k]['pct']:.0f}%", va="center", fontsize=8)
a.set_xlabel("% of EQ dark families with no detected homolog"); a.set_xlim(0, 100)
a.axvline(50, ls=":", color="grey", lw=0.8)
a.set_title("(a) No detected homolog per named catalogue", loc="left", fontsize=11, weight="bold")
a.barh([], [], color=COL["global"], label="global catalogue")
a.barh([], [], color=COL["desert"], label="other desert")
a.legend(loc="upper right", fontsize=8, frameon=False)

# (b) no-hit fraction vs depth
pts = [(META[k][2], ROWS[k]["pct"], k) for k in order if META[k][2]]
for n, pct, k in pts:
    b.scatter(n, pct, s=60, color=COL[META[k][1]], edgecolor="k", lw=0.4, zorder=5)
    b.annotate(META[k][0].split(" (")[0], (n, pct), fontsize=7, xytext=(5, 3), textcoords="offset points")
b.set_xscale("log"); b.set_xlabel("catalogue depth (# target proteins, log)")
b.set_ylabel("% EQ dark-family no-hits"); b.set_ylim(0, 100)
b.set_title("(b) No-hit fraction depends on named reference", loc="left", fontsize=11, weight="bold")
b.axhline(min(p[1] for p in pts), ls="--", color="grey", lw=0.8,
          label=f"floor ≈ {min(p[1] for p in pts):.0f}% (deepest catalogue)")
b.legend(loc="upper right", fontsize=8, frameon=False)

fig.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(os.path.join(FIG, f"fig_darkmatter_global_specificity.{ext}"), bbox_inches="tight", dpi=200)
print("wrote fig_darkmatter_global_specificity ; catalogues:", ", ".join(order))
print("  " + " | ".join(f"{k}={ROWS[k]['pct']:.1f}%" for k in order))
