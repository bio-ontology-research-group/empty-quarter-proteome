#!/usr/bin/env python3
"""Section 3.1 figure — the two controls behind the absence claim.

DESIGN DECISION (2026-08-17). The per-catalogue absence percentages live in tab:specificity, not here: a
table presents nine rows of exact numbers better than bars, and an earlier draft of this figure duplicated
that table. This figure therefore carries only what a table cannot show, one panel per control:

  (a) observed absence against the absence expected from catalogue size alone (e-value-corrected
      rarefaction of MGnify90). "Below the diagonal" makes content-beyond-depth legible at a glance and is
      otherwise hard to convey. MGnify90 sits on the diagonal by construction, which is the internal check.
      Atacama is drawn open: its control gap of -29.97 shows a single 35-sample study samples protein space
      too narrowly to beat a random draw even for conserved proteins, so its position says nothing about
      desert content.
  (b) median site occupancy of the reference-absent families against those with a homolog, within length
      strata. Shows "widespread, and marginally less so, never confined" as a consistent pattern across
      five strata rather than as a pair of medians in prose. Drawn as dumbbells because the gap is the
      message. The x axis starts at zero deliberately: the gap is small and truncating the axis would
      overstate it.

This supersedes make_global_specificity_figure.py (23,382-orthogroup basis) and replaces the cut
fig:endemism, whose other two panels showed a dropped analysis and a duplicate of the table.

Colour encodes the query set, never the catalogue type: focal set #4477aa, comparator #cc6677. Checked for
colour-vision separation before use (OKLab dE: normal 21.3, deuteranopia and protanopia 18.5, against a >=8
CVD floor and >=15 normal-vision floor). A grey comparator was rejected: deuteranopia dE 4.4.

Reads results/specificity_clean/{content_vs_depth,b6_absence_occupancy_by_length}.tsv.
Writes figures/fig_specificity_controls.{pdf,png}.
"""
import os, csv
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.join(HERE, "..", "..")
SRC = os.path.join(REPO, "results", "specificity_clean")
FIG = os.path.join(REPO, "figures")
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})

from palette import C_FOCAL, C_COMP
INK, MUTED, GRID = "#222222", "#666666", "#dddddd"
N_SITES = 61

PRETTY = {
    "atacama": "Atacama", "desert_autotroph": "Desert autotrophy", "arid": "Arid soils",
    "antarctic_colddesert": "Antarctic cold desert", "mojave": "Mojave",
    "deserts_combined": "Five deserts, combined", "gem": "GEM", "spire": "SPIRE representatives",
    "mgnify90": "MGnify90",
}
STRATUM = {"lt60aa": "< 60", "60-99": "60–99", "100-149": "100–149",
           "150-249": "150–249", "ge250aa": "≥ 250"}

def read(fn):
    p = os.path.join(SRC, fn)
    if not os.path.exists(p):
        raise SystemExit("missing input: " + p)
    with open(p) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))

gapr = {r["catalogue"]: (float(r["dark_obs"]), float(r["dark_exp_size_only"])) for r in read("content_vs_depth.tsv")}
occ = read("b6_absence_occupancy_by_length.tsv")

fig, (a, b) = plt.subplots(1, 2, figsize=(11.6, 5.0))

# ---- (a) observed vs size-matched expectation -------------------------------------------------------
a.plot([30, 105], [30, 105], ls="--", lw=0.9, color="#999999", zorder=1)
a.text(37.5, 33.0, "no content beyond size", rotation=38, fontsize=7.5, color=MUTED,
       ha="left", va="center", rotation_mode="anchor")
OVERRIDE = {  # crowded middle; defaults ran a label across a neighbouring marker (found by rendering)
    "antarctic_colddesert": (-8, 10, "right", "bottom"),
    "mojave": (9, -9, "left", "top"),
}
for k, (obs, exp) in sorted(gapr.items(), key=lambda kv: kv[1][1]):
    open_marker = (k == "atacama")
    a.scatter(exp, obs, s=64, zorder=5, facecolor="white" if open_marker else C_FOCAL,
              edgecolor=C_FOCAL, lw=1.4 if open_marker else 0.5)
    if k in OVERRIDE:
        dx, dy, ha, va = OVERRIDE[k]
    else:
        right = exp > 74
        dx, dy, ha, va = (-8 if right else 8), -9, ("right" if right else "left"), "top"
    a.annotate(PRETTY.get(k, k), (exp, obs), textcoords="offset points", xytext=(dx, dy),
               fontsize=7.5, color=INK, ha=ha, va=va)
a.set_xlabel("absence expected from catalogue size alone (%)")
a.set_ylabel("absence observed (%)")
a.set_xlim(33, 100); a.set_ylim(28, 104)
a.grid(color=GRID, lw=0.6, zorder=0); a.set_axisbelow(True)
a.set_title("a  Size does not explain the absence", loc="left", fontsize=10)
a.scatter([], [], s=64, facecolor="white", edgecolor=C_FOCAL, lw=1.4,
          label="narrow single study, not interpretable")
a.legend(loc="upper left", frameon=False, fontsize=8)

# ---- (b) occupancy dumbbells within length strata ---------------------------------------------------
ys = list(range(len(occ)))[::-1]          # shortest stratum at the top
for y, r in zip(ys, occ):
    xa, xh = float(r["median_sites_absent"]), float(r["median_sites_has_hit"])
    b.plot([xa, xh], [y, y], color="#bbbbbb", lw=2.0, zorder=2, solid_capstyle="round")
    b.scatter([xh], [y], s=58, color=C_COMP, edgecolor="white", lw=0.8, zorder=4)
    b.scatter([xa], [y], s=58, color=C_FOCAL, edgecolor="white", lw=0.8, zorder=5)
    b.text(min(xa, xh) - 1.6, y, f'{xa:.0f}', ha="right", va="center", fontsize=8, color=INK)
    b.text(max(xa, xh) + 1.6, y, f'{xh:.0f}', ha="left", va="center", fontsize=8, color=MUTED)
b.axvline(N_SITES, ls=":", lw=0.9, color="#999999", zorder=1)
b.text(N_SITES - 1.0, 1.6, "all 61 sites", rotation=90, fontsize=7.5, color=MUTED,
       ha="right", va="center")
b.set_yticks(ys)
b.set_yticklabels([f'{STRATUM.get(r["len_stratum"], r["len_stratum"])}\n'
                  f'n={int(r["n_absent"]):,} / {int(r["n_has_hit"]):,}' for r in occ], fontsize=8)
b.set_ylabel("family representative length (aa)")
b.set_xlabel("median number of sites occupied (of 61)")
b.set_xlim(0, 66); b.set_ylim(-0.7, len(occ) - 0.3)
b.grid(axis="x", color=GRID, lw=0.6, zorder=0); b.set_axisbelow(True)
b.set_title("b  Nor is local rarity", loc="left", fontsize=10)
b.scatter([], [], s=58, color=C_FOCAL, label="absent from every reference")
b.scatter([], [], s=58, color=C_COMP, label="has a homolog")
b.legend(loc="lower left", frameon=False, fontsize=8)

fig.tight_layout()
os.makedirs(FIG, exist_ok=True)
for ext in ("pdf", "png"):
    fig.savefig(os.path.join(FIG, f"fig_specificity_controls.{ext}"), bbox_inches="tight", dpi=200)
print(f"wrote figures/fig_specificity_controls.{{pdf,png}} — {len(gapr)} catalogues, {len(occ)} length strata")
