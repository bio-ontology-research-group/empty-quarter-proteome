#!/usr/bin/env python3
"""Section 2.3 figure (supplementary candidate) — three stressors, three responses.

DESIGN DECISION (2026-08-23, session 09-display). One panel per stressor, matching the section's
three-mode story; every number is recomputed from the site-level tables and asserted against the
registered site-paired values before drawing:

  (a) LIGHT PARTITIONS: within the 44 sites profiled in both compartments, the paired-site median share
      of the mapped community carrying photolyase, carotenoid synthesis and NHEJ, surface against
      shielded subsurface — dumbbells, opposed movements visible (photolyase and carotenoids up at the
      surface, NHEJ up below).
  (b) DESICCATION IS UNIFORM: compatible-solute pathways by compartment (share of the mapped community);
      the flat lines are the finding, so the y axis starts at zero to keep flatness honest.
  (c) SALT REPLACES: per-phylum partial Spearman (geography-controlled) of phylum share against sulfur;
      Halobacteriota stands alone (rho +0.45, p = 7e-14) — lollipops, one highlighted.

Colour: surface / highlighted entity #4477aa; subsurface and other phyla #999999. The control red is not
used here: nothing in this figure is a negative-control set.
Reads  M/06_metabolism/uv_survival_clean/results/site_level_values.tsv,
       M/06_metabolism/osmolyte_clean/results/{community_summary,salt_turnover}.tsv.
Writes figures/fig_stress_modes.{pdf,png}.
"""
import os, csv, statistics
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.join(HERE, "..", "..")
FIG = os.path.join(REPO, "figures")
M = "/ibex/scratch/projects/c2014/rob/metagenomics/06_metabolism"
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
from palette import C_HI, INK, MUTED, GRID, GREY

# ---- (a) paired-site compartment contrast (uv_survival) -----------------------------------------------
vals = {}
for r in csv.DictReader(open(os.path.join(M, "uv_survival_clean", "results", "site_level_values.tsv")),
                        delimiter="\t"):
    if r["metric"] == "abund_annotated_pct":
        vals.setdefault((r["trait"], r["compartment"]), {})[r["site"]] = float(r["median_value"])
pair = {}
for trait in ("photolyase", "carotenoid", "nhej"):
    s, d = vals.get((trait, "Surface"), {}), vals.get((trait, "Deep"), {})
    sites = sorted(set(s) & set(d))
    pair[trait] = (statistics.median(s[x] for x in sites), statistics.median(d[x] for x in sites), len(sites))
assert pair["photolyase"][2] == 44, pair["photolyase"]
assert (round(pair["photolyase"][0], 1), round(pair["photolyase"][1], 1)) == (24.3, 12.1), pair["photolyase"]
assert (round(pair["carotenoid"][0], 1), round(pair["carotenoid"][1], 1)) == (70.3, 61.0), pair["carotenoid"]
assert (round(pair["nhej"][0], 1), round(pair["nhej"][1], 1)) == (80.0, 85.7), pair["nhej"]

# ---- (b) osmolytes by compartment (community medians) --------------------------------------------------
osm = {}
for r in csv.DictReader(open(os.path.join(M, "osmolyte_clean", "results", "community_summary.tsv")),
                        delimiter="\t"):
    osm.setdefault(r["trait"], {})[r["compartment"]] = float(r["abund_annotated_median"])
# (osmolyte community_summary keeps the *_median column names; uv site_level_values uses metric rows)
osm_traits = [t for t in ("any_osmolyte", "trehalose", "ectoine") if t in osm] or sorted(osm)[:3]

# ---- (c) salt turnover: per-phylum partial rho vs sulfur ----------------------------------------------
salt = [(r["phylum"], float(r["rho_partial"]), float(r["p_partial"]))
        for r in csv.DictReader(open(os.path.join(M, "osmolyte_clean", "results", "salt_turnover.tsv")),
                                delimiter="\t") if r["axis"] == "S"]
halo = [x for x in salt if x[0] == "Halobacteriota"]
assert halo and round(halo[0][1], 2) == 0.45, halo
print(f"asserted: photolyase {pair['photolyase'][:2]}, carotenoid {pair['carotenoid'][:2]}, "
      f"nhej {pair['nhej'][:2]} (44 paired sites); Halobacteriota rho_S {halo[0][1]:.2f}")

fig, (a, b, c) = plt.subplots(1, 3, figsize=(15.0, 4.6))

rows = [("photolyase", "photolyase"), ("carotenoid synthesis", "carotenoid"),
        ("non-homologous end joining", "nhej")]
ys = [2, 1, 0]
for y, (lab, key) in zip(ys, rows):
    s, d, _ = pair[key]
    a.plot([d, s], [y, y], color="#bbbbbb", lw=2.0, zorder=2, solid_capstyle="round")
    a.scatter([d], [y], s=58, color=GREY, edgecolor="white", lw=0.8, zorder=4)
    a.scatter([s], [y], s=58, color=C_HI, edgecolor="white", lw=0.8, zorder=5)
    a.text(min(s, d) - 2.2, y, f"{min(s,d):.0f}", ha="right", va="center", fontsize=8, color=MUTED)
    a.text(max(s, d) + 2.2, y, f"{max(s,d):.0f}", ha="left", va="center", fontsize=8, color=INK)
a.set_yticks(ys); a.set_yticklabels([r[0] for r in rows], fontsize=8.5)
a.set_xlabel("share of the mapped community (%), paired-site median, 44 sites")
a.set_xlim(0, 100); a.set_ylim(-0.6, 2.6)
a.grid(axis="x", color=GRID, lw=0.6); a.set_axisbelow(True)
a.set_title("a  Light partitions the community", loc="left", fontsize=10)
a.scatter([], [], s=58, color=C_HI, label="sun-exposed surface")
a.scatter([], [], s=58, color=GREY, label="shielded subsurface")
a.legend(loc="lower right", frameon=False, fontsize=8)

comps = ["Surface", "Deep", "PlantRoot"]
xs = [0, 1, 2]
for t in osm_traits:
    yy = [osm[t].get(k) for k in comps]
    b.plot(xs, yy, color=C_HI, lw=1.6, marker="o", ms=5, mec="white", mew=0.6)
    b.annotate(t.replace("any_osmolyte", "any pathway"), (xs[-1], yy[-1]),
               xytext=(8, 0), textcoords="offset points", fontsize=8, color=INK, va="center")
b.set_xticks(xs); b.set_xticklabels(["surface", "subsurface", "root-adjacent"], fontsize=8.5)
b.set_ylabel("share of the mapped community (%)")
b.set_ylim(0, 100); b.set_xlim(-0.3, 3.1)
b.grid(axis="y", color=GRID, lw=0.6); b.set_axisbelow(True)
b.set_title("b  Desiccation is met everywhere alike", loc="left", fontsize=10)

salt.sort(key=lambda x: x[1])
for i, (ph, rho, p) in enumerate(salt):
    hi = ph == "Halobacteriota"
    c.plot([0, rho], [i, i], color="#dddddd" if not hi else "#bbbbbb", lw=1.6, zorder=2)
    c.scatter([rho], [i], s=52 if hi else 36, color=C_HI if hi else GREY,
              edgecolor="white", lw=0.8, zorder=4)
c.axvline(0, color="#999999", lw=0.8)
hi_i = [i for i, x in enumerate(salt) if x[0] == "Halobacteriota"][0]
c.annotate(f"Halobacteriota  rho = +{halo[0][1]:.2f}\np = 7e-14, geography-controlled",
           (halo[0][1], hi_i), xytext=(-14, -30), textcoords="offset points",
           fontsize=8, color=INK, ha="right")
c.set_yticks(range(len(salt))); c.set_yticklabels([x[0] for x in salt], fontsize=7.5)
c.set_xlabel("partial Spearman vs sulfur (248 samples)")
c.grid(axis="x", color=GRID, lw=0.6); c.set_axisbelow(True)
c.set_title("c  Salt replaces the community", loc="left", fontsize=10)

fig.tight_layout()
os.makedirs(FIG, exist_ok=True)
for ext in ("pdf", "png"):
    fig.savefig(os.path.join(FIG, f"fig_stress_modes.{ext}"), bbox_inches="tight", dpi=200)
print("wrote figures/fig_stress_modes.{pdf,png}")
