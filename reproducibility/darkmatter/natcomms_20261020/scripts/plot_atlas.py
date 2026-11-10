#!/usr/bin/env python3
"""Render manuscript figures from compact, retained numerical results."""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

P = Path(__file__).resolve().parents[2]
R = Path(__file__).resolve().parents[1] / "results"
parser = argparse.ArgumentParser()
parser.add_argument("--out", type=Path, required=True)
args = parser.parse_args()
args.out.mkdir(exist_ok=True, parents=True)
plt.rcParams.update({"font.size": 9, "axes.spines.top": False,
                     "axes.spines.right": False, "pdf.fonttype": 42,
                     "axes.labelcolor": "#253444", "text.color": "#253444"})
TEAL, GREY, ORANGE = "#187f83", "#73849b", "#d18948"


def save(fig, name):
    fig.savefig(args.out / f"{name}.pdf", bbox_inches="tight")
    fig.savefig(args.out / f"{name}.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


def panel(ax, letter, title):
    ax.set_title(title, loc="left", fontsize=10, pad=12)
    ax.text(-.12, 1.06, letter, transform=ax.transAxes, weight="bold", fontsize=12)


composition = pd.read_csv(P / "inputs/06_darkmatter/lineage_darkness_clean/results/v2_dark_composition_by_phylum.tsv.gz", sep="\t")
composition = composition.nlargest(6, "dark_family_share").sort_values("dark_family_share")
sites = pd.read_csv(R / "ecology_site_contrasts.tsv", sep="\t")
sites = sites[(sites.metric == "thermo_read_percent") & (sites.contrast == "Deep-Surface")]
metabolism = pd.read_csv(P / "round4_corrections/inputs/v2_thermo_metabolism_summary.tsv.gz", sep="\t")
metabolism = metabolism[metabolism.group == "thermo_237"].set_index("feature")
features = [
    ("M00001", "Glycolysis module"),
    ("M00012", "Glyoxylate-cycle module"),
    ("branched_chain_aa_ABC", "Branched-chain amino-acid transport"),
    ("peptide_nickel_ABC", "Peptide/nickel ABC transport"),
    ("oligopeptide_POT", "Oligopeptide transport"),
    ("sugar_ABC_msm", "Selected sugar ABC transport"),
]
assert len(sites) == 37
fig = plt.figure(figsize=(7.6, 6.7), layout="constrained")
gs = fig.add_gridspec(2, 2, width_ratios=[1.25, 1], height_ratios=[1.2, 1])
ax = fig.add_subplot(gs[0, 0])
y = np.arange(len(composition))
ax.barh(y + .18, 100 * composition.species_share, .32, color=GREY, label="Species")
ax.barh(y - .18, 100 * composition.dark_family_share, .32, color=TEAL, label="Dark families")
ax.set(yticks=y, yticklabels=composition.phylum, xlabel="Share of catalogue (%)", xlim=(0, 40))
ax.tick_params(axis="y", labelsize=8)
ax.legend(frameon=False, loc="lower right", fontsize=8)
panel(ax, "a", "Contribution to family diversity")
ax = fig.add_subplot(gs[0, 1])
for _, row in sites.iterrows():
    ax.plot([0, 1], [row.Surface, row.Deep], color=GREY, alpha=.40, lw=.8)
ax.scatter(np.zeros(len(sites)), sites.Surface, color=GREY, s=13, zorder=3)
ax.scatter(np.ones(len(sites)), sites.Deep, color=TEAL, s=13, zorder=3)
ax.set(xticks=[0, 1], xticklabels=["Surface", "Subsurface"], xlim=(-.2, 1.2),
       ylabel="Thermoplasmatota (% of all reads)")
ax.text(.03, .96, "37 sites · 69 paired visits", transform=ax.transAxes, va="top", fontsize=8)
panel(ax, "b", "Within-visit depth contrast")
ax = fig.add_subplot(gs[1, :])
values = [100 * metabolism.loc[k, "share"] for k, _ in features]
y = np.arange(len(features))
ax.barh(y, values, color=TEAL, height=.55)
ax.set(yticks=y, yticklabels=[s for _, s in features], xlabel="Representatives with marker set (%)",
       xlim=(0, 100))
ax.invert_yaxis()
for pos, value in zip(y, values):
    ax.text(value + 1.5, pos, f"{value:.0f}", va="center", fontsize=8)
panel(ax, "c", "Encoded carbon metabolism and transport (237 genomes)")
save(fig, "fig_archaeal_ecology")

matched = pd.read_csv(P / "round7_final/h3/matched_full_census_summary.tsv", sep="\t")
matched = matched[matched.run_label == "harmonized_seed20261014"].set_index("cohort")
depth = pd.read_csv(P / "inputs/06_darkmatter/specificity_clean/results/content_vs_depth.tsv.gz", sep="\t")
occupancy = pd.read_csv(P / "inputs/06_darkmatter/specificity_clean/results/b6_absence_occupancy_by_length.tsv.gz", sep="\t")
fig, axes = plt.subplots(1, 3, figsize=(8.2, 3.3), layout="constrained")
ax = axes[0]
for i, (field, label, colour) in enumerate([("smag_hit_pct", "Soil (SMAG)", TEAL),
                                          ("gem_hit_pct", "Global (GEM)", GREY)]):
    ax.bar(np.arange(2) + (i - .5) * .32,
           matched.loc[["Dark", "Control"], field], width=.30, color=colour, label=label)
ax.set(xticks=[0, 1], xticklabels=["Dark", "Control"], ylabel="Queries with a hit (%)", ylim=(0, 105))
ax.legend(frameon=False, fontsize=7, loc="upper left", bbox_to_anchor=(0, .78))
panel(ax, "a", "Matched reference depth")
ax = axes[1]
for prefix, colour, label in [("dark", TEAL, "Dark"), ("ctrl", GREY, "Control")]:
    ax.scatter(depth[prefix + "_exp_size_only"], depth[prefix + "_obs"], s=23, color=colour, label=label)
ax.plot([0, 100], [0, 100], color=".65", lw=.8)
for name, label, offset in [("deserts_combined", "Deserts", (5, -10)),
                             ("spire", "SPIRE", (-29, -14)),
                             ("mgnify90", "MGnify90", (-44, 9))]:
    row = depth.set_index("catalogue").loc[name]
    ax.annotate(label, (row.dark_exp_size_only, row.dark_obs), xytext=offset,
                textcoords="offset points", fontsize=6.5, color=TEAL)
ax.set(xlabel="Expected absence (%)", ylabel="Observed absence (%)", xlim=(0, 100), ylim=(0, 100))
ax.legend(frameon=False, fontsize=7)
panel(ax, "b", "Reference composition")
ax = axes[2]
for col, colour, label in [("median_sites_absent", TEAL, "No reference hit"),
                           ("median_sites_has_hit", GREY, "Reference hit")]:
    ax.plot(np.arange(len(occupancy)), occupancy[col], "o-", color=colour, ms=4, label=label)
length_labels = {"lt60aa": "<60", "60-99": "60–99", "100-149": "100–149",
                 "150-249": "150–249", "ge250aa": "≥250"}
ax.set(xticks=np.arange(len(occupancy)), xticklabels=[length_labels.get(x, x) for x in occupancy.len_stratum],
       xlabel="Representative length (aa)", ylabel="Median occupied locations", ylim=(0, 61))
ax.tick_params(axis="x", rotation=35, labelsize=7)
ax.legend(frameon=False, fontsize=7, loc="lower left")
panel(ax, "c", "Geographic recurrence")
save(fig, "fig_reference_atlas")
