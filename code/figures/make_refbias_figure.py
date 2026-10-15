#!/usr/bin/env python3
"""Section 2.1 figure (supplementary candidate) — the predictor works on its own terms; the ceiling is
in the reference.

DESIGN DECISION (2026-08-23, session 09-display). Two panels on the 448 matched samples / 5,438 catalogue:

  (a) THE PREDICTOR WORKS: per-sample ECDFs of KO presence-recall and rank agreement (Spearman on shared
      KOs). Medians 86.0% and 0.839 are marked. This panel concedes the predictor's competence, which is
      what makes panel (b) a statement about the reference rather than the tool.
  (b) THE TWO CEILINGS: proportion bars. Of protein-coding genes, 55.6% carry no KO at all; of the
      11,374 KOs the genomes do encode, 22.1% (2,517) lie outside the predictor's emittable family set.
      The unreachable share is the focal colour; the annotated/reachable remainder is grey.

Colour: unreachable/dark share #4477aa; reachable remainder #999999; no control set in this figure.
Reads  M/06_reference_bias/refbias_clean/results/{per_sample_comparison,ceiling_summary,ko_universe}.tsv.
Writes figures/fig_refbias_ceiling.{pdf,png}.
Asserted against registered §2.1 numbers: recall median 0.860, Spearman 0.839, NSTI 0.159,
no-KO 55.6% (Bakta denominator), 2,517/11,374 outside the reference, per-sample median 14.0%.
"""
import os, csv, statistics
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.join(HERE, "..", "..")
FIG = os.path.join(REPO, "figures")
WB = "/ibex/scratch/projects/c2014/rob/metagenomics/06_reference_bias/refbias_clean/results"
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
from palette import C_DARK, INK, MUTED, GRID, GREY

ps = list(csv.DictReader(open(os.path.join(WB, "per_sample_comparison.tsv")), delimiter="\t"))
recall = sorted(float(r["recall_of_mag"]) for r in ps)
rho = sorted(float(r["spearman_shared"]) for r in ps)
nsti = statistics.median(float(r["weighted_nsti"]) for r in ps)
outside = statistics.median(float(r["frac_mag_outside_reference"]) for r in ps)
assert len(ps) == 448 and round(statistics.median(recall), 2) == 0.86 \
    and round(statistics.median(rho), 3) == 0.839 and round(nsti, 3) == 0.159 \
    and round(outside, 2) == 0.14, (len(ps), statistics.median(recall), statistics.median(rho), nsti, outside)
ceil = {(r["cohort"], r["denominator"]): float(r["frac_no_ko"])
        for r in csv.DictReader(open(os.path.join(WB, "ceiling_summary.tsv")), delimiter="\t")}
no_ko = ceil[("all", "bakta")]
ko = {r["quantity"]: int(r["value"]) for r in csv.DictReader(open(os.path.join(WB, "ko_universe.tsv")), delimiter="\t")}
assert round(100 * no_ko, 1) == 55.6 and ko["mag_kos_outside_picrust_reference"] == 2517 \
    and ko["mag_kos_detected_in_cohort"] == 11374
print(f"asserted: recall {statistics.median(recall):.3f}, rho {statistics.median(rho):.3f}, "
      f"nsti {nsti:.3f}, no-KO {100*no_ko:.1f}%, outside {ko['mag_kos_outside_picrust_reference']}/11,374")

fig, (a, b) = plt.subplots(1, 2, figsize=(11.6, 4.4))

n = len(recall)
for xs, col, lab, dx, dy, ha in ((recall, C_DARK, "KO presence-recall", 10, -18, "left"),
                                 (rho, GREY, "rank agreement (Spearman, shared KOs)", -8, 10, "right")):
    a.step(xs, [i / n for i in range(1, n + 1)], where="post", color=col, lw=1.8)
    med = statistics.median(xs)
    a.scatter([med], [0.5], s=40, color=col, zorder=5, edgecolor="white", lw=0.8)
    a.annotate(f"median {med:.3f}\n{lab}", (med, 0.5), xytext=(dx, dy),
               textcoords="offset points", fontsize=8, color=col if col == C_DARK else MUTED, ha=ha)
a.text(0.03, 0.95, f"448 matched samples\nmedian weighted NSTI {nsti:.3f}",
       transform=a.transAxes, fontsize=8, color=MUTED, va="top")
a.set_xlabel("per-sample agreement between prediction and encoded content")
a.set_ylabel("cumulative fraction of samples")
a.set_xlim(0.5, 1.0); a.set_ylim(0, 1.02)
a.grid(color=GRID, lw=0.6); a.set_axisbelow(True)
a.set_title("a  Read on its own terms, the predictor works", loc="left", fontsize=10)

rows = [("protein-coding genes\n(19.2M, Bakta annotation)", 100 * no_ko, "55.6% carry no KO"),
        ("KO families the genomes encode\n(11,374 in this cohort)",
         100 * ko["mag_kos_outside_picrust_reference"] / ko["mag_kos_detected_in_cohort"],
         "22.1% outside the predictor's set")]
for y, (lab, pct, note) in zip((1, 0), rows):
    b.barh([y], [pct], height=0.42, color=C_DARK, lw=0)
    b.barh([y], [100 - pct], left=[pct + 0.35], height=0.42, color="#e2e2e2", lw=0)
    b.text(pct + 2.5, y + 0.31, note, fontsize=8.5, color=INK)
b.set_yticks([1, 0]); b.set_yticklabels([r[0] for r in rows], fontsize=8.5)
b.set_xlabel("share unreachable by any orthologue-based predictor (%)")
b.set_xlim(0, 102); b.set_ylim(-0.55, 1.75)
b.grid(axis="x", color=GRID, lw=0.6); b.set_axisbelow(True)
b.set_title("b  What no placement accuracy can reach", loc="left", fontsize=10)
b.text(0.97, 0.06, "median 14.0% of any one sample's\nencoded KOs are outside the set",
       transform=b.transAxes, fontsize=8, color=MUTED, ha="right")

fig.tight_layout()
os.makedirs(FIG, exist_ok=True)
for ext in ("pdf", "png"):
    fig.savefig(os.path.join(FIG, f"fig_refbias_ceiling.{ext}"), bbox_inches="tight", dpi=200)
print("wrote figures/fig_refbias_ceiling.{pdf,png}")
