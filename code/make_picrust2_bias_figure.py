#!/usr/bin/env python3
"""Figure F5 — reference bias: PICRUSt2 (16S) vs MAG-encoded gene content. Writes pdf+png to ../figures/.
Panels:
  (a) PICRUSt2/MAG KO-richness ratios under the shared-relative-abundance and
      asymmetric MAG-any-presence denominator rules.
  (b) Across-sample abundance correlation for encoded signature functions.
  (c) KO-space ceiling: per-MAG fraction of proteins with NO KO, stratified by GTDB novelty (darker = more novel).
  (d) Weighted NSTI and KO presence-recall across matched samples.
Run: python3 code/make_picrust2_bias_figure.py
"""
import os, csv
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.join(HERE, "..")
PB = os.path.join(REPO, "results", "picrust2_bias")
NOVELTY = os.path.join(REPO, "reproducibility", "historical_mag_novelty_class.tsv")
FIG = os.path.join(REPO, "figures")
os.makedirs(FIG, exist_ok=True)
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})

C = {"core": "#4878a8", "marq": "#d1495b", "fp": "#e08a1e", "novel": "#7a5195"}

# ---------- load ----------
conc = list(csv.DictReader(open(os.path.join(PB, "per_sample_concordance.tsv")), delimiter="\t"))
nsti = np.array([float(r["nsti"]) for r in conc])
recall = np.array([float(r["recall"]) for r in conc])
threshold_rows = list(csv.DictReader(
    open(os.path.join(PB, "rebuild", "threshold_sweep.tsv")), delimiter="\t"
))

# GTDB novelty per MAG (panel c).  This compact historical-R226 snapshot is
# committed and checksum-pinned; figure regeneration must not depend on a
# workstation-only workbench symlink.
tax = {
    r["mag"]: r["novelty_class"].replace(" ", "\n", 1)
    for r in csv.DictReader(open(NOVELTY), delimiter="\t")
}
dark = list(csv.DictReader(open(os.path.join(PB, "dark_fraction.tsv")), delimiter="\t"))
ORDER = ["Known\nspecies", "Novel\nspecies", "Novel\ngenus", "Novel\nfamily+"]
expected_classes = set(ORDER)
unexpected_classes = sorted(set(tax.values()) - expected_classes)
missing_taxonomy = sorted({r["mag"] for r in dark} - set(tax))
if unexpected_classes:
    raise ValueError(f"unexpected novelty classes: {unexpected_classes}")
if missing_taxonomy:
    raise ValueError(f"historical taxonomy snapshot misses {len(missing_taxonomy)} MAGs")
by_nov = {k: [] for k in ORDER}
for r in dark:
    by_nov[tax[r["mag"]]].append(float(r["frac_no_ko"]))
overall_med = np.median([float(r["frac_no_ko"]) for r in dark])

# signature functions (panel b)
marq = list(csv.DictReader(open(os.path.join(PB, "marquee_functions.tsv")), delimiter="\t"))
MSHOW = ["coxL", "coxM", "coxS", "hyaB", "ectC", "ectA", "ectB", "otsA", "otsB", "phr"]
marq = [m for m in marq if m["gene"] in MSHOW]
marq.sort(key=lambda m: float(m["mag_prevalence"]))

# ---------- figure ----------
fig, axes = plt.subplots(2, 2, figsize=(11, 9))
(a, b), (c, d) = axes

# (a) threshold and denominator-rule sensitivity of the richness ratio
threshold_labels = [row["threshold"] for row in threshold_rows]
threshold_labels = ["0" if value == "0" else value.replace("e-0", "e-") for value in threshold_labels]
shared_ratios = [float(row["ratio_picrust_over_mag_relab"]) for row in threshold_rows]
asymmetric_ratios = [float(row["ratio_picrust_over_mag_presence"]) for row in threshold_rows]
xv = np.arange(len(threshold_rows))
a.plot(xv, shared_ratios, marker="o", color=C["fp"], lw=2,
       label="same relative-abundance threshold")
a.plot(xv, asymmetric_ratios, marker="s", color=C["core"], lw=2,
       label="PIC threshold / fixed MAG any-presence")
a.axhline(1, color="grey", ls="--", lw=1, label="parity")
a.annotate(f"{shared_ratios[0]:.2f}×", (xv[0], shared_ratios[0]), xytext=(5, 5),
           textcoords="offset points", ha="left", fontsize=7)
a.annotate(f"{asymmetric_ratios[0]:.2f}×", (xv[0], asymmetric_ratios[0]), xytext=(5, 5),
           textcoords="offset points", ha="left", fontsize=7)
a.set_xticks(xv); a.set_xticklabels(threshold_labels, rotation=35, ha="right")
a.set_xlabel("PICRUSt2 detection threshold")
a.set_ylabel("predicted / encoded KO-richness ratio")
a.set_yscale("log")
a.set_ylim(min(asymmetric_ratios) * 0.65, max(shared_ratios) * 1.45)
a.set_title("(a) Denominator rule dominates KO-richness ratios", loc="left", fontsize=11, weight="bold")
a.legend(loc="lower left", fontsize=7, frameon=False)

# (b) per-function across-sample abundance correlation: can PICRUSt2 rank WHICH samples are enriched?
yvb = np.arange(len(marq))
rhos = [float(m["spearman_abund"]) if m["spearman_abund"] not in ("nan", "") else 0.0 for m in marq]
b.barh(yvb, rhos, color=C["marq"], edgecolor="k", lw=0.3)
b.axvline(0, color="k", lw=0.8)
b.set_yticks(yvb); b.set_yticklabels([m["gene"] for m in marq], fontsize=8)
b.set_xlabel("across-sample abundance Spearman ρ\n(PICRUSt2 vs MAG, per function)")
b.set_xlim(-0.3, 0.5)
b.set_title("(b) Signature-function agreement is heterogeneous", loc="left", fontsize=11, weight="bold")

# (c) dark fraction by novelty
data = [np.array(by_nov[k]) * 100 for k in ORDER]
labels = [f"{k}\n(n={len(by_nov[k]):,})" for k in ORDER]
parts = c.violinplot(data, showmedians=True, showextrema=False)
for pc in parts["bodies"]:
    pc.set_facecolor(C["novel"]); pc.set_alpha(0.55); pc.set_edgecolor("k"); pc.set_linewidth(0.4)
parts["cmedians"].set_color("k")
c.axhline(overall_med * 100, ls="--", color="grey", lw=1, label=f"overall median {overall_med*100:.0f}%")
c.set_xticks(range(1, len(ORDER) + 1)); c.set_xticklabels(labels, fontsize=8)
c.set_ylabel("% of MAG proteins with NO KO"); c.set_ylim(0, 100)
c.set_title("(c) The KO-space ceiling, by GTDB novelty", loc="left", fontsize=11, weight="bold")
c.legend(loc="lower right", fontsize=8, frameon=False)

# (d) placement and KO-core recovery
boxes = d.boxplot([nsti, recall], tick_labels=["weighted NSTI", "KO presence-recall"],
                  patch_artist=True, showfliers=False)
for patch, colour in zip(boxes["boxes"], [C["core"], C["novel"]]):
    patch.set_facecolor(colour)
    patch.set_alpha(0.65)
for median in boxes["medians"]:
    median.set_color("black")
d.set_ylim(0, 1)
d.set_ylabel("per-sample value")
d.text(1, np.median(nsti) + 0.03, f"median {np.median(nsti):.2f}", ha="center", fontsize=8)
d.text(2, np.median(recall) - 0.07, f"median {np.median(recall)*100:.0f}%", ha="center", fontsize=8)
d.set_title("(d) Placement and KO-core recovery remain high", loc="left", fontsize=11, weight="bold")

fig.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(os.path.join(FIG, f"fig_picrust2_bias.{ext}"), bbox_inches="tight", dpi=200)
print(f"wrote fig_picrust2_bias  matched={len(conc)}")
