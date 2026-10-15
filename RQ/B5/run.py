#!/usr/bin/env python3
"""RQ-B5: How is the dark proteome partitioned by depth/compartment?

Partitions ~23k dark-protein families across three sampling compartments
(Surface / Deep / PlantRoot-rhizosphere) using per-sample CoverM abundance.
Per-family Kruskal-Wallis across the 3 compartments + BH-FDR; significant
families are classified by the compartment with the highest mean abundance.

Run from repo root:  python3 RQ/B5/run.py
Stdlib + scipy + matplotlib. Streams the 92 MB abundance matrix line-by-line.
"""
import os
import csv
import math
import sys
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import kruskal

# ---- resolve paths relative to repo root -------------------------------------
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
ECO = os.path.join(ROOT, "results", "darkmatter_ecology")
ENV = os.path.join(ECO, "env_matrix.tsv")
ABUND = os.path.join(ECO, "darkfam_abundance.tsv")
OUTDIR = HERE

COMPARTMENTS = ["Surface", "Deep", "PlantRoot"]
PRESENT = 0.0  # abundance > 0 counts as present

# ---- 1. sample -> compartment ------------------------------------------------
samp2comp = {}
with open(ENV) as fh:
    r = csv.DictReader(fh, delimiter="\t")
    for row in r:
        c = row["compartment"].strip()
        if c in COMPARTMENTS:
            samp2comp[row["sample"]] = c

n_per_comp = defaultdict(int)
for c in samp2comp.values():
    n_per_comp[c] += 1
print("[1] samples per compartment (env_matrix):")
for c in COMPARTMENTS:
    print(f"      {c:10s} {n_per_comp[c]}")
print(f"      TOTAL      {len(samp2comp)}")

# ---- 2. stream abundance, KW per family --------------------------------------
with open(ABUND) as fh:
    header = fh.readline().rstrip("\n").split("\t")
sample_cols = header[1:]
# index of abundance columns that map to a usable compartment
col_idx_by_comp = {c: [] for c in COMPARTMENTS}
for j, s in enumerate(sample_cols, start=1):  # +1 because col0 is family_id
    c = samp2comp.get(s)
    if c:
        col_idx_by_comp[c].append(j)
used_cols = sum(len(v) for v in col_idx_by_comp.values())
print(f"\n[2] abundance sample columns: {len(sample_cols)}; "
      f"mapped to a compartment: {used_cols}")

results = []  # per family: dict
n_total = 0
with open(ABUND) as fh:
    fh.readline()  # skip header
    for line in fh:
        parts = line.rstrip("\n").split("\t")
        fam = parts[0]
        groups = []
        means = {}
        present = {}
        n_obs = {}
        for c in COMPARTMENTS:
            vals = []
            for j in col_idx_by_comp[c]:
                v = parts[j]
                if v == "":
                    continue
                vals.append(float(v))
            groups.append(vals)
            n = len(vals)
            n_obs[c] = n
            means[c] = (sum(vals) / n) if n else 0.0
            present[c] = (sum(1 for x in vals if x > PRESENT) / n) if n else 0.0
        n_total += 1
        # require data in all three groups and some variation
        nonzero_groups = [g for g in groups if any(x > 0 for x in g)]
        if len(nonzero_groups) < 2:
            stat, p = float("nan"), float("nan")
        else:
            try:
                stat, p = kruskal(*groups)
            except ValueError:
                stat, p = float("nan"), float("nan")
        results.append({
            "family_id": fam,
            "mean_Surface": means["Surface"],
            "mean_Deep": means["Deep"],
            "mean_PlantRoot": means["PlantRoot"],
            "prev_Surface": present["Surface"],
            "prev_Deep": present["Deep"],
            "prev_PlantRoot": present["PlantRoot"],
            "kw_H": stat,
            "p": p,
        })

print(f"      families processed: {n_total}")

# ---- BH-FDR over families with a finite p ------------------------------------
tested = [r for r in results if not math.isnan(r["p"])]
tested.sort(key=lambda r: r["p"])
m = len(tested)
for i, r in enumerate(tested, start=1):
    r["q"] = min(1.0, r["p"] * m / i)
# enforce monotonicity
qmin = 1.0
for r in reversed(tested):
    qmin = min(qmin, r["q"])
    r["q"] = qmin
for r in results:
    if math.isnan(r["p"]):
        r["q"] = float("nan")

QCUT = 0.05
sig = [r for r in tested if r["q"] < QCUT]
# classify each significant family by max-mean compartment
enriched = defaultdict(list)
for r in sig:
    mvals = {c: r[f"mean_{c}"] for c in COMPARTMENTS}
    top = max(mvals, key=mvals.get)
    r["enriched_in"] = top
    enriched[top].append(r)

print(f"\n[3] Kruskal-Wallis (3 compartments) + BH-FDR:")
print(f"      families tested (finite p): {m}")
print(f"      significant q<{QCUT}: {len(sig)} "
      f"({100*len(sig)/n_total:.1f}% of {n_total})")
for c in COMPARTMENTS:
    print(f"      enriched in {c:10s}: {len(enriched[c])}")
broad = n_total - len(sig)
print(f"      NOT compartment-specific (broadly distributed): {broad} "
      f"({100*broad/n_total:.1f}%)")

# ---- write per-family table --------------------------------------------------
tab = os.path.join(OUTDIR, "compartment_partition.tsv")
cols = ["family_id", "mean_Surface", "mean_Deep", "mean_PlantRoot",
        "prev_Surface", "prev_Deep", "prev_PlantRoot",
        "kw_H", "p", "q", "enriched_in"]
with open(tab, "w") as out:
    out.write("\t".join(cols) + "\n")
    for r in results:
        r.setdefault("enriched_in", "")
        if "q" not in r:
            r["q"] = float("nan")
        out.write("\t".join(
            (f"{r[c]:.6g}" if isinstance(r.get(c), float) else str(r.get(c, "")))
            for c in cols) + "\n")
print(f"\n[tables] {tab}")

# summary tsv
summ = os.path.join(OUTDIR, "compartment_summary.tsv")
with open(summ, "w") as out:
    out.write("metric\tvalue\n")
    for c in COMPARTMENTS:
        out.write(f"n_samples_{c}\t{n_per_comp[c]}\n")
    out.write(f"n_families_total\t{n_total}\n")
    out.write(f"n_families_tested\t{m}\n")
    out.write(f"n_sig_q<0.05\t{len(sig)}\n")
    for c in COMPARTMENTS:
        out.write(f"n_enriched_{c}\t{len(enriched[c])}\n")
    out.write(f"n_broad\t{broad}\n")
    out.write(f"permanova_compartment_R2\t0.020\n")
print(f"[tables] {summ}")

# ---- 4. figure ---------------------------------------------------------------
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.4))

# panel A: enriched families per compartment
order = ["Surface", "Deep", "PlantRoot"]
labels = ["Surface", "Deep", "Rhizosphere\n(PlantRoot)"]
counts = [len(enriched[c]) for c in order]
colors = ["#e8a33d", "#7b5e3b", "#3d8c5a"]
bars = ax1.bar(labels, counts, color=colors, edgecolor="black", linewidth=0.6)
for b, v in zip(bars, counts):
    ax1.text(b.get_x() + b.get_width() / 2, v, str(v),
             ha="center", va="bottom", fontsize=10)
ax1.set_ylabel("# compartment-enriched dark families")
ax1.set_title(f"A  Compartment-specific dark families\n"
              f"({len(sig)} / {n_total} = {100*len(sig)/n_total:.0f}% KW q<0.05)")
ax1.spines["top"].set_visible(False)
ax1.spines["right"].set_visible(False)

# panel B: example abundance split for the most significant family per compartment
# pick the strongest (lowest q) enriched family in each compartment
examples = []
for c in order:
    fams = sorted(enriched[c], key=lambda r: r["q"])
    if fams:
        examples.append(fams[0])
x = list(range(len(order)))
width = 0.25
for k, r in enumerate(examples):
    yvals = [r[f"mean_{c}"] for c in order]
    short = r["family_id"].split("::")[0][:14]
    ax2.bar([xi + (k - 1) * width for xi in x], yvals, width,
            color=colors[k], edgecolor="black", linewidth=0.4,
            label=f"{r['enriched_in']}: {short}")
ax2.set_xticks(x)
ax2.set_xticklabels(["Surface", "Deep", "Rhizo"])
ax2.set_ylabel("mean CoverM abundance")
ax2.set_title("B  Example enriched families\n(strongest per compartment)")
ax2.legend(fontsize=7, frameon=False)
ax2.spines["top"].set_visible(False)
ax2.spines["right"].set_visible(False)

fig.tight_layout()
pdf = os.path.join(OUTDIR, "fig_compartment.pdf")
png = os.path.join(OUTDIR, "fig_compartment.png")
fig.savefig(pdf)
fig.savefig(png, dpi=150)
print(f"\n[figure] {pdf}")
print(f"[figure] {png}")
print("\nDONE.")
