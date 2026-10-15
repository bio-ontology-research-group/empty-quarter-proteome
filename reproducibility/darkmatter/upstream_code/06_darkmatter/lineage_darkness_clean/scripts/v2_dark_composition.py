#!/usr/bin/env python
"""ROUND TWO, GOAL C: dark composition versus community composition, by phylum.

community  = share of the 5,438 species reps, and median per-sample share of mapped
             relative abundance (W/10_thermal/abundance/rep_relative_abundance.tsv,
             each sample column renormalised to sum 1 over the 5,438 reps)
dark       = share of the single-phylum attributable dark families
             (W/06_darkmatter/darkfam_basis_clean/results/fam_phylum.tsv), and
             share of all dark proteins (results/v2_rep_dark.tsv, native key)

Writes results/v2_dark_composition_by_phylum.tsv and
figures/v2_dark_vs_community_composition.png.
"""
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.abspath(os.path.join(HERE, ".."))
RES = os.path.join(BASE, "results")
FIG = os.path.join(BASE, "figures")
sys.path.insert(0, HERE)
import common as eq

W = eq.W
F_FAM = W + "/06_darkmatter/darkfam_basis_clean/results/fam_phylum.tsv"
F_ABU = W + "/10_thermal/abundance/rep_relative_abundance.tsv"

lines = []


def say(s=""):
    print(s)
    lines.append(s)


# ---- per-rep phylum + native protein counts
reps = pd.read_csv(os.path.join(RES, "v2_reps.tsv"), sep="\t")[
    ["genome", "phylum", "n_proteins", "n_dark"]]
say("GOAL C. DARK COMPOSITION VERSUS COMMUNITY COMPOSITION")
say("=====================================================")
say("species reps %d   phyla %d   proteins %d   dark proteins %d"
    % (len(reps), reps.phylum.nunique(), int(reps.n_proteins.sum()), int(reps.n_dark.sum())))

sp = reps.groupby("phylum").agg(n_species=("genome", "size"),
                                proteins=("n_proteins", "sum"),
                                dark_proteins=("n_dark", "sum")).reset_index()
sp["species_share"] = sp.n_species / sp.n_species.sum()
sp["protein_share"] = sp.proteins / sp.proteins.sum()
sp["dark_protein_share"] = sp.dark_proteins / sp.dark_proteins.sum()

# ---- dark families, single-phylum only
fam = pd.read_csv(F_FAM, sep="\t")
single = fam[fam.n_phyla == 1]
say("dark families attributable %d, single-phylum %d (%.1f%%) - the denominator used"
    % (len(fam), len(single), 100.0 * len(single) / len(fam)))
fc = single.phyla.value_counts().rename_axis("phylum").reset_index(name="dark_families")
fc["dark_family_share"] = fc.dark_families / fc.dark_families.sum()

# ---- abundance, per-sample share of mapped, median across samples
ab = pd.read_csv(F_ABU, sep="\t", index_col=0)
tot = ab.sum(axis=0)
drop = list(tot[tot <= 0].index)
ab = ab.drop(columns=drop)
say("abundance matrix %d reps x %d samples (%d all-zero sample columns dropped: %s)"
    % (ab.shape[0], ab.shape[1], len(drop), ", ".join(drop) if drop else "none"))
ab = ab / ab.sum(axis=0)
ph = reps.set_index("genome")["phylum"].reindex(ab.index)
say("abundance rows without a phylum: %d" % int(ph.isna().sum()))
byph = ab.groupby(ph.values).sum()
med = byph.median(axis=1).rename("median_abundance_share").rename_axis("phylum").reset_index()

# ---- assemble
t = sp.merge(fc, on="phylum", how="outer").merge(med, on="phylum", how="outer")
for c in ["n_species", "dark_families", "proteins", "dark_proteins"]:
    t[c] = t[c].fillna(0).astype(int)
for c in ["species_share", "protein_share", "dark_protein_share",
          "dark_family_share", "median_abundance_share"]:
    t[c] = t[c].fillna(0.0)
t["enrichment_families"] = np.where(t.species_share > 0, t.dark_family_share / t.species_share, np.nan)
t["enrichment_proteins"] = np.where(t.protein_share > 0, t.dark_protein_share / t.protein_share, np.nan)
t = t[["phylum", "n_species", "species_share", "median_abundance_share",
       "dark_family_share", "protein_share", "dark_protein_share",
       "enrichment_families", "enrichment_proteins"]]
t = t.sort_values("n_species", ascending=False)
t.round(6).to_csv(os.path.join(RES, "v2_dark_composition_by_phylum.tsv"), sep="\t", index=False)

say()
say("PHYLA WITH >= 20 SPECIES REPS, sorted by enrichment_families")
big = t[t.n_species >= 20].sort_values("enrichment_families", ascending=False)
say(big.round(4).to_string(index=False))
say()
say("TOP 10 PHYLA BY SPECIES COUNT")
say(t.head(10).round(4).to_string(index=False))
say()
fam_only = t[(t.n_species == 0) & (t.dark_family_share > 0)]
if len(fam_only):
    say("phyla with dark families but no species rep in the 5,438: %s"
        % ", ".join(fam_only.phylum))

# ---- figure
b = big.iloc[::-1]
y = np.arange(len(b))
h = 0.38
fig, ax = plt.subplots(figsize=(10, max(5, 0.42 * len(b))))
ax.barh(y + h / 2, 100 * b.species_share, height=h, color="#6b8fb4", label="community (share of the 5,438 species reps)")
ax.barh(y - h / 2, 100 * b.dark_family_share, height=h, color="#c1543a", label="dark (share of single-phylum dark families)")
ax.set_yticks(y)
ax.set_yticklabels(b.phylum, fontsize=9)
ax.set_xlabel("share of the total (%)")
ax.set_title("Dark-family composition versus community composition, phyla with >=20 species reps")
xmax = 100 * max(b.species_share.max(), b.dark_family_share.max())
ax.set_xlim(0, xmax * 1.22)
for i, (_, r) in enumerate(b.iterrows()):
    ax.text(100 * max(r.species_share, r.dark_family_share) + xmax * 0.015, i,
            "%.2fx" % r.enrichment_families, va="center", fontsize=8, color="0.25")
ax.legend(fontsize=8, loc="lower right")
ax.grid(axis="x", lw=0.4, alpha=0.4)
ax.set_axisbelow(True)
fig.tight_layout()
fig.savefig(os.path.join(FIG, "v2_dark_vs_community_composition.png"), dpi=150)
plt.close(fig)

with open(os.path.join(RES, "v2_composition_report.txt"), "w") as fh:
    fh.write("\n".join(lines) + "\n")
print("\nwrote results/v2_dark_composition_by_phylum.tsv, results/v2_composition_report.txt "
      "and figures/v2_dark_vs_community_composition.png")
