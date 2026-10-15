#!/usr/bin/env python3
"""
RQ-E2: What hypotheses emerge from the historical dark-family dN/dS>1 screen tail?

Run from repo root:  python3 RQ/E2/run.py
Inputs:  results/darkmatter_ecology/positive_selection/{positive_families,struct_function_confident,
         struct_master,boltz_affinity}.tsv
Outputs: RQ/E2/table_*.tsv, RQ/E2/fig_posel_function.{pdf,png}
Stdlib + matplotlib only. REAL data; nothing invented.
"""
import csv, os, statistics
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

IN = "results/darkmatter_ecology/positive_selection"
OUT = "RQ/E2"
os.makedirs(OUT, exist_ok=True)


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


# ---------------------------------------------------------------- load
fam = read_tsv(os.path.join(IN, "positive_families.tsv"))
master = read_tsv(os.path.join(IN, "struct_master.tsv"))
conf = read_tsv(os.path.join(IN, "struct_function_confident.tsv"))
boltz = read_tsv(os.path.join(IN, "boltz_affinity.tsv"))

# ================================================================ TASK 1
n_fam = len(fam)
robust = [r for r in fam if r["robust"] == "True"]
dnds = [float(r["dnds"]) for r in fam]
members = [int(r["n_members"]) for r in fam]
nphyla = [int(r["n_phyla"]) for r in fam]
cconf = [float(r["Cconf"]) for r in fam if r["Cconf"] not in ("", "NA")]
env_track = [r for r in fam if r["env"].strip() != ""]

phylum_ct = {}
for r in fam:
    phylum_ct[r["dom_phylum"]] = phylum_ct.get(r["dom_phylum"], 0) + 1
phylum_ct = dict(sorted(phylum_ct.items(), key=lambda kv: -kv[1]))

with open(os.path.join(OUT, "table_positive_families_summary.tsv"), "w") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["metric", "value"])
    w.writerow(["n_dark_families_dNdS_gt1", n_fam])
    w.writerow(["n_robust (mean_dS>0.05 & n_pairs>=3)", len(robust)])
    w.writerow(["dNdS_max", f"{max(dnds):.3f}"])
    w.writerow(["dNdS_median", f"{statistics.median(dnds):.3f}"])
    w.writerow(["dNdS_min", f"{min(dnds):.3f}"])
    w.writerow(["n_members_median", int(statistics.median(members))])
    w.writerow(["n_members_max", max(members)])
    w.writerow(["Cconf_median (calibrated EC confidence)", f"{statistics.median(cconf):.3f}"])
    w.writerow(["n_env_tracking (indicator gradient)", len(env_track)])
    w.writerow(["n_phyla_represented", len(phylum_ct)])
    for p, c in phylum_ct.items():
        w.writerow([f"phylum::{p}", c])

# ================================================================ TASK 2
# Descriptive functional categorisation of the confident structural homologs.
# Category assignment is by the FoldSeek target's known function (manual, from fs_function string).
CATEGORY = {
    "V5PRr1__self_SemiBin_543::V5PRr1_235427_14": ("Uncharacterised conserved (DUF/orphan)", "other"),
    "F2Sr2__self_SemiBin_480::F2Sr2_109584_65": ("Toxin-like peptide (tamapin scorpion toxin fold)", "recognition"),
    "F38PRr1__self_SemiBin_216::F38PRr1_858279_42": ("Uncharacterised conserved (DUF/orphan)", "other"),
    "36PRr1__single_SemiBin_29_sub::36PRr1_1482479_11": ("Uncharacterised conserved (DUF/orphan)", "other"),
    "58PRr3__self_SemiBin_361::58PRr3_401505_1": ("Carbohydrate-binding lectin (Photorhabdus PHL)", "recognition"),
    "V12Dr3__self_SemiBin_27_sub::V12Dr3_1574070_5": ("DNA-binding regulator (Arg repressor)", "recognition"),
    "S11PRr1__self_SemiBin_3255::S11PRr1_690820_66": ("Solute uptake (ABC transporter NosDFY)", "recognition"),
    "F4PRr2__self_SemiBin_1826::F4PRr2_783354_7": ("Lipid/cell-surface carrier (LprG)", "recognition"),
}
master_by_fam = {r["family"]: r for r in master}

rows_cat = []
for r in conf:
    fam_id = r["family"]
    cat, group = CATEGORY.get(fam_id, ("Unclassified", "other"))
    m = master_by_fam.get(fam_id, {})
    rows_cat.append({
        "family": fam_id,
        "n_members": r["n_members"],
        "dnds": r["dnds"],
        "plddt": r["plddt"],
        "pocket_prob": r["pocket_prob"],
        "n_pos_sites": r["n_pos_sites"],
        "fs_tm": r["fs_tm"],
        "fs_target": r["fs_target"],
        "function_category": cat,
        "functional_hypothesis_group": group,
        "structural_homolog": r["function"],
    })

with open(os.path.join(OUT, "table_functional_categories.tsv"), "w") as fh:
    cols = ["family", "n_members", "dnds", "plddt", "pocket_prob", "n_pos_sites",
            "fs_tm", "fs_target", "function_category", "functional_hypothesis_group", "structural_homolog"]
    w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
    w.writeheader()
    for r in rows_cat:
        w.writerow(r)

# structural-decodability tiers across the full historical dN/dS>1 screen tail
n_homolog = sum(1 for r in master if r["has_struct_homolog"] == "1")
n_confident = sum(1 for r in master if r["fs_confident"] == "1")
n_dark = len(master) - n_homolog
n_homolog_lowconf = n_homolog - n_confident
n_recognition = sum(1 for r in rows_cat if r["functional_hypothesis_group"] == "recognition")
n_other = n_confident - n_recognition

# Nominal FUBAR-site overlap with predicted pockets, recomputed from master.
# This tally is descriptive; it is not an interface-enrichment test.
pocket_fam = [r for r in master if float(r["pocket_prob"]) >= 0.5]
sites_total = sum(int(r["n_pos_sites"]) for r in pocket_fam)
in_pocket = sum(int(r["pos_in_pocket"]) for r in pocket_fam)

with open(os.path.join(OUT, "table_structural_tiers.tsv"), "w") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["tier", "n_families"])
    w.writerow(["total_dnds_gt1_screen_tail", len(master)])
    w.writerow(["structurally_dark (no PDB/SwissProt/AFDB homolog)", n_dark])
    w.writerow(["structural_homolog_but_not_confident", n_homolog_lowconf])
    w.writerow(["confident_homolog (prob>=.9,TM>=.5,pLDDT>=.5)", n_confident])
    w.writerow(["  of confident: recognition/binding/surface/transport hypotheses", n_recognition])
    w.writerow(["  of confident: uncharacterised/other", n_other])
    w.writerow(["families with predicted pocket", len(pocket_fam)])
    w.writerow(["nominal FUBAR-positive sites in predicted-pocket families", sites_total])
    w.writerow(["nominal FUBAR-positive sites inside predicted pockets", in_pocket])

# ================================================================ TASK 3 (Boltz)
with open(os.path.join(OUT, "table_boltz_hypotheses.tsv"), "w") as fh:
    cols = list(boltz[0].keys())
    w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
    w.writeheader()
    for r in boltz:
        w.writerow(r)

# ================================================================ FIGURE
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.6))

# Panel A: structural-decodability tiers of the 148
labels_a = ["Structurally dark\n(no homolog)",
            "Homolog,\nlow confidence",
            "Confident\nstructural homolog"]
vals_a = [n_dark, n_homolog_lowconf, n_confident]
colors_a = ["#444444", "#c0a060", "#2a7fb8"]
bars = ax1.bar(labels_a, vals_a, color=colors_a, edgecolor="black", linewidth=0.6)
for b, v in zip(bars, vals_a):
    ax1.text(b.get_x() + b.get_width() / 2, v + 1.5, str(v), ha="center", fontsize=11, fontweight="bold")
ax1.set_ylabel("dN/dS>1 screen-tail families")
ax1.set_title(f"A. Structural decodability of the {len(master)} dN/dS>1 families", fontsize=10.5)
ax1.set_ylim(0, max(vals_a) * 1.18)
ax1.tick_params(axis="x", labelsize=8.5)

# Panel B: descriptive functional category of the confident set
cat_ct = {}
cat_grp = {}
for r in rows_cat:
    cat_ct[r["function_category"]] = cat_ct.get(r["function_category"], 0) + 1
    cat_grp[r["function_category"]] = r["functional_hypothesis_group"]
items = sorted(cat_ct.items(), key=lambda kv: (cat_grp[kv[0]] != "recognition", -kv[1]))
clabels = [k for k, _ in items]
cvals = [v for _, v in items]
ccolors = ["#d1495b" if cat_grp[k] == "recognition" else "#8d99ae" for k in clabels]
ypos = range(len(clabels))
ax2.barh(list(ypos), cvals, color=ccolors, edgecolor="black", linewidth=0.6)
ax2.set_yticks(list(ypos))
ax2.set_yticklabels(clabels, fontsize=8)
ax2.invert_yaxis()
ax2.set_xlabel("Families (confident structural homolog)")
ax2.set_title(f"B. Structural annotations include recognition/binding hypotheses ({n_recognition}/{n_confident})", fontsize=9)
ax2.set_xlim(0, max(cvals) + 0.8)
for i, v in enumerate(cvals):
    ax2.text(v + 0.05, i, str(v), va="center", fontsize=9)
legend = [Patch(facecolor="#d1495b", edgecolor="black", label="recognition / binding / surface / transport hypothesis"),
          Patch(facecolor="#8d99ae", edgecolor="black", label="uncharacterised / other")]
ax2.legend(handles=legend, fontsize=7.5, loc="upper right", frameon=False)

fig.suptitle("Structural annotations in the historical dark-family dN/dS>1 screen",
             fontsize=12, fontweight="bold")
fig.tight_layout(rect=(0, 0, 1, 0.95))
fig.savefig(os.path.join(OUT, "fig_posel_function.pdf"))
fig.savefig(os.path.join(OUT, "fig_posel_function.png"), dpi=200)

# ---------------------------------------------------------------- console
print(f"[E2] {n_fam} dark families dN/dS>1 ({len(robust)} robust); {len(env_track)} env-tracking")
print(f"[E2] structural tiers: dark={n_dark} lowconf={n_homolog_lowconf} confident={n_confident} "
      f"(recognition={n_recognition}/{n_confident})")
print(f"[E2] pocket test: {in_pocket} in-pocket of {sites_total} sites in {len(pocket_fam)} families")
print(f"[E2] phyla: {phylum_ct}")
print("[E2] wrote tables + fig_posel_function.{pdf,png}")
