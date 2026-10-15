#!/usr/bin/env python3
"""Step 3 -- pathogen census on the CURRENT 5,438-representative catalogue.

Reuses, verbatim, the genus and species panels defined by the PI's
06_resistome_safety/code/pathogen_census.py (imported at runtime, not copied) and
re-runs the same two-resolution count on the refined R232 taxonomy of the 5,438
dRep representatives. The historical run (PATHOGENS.md) was on the 15,469 MAG set:
7 MAGs / 0.05% named pathogen species. On a species basis the statistic RESTATES
rather than reproduces -- the natural unit is "N named pathogen species of 5,438
species" -- so both framings are written out.

Out (results/):
  pathogen_species_5438.tsv   named pathogen species detected among the reps
  pathogen_genus_5438.tsv     genus-level environmental-relative census
  pathogen_census_summary.tsv headline scalars

Run: python3 scripts/amr_3_pathogens.py   (taxonomy lookup only; seconds)
"""
import csv, importlib.util, os
from collections import Counter, defaultdict

M = "/ibex/scratch/projects/c2014/rob/metagenomics"
B = f"{M}/06_metabolism/resistome_clean"
RES = f"{B}/results"
GTDB = [f"{M}/05_refine/_gtdbtk_r232/gtdbtk.bac120.summary.tsv",
        f"{M}/05_refine/_gtdbtk_r232/gtdbtk.ar53.summary.tsv"]

# ---- import the PI's panels verbatim -----------------------------------------
spec = importlib.util.spec_from_file_location(
    "pathogen_census", f"{M}/06_resistome_safety/code/pathogen_census.py")
pc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pc)
HA, PL, SP = pc.HUMAN_ANIMAL_GENERA, pc.PLANT_GENERA, pc.PATHOGEN_SPECIES
print(f"[pathogens] panels: human/animal genera={len(HA)}, plant genera={len(PL)}, "
      f"named species={len(SP)}")

def rank(lin, pre):
    for t in lin.split(";"):
        t = t.strip()
        if t.startswith(pre):
            return t[3:]
    return ""

total = 0
genus_hits = defaultdict(lambda: defaultdict(Counter))
species_hits = Counter()
cat_mags = defaultdict(set)
for p in GTDB:
    for r in csv.DictReader(open(p), delimiter="\t"):
        total += 1
        mag, lin = r["user_genome"], r["classification"]
        g, s = rank(lin, "g__"), rank(lin, "s__")
        cat = "human_animal" if g in HA else ("plant" if g in PL else None)
        if cat:
            genus_hits[cat][g][s or "(unnamed sp.)"] += 1
            cat_mags[cat].add(mag)
        if s in SP:
            species_hits[s] += 1
assert total == 5438

with open(f"{RES}/pathogen_species_5438.tsv", "w") as o:
    o.write("named_pathogen_species\tn_species_reps\n")
    for s, n in species_hits.most_common():
        o.write(f"{s}\t{n}\n")
with open(f"{RES}/pathogen_genus_5438.tsv", "w") as o:
    o.write("category\tgenus\tn_species_reps\tnamed_pathogen_species_present\n")
    for cat in ("human_animal", "plant"):
        for g, spc in sorted(genus_hits[cat].items(), key=lambda kv: -sum(kv[1].values())):
            patho = sorted(s for s in spc if s in SP)
            o.write(f"{cat}\t{g}\t{sum(spc.values())}\t{';'.join(patho) or '-'}\n")

n_ha, n_pl = len(cat_mags["human_animal"]), len(cat_mags["plant"])
n_named_reps = sum(species_hits.values())
with open(f"{RES}/pathogen_census_summary.tsv", "w") as o:
    o.write("metric\tvalue\n")
    for k, v in [("n_species", total),
                 ("n_named_pathogen_species", len(species_hits)),
                 ("n_reps_named_pathogen", n_named_reps),
                 ("pct_reps_named_pathogen", f"{100*n_named_reps/total:.3f}"),
                 ("n_reps_pathogen_containing_genus_human_animal", n_ha),
                 ("pct_human_animal_genus", f"{100*n_ha/total:.2f}"),
                 ("n_reps_pathogen_containing_genus_plant", n_pl),
                 ("pct_plant_genus", f"{100*n_pl/total:.2f}")]:
        o.write(f"{k}\t{v}\n")
print(f"[pathogens] named pathogen species among 5,438 reps: {dict(species_hits)}")
print(f"[pathogens] genus-level environmental relatives: human/animal {n_ha} "
      f"({100*n_ha/total:.1f}%), plant {n_pl} ({100*n_pl/total:.1f}%)")
print("[pathogens] done")
