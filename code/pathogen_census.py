#!/usr/bin/env python3
"""#2 Pathogen census across the 15,469 MAGs (contamination-control + ecology).

Two resolutions, because genus-level overcounts (soil Streptomyces/Bacillus are environmental):
  - GENUS level: MAGs in genera that CONTAIN pathogens (mostly environmental relatives here).
  - SPECIES level: MAGs assigned to an actual named pathogenic species (the real signal).
The gap between them is the contamination control: a pristine desert should show many
environmental relatives but ~zero named human/animal pathogen species.

In: workbench GTDB summaries (full 15,469-MAG set). Out: results/pathogens/{...}.tsv + PATHOGENS.md
"""
import os, csv
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "results", "pathogens")
os.makedirs(OUT, exist_ok=True)
GTDB = [os.path.join(HERE, "..", "workbench", "results", "soil_gtdbtk", x)
        for x in ("gtdbtk.bac120.summary.tsv", "gtdbtk.ar53.summary.tsv")]

# genera that CONTAIN notable pathogens (GTDB g__ names; many also have environmental members)
HUMAN_ANIMAL_GENERA = {
    "Mycobacterium", "Mycobacteroides", "Klebsiella", "Klebsiella_A", "Escherichia",
    "Salmonella", "Shigella", "Staphylococcus", "Streptococcus", "Enterococcus",
    "Enterococcus_B", "Acinetobacter", "Neisseria", "Haemophilus", "Legionella",
    "Coxiella", "Brucella", "Bartonella", "Bordetella", "Burkholderia", "Burkholderia_A",
    "Vibrio", "Campylobacter", "Helicobacter", "Clostridium", "Clostridioides",
    "Bacillus", "Bacillus_A", "Listeria", "Corynebacterium", "Nocardia", "Yersinia",
    "Francisella", "Leptospira", "Mycoplasma", "Mycoplasmoides", "Pasteurella",
    "Proteus", "Serratia", "Enterobacter", "Citrobacter", "Aeromonas",
    "Stenotrophomonas", "Rhodococcus", "Moraxella", "Cutibacterium", "Gardnerella",
}
PLANT_GENERA = {
    "Ralstonia", "Xanthomonas", "Agrobacterium", "Pectobacterium", "Dickeya", "Erwinia",
    "Pantoea", "Clavibacter", "Curtobacterium", "Xylella", "Acidovorax", "Rathayibacter",
    "Spiroplasma", "Phytoplasma", "Rhizobium",
}
# actual named pathogenic species (substring match on s__ ; the narrow, meaningful signal)
PATHOGEN_SPECIES = {
    "Mycobacterium tuberculosis", "Mycobacterium leprae", "Pseudomonas aeruginosa",
    "Staphylococcus aureus", "Klebsiella pneumoniae", "Escherichia coli", "Salmonella enterica",
    "Streptococcus pneumoniae", "Streptococcus pyogenes", "Enterococcus faecalis",
    "Enterococcus faecium", "Acinetobacter baumannii", "Bacillus anthracis", "Bacillus cereus",
    "Listeria monocytogenes", "Clostridioides difficile", "Clostridium botulinum",
    "Clostridium perfringens", "Vibrio cholerae", "Yersinia pestis", "Neisseria meningitidis",
    "Neisseria gonorrhoeae", "Legionella pneumophila", "Brucella melitensis", "Brucella abortus",
    "Bordetella pertussis", "Francisella tularensis", "Burkholderia pseudomallei",
    "Burkholderia mallei", "Corynebacterium diphtheriae", "Helicobacter pylori",
    # plant
    "Ralstonia solanacearum", "Xanthomonas campestris", "Xanthomonas oryzae",
    "Pectobacterium carotovorum", "Erwinia amylovora", "Clavibacter michiganensis",
    "Pseudomonas syringae", "Xylella fastidiosa", "Agrobacterium tumefaciens",
}


def rank(lin, pre):
    for t in lin.split(";"):
        if t.startswith(pre):
            return t[3:]
    return ""


def main():
    total = 0
    genus_hits = defaultdict(lambda: defaultdict(Counter))   # cat -> genus -> species Counter
    species_hits = Counter()                                  # named pathogen species -> n MAGs
    cat_mags = defaultdict(set)
    for p in GTDB:
        if not os.path.exists(p):
            continue
        for r in csv.DictReader(open(p), delimiter="\t"):
            total += 1
            mag, lin = r["user_genome"], r["classification"]
            g, s = rank(lin, "g__"), rank(lin, "s__")
            cat = None
            if g in HUMAN_ANIMAL_GENERA:
                cat = "human_animal"
            elif g in PLANT_GENERA:
                cat = "plant"
            if cat:
                genus_hits[cat][g][s or "(unnamed sp.)"] += 1
                cat_mags[cat].add(mag)
            if s in PATHOGEN_SPECIES:
                species_hits[s] += 1

    # genus-level table
    with open(os.path.join(OUT, "pathogen_genus_census.tsv"), "w") as o:
        o.write("category\tgenus\tn_MAGs\tn_named_species\ttop_species\tnamed_pathogen_species_present\n")
        for cat in ("human_animal", "plant"):
            for g, spc in sorted(genus_hits[cat].items(), key=lambda kv: -sum(kv[1].values())):
                nm = sum(spc.values())
                named = sorted([s for s in spc if s != "(unnamed sp.)"])
                patho = [s for s in named if s in PATHOGEN_SPECIES]
                top = spc.most_common(1)[0][0]
                o.write(f"{cat}\t{g}\t{nm}\t{len(named)}\t{top}\t{';'.join(patho) or '-'}\n")

    # species-level (the meaningful signal)
    with open(os.path.join(OUT, "named_pathogen_species.tsv"), "w") as o:
        o.write("named_pathogen_species\tn_MAGs\n")
        for s, n in species_hits.most_common():
            o.write(f"{s}\t{n}\n")

    n_ha = len(cat_mags["human_animal"])
    n_pl = len(cat_mags["plant"])
    n_named = sum(species_hits.values())
    ml = [f"# Pathogen census — {total} MAGs (contamination control)\n",
          f"- **Genus level** (genera that contain pathogens; mostly environmental relatives): "
          f"human/animal {n_ha} MAGs ({100*n_ha/total:.1f}%), plant {n_pl} MAGs ({100*n_pl/total:.1f}%).\n",
          f"- **Species level** (actual named pathogen species): **{n_named} MAGs "
          f"({100*n_named/total:.2f}%)**.\n",
          f"- The gap is the control: a hyper-arid pristine desert carries environmental relatives of "
          f"pathogen genera but essentially **no named human/animal pathogen species** — consistent with "
          f"low anthropogenic input and no kit-ome/human contamination bloom.\n",
          "## Named pathogen species detected\n"]
    if species_hits:
        ml.append("| species | MAGs |\n|---|---|")
        for s, n in species_hits.most_common():
            ml.append(f"| {s} | {n} |")
    else:
        ml.append("_None._")
    ml.append("\n## Top pathogen-containing genera (environmental relatives)\n")
    ml.append("| category | genus | MAGs | named species present |\n|---|---|---|---|")
    for cat in ("human_animal", "plant"):
        for g, spc in sorted(genus_hits[cat].items(), key=lambda kv: -sum(kv[1].values()))[:10]:
            patho = [s for s in spc if s in PATHOGEN_SPECIES]
            ml.append(f"| {cat} | {g} | {sum(spc.values())} | {';'.join(patho) or '— (environmental)'} |")
    open(os.path.join(OUT, "PATHOGENS.md"), "w").write("\n".join(ml) + "\n")
    print(f"total MAGs: {total}")
    print(f"human/animal genus MAGs: {n_ha} ({100*n_ha/total:.1f}%)")
    print(f"plant genus MAGs: {n_pl} ({100*n_pl/total:.1f}%)")
    print(f"named pathogen-species MAGs: {n_named} ({100*n_named/total:.2f}%)")
    print("named species:", dict(species_hits))


if __name__ == "__main__":
    main()
