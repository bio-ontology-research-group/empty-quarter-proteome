#!/usr/bin/env python3
"""Reproduce the registered pathogen census on the 5,438 r232 representatives, and inventory
the sequenced samples that are not in the canonical 551.

The census script in W/06_resistome_safety/code/pathogen_census.py matches GTDB s__ strings by
exact equality against an unsuffixed NCBI-style list. GTDB appends genus suffixes (_A, _B), so any
pathogen whose GTDB genus is suffixed can never match. This recomputes the census as written, then
recounts with a suffix-blind matcher to show exactly which species the strict version misses.
"""
import csv
import os
import glob
import pandas as pd

W = "/ibex/scratch/projects/c2014/rob/metagenomics"
D = "/ibex/scratch/projects/c2014/EmptyQuarter_Data"
K = "/ibex/scratch/projects/c2014/alelopezv/eqmgc/validation/kraken2/r1_list.txt"
RES = "/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/contamination_clean/results"

PATHOGEN_SPECIES = {  # verbatim from 06_resistome_safety/code/pathogen_census.py
    "Mycobacterium tuberculosis", "Mycobacterium leprae", "Pseudomonas aeruginosa",
    "Staphylococcus aureus", "Klebsiella pneumoniae", "Escherichia coli", "Salmonella enterica",
    "Streptococcus pneumoniae", "Streptococcus pyogenes", "Enterococcus faecalis",
    "Enterococcus faecium", "Acinetobacter baumannii", "Bacillus anthracis", "Bacillus cereus",
    "Listeria monocytogenes", "Clostridioides difficile", "Clostridium botulinum",
    "Clostridium perfringens", "Vibrio cholerae", "Yersinia pestis", "Neisseria meningitidis",
    "Neisseria gonorrhoeae", "Legionella pneumophila", "Brucella melitensis", "Brucella abortus",
    "Bordetella pertussis", "Francisella tularensis", "Burkholderia pseudomallei",
    "Burkholderia mallei", "Corynebacterium diphtheriae", "Helicobacter pylori",
    "Ralstonia solanacearum", "Xanthomonas campestris", "Xanthomonas oryzae",
    "Pectobacterium carotovorum", "Erwinia amylovora", "Clavibacter michiganensis",
    "Pseudomonas syringae", "Xylella fastidiosa", "Agrobacterium tumefaciens",
}


def base(name):
    """Strip GTDB genus and species suffixes: 'Bacillus_A anthracis' -> 'Bacillus anthracis'."""
    if " " not in name:
        return name
    g, s = name.split(" ", 1)
    return f"{g.split('_')[0]} {s.split('_')[0]}"


rows = []
for f in ("gtdbtk.bac120.summary.tsv", "gtdbtk.ar53.summary.tsv"):
    for r in csv.DictReader(open(os.path.join(W, "05_refine/_gtdbtk_r232", f)), delimiter="\t"):
        sp = next((t[3:] for t in r["classification"].split(";") if t.startswith("s__")), "")
        rows.append({"genome": r["user_genome"], "species": sp,
                     "strict": sp in PATHOGEN_SPECIES,
                     "suffix_blind": base(sp) in PATHOGEN_SPECIES})
g = pd.DataFrame(rows)
n = len(g)
print(f"representatives scanned: {n}")
print(f"STRICT (as registered): {g.strict.sum()} MAGs = {100*g.strict.sum()/n:.2f}% "
      f"-> {sorted(g.loc[g.strict, 'species'].unique())}")
print(f"SUFFIX-BLIND:           {g.suffix_blind.sum()} MAGs = {100*g.suffix_blind.sum()/n:.2f}% "
      f"-> {sorted(g.loc[g.suffix_blind, 'species'].unique())}")
missed = g[g.suffix_blind & ~g.strict]
print("\nMissed purely by GTDB suffix naming:")
print(missed[["genome", "species"]].to_string(index=False))
g[g.suffix_blind | g.strict].to_csv(os.path.join(RES, "pathogen_census_5438.tsv"),
                                    sep="\t", index=False)

# ---------------------------------------------------------------- controls inventory
canon = set(pd.read_csv(os.path.join(W, "metadata/sample_metadata.tsv"), sep="\t")["sample"])
kr = {os.path.basename(x).replace("_R1.fastq.gz", "").strip()
      for x in open(K) if x.strip()}
reads_dir = os.path.join(D, "soil/wgs_assemblies/processed_reads")
on_disk = {os.path.basename(p).split("_R1.fastq.gz")[0]
           for p in glob.glob(os.path.join(reads_dir, "*_R1.fastq.gz"))}
inv = []
for s in sorted(on_disk - canon):
    inv.append({
        "sample_id": s,
        "in_kraken_554": s in kr,
        "reads_R1_bytes": os.path.getsize(os.path.join(reads_dir, f"{s}_R1.fastq.gz")),
        "has_assembly": os.path.exists(os.path.join(D, "soil/wgs_assemblies/assemblies", s,
                                                    f"{s}_contigs.fa")),
        "in_abundance_matrix": False,
    })
inv = pd.DataFrame(inv)
inv.to_csv(os.path.join(RES, "noncanonical_sequenced_samples.tsv"), sep="\t", index=False)
print(f"\nsequenced soil libraries not in the canonical 551: {len(inv)}")
print(f"of these, in the 554-sample kraken list: {sorted(inv.loc[inv.in_kraken_554,'sample_id'])}")
