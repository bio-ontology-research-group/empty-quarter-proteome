#!/usr/bin/env python3
"""Split the harmonised protein FASTA into desert (EQ) queries and external targets,
and write the genome -> target-group map used to reduce the search output."""
import argparse
from pathlib import Path
import pandas as pd

EXCLUDED = {"External_Comparators_CrossResource_Duplicate_54",
            "SPIRE_ProGenomes3_SpecI_Excluded_148"}
SHORT = {"EQ_Thermoplasmatota_237": "eq", "External_Comparators_Deduplicated_247": "external247",
         "Maeke_2025_Penumbrarchaeia_35": "maeke35", "SMAG_Thermoplasmatota_354": "smag354",
         "SPIRE_Thermoplasmatota_1084_DeNovo": "spire1084"}
ap = argparse.ArgumentParser()
ap.add_argument("--root", type=Path, required=True)
ap.add_argument("--out", type=Path, required=True)
a = ap.parse_args()
gm = pd.read_csv(a.root / "metadata/harmonized_genome_manifest.tsv", sep="\t", dtype=str)
gm = gm[~gm.cohort.isin(EXCLUDED)].copy()
assert len(gm) == 1957 and gm.genome_id.is_unique
gm["cohort_short"] = gm.cohort.map(SHORT)
assert not gm.genome_id.str.contains("::").any()
gm["group"] = gm.cohort_short + "|" + gm.tax_order
gm[["genome_id", "group"]].to_csv(a.out / "genome_groups.tsv", sep="\t", index=False, header=False)
coh = dict(zip(gm.genome_id, gm.cohort_short))
n = {"eq": 0, "ext": 0}
with open(a.root / "genecall/harmonized_proteins.faa") as fh, \
        open(a.out / "eq_proteins.faa", "w") as fe, open(a.out / "external_proteins.faa", "w") as fx:
    dest = None
    for line in fh:
        if line[0] == ">":
            g = line[1:].split()[0].split("::")[0]
            c = coh[g]  # KeyError if the FASTA holds a genome outside the 1,957
            if c == "eq":
                dest = fe; n["eq"] += 1
            else:
                dest = fx; n["ext"] += 1
        dest.write(line)
print("proteins", n, sum(n.values()))
assert sum(n.values()) == 3617011
