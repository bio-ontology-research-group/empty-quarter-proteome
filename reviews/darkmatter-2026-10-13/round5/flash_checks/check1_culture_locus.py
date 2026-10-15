#!/usr/bin/env python3
"""Check 1: Exploratory analysis of culture vs MAG gene calling discrepancies.
Investigates whether DNA-supported culture no-CDS discrepancies persist in
complete, long (>=100 aa), and recurrent family calls, and across individual pairs.
"""
import json
from pathlib import Path
import pandas as pd
from scipy.stats import fisher_exact

BASE = Path(__file__).resolve().parents[4]
DM = BASE / "reproducibility" / "darkmatter"

cigar = pd.read_csv(DM / "culture_nucleotide_cigar" / "per_gene.tsv", sep="\t")
locus = pd.read_csv(DM / "culture_locus_all" / "per_gene.tsv", sep="\t")
fam = pd.read_csv(DM / "culture_family_membership" / "per_gene.tsv", sep="\t")

# Merge
m = cigar.merge(locus[["protein", "classification", "overlapping_cds_hits"]], on="protein", how="left")
m = m.merge(fam[["protein", "cluster90", "admitted_family"]], on="protein", how="left")

el = m[m.eligible == 1].copy()
dna = el[el.dna_supported == 1].copy()

def classify_stats(df):
    vc = df["classification"].value_counts().to_dict()
    total = len(df)
    return {
        "total": total,
        "no_major_cds": vc.get("no_major_CDS_overlap", 0),
        "no_major_cds_pct": vc.get("no_major_CDS_overlap", 0) / total * 100 if total else 0.0,
        "same_orientation_start_phase": vc.get("same_orientation_start_phase", 0),
        "same_orientation_start_phase_pct": vc.get("same_orientation_start_phase", 0) / total * 100 if total else 0.0,
        "other_frame_or_orientation": vc.get("other_frame_or_orientation", 0),
        "other_frame_or_orientation_pct": vc.get("other_frame_or_orientation", 0) / total * 100 if total else 0.0,
    }

# Subsets
dark_dna = dna[dna.dark == 1]
ann_dna = dna[dna.dark == 0]

dark_comp = dark_dna[(dark_dna.partial == "00") | (dark_dna.partial == 0)]
ann_comp = ann_dna[(ann_dna.partial == "00") | (ann_dna.partial == 0)]

dark_comp_long = dark_dna[((dark_dna.partial == "00") | (dark_dna.partial == 0)) & (dark_dna.nt >= 300)]
ann_comp_long = ann_dna[((ann_dna.partial == "00") | (ann_dna.partial == 0)) & (ann_dna.nt >= 300)]

dark_rec = dark_dna[dark_dna.admitted_family.notna() & (dark_dna.admitted_family != "")]

# Fisher exact test for complete long:
# Contingency table: [[dark_no_cds, dark_cds_present], [ann_no_cds, ann_cds_present]]
d_no = (dark_comp_long.classification == "no_major_CDS_overlap").sum()
d_yes = len(dark_comp_long) - d_no
a_no = (ann_comp_long.classification == "no_major_CDS_overlap").sum()
a_yes = len(ann_comp_long) - a_no

table = [[int(d_no), int(d_yes)], [int(a_no), int(a_yes)]]
odds_ratio, p_val = fisher_exact(table)

# Per pair
per_pair = dark_dna.groupby("rep").agg(
    n_dna_dark=("protein", "count"),
    n_no_cds=("classification", lambda s: (s == "no_major_CDS_overlap").sum()),
    n_same_phase=("classification", lambda s: (s == "same_orientation_start_phase").sum()),
    n_other=("classification", lambda s: (s == "other_frame_or_orientation").sum())
).reset_index()
per_pair["frac_no_cds"] = per_pair["n_no_cds"] / per_pair["n_dna_dark"]

results = {
    "all_dna_supported_dark": classify_stats(dark_dna),
    "complete_dna_supported_dark": classify_stats(dark_comp),
    "complete_long_dna_supported_dark": classify_stats(dark_comp_long),
    "recurrent_family_dna_supported_dark": classify_stats(dark_rec),
    "all_dna_supported_annotated": classify_stats(ann_dna),
    "complete_long_dna_supported_annotated": classify_stats(ann_comp_long),
    "complete_long_contingency_table": table,
    "complete_long_fisher_odds_ratio": float(odds_ratio),
    "complete_long_fisher_pvalue": float(p_val),
    "per_pair_summary": {
        "n_pairs": len(per_pair),
        "pairs_with_no_cds": int((per_pair.n_no_cds > 0).sum()),
        "median_frac_no_cds": float(per_pair.frac_no_cds.median()),
        "min_frac_no_cds": float(per_pair.frac_no_cds.min()),
        "max_frac_no_cds": float(per_pair.frac_no_cds.max()),
    }
}

out_dir = Path(__file__).resolve().parent
with open(out_dir / "check1_culture_locus.json", "w") as f:
    json.dump(results, f, indent=2)

per_pair.to_csv(out_dir / "check1_culture_locus_per_pair.tsv", sep="\t", index=False)

print(json.dumps(results, indent=2))
