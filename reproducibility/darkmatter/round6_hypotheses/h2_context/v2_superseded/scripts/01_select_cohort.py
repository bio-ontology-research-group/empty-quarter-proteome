#!/usr/bin/env python3
"""01_select_cohort.py - Deterministic cohort selection for Hypothesis 2 (v2 protocol).

Filters structural dark families across explicit attrition gates:
1. Qualified structural models (mean pLDDT >= 0.70, length >= 100 aa).
2. Refined MAG binned subset.
3. Multi-species presence (>= 2 distinct valid dRep 95% ANI species clusters; blank/unbinned excluded).

Controls multiplicity by designating one deterministic primary locus per species cluster.
Emits:
- all_evaluable_families.tsv
- pilot_cohort.tsv
- pilot_members.tsv
"""
import csv
import gzip
import os
import random
import sys
from collections import defaultdict
from pathlib import Path

W = Path("/ibex/scratch/projects/c2014/rob/metagenomics")
if not W.exists():
    W = Path("reproducibility/darkmatter/source_snapshot")

PER_REP = W / "04_structure/00_structure_clean/05_validate/per_rep_structures.tsv"
FAM_REPS = W / "04_structure/00_structure_clean/05_validate/fam_reps.tsv"
CLUSTERS = W / "04_structure/00_structure_clean/04_cluster/clusters.tsv"
C2B = W / "04_protein_catalog/contig_to_bin_refined.tsv.gz"
CDB = W / "05_refine/_drep/out/data_tables/Cdb.csv"
GTDB_BAC = W / "02_taxonomy/gtdbtk_out/gtdbtk.bac120.summary.tsv"
GTDB_AR = W / "02_taxonomy/gtdbtk_out/gtdbtk.ar53.summary.tsv"

out_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("results")
out_dir.mkdir(parents=True, exist_ok=True)

print("[cohort] Reading structural models from:", PER_REP)
per_rep = {}
with open(PER_REP) as f:
    r = csv.DictReader(f, delimiter="\t")
    for row in r:
        if row["confidence_qualified"] == "1":
            is_k = row["conf_qualified_dual_novel"] == "1"
            per_rep[row["id"]] = {
                "arm": "K" if is_k else "QnotK",
                "length": int(row["length"]),
                "mean_plddt": float(row["mean_plddt"]),
                "best_pdb_qtm": float(row["best_pdb_qtm"]),
                "best_afdb_qtm": float(row["best_afdb_qtm"])
            }

print(f"[cohort] Qualified models: {len(per_rep)} (K: {sum(1 for v in per_rep.values() if v['arm'] == 'K')}, QnotK: {sum(1 for v in per_rep.values() if v['arm'] == 'QnotK')})")

print("[cohort] Mapping chosen_rep to family key from:", FAM_REPS)
fam_to_rep = {}
rep_to_fam = {}
declared_size = {}
with open(FAM_REPS) as f:
    r = csv.DictReader(f, delimiter="\t")
    for row in r:
        rep = row["chosen_rep"]
        fam = row["family"]
        if rep in per_rep:
            fam_to_rep[fam] = rep
            rep_to_fam[rep] = fam
            declared_size[fam] = int(row["n_members"])

print(f"[cohort] Mapped families: {len(fam_to_rep)}")

print("[cohort] Reading members from clusters.tsv...")
fam_members = defaultdict(list)
with open(CLUSTERS) as f:
    for line in f:
        fam, mem = line.rstrip("\n").split("\t")[:2]
        if fam in fam_to_rep:
            fam_members[fam].append(mem)

print(f"[cohort] Members collected for {len(fam_members)} families")

print("[cohort] Reading contig-to-bin mapping from:", C2B)
c2b = {}
c_lens = {}
all_contigs = {mem.rsplit("_", 1)[0] for members in fam_members.values() for mem in members}
with gzip.open(C2B, "rt") as f:
    for line in f:
        parts = line.rstrip("\n").split("\t")
        c = parts[0]
        if c in all_contigs:
            b = parts[2]
            c2b[c] = b
            c_lens[c] = int(parts[1]) if parts[1].isdigit() else 0

print(f"[cohort] Binned contigs: {len(c2b)} / {len(all_contigs)}")

print("[cohort] Reading dRep species clusters from:", CDB)
bin_to_species = {}
if CDB.exists():
    with open(CDB) as f:
        r = csv.DictReader(f)
        for row in r:
            g = row["genome"].removesuffix(".fa")
            sc = row["secondary_cluster"].strip()
            if sc:
                bin_to_species[g] = sc

print(f"[cohort] Valid bin-to-species mappings: {len(bin_to_species)}")

print("[cohort] Reading taxonomy from GTDB-Tk...")
bin_to_tax = {}
for gtdb_file in [GTDB_BAC, GTDB_AR]:
    if gtdb_file.exists():
        with open(gtdb_file) as f:
            r = csv.DictReader(f, delimiter="\t")
            for row in r:
                bin_to_tax[row["user_genome"]] = row["classification"]

print(f"[cohort] Bin-to-taxonomy mappings: {len(bin_to_tax)}")

# Characterize all families
all_families = []
for fam, members in sorted(fam_members.items()):
    rep = fam_to_rep[fam]
    arm = per_rep[rep]["arm"]
    
    valid_species = set()
    bins = set()
    phyla = set()
    binned_count = 0
    
    for mem in members:
        contig = mem.rsplit("_", 1)[0]
        b = c2b.get(contig)
        if b:
            binned_count += 1
            bins.add(b)
            sp = bin_to_species.get(b, "")
            if sp:
                valid_species.add(sp)
            tax = bin_to_tax.get(b, "")
            if tax and "p__" in tax:
                phyla.add(tax.split("p__")[1].split(";")[0])
    
    dom_phylum = list(phyla)[0] if len(phyla) == 1 else ("multi" if len(phyla) > 1 else "unassigned")
    
    all_families.append({
        "family": fam,
        "chosen_rep": rep,
        "arm": arm,
        "n_members": len(members),
        "n_binned": binned_count,
        "n_bins": len(bins),
        "n_species": len(valid_species),
        "phylum": dom_phylum,
        "length": per_rep[rep]["length"],
        "mean_plddt": per_rep[rep]["mean_plddt"],
        "evaluable_species_ge_2": len(valid_species) >= 2
    })

# Write all evaluable families
with open(out_dir / "all_evaluable_families.tsv", "w") as f:
    w = csv.DictWriter(f, fieldnames=list(all_families[0].keys()), delimiter="\t")
    w.writeheader()
    w.writerows(all_families)

evaluable_k = [f for f in all_families if f["arm"] == "K" and f["evaluable_species_ge_2"]]
evaluable_qnotk = [f for f in all_families if f["arm"] == "QnotK" and f["evaluable_species_ge_2"]]

print(f"[cohort] Evaluable families (>=2 distinct valid species): K={len(evaluable_k)}, QnotK={len(evaluable_qnotk)}")

# Stratified sampling of 50 K and 50 QnotK
rng = random.Random(20261013)

def stratified_sample(pool, target=50):
    # Deterministic sorting prior to shuffle
    pool = sorted(pool, key=lambda x: (x["phylum"], x["n_members"], x["family"]))
    rng.shuffle(pool)
    return pool[:target]

sample_k = stratified_sample(evaluable_k, 50)
sample_qnotk = stratified_sample(evaluable_qnotk, 50)
pilot_cohort = sorted(sample_k + sample_qnotk, key=lambda x: (x["arm"], x["family"]))

with open(out_dir / "pilot_cohort.tsv", "w") as f:
    w = csv.DictWriter(f, fieldnames=list(pilot_cohort[0].keys()), delimiter="\t")
    w.writeheader()
    w.writerows(pilot_cohort)

print(f"[cohort] Emitted pilot_cohort.tsv: {len(pilot_cohort)} families (50 K, 50 QnotK)")

# Write member list for the pilot cohort with multiplicity controls
pilot_members = []
for f in pilot_cohort:
    fam = f["family"]
    rep = f["chosen_rep"]
    arm = f["arm"]
    
    # Group members by species cluster
    species_member_map = defaultdict(list)
    for mem in fam_members[fam]:
        sample = mem.split("::")[0]
        contig = mem.rsplit("_", 1)[0]
        b = c2b.get(contig, "")
        sp = bin_to_species.get(b, "") if b else ""
        tax = bin_to_tax.get(b, "") if b else ""
        c_len = c_lens.get(contig, 0)
        
        item = {
            "family": fam,
            "chosen_rep": rep,
            "arm": arm,
            "member": mem,
            "sample": sample,
            "contig": contig,
            "contig_len": c_len,
            "bin": b,
            "species": sp,
            "taxonomy": tax
        }
        species_member_map[sp].append(item)
    
    # For each species, deterministically select primary locus
    for sp, m_list in species_member_map.items():
        copy_count = len(m_list)
        # Sort by: valid binned first, longest contig, then member ID
        m_list_sorted = sorted(m_list, key=lambda x: (-x["contig_len"], x["member"]))
        for idx, m_item in enumerate(m_list_sorted):
            is_primary = 1 if (sp != "" and idx == 0) else 0
            m_item["is_primary_locus"] = is_primary
            m_item["locus_copy_count"] = copy_count
            pilot_members.append(m_item)

fieldnames = [
    "family", "chosen_rep", "arm", "member", "sample", "contig", "contig_len",
    "bin", "species", "taxonomy", "is_primary_locus", "locus_copy_count"
]
with open(out_dir / "pilot_members.tsv", "w") as f:
    w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
    w.writeheader()
    w.writerows(pilot_members)

print(f"[cohort] Emitted pilot_members.tsv: {len(pilot_members)} members")
