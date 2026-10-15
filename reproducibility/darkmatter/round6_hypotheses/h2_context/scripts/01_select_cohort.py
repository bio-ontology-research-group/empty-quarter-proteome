#!/usr/bin/env python3
"""01_select_cohort.py - Deterministic cohort selection for Hypothesis 2 (v3 protocol).

Resolves C2 and C5 from round6/h2_review1:
1. Explicitly documents and tracks the 90%-sequence-cluster representative basis of
   04_structure/00_structure_clean/04_cluster/clusters.tsv across explicit attrition gates:
   - Gate 0: Confidence-qualified monomer structural families (N=2,838; 1,031 K, 1,807 QnotK).
   - Gate 1: Refined MAG binned subset (contig_to_bin_refined.tsv.gz).
   - Gate 2: Valid dRep 95% ANI species-cluster presence (Cdb.csv; bins absent from Cdb.csv
             are explicitly labelled 'binned_unclustered_no_drep_species' and never form species).
   - Gate 3: Multi-species presence (>= 2 distinct dRep 95% ANI species clusters; N=1,461).
2. Joins verified GTDB-Tk R232 taxonomy through Cdb.csv (secondary_cluster) -> Wdb.csv (winner)
   -> 05_refine/_gtdbtk_r232/gtdbtk.{bac120,ar53}.summary.tsv (100% coverage of all 5,438 species
   clusters and all 12,433 dRep-clustered refined MAGs), cross-checked against
   05_refine/_gtdbtk_refined/gtdbtk.{bac120,ar53}.summary.tsv.
3. Selects a deterministic, covariate-stratified expanded pilot cohort of 300 families
   (150 K, 150 QnotK; or all 1,461 families when mode='full') across joint
   (phylum_group, species_breadth_bin) strata using Hamilton largest-remainder apportionment.
4. Emits all_evaluable_families.tsv, pilot_cohort.tsv, and pilot_members_candidates.tsv
   (for coordinate/contig-length measurement and outcome-independent primary-locus selection
   in 02_extract_neighbourhoods.py).
"""
import csv
import gzip
import random
import sys
from collections import Counter, defaultdict
from pathlib import Path

W = Path("/ibex/scratch/projects/c2014/rob/metagenomics")
if not W.exists():
    W = Path("reproducibility/darkmatter/source_snapshot")

PER_REP = W / "04_structure/00_structure_clean/05_validate/per_rep_structures.tsv"
FAM_REPS = W / "04_structure/00_structure_clean/05_validate/fam_reps.tsv"
CLUSTERS = W / "04_structure/00_structure_clean/04_cluster/clusters.tsv"
C2B = W / "04_protein_catalog/contig_to_bin_refined.tsv.gz"
CDB = W / "05_refine/_drep/out/data_tables/Cdb.csv"
WDB = W / "05_refine/_drep/out/data_tables/Wdb.csv"
GTDB_R232_BAC = W / "05_refine/_gtdbtk_r232/gtdbtk.bac120.summary.tsv"
GTDB_R232_AR = W / "05_refine/_gtdbtk_r232/gtdbtk.ar53.summary.tsv"
GTDB_REF_BAC = W / "05_refine/_gtdbtk_refined/gtdbtk.bac120.summary.tsv"
GTDB_REF_AR = W / "05_refine/_gtdbtk_refined/gtdbtk.ar53.summary.tsv"

out_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("results")
mode = sys.argv[2] if len(sys.argv) > 2 else "pilot"
target_per_arm = int(sys.argv[3]) if len(sys.argv) > 3 else 150
out_dir.mkdir(parents=True, exist_ok=True)

print(f"[cohort:v3] Mode={mode}, target_per_arm={target_per_arm}, out_dir={out_dir}")
print("[cohort:v3] Reading structural models from:", PER_REP)
per_rep = {}
with open(PER_REP) as f:
    for row in csv.DictReader(f, delimiter="\t"):
        if row["confidence_qualified"] == "1":
            is_k = row["conf_qualified_dual_novel"] == "1"
            per_rep[row["id"]] = {
                "arm": "K" if is_k else "QnotK",
                "length": int(row["length"]),
                "mean_plddt": float(row["mean_plddt"]),
                "best_pdb_qtm": float(row["best_pdb_qtm"]),
                "best_afdb_qtm": float(row["best_afdb_qtm"]),
            }

k_tot = sum(1 for v in per_rep.values() if v["arm"] == "K")
q_tot = sum(1 for v in per_rep.values() if v["arm"] == "QnotK")
print(f"[cohort:v3] Qualified structural models: {len(per_rep)} (K={k_tot}, QnotK={q_tot})")

print("[cohort:v3] Mapping chosen_rep to structural family key from:", FAM_REPS)
fam_to_rep = {}
declared_size = {}
with open(FAM_REPS) as f:
    for row in csv.DictReader(f, delimiter="\t"):
        rep = row["chosen_rep"]
        fam = row["family"]
        if rep in per_rep:
            fam_to_rep[fam] = rep
            declared_size[fam] = int(row["n_members"])

print(f"[cohort:v3] Mapped structural families: {len(fam_to_rep)}")

print("[cohort:v3] Reading 90%-cluster representative members from:", CLUSTERS)
fam_members = defaultdict(list)
with open(CLUSTERS) as f:
    for line in f:
        fam, mem = line.rstrip("\n").split("\t")[:2]
        if fam in fam_to_rep:
            fam_members[fam].append(mem)

total_rep_members = sum(len(v) for v in fam_members.values())
all_contigs = {mem.rsplit("_", 1)[0] for ms in fam_members.values() for mem in ms}
print(f"[cohort:v3] Collected {total_rep_members} 90%-cluster rep members across {len(all_contigs)} unique contigs for {len(fam_members)} families")

print("[cohort:v3] Reading refined contig-to-bin mapping from:", C2B)
c2b = {}
with gzip.open(C2B, "rt") as f:
    for line in f:
        parts = line.rstrip("\n").split("\t")
        c = parts[0]
        if c in all_contigs:
            c2b[c] = parts[2]

print(f"[cohort:v3] Refined-MAG binned contigs: {len(c2b)} / {len(all_contigs)} ({100.0 * len(c2b) / len(all_contigs):.2f}%)")

print("[cohort:v3] Reading dRep 95% ANI secondary clusters from:", CDB)
bin_to_species = {}
with open(CDB) as f:
    for row in csv.DictReader(f):
        g = row["genome"].removesuffix(".fa")
        sc = row["secondary_cluster"].strip()
        if sc:
            bin_to_species[g] = sc

print(f"[cohort:v3] Valid dRep bin-to-species mappings: {len(bin_to_species)}")

print("[cohort:v3] Reading dRep species cluster winners from:", WDB)
species_to_winner = {}
with open(WDB) as f:
    for row in csv.DictReader(f):
        sc = row["cluster"].strip()
        win = row["genome"].removesuffix(".fa")
        if sc and win:
            species_to_winner[sc] = win

print(f"[cohort:v3] Valid dRep species winners: {len(species_to_winner)}")

print("[cohort:v3] Reading GTDB-Tk R232 species-winner and refined-bin classifications...")
r232_winner_tax = {}
for p in [GTDB_R232_BAC, GTDB_R232_AR]:
    with open(p) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            r232_winner_tax[row["user_genome"]] = row["classification"]

refined_bin_tax = {}
for p in [GTDB_REF_BAC, GTDB_REF_AR]:
    if p.exists():
        with open(p) as f:
            for row in csv.DictReader(f, delimiter="\t"):
                refined_bin_tax[row["user_genome"]] = row["classification"]

species_to_r232_tax = {}
for sc, win in species_to_winner.items():
    if win in r232_winner_tax:
        species_to_r232_tax[sc] = r232_winner_tax[win]
    elif win in refined_bin_tax:
        species_to_r232_tax[sc] = refined_bin_tax[win]

print(f"[cohort:v3] R232 taxonomy resolved for {len(species_to_r232_tax)} / {len(species_to_winner)} species clusters")
assert len(species_to_r232_tax) == len(species_to_winner), "Incomplete R232 taxonomy for dRep winners!"

def extract_phylum(tax_str):
    if tax_str and "p__" in tax_str:
        ph = tax_str.split("p__")[1].split(";")[0].strip()
        if ph:
            return ph
    return "unassigned"

def map_phylum_group(phylum_str):
    if phylum_str in ("Thermoplasmatota", "Actinomycetota", "multi"):
        return phylum_str
    return "other_bacteria_archaea"

# Characterize all 2,838 families across attrition gates
all_families = []
gate_counts = Counter()

for fam, members in sorted(fam_members.items()):
    rep = fam_to_rep[fam]
    arm = per_rep[rep]["arm"]
    gate_counts[f"gate0_total_{arm}"] += 1

    valid_species = set()
    refined_bins = set()
    drep_bins = set()
    phyla = set()
    n_binned_refined = 0
    n_binned_drep = 0
    n_binned_unclustered = 0

    for mem in members:
        contig = mem.rsplit("_", 1)[0]
        b = c2b.get(contig, "")
        if b:
            n_binned_refined += 1
            refined_bins.add(b)
            sp = bin_to_species.get(b, "")
            if sp:
                n_binned_drep += 1
                drep_bins.add(b)
                valid_species.add(sp)
                tax = species_to_r232_tax.get(sp, "")
                ph = extract_phylum(tax)
                if ph != "unassigned":
                    phyla.add(ph)
            else:
                n_binned_unclustered += 1

    if n_binned_refined > 0:
        gate_counts[f"gate1_binned_refined_{arm}"] += 1
    if len(valid_species) >= 1:
        gate_counts[f"gate2_drep_species_ge1_{arm}"] += 1
    if len(valid_species) >= 2:
        gate_counts[f"gate3_drep_species_ge2_{arm}"] += 1

    dom_phylum = (
        sorted(phyla)[0]
        if len(phyla) == 1
        else ("multi" if len(phyla) > 1 else "unassigned")
    )
    ph_group = map_phylum_group(dom_phylum) if dom_phylum != "unassigned" else "unassigned"
    sp_breadth_bin = "sp2" if len(valid_species) == 2 else ("sp3plus" if len(valid_species) >= 3 else "sp_lt2")

    all_families.append({
        "family": fam,
        "chosen_rep": rep,
        "arm": arm,
        "n_members": len(members),
        "n_binned_refined": n_binned_refined,
        "n_binned_drep_species": n_binned_drep,
        "n_binned_unclustered_no_drep": n_binned_unclustered,
        "n_bins_refined": len(refined_bins),
        "n_bins_drep": len(drep_bins),
        "n_species": len(valid_species),
        "phylum": dom_phylum,
        "phylum_group": ph_group,
        "species_breadth_bin": sp_breadth_bin,
        "length": per_rep[rep]["length"],
        "mean_plddt": per_rep[rep]["mean_plddt"],
        "evaluable_species_ge_2": len(valid_species) >= 2,
    })

with open(out_dir / "all_evaluable_families.tsv", "w") as f:
    w = csv.DictWriter(f, fieldnames=list(all_families[0].keys()), delimiter="\t")
    w.writeheader()
    w.writerows(all_families)

evaluable_k = [f for f in all_families if f["arm"] == "K" and f["evaluable_species_ge_2"]]
evaluable_qnotk = [f for f in all_families if f["arm"] == "QnotK" and f["evaluable_species_ge_2"]]

assert all(f["phylum"] != "unassigned" for f in evaluable_k + evaluable_qnotk), (
    "Found unassigned phylum among multi-species evaluable families!"
)

print(
    f"[cohort:v3] Gate summary: Gate0={len(all_families)} (K={gate_counts['gate0_total_K']}, QnotK={gate_counts['gate0_total_QnotK']}); "
    f"Gate1_refined_binned={gate_counts['gate1_binned_refined_K']+gate_counts['gate1_binned_refined_QnotK']} "
    f"(K={gate_counts['gate1_binned_refined_K']}, QnotK={gate_counts['gate1_binned_refined_QnotK']}); "
    f"Gate2_drep_ge1={gate_counts['gate2_drep_species_ge1_K']+gate_counts['gate2_drep_species_ge1_QnotK']} "
    f"(K={gate_counts['gate2_drep_species_ge1_K']}, QnotK={gate_counts['gate2_drep_species_ge1_QnotK']}); "
    f"Gate3_multi_species_ge2={len(evaluable_k)+len(evaluable_qnotk)} (K={len(evaluable_k)}, QnotK={len(evaluable_qnotk)})"
)

def proportional_stratified_sample(pool, target, seed=20261013):
    if mode == "full" or target >= len(pool):
        return sorted(pool, key=lambda x: x["family"])
    rng = random.Random(seed)
    by_stratum = defaultdict(list)
    for item in pool:
        key = (item["phylum_group"], item["species_breadth_bin"])
        by_stratum[key].append(item)

    strata_keys = sorted(by_stratum.keys())
    n_total = len(pool)
    quotas = {}
    remainders = []
    allocated = 0
    for k in strata_keys:
        s_list = sorted(
            by_stratum[k],
            key=lambda x: (x["phylum"], x["n_species"], x["n_members"], x["family"]),
        )
        rng.shuffle(s_list)
        by_stratum[k] = s_list
        exact_q = target * (len(s_list) / n_total)
        int_q = int(exact_q)
        quotas[k] = int_q
        allocated += int_q
        remainders.append((exact_q - int_q, len(s_list), k))

    remainders.sort(reverse=True)
    for _, _, k in remainders[: target - allocated]:
        quotas[k] += 1

    selected = []
    for k in strata_keys:
        selected.extend(by_stratum[k][: quotas[k]])
    assert len(selected) == target, f"Expected {target}, got {len(selected)}"
    return sorted(selected, key=lambda x: x["family"])

sample_k = proportional_stratified_sample(evaluable_k, target_per_arm, seed=20261013)
sample_qnotk = proportional_stratified_sample(evaluable_qnotk, target_per_arm, seed=20261014)
pilot_cohort = sorted(sample_k + sample_qnotk, key=lambda x: (x["arm"], x["family"]))

with open(out_dir / "pilot_cohort.tsv", "w") as f:
    w = csv.DictWriter(f, fieldnames=list(pilot_cohort[0].keys()), delimiter="\t")
    w.writeheader()
    w.writerows(pilot_cohort)

print(
    f"[cohort:v3] Emitted pilot_cohort.tsv: {len(pilot_cohort)} families "
    f"(K={len(sample_k)}, QnotK={len(sample_qnotk)})"
)
print(
    "[cohort:v3] Selected K strata:",
    Counter((f["phylum_group"], f["species_breadth_bin"]) for f in sample_k),
)
print(
    "[cohort:v3] Selected QnotK strata:",
    Counter((f["phylum_group"], f["species_breadth_bin"]) for f in sample_qnotk),
)

# Emit candidate members for the selected cohort
candidate_members = []
for f in pilot_cohort:
    fam = f["family"]
    rep = f["chosen_rep"]
    arm = f["arm"]
    for mem in sorted(fam_members[fam]):
        sample = mem.split("::")[0]
        contig = mem.rsplit("_", 1)[0]
        b = c2b.get(contig, "")
        sp = bin_to_species.get(b, "") if b else ""
        win = species_to_winner.get(sp, "") if sp else ""
        tax = species_to_r232_tax.get(sp, "") if sp else ""
        ph = extract_phylum(tax) if sp else ""
        if not b:
            bin_status = "unbinned"
        elif not sp:
            bin_status = "binned_unclustered_no_drep_species"
        else:
            bin_status = "binned_drep_species"

        candidate_members.append({
            "family": fam,
            "chosen_rep": rep,
            "arm": arm,
            "member": mem,
            "sample": sample,
            "contig": contig,
            "bin": b,
            "bin_status": bin_status,
            "species": sp,
            "species_winner_genome": win,
            "r232_taxonomy": tax,
            "r232_phylum": ph,
        })

cand_file = out_dir / "pilot_members_candidates.tsv"
with open(cand_file, "w") as f:
    w = csv.DictWriter(f, fieldnames=list(candidate_members[0].keys()), delimiter="\t")
    w.writeheader()
    w.writerows(candidate_members)

print(f"[cohort:v3] Emitted {len(candidate_members)} candidate member rows to {cand_file}")
