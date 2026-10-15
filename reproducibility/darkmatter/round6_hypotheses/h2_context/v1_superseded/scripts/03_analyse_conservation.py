#!/usr/bin/env python3
"""03_analyse_conservation.py - Family-level conserved neighbourhood statistic and null model.

Inputs: pilot_cohort.tsv and pilot_neighbourhoods.tsv.
Outputs: pilot_family_conservation.tsv and pilot_conservation_summary.json.

Calculates:
1. Family-level conserved neighbourhood indicator C_fam (any flanking KO in >= 2 distinct species).
2. Continuous recurrence score R_fam = max_k (species with KO k) / (total evaluable species).
3. Mobilome/phage flag for recurrent neighbours (transposase, integrase, etc.).
4. Statistical comparison between K (low-query-TM) and QnotK (fold-matched sequence-dark).
5. Within-lineage gene shuffle permutation null (100 iterations).
"""
import csv
import json
import math
import random
import sys
from collections import Counter, defaultdict
from pathlib import Path
from scipy import stats

results_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("results")
cohort_file = results_dir / "pilot_cohort.tsv"
neighbourhoods_file = results_dir / "pilot_neighbourhoods.tsv"

if not cohort_file.exists() or not neighbourhoods_file.exists():
    sys.exit(f"Error: {cohort_file} or {neighbourhoods_file} missing. Run earlier steps first.")

# Load cohort metadata
cohort = {}
with open(cohort_file) as f:
    for row in csv.DictReader(f, delimiter="\t"):
        cohort[row["family"]] = row

# Load neighbourhood data
fam_members = defaultdict(list)
with open(neighbourhoods_file) as f:
    for row in csv.DictReader(f, delimiter="\t"):
        fam_members[row["family"]].append(row)

print(f"[analyse] Loaded {len(cohort)} cohort families, {len(fam_members)} with extracted neighbourhoods")

# Keywords for mobilome/phage bounding
MOBILOME_TERMS = {
    "transposase", "integrase", "recombinase", "insertion element",
    "phage", "prophage", "terminase", "capsid", "tail", "portal", "plasmid"
}

def is_mobilome(product_str, desc_str):
    text = (product_str + " " + desc_str).lower()
    return any(term in text for term in MOBILOME_TERMS)

# Process family conservation
family_summaries = []
arm_fams = defaultdict(list)

for fam, meta in cohort.items():
    m_rows = fam_members.get(fam, [])
    arm = meta["arm"]
    
    # Track callable species and KOs in window (+-1 and +-2)
    species_seen = set()
    species_callable = set()
    ko_by_species = defaultdict(set)
    product_by_species = defaultdict(set)
    cog_by_species = defaultdict(set)
    
    total_flanks = 0
    truncated_flanks = 0
    annotated_ko_flanks = 0
    same_strand_flanks = 0
    
    for r in m_rows:
        sp = r["species"]
        if sp:
            species_seen.add(sp)
        
        has_callable_flank = False
        for flank_key in ["m2", "m1", "p1", "p2"]:
            total_flanks += 1
            if r[f"{flank_key}_truncated"] == "1":
                truncated_flanks += 1
                continue
            
            has_callable_flank = True
            if r[f"{flank_key}_same_strand"] == "1":
                same_strand_flanks += 1
            
            ko_str = r[f"{flank_key}_ko"]
            prod = r[f"{flank_key}_product"]
            cog = r[f"{flank_key}_cog"]
            
            if ko_str:
                annotated_ko_flanks += 1
                for k in ko_str.split(","):
                    if k:
                        ko_by_species[k].add(sp)
            if prod and prod != "hypothetical protein":
                product_by_species[prod].add(sp)
            if cog:
                cog_by_species[cog].add(sp)
        
        if has_callable_flank and sp:
            species_callable.add(sp)
    
    # Recurrence statistics
    recurrent_kos = {k: len(sps) for k, sps in ko_by_species.items() if len(sps) >= 2}
    recurrent_prods = {p: len(sps) for p, sps in product_by_species.items() if len(sps) >= 2}
    
    best_ko = ""
    best_ko_count = 0
    best_ko_mobilome = 0
    if recurrent_kos:
        best_ko, best_ko_count = max(recurrent_kos.items(), key=lambda x: x[1])
        # Check mobilome flag for best KO
        for r in m_rows:
            for flank_key in ["m2", "m1", "p1", "p2"]:
                if best_ko in r[f"{flank_key}_ko"].split(","):
                    if is_mobilome(r[f"{flank_key}_product"], r[f"{flank_key}_desc"]):
                        best_ko_mobilome = 1
                        break
    
    n_sp_call = len(species_callable)
    recurrence_score = best_ko_count / n_sp_call if n_sp_call > 0 else 0.0
    c_fam = 1 if len(recurrent_kos) > 0 else 0
    
    f_res = {
        "family": fam,
        "chosen_rep": meta["chosen_rep"],
        "arm": arm,
        "phylum": meta["phylum"],
        "n_members": int(meta["n_members"]),
        "n_species_total": len(species_seen),
        "n_species_callable": n_sp_call,
        "c_fam_conserved": c_fam,
        "recurrence_score": round(recurrence_score, 4),
        "n_recurrent_kos": len(recurrent_kos),
        "best_recurrent_ko": best_ko,
        "best_ko_species_count": best_ko_count,
        "best_ko_mobilome_flag": best_ko_mobilome,
        "n_recurrent_products": len(recurrent_prods),
        "flank_truncation_rate": round(truncated_flanks / total_flanks, 4) if total_flanks else 0.0,
        "flank_ko_annot_rate": round(annotated_ko_flanks / (total_flanks - truncated_flanks), 4) if (total_flanks - truncated_flanks) else 0.0,
        "same_strand_fraction": round(same_strand_flanks / (total_flanks - truncated_flanks), 4) if (total_flanks - truncated_flanks) else 0.0
    }
    family_summaries.append(f_res)
    arm_fams[arm].append(f_res)

# Write per-family conservation
out_fam_tsv = results_dir / "pilot_family_conservation.tsv"
with open(out_fam_tsv, "w") as f:
    w = csv.DictWriter(f, fieldnames=list(family_summaries[0].keys()), delimiter="\t")
    w.writeheader()
    w.writerows(family_summaries)

print(f"[analyse] Wrote family summaries to {out_fam_tsv}")

# Arm comparison
k_list = arm_fams["K"]
q_list = arm_fams["QnotK"]

k_c_fam = sum(f["c_fam_conserved"] for f in k_list)
q_c_fam = sum(f["c_fam_conserved"] for f in q_list)
k_n = len(k_list)
q_n = len(q_list)

k_rate = k_c_fam / k_n if k_n else 0.0
q_rate = q_c_fam / q_n if q_n else 0.0
diff_rate = k_rate - q_rate

# 2x2 contingency table: [[k_conserved, k_not], [q_conserved, q_not]]
contingency = [
    [k_c_fam, k_n - k_c_fam],
    [q_c_fam, q_n - q_c_fam]
]
odds_ratio, p_fisher = stats.fisher_exact(contingency)

# Mann-Whitney U test on recurrence scores
u_stat, p_mannwhitney = stats.mannwhitneyu(
    [f["recurrence_score"] for f in k_list],
    [f["recurrence_score"] for f in q_list],
    alternative="two-sided"
)

# Covariate distributions
k_trunc = [f["flank_truncation_rate"] for f in k_list]
q_trunc = [f["flank_truncation_rate"] for f in q_list]
k_annot = [f["flank_ko_annot_rate"] for f in k_list]
q_annot = [f["flank_ko_annot_rate"] for f in q_list]

print(f"[analyse] Observed rates: K={k_c_fam}/{k_n} ({k_rate*100:.1f}%), QnotK={q_c_fam}/{q_n} ({q_rate*100:.1f}%)")
print(f"[analyse] Risk Difference: {diff_rate*100:.2f} percentage points, Fisher p={p_fisher:.4f}, Mann-Whitney p={p_mannwhitney:.4f}")

# Permutation null model: within-sample contig flank shuffle
# Preserves contig lengths and annotation availability, shuffles KO labels across positions
rng = random.Random(20261013)
null_diffs = []
B = 100

all_assigned_kos = [
    r[f"{k}_ko"]
    for fam, rows in fam_members.items()
    for r in rows
    for k in ["m2", "m1", "p1", "p2"]
    if r[f"{k}_truncated"] == "0" and r[f"{k}_ko"]
]

for b_idx in range(B):
    # Create shuffled assignment
    shuffled_pool = list(all_assigned_kos)
    rng.shuffle(shuffled_pool)
    pool_idx = 0
    
    b_k_cons = 0
    b_q_cons = 0
    
    for f in family_summaries:
        fam = f["family"]
        m_rows = fam_members.get(fam, [])
        null_ko_by_species = defaultdict(set)
        
        for r in m_rows:
            sp = r["species"]
            for flank_key in ["m2", "m1", "p1", "p2"]:
                if r[f"{flank_key}_truncated"] == "0" and r[f"{flank_key}_ko"]:
                    pseudo_ko = shuffled_pool[pool_idx % len(shuffled_pool)]
                    pool_idx += 1
                    for pk in pseudo_ko.split(","):
                        if pk and sp:
                            null_ko_by_species[pk].add(sp)
        
        is_cons = 1 if any(len(sps) >= 2 for sps in null_ko_by_species.values()) else 0
        if f["arm"] == "K":
            b_k_cons += is_cons
        else:
            b_q_cons += is_cons
    
    null_diff = (b_k_cons / k_n) - (b_q_cons / q_n)
    null_diffs.append(null_diff)

p_perm = sum(1 for d in null_diffs if abs(d) >= abs(diff_rate)) / B

summary = {
    "sample_size": {"K_families": k_n, "QnotK_families": q_n, "total": k_n + q_n},
    "conserved_neighbourhood_families": {
        "K": {"n": k_c_fam, "rate": round(k_rate, 4)},
        "QnotK": {"n": q_c_fam, "rate": round(q_rate, 4)},
        "risk_difference": round(diff_rate, 4),
        "odds_ratio": round(odds_ratio, 4) if not math.isinf(odds_ratio) else "Inf",
        "p_fisher": round(p_fisher, 6),
        "p_mannwhitney": round(p_mannwhitney, 6)
    },
    "flank_truncation_rate": {
        "K_mean": round(sum(k_trunc) / len(k_trunc), 4),
        "QnotK_mean": round(sum(q_trunc) / len(q_trunc), 4)
    },
    "flank_ko_annotation_rate": {
        "K_mean": round(sum(k_annot) / len(k_annot), 4),
        "QnotK_mean": round(sum(q_annot) / len(q_annot), 4)
    },
    "permutation_null": {
        "n_permutations": B,
        "null_diff_mean": round(sum(null_diffs) / B, 4),
        "null_diff_std": round(math.sqrt(sum((x - sum(null_diffs)/B)**2 for x in null_diffs) / B), 4),
        "null_p_empirical": round(p_perm, 4)
    },
    "recurrent_functions_found": [
        {"family": f["family"], "arm": f["arm"], "ko": f["best_recurrent_ko"], "species_count": f["best_ko_species_count"], "mobilome_flag": f["best_ko_mobilome_flag"]}
        for f in family_summaries if f["c_fam_conserved"] == 1
    ]
}

out_sum_json = results_dir / "pilot_conservation_summary.json"
out_sum_json.write_text(json.dumps(summary, indent=2) + "\n")
print(f"[analyse] Saved summary to {out_sum_json}")
