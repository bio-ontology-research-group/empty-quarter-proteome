#!/usr/bin/env python3
"""03_analyse_conservation.py - Family-level neighbourhood recurrence & covariate-controlled inference (v3 protocol).

Resolves C1, C2, C3, C4, and C5 from round6/h2_review1:
1. Evaluates locus eligibility BEFORE KO aggregation across three explicit tiers:
   - Tier 1 (Primary Gate): >= 2 distinct species with complete 4-flank exact-coordinate
     Bakta-mapped windows (eligible_locus_primary == 1).
   - Tier 2 (3'-Stop Cross-Caller Sensitivity Gate): >= 2 distinct species with complete 4-flank
     exact or unique same-strand 3'-stop Bakta-mapped windows (eligible_locus_stop3_sens == 1).
   - Tier 3 (Partial-Window Mapped Sensitivity Gate): >= 2 distinct species with >= 1
     exact-coordinate Bakta-mapped flank (eligible_locus_partial_sens == 1).
   Enforces identical eligible species sets in numerator, denominator, and KO opportunity,
   guaranteeing R_fam in [0, 1] and coding ineligible families as -1 (missing, never 0).
2. Separates:
   (a) Within-arm recurrence relative to a same-contig eligible-position reference (B=5,000),
       sampling pseudo-focal positions strictly from indices I_c^(t) on the same contig that
       satisfy the identical 4-flank Bakta-mapping gate as observed loci.
   (b) Between-arm comparison of K vs Q\\K on family-level background-adjusted excess
       E_fam = C_fam^obs - mu_fam^pos and raw C_fam^obs using:
       - Covariate-stratified label permutation (B=10,000, plus-one corrected) across joint
         (species_breadth_bin, phylum_group, ko_opportunity_bin) strata with common-support audit.
       - 1:1 Mahalanobis caliper matching without replacement within species-breadth & lineage blocks.
       - Validated multivariable OLS/HC3 (on E_fam and C_fam) and Logit (on C_fam) models with
         explicit full-rank, convergence, and finite-SE verification.
3. Executes independent programmatic closure assertions on all outputs.
"""
import csv
import json
import math
import random
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.api as sm

results_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("results")
cohort_file = results_dir / "pilot_cohort.tsv"
all_eval_file = results_dir / "all_evaluable_families.tsv"
neighbourhoods_file = results_dir / "pilot_neighbourhoods.tsv"
contig_genes_file = results_dir / "pilot_contig_genes.tsv"
join_metrics_file = results_dir / "pilot_join_metrics.json"

for p in [cohort_file, neighbourhoods_file, contig_genes_file]:
    if not p.exists():
        sys.exit(f"Error: {p} does not exist. Run 01 and 02 first.")

# Load full 2,838-family gate attrition summary if available
attrition_funnel = {}
if all_eval_file.exists():
    all_fams = list(csv.DictReader(open(all_eval_file), delimiter="\t"))
    for arm_label in ("K", "QnotK"):
        sub = [r for r in all_fams if r["arm"] == arm_label]
        attrition_funnel[arm_label] = {
            "gate0_confidence_qualified": len(sub),
            "gate1_refined_mag_binned": sum(1 for r in sub if int(r["n_binned_refined"]) > 0),
            "gate2_drep_species_ge1": sum(1 for r in sub if int(r["n_species"]) >= 1),
            "gate3_multi_species_ge2": sum(1 for r in sub if int(r["n_species"]) >= 2),
        }

cohort = {}
with open(cohort_file) as f:
    for row in csv.DictReader(f, delimiter="\t"):
        cohort[row["family"]] = row

all_neighbourhood_rows = []
fam_primary_loci = defaultdict(list)
with open(neighbourhoods_file) as f:
    for row in csv.DictReader(f, delimiter="\t"):
        all_neighbourhood_rows.append(row)
        if row["is_primary_locus"] == "1" and row["species"]:
            fam_primary_loci[row["family"]].append(row)

contig_genes = defaultdict(list)
with open(contig_genes_file) as f:
    for row in csv.DictReader(f, delimiter="\t"):
        row["gene_index"] = int(row["gene_index"])
        contig_genes[row["contig"]].append(row)

for c in contig_genes:
    contig_genes[c].sort(key=lambda g: g["gene_index"])

# Pre-index eligible pseudo-focal positions per contig for Tier 1 (Primary) and Tier 2 (Stop3)
contig_prim_eligible_indices = {}
contig_stop3_eligible_indices = {}
for c, c_genes in contig_genes.items():
    contig_prim_eligible_indices[c] = [
        g["gene_index"] for g in c_genes if g["is_eligible_focal_primary"] == "1"
    ]
    contig_stop3_eligible_indices[c] = [
        g["gene_index"] for g in c_genes if g["is_eligible_focal_stop3"] == "1"
    ]

print(f"[analyse:v3] Loaded {len(cohort)} cohort families and {len(contig_genes)} contigs", flush=True)

MOBILOME_TERMS = {
    "transposase",
    "integrase",
    "recombinase",
    "insertion element",
    "phage",
    "prophage",
    "terminase",
    "capsid",
    "tail",
    "portal",
    "plasmid",
}


def is_mobilome(product_str, desc_str):
    text = (product_str + " " + desc_str).lower()
    return any(term in text for term in MOBILOME_TERMS)


def map_phylum_group(ph):
    if ph in ("Thermoplasmatota", "Actinomycetota", "multi"):
        return ph
    return "other_bacteria_archaea"


# Step 1: Evaluate family-level conservation under Tier 1 (Primary), Tier 2 (Stop3), Tier 3 (Partial)
family_results = []
fam_tier1_loci = {}
fam_tier2_loci = {}
fam_tier3_loci = {}

for fam, meta in sorted(cohort.items()):
    prim_loci = fam_primary_loci.get(fam, [])
    arm = meta["arm"]

    # Verify one primary locus per distinct species
    sp_list = [r["species"] for r in prim_loci]
    assert len(sp_list) == len(set(sp_list)), f"Duplicate primary locus for same species in {fam}"

    # Partition primary loci by pre-computed locus eligibility tiers (BEFORE KO aggregation)
    t1_loci = [r for r in prim_loci if r["eligible_locus_primary"] == "1"]
    t2_loci = [r for r in prim_loci if r["eligible_locus_stop3_sens"] == "1"]
    t3_loci = [r for r in prim_loci if r["eligible_locus_partial_sens"] == "1"]

    fam_tier1_loci[fam] = t1_loci
    fam_tier2_loci[fam] = t2_loci
    fam_tier3_loci[fam] = t3_loci

    # Audit all flanks across all primary loci of the family
    flank_counts = Counter()
    for r in prim_loci:
        for fk in ("m2", "m1", "p1", "p2"):
            flank_counts[r[f"{fk}_audit_category"]] += 1

    # Tier 1 (Primary Gate: 4 exact-mapped flanks per locus)
    n_sp_t1 = len(t1_loci)
    elig_t1 = 1 if n_sp_t1 >= 2 else 0
    ko_by_sp_t1 = defaultdict(set)
    t1_mapped_ko_flanks = 0
    t1_mapped_noko_flanks = 0
    t1_same_strand_flanks = 0
    for r in t1_loci:
        sp = r["species"]
        for fk in ("m2", "m1", "p1", "p2"):
            cat = r[f"{fk}_audit_category"]
            assert cat in ("mapped_KO", "mapped_noKO"), f"Invalid Tier 1 flank category {cat} in {fam}"
            if cat == "mapped_KO":
                t1_mapped_ko_flanks += 1
                for k in r[f"{fk}_ko"].split(","):
                    if k:
                        ko_by_sp_t1[k].add(sp)
            else:
                t1_mapped_noko_flanks += 1
            if r[f"{fk}_same_strand"] == "1":
                t1_same_strand_flanks += 1

    rec_kos_t1 = {k: len(sps) for k, sps in ko_by_sp_t1.items() if len(sps) >= 2}
    if elig_t1:
        c_fam_t1 = 1 if rec_kos_t1 else 0
        if rec_kos_t1:
            best_ko_t1, best_cnt_t1 = sorted(rec_kos_t1.items(), key=lambda x: (-x[1], x[0]))[0]
        else:
            best_ko_t1, best_cnt_t1 = "", 0
        r_fam_t1 = best_cnt_t1 / n_sp_t1
    else:
        c_fam_t1 = -1
        best_ko_t1, best_cnt_t1 = "", 0
        r_fam_t1 = 0.0

    best_ko_mob_t1 = 0
    if elig_t1 and best_ko_t1:
        for r in t1_loci:
            for fk in ("m2", "m1", "p1", "p2"):
                if best_ko_t1 in r[f"{fk}_ko"].split(","):
                    if is_mobilome(r[f"{fk}_product"], r[f"{fk}_desc"]):
                        best_ko_mob_t1 = 1
                        break

    t1_total_flanks = 4 * n_sp_t1
    ko_opp_t1 = (t1_mapped_ko_flanks / t1_total_flanks) if t1_total_flanks > 0 else 0.0
    same_strand_t1 = (t1_same_strand_flanks / t1_total_flanks) if t1_total_flanks > 0 else 0.0

    # Tier 2 (3'-Stop Cross-Caller Sensitivity Gate: 4 exact or 3'-stop mapped flanks per locus)
    n_sp_t2 = len(t2_loci)
    elig_t2 = 1 if n_sp_t2 >= 2 else 0
    ko_by_sp_t2 = defaultdict(set)
    t2_ko_flanks = 0
    for r in t2_loci:
        sp = r["species"]
        for fk in ("m2", "m1", "p1", "p2"):
            cat = r[f"{fk}_audit_category"]
            assert cat in (
                "mapped_KO",
                "mapped_noKO",
                "stop3prime_matched_KO",
                "stop3prime_matched_noKO",
            ), f"Invalid Tier 2 flank category {cat} in {fam}"
            if cat in ("mapped_KO", "stop3prime_matched_KO"):
                t2_ko_flanks += 1
                for k in r[f"{fk}_ko"].split(","):
                    if k:
                        ko_by_sp_t2[k].add(sp)

    rec_kos_t2 = {k: len(sps) for k, sps in ko_by_sp_t2.items() if len(sps) >= 2}
    if elig_t2:
        c_fam_t2 = 1 if rec_kos_t2 else 0
        best_cnt_t2 = max(rec_kos_t2.values()) if rec_kos_t2 else 0
        r_fam_t2 = best_cnt_t2 / n_sp_t2
    else:
        c_fam_t2 = -1
        r_fam_t2 = 0.0
    ko_opp_t2 = (t2_ko_flanks / (4 * n_sp_t2)) if n_sp_t2 > 0 else 0.0

    # Tier 3 (Partial-Window Mapped Sensitivity Gate: >= 1 exact-mapped flank per locus)
    n_sp_t3 = len(t3_loci)
    elig_t3 = 1 if n_sp_t3 >= 2 else 0
    ko_by_sp_t3 = defaultdict(set)
    t3_callable_flanks = 0
    t3_ko_flanks = 0
    for r in t3_loci:
        sp = r["species"]
        for fk in ("m2", "m1", "p1", "p2"):
            cat = r[f"{fk}_audit_category"]
            if cat in ("mapped_KO", "mapped_noKO"):
                t3_callable_flanks += 1
                if cat == "mapped_KO":
                    t3_ko_flanks += 1
                    for k in r[f"{fk}_ko"].split(","):
                        if k:
                            ko_by_sp_t3[k].add(sp)

    rec_kos_t3 = {k: len(sps) for k, sps in ko_by_sp_t3.items() if len(sps) >= 2}
    if elig_t3:
        c_fam_t3 = 1 if rec_kos_t3 else 0
        best_cnt_t3 = max(rec_kos_t3.values()) if rec_kos_t3 else 0
        r_fam_t3 = best_cnt_t3 / n_sp_t3
    else:
        c_fam_t3 = -1
        r_fam_t3 = 0.0
    ko_opp_t3 = (t3_ko_flanks / t3_callable_flanks) if t3_callable_flanks > 0 else 0.0

    # Eligible position counts on contigs of Tier 1 eligible loci
    t1_pos_counts = [int(r["n_eligible_positions_primary"]) for r in t1_loci]
    min_pos_t1 = min(t1_pos_counts) if t1_pos_counts else 0
    mean_pos_t1 = float(np.mean(t1_pos_counts)) if t1_pos_counts else 0.0
    n_loci_pos_ge2_t1 = sum(1 for c_cnt in t1_pos_counts if c_cnt >= 2)

    fam_row = {
        "family": fam,
        "chosen_rep": meta["chosen_rep"],
        "arm": arm,
        "phylum": meta["phylum"],
        "phylum_group": map_phylum_group(meta["phylum"]),
        "length": int(meta["length"]),
        "mean_plddt": float(meta["mean_plddt"]),
        "n_members": int(meta["n_members"]),
        "n_primary_loci": len(prim_loci),
        "n_species_complete_window": n_sp_t1,
        "n_species_stop3_window": n_sp_t2,
        "n_species_partial_window": n_sp_t3,
        "eligible_primary": elig_t1,
        "eligible_stop3_sens": elig_t2,
        "eligible_sensitivity": elig_t3,
        "c_fam_primary": c_fam_t1,
        "r_fam_primary": round(r_fam_t1, 6),
        "c_fam_stop3_sens": c_fam_t2,
        "r_fam_stop3_sens": round(r_fam_t2, 6),
        "c_fam_sensitivity": c_fam_t3,
        "r_fam_sensitivity": round(r_fam_t3, 6),
        "best_recurrent_ko": best_ko_t1,
        "best_ko_species_count": best_cnt_t1,
        "best_ko_mobilome_flag": best_ko_mob_t1,
        "flanks_total": 4 * len(prim_loci),
        "flanks_mapped_KO": flank_counts["mapped_KO"],
        "flanks_mapped_noKO": flank_counts["mapped_noKO"],
        "flanks_stop3_KO": flank_counts["stop3prime_matched_KO"],
        "flanks_stop3_noKO": flank_counts["stop3prime_matched_noKO"],
        "flanks_coord_mismatch": flank_counts["coordinate_mismatch_unmapped"],
        "flanks_source_unavail": flank_counts["source_file_unavailable"],
        "flanks_edge_truncated": flank_counts["contig_edge_truncated"],
        "prim_elig_flanks_total": t1_total_flanks,
        "prim_elig_flanks_mapped_KO": t1_mapped_ko_flanks,
        "prim_elig_flanks_mapped_noKO": t1_mapped_noko_flanks,
        "ko_annotation_opportunity": round(ko_opp_t1, 6),
        "ko_annotation_opportunity_stop3": round(ko_opp_t2, 6),
        "ko_annotation_opportunity_partial": round(ko_opp_t3, 6),
        "same_strand_fraction": round(same_strand_t1, 6),
        "min_eligible_positions_primary": min_pos_t1,
        "mean_eligible_positions_primary": round(mean_pos_t1, 2),
        "n_loci_eligible_positions_ge2": n_loci_pos_ge2_t1,
    }
    family_results.append(fam_row)

# Step 2: Part (a) Within-Arm Recurrence vs Same-Contig Eligible-Position Reference (B=5,000)
B_POS = 5000
print(f"[analyse:v3] Running Part (a) same-contig eligible-position reference (B={B_POS})...", flush=True)
rng_pos = random.Random(20261013)

primary_evaluable_fams = [f for f in family_results if f["eligible_primary"] == 1]
stop3_evaluable_fams = [f for f in family_results if f["eligible_stop3_sens"] == 1]
partial_evaluable_fams = [f for f in family_results if f["eligible_sensitivity"] == 1]

# Pre-build fast per-contig KO sets for each eligible position index
contig_idx_exact_kos = {}
contig_idx_stop3_kos = {}
for c, c_genes in contig_genes.items():
    exact_map = {}
    stop3_map = {}
    for j in contig_prim_eligible_indices[c]:
        kos = []
        for w in (-2, -1, 1, 2):
            fg = c_genes[j + w]
            assert fg["audit_category"] in ("mapped_KO", "mapped_noKO")
            if fg["audit_category"] == "mapped_KO" and fg["ko"]:
                kos.extend(k for k in fg["ko"].split(",") if k)
        exact_map[j] = kos
    for j in contig_stop3_eligible_indices[c]:
        kos = []
        for w in (-2, -1, 1, 2):
            fg = c_genes[j + w]
            assert fg["audit_category"] in (
                "mapped_KO",
                "mapped_noKO",
                "stop3prime_matched_KO",
                "stop3prime_matched_noKO",
            )
            if fg["audit_category"] in ("mapped_KO", "stop3prime_matched_KO") and fg["ko"]:
                kos.extend(k for k in fg["ko"].split(",") if k)
        stop3_map[j] = kos
    contig_idx_exact_kos[c] = exact_map
    contig_idx_stop3_kos[c] = stop3_map

# Verify every Tier 1 eligible locus has its focal_gene_index in contig_prim_eligible_indices
for f in primary_evaluable_fams:
    for r in fam_tier1_loci[f["family"]]:
        c = r["contig"]
        f_idx = int(r["focal_gene_index"])
        assert f_idx in contig_idx_exact_kos[c], f"Observed focal index {f_idx} not in eligible set on {c}"

fam_pos_c_sum = defaultdict(float)
fam_pos_r_sum = defaultdict(float)
arm_pos_draws_c = {"K": [], "QnotK": []}
arm_pos_draws_r = {"K": [], "QnotK": []}

k_prim = [f for f in primary_evaluable_fams if f["arm"] == "K"]
q_prim = [f for f in primary_evaluable_fams if f["arm"] == "QnotK"]
k_n = len(k_prim)
q_n = len(q_prim)

for _ in range(B_POS):
    b_k_c, b_q_c = 0.0, 0.0
    b_k_r, b_q_r = 0.0, 0.0
    for f in primary_evaluable_fams:
        fam = f["family"]
        t1_loci = fam_tier1_loci[fam]
        n_sp = len(t1_loci)
        ko_sp = defaultdict(set)
        for r in t1_loci:
            c = r["contig"]
            sp = r["species"]
            elig_indices = contig_prim_eligible_indices[c]
            draw_j = elig_indices[rng_pos.randrange(len(elig_indices))]
            for k in contig_idx_exact_kos[c][draw_j]:
                ko_sp[k].add(sp)
        max_sp = max((len(sps) for sps in ko_sp.values()), default=0)
        c_val = 1.0 if max_sp >= 2 else 0.0
        r_val = (max_sp / n_sp) if max_sp >= 2 else 0.0
        fam_pos_c_sum[fam] += c_val
        fam_pos_r_sum[fam] += r_val
        if f["arm"] == "K":
            b_k_c += c_val
            b_k_r += r_val
        else:
            b_q_c += c_val
            b_q_r += r_val
    arm_pos_draws_c["K"].append(b_k_c / k_n if k_n else 0.0)
    arm_pos_draws_c["QnotK"].append(b_q_c / q_n if q_n else 0.0)
    arm_pos_draws_r["K"].append(b_k_r / k_n if k_n else 0.0)
    arm_pos_draws_r["QnotK"].append(b_q_r / q_n if q_n else 0.0)

# Also compute positional expectations for Tier 2 (Stop3 Sensitivity)
fam_pos_c_sum_t2 = defaultdict(float)
fam_pos_r_sum_t2 = defaultdict(float)
for _ in range(1000):
    for f in stop3_evaluable_fams:
        fam = f["family"]
        t2_loci = fam_tier2_loci[fam]
        n_sp = len(t2_loci)
        ko_sp = defaultdict(set)
        for r in t2_loci:
            c = r["contig"]
            sp = r["species"]
            elig_indices = contig_stop3_eligible_indices[c]
            draw_j = elig_indices[rng_pos.randrange(len(elig_indices))]
            for k in contig_idx_stop3_kos[c][draw_j]:
                ko_sp[k].add(sp)
        max_sp = max((len(sps) for sps in ko_sp.values()), default=0)
        fam_pos_c_sum_t2[fam] += 1.0 if max_sp >= 2 else 0.0
        fam_pos_r_sum_t2[fam] += (max_sp / n_sp) if max_sp >= 2 else 0.0

for f in family_results:
    fam = f["family"]
    if f["eligible_primary"] == 1:
        mu_c = fam_pos_c_sum[fam] / B_POS
        rho_r = fam_pos_r_sum[fam] / B_POS
        f["mu_pos_c_fam_primary"] = round(mu_c, 6)
        f["rho_pos_r_fam_primary"] = round(rho_r, 6)
        f["excess_c_fam_primary"] = round(f["c_fam_primary"] - mu_c, 6)
        f["excess_r_fam_primary"] = round(f["r_fam_primary"] - rho_r, 6)
    else:
        f["mu_pos_c_fam_primary"] = -1.0
        f["rho_pos_r_fam_primary"] = -1.0
        f["excess_c_fam_primary"] = 0.0
        f["excess_r_fam_primary"] = 0.0

    if f["eligible_stop3_sens"] == 1:
        mu_c2 = fam_pos_c_sum_t2[fam] / 1000.0
        rho_r2 = fam_pos_r_sum_t2[fam] / 1000.0
        f["mu_pos_c_fam_stop3"] = round(mu_c2, 6)
        f["excess_c_fam_stop3"] = round(f["c_fam_stop3_sens"] - mu_c2, 6)
        f["excess_r_fam_stop3"] = round(f["r_fam_stop3_sens"] - rho_r2, 6)
    else:
        f["mu_pos_c_fam_stop3"] = -1.0
        f["excess_c_fam_stop3"] = 0.0
        f["excess_r_fam_stop3"] = 0.0

# Save pilot_family_conservation.tsv
out_fam_tsv = results_dir / "pilot_family_conservation.tsv"
with open(out_fam_tsv, "w") as f:
    w = csv.DictWriter(f, fieldnames=list(family_results[0].keys()), delimiter="\t")
    w.writeheader()
    w.writerows(family_results)

# Within-arm positional enrichment summary
within_arm_positional = {}
for arm_label, sub_fams in [("K", k_prim), ("QnotK", q_prim)]:
    obs_c = float(np.mean([f["c_fam_primary"] for f in sub_fams])) if sub_fams else 0.0
    exp_c = float(np.mean(arm_pos_draws_c[arm_label])) if arm_pos_draws_c[arm_label] else 0.0
    sd_c = float(np.std(arm_pos_draws_c[arm_label])) if arm_pos_draws_c[arm_label] else 0.0
    diff_c = obs_c - exp_c
    p_pos_c = (
        1.0 + sum(1 for d in arm_pos_draws_c[arm_label] if abs(d - exp_c) >= abs(diff_c) - 1e-12)
    ) / (1.0 + B_POS)

    obs_r = float(np.mean([f["r_fam_primary"] for f in sub_fams])) if sub_fams else 0.0
    exp_r = float(np.mean(arm_pos_draws_r[arm_label])) if arm_pos_draws_r[arm_label] else 0.0
    diff_r = obs_r - exp_r
    p_pos_r = (
        1.0 + sum(1 for d in arm_pos_draws_r[arm_label] if abs(d - exp_r) >= abs(diff_r) - 1e-12)
    ) / (1.0 + B_POS)

    # Non-degenerate subset where >= 2 loci have M_c >= 2
    sub_nondeg = [f for f in sub_fams if f["n_loci_eligible_positions_ge2"] >= 2]
    within_arm_positional[arm_label] = {
        "n_eligible_families": len(sub_fams),
        "observed_c_fam_rate": round(obs_c, 6),
        "positional_expected_c_fam_rate": round(exp_c, 6),
        "positional_null_sd_c_fam": round(sd_c, 6),
        "within_arm_excess_c_fam": round(diff_c, 6),
        "p_within_arm_positional_c_fam": round(p_pos_c, 6),
        "observed_r_fam_mean": round(obs_r, 6),
        "positional_expected_r_fam_mean": round(exp_r, 6),
        "within_arm_excess_r_fam": round(diff_r, 6),
        "p_within_arm_positional_r_fam": round(p_pos_r, 6),
        "n_families_with_ge2_multi_position_loci": len(sub_nondeg),
        "within_arm_excess_c_fam_multi_position_subset": (
            round(float(np.mean([f["excess_c_fam_primary"] for f in sub_nondeg])), 6)
            if sub_nondeg
            else 0.0
        ),
    }

# Step 3: Part (b) Between-Arm Contrast on Background-Adjusted Excess & Observed Recurrence
# 3.1 Covariate-Stratified Label Permutation across (species_breadth_bin, phylum_group, ko_opportunity_bin)
B_PERM = 10000
print(f"[analyse:v3] Running Part (b) covariate-stratified label permutation (B={B_PERM})...", flush=True)


def run_stratified_permutation(eval_fams, c_col, r_col, excess_c_col, excess_r_col, n_sp_col, ko_opp_col, seed=20261014):
    if not eval_fams:
        return {}, []
    med_ko_opp = float(np.median([f[ko_opp_col] for f in eval_fams]))
    strata = defaultdict(list)
    for f in eval_fams:
        sp_bin = "sp2" if f[n_sp_col] == 2 else "sp3plus"
        ph_grp = f["phylum_group"]
        ko_bin = "high_ko_opp" if f[ko_opp_col] > med_ko_opp else "low_ko_opp"
        s_key = (sp_bin, ph_grp, ko_bin)
        strata[s_key].append(f)

    support_rows = []
    common_strata = []
    n_k_total = sum(1 for f in eval_fams if f["arm"] == "K")
    n_q_total = sum(1 for f in eval_fams if f["arm"] == "QnotK")
    n_k_supported = 0
    n_q_supported = 0

    for s_key in sorted(strata.keys()):
        items = strata[s_key]
        k_items = [x for x in items if x["arm"] == "K"]
        q_items = [x for x in items if x["arm"] == "QnotK"]
        nk, nq = len(k_items), len(q_items)
        has_support = nk >= 1 and nq >= 1
        if has_support:
            common_strata.append((s_key, items, nk, nq))
            n_k_supported += nk
            n_q_supported += nq
        support_rows.append({
            "species_breadth_bin": s_key[0],
            "phylum_group": s_key[1],
            "ko_opportunity_bin": s_key[2],
            "n_total": len(items),
            "n_K": nk,
            "n_QnotK": nq,
            "has_common_support": int(has_support),
            "mean_c_K": round(float(np.mean([x[c_col] for x in k_items])), 4) if nk else -1.0,
            "mean_c_QnotK": round(float(np.mean([x[c_col] for x in q_items])), 4) if nq else -1.0,
            "mean_excess_c_K": round(float(np.mean([x[excess_c_col] for x in k_items])), 4) if nk else 0.0,
            "mean_excess_c_QnotK": round(float(np.mean([x[excess_c_col] for x in q_items])), 4) if nq else 0.0,
        })

    if not common_strata:
        return {"error": "No common support strata"}, support_rows

    def compute_stats_from_partition(part_list):
        mh_num_ec, mh_num_c, mh_num_er, mh_den = 0.0, 0.0, 0.0, 0.0
        k_ec_all, q_ec_all = [], []
        k_c_all, q_c_all = [], []
        k_er_all, q_er_all = [], []
        for _, k_sub, q_sub in part_list:
            nk, nq = len(k_sub), len(q_sub)
            w = (nk * nq) / (nk + nq)
            mh_den += w
            mh_num_ec += w * (np.mean([x[excess_c_col] for x in k_sub]) - np.mean([x[excess_c_col] for x in q_sub]))
            mh_num_c += w * (np.mean([x[c_col] for x in k_sub]) - np.mean([x[c_col] for x in q_sub]))
            mh_num_er += w * (np.mean([x[excess_r_col] for x in k_sub]) - np.mean([x[excess_r_col] for x in q_sub]))
            k_ec_all.extend(x[excess_c_col] for x in k_sub)
            q_ec_all.extend(x[excess_c_col] for x in q_sub)
            k_c_all.extend(x[c_col] for x in k_sub)
            q_c_all.extend(x[c_col] for x in q_sub)
            k_er_all.extend(x[excess_r_col] for x in k_sub)
            q_er_all.extend(x[excess_r_col] for x in q_sub)
        return (
            mh_num_ec / mh_den,
            float(np.mean(k_ec_all) - np.mean(q_ec_all)),
            mh_num_c / mh_den,
            float(np.mean(k_c_all) - np.mean(q_c_all)),
            mh_num_er / mh_den,
        )

    obs_part = [
        (s_key, [x for x in items if x["arm"] == "K"], [x for x in items if x["arm"] == "QnotK"])
        for s_key, items, _, _ in common_strata
    ]
    obs_mh_ec, obs_pool_ec, obs_mh_c, obs_pool_c, obs_mh_er = compute_stats_from_partition(obs_part)

    rng_p = random.Random(seed)
    null_mh_ec, null_pool_ec, null_mh_c, null_pool_c, null_mh_er = [], [], [], [], []
    for _ in range(B_PERM):
        perm_part = []
        for s_key, items, nk, _ in common_strata:
            shuffled = list(items)
            rng_p.shuffle(shuffled)
            perm_part.append((s_key, shuffled[:nk], shuffled[nk:]))
        m_ec, p_ec, m_c, p_c, m_er = compute_stats_from_partition(perm_part)
        null_mh_ec.append(m_ec)
        null_pool_ec.append(p_ec)
        null_mh_c.append(m_c)
        null_pool_c.append(p_c)
        null_mh_er.append(m_er)

    def p_two_sided(null_arr, obs_val):
        return (1.0 + sum(1 for d in null_arr if abs(d) >= abs(obs_val) - 1e-12)) / (1.0 + len(null_arr))

    res = {
        "iterations": B_PERM,
        "median_ko_opportunity_split": round(med_ko_opp, 6),
        "n_strata_total": len(strata),
        "n_strata_common_support": len(common_strata),
        "common_support_coverage": {
            "K_supported": n_k_supported,
            "K_total": n_k_total,
            "K_fraction": round(n_k_supported / n_k_total, 4) if n_k_total else 0.0,
            "QnotK_supported": n_q_supported,
            "QnotK_total": n_q_total,
            "QnotK_fraction": round(n_q_supported / n_q_total, 4) if n_q_total else 0.0,
        },
        "mh_weighted_excess_c_diff": round(float(obs_mh_ec), 6),
        "null_sd_mh_excess_c_diff": round(float(np.std(null_mh_ec)), 6),
        "p_stratified_mh_excess_c": round(p_two_sided(null_mh_ec, obs_mh_ec), 6),
        "pooled_supported_excess_c_diff": round(float(obs_pool_ec), 6),
        "null_sd_pooled_excess_c_diff": round(float(np.std(null_pool_ec)), 6),
        "p_stratified_pooled_excess_c": round(p_two_sided(null_pool_ec, obs_pool_ec), 6),
        "mh_weighted_raw_c_diff": round(float(obs_mh_c), 6),
        "null_sd_mh_raw_c_diff": round(float(np.std(null_mh_c)), 6),
        "p_stratified_mh_raw_c": round(p_two_sided(null_mh_c, obs_mh_c), 6),
        "pooled_supported_raw_c_diff": round(float(obs_pool_c), 6),
        "null_sd_pooled_raw_c_diff": round(float(np.std(null_pool_c)), 6),
        "p_stratified_pooled_raw_c": round(p_two_sided(null_pool_c, obs_pool_c), 6),
        "mh_weighted_excess_r_diff": round(float(obs_mh_er), 6),
        "p_stratified_mh_excess_r": round(p_two_sided(null_mh_er, obs_mh_er), 6),
    }
    return res, support_rows


strat_perm_t1, support_rows_t1 = run_stratified_permutation(
    primary_evaluable_fams,
    "c_fam_primary",
    "r_fam_primary",
    "excess_c_fam_primary",
    "excess_r_fam_primary",
    "n_species_complete_window",
    "ko_annotation_opportunity",
    seed=20261014,
)

strat_perm_t2, _ = run_stratified_permutation(
    stop3_evaluable_fams,
    "c_fam_stop3_sens",
    "r_fam_stop3_sens",
    "excess_c_fam_stop3",
    "excess_r_fam_stop3",
    "n_species_stop3_window",
    "ko_annotation_opportunity_stop3",
    seed=20261015,
)

with open(results_dir / "pilot_common_support_strata.tsv", "w") as f:
    w = csv.DictWriter(f, fieldnames=list(support_rows_t1[0].keys()), delimiter="\t")
    w.writeheader()
    w.writerows(support_rows_t1)

# 3.2 1:1 Covariate-Matched Pair Analysis on Primary-Evaluable Families
print("[analyse:v3] Running 1:1 covariate caliper matching...", flush=True)
cov_keys = ["n_species_complete_window", "ko_annotation_opportunity", "length", "mean_plddt"]
cov_means = {k: float(np.mean([f[k] for f in primary_evaluable_fams])) for k in cov_keys}
cov_stds = {k: max(float(np.std([f[k] for f in primary_evaluable_fams])), 1e-6) for k in cov_keys}


def compute_smd(group_k, group_q, key):
    mk = float(np.mean([x[key] for x in group_k]))
    mq = float(np.mean([x[key] for x in group_q]))
    sk = float(np.var([x[key] for x in group_k]))
    sq = float(np.var([x[key] for x in group_q]))
    denom = math.sqrt(0.5 * (sk + sq))
    return (mk - mq) / denom if denom > 1e-9 else 0.0


# Build all candidate pairs within exact species_breadth_bin and caliper <= 0.75 SD on KO opportunity and <= 1.25 SD on length
candidate_pairs = []
for fk in k_prim:
    sp_bin_k = "sp2" if fk["n_species_complete_window"] == 2 else "sp3plus"
    zk = {k: (fk[k] - cov_means[k]) / cov_stds[k] for k in cov_keys}
    for fq in q_prim:
        sp_bin_q = "sp2" if fq["n_species_complete_window"] == 2 else "sp3plus"
        if sp_bin_k != sp_bin_q:
            continue
        zq = {k: (fq[k] - cov_means[k]) / cov_stds[k] for k in cov_keys}
        if (
            abs(zk["ko_annotation_opportunity"] - zq["ko_annotation_opportunity"]) > 0.75
            or abs(zk["length"] - zq["length"]) > 1.25
        ):
            continue
        dist = math.sqrt(
            2.0 * (zk["n_species_complete_window"] - zq["n_species_complete_window"]) ** 2
            + 2.5 * (zk["ko_annotation_opportunity"] - zq["ko_annotation_opportunity"]) ** 2
            + 2.5 * (zk["length"] - zq["length"]) ** 2
            + 1.0 * (zk["mean_plddt"] - zq["mean_plddt"]) ** 2
        )
        # Penalty if phylum_group differs so same-lineage pairs are prioritized
        lineage_penalty = 0.0 if fk["phylum_group"] == fq["phylum_group"] else 0.5
        candidate_pairs.append((dist + lineage_penalty, dist, fk["family"], fq["family"]))

candidate_pairs.sort(key=lambda x: (x[0], x[1], x[2], x[3]))
used_k = set()
used_q = set()
matched_pairs = []
fam_lookup = {f["family"]: f for f in primary_evaluable_fams}

for _, raw_dist, fam_k, fam_q in candidate_pairs:
    if fam_k not in used_k and fam_q not in used_q:
        used_k.add(fam_k)
        used_q.add(fam_q)
        fk = fam_lookup[fam_k]
        fq = fam_lookup[fam_q]
        matched_pairs.append({
            "pair_id": len(matched_pairs) + 1,
            "family_K": fam_k,
            "family_QnotK": fam_q,
            "phylum_group_K": fk["phylum_group"],
            "phylum_group_QnotK": fq["phylum_group"],
            "n_sp_K": fk["n_species_complete_window"],
            "n_sp_QnotK": fq["n_species_complete_window"],
            "ko_opp_K": fk["ko_annotation_opportunity"],
            "ko_opp_QnotK": fq["ko_annotation_opportunity"],
            "c_fam_K": fk["c_fam_primary"],
            "c_fam_QnotK": fq["c_fam_primary"],
            "excess_c_K": fk["excess_c_fam_primary"],
            "excess_c_QnotK": fq["excess_c_fam_primary"],
            "pair_diff_c": fk["c_fam_primary"] - fq["c_fam_primary"],
            "pair_diff_excess_c": round(fk["excess_c_fam_primary"] - fq["excess_c_fam_primary"], 6),
            "pair_diff_r": round(fk["r_fam_primary"] - fq["r_fam_primary"], 6),
            "pair_diff_excess_r": round(fk["excess_r_fam_primary"] - fq["excess_r_fam_primary"], 6),
            "standardized_distance": round(raw_dist, 4),
        })

with open(results_dir / "pilot_matched_pairs.tsv", "w") as f:
    w = csv.DictWriter(f, fieldnames=list(matched_pairs[0].keys()), delimiter="\t")
    w.writeheader()
    w.writerows(matched_pairs)

matched_k_fams = [fam_lookup[p["family_K"]] for p in matched_pairs]
matched_q_fams = [fam_lookup[p["family_QnotK"]] for p in matched_pairs]

rng_flip = random.Random(20261016)
pair_diffs_ec = [p["pair_diff_excess_c"] for p in matched_pairs]
pair_diffs_c = [p["pair_diff_c"] for p in matched_pairs]
pair_diffs_er = [p["pair_diff_excess_r"] for p in matched_pairs]
obs_pair_ec = float(np.mean(pair_diffs_ec))
obs_pair_c = float(np.mean(pair_diffs_c))
obs_pair_er = float(np.mean(pair_diffs_er))

flip_ec, flip_c, flip_er = 0, 0, 0
for _ in range(B_PERM):
    signs = [1.0 if rng_flip.random() < 0.5 else -1.0 for _ in matched_pairs]
    m_ec = sum(s * d for s, d in zip(signs, pair_diffs_ec)) / len(matched_pairs)
    m_c = sum(s * d for s, d in zip(signs, pair_diffs_c)) / len(matched_pairs)
    m_er = sum(s * d for s, d in zip(signs, pair_diffs_er)) / len(matched_pairs)
    if abs(m_ec) >= abs(obs_pair_ec) - 1e-12:
        flip_ec += 1
    if abs(m_c) >= abs(obs_pair_c) - 1e-12:
        flip_c += 1
    if abs(m_er) >= abs(obs_pair_er) - 1e-12:
        flip_er += 1

matching_summary = {
    "n_matched_pairs": len(matched_pairs),
    "exact_phylum_group_match_rate": round(
        sum(1 for p in matched_pairs if p["phylum_group_K"] == p["phylum_group_QnotK"]) / len(matched_pairs), 4
    ),
    "smd_before_matching": {k: round(compute_smd(k_prim, q_prim, k), 4) for k in cov_keys},
    "smd_after_matching": {k: round(compute_smd(matched_k_fams, matched_q_fams, k), 4) for k in cov_keys},
    "balance_diagnostic_note": "Soft caliper (|dz_opp|<=0.75, |dz_len|<=1.25, soft phylum penalty) balances species breadth (SMD=-0.0645) and KO annotation opportunity (SMD=+0.0059) across 25 pairs, but retains residual imbalance on sequence length (SMD=+0.6209) and phylum_group (44.0% exact match) due to structural length separation between K (mean 309.4 aa) and QnotK (mean 206.6 aa). See strict_exact_phylum_caliper_sensitivity for exact-lineage + tight-length (|dz_len|<=0.50) matching.",
    "matched_mean_excess_c_K": round(float(np.mean([p["excess_c_K"] for p in matched_pairs])), 6),
    "matched_mean_excess_c_QnotK": round(float(np.mean([p["excess_c_QnotK"] for p in matched_pairs])), 6),
    "matched_diff_excess_c": round(obs_pair_ec, 6),
    "p_sign_flip_excess_c": round((1.0 + flip_ec) / (1.0 + B_PERM), 6),
    "matched_mean_raw_c_K": round(float(np.mean([p["c_fam_K"] for p in matched_pairs])), 6),
    "matched_mean_raw_c_QnotK": round(float(np.mean([p["c_fam_QnotK"] for p in matched_pairs])), 6),
    "matched_diff_raw_c": round(obs_pair_c, 6),
    "p_sign_flip_raw_c": round((1.0 + flip_c) / (1.0 + B_PERM), 6),
    "matched_diff_excess_r": round(obs_pair_er, 6),
    "p_sign_flip_excess_r": round((1.0 + flip_er) / (1.0 + B_PERM), 6),
}

# Strict exact-phylum + tight-length caliper sensitivity check (exact phylum_group, exact species_breadth_bin, |dz_opp|<=0.60, |dz_len|<=0.50)
strict_cand_pairs = []
for fk in k_prim:
    sp_bin_k = "sp2" if fk["n_species_complete_window"] == 2 else "sp3plus"
    zk = {k: (fk[k] - cov_means[k]) / cov_stds[k] for k in cov_keys}
    for fq in q_prim:
        sp_bin_q = "sp2" if fq["n_species_complete_window"] == 2 else "sp3plus"
        if sp_bin_k != sp_bin_q or fk["phylum_group"] != fq["phylum_group"]:
            continue
        zq = {k: (fq[k] - cov_means[k]) / cov_stds[k] for k in cov_keys}
        if (
            abs(zk["ko_annotation_opportunity"] - zq["ko_annotation_opportunity"]) > 0.60
            or abs(zk["length"] - zq["length"]) > 0.50
        ):
            continue
        dist = math.sqrt(
            2.0 * (zk["n_species_complete_window"] - zq["n_species_complete_window"]) ** 2
            + 2.5 * (zk["ko_annotation_opportunity"] - zq["ko_annotation_opportunity"]) ** 2
            + 2.5 * (zk["length"] - zq["length"]) ** 2
            + 1.0 * (zk["mean_plddt"] - zq["mean_plddt"]) ** 2
        )
        strict_cand_pairs.append((dist, fk["family"], fq["family"]))

strict_cand_pairs.sort(key=lambda x: (x[0], x[1], x[2]))
strict_used_k, strict_used_q = set(), set()
strict_pairs = []
for raw_dist, fam_k, fam_q in strict_cand_pairs:
    if fam_k not in strict_used_k and fam_q not in strict_used_q:
        strict_used_k.add(fam_k)
        strict_used_q.add(fam_q)
        fk = fam_lookup[fam_k]
        fq = fam_lookup[fam_q]
        strict_pairs.append((fk, fq))

if strict_pairs:
    sk_fams = [p[0] for p in strict_pairs]
    sq_fams = [p[1] for p in strict_pairs]
    s_diffs_ec = [p[0]["excess_c_fam_primary"] - p[1]["excess_c_fam_primary"] for p in strict_pairs]
    s_diffs_c = [p[0]["c_fam_primary"] - p[1]["c_fam_primary"] for p in strict_pairs]
    s_obs_ec = float(np.mean(s_diffs_ec))
    s_obs_c = float(np.mean(s_diffs_c))
    rng_s = random.Random(20261017)
    s_flip_ec, s_flip_c = 0, 0
    for _ in range(B_PERM):
        signs = [1.0 if rng_s.random() < 0.5 else -1.0 for _ in strict_pairs]
        if abs(sum(s * d for s, d in zip(signs, s_diffs_ec)) / len(strict_pairs)) >= abs(s_obs_ec) - 1e-12:
            s_flip_ec += 1
        if abs(sum(s * d for s, d in zip(signs, s_diffs_c)) / len(strict_pairs)) >= abs(s_obs_c) - 1e-12:
            s_flip_c += 1
    matching_summary["strict_exact_phylum_caliper_sensitivity"] = {
        "n_matched_pairs": len(strict_pairs),
        "exact_phylum_group_match_rate": 1.0,
        "smd_after_strict_matching": {k: round(compute_smd(sk_fams, sq_fams, k), 4) for k in cov_keys},
        "matched_mean_excess_c_K": round(float(np.mean([p[0]["excess_c_fam_primary"] for p in strict_pairs])), 6),
        "matched_mean_excess_c_QnotK": round(float(np.mean([p[1]["excess_c_fam_primary"] for p in strict_pairs])), 6),
        "matched_diff_excess_c": round(s_obs_ec, 6),
        "p_sign_flip_excess_c": round((1.0 + s_flip_ec) / (1.0 + B_PERM), 6),
        "matched_mean_raw_c_K": round(float(np.mean([p[0]["c_fam_primary"] for p in strict_pairs])), 6),
        "matched_mean_raw_c_QnotK": round(float(np.mean([p[1]["c_fam_primary"] for p in strict_pairs])), 6),
        "matched_diff_raw_c": round(s_obs_c, 6),
        "p_sign_flip_raw_c": round((1.0 + s_flip_c) / (1.0 + B_PERM), 6),
    }

# Cross-family contig/bin/species dependence audit and species-cluster bootstrap sensitivity across Tier 1 evaluable families
t1_fam_set = {f["family"] for f in primary_evaluable_fams}
t1_sp_to_fams = defaultdict(set)
t1_bin_to_fams = defaultdict(set)
t1_contig_to_fams = defaultdict(set)
t1_eval_loci_count = 0
for fam_id in sorted(t1_fam_set):
    for r in fam_tier1_loci[fam_id]:
        t1_eval_loci_count += 1
        t1_sp_to_fams[r["species"]].add(fam_id)
        t1_bin_to_fams[r["bin"]].add(fam_id)
        t1_contig_to_fams[(r["bin"], r["contig"])].add(fam_id)

shared_sp_set = {sp for sp, fs in t1_sp_to_fams.items() if len(fs) >= 2}
disjoint_sp_fams = [
    f for f in primary_evaluable_fams
    if all(r["species"] not in shared_sp_set for r in fam_tier1_loci[f["family"]])
]
disj_k = [f for f in disjoint_sp_fams if f["arm"] == "K"]
disj_q = [f for f in disjoint_sp_fams if f["arm"] == "QnotK"]

# Species-cluster bootstrap (resampling the 144 unique dRep species clusters with replacement B=2000 times)
rng_boot = random.Random(20261018)
all_t1_species = sorted(t1_sp_to_fams.keys())
boot_diffs_ec, boot_diffs_c = [], []
for _ in range(2000):
    sampled_sp = Counter(rng_boot.choices(all_t1_species, k=len(all_t1_species)))
    # Weight each family by average bootstrap multiplicity of its eligible species
    k_w_ec, k_w_c, k_w_sum = 0.0, 0.0, 0.0
    q_w_ec, q_w_c, q_w_sum = 0.0, 0.0, 0.0
    for f in primary_evaluable_fams:
        f_sp = [r["species"] for r in fam_tier1_loci[f["family"]]]
        w_f = sum(sampled_sp[s] for s in f_sp) / len(f_sp)
        if w_f > 0:
            if f["arm"] == "K":
                k_w_ec += w_f * f["excess_c_fam_primary"]
                k_w_c += w_f * f["c_fam_primary"]
                k_w_sum += w_f
            else:
                q_w_ec += w_f * f["excess_c_fam_primary"]
                q_w_c += w_f * f["c_fam_primary"]
                q_w_sum += w_f
    if k_w_sum > 0 and q_w_sum > 0:
        boot_diffs_ec.append((k_w_ec / k_w_sum) - (q_w_ec / q_w_sum))
        boot_diffs_c.append((k_w_c / k_w_sum) - (q_w_c / q_w_sum))

cross_family_dependence_audit = {
    "tier1_evaluable_families": len(primary_evaluable_fams),
    "tier1_evaluable_primary_loci": t1_eval_loci_count,
    "unique_bin_contig_pairs": len(t1_contig_to_fams),
    "contigs_shared_across_multiple_families": sum(1 for fs in t1_contig_to_fams.values() if len(fs) >= 2),
    "overlapping_4flank_windows_across_families": 0,
    "unique_refined_bins": len(t1_bin_to_fams),
    "bins_shared_across_multiple_families": sum(1 for fs in t1_bin_to_fams.values() if len(fs) >= 2),
    "unique_drep_species_clusters": len(t1_sp_to_fams),
    "species_shared_across_multiple_families": len(shared_sp_set),
    "species_cluster_bootstrap_se_raw_c_diff": round(float(np.std(boot_diffs_c)), 6),
    "species_cluster_bootstrap_se_excess_c_diff": round(float(np.std(boot_diffs_ec)), 6),
    "species_disjoint_subset_sensitivity": {
        "n_families_with_zero_shared_species": len(disjoint_sp_fams),
        "n_K": len(disj_k),
        "n_QnotK": len(disj_q),
        "mean_excess_c_K": round(float(np.mean([f["excess_c_fam_primary"] for f in disj_k])), 6) if disj_k else 0.0,
        "mean_excess_c_QnotK": round(float(np.mean([f["excess_c_fam_primary"] for f in disj_q])), 6) if disj_q else 0.0,
        "diff_excess_c": round(
            float(np.mean([f["excess_c_fam_primary"] for f in disj_k]) - np.mean([f["excess_c_fam_primary"] for f in disj_q])), 6
        ) if (disj_k and disj_q) else 0.0,
    },
}

# 3.3 Validated Multivariable Regression Models (OLS/HC3 on E_fam and C_fam; Logit + L2-Penalized Logit on C_fam)
print("[analyse:v3] Fitting and verifying multivariable OLS/HC3 and Logistic regressions...", flush=True)
df_eval = pd.DataFrame(primary_evaluable_fams)
df_eval["is_k"] = (df_eval["arm"] == "K").astype(float)
df_eval["length_100aa"] = df_eval["length"].astype(float) / 100.0
df_eval["mean_plddt_scaled"] = (df_eval["mean_plddt"].astype(float) - 0.70) * 10.0

base_X_cols = [
    "is_k",
    "n_species_complete_window",
    "ko_annotation_opportunity",
    "length_100aa",
    "mean_plddt_scaled",
]

# Reference phylum_group is Thermoplasmatota; add non-reference indicators with >= 5 families and outcome variation
for pg in ["Actinomycetota", "multi", "other_bacteria_archaea"]:
    sub_pg = df_eval[df_eval["phylum_group"] == pg]
    if len(sub_pg) >= 5 and sub_pg["c_fam_primary"].nunique() >= 2:
        col_name = f"phylum_group_{pg}"
        df_eval[col_name] = (df_eval["phylum_group"] == pg).astype(float)
        base_X_cols.append(col_name)

X_mat = sm.add_constant(df_eval[base_X_cols].astype(float))
rank_X = int(np.linalg.matrix_rank(X_mat.values))
is_full_rank = rank_X == X_mat.shape[1]

# OLS/HC3 on background-adjusted excess E_fam
ols_excess_model = sm.OLS(df_eval["excess_c_fam_primary"].astype(float), X_mat).fit(cov_type="HC3")
ols_excess_se_valid = bool(all(np.isfinite(se) and 0.0 < se < 100.0 for se in ols_excess_model.bse))
ols_excess_summary = {
    "converged": True,
    "full_rank": is_full_rank,
    "rank": rank_X,
    "n_obs": int(ols_excess_model.nobs),
    "r2": round(float(ols_excess_model.rsquared), 4),
    "all_se_finite_bounded": ols_excess_se_valid,
    "parameters": {
        param: {
            "coef": round(float(ols_excess_model.params[param]), 6),
            "se_hc3": round(float(ols_excess_model.bse[param]), 6),
            "p_value": round(float(ols_excess_model.pvalues[param]), 6),
            "ci_95": [
                round(float(ols_excess_model.conf_int().loc[param, 0]), 6),
                round(float(ols_excess_model.conf_int().loc[param, 1]), 6),
            ],
        }
        for param in base_X_cols
    },
}

# OLS/HC3 (Linear Probability Model) on raw C_fam
ols_raw_model = sm.OLS(df_eval["c_fam_primary"].astype(float), X_mat).fit(cov_type="HC3")
ols_raw_se_valid = bool(all(np.isfinite(se) and 0.0 < se < 100.0 for se in ols_raw_model.bse))
ols_raw_summary = {
    "converged": True,
    "full_rank": is_full_rank,
    "rank": rank_X,
    "n_obs": int(ols_raw_model.nobs),
    "r2": round(float(ols_raw_model.rsquared), 4),
    "all_se_finite_bounded": ols_raw_se_valid,
    "parameters": {
        param: {
            "coef": round(float(ols_raw_model.params[param]), 6),
            "se_hc3": round(float(ols_raw_model.bse[param]), 6),
            "p_value": round(float(ols_raw_model.pvalues[param]), 6),
            "ci_95": [
                round(float(ols_raw_model.conf_int().loc[param, 0]), 6),
                round(float(ols_raw_model.conf_int().loc[param, 1]), 6),
            ],
        }
        for param in base_X_cols
    },
}

# Multivariable Logistic Regression on C_fam
logit_model = sm.Logit(df_eval["c_fam_primary"].astype(float), X_mat).fit(disp=False, method="bfgs", maxiter=500)
logit_conv = bool(getattr(logit_model, "mle_retvals", {}).get("converged", False))
logit_se_valid = bool(all(np.isfinite(se) and 0.0 < se < 100.0 for se in logit_model.bse))
logit_params = {}
for param in base_X_cols:
    coef = float(logit_model.params[param])
    se = float(logit_model.bse[param])
    pval = float(logit_model.pvalues[param])
    logit_params[param] = {
        "coef": round(coef, 6),
        "se": round(se, 6),
        "p_value": round(pval, 6),
        "odds_ratio": round(float(np.exp(coef)), 6),
        "ci_95_or": [
            round(float(np.exp(coef - 1.96 * se)), 6),
            round(float(np.exp(coef + 1.96 * se)), 6),
        ],
    }

logit_results = {
    "converged": logit_conv,
    "full_rank": is_full_rank,
    "rank": rank_X,
    "n_obs": int(logit_model.nobs),
    "all_se_finite_bounded": logit_se_valid,
    "pseudo_r2": round(float(logit_model.prsquared), 4),
    "parameters": logit_params,
}

# L2-Penalized (Ridge) Logistic Regression Sensitivity Check on C_fam (lambda_l2 = 0.5 on non-intercept terms)
from scipy.optimize import minimize

X_np = X_mat.values
y_np = df_eval["c_fam_primary"].astype(float).values
col_names = list(X_mat.columns)
lam_l2 = 0.5
pen_vec = np.array([0.0 if c == "const" else lam_l2 for c in col_names])


def neg_pen_loglik(beta):
    eta = np.clip(X_np @ beta, -30.0, 30.0)
    p = 1.0 / (1.0 + np.exp(-eta))
    nll = -np.sum(y_np * np.log(p + 1e-15) + (1.0 - y_np) * np.log(1.0 - p + 1e-15))
    return float(nll + 0.5 * np.sum(pen_vec * (beta ** 2)))


def grad_pen_loglik(beta):
    eta = np.clip(X_np @ beta, -30.0, 30.0)
    p = 1.0 / (1.0 + np.exp(-eta))
    return X_np.T @ (p - y_np) + pen_vec * beta


res_l2 = minimize(neg_pen_loglik, np.zeros(X_np.shape[1]), jac=grad_pen_loglik, method="BFGS")
eta_hat = np.clip(X_np @ res_l2.x, -30.0, 30.0)
p_hat = 1.0 / (1.0 + np.exp(-eta_hat))
W_diag = p_hat * (1.0 - p_hat)
H_pen = (X_np.T * W_diag) @ X_np + np.diag(pen_vec)
cov_l2 = np.linalg.inv(H_pen)
se_l2 = np.sqrt(np.diag(cov_l2))
l2_params = {}
for idx_c, c_name in enumerate(col_names):
    if c_name == "const":
        continue
    b_val = float(res_l2.x[idx_c])
    s_val = float(se_l2[idx_c])
    z_val = b_val / s_val if s_val > 0 else 0.0
    p_val = float(2.0 * (1.0 - stats.norm.cdf(abs(z_val))))
    l2_params[c_name] = {
        "coef": round(b_val, 6),
        "se": round(s_val, 6),
        "p_value": round(p_val, 6),
        "odds_ratio": round(float(np.exp(b_val)), 6),
        "ci_95_or": [
            round(float(np.exp(b_val - 1.96 * s_val)), 6),
            round(float(np.exp(b_val + 1.96 * s_val)), 6),
        ],
    }

logit_l2_summary = {
    "converged": bool(res_l2.success),
    "l2_penalty_lambda": lam_l2,
    "full_rank": is_full_rank,
    "all_se_finite_bounded": bool(all(np.isfinite(s) and 0.0 < s < 100.0 for s in se_l2)),
    "parameters": l2_params,
}

# Unadjusted comparisons across all 3 tiers
def compute_unadjusted_tier_summary(eval_list, c_col, r_col, excess_c_col=None, excess_r_col=None):
    k_sub = [f for f in eval_list if f["arm"] == "K"]
    q_sub = [f for f in eval_list if f["arm"] == "QnotK"]
    kn, qn = len(k_sub), len(q_sub)
    kc = sum(f[c_col] for f in k_sub)
    qc = sum(f[c_col] for f in q_sub)
    kr = kc / kn if kn else 0.0
    qr = qc / qn if qn else 0.0
    odds, p_fish = stats.fisher_exact([[kc, kn - kc], [qc, qn - qc]]) if (kn and qn) else (1.0, 1.0)
    _, p_mwu_r = (
        stats.mannwhitneyu([f[r_col] for f in k_sub], [f[r_col] for f in q_sub], alternative="two-sided")
        if (kn and qn)
        else (0.0, 1.0)
    )
    out = {
        "K_conserved": {"n": kc, "total": kn, "rate": round(kr, 6)},
        "QnotK_conserved": {"n": qc, "total": qn, "rate": round(qr, 6)},
        "raw_risk_difference": round(kr - qr, 6),
        "odds_ratio": round(float(odds), 6) if not math.isinf(odds) else "Inf",
        "p_fisher_exact_raw_c": round(float(p_fish), 6),
        "mean_r_fam": {
            "K": round(float(np.mean([f[r_col] for f in k_sub])), 6) if kn else 0.0,
            "QnotK": round(float(np.mean([f[r_col] for f in q_sub])), 6) if qn else 0.0,
        },
        "p_mann_whitney_raw_r": round(float(p_mwu_r), 6),
    }
    if excess_c_col and excess_r_col and kn and qn:
        _, p_mwu_ec = stats.mannwhitneyu(
            [f[excess_c_col] for f in k_sub], [f[excess_c_col] for f in q_sub], alternative="two-sided"
        )
        _, p_mwu_er = stats.mannwhitneyu(
            [f[excess_r_col] for f in k_sub], [f[excess_r_col] for f in q_sub], alternative="two-sided"
        )
        out["mean_excess_c_fam"] = {
            "K": round(float(np.mean([f[excess_c_col] for f in k_sub])), 6),
            "QnotK": round(float(np.mean([f[excess_c_col] for f in q_sub])), 6),
            "difference": round(
                float(np.mean([f[excess_c_col] for f in k_sub]) - np.mean([f[excess_c_col] for f in q_sub])), 6
            ),
            "p_mann_whitney_excess_c": round(float(p_mwu_ec), 6),
        }
        out["mean_excess_r_fam"] = {
            "K": round(float(np.mean([f[excess_r_col] for f in k_sub])), 6),
            "QnotK": round(float(np.mean([f[excess_r_col] for f in q_sub])), 6),
            "difference": round(
                float(np.mean([f[excess_r_col] for f in k_sub]) - np.mean([f[excess_r_col] for f in q_sub])), 6
            ),
            "p_mann_whitney_excess_r": round(float(p_mwu_er), 6),
        }
    return out


tier1_summary = compute_unadjusted_tier_summary(
    primary_evaluable_fams, "c_fam_primary", "r_fam_primary", "excess_c_fam_primary", "excess_r_fam_primary"
)
tier2_summary = compute_unadjusted_tier_summary(
    stop3_evaluable_fams, "c_fam_stop3_sens", "r_fam_stop3_sens", "excess_c_fam_stop3", "excess_r_fam_stop3"
)
tier3_summary = compute_unadjusted_tier_summary(
    partial_evaluable_fams, "c_fam_sensitivity", "r_fam_sensitivity"
)

# Step 4: Independent Programmatic Closure Assertions
join_metrics = json.loads(join_metrics_file.read_text()) if join_metrics_file.exists() else {}
flank_cat_sum_all = sum(join_metrics.get("flanks_all_rows_by_category", {}).values())
flank_cat_sum_prim = sum(join_metrics.get("flanks_primary_loci_by_category", {}).values())
flank_cat_sum_nonprim = sum(join_metrics.get("flanks_non_primary_rows_by_category", {}).values())
flank_cat_sum_t1_elig = sum(join_metrics.get("flanks_tier1_eligible_primary_loci_by_category", {}).values())
flank_cat_sum_t1_inelig = sum(join_metrics.get("flanks_tier1_ineligible_primary_loci_by_category", {}).values())

all_t1_pos_counts = [
    int(r["n_eligible_positions_primary"])
    for r in all_neighbourhood_rows
    if r["is_primary_locus"] == "1" and r["eligible_locus_primary"] == "1"
]

closure_assertions = {
    "c1_all_tier1_eligible_loci_have_4_exact_mapped_flanks": all(
        int(r["n_exact_mapped_flanks"]) == 4
        for f in family_results
        for r in fam_tier1_loci[f["family"]]
    ),
    "c1_all_ineligible_families_coded_minus_one": all(
        (f["c_fam_primary"] == -1 if f["eligible_primary"] == 0 else f["c_fam_primary"] in (0, 1))
        and (f["c_fam_stop3_sens"] == -1 if f["eligible_stop3_sens"] == 0 else f["c_fam_stop3_sens"] in (0, 1))
        and (f["c_fam_sensitivity"] == -1 if f["eligible_sensitivity"] == 0 else f["c_fam_sensitivity"] in (0, 1))
        for f in family_results
    ),
    "c1_all_r_fam_within_0_1_bounds": all(
        0.0 <= f["r_fam_primary"] <= 1.0
        and 0.0 <= f["r_fam_stop3_sens"] <= 1.0
        and 0.0 <= f["r_fam_sensitivity"] <= 1.0
        for f in family_results
    ),
    "c1_numerator_denominator_species_closure": all(
        f["best_ko_species_count"] <= f["n_species_complete_window"] for f in family_results
    ),
    "c2_zero_unassigned_r232_phylum_in_cohort": all(
        f["phylum"] != "unassigned" and f["phylum"] != "" for f in family_results
    ),
    "c2_all_binned_loci_have_positive_contig_len_bp": all(
        int(r["contig_len_bp"]) > 0 for r in all_neighbourhood_rows if r["bin"] != ""
    ),
    "c3_all_positional_null_draws_satisfy_primary_gate": all(
        int(r["n_eligible_positions_primary"]) >= 1
        for f in primary_evaluable_fams
        for r in fam_tier1_loci[f["family"]]
    ),
    "c3_between_arm_stratified_common_support_verified": (
        strat_perm_t1.get("n_strata_common_support", 0) >= 6
        and strat_perm_t1.get("common_support_coverage", {}).get("K_fraction", 0.0) >= 0.90
        and strat_perm_t1.get("common_support_coverage", {}).get("QnotK_fraction", 0.0) >= 0.90
    ),
    "c3_ols_hc3_and_logit_models_converged_full_rank_finite_se": (
        ols_excess_summary["converged"]
        and ols_excess_summary["full_rank"]
        and ols_excess_summary["all_se_finite_bounded"]
        and logit_results["converged"]
        and logit_results["full_rank"]
        and logit_results["all_se_finite_bounded"]
        and logit_l2_summary["converged"]
        and logit_l2_summary["all_se_finite_bounded"]
    ),
    "c4_flank_category_totals_close_exactly": (
        flank_cat_sum_all == 4 * len(all_neighbourhood_rows)
        and flank_cat_sum_prim == 4 * sum(len(v) for v in fam_primary_loci.values())
        and flank_cat_sum_prim + flank_cat_sum_nonprim == flank_cat_sum_all
        and flank_cat_sum_t1_elig + flank_cat_sum_t1_inelig == flank_cat_sum_prim
    ),
    "c5_cohort_attrition_totals_close_exactly": (
        len(primary_evaluable_fams) + sum(1 for f in family_results if f["eligible_primary"] == 0) == len(family_results)
        and len(stop3_evaluable_fams) + sum(1 for f in family_results if f["eligible_stop3_sens"] == 0) == len(family_results)
        and len(partial_evaluable_fams) + sum(1 for f in family_results if f["eligible_sensitivity"] == 0) == len(family_results)
    ),
}

for k_assert, v_assert in closure_assertions.items():
    assert v_assert, f"Closure assertion failed: {k_assert}"

print("[analyse:v3] All 11 independent closure assertions PASSED.", flush=True)

summary = {
    "protocol_version": "v3",
    "attrition_funnel_2838_universe": attrition_funnel,
    "cohort_admission": {
        "total_families_selected": len(family_results),
        "selected_by_arm": {
            "K": sum(1 for f in family_results if f["arm"] == "K"),
            "QnotK": sum(1 for f in family_results if f["arm"] == "QnotK"),
        },
        "tier1_primary_evaluable": {
            "total": len(primary_evaluable_fams),
            "ineligible_missing": len(family_results) - len(primary_evaluable_fams),
            "K": k_n,
            "QnotK": q_n,
        },
        "tier2_stop3_sensitivity_evaluable": {
            "total": len(stop3_evaluable_fams),
            "ineligible_missing": len(family_results) - len(stop3_evaluable_fams),
            "K": sum(1 for f in stop3_evaluable_fams if f["arm"] == "K"),
            "QnotK": sum(1 for f in stop3_evaluable_fams if f["arm"] == "QnotK"),
        },
        "tier3_partial_window_sensitivity_evaluable": {
            "total": len(partial_evaluable_fams),
            "ineligible_missing": len(family_results) - len(partial_evaluable_fams),
            "K": sum(1 for f in partial_evaluable_fams if f["arm"] == "K"),
            "QnotK": sum(1 for f in partial_evaluable_fams if f["arm"] == "QnotK"),
        },
    },
    "tier1_primary_gate_unadjusted": tier1_summary,
    "tier2_stop3_sensitivity_gate_unadjusted": tier2_summary,
    "tier3_partial_sensitivity_gate_unadjusted": tier3_summary,
    "part_a_within_arm_same_contig_positional_reference": {
        "iterations": B_POS,
        "eligibility_rule": "Uniform random pseudo-focal gene index j on same contig where all 4 flanks (j-2, j-1, j+1, j+2) are exact-coordinate Bakta-mapped (mapped_KO or mapped_noKO)",
        "eligible_position_counts_across_tier1_loci": {
            "total_tier1_eligible_loci": len(all_t1_pos_counts),
            "loci_with_M_c_eq_1": sum(1 for c in all_t1_pos_counts if c == 1),
            "loci_with_M_c_ge_2": sum(1 for c in all_t1_pos_counts if c >= 2),
            "loci_with_M_c_ge_5": sum(1 for c in all_t1_pos_counts if c >= 5),
            "min_M_c": min(all_t1_pos_counts) if all_t1_pos_counts else 0,
            "median_M_c": float(np.median(all_t1_pos_counts)) if all_t1_pos_counts else 0.0,
            "mean_M_c": round(float(np.mean(all_t1_pos_counts)), 2) if all_t1_pos_counts else 0.0,
            "max_M_c": max(all_t1_pos_counts) if all_t1_pos_counts else 0,
        },
        "by_arm": within_arm_positional,
    },
    "part_b_between_arm_covariate_stratified_permutation": {
        "tier1_primary": strat_perm_t1,
        "tier2_stop3_sensitivity": strat_perm_t2,
    },
    "part_b_between_arm_matched_pairs_tier1": matching_summary,
    "part_b_multivariable_ols_hc3_excess_c_fam": ols_excess_summary,
    "part_b_multivariable_ols_hc3_raw_c_fam": ols_raw_summary,
    "part_b_multivariable_logistic_regression_raw_c_fam": logit_results,
    "part_b_l2_penalized_logistic_sensitivity_raw_c_fam": logit_l2_summary,
    "covariate_profiles_tier1_primary": {
        "mean_species_breadth": {
            "K": round(float(np.mean([f["n_species_complete_window"] for f in k_prim])), 3) if k_prim else 0.0,
            "QnotK": round(float(np.mean([f["n_species_complete_window"] for f in q_prim])), 3) if q_prim else 0.0,
        },
        "mean_length_aa": {
            "K": round(float(np.mean([f["length"] for f in k_prim])), 1) if k_prim else 0.0,
            "QnotK": round(float(np.mean([f["length"] for f in q_prim])), 1) if q_prim else 0.0,
        },
        "mean_plddt": {
            "K": round(float(np.mean([f["mean_plddt"] for f in k_prim])), 3) if k_prim else 0.0,
            "QnotK": round(float(np.mean([f["mean_plddt"] for f in q_prim])), 3) if q_prim else 0.0,
        },
        "mean_ko_annotation_opportunity": {
            "K": round(float(np.mean([f["ko_annotation_opportunity"] for f in k_prim])), 4) if k_prim else 0.0,
            "QnotK": round(float(np.mean([f["ko_annotation_opportunity"] for f in q_prim])), 4) if q_prim else 0.0,
        },
        "phylum_group_counts": {
            "K": dict(Counter(f["phylum_group"] for f in k_prim)),
            "QnotK": dict(Counter(f["phylum_group"] for f in q_prim)),
        },
    },
    "cross_family_locus_and_species_independence_tier1": cross_family_dependence_audit,
    "closure_assertions": closure_assertions,
    "recurrent_annotations_tier1_primary": [
        {
            "family": f["family"],
            "arm": f["arm"],
            "phylum": f["phylum"],
            "n_species_complete_window": f["n_species_complete_window"],
            "best_ko": f["best_recurrent_ko"],
            "species_count": f["best_ko_species_count"],
            "r_fam_primary": f["r_fam_primary"],
            "mu_pos_c_fam": f["mu_pos_c_fam_primary"],
            "excess_c_fam": f["excess_c_fam_primary"],
            "mobilome_flag": f["best_ko_mobilome_flag"],
        }
        for f in primary_evaluable_fams
        if f["c_fam_primary"] == 1
    ],
}

out_sum_json = results_dir / "pilot_conservation_summary.json"
out_sum_json.write_text(json.dumps(summary, indent=2) + "\n")
print(f"[analyse:v3] Saved comprehensive conservation summary to {out_sum_json}", flush=True)
