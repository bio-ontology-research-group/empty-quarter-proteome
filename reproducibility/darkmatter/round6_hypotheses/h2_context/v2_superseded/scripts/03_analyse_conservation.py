#!/usr/bin/env python3
"""03_analyse_conservation.py - Rigorous family-level conservation analysis (v2 protocol).

Implements:
1. Strict eligibility gates:
   - Primary: >= 2 distinct species with complete four-flank windows (ineligible -> missing, not negative).
   - Sensitivity: >= 2 distinct species with >= 1 callable flank.
2. Multiplicity control: evaluates only designated primary locus per species.
3. Family-level estimands: C_fam (binary KO recurrence) and R_fam (continuous recurrence score).
4. Host/contig-preserving random focal-position permutation null (B=1000):
   - Uniform pseudo-focal sampling on the same contig.
   - Flanks inherit true contig geometry and annotations.
   - Monte Carlo plus-one corrected p-value.
5. Multivariable logistic regression adjusting for species breadth, lineage, length, pLDDT, and annotation opportunity.
6. Mobilome / phage overlap audit.
Emits:
- pilot_family_conservation.tsv
- pilot_conservation_summary.json
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
neighbourhoods_file = results_dir / "pilot_neighbourhoods.tsv"
contig_genes_file = results_dir / "pilot_contig_genes.tsv"

for p in [cohort_file, neighbourhoods_file, contig_genes_file]:
    if not p.exists():
        sys.exit(f"Error: {p} does not exist. Run earlier steps first.")

# Load cohort metadata
cohort = {}
with open(cohort_file) as f:
    for row in csv.DictReader(f, delimiter="\t"):
        cohort[row["family"]] = row

# Load neighbourhood rows (filter to primary loci only)
fam_members_all = defaultdict(list)
fam_primary_loci = defaultdict(list)
with open(neighbourhoods_file) as f:
    for row in csv.DictReader(f, delimiter="\t"):
        fam = row["family"]
        fam_members_all[fam].append(row)
        if row["is_primary_locus"] == "1" and row["species"]:
            fam_primary_loci[fam].append(row)

# Load contig genes for host/contig-preserving null
contig_genes = defaultdict(list)
with open(contig_genes_file) as f:
    for row in csv.DictReader(f, delimiter="\t"):
        contig_genes[row["contig"]].append(row)

print(f"[analyse] Loaded {len(cohort)} cohort families, {len(contig_genes)} member contigs")

# Mobilome / phage terms for companion overlap bounding
MOBILOME_TERMS = {
    "transposase", "integrase", "recombinase", "insertion element",
    "phage", "prophage", "terminase", "capsid", "tail", "portal", "plasmid"
}

def is_mobilome(product_str, desc_str):
    text = (product_str + " " + desc_str).lower()
    return any(term in text for term in MOBILOME_TERMS)

# Evaluate family conservation under Primary (complete window) and Sensitivity gates
family_results = []
primary_evaluable_fams = []
sensitivity_evaluable_fams = []

for fam, meta in sorted(cohort.items()):
    prim_loci = fam_primary_loci.get(fam, [])
    arm = meta["arm"]
    
    # Track loci per species
    complete_window_species = set()
    partial_window_species = set()
    
    # Store KOs by species for both gates
    ko_by_sp_primary = defaultdict(set)
    ko_by_sp_sensitivity = defaultdict(set)
    prod_by_sp_primary = defaultdict(set)
    
    total_flanks = 0
    edge_truncated_flanks = 0
    mapped_ko_flanks = 0
    mapped_noko_flanks = 0
    coord_mismatch_flanks = 0
    source_unavail_flanks = 0
    same_strand_flanks = 0
    
    for r in prim_loci:
        sp = r["species"]
        is_complete = True
        has_any_callable = False
        
        for flank_key in ["m2", "m1", "p1", "p2"]:
            total_flanks += 1
            audit = r[f"{flank_key}_audit_category"]
            if audit == "contig_edge_truncated":
                edge_truncated_flanks += 1
                is_complete = False
            elif audit == "source_file_unavailable":
                source_unavail_flanks += 1
                has_any_callable = True
            elif audit == "coordinate_mismatch":
                coord_mismatch_flanks += 1
                has_any_callable = True
            elif audit == "mapped_noKO":
                mapped_noko_flanks += 1
                has_any_callable = True
            elif audit == "mapped_KO":
                mapped_ko_flanks += 1
                has_any_callable = True
                
            if r[f"{flank_key}_same_strand"] == "1":
                same_strand_flanks += 1
                
            ko_str = r[f"{flank_key}_ko"]
            prod = r[f"{flank_key}_product"]
            if ko_str:
                for k in ko_str.split(","):
                    if k:
                        ko_by_sp_sensitivity[k].add(sp)
                        if r[f"{flank_key}_truncated"] == "0":
                            ko_by_sp_primary[k].add(sp)
            if prod and prod not in ["hypothetical protein", "CONTIG_EDGE", "UNAVAILABLE", "COORD_MISMATCH"]:
                prod_by_sp_primary[prod].add(sp)
                
        if is_complete:
            complete_window_species.add(sp)
        if has_any_callable:
            partial_window_species.add(sp)
            
    n_sp_complete = len(complete_window_species)
    n_sp_partial = len(partial_window_species)
    
    eligible_primary = 1 if n_sp_complete >= 2 else 0
    eligible_sensitivity = 1 if n_sp_partial >= 2 else 0
    
    # Primary gate recurrence
    recurrent_kos_primary = {k: len(sps) for k, sps in ko_by_sp_primary.items() if len(sps) >= 2}
    c_fam_primary = 1 if (eligible_primary and len(recurrent_kos_primary) > 0) else (0 if eligible_primary else -1)
    
    best_ko_primary = ""
    best_ko_sp_count = 0
    best_ko_mobilome = 0
    if recurrent_kos_primary and eligible_primary:
        best_ko_primary, best_ko_sp_count = max(recurrent_kos_primary.items(), key=lambda x: x[1])
        for r in prim_loci:
            for flank_key in ["m2", "m1", "p1", "p2"]:
                if best_ko_primary in r[f"{flank_key}_ko"].split(","):
                    if is_mobilome(r[f"{flank_key}_product"], r[f"{flank_key}_desc"]):
                        best_ko_mobilome = 1
                        break
                        
    r_fam_primary = (best_ko_sp_count / n_sp_complete) if (eligible_primary and n_sp_complete > 0) else 0.0
    
    # Sensitivity gate recurrence
    recurrent_kos_sens = {k: len(sps) for k, sps in ko_by_sp_sensitivity.items() if len(sps) >= 2}
    c_fam_sensitivity = 1 if (eligible_sensitivity and len(recurrent_kos_sens) > 0) else (0 if eligible_sensitivity else -1)
    
    callable_flanks = total_flanks - edge_truncated_flanks
    ko_annot_rate = mapped_ko_flanks / callable_flanks if callable_flanks > 0 else 0.0
    
    fam_row = {
        "family": fam,
        "chosen_rep": meta["chosen_rep"],
        "arm": arm,
        "phylum": meta["phylum"],
        "length": int(meta["length"]),
        "mean_plddt": float(meta["mean_plddt"]),
        "n_members": int(meta["n_members"]),
        "n_primary_loci": len(prim_loci),
        "n_species_complete_window": n_sp_complete,
        "n_species_partial_window": n_sp_partial,
        "eligible_primary": eligible_primary,
        "eligible_sensitivity": eligible_sensitivity,
        "c_fam_primary": c_fam_primary,
        "c_fam_sensitivity": c_fam_sensitivity,
        "r_fam_primary": round(r_fam_primary, 4),
        "best_recurrent_ko": best_ko_primary,
        "best_ko_species_count": best_ko_sp_count,
        "best_ko_mobilome_flag": best_ko_mobilome,
        "flanks_total": total_flanks,
        "flanks_edge_truncated": edge_truncated_flanks,
        "flanks_mapped_KO": mapped_ko_flanks,
        "flanks_mapped_noKO": mapped_noko_flanks,
        "flanks_coord_mismatch": coord_mismatch_flanks,
        "flanks_source_unavail": source_unavail_flanks,
        "ko_annotation_opportunity": round(ko_annot_rate, 4),
        "same_strand_fraction": round(same_strand_flanks / callable_flanks, 4) if callable_flanks > 0 else 0.0
    }
    family_results.append(fam_row)
    if eligible_primary:
        primary_evaluable_fams.append(fam_row)
    if eligible_sensitivity:
        sensitivity_evaluable_fams.append(fam_row)

# Save per-family results
out_fam_tsv = results_dir / "pilot_family_conservation.tsv"
with open(out_fam_tsv, "w") as f:
    w = csv.DictWriter(f, fieldnames=list(family_results[0].keys()), delimiter="\t")
    w.writeheader()
    w.writerows(family_results)

print(f"[analyse] Saved family conservation metrics to {out_fam_tsv}")
print(f"[analyse] Primary evaluable families (>=2 complete species): {len(primary_evaluable_fams)} / {len(family_results)}")
print(f"[analyse] Sensitivity evaluable families (>=2 partial species): {len(sensitivity_evaluable_fams)} / {len(family_results)}")

# Statistical comparison for Primary Gate
k_prim = [f for f in primary_evaluable_fams if f["arm"] == "K"]
q_prim = [f for f in primary_evaluable_fams if f["arm"] == "QnotK"]

k_prim_cons = sum(f["c_fam_primary"] for f in k_prim)
q_prim_cons = sum(f["c_fam_primary"] for f in q_prim)
k_n = len(k_prim)
q_n = len(q_prim)

k_rate = k_prim_cons / k_n if k_n else 0.0
q_rate = q_prim_cons / q_n if q_n else 0.0
diff_rate = k_rate - q_rate

contingency_prim = [
    [k_prim_cons, k_n - k_prim_cons],
    [q_prim_cons, q_n - q_prim_cons]
]
odds_ratio_prim, p_fisher_prim = stats.fisher_exact(contingency_prim)

# Mann-Whitney on recurrence scores
u_stat, p_mannwhitney_prim = stats.mannwhitneyu(
    [f["r_fam_primary"] for f in k_prim] if k_prim else [0],
    [f["r_fam_primary"] for f in q_prim] if q_prim else [0],
    alternative="two-sided"
)

# Host/contig-preserving random focal-position permutation null
# Preserves contig architecture, gene count, edge truncation, and real KOs on each contig
print("[analyse] Running host/contig-preserving permutation null (B=1000)...", flush=True)
rng = random.Random(20261013)
B = 1000
null_diffs = []

# Pre-index contig CDSs by contig
contig_cds_lookup = {c: genes for c, genes in contig_genes.items()}

for b_idx in range(B):
    b_k_cons = 0
    b_q_cons = 0
    
    for f in primary_evaluable_fams:
        fam = f["family"]
        prim_loci = fam_primary_loci[fam]
        
        null_ko_by_sp = defaultdict(set)
        
        for r in prim_loci:
            sp = r["species"]
            c = r["contig"]
            c_cds = contig_cds_lookup.get(c, [])
            n_cds = len(c_cds)
            if n_cds == 0:
                continue
                
            # Sample random pseudo-focal locus uniformly on the SAME contig
            pseudo_idx = rng.randint(0, n_cds - 1)
            
            # Extract flanks under the same window definition
            for offset in [-2, -1, 1, 2]:
                target_idx = pseudo_idx + offset
                if 0 <= target_idx < n_cds:
                    target_cds = c_cds[target_idx]
                    ko_val = target_cds["ko"]
                    if ko_val:
                        for pk in ko_val.split(","):
                            if pk:
                                null_ko_by_sp[pk].add(sp)
                                
        has_recurrent = 1 if any(len(sps) >= 2 for sps in null_ko_by_sp.values()) else 0
        if f["arm"] == "K":
            b_k_cons += has_recurrent
        else:
            b_q_cons += has_recurrent
            
    b_k_rate = b_k_cons / k_n if k_n else 0.0
    b_q_rate = b_q_cons / q_n if q_n else 0.0
    null_diffs.append(b_k_rate - b_q_rate)

# Monte Carlo plus-one corrected p-value
p_perm_plus_one = (1.0 + sum(1 for d in null_diffs if abs(d) >= abs(diff_rate))) / (1.0 + B)
null_diff_mean = float(np.mean(null_diffs))
null_diff_std = float(np.std(null_diffs))

print(f"[analyse] Permutation null completed: mean={null_diff_mean:.4f}, std={null_diff_std:.4f}, p_plus_one={p_perm_plus_one:.4f}")

# Multivariable Logistic Regression
# Models logit(P(C_fam = 1)) ~ arm_K + n_species + n_members + length + mean_plddt + ko_opportunity + phylum
logit_results = {}
try:
    df_eval = pd.DataFrame(primary_evaluable_fams)
    df_eval["is_k"] = (df_eval["arm"] == "K").astype(int)
    
    # Feature matrix
    X_cols = ["is_k", "n_species_complete_window", "n_members", "length", "mean_plddt", "ko_annotation_opportunity"]
    
    # Add dominant phylum indicator if variation exists
    phylum_counts = df_eval["phylum"].value_counts()
    frequent_phyla = phylum_counts[phylum_counts >= 5].index.tolist()
    for ph in frequent_phyla:
        col_name = f"phylum_{ph}"
        df_eval[col_name] = (df_eval["phylum"] == ph).astype(int)
        X_cols.append(col_name)
        
    y = df_eval["c_fam_primary"]
    X = sm.add_constant(df_eval[X_cols])
    
    if y.nunique() < 2:
        logit_results = {"converged": False, "error": "Outcome y has no variation in primary evaluable subset"}
    else:
        try:
            model = sm.Logit(y, X).fit(disp=False, method="newton", maxiter=100)
        except Exception:
            model = sm.Logit(y, X).fit_regularized(disp=False, alpha=1.0)
            
        logit_summary = {}
        for param in X_cols:
            coef = float(model.params[param]) if param in model.params else 0.0
            se = float(model.bse[param]) if hasattr(model, "bse") and param in model.bse else 0.0
            pval = float(model.pvalues[param]) if hasattr(model, "pvalues") and param in model.pvalues else 1.0
            odds = float(np.exp(coef)) if abs(coef) < 700 else "Inf"
            ci_lower = float(np.exp(coef - 1.96 * se)) if se > 0 and abs(coef - 1.96 * se) < 700 else 0.0
            ci_upper = float(np.exp(coef + 1.96 * se)) if se > 0 and abs(coef + 1.96 * se) < 700 else "Inf"
            logit_summary[param] = {
                "coef": round(coef, 4),
                "se": round(se, 4) if se > 0 else "N/A",
                "p_value": round(pval, 4),
                "odds_ratio": round(odds, 4) if isinstance(odds, float) else odds,
                "ci_95": [round(ci_lower, 4) if isinstance(ci_lower, float) else ci_lower, round(ci_upper, 4) if isinstance(ci_upper, float) else ci_upper]
            }
        logit_results = {
            "converged": bool(getattr(model, "mle_retvals", {}).get("converged", True)),
            "pseudo_r2": round(float(getattr(model, "prsquared", 0.0)), 4),
            "parameters": logit_summary
        }
except Exception as e:
    logit_results = {"converged": False, "error": str(e)}

# Sensitivity Gate Comparison
k_sens = [f for f in sensitivity_evaluable_fams if f["arm"] == "K"]
q_sens = [f for f in sensitivity_evaluable_fams if f["arm"] == "QnotK"]
k_sens_cons = sum(f["c_fam_sensitivity"] for f in k_sens)
q_sens_cons = sum(f["c_fam_sensitivity"] for f in q_sens)
k_sens_n = len(k_sens)
q_sens_n = len(q_sens)
k_sens_rate = k_sens_cons / k_sens_n if k_sens_n else 0.0
q_sens_rate = q_sens_cons / q_sens_n if q_sens_n else 0.0

contingency_sens = [
    [k_sens_cons, k_sens_n - k_sens_cons],
    [q_sens_cons, q_sens_n - q_sens_cons]
]
odds_ratio_sens, p_fisher_sens = stats.fisher_exact(contingency_sens)

# Compile comprehensive summary JSON
summary = {
    "cohort_admission": {
        "total_families_selected": len(cohort),
        "total_primary_evaluable": len(primary_evaluable_fams),
        "primary_attrition_ineligible": len(cohort) - len(primary_evaluable_fams),
        "primary_evaluable_by_arm": {"K": k_n, "QnotK": q_n},
        "total_sensitivity_evaluable": len(sensitivity_evaluable_fams),
        "sensitivity_evaluable_by_arm": {"K": k_sens_n, "QnotK": q_sens_n}
    },
    "primary_gate_results": {
        "definition": ">=2 distinct species with complete four-flank observable windows",
        "K_conserved": {"n": k_prim_cons, "total": k_n, "rate": round(k_rate, 4)},
        "QnotK_conserved": {"n": q_prim_cons, "total": q_n, "rate": round(q_rate, 4)},
        "risk_difference": round(diff_rate, 4),
        "odds_ratio": round(odds_ratio_prim, 4) if not math.isinf(odds_ratio_prim) else "Inf",
        "p_fisher_exact": round(p_fisher_prim, 6),
        "p_mann_whitney_r_fam": round(p_mannwhitney_prim, 6)
    },
    "sensitivity_gate_results": {
        "definition": ">=2 distinct species with >=1 callable flank",
        "K_conserved": {"n": k_sens_cons, "total": k_sens_n, "rate": round(k_sens_rate, 4)},
        "QnotK_conserved": {"n": q_sens_cons, "total": q_sens_n, "rate": round(q_sens_rate, 4)},
        "risk_difference": round(k_sens_rate - q_sens_rate, 4),
        "odds_ratio": round(odds_ratio_sens, 4) if not math.isinf(odds_ratio_sens) else "Inf",
        "p_fisher_exact": round(p_fisher_sens, 6)
    },
    "host_contig_preserving_permutation_null": {
        "iterations": B,
        "sampling_unit": "uniform random pseudo-focal locus on same contig",
        "null_diff_mean": round(null_diff_mean, 4),
        "null_diff_std": round(null_diff_std, 4),
        "p_permutation_plus_one": round(p_perm_plus_one, 4)
    },
    "multivariable_logistic_regression": logit_results,
    "covariate_profiles_primary": {
        "mean_species_breadth": {
            "K": round(float(np.mean([f["n_species_complete_window"] for f in k_prim])), 2) if k_prim else 0.0,
            "QnotK": round(float(np.mean([f["n_species_complete_window"] for f in q_prim])), 2) if q_prim else 0.0
        },
        "mean_length_aa": {
            "K": round(float(np.mean([f["length"] for f in k_prim])), 1) if k_prim else 0.0,
            "QnotK": round(float(np.mean([f["length"] for f in q_prim])), 1) if q_prim else 0.0
        },
        "mean_plddt": {
            "K": round(float(np.mean([f["mean_plddt"] for f in k_prim])), 3) if k_prim else 0.0,
            "QnotK": round(float(np.mean([f["mean_plddt"] for f in q_prim])), 3) if q_prim else 0.0
        },
        "mean_ko_annotation_opportunity": {
            "K": round(float(np.mean([f["ko_annotation_opportunity"] for f in k_prim])), 3) if k_prim else 0.0,
            "QnotK": round(float(np.mean([f["ko_annotation_opportunity"] for f in q_prim])), 3) if q_prim else 0.0
        }
    },
    "recurrent_annotations_primary": [
        {
            "family": f["family"],
            "arm": f["arm"],
            "best_ko": f["best_recurrent_ko"],
            "species_count": f["best_ko_species_count"],
            "mobilome_flag": f["best_ko_mobilome_flag"]
        }
        for f in primary_evaluable_fams if f["c_fam_primary"] == 1
    ]
}

out_sum_json = results_dir / "pilot_conservation_summary.json"
out_sum_json.write_text(json.dumps(summary, indent=2) + "\n")
print(f"[analyse] Saved comprehensive conservation summary to {out_sum_json}")
