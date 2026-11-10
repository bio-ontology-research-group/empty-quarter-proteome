#!/usr/bin/env python3
"""Stage 3 Collector of `exp1_metabolic`: Aggregates the 8 shards across all 484 deduplicated genomes
(237 EQ MAGs + 247 external comparators), validates exact ID coverage and shard completeness,
enforces NA on invalid/infeasible models, performs descriptive unadjusted and lineage-matched comparisons
under both canonical (`refined_representatives.tsv`) and legacy taxonomy, and emits supported candidate
deficiencies for Stage 4 sequence/synteny verification (`EXP1-ACCEPTANCE-AND-COLLECTOR-REVIEW.txt` Item 3 & 4).
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy import stats

CODE_DIR = Path(__file__).resolve().parent
if str(CODE_DIR) not in sys.path:
    sys.path.insert(0, str(CODE_DIR))

from exp1_sbml_highs_engine import (
    BIO_B2_OMIT_7_KEEP_F420,
    SCEN_S4_OXIC_DESERT,
)

REVALIDATED_FEASIBLE_STATUSES = {
    "optimal",
    "optimal_after_cycle_qc_filtering",
    "optimal_draft_no_gapfill_needed",
    "heuristic_incumbent_bounded",
    "post_pruning_cost_changed_not_proven_optimal",
}

MISSING_TAXONOMY_PLACEHOLDERS = {
    "", "none", "nan", "f__", "o__", "g__", "unclassified", "na", "null"
}


def compute_mannwhitney_and_spearman(
    df: pd.DataFrame,
    metric_col: str,
) -> dict[str, Any]:
    """Compute descriptive unadjusted EQ vs External Mann-Whitney U, exploratory completeness correlation,
    and HQ-restricted sensitivity analysis (`EXP1-ACCEPTANCE-AND-COLLECTOR-REVIEW.txt` Item 3).
    Note: Unadjusted Mann-Whitney and completeness Spearman are descriptive/unadjusted; HQ restriction
    is sensitivity analysis only. Invalid/infeasible models are strictly excluded from numerical medians/MW.
    """
    sub = df[df[metric_col].notna()].copy()
    eq_sub = sub[sub["cohort"] == "EQ_desert"]
    ext_sub = sub[sub["cohort"] == "external_comparator"]
    eq_vals = eq_sub[metric_col].to_numpy(dtype=float)
    ext_vals = ext_sub[metric_col].to_numpy(dtype=float)

    # Denominator and status accounting across all evaluated models in group
    all_eq = df[df["cohort"] == "EQ_desert"]
    all_ext = df[df["cohort"] == "external_comparator"]
    n_eq_total = len(all_eq)
    n_ext_total = len(all_ext)
    eq_status_counts = all_eq["solver_status"].value_counts().to_dict() if not all_eq.empty else {}
    ext_status_counts = all_ext["solver_status"].value_counts().to_dict() if not all_ext.empty else {}

    if len(eq_vals) > 0 and len(ext_vals) > 0:
        u_stat, p_val = stats.mannwhitneyu(eq_vals, ext_vals, alternative="two-sided")
        eq_med = float(np.median(eq_vals))
        ext_med = float(np.median(ext_vals))
    else:
        u_stat, p_val, eq_med, ext_med = math.nan, math.nan, math.nan, math.nan

    # Descriptive Spearman correlation with completeness within EQ_desert
    eq_comp = sub[(sub["cohort"] == "EQ_desert") & (sub["completeness"].notna())]
    if len(eq_comp) >= 5:
        rho_comp, p_comp = stats.spearmanr(eq_comp["completeness"], eq_comp[metric_col])
    else:
        rho_comp, p_comp = math.nan, math.nan

    # HQ sensitivity analysis (completeness >= 90%, contamination <= 5%)
    all_hq_eq = all_eq[(all_eq["completeness"] >= 90.0) & (all_eq["contamination"] <= 5.0)]
    all_hq_ext = all_ext[(all_ext["completeness"] >= 90.0) & (all_ext["contamination"] <= 5.0)]
    hq_sub = sub[(sub["completeness"] >= 90.0) & (sub["contamination"] <= 5.0)]
    hq_eq_vals = hq_sub[hq_sub["cohort"] == "EQ_desert"][metric_col].to_numpy(dtype=float)
    hq_ext_vals = hq_sub[hq_sub["cohort"] == "external_comparator"][metric_col].to_numpy(dtype=float)
    if len(hq_eq_vals) > 0 and len(hq_ext_vals) > 0:
        _, hq_p = stats.mannwhitneyu(hq_eq_vals, hq_ext_vals, alternative="two-sided")
        hq_eq_med = float(np.median(hq_eq_vals))
        hq_ext_med = float(np.median(hq_ext_vals))
    else:
        hq_p, hq_eq_med, hq_ext_med = math.nan, math.nan, math.nan

    return {
        "metric": metric_col,
        "comparison_type": "descriptive_unadjusted_with_sensitivity",
        "n_eq_total_evaluated": n_eq_total,
        "n_eq_valid_accepted_models": len(eq_vals),
        "n_external_total_evaluated": n_ext_total,
        "n_external_valid_accepted_models": len(ext_vals),
        "eq_median_unadjusted": eq_med,
        "external_median_unadjusted": ext_med,
        "unadjusted_mannwhitney_u": float(u_stat) if not math.isnan(u_stat) else math.nan,
        "unadjusted_mannwhitney_p_exploratory": float(p_val) if not math.isnan(p_val) else math.nan,
        "eq_completeness_spearman_rho_descriptive": float(rho_comp) if not math.isnan(rho_comp) else math.nan,
        "eq_completeness_spearman_p_exploratory": float(p_comp) if not math.isnan(p_comp) else math.nan,
        "n_hq_sensitivity_eq_total": len(all_hq_eq),
        "n_hq_sensitivity_eq_valid": len(hq_eq_vals),
        "n_hq_sensitivity_external_total": len(all_hq_ext),
        "n_hq_sensitivity_external_valid": len(hq_ext_vals),
        "hq_sensitivity_eq_median": hq_eq_med,
        "hq_sensitivity_external_median": hq_ext_med,
        "hq_sensitivity_mannwhitney_p_exploratory": float(hq_p) if not math.isnan(hq_p) else math.nan,
        "eq_status_counts_json": json.dumps(eq_status_counts),
        "external_status_counts_json": json.dumps(ext_status_counts),
    }


def main() -> None:
    ap = argparse.ArgumentParser(description="Stage 3 Collector: Aggregate 484-genome shards and emit Stage 4 targets.")
    ap.add_argument("--inventory-484-tsv", type=Path, required=True)
    ap.add_argument("--shards-dir", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--num-tasks", type=int, default=8)
    args = ap.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    inv_df = pd.read_csv(args.inventory_484_tsv, sep="\t")

    summary_columns = [
        "genome", "cohort", "canonical_order", "canonical_family", "canonical_genus",
        "completeness", "contamination", "scenario_id", "redox_regime", "biomass_recipe_id",
        "solver_status", "optimizer_status", "raw_solver_status", "enumeration_completeness",
        "pool_exhaustion_type", "cost_window_anchor_type", "draft_already_grows", "draft_growth_rate",
        "min_total_cost", "max_ensemble_cost", "first_solve_mip_gap", "first_solve_dual_bound",
        "n_solutions_attempted", "n_optimization_feasible_solutions", "n_cycle_qc_rejected_solutions",
        "n_distinct_valid_solutions", "fewer_than_target_reason", "min_added_rxns_count",
        "max_added_rxns_count", "n_rxns_observed_in_all_enumerated_models", "n_rxns_essential_in_declared_supernetwork"
    ]
    rf_columns = [
        "genome", "cohort", "canonical_order", "canonical_family", "canonical_genus",
        "completeness", "contamination", "scenario_id", "rxn_id", "base_id", "rxn_name",
        "is_transport", "candidate_weight", "prior_evidence", "freq_B0_baseline",
        "freq_B2_omit_7_keep_f420", "freq_difference_B0_minus_B2",
        "biomass_sensitivity_classification", "stability_tier_B0", "stability_tier_B2",
        "supernetwork_ko_necessity_B0", "supernetwork_ko_necessity_B2"
    ]
    wkf_columns = [
        "genome", "cohort", "scenario_id", "biomass_recipe_id", "record_type", "target_id",
        "target_name", "essential_in_native", "essential_in_desert", "growth_with_target",
        "growth_without_target", "fva_min_flux", "fva_max_flux"
    ]

    sum_dfs: list[pd.DataFrame] = []
    ledger_dfs: list[pd.DataFrame] = []
    rf_dfs: list[pd.DataFrame] = []
    wkf_dfs: list[pd.DataFrame] = []
    supp_dfs: list[pd.DataFrame] = []
    shard_receipts: list[dict[str, Any]] = []

    for tid in range(args.num_tasks):
        tag = f"shard_{tid:02d}"
        s_p = args.shards_dir / f"exp1_{tag}_summary.tsv"
        l_p = args.shards_dir / f"exp1_{tag}_solution_ledger.tsv"
        r_p = args.shards_dir / f"exp1_{tag}_rxn_freq.tsv"
        w_p = args.shards_dir / f"exp1_{tag}_wkf.tsv"
        u_p = args.shards_dir / f"exp1_{tag}_supplement.tsv"
        if not s_p.exists():
            raise FileNotFoundError(f"Missing shard summary file: {s_p}")
        sdf = pd.read_csv(s_p, sep="\t")
        sum_dfs.append(sdf)
        n_ledger = 0
        if l_p.exists() and l_p.stat().st_size > 10:
            ldf = pd.read_csv(l_p, sep="\t")
            ledger_dfs.append(ldf)
            n_ledger = len(ldf)
        n_rf = 0
        if r_p.exists() and r_p.stat().st_size > 10:
            rdf = pd.read_csv(r_p, sep="\t")
            rf_dfs.append(rdf)
            n_rf = len(rdf)
        n_wkf = 0
        if w_p.exists() and w_p.stat().st_size > 10:
            wdf = pd.read_csv(w_p, sep="\t")
            wkf_dfs.append(wdf)
            n_wkf = len(wdf)
        n_supp = 0
        if u_p.exists() and u_p.stat().st_size > 10:
            udf = pd.read_csv(u_p, sep="\t")
            supp_dfs.append(udf)
            n_supp = len(udf)

        shard_receipts.append(
            {
                "shard_id": tid,
                "n_summary_rows": len(sdf),
                "n_unique_genomes": int(sdf["genome"].nunique()) if not sdf.empty else 0,
                "n_solution_ledger_rows": n_ledger,
                "n_rxn_freq_rows": n_rf,
                "n_wkf_rows": n_wkf,
                "n_supplement_rows": n_supp,
            }
        )

    # Schema-valid concatenation
    full_summary_df = (
        pd.concat(sum_dfs, ignore_index=True) if sum_dfs else pd.DataFrame(columns=summary_columns)
    ).sort_values(["cohort", "genome", "scenario_id", "biomass_recipe_id"]).reset_index(drop=True)

    full_rf_df = (
        pd.concat(rf_dfs, ignore_index=True) if rf_dfs else pd.DataFrame(columns=rf_columns)
    ).sort_values(["cohort", "genome", "scenario_id", "rxn_id"]).reset_index(drop=True)

    full_wkf_df = (
        pd.concat(wkf_dfs, ignore_index=True) if wkf_dfs else pd.DataFrame(columns=wkf_columns)
    ).sort_values(["cohort", "genome", "scenario_id", "biomass_recipe_id", "record_type", "target_id"]).reset_index(drop=True)

    full_supp_df = (
        pd.concat(supp_dfs, ignore_index=True) if supp_dfs else pd.DataFrame()
    )

    full_ledger_df = (
        pd.concat(ledger_dfs, ignore_index=True) if ledger_dfs else pd.DataFrame()
    )
    if not full_ledger_df.empty:
        full_ledger_df = full_ledger_df.sort_values(
            ["cohort", "genome", "scenario_id", "biomass_recipe_id", "solution_index"]
        ).reset_index(drop=True)
        full_ledger_df.to_csv(args.out_dir / "exp1_full484_per_solution_ledger.tsv", sep="\t", index=False)

        # Fail-closed ledger correspondence audit (EXP1-LEDGER-FINAL-ACCEPTANCE.txt)
        acc_ledger = full_ledger_df[full_ledger_df["validity_disposition"] == "accepted_feasible"]
        acc_counts = acc_ledger.groupby(["genome", "scenario_id", "biomass_recipe_id"]).size().to_dict() if not acc_ledger.empty else {}
        for _, srow in full_summary_df.iterrows():
            k = (str(srow["genome"]), str(srow["scenario_id"]), str(srow["biomass_recipe_id"]))
            exp_count = int(srow["n_distinct_valid_solutions"]) if srow["solver_status"] in REVALIDATED_FEASIBLE_STATUSES and pd.notna(srow["n_distinct_valid_solutions"]) else 0
            act_count = acc_counts.get(k, 0)
            if exp_count != act_count:
                raise RuntimeError(f"Full484 ledger count mismatch for {k}: summary claims {exp_count}, ledger has {act_count}")

    # 1. Exact ID Coverage Verification (`EXP1-ACCEPTANCE-AND-COLLECTOR-REVIEW.txt` Item 3)
    expected_eligible_genomes = set(
        inv_df[(inv_df["has_draft_sbml"] == True) & (inv_df["qc_error"].fillna("") == "")]["genome"].astype(str)
    )
    collected_genomes = set(full_summary_df["genome"].astype(str)) if not full_summary_df.empty else set()
    missing_genomes = expected_eligible_genomes - collected_genomes
    unexpected_genomes = collected_genomes - expected_eligible_genomes

    if missing_genomes or unexpected_genomes:
        raise RuntimeError(
            f"Stage 3 Collector exact genome ID mismatch! "
            f"Expected {len(expected_eligible_genomes)} eligible genomes, got {len(collected_genomes)}. "
            f"Missing ({len(missing_genomes)}): {sorted(missing_genomes)[:5]}, "
            f"Unexpected ({len(unexpected_genomes)}): {sorted(unexpected_genomes)[:5]}"
        )

    # Check for duplicate evaluations
    dup_mask = full_summary_df.duplicated(subset=["genome", "scenario_id", "biomass_recipe_id"])
    if dup_mask.any():
        n_dups = int(dup_mask.sum())
        raise RuntimeError(f"Duplicate evaluations detected in collector: {n_dups} duplicate rows found")

    # 2. Strict NA assignment for invalid/infeasible models (`EXP1-ACCEPTANCE-AND-COLLECTOR-REVIEW.txt` Item 3)
    invalid_mask = (
        (~full_summary_df["solver_status"].isin(REVALIDATED_FEASIBLE_STATUSES))
        | (full_summary_df["n_distinct_valid_solutions"] == 0)
    )
    metric_cols_to_na = [
        "min_added_rxns_count",
        "max_added_rxns_count",
        "n_rxns_observed_in_all_enumerated_models",
        "n_rxns_essential_in_declared_supernetwork",
        "min_total_cost",
        "max_ensemble_cost",
    ]
    for c in metric_cols_to_na:
        if c in full_summary_df.columns:
            full_summary_df.loc[invalid_mask, c] = np.nan

    # Write receipts
    with open(args.out_dir / "exp1_stage3_shard_collection_receipts.json", "w", encoding="utf-8") as fh:
        json.dump(
            {
                "n_shards": args.num_tasks,
                "n_total_expected_eligible_genomes": len(expected_eligible_genomes),
                "n_total_collected_genomes": len(collected_genomes),
                "exact_genome_ids_matched": True,
                "no_duplicate_evaluations": True,
                "shards": shard_receipts,
            },
            fh,
            indent=2,
        )

    full_summary_df.to_csv(args.out_dir / "exp1_full484_ensemble_summary.tsv.gz", sep="\t", index=False, compression="gzip")
    full_rf_df.to_csv(
        args.out_dir / "exp1_full484_reaction_stability_and_biomass_sensitivity.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    full_wkf_df.to_csv(
        args.out_dir / "exp1_full484_withdrawal_knockout_fva.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    if not full_supp_df.empty:
        full_supp_df.to_csv(args.out_dir / "exp1_full484_supplement_rescue.tsv", sep="\t", index=False)

    # 3. Cohort & Quality Comparisons across Scenarios and Biomass Perturbations
    comp_rows: list[dict[str, Any]] = []
    for (scen_id, bio_id), grp in full_summary_df.groupby(["scenario_id", "biomass_recipe_id"], dropna=False):
        for col in (
            "min_total_cost",
            "min_added_rxns_count",
            "n_rxns_observed_in_all_enumerated_models",
            "n_rxns_essential_in_declared_supernetwork",
        ):
            st_row = compute_mannwhitney_and_spearman(grp, col)
            st_row["scenario_id"] = scen_id
            st_row["biomass_recipe_id"] = bio_id
            comp_rows.append(st_row)

    pd.DataFrame(comp_rows).to_csv(
        args.out_dir / "exp1_full484_cohort_and_quality_comparisons.tsv", sep="\t", index=False
    )

    # 4. Actual canonical order/family/genus descriptive comparisons (not pooling missing taxonomy placeholders)
    #    and legacy family sensitivity comparison (`EXP1-ACCEPTANCE-AND-COLLECTOR-REVIEW.txt` Item 3).
    leg_map = dict(zip(inv_df["genome"].astype(str), inv_df["legacy_classification"].astype(str)))
    s4_b2 = full_summary_df[
        (full_summary_df["scenario_id"] == SCEN_S4_OXIC_DESERT)
        & (full_summary_df["biomass_recipe_id"] == BIO_B2_OMIT_7_KEEP_F420)
    ].copy()
    s4_b2["legacy_classification"] = s4_b2["genome"].map(leg_map)
    s4_b2["legacy_family"] = s4_b2["legacy_classification"].apply(
        lambda s: next((p for p in str(s).split(";") if p.startswith("f__")), "f__")
    )

    lineage_comparisons: list[tuple[str, str, str]] = [
        ("canonical_gtdb_r226_r232", "order", "canonical_order"),
        ("canonical_gtdb_r226_r232", "family", "canonical_family"),
        ("canonical_gtdb_r226_r232", "genus", "canonical_genus"),
        ("legacy_gtdb_sensitivity", "family", "legacy_family"),
    ]

    lineage_rows: list[dict[str, Any]] = []
    for tax_scheme, rank, col_name in lineage_comparisons:
        if col_name not in s4_b2.columns:
            continue
        for taxon_val, tgrp in s4_b2.groupby(col_name, dropna=False):
            tax_str = str(taxon_val).strip()
            # Do NOT pool missing taxonomy placeholders into shared lineages
            is_placeholder = (
                tax_str.lower() in MISSING_TAXONOMY_PLACEHOLDERS
                or tax_str.endswith("__")
                or tax_str == ""
            )

            all_eq = tgrp[tgrp["cohort"] == "EQ_desert"]
            all_ext = tgrp[tgrp["cohort"] == "external_comparator"]
            n_eq_total = len(all_eq)
            n_ext_total = len(all_ext)
            if n_eq_total == 0 and n_ext_total == 0:
                continue

            eq_valid_added = all_eq["min_added_rxns_count"].dropna().to_numpy(dtype=float)
            ext_valid_added = all_ext["min_added_rxns_count"].dropna().to_numpy(dtype=float)
            eq_valid_ess = all_eq["n_rxns_essential_in_declared_supernetwork"].dropna().to_numpy(dtype=float)
            ext_valid_ess = all_ext["n_rxns_essential_in_declared_supernetwork"].dropna().to_numpy(dtype=float)

            is_shared = (n_eq_total > 0 and n_ext_total > 0 and not is_placeholder)

            if is_shared and len(eq_valid_added) > 0 and len(ext_valid_added) > 0:
                _, p_added = stats.mannwhitneyu(eq_valid_added, ext_valid_added, alternative="two-sided")
            else:
                p_added = math.nan

            lineage_rows.append(
                {
                    "taxonomy_scheme": tax_scheme,
                    "taxonomic_rank": rank,
                    "taxon_name": tax_str,
                    "is_missing_taxonomy_placeholder": is_placeholder,
                    "is_shared_lineage": is_shared,
                    "n_eq_genomes_total": n_eq_total,
                    "n_eq_genomes_valid_models": len(eq_valid_added),
                    "n_external_genomes_total": n_ext_total,
                    "n_external_genomes_valid_models": len(ext_valid_added),
                    "eq_median_min_added_rxns_B2_S4": float(np.median(eq_valid_added)) if len(eq_valid_added) > 0 else math.nan,
                    "external_median_min_added_rxns_B2_S4": float(np.median(ext_valid_added)) if len(ext_valid_added) > 0 else math.nan,
                    "exploratory_mannwhitney_p_added_rxns": float(p_added) if not math.isnan(p_added) else math.nan,
                    "eq_median_essential_supernetwork_rxns_B2_S4": float(np.median(eq_valid_ess)) if len(eq_valid_ess) > 0 else math.nan,
                    "external_median_essential_supernetwork_rxns_B2_S4": float(np.median(ext_valid_ess)) if len(ext_valid_ess) > 0 else math.nan,
                }
            )

    pd.DataFrame(lineage_rows).to_csv(
        args.out_dir / "exp1_full484_lineage_matched_comparisons_order_family_genus.tsv",
        sep="\t",
        index=False,
    )

    # 5. Emit Supported Candidate Deficiencies for Stage 4 Sequence/Synteny Verification (`EXP1-ACCEPTANCE-AND-COLLECTOR-REVIEW.txt` Item 4)
    #    Criteria:
    #    - High-quality / near-complete EQ MAGs (completeness >= 85.0, contamination <= 5.0),
    #    - Directly verified scenario/biomass results for asserted robustness:
    #      'robust_across_all_tested_biomass_perturbations' marked directly verified;
    #      incomplete sensitivity marked exploratory and explicitly segregated.
    #    - Preserves scenario_id condition in candidate summary table!
    #    - Termed 'sequence_synteny_followup_selection_not_functional_proof' (NOT universal necessity).
    eq_hq_rf = full_rf_df[
        (full_rf_df["cohort"] == "EQ_desert")
        & (full_rf_df["completeness"] >= 85.0)
        & (full_rf_df["contamination"] <= 5.0)
        & (full_rf_df["freq_B2_omit_7_keep_f420"] >= 0.80)
        & (full_rf_df["supernetwork_ko_necessity_B2"] == "model_specific_essential_in_declared_supernetwork")
    ].copy()

    if not eq_hq_rf.empty:
        is_directly_verified_robust = (
            eq_hq_rf["biomass_sensitivity_classification"] == "robust_across_all_tested_biomass_perturbations"
        )
        eq_hq_rf["followup_selection_tier"] = np.where(
            is_directly_verified_robust,
            "directly_verified_robust_candidate",
            "exploratory_incomplete_biomass_sensitivity_candidate",
        )
        eq_hq_rf["epistemic_status"] = "sequence_synteny_followup_selection_not_functional_proof"

        cand_def = eq_hq_rf[
            ~eq_hq_rf["biomass_sensitivity_classification"].str.contains("vanishes_when_7_non_f420", na=False)
        ].copy()

        cand_def_summary = (
            cand_def.groupby(
                [
                    "scenario_id",
                    "rxn_id",
                    "base_id",
                    "rxn_name",
                    "is_transport",
                    "biomass_sensitivity_classification",
                    "followup_selection_tier",
                    "epistemic_status",
                ],
                dropna=False,
            )
            .agg(
                n_hq_eq_genomes_requiring=("genome", "nunique"),
                n_hq_eq_families_requiring=("canonical_family", "nunique"),
                median_freq_B2=("freq_B2_omit_7_keep_f420", "median"),
                median_candidate_weight=("candidate_weight", "median"),
                example_genomes_csv=("genome", lambda s: ",".join(sorted(set(s.astype(str)))[:10])),
                families_csv=("canonical_family", lambda s: ",".join(sorted(set(s.astype(str))))),
            )
            .reset_index()
            .sort_values(
                ["followup_selection_tier", "n_hq_eq_genomes_requiring", "n_hq_eq_families_requiring"],
                ascending=[True, False, False],
            )
        )
    else:
        cand_def = pd.DataFrame()
        cand_def_summary = pd.DataFrame(
            columns=[
                "scenario_id",
                "rxn_id",
                "base_id",
                "rxn_name",
                "is_transport",
                "biomass_sensitivity_classification",
                "followup_selection_tier",
                "epistemic_status",
                "n_hq_eq_genomes_requiring",
                "n_hq_eq_families_requiring",
                "median_freq_B2",
                "median_candidate_weight",
                "example_genomes_csv",
                "families_csv",
            ]
        )

    cand_def.to_csv(args.out_dir / "exp1_stage4_supported_candidate_deficiencies_instances.tsv", sep="\t", index=False)
    cand_def_summary.to_csv(args.out_dir / "exp1_stage4_supported_candidate_deficiencies_summary.tsv", sep="\t", index=False)

    print(
        f"[Stage 3 Collector Complete] Aggregated {len(collected_genomes)} genomes across {args.num_tasks} shards, "
        f"{len(cand_def)} supported HQ EQ deficiency instances across {len(cand_def_summary)} reactions."
    )


if __name__ == "__main__":
    main()
