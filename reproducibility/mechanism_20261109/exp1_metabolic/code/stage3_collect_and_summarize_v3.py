#!/usr/bin/env python3
"""Stage 3 Collector (v3) of `exp1_metabolic`:
1. Aggregates the 8 shards across all 484 deduplicated genomes (237 EQ MAGs + 247 external comparators).
2. Validates exact ID coverage, shard completeness, and fail-closed per-solution ledger correspondence.
3. Strict NA assignment for invalid/infeasible models in summary metrics.
4. Corrected unknown-as-zero biomass sensitivity & reaction-frequency repair (`EXP1-FINAL-CONTROLS-AND-COLLECTOR.txt` Item 2):
   - Joins exact (genome, scenario, biomass) summary status and valid-solution counts to reaction rows.
   - If an ensemble failed or is infeasible, frequency is NaN and status is explicit (never 0.0 or "absent").
   - Invalid/no-baseline sensitivity is strictly NaN.
   - Biomass sensitivity classification recomputed: observed across tested recipes is labeled
     `observed_across_all_tested_biomass_perturbations` (NOT `directly_verified_robust_candidate`).
   - B2 supernetwork knockout is labeled as a conditional B2 requirement (`conditional_b2_essential_candidate`).
5. Integrates same-medium 7 HQ external masking controls (`exp1_masking_controls_7hq.tsv`).
6. Preserves schema-valid zero eligible outputs for Stage 4 sequence/synteny verification.
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
    BIO_B0_BASELINE,
    BIO_B1_OMIT_8,
    BIO_B2_OMIT_7_KEEP_F420,
    BIO_B3_OMIT_LIPIDS,
    BIO_B4_CORE_SKELETON,
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
    """
    sub = df[df[metric_col].notna()].copy()
    eq_sub = sub[sub["cohort"] == "EQ_desert"]
    ext_sub = sub[sub["cohort"] == "external_comparator"]
    eq_vals = eq_sub[metric_col].to_numpy(dtype=float)
    ext_vals = ext_sub[metric_col].to_numpy(dtype=float)

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

    eq_comp = sub[(sub["cohort"] == "EQ_desert") & (sub["completeness"].notna())]
    if len(eq_comp) >= 5:
        rho_comp, p_comp = stats.spearmanr(eq_comp["completeness"], eq_comp[metric_col])
    else:
        rho_comp, p_comp = math.nan, math.nan

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


def repair_reaction_frequency_and_biomass_sensitivity(
    rf_df: pd.DataFrame,
    summary_df: pd.DataFrame,
) -> pd.DataFrame:
    """Repair reaction-frequency table by joining exact (genome, scenario, biomass) summary status
    and valid-solution counts (`EXP1-FINAL-CONTROLS-AND-COLLECTOR.txt` Item 2).
    """
    if rf_df.empty:
        return rf_df.copy()

    repaired_rf = rf_df.copy()

    # Build exact status and valid-solution lookup: (genome, scenario_id, biomass_recipe_id) -> (status, n_valid)
    status_map: dict[tuple[str, str, str], tuple[str, int]] = {}
    for _, srow in summary_df.iterrows():
        k = (str(srow["genome"]), str(srow["scenario_id"]), str(srow["biomass_recipe_id"]))
        st = str(srow["solver_status"])
        nsol = int(srow["n_distinct_valid_solutions"]) if pd.notna(srow.get("n_distinct_valid_solutions")) else 0
        status_map[k] = (st, nsol)

    recipes = [
        ("freq_B0_template_baseline", "stability_tier_B0", "supernetwork_ko_necessity_B0", BIO_B0_BASELINE),
        ("freq_B1_omit_all8_cofactors", None, None, BIO_B1_OMIT_8),
        ("freq_B2_omit_7_keep_f420", "stability_tier_B2", "supernetwork_ko_necessity_B2", BIO_B2_OMIT_7_KEEP_F420),
        ("freq_B3_omit_surface_lipids", None, None, BIO_B3_OMIT_LIPIDS),
        ("freq_B4_core_skeleton", None, None, BIO_B4_CORE_SKELETON),
    ]

    new_rows: list[dict[str, Any]] = []
    for _, r in repaired_rf.iterrows():
        row_dict = dict(r)
        g = str(r["genome"])
        scen = str(r["scenario_id"])

        # Determine validity of each biomass perturbation recipe
        validity: dict[str, bool] = {}
        statuses: dict[str, str] = {}
        for freq_col, stab_col, ko_col, bio_id in recipes:
            st, nsol = status_map.get((g, scen, bio_id), ("not_evaluated", 0))
            is_valid = bool(st in REVALIDATED_FEASIBLE_STATUSES and nsol > 0)
            validity[bio_id] = is_valid
            statuses[bio_id] = st

            if not is_valid:
                # Strictly NA frequency when ensemble failed or is infeasible
                row_dict[freq_col] = np.nan
                if stab_col and stab_col in row_dict:
                    row_dict[stab_col] = f"unsolved_or_infeasible_{st}"
                if ko_col and ko_col in row_dict:
                    row_dict[ko_col] = f"unsolved_or_infeasible_{st}"
            else:
                # Honest reporting: zero frequency in finite heuristic ensemble means not observed, not proven absent
                f_val = row_dict.get(freq_col)
                if pd.notna(f_val) and float(f_val) == 0.0:
                    if stab_col and stab_col in row_dict and row_dict[stab_col] == "absent":
                        row_dict[stab_col] = "not_observed_in_finite_ensemble"
                    if ko_col and ko_col in row_dict and str(row_dict[ko_col]).startswith("not_in_"):
                        row_dict[ko_col] = "not_observed_in_finite_ensemble"

        # Difference B0 - B2
        if validity[BIO_B0_BASELINE] and validity[BIO_B2_OMIT_7_KEEP_F420]:
            f0 = float(row_dict["freq_B0_template_baseline"])
            f2 = float(row_dict["freq_B2_omit_7_keep_f420"])
            row_dict["freq_difference_B0_minus_B2"] = round(f0 - f2, 6)
        else:
            row_dict["freq_difference_B0_minus_B2"] = np.nan

        # Biomass sensitivity classification
        if not validity[BIO_B0_BASELINE]:
            bio_class = f"unknown_biomass_sensitivity__baseline_B0_{statuses[BIO_B0_BASELINE]}"
        elif not validity[BIO_B2_OMIT_7_KEEP_F420]:
            bio_class = f"unknown_biomass_sensitivity__B2_{statuses[BIO_B2_OMIT_7_KEEP_F420]}"
        else:
            # Check other perturbations if scheduled (S3..S6)
            other_valid = all(validity[b] for b in (BIO_B1_OMIT_8, BIO_B3_OMIT_LIPIDS, BIO_B4_CORE_SKELETON))
            if not other_valid:
                # Some perturbation failed: cannot confirm robust across all assumptions
                in_b0 = (float(row_dict.get("freq_B0_template_baseline", 0.0) or 0.0) > 0.0)
                in_b2 = (float(row_dict.get("freq_B2_omit_7_keep_f420", 0.0) or 0.0) > 0.0)
                if in_b0 and not in_b2:
                    bio_class = "biomass_assumption_sensitive__vanishes_when_7_non_f420_template_cofactors_omitted"
                else:
                    bio_class = "unknown_biomass_sensitivity__incomplete_perturbations"
            else:
                f0 = float(row_dict["freq_B0_template_baseline"])
                f1 = float(row_dict["freq_B1_omit_all8_cofactors"])
                f2 = float(row_dict["freq_B2_omit_7_keep_f420"])
                f3 = float(row_dict["freq_B3_omit_surface_lipids"])
                f4 = float(row_dict["freq_B4_core_skeleton"])

                in_b0 = (f0 > 0.0)
                in_b1 = (f1 > 0.0)
                in_b2 = (f2 > 0.0)
                in_b3 = (f3 > 0.0)
                in_b4 = (f4 > 0.0)

                if in_b0 and in_b1 and in_b2 and in_b3 and in_b4:
                    bio_class = "observed_across_all_tested_biomass_perturbations"
                elif in_b0 and not (in_b1 and in_b2 and in_b3 and in_b4):
                    if not in_b2:
                        bio_class = "biomass_assumption_sensitive__vanishes_when_7_non_f420_template_cofactors_omitted"
                    elif in_b2 and not in_b1:
                        bio_class = "biomass_assumption_sensitive__vanishes_when_f420_omitted"
                    elif in_b1 and not in_b3:
                        bio_class = "biomass_assumption_sensitive__vanishes_when_specialized_lipids_or_homospermidine_omitted"
                    else:
                        bio_class = "biomass_assumption_sensitive__vanishes_in_core_skeleton_perturbation"
                else:
                    bio_class = "not_in_baseline_b0_ensemble"

        row_dict["biomass_sensitivity_classification"] = bio_class
        new_rows.append(row_dict)

    return pd.DataFrame(new_rows)


def main() -> None:
    ap = argparse.ArgumentParser(description="Stage 3 Collector (v3): Aggregate 484-genome shards with repaired sensitivity.")
    ap.add_argument("--inventory-484-tsv", type=Path, required=True)
    ap.add_argument("--shards-dir", type=Path, required=True)
    ap.add_argument("--masking-controls-tsv", type=Path, default=None)
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
        "is_transport", "candidate_weight", "prior_evidence", "freq_B0_template_baseline",
        "freq_B1_omit_all8_cofactors", "freq_B2_omit_7_keep_f420", "freq_B3_omit_surface_lipids",
        "freq_B4_core_skeleton", "freq_difference_B0_minus_B2",
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

    # Concatenate shards
    full_summary_df = (
        pd.concat(sum_dfs, ignore_index=True) if sum_dfs else pd.DataFrame(columns=summary_columns)
    ).sort_values(["cohort", "genome", "scenario_id", "biomass_recipe_id"]).reset_index(drop=True)

    raw_rf_df = (
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

    # 1. Exact ID Coverage Verification
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

    dup_mask = full_summary_df.duplicated(subset=["genome", "scenario_id", "biomass_recipe_id"])
    if dup_mask.any():
        raise RuntimeError(f"Duplicate evaluations detected in collector: {int(dup_mask.sum())} duplicate rows")

    # 2. Strict NA assignment for invalid/infeasible models in summary metrics
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

    # 3. Fail-Closed Per-Solution Ledger Correspondence Audit
    if not full_ledger_df.empty:
        full_ledger_df = full_ledger_df.sort_values(
            ["cohort", "genome", "scenario_id", "biomass_recipe_id", "solution_index"]
        ).reset_index(drop=True)
        full_ledger_df.to_csv(args.out_dir / "exp1_full484_per_solution_ledger.tsv", sep="\t", index=False)

        acc_ledger = full_ledger_df[full_ledger_df["validity_disposition"] == "accepted_feasible"]
        acc_counts = acc_ledger.groupby(["genome", "scenario_id", "biomass_recipe_id"]).size().to_dict() if not acc_ledger.empty else {}
        for _, srow in full_summary_df.iterrows():
            k = (str(srow["genome"]), str(srow["scenario_id"]), str(srow["biomass_recipe_id"]))
            exp_count = int(srow["n_distinct_valid_solutions"]) if srow["solver_status"] in REVALIDATED_FEASIBLE_STATUSES and pd.notna(srow["n_distinct_valid_solutions"]) else 0
            act_count = acc_counts.get(k, 0)
            if exp_count != act_count:
                raise RuntimeError(f"Full484 ledger count mismatch for {k}: summary claims {exp_count}, ledger has {act_count}")

    # 4. Correct Reaction-Frequency and Biomass Sensitivity Table (`EXP1-FINAL-CONTROLS-AND-COLLECTOR.txt` Item 2)
    full_rf_df = repair_reaction_frequency_and_biomass_sensitivity(raw_rf_df, full_summary_df)

    # 5. Masking Controls Verification
    masking_summary: dict[str, Any] = {}
    if args.masking_controls_tsv and args.masking_controls_tsv.exists():
        mask_df = pd.read_csv(args.masking_controls_tsv, sep="\t")
        mask_dest = args.out_dir / args.masking_controls_tsv.name
        mask_df.to_csv(mask_dest, sep="\t", index=False)
        rec_col = (
            "masked_rxn_recovery_rate_under_medium"
            if "masked_rxn_recovery_rate_under_medium" in mask_df.columns
            else "true_masked_rxn_recovery_rate"
        )
        eval_mask = mask_df["is_evaluable"] == True if "is_evaluable" in mask_df.columns else pd.Series(False, index=mask_df.index)
        mean_rec = float(mask_df[eval_mask][rec_col].dropna().mean()) if (
            eval_mask.any() and rec_col in mask_df.columns
        ) else math.nan
        masking_summary = {
            "masking_controls_present": True,
            "masking_controls_file": args.masking_controls_tsv.name,
            "n_masking_rows": len(mask_df),
            "n_evaluable_masking_rows": int(eval_mask.sum()),
            "mean_recovery_rate": mean_rec,
        }
        # Also copy companion ledger if present
        led_src = args.masking_controls_tsv.parent / (args.masking_controls_tsv.stem + "_per_solution_ledger.tsv")
        if led_src.exists():
            led_df = pd.read_csv(led_src, sep="\t")
            led_df.to_csv(args.out_dir / led_src.name, sep="\t", index=False)
            masking_summary["ledger_file"] = led_src.name
            masking_summary["n_ledger_rows"] = len(led_df)
    else:
        masking_summary = {"masking_controls_present": False}

    # Write primary full484 outputs
    full_summary_df.to_csv(args.out_dir / "exp1_full484_ensemble_summary.tsv.gz", sep="\t", index=False, compression="gzip")
    full_rf_df.to_csv(args.out_dir / "exp1_full484_reaction_stability_and_biomass_sensitivity.tsv.gz", sep="\t", index=False, compression="gzip")
    full_wkf_df.to_csv(args.out_dir / "exp1_full484_withdrawal_knockout_fva.tsv.gz", sep="\t", index=False, compression="gzip")
    if not full_supp_df.empty:
        full_supp_df.to_csv(args.out_dir / "exp1_full484_supplement_rescue.tsv", sep="\t", index=False)

    # 6. Descriptive Comparisons
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

    pd.DataFrame(comp_rows).to_csv(args.out_dir / "exp1_full484_cohort_and_quality_comparisons.tsv", sep="\t", index=False)

    # 7. Lineage matched comparisons
    leg_map = dict(zip(inv_df["genome"].astype(str), inv_df["legacy_classification"].astype(str)))
    s4_b2 = full_summary_df[
        (full_summary_df["scenario_id"] == SCEN_S4_OXIC_DESERT)
        & (full_summary_df["biomass_recipe_id"] == BIO_B2_OMIT_7_KEEP_F420)
    ].copy()
    s4_b2["legacy_classification"] = s4_b2["genome"].map(leg_map)
    s4_b2["legacy_family"] = s4_b2["legacy_classification"].apply(
        lambda s: next((p for p in str(s).split(";") if p.startswith("f__")), "f__")
    )

    lineage_comparisons = [
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

    # 8. Emit Supported Candidate Deficiencies for Stage 4 Sequence/Synteny Verification
    #    CRITICAL HONEST TIERS (`EXP1-FINAL-CONTROLS-AND-COLLECTOR.txt` Item 2):
    #    - High-quality EQ MAGs (completeness >= 85.0, contamination <= 5.0)
    #    - Positive frequency in B2 (freq >= 0.80)
    #    - Model-specific essential in B2 declared supernetwork
    #    - Valid evaluated B2 ensemble (not NaN)
    #    - Does NOT vanish when 7 non-F420 template cofactors omitted
    #    - Does NOT have unknown sensitivity due to unsolved/infeasible models
    #    - Labels tier honestly:
    #      * `conditional_b0_b2_essential_candidate` if KO necessity verified in both B0 and B2 AND observed across perturbations
    #      * `conditional_b2_essential_candidate` if KO necessity established only under B2
    #      * NEVER `directly_verified_robust_candidate` without multi-recipe necessity proof
    #    - Epistemic status: `conditional_b2_supernetwork_requirement_sequence_synteny_followup_not_functional_proof`
    b2_col = "freq_B2_omit_7_keep_f420"
    b2_ko_col = "supernetwork_ko_necessity_B2"
    b0_ko_col = "supernetwork_ko_necessity_B0"

    eq_hq_rf = full_rf_df[
        (full_rf_df["cohort"] == "EQ_desert")
        & (full_rf_df["completeness"] >= 85.0)
        & (full_rf_df["contamination"] <= 5.0)
        & (full_rf_df[b2_col].notna())
        & (full_rf_df[b2_col] >= 0.80)
        & (full_rf_df[b2_ko_col] == "model_specific_essential_in_declared_supernetwork")
    ].copy()

    if not eq_hq_rf.empty:
        # Exclude artifacts and unknowns
        non_artifact = (
            ~eq_hq_rf["biomass_sensitivity_classification"].str.contains("vanishes_when_7_non_f420", na=False)
            & ~eq_hq_rf["biomass_sensitivity_classification"].str.startswith("unknown_biomass_sensitivity", na=False)
        )
        cand_def = eq_hq_rf[non_artifact].copy()

        if not cand_def.empty:
            has_b0_ko = (cand_def[b0_ko_col] == "model_specific_essential_in_declared_supernetwork")
            is_obs_all = (cand_def["biomass_sensitivity_classification"] == "observed_across_all_tested_biomass_perturbations")
            cand_def["followup_selection_tier"] = np.where(
                has_b0_ko & is_obs_all,
                "conditional_b0_b2_essential_candidate",
                "conditional_b2_essential_candidate",
            )
            cand_def["epistemic_status"] = "conditional_b2_supernetwork_requirement_sequence_synteny_followup_not_functional_proof"

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
                    median_freq_B2=(b2_col, "median"),
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
                    "scenario_id", "rxn_id", "base_id", "rxn_name", "is_transport",
                    "biomass_sensitivity_classification", "followup_selection_tier", "epistemic_status",
                    "n_hq_eq_genomes_requiring", "n_hq_eq_families_requiring", "median_freq_B2",
                    "median_candidate_weight", "example_genomes_csv", "families_csv",
                ]
            )
    else:
        cand_def = pd.DataFrame()
        cand_def_summary = pd.DataFrame(
            columns=[
                "scenario_id", "rxn_id", "base_id", "rxn_name", "is_transport",
                "biomass_sensitivity_classification", "followup_selection_tier", "epistemic_status",
                "n_hq_eq_genomes_requiring", "n_hq_eq_families_requiring", "median_freq_B2",
                "median_candidate_weight", "example_genomes_csv", "families_csv",
            ]
        )

    cand_def.to_csv(args.out_dir / "exp1_stage4_supported_candidate_deficiencies_instances.tsv", sep="\t", index=False)
    cand_def_summary.to_csv(args.out_dir / "exp1_stage4_supported_candidate_deficiencies_summary.tsv", sep="\t", index=False)

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
                "masking_controls": masking_summary,
                "n_supported_candidate_instances": len(cand_def),
                "n_supported_candidate_reactions": len(cand_def_summary),
            },
            fh,
            indent=2,
        )

    print(
        f"[Stage 3 Collector v3 Complete] Aggregated {len(collected_genomes)} genomes across {args.num_tasks} shards, "
        f"{len(cand_def)} supported HQ EQ deficiency instances across {len(cand_def_summary)} reactions."
    )


if __name__ == "__main__":
    main()
