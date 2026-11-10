#!/usr/bin/env python3
"""Stage 2 Focused Release Validation Script.
Executes the bounded real focused release validation set on Shaheen per `EXP1-FULL-RELEASE-DECISION.txt`:
1. The 68 problematic scenario cases that failed in run_15873653 (64 numerical coupling violations + 4 fba verification failures).
2. Representative previously-feasible cases (optimal and heuristic_incumbent_bounded).
3. Representative proven-infeasible cases (infeasible and infeasible_under_fixed_medium).

Uses original native constraints and corrected residual/bound checks in `solve_fba_highs`.
Records both aggregate summary rows and compact per-solution ledger rows (`exp1_focused_validation_per_solution_ledger.tsv`).
Enforces strict per-solution correctness verification: zero invalid accepted solutions, explicit unknown/infeasible outcomes.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from exp1_sbml_highs_engine import (
    BIOMASS_PERTURBATION_RECIPES,
    NUTRIENT_SCENARIOS,
    OBSERVED_COHORT_538_POOL_SCOPE,
    SCEN_S0_NATIVE,
    SCEN_S3_ANOXIC_DESERT,
    SCEN_S4_OXIC_DESERT,
    SCEN_S5_ANOXIC_DEFINED,
    SCEN_S6_OXIC_DEFINED,
    enumerate_alternative_gapfill_ensemble,
)
from stage2_pilot24_ensemble import (
    REQUIRED_VERIFICATION_COLUMNS,
    build_genome_candidate_weights,
    load_genome_archive_payload,
    validate_pilot24_scientific_gate,
)

REVALIDATED_FEASIBLE_STATUSES = {
    "optimal",
    "optimal_after_cycle_qc_filtering",
    "optimal_draft_no_gapfill_needed",
    "heuristic_incumbent_bounded",
    "post_pruning_cost_changed_not_proven_optimal",
}
PROVEN_INFEASIBLE_STATUSES = {
    "infeasible",
    "infeasible_under_fixed_medium",
}
TIMEOUT_STATUSES = {
    "time_limit_no_incumbent",
    "time_limit_heuristic_incumbent",
}
NUMERICAL_FAILURE_STATUSES = {
    "numerical_coupling_violation",
    "fba_verification_failed",
    "numerical_error",
}
KNOWN_STATUSES = REVALIDATED_FEASIBLE_STATUSES | PROVEN_INFEASIBLE_STATUSES | TIMEOUT_STATUSES | NUMERICAL_FAILURE_STATUSES


def run_single_validation_task(task: dict[str, Any]) -> dict[str, Any]:
    genome = task["genome"]
    cohort = task.get("cohort", "unknown")
    canon_order = task.get("canonical_order", "")
    canon_family = task.get("canonical_family", "")
    scen_id = task["scenario_id"]
    bio_id = task["biomass_recipe_id"]
    old_status = task.get("old_status", "unknown")
    case_type = task.get("case_type", "unknown")
    gapseq_out_dir = Path(task["gapseq_out_dir"])
    supp_models_dir = Path(task["supp_models_dir"])
    universal_rxns = task["universal_rxns"]
    milp_time_limit_s = float(task.get("milp_time_limit_s", 10.0))

    payload = load_genome_archive_payload(genome, gapseq_out_dir, supp_models_dir)
    draft_mod = payload["draft_model"]
    cand_weights = build_genome_candidate_weights(universal_rxns, payload["rxn_evi"], payload["trans_evi"])

    scenarios = dict(NUTRIENT_SCENARIOS)
    if payload.get("desert_med"):
        d_med = payload["desert_med"]
        scenarios[SCEN_S3_ANOXIC_DESERT] = {
            "description": "Per-genome gapseq desert medium (anoxic)",
            "redox_regime": "anoxic",
            "compounds": {cid: (cname, u) for cid, (cname, u) in d_med.items() if cid != "cpd00007"},
        }
        scenarios[SCEN_S4_OXIC_DESERT] = {
            "description": "Per-genome gapseq desert medium (oxic with O2)",
            "redox_regime": "oxic",
            "compounds": {**d_med, "cpd00007": ("O2", 10.0)},
        }
    if payload.get("native_med"):
        scenarios[SCEN_S0_NATIVE] = {
            "description": "Per-genome gapseq predicted native medium",
            "redox_regime": "oxic" if "cpd00007" in payload["native_med"] else "anoxic",
            "compounds": payload["native_med"],
        }

    scen_info = scenarios[scen_id]
    scen_cpds = scen_info["compounds"]

    n_sol_cap = 20 if bio_id in {"B0_methanosarcina_template_baseline", "B2_perturbation_omit_7_non_f420_template_cofactors_keep_f420"} else 10

    t0 = time.time()
    ens = enumerate_alternative_gapfill_ensemble(
        draft_model=draft_mod,
        candidate_reactions=universal_rxns,
        candidate_weights=cand_weights,
        scenario_compounds=scen_cpds,
        biomass_recipe_id=bio_id,
        min_growth=0.01,
        n_Target_solutions=n_sol_cap,
        cost_slack_frac=0.25,
        candidate_pool_scope=OBSERVED_COHORT_538_POOL_SCOPE,
        time_limit_per_milp_s=milp_time_limit_s,
    )
    dt = time.time() - t0

    sols = ens.get("solutions", [])
    opt_sols = ens.get("optimization_feasible_solutions", [])
    solver_st = ens["solver_status"]
    is_reval_feasible = solver_st in REVALIDATED_FEASIBLE_STATUSES

    if sols:
        n_added_min = float(min(len(s["added_rxns"]) for s in sols))
        n_added_max = float(max(len(s["added_rxns"]) for s in sols))
        n_obs_all = float(sum(
            1 for rf in ens["reaction_frequencies"].values()
            if rf["stability_tier"] == "observed_in_all_enumerated_models"
        ))
        n_super_essential = float(sum(
            1 for rf in ens["reaction_frequencies"].values()
            if rf["supernetwork_ko_necessity"] == "model_specific_essential_in_declared_supernetwork"
        ))
    else:
        n_added_min = math.nan
        n_added_max = math.nan
        n_obs_all = math.nan
        n_super_essential = math.nan

    summary_row = {
        "genome": genome,
        "cohort": cohort,
        "scenario_id": scen_id,
        "biomass_recipe_id": bio_id,
        "case_type": case_type,
        "old_status": old_status,
        "solver_status": solver_st,
        "optimizer_status": ens["optimizer_status"],
        "raw_solver_status": ens["raw_solver_status"],
        "draft_already_grows": ens["draft_already_grows"],
        "draft_growth_rate": ens["draft_growth_rate"],
        "min_total_cost": ens["min_total_cost"],
        "max_ensemble_cost": ens["max_ensemble_cost"],
        "first_solve_mip_gap": ens["first_solve_mip_gap"],
        "first_solve_dual_bound": ens["first_solve_dual_bound"],
        "n_solutions_attempted": ens["n_solutions_attempted"],
        "n_optimization_feasible_solutions": ens["n_optimization_feasible_solutions"],
        "n_cycle_qc_rejected_solutions": ens["n_cycle_qc_rejected_solutions"],
        "n_distinct_valid_solutions": ens["n_distinct_valid_solutions"],
        "fewer_than_target_reason": ens["fewer_than_target_reason"],
        "min_added_rxns_count": n_added_min,
        "max_added_rxns_count": n_added_max,
        "n_rxns_observed_in_all_enumerated_models": n_obs_all,
        "n_rxns_essential_in_declared_supernetwork": n_super_essential,
        "max_solution_residual": ens.get("max_solution_residual", math.nan),
        "max_solution_bound_viol": ens.get("max_solution_bound_viol", math.nan),
        "min_solution_verified_growth": ens.get("min_solution_verified_growth", math.nan),
        "runtime_seconds": dt,
    }

    ledger_rows: list[dict[str, Any]] = []
    if sols:
        for s in sols:
            r_ids = sorted(s.get("added_rxns", []))
            res_val = float(s.get("fba_check_max_residual", math.nan))
            bnd_val = float(s.get("fba_check_max_bound_viol", math.nan))
            gr_val = float(s.get("verified_fba_growth", math.nan))
            cyc_st = str(s.get("cycle_qc_status", "unknown"))
            spur_atp = s.get("post_pruning_egc_has_spurious_atp_cycle", None)
            science_eligible = s.get("science_qc_eligible", None)

            # Explicit cycle-QC pass: must be explicitly evaluated and pass with no spurious cycle,
            # or legitimate documented case where ATP species is absent (not required)
            valid_cycle_pass = (cyc_st in {"passed_closed_uptake_egc_check", "passed_no_free_atp_flux"} and spur_atp is False)
            legitimate_no_atp = (cyc_st in {"evaluated_no_atp_species", "not_required_no_atp_species"} and spur_atp is not True)
            cycle_qc_ok = bool((valid_cycle_pass or legitimate_no_atp) and (science_eligible is not False))

            checks_ok = bool(
                math.isfinite(res_val) and 0.0 <= res_val <= 1e-5
                and math.isfinite(bnd_val) and 0.0 <= bnd_val <= 1e-5
                and math.isfinite(gr_val) and gr_val >= 0.00999
                and cycle_qc_ok
            )
            disposition = (
                "accepted_feasible" if (is_reval_feasible and checks_ok)
                else ("rejected_cycle_qc" if (is_reval_feasible and not cycle_qc_ok)
                else "unverified_solution")
            )
            ledger_rows.append({
                "genome": genome,
                "cohort": cohort,
                "scenario_id": scen_id,
                "biomass_recipe_id": bio_id,
                "case_type": case_type,
                "solution_index": int(s.get("solution_index", 0)),
                "validity_disposition": disposition,
                "solver_status": solver_st,
                "selected_reaction_ids": ";".join(r_ids),
                "n_added_rxns": len(r_ids),
                "total_cost": float(s.get("total_cost", math.nan)),
                "raw_milp_objective": float(s.get("raw_milp_objective", math.nan)),
                "mip_gap": float(s.get("mip_gap", math.nan)),
                "dual_bound": float(s.get("dual_bound", math.nan)),
                "verified_growth": gr_val,
                "stoichiometric_max_residual": res_val,
                "max_bound_violation": bnd_val,
                "residual_bounds_growth_verified": checks_ok,
                "cycle_qc_status": cyc_st,
                "spurious_atp_cycle": spur_atp if spur_atp is not None else math.nan,
                "explicit_reason": str(ens.get("fewer_than_target_reason", "")),
                "code_version": "exp1_sbml_highs_engine_v2_exact_ledger",
            })
    elif opt_sols:
        for s in opt_sols:
            r_ids = sorted(s.get("added_rxns", []))
            res_val = float(s.get("fba_check_max_residual", math.nan))
            bnd_val = float(s.get("fba_check_max_bound_viol", math.nan))
            gr_val = float(s.get("verified_fba_growth", math.nan))
            cyc_st = str(s.get("cycle_qc_status", "rejected"))
            spur_atp = s.get("post_pruning_egc_has_spurious_atp_cycle", True)
            ledger_rows.append({
                "genome": genome,
                "cohort": cohort,
                "scenario_id": scen_id,
                "biomass_recipe_id": bio_id,
                "case_type": case_type,
                "solution_index": int(s.get("optimization_solution_index", -1)),
                "validity_disposition": "rejected_cycle_qc",
                "solver_status": solver_st,
                "selected_reaction_ids": ";".join(r_ids),
                "n_added_rxns": len(r_ids),
                "total_cost": float(s.get("total_cost", math.nan)),
                "raw_milp_objective": float(s.get("raw_milp_objective", math.nan)),
                "mip_gap": float(s.get("mip_gap", math.nan)),
                "dual_bound": float(s.get("dual_bound", math.nan)),
                "verified_growth": gr_val,
                "stoichiometric_max_residual": res_val,
                "max_bound_violation": bnd_val,
                "residual_bounds_growth_verified": False,
                "cycle_qc_status": cyc_st,
                "spurious_atp_cycle": spur_atp if spur_atp is not None else True,
                "explicit_reason": str(ens.get("fewer_than_target_reason", "all_solutions_rejected_by_cycle_qc")),
                "code_version": "exp1_sbml_highs_engine_v2_exact_ledger",
            })
    else:
        ledger_rows.append({
            "genome": genome,
            "cohort": cohort,
            "scenario_id": scen_id,
            "biomass_recipe_id": bio_id,
            "case_type": case_type,
            "solution_index": -1,
            "validity_disposition": f"unresolved_{solver_st}",
            "solver_status": solver_st,
            "selected_reaction_ids": "",
            "n_added_rxns": math.nan,
            "total_cost": math.nan,
            "raw_milp_objective": math.nan,
            "mip_gap": float(ens.get("first_solve_mip_gap", math.nan)),
            "dual_bound": float(ens.get("first_solve_dual_bound", math.nan)),
            "verified_growth": math.nan,
            "stoichiometric_max_residual": math.nan,
            "max_bound_violation": math.nan,
            "residual_bounds_growth_verified": False,
            "cycle_qc_status": "none_evaluated",
            "spurious_atp_cycle": math.nan,
            "explicit_reason": str(ens.get("fewer_than_target_reason", solver_st)),
            "code_version": "exp1_sbml_highs_engine_v2_exact_ledger",
        })

    return {
        "summary_row": summary_row,
        "ledger_rows": ledger_rows,
    }


def validate_focused_release_results(
    summary_df: pd.DataFrame,
    ledger_df: pd.DataFrame,
    expected_n_cases: int | None = None,
) -> dict[str, Any]:
    """Audit focused release validation results (`EXP1-FRESH-HANDOFF.txt` & `EXP1-LEDGER-FINAL-ACCEPTANCE.txt`):
    1. Zero invalid accepted solutions: every accepted row must have finite verified residual <= 1e-5,
       bound violation <= 1e-5, verified growth >= 0.00999, non-negative cost, valid added rxns count, distinct solutions > 0.
    2. Non-accepted outcomes must have NA added counts and NA costs, with zero valid solutions.
    3. Ground-truth controls and historical comparisons:
       - Positive controls (representative_feasible or synthetic_positive) must solve feasible and verified.
       - Synthetic known-impossible negative control (if present) must be proven infeasible.
       - Historical comparison cases: old numerical solver labels are comparison only; models under fixed medium remained
         proven infeasible, while models under native predicted medium that solved with HiGHS are verified witnesses.
    4. Per-solution ledger matches summary:
       - Ledger must NOT be empty when summary is non-empty.
       - Exact task keys (genome, scenario_id, biomass_recipe_id) match summary.
       - Accepted row count per task in ledger exactly matches n_distinct_valid_solutions in summary.
       - All accepted rows satisfy finite residual, bounds, growth, and explicit cycle-QC pass.
       - Non-accepted rows in ledger have residual_bounds_growth_verified == False.
    """
    n_cases = len(summary_df)
    n_cases_ok = bool(n_cases == expected_n_cases) if expected_n_cases is not None else bool(n_cases > 0)

    accepted_mask = summary_df["solver_status"].isin(REVALIDATED_FEASIBLE_STATUSES)
    acc_df = summary_df[accepted_mask]
    non_acc_df = summary_df[~accepted_mask]

    # Required verification columns check
    missing_cols = [c for c in REQUIRED_VERIFICATION_COLUMNS if c not in summary_df.columns]
    cols_ok = len(missing_cols) == 0

    acc_valid = True
    if not acc_df.empty and cols_ok:
        d_ok = (acc_df["n_distinct_valid_solutions"].notna()) & (acc_df["n_distinct_valid_solutions"] > 0)
        a_ok = (acc_df["min_added_rxns_count"].notna()) & (acc_df["min_added_rxns_count"] >= 0)
        c_ok = (acc_df["min_total_cost"].notna()) & (acc_df["min_total_cost"] >= 0.0)
        g_ok = (acc_df["min_solution_verified_growth"].notna()) & np.isfinite(acc_df["min_solution_verified_growth"]) & (acc_df["min_solution_verified_growth"] >= 0.00999)
        r_ok = (acc_df["max_solution_residual"].notna()) & np.isfinite(acc_df["max_solution_residual"]) & (acc_df["max_solution_residual"] <= 1e-5)
        b_ok = (acc_df["max_solution_bound_viol"].notna()) & np.isfinite(acc_df["max_solution_bound_viol"]) & (acc_df["max_solution_bound_viol"] <= 1e-5)
        acc_valid = bool(d_ok.all() and a_ok.all() and c_ok.all() and g_ok.all() and r_ok.all() and b_ok.all())
    elif not acc_df.empty and not cols_ok:
        acc_valid = False

    non_acc_valid = True
    if not non_acc_df.empty:
        na_dist = bool((non_acc_df["n_distinct_valid_solutions"].fillna(0) == 0).all()) if "n_distinct_valid_solutions" in non_acc_df.columns else True
        na_add = bool(non_acc_df["min_added_rxns_count"].isna().all()) if "min_added_rxns_count" in non_acc_df.columns else True
        na_cost = bool(non_acc_df["min_total_cost"].isna().all()) if "min_total_cost" in non_acc_df.columns else True
        non_acc_valid = bool(na_dist and na_add and na_cost)

    # Per-solution ledger integrity and correspondence
    ledger_valid = True
    if ledger_df is None or (ledger_df.empty and not summary_df.empty):
        ledger_valid = False
    elif not ledger_df.empty:
        # Key correspondence
        sum_keys = set(zip(summary_df["genome"].astype(str), summary_df["scenario_id"].astype(str), summary_df["biomass_recipe_id"].astype(str)))
        led_keys = set(zip(ledger_df["genome"].astype(str), ledger_df["scenario_id"].astype(str), ledger_df["biomass_recipe_id"].astype(str)))
        keys_match = bool(sum_keys == led_keys)

        # Count correspondence
        acc_l = ledger_df[ledger_df["validity_disposition"] == "accepted_feasible"]
        acc_counts = acc_l.groupby(["genome", "scenario_id", "biomass_recipe_id"]).size().to_dict() if not acc_l.empty else {}
        counts_match = True
        for _, srow in summary_df.iterrows():
            k = (str(srow["genome"]), str(srow["scenario_id"]), str(srow["biomass_recipe_id"]))
            exp_count = int(srow["n_distinct_valid_solutions"]) if srow["solver_status"] in REVALIDATED_FEASIBLE_STATUSES and pd.notna(srow["n_distinct_valid_solutions"]) else 0
            act_count = acc_counts.get(k, 0)
            if exp_count != act_count:
                counts_match = False
                break

        # Accepted rows verification
        acc_rows_valid = True
        if not acc_l.empty:
            l_res = (acc_l["stoichiometric_max_residual"].notna()) & np.isfinite(acc_l["stoichiometric_max_residual"]) & (acc_l["stoichiometric_max_residual"] >= 0.0) & (acc_l["stoichiometric_max_residual"] <= 1e-5)
            l_bnd = (acc_l["max_bound_violation"].notna()) & np.isfinite(acc_l["max_bound_violation"]) & (acc_l["max_bound_violation"] >= 0.0) & (acc_l["max_bound_violation"] <= 1e-5)
            l_gr = (acc_l["verified_growth"].notna()) & np.isfinite(acc_l["verified_growth"]) & (acc_l["verified_growth"] >= 0.00999)
            l_v = acc_l["residual_bounds_growth_verified"] == True
            l_cyc = acc_l["cycle_qc_status"].isin({"passed_closed_uptake_egc_check", "passed_no_free_atp_flux", "evaluated_no_atp_species", "not_required_no_atp_species"})
            l_spur = acc_l["spurious_atp_cycle"] != True
            acc_rows_valid = bool(l_res.all() and l_bnd.all() and l_gr.all() and l_v.all() and l_cyc.all() and l_spur.all())
        elif len(acc_df) > 0:
            acc_rows_valid = False

        # Non-accepted rows verification
        non_acc_l = ledger_df[ledger_df["validity_disposition"] != "accepted_feasible"]
        non_acc_rows_valid = True
        if not non_acc_l.empty:
            non_acc_rows_valid = bool((non_acc_l["residual_bounds_growth_verified"] == False).all())

        ledger_valid = bool(keys_match and counts_match and acc_rows_valid and non_acc_rows_valid)

    # Positive controls
    pos_ctrl = summary_df[summary_df["case_type"].isin(["representative_feasible", "synthetic_positive_control"])]
    pos_ctrl_ok = bool(not pos_ctrl.empty and pos_ctrl["solver_status"].isin(REVALIDATED_FEASIBLE_STATUSES).all())

    # Negative controls and historical comparisons
    # Item 5 of EXP1-LEDGER-FINAL-ACCEPTANCE.txt:
    # Synthetic known-impossible negative control provides ground truth;
    # historical infeasible labels provide comparison only (old solver labels are not ground truth negatives).
    syn_neg = summary_df[summary_df["case_type"] == "known_impossible_negative_control"]
    syn_neg_ok = bool(syn_neg["solver_status"].isin(PROVEN_INFEASIBLE_STATUSES).all()) if not syn_neg.empty else True

    hist_infeas = summary_df[summary_df["case_type"].isin(["representative_infeasible", "historical_infeasible_comparison"])]
    if not hist_infeas.empty:
        # Each historical case must be either proven infeasible under fixed medium or verified feasible under native medium
        fixed_med_cases = hist_infeas[hist_infeas["scenario_id"] != "S0_gapseq_native_predicted"]
        native_cases = hist_infeas[hist_infeas["scenario_id"] == "S0_gapseq_native_predicted"]
        fixed_ok = bool(fixed_med_cases["solver_status"].isin(PROVEN_INFEASIBLE_STATUSES).all()) if not fixed_med_cases.empty else True
        native_ok = bool(native_cases["solver_status"].isin(REVALIDATED_FEASIBLE_STATUSES).all()) if not native_cases.empty else True
        hist_ok = bool(fixed_ok and native_ok)
    else:
        hist_ok = True

    neg_ctrl_ok = bool(syn_neg_ok and hist_ok)

    # Statuses honest
    actual_st = set(summary_df["solver_status"].unique())
    statuses_ok = bool(actual_st.issubset(KNOWN_STATUSES) and len(actual_st) > 0)

    checks = {
        "expected_case_count_matched": n_cases_ok,
        "required_verification_columns_present": cols_ok,
        "all_accepted_solutions_strictly_verified": acc_valid,
        "non_accepted_models_clean_na_costs_and_counts": non_acc_valid,
        "per_solution_ledger_integrity_verified": ledger_valid,
        "representative_positive_controls_feasible": pos_ctrl_ok,
        "negative_controls_and_historical_comparisons_sound": neg_ctrl_ok,
        "statuses_honest_and_known": statuses_ok,
    }
    validation_passed = all(checks.values())

    status_counts = summary_df["solver_status"].value_counts().to_dict()
    case_breakdowns = summary_df.groupby(["case_type", "solver_status"]).size().unstack(fill_value=0).to_dict(orient="index")

    return {
        "validation_passed": validation_passed,
        "n_cases_evaluated": len(summary_df),
        "n_accepted_feasible": int(len(acc_df)),
        "n_non_accepted": int(len(non_acc_df)),
        "status_counts": status_counts,
        "case_breakdowns": case_breakdowns,
        "checks": checks,
    }


def main() -> None:
    ap = argparse.ArgumentParser(description="Focused Release Validation Job on Shaheen.")
    ap.add_argument("--old-summary-tsv", type=Path, required=True)
    ap.add_argument("--pool-pkl", type=Path, required=True)
    ap.add_argument("--gapseq-out-dir", type=Path, required=True)
    ap.add_argument("--supp-models-dir", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--milp-time-limit-s", type=float, default=10.0)
    ap.add_argument("--workers", type=int, default=32)
    args = ap.parse_args()

    t0 = time.time()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    old_df = pd.read_csv(args.old_summary_tsv, sep="\t")
    print(f"[{time.strftime('%X')}] Loaded {len(old_df)} reference summary rows from {args.old_summary_tsv}")

    # 1. Select the 68 numerical failures (64 numerical_coupling_violation + 4 fba_verification_failed)
    failed_rows = old_df[old_df["solver_status"].isin(["numerical_coupling_violation", "fba_verification_failed"])].copy()
    failed_rows["case_type"] = "previously_problematic_numerical_failure"

    # 2. Select 8 representative previously-feasible cases (4 optimal, 4 heuristic_incumbent_bounded)
    feas_rows = old_df[old_df["solver_status"].isin(["optimal", "heuristic_incumbent_bounded"])].copy()
    rep_feas = pd.concat([
        feas_rows[feas_rows["solver_status"] == "optimal"].drop_duplicates("genome").head(4),
        feas_rows[feas_rows["solver_status"] == "heuristic_incumbent_bounded"].drop_duplicates("genome").head(4),
    ], ignore_index=True)
    rep_feas["case_type"] = "representative_feasible"

    # 3. Select 8 representative proven-infeasible cases (4 infeasible, 4 infeasible_under_fixed_medium)
    infeas_rows = old_df[old_df["solver_status"].isin(["infeasible", "infeasible_under_fixed_medium"])].copy()
    rep_infeas = pd.concat([
        infeas_rows[infeas_rows["solver_status"] == "infeasible"].drop_duplicates("genome").head(4),
        infeas_rows[infeas_rows["solver_status"] == "infeasible_under_fixed_medium"].drop_duplicates("genome").head(4),
    ], ignore_index=True)
    rep_infeas["case_type"] = "representative_infeasible"

    selected_tasks_df = pd.concat([failed_rows, rep_feas, rep_infeas], ignore_index=True).drop_duplicates(
        subset=["genome", "scenario_id", "biomass_recipe_id"]
    ).reset_index(drop=True)

    print(f"[{time.strftime('%X')}] Selected {len(selected_tasks_df)} validation cases:")
    print(selected_tasks_df.groupby(["case_type", "solver_status"]).size().to_string())

    import pickle
    with open(args.pool_pkl, "rb") as fh:
        pool_data = pickle.load(fh)
    universal_rxns = pool_data["reactions"]

    tasks: list[dict[str, Any]] = []
    for _, r in selected_tasks_df.iterrows():
        tasks.append({
            "genome": str(r["genome"]),
            "cohort": str(r.get("cohort", "unknown")),
            "canonical_order": str(r.get("canonical_order", "")),
            "canonical_family": str(r.get("canonical_family", "")),
            "scenario_id": str(r["scenario_id"]),
            "biomass_recipe_id": str(r["biomass_recipe_id"]),
            "old_status": str(r["solver_status"]),
            "case_type": str(r["case_type"]),
            "gapseq_out_dir": str(args.gapseq_out_dir),
            "supp_models_dir": str(args.supp_models_dir),
            "universal_rxns": universal_rxns,
            "milp_time_limit_s": args.milp_time_limit_s,
        })

    print(f"[{time.strftime('%X')}] Launching {len(tasks)} validation tasks across {args.workers} worker processes...")
    all_summary: list[dict[str, Any]] = []
    all_ledger: list[dict[str, Any]] = []

    with ProcessPoolExecutor(max_workers=min(args.workers, max(1, len(tasks)))) as pool:
        futures = {pool.submit(run_single_validation_task, t): t for t in tasks}
        for fut in as_completed(futures):
            res = fut.result()
            all_summary.append(res["summary_row"])
            all_ledger.extend(res["ledger_rows"])
            s_row = res["summary_row"]
            print(f"  [DONE] {s_row['genome']} {s_row['scenario_id']} {s_row['biomass_recipe_id']}: {s_row['old_status']} -> {s_row['solver_status']} (residual: {s_row['max_solution_residual']:.1e}, runtime: {s_row['runtime_seconds']:.1f}s)", flush=True)

    summary_df = pd.DataFrame(all_summary).sort_values(["case_type", "genome", "scenario_id", "biomass_recipe_id"]).reset_index(drop=True)
    ledger_df = pd.DataFrame(all_ledger).sort_values(["case_type", "genome", "scenario_id", "biomass_recipe_id", "solution_index"]).reset_index(drop=True)

    summary_tsv = args.out_dir / "exp1_focused_validation_summary.tsv"
    ledger_tsv = args.out_dir / "exp1_focused_validation_per_solution_ledger.tsv"
    report_json = args.out_dir / "exp1_focused_validation_report.json"

    summary_df.to_csv(summary_tsv, sep="\t", index=False)
    ledger_df.to_csv(ledger_tsv, sep="\t", index=False)

    report = validate_focused_release_results(summary_df, ledger_df)
    report["elapsed_seconds"] = round(time.time() - t0, 2)

    with open(report_json, "w", encoding="utf-8") as fh:
        json.dump(report, fh, indent=2, sort_keys=True)

    print(f"\n[{time.strftime('%X')}] Focused Release Validation Report:")
    print(json.dumps(report, indent=2))

    if not report["validation_passed"]:
        print(f"[ERROR] Focused release validation FAILED checks: {report['checks']}", file=sys.stderr)
        sys.exit(2)

    print(f"\n[SUCCESS] Focused release validation PASSED across all {len(summary_df)} cases in {report['elapsed_seconds']}s.")


if __name__ == "__main__":
    main()
