#!/usr/bin/env python3
"""Stage 2 Focused Diagnostic Retry Script.
Executes the exact 68 combinations that failed with numerical_coupling_violation or
fba_verification_failed in run_15873653 to verify that tightened HiGHS tolerances (1e-9)
and native reaction bounds eliminate numerical violations without arbitrary flux caps.
Per EXP1-QC-DESIGN-CORRECTION.txt: Runs in Slurm, records differentiated rates.
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
    build_genome_candidate_weights,
    load_genome_archive_payload,
)


def run_single_retry_task(task: dict[str, Any]) -> dict[str, Any]:
    genome = task["genome"]
    scen_id = task["scenario_id"]
    bio_id = task["biomass_recipe_id"]
    old_status = task["old_status"]
    gapseq_out_dir = Path(task["gapseq_out_dir"])
    supp_models_dir = Path(task["supp_models_dir"])
    universal_rxns = task["universal_rxns"]
    milp_time_limit_s = float(task.get("milp_time_limit_s", 10.0))

    payload = load_genome_archive_payload(genome, gapseq_out_dir, supp_models_dir)
    draft_mod = payload["draft_model"]
    cand_weights = build_genome_candidate_weights(universal_rxns, payload["rxn_evi"], payload["trans_evi"])

    # Determine scenario compounds
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
    v_growth = sols[0]["verified_fba_growth"] if sols else math.nan

    return {
        "genome": genome,
        "scenario_id": scen_id,
        "biomass_recipe_id": bio_id,
        "old_status": old_status,
        "new_solver_status": ens["solver_status"],
        "new_optimizer_status": ens.get("optimizer_status", "unknown"),
        "raw_solver_status": ens.get("raw_solver_status", -1),
        "fewer_than_target_reason": ens.get("fewer_than_target_reason", ""),
        "n_distinct_valid_solutions": len(sols),
        "n_optimization_feasible": ens.get("n_optimization_feasible_solutions", 0),
        "n_cycle_qc_rejected": ens.get("n_cycle_qc_rejected_solutions", 0),
        "verified_fba_growth": v_growth,
        "first_solve_dual_bound": ens.get("first_solve_dual_bound", math.nan),
        "first_solve_mip_gap": ens.get("first_solve_mip_gap", math.nan),
        "runtime_seconds": dt,
    }


def main() -> None:
    ap = argparse.ArgumentParser(description="Focused diagnostic retry for exp1 68 numerical failure rows.")
    ap.add_argument("--old-summary-tsv", type=Path, required=True)
    ap.add_argument("--pool-pkl", type=Path, required=True)
    ap.add_argument("--gapseq-out-dir", type=Path, required=True)
    ap.add_argument("--supp-models-dir", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--milp-time-limit-s", type=float, default=10.0)
    ap.add_argument("--workers", type=int, default=32)
    args = ap.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    old_df = pd.read_csv(args.old_summary_tsv, sep="\t")
    failed_rows = old_df[old_df["solver_status"].isin(["numerical_coupling_violation", "fba_verification_failed"])].copy()

    print(f"[{time.strftime('%X')}] Loaded {len(failed_rows)} target failing rows from {args.old_summary_tsv}")
    print(f"Status breakdown: {failed_rows['solver_status'].value_counts().to_dict()}")

    import pickle
    with open(args.pool_pkl, "rb") as fh:
        pool_data = pickle.load(fh)
    universal_rxns = pool_data["reactions"]

    tasks = []
    for _, r in failed_rows.iterrows():
        tasks.append({
            "genome": r["genome"],
            "scenario_id": r["scenario_id"],
            "biomass_recipe_id": r["biomass_recipe_id"],
            "old_status": r["solver_status"],
            "gapseq_out_dir": str(args.gapseq_out_dir),
            "supp_models_dir": str(args.supp_models_dir),
            "universal_rxns": universal_rxns,
            "milp_time_limit_s": args.milp_time_limit_s,
        })

    print(f"[{time.strftime('%X')}] Launching {len(tasks)} retry tasks across {args.workers} worker processes...")
    results = []
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(run_single_retry_task, t): t for t in tasks}
        for fut in as_completed(futures):
            res = fut.result()
            results.append(res)
            print(f"  [DONE] {res['genome']} {res['scenario_id']} {res['biomass_recipe_id']}: {res['old_status']} -> {res['new_solver_status']} (runtime: {res['runtime_seconds']:.1f}s)")

    res_df = pd.DataFrame(results)
    summary_path = args.out_dir / "exp1_diagnostic_retry_summary.tsv"
    res_df.to_csv(summary_path, sep="\t", index=False)
    print(f"\nWrote summary to {summary_path}")

    # Compute differentiated rates
    new_status_counts = res_df["new_solver_status"].value_counts().to_dict()
    revalidated_feasible = sum(new_status_counts.get(s, 0) for s in ["optimal", "heuristic_incumbent_bounded", "optimal_after_cycle_qc_filtering"])
    proven_infeasible = sum(new_status_counts.get(s, 0) for s in ["infeasible", "infeasible_under_fixed_medium"])
    numerical_failures = sum(new_status_counts.get(s, 0) for s in ["numerical_coupling_violation", "fba_verification_failed", "numerical_error"])
    timeouts = sum(new_status_counts.get(s, 0) for s in ["time_limit_no_incumbent", "time_limit_heuristic_incumbent"])

    report = {
        "n_retried": len(res_df),
        "status_distribution_old": failed_rows["solver_status"].value_counts().to_dict(),
        "status_distribution_new": new_status_counts,
        "n_revalidated_feasible": revalidated_feasible,
        "rate_revalidated_feasible": revalidated_feasible / len(res_df),
        "n_proven_infeasible": proven_infeasible,
        "rate_proven_infeasible": proven_infeasible / len(res_df),
        "n_numerical_failures": numerical_failures,
        "rate_numerical_failures": numerical_failures / len(res_df),
        "n_timeouts": timeouts,
        "rate_timeouts": timeouts / len(res_df),
        "remediation_success": numerical_failures == 0,
    }

    report_path = args.out_dir / "exp1_diagnostic_retry_report.json"
    with open(report_path, "w") as fh:
        json.dump(report, fh, indent=2)

    print("\n" + "=" * 60)
    print("EXP1 DIAGNOSTIC RETRY RESULTS:")
    print(f"Total retried: {len(res_df)}")
    print(f"New Status Distribution: {new_status_counts}")
    print(f"Revalidated Feasible: {revalidated_feasible} ({revalidated_feasible / len(res_df) * 100:.1f}%)")
    print(f"Proven Infeasible: {proven_infeasible} ({proven_infeasible / len(res_df) * 100:.1f}%)")
    print(f"Numerical Failures: {numerical_failures} ({numerical_failures / len(res_df) * 100:.1f}%)")
    print(f"Timeouts: {timeouts} ({timeouts / len(res_df) * 100:.1f}%)")
    print("=" * 60)


if __name__ == "__main__":
    main()
