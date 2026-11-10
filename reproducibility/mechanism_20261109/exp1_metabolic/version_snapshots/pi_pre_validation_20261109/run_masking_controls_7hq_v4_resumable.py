#!/usr/bin/env python3
"""Resumable, Checkpointed Masking Controls Runner for the 7 High-Quality External Reference Genomes (v4).
Prepared per `EXP1-V4-REVIEW-CORRECTIONS.txt` and `EXP1-CONFIRMED-FORK-SPIN-REPAIR.txt`:
1. Uses multiprocessing 'spawn' context exclusively to prevent HiGHS thread pool deadlocks across fork.
2. Parallelizes baseline solves (Step 1) across workers instead of serial execution.
3. Checkpoints baselines, media provenance, masked runs, and per-solution ledger rows atomically
   with an exact SHA-256 commit marker JSON to guarantee cross-file consistency.
4. Resumes ONLY when input hash, solver config hash, engine code SHA-256, and per-case ledger rows
   match the commit marker; otherwise recomputes explicitly.
5. Strictly consumes exact engine solution keys (added_rxns, total_cost, raw_milp_objective,
   verified_fba_growth, fba_check_max_residual, fba_check_max_bound_viol, mip_gap, dual_bound,
   post_pruning_egc_has_spurious_atp_cycle, cycle_qc_status, science_qc_eligible) and rejects missing fields.
6. Derives validity_disposition honestly without invented QC passes.
7. Preserves fresh v4 output namespace and isolates old v3 outputs.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import multiprocessing as mp
import os
import pickle
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

CODE_DIR = Path(__file__).resolve().parent
if str(CODE_DIR) not in sys.path:
    sys.path.insert(0, str(CODE_DIR))

# Import from versioned v4 engine with explicit threads=1 across all HiGHS invocations
try:
    from exp1_sbml_highs_engine_v4 import (
        BIO_B2_OMIT_7_KEEP_F420,
        NUTRIENT_SCENARIOS,
        SCEN_S4_OXIC_DESERT,
        SBMLModelData,
        apply_artificial_gene_masking,
        enumerate_alternative_gapfill_ensemble,
    )
except ImportError:
    from exp1_sbml_highs_engine import (
        BIO_B2_OMIT_7_KEEP_F420,
        NUTRIENT_SCENARIOS,
        SCEN_S4_OXIC_DESERT,
        SBMLModelData,
        apply_artificial_gene_masking,
        enumerate_alternative_gapfill_ensemble,
    )

from stage2_pilot24_ensemble import (
    build_genome_candidate_weights,
    load_genome_archive_payload,
)

EXACT_7_HQ_EXTERNAL_IDS = [
    "MGYG000519815",
    "MGYG000520046",
    "MGYG000519748",
    "MGYG000510876",
    "MGYG000495096",
    "MGYG000526805",
    "MGYG000523080",
]

REVALIDATED_FEASIBLE_STATUSES = {
    "optimal",
    "optimal_after_cycle_qc_filtering",
    "optimal_draft_no_gapfill_needed",
    "heuristic_incumbent_bounded",
    "post_pruning_cost_changed_not_proven_optimal",
}

REQUIRED_ENGINE_SOL_KEYS = (
    "added_rxns",
    "total_cost",
    "raw_milp_objective",
    "verified_fba_growth",
    "fba_check_max_residual",
    "fba_check_max_bound_viol",
    "mip_gap",
    "dual_bound",
    "post_pruning_egc_has_spurious_atp_cycle",
    "cycle_qc_status",
    "science_qc_eligible",
)


def derive_stable_masking_seed(genome: str, mask_fraction: float, replicate: int) -> tuple[int, str]:
    """Derive a stable, architecture- and process-independent 32-bit integer seed via SHA-256."""
    seed_str = f"exp1_masking_{genome}_{int(round(mask_fraction * 100))}_{replicate}"
    seed_val = int(hashlib.sha256(seed_str.encode("utf-8")).hexdigest()[:8], 16)
    return seed_val, seed_str


def compute_file_sha256(path: Path) -> str:
    """Compute sha256 checksum of an existing file."""
    if not path.exists():
        return ""
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while chunk := fh.read(65536):
            h.update(chunk)
    return h.hexdigest()


def serialize_and_hash_medium(compounds: dict[str, tuple[str, float]]) -> tuple[str, str]:
    """Canonical JSON serialization and SHA-256 of medium compounds without arbitrary rounding."""
    sorted_items = {cid: [str(cname), float(ub)] for cid, (cname, ub) in sorted(compounds.items())}
    ser_json = json.dumps(sorted_items, sort_keys=True)
    med_hash = hashlib.sha256(ser_json.encode("utf-8")).hexdigest()
    return ser_json, med_hash


def compute_exact_input_provenance_hash(
    archive_sha256: str,
    draft_model_id: str,
    medium_hash: str,
    candidate_pool_sha256: str,
    biomass_recipe_id: str,
) -> str:
    """Compute exact input provenance hash from raw source archive SHA, model ID, medium hash,
    candidate pool SHA, and biomass recipe ID.
    """
    h = hashlib.sha256()
    h.update(archive_sha256.encode("utf-8"))
    h.update(draft_model_id.encode("utf-8"))
    h.update(medium_hash.encode("utf-8"))
    h.update(candidate_pool_sha256.encode("utf-8"))
    h.update(biomass_recipe_id.encode("utf-8"))
    return h.hexdigest()


def compute_case_config_hash(
    input_hash: str,
    solver_config_hash: str,
    engine_sha256: str,
    condition_id: str,
) -> str:
    """Compute exact per-case configuration hash combining inputs, solver parameters, and engine code."""
    h = hashlib.sha256()
    h.update(input_hash.encode("utf-8"))
    h.update(solver_config_hash.encode("utf-8"))
    h.update(engine_sha256.encode("utf-8"))
    h.update(condition_id.encode("utf-8"))
    return h.hexdigest()


def build_solution_ledger_row(
    sol: dict[str, Any],
    condition_type: str,
    genome: str,
    cohort: str,
    canonical_order: str,
    canonical_family: str,
    mask_fraction: float,
    replicate: int,
    solution_index: int,
    solver_status: str,
    code_version: str = "v4_masking_controls_resumable",
) -> dict[str, Any]:
    """Build a validated per-solution ledger row from an actual engine solution dict.
    Strictly validates that all REQUIRED_ENGINE_SOL_KEYS are present and rejects missing fields.
    Derives validity_disposition honestly without invented QC passes.
    """
    missing = [k for k in REQUIRED_ENGINE_SOL_KEYS if k not in sol]
    if missing:
        raise KeyError(f"Engine solution dictionary missing required key(s): {missing}")

    added_rxns = sol["added_rxns"]
    if not isinstance(added_rxns, (list, set, tuple)):
        raise TypeError(f"added_rxns must be a list/set/tuple, got {type(added_rxns)}")
    rids_str = ",".join(sorted(added_rxns))
    n_added = len(added_rxns)

    total_cost = float(sol["total_cost"])
    raw_obj = float(sol["raw_milp_objective"])
    raw_gap = sol["mip_gap"]
    mip_gap = float(raw_gap) if (raw_gap is not None and not (isinstance(raw_gap, float) and math.isnan(raw_gap))) else math.nan
    raw_db = sol["dual_bound"]
    dual_bnd = float(raw_db) if (raw_db is not None and not (isinstance(raw_db, float) and math.isnan(raw_db))) else math.nan

    v_growth = float(sol["verified_fba_growth"])
    res_val = float(sol["fba_check_max_residual"])
    bnd_val = float(sol["fba_check_max_bound_viol"])
    has_spurious = bool(sol["post_pruning_egc_has_spurious_atp_cycle"])
    cyc_status = str(sol["cycle_qc_status"])
    sci_eligible = bool(sol["science_qc_eligible"])

    residual_bounds_growth_verified = bool(
        math.isfinite(res_val) and res_val <= 1e-5
        and math.isfinite(bnd_val) and bnd_val <= 1e-5
        and math.isfinite(v_growth) and v_growth >= 0.00999
    )

    if (
        sci_eligible
        and not has_spurious
        and residual_bounds_growth_verified
    ):
        validity_disposition = "accepted_feasible"
        explicit_reason = ""
    else:
        if has_spurious:
            validity_disposition = "rejected_spurious_atp_cycle"
            explicit_reason = f"Post-pruning EGC check detected spurious ATP cycle ({cyc_status})"
        elif not math.isfinite(v_growth) or v_growth < 0.00999:
            validity_disposition = "rejected_insufficient_growth"
            explicit_reason = f"Verified FBA growth {v_growth} < 0.01"
        elif not math.isfinite(res_val) or res_val > 1e-5:
            validity_disposition = "rejected_residual_violation"
            explicit_reason = f"FBA check max residual {res_val} > 1e-5"
        elif not math.isfinite(bnd_val) or bnd_val > 1e-5:
            validity_disposition = "rejected_bound_violation"
            explicit_reason = f"FBA check max bound viol {bnd_val} > 1e-5"
        else:
            validity_disposition = "rejected_cycle_qc_ineligible"
            explicit_reason = f"Solution rejected by cycle QC: {cyc_status}"

    return {
        "condition_type": condition_type,
        "genome": genome,
        "cohort": cohort,
        "canonical_order": canonical_order,
        "canonical_family": canonical_family,
        "scenario_id": SCEN_S4_OXIC_DESERT,
        "biomass_recipe_id": BIO_B2_OMIT_7_KEEP_F420,
        "mask_fraction": mask_fraction,
        "replicate": replicate,
        "solution_index": solution_index,
        "validity_disposition": validity_disposition,
        "solver_status": solver_status,
        "selected_reaction_ids": rids_str,
        "n_added_rxns": n_added,
        "total_cost": total_cost,
        "raw_milp_objective": raw_obj,
        "mip_gap": mip_gap,
        "dual_bound": dual_bnd,
        "verified_growth": v_growth,
        "stoichiometric_max_residual": res_val,
        "max_bound_violation": bnd_val,
        "residual_bounds_growth_verified": residual_bounds_growth_verified,
        "cycle_qc_status": cyc_status,
        "spurious_atp_cycle": has_spurious,
        "explicit_reason": explicit_reason,
        "code_version": code_version,
    }


def build_empty_solver_ledger_row(
    condition_type: str,
    genome: str,
    cohort: str,
    canonical_order: str,
    canonical_family: str,
    mask_fraction: float,
    replicate: int,
    solver_status: str,
    code_version: str = "v4_masking_controls_resumable",
) -> dict[str, Any]:
    """Build a ledger row representing a failed or infeasible solver outcome."""
    return {
        "condition_type": condition_type,
        "genome": genome,
        "cohort": cohort,
        "canonical_order": canonical_order,
        "canonical_family": canonical_family,
        "scenario_id": SCEN_S4_OXIC_DESERT,
        "biomass_recipe_id": BIO_B2_OMIT_7_KEEP_F420,
        "mask_fraction": mask_fraction,
        "replicate": replicate,
        "solution_index": 0,
        "validity_disposition": f"rejected_{solver_status}",
        "solver_status": solver_status,
        "selected_reaction_ids": "",
        "n_added_rxns": 0,
        "total_cost": math.nan,
        "raw_milp_objective": math.nan,
        "mip_gap": math.nan,
        "dual_bound": math.nan,
        "verified_growth": math.nan,
        "stoichiometric_max_residual": math.nan,
        "max_bound_violation": math.nan,
        "residual_bounds_growth_verified": False,
        "cycle_qc_status": "not_applicable_solver_failed",
        "spurious_atp_cycle": math.nan,
        "explicit_reason": f"Solver returned {solver_status}",
        "code_version": code_version,
    }


def evaluate_single_baseline_case(task: dict[str, Any]) -> dict[str, Any]:
    """Evaluate one matched baseline (0% mask) under per-genome S4 medium and B2 biomass."""
    genome = task["genome"]
    cohort = task["cohort"]
    canon_order = task.get("canonical_order", "")
    canon_family = task.get("canonical_family", "")
    completeness = task["completeness"]
    contamination = task["contamination"]
    draft_mod: SBMLModelData = task["draft_model"]
    s4_compounds = task["s4_compounds"]
    medium_hash = task["medium_hash"]
    archive_sha = task["archive_sha256"]
    cand_rxns = task["candidate_reactions"]
    cand_weights = task["candidate_weights"]
    milp_time_limit_s = task["milp_time_limit_s"]
    input_hash = task["input_hash"]
    case_config_hash = task.get("case_config_hash", "")

    t0 = time.time()
    base_ens = enumerate_alternative_gapfill_ensemble(
        draft_model=draft_mod,
        candidate_reactions=cand_rxns,
        candidate_weights=cand_weights,
        scenario_compounds=s4_compounds,
        biomass_recipe_id=BIO_B2_OMIT_7_KEEP_F420,
        min_growth=0.01,
        n_Target_solutions=8,
        cost_slack_frac=0.25,
        candidate_pool_scope="artificial_gene_masking_control_baseline",
        time_limit_per_milp_s=milp_time_limit_s,
    )
    dt = round(time.time() - t0, 2)
    st = str(base_ens["solver_status"])
    nsol = int(base_ens.get("n_distinct_valid_solutions", 0))
    base_valid = bool(st in REVALIDATED_FEASIBLE_STATUSES and nsol > 0)
    gapfill_rxns = sorted(base_ens.get("reaction_frequencies", {}).keys()) if base_valid else []

    summary_row = {
        "condition_type": "baseline",
        "genome": genome,
        "cohort": cohort,
        "canonical_order": canon_order,
        "canonical_family": canon_family,
        "completeness": completeness,
        "contamination": contamination,
        "scenario_id": SCEN_S4_OXIC_DESERT,
        "biomass_recipe_id": BIO_B2_OMIT_7_KEEP_F420,
        "mask_fraction": 0.0,
        "replicate": 0,
        "medium_hash": medium_hash,
        "archive_sha256": archive_sha,
        "input_hash": input_hash,
        "case_config_hash": case_config_hash,
        "solver_status": st,
        "n_valid_solutions": nsol,
        "n_gapfill_rxns": len(gapfill_rxns) if base_valid else math.nan,
        "gapfill_reaction_ids": ",".join(gapfill_rxns),
        "runtime_seconds": dt,
    }

    ledger_rows: list[dict[str, Any]] = []
    # Collect all solutions (both cycle-QC-eligible and optimization-feasible)
    ens_solutions = base_ens.get("optimization_feasible_solutions", base_ens.get("solutions", []))
    if ens_solutions and base_valid:
        for s_idx, sol in enumerate(ens_solutions):
            ledger_rows.append(
                build_solution_ledger_row(
                    sol=sol,
                    condition_type="baseline",
                    genome=genome,
                    cohort=cohort,
                    canonical_order=canon_order,
                    canonical_family=canon_family,
                    mask_fraction=0.0,
                    replicate=0,
                    solution_index=s_idx,
                    solver_status=st,
                )
            )
    else:
        ledger_rows.append(
            build_empty_solver_ledger_row(
                condition_type="baseline",
                genome=genome,
                cohort=cohort,
                canonical_order=canon_order,
                canonical_family=canon_family,
                mask_fraction=0.0,
                replicate=0,
                solver_status=st,
            )
        )

    return {
        "genome": genome,
        "status": st,
        "n_valid_solutions": nsol,
        "gapfill_rxns": gapfill_rxns,
        "summary_row": summary_row,
        "ledger_rows": ledger_rows,
        "runtime_seconds": dt,
    }


def evaluate_single_mask_case(task: dict[str, Any]) -> dict[str, Any]:
    """Evaluate one masking replicate under the exact same S4 medium as baseline."""
    genome = task["genome"]
    cohort = task["cohort"]
    canon_order = task.get("canonical_order", "")
    canon_family = task.get("canonical_family", "")
    completeness = task["completeness"]
    contamination = task["contamination"]
    mask_frac = task["mask_fraction"]
    rep = task["replicate"]
    draft_mod: SBMLModelData = task["draft_model"]
    s4_compounds = task["s4_compounds"]
    medium_hash = task["medium_hash"]
    cand_rxns = task["candidate_reactions"]
    cand_weights = task["candidate_weights"]
    milp_time_limit_s = task["milp_time_limit_s"]
    input_hash = task["input_hash"]
    case_config_hash = task.get("case_config_hash", "")
    base_status = task["baseline_status"]
    base_nsol = task["baseline_n_valid_solutions"]
    base_rxns = set(task["baseline_gapfill_rxns"])

    seed_val, seed_str = derive_stable_masking_seed(genome, mask_frac, rep)
    rng = np.random.default_rng(seed_val)

    masked_draft, mask_meta = apply_artificial_gene_masking(draft_mod, mask_fraction=mask_frac, rng=rng)
    masked_weights = dict(cand_weights)
    for rem_rid in mask_meta["removed_rxn_ids"]:
        masked_weights[rem_rid] = 100.0

    t0 = time.time()
    m_ens = enumerate_alternative_gapfill_ensemble(
        draft_model=masked_draft,
        candidate_reactions=cand_rxns,
        candidate_weights=masked_weights,
        scenario_compounds=s4_compounds,
        biomass_recipe_id=BIO_B2_OMIT_7_KEEP_F420,
        min_growth=0.01,
        n_Target_solutions=8,
        cost_slack_frac=0.25,
        candidate_pool_scope="artificial_gene_masking_control",
        time_limit_per_milp_s=milp_time_limit_s,
    )
    runtime_s = round(time.time() - t0, 2)

    masked_status = str(m_ens["solver_status"])
    masked_nsol = int(m_ens.get("n_distinct_valid_solutions", 0))
    base_valid = bool(base_status in REVALIDATED_FEASIBLE_STATUSES and base_nsol > 0)
    masked_valid = bool(masked_status in REVALIDATED_FEASIBLE_STATUSES and masked_nsol > 0)

    if masked_valid:
        m_rxns = set(m_ens["reaction_frequencies"].keys())
        n_masked_rxns: float = float(len(m_rxns))
    else:
        m_rxns = set()
        n_masked_rxns = math.nan

    if base_valid:
        n_base_rxns: float = float(len(base_rxns))
    else:
        n_base_rxns = math.nan

    is_evaluable = bool(base_valid and masked_valid)
    removed_set = set(mask_meta["removed_rxn_ids"])

    if is_evaluable:
        newly_induced = m_rxns - base_rxns
        recovered_true_masked = newly_induced & removed_set
        alternative_candidate_rxns = newly_induced - removed_set

        n_newly_induced: float = float(len(newly_induced))
        n_recovered: float = float(len(recovered_true_masked))
        n_alternative: float = float(len(alternative_candidate_rxns))

        recovery_rate: float = (n_recovered / len(removed_set)) if len(removed_set) > 0 else math.nan
        alternative_fraction: float = (n_alternative / n_newly_induced) if n_newly_induced > 0 else (0.0 if len(newly_induced) == 0 else math.nan)
        eval_status = "evaluable_valid_ensembles"
    else:
        n_newly_induced = math.nan
        n_recovered = math.nan
        n_alternative = math.nan
        recovery_rate = math.nan
        alternative_fraction = math.nan
        if not base_valid:
            eval_status = f"non_evaluable_invalid_baseline_{base_status}"
        else:
            eval_status = f"non_evaluable_invalid_masked_{masked_status}"

    masked_genes_str = ",".join(sorted(mask_meta.get("masked_gene_ids", [])))
    removed_rxns_str = ",".join(sorted(mask_meta.get("removed_rxn_ids", [])))

    summary_row = {
        "condition_type": "masked",
        "genome": genome,
        "cohort": cohort,
        "canonical_order": canon_order,
        "canonical_family": canon_family,
        "completeness": completeness,
        "contamination": contamination,
        "scenario_id": SCEN_S4_OXIC_DESERT,
        "biomass_recipe_id": BIO_B2_OMIT_7_KEEP_F420,
        "mask_fraction": mask_frac,
        "replicate": rep,
        "seed_val": seed_val,
        "seed_str": seed_str,
        "medium_hash": medium_hash,
        "input_hash": input_hash,
        "case_config_hash": case_config_hash,
        "n_total_genes": mask_meta["n_total_genes"],
        "n_masked_genes": mask_meta["n_masked_genes"],
        "masked_gene_ids": masked_genes_str,
        "n_removed_gpr_reactions": mask_meta["n_removed_gpr_reactions"],
        "removed_rxn_ids": removed_rxns_str,
        "baseline_solver_status": base_status,
        "baseline_n_valid_solutions": base_nsol,
        "n_baseline_gapfill_rxns": n_base_rxns,
        "masked_solver_status": masked_status,
        "masked_n_valid_solutions": masked_nsol,
        "n_masked_gapfill_rxns": n_masked_rxns,
        "is_evaluable": is_evaluable,
        "evaluation_status": eval_status,
        "n_newly_induced_gapfill_rxns": n_newly_induced,
        "n_true_masked_rxns_recovered_in_gapfill": n_recovered,
        "n_alternative_candidate_rxns_induced": n_alternative,
        "masked_rxn_recovery_rate_under_medium": recovery_rate,
        "alternative_candidate_fraction_under_medium": alternative_fraction,
        "runtime_seconds": runtime_s,
    }

    ledger_rows: list[dict[str, Any]] = []
    ens_solutions = m_ens.get("optimization_feasible_solutions", m_ens.get("solutions", []))
    if ens_solutions and masked_valid:
        for s_idx, sol in enumerate(ens_solutions):
            ledger_rows.append(
                build_solution_ledger_row(
                    sol=sol,
                    condition_type="masked",
                    genome=genome,
                    cohort=cohort,
                    canonical_order=canon_order,
                    canonical_family=canon_family,
                    mask_fraction=mask_frac,
                    replicate=rep,
                    solution_index=s_idx,
                    solver_status=masked_status,
                )
            )
    else:
        ledger_rows.append(
            build_empty_solver_ledger_row(
                condition_type="masked",
                genome=genome,
                cohort=cohort,
                canonical_order=canon_order,
                canonical_family=canon_family,
                mask_fraction=mask_frac,
                replicate=rep,
                solver_status=masked_status,
            )
        )

    return {
        "summary_row": summary_row,
        "ledger_rows": ledger_rows,
    }


def validate_masking_controls_results(
    masking_df: pd.DataFrame,
    expected_genomes: set[str],
    ledger_df: pd.DataFrame | None = None,
    media_provenance: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Validate completeness, integrity, and scientific consistency of masking controls."""
    actual_genomes = set(masking_df["genome"].astype(str).unique()) if not masking_df.empty else set()
    genomes_ok = bool(actual_genomes == expected_genomes)

    expected_tuples = {
        (g, frac, rep)
        for g in expected_genomes
        for frac in (0.10, 0.20, 0.30)
        for rep in (0, 1, 2)
    }
    actual_tuples = set(zip(
        masking_df["genome"].astype(str),
        masking_df["mask_fraction"].astype(float).round(2),
        masking_df["replicate"].astype(int),
    )) if not masking_df.empty else set()
    tuples_exact_match = bool(actual_tuples == expected_tuples)

    seeds_derived_ok = True
    if not masking_df.empty:
        for _, row in masking_df.iterrows():
            exp_seed, exp_str = derive_stable_masking_seed(str(row["genome"]), float(row["mask_fraction"]), int(row["replicate"]))
            if int(row["seed_val"]) != exp_seed or str(row["seed_str"]) != exp_str:
                seeds_derived_ok = False
                break

    same_scen = bool((masking_df["scenario_id"] == SCEN_S4_OXIC_DESERT).all()) if not masking_df.empty else False
    same_bio = bool((masking_df["biomass_recipe_id"] == BIO_B2_OMIT_7_KEEP_F420).all()) if not masking_df.empty else False
    same_medium_biomass_ok = bool(same_scen and same_bio)

    med_hash_ok = True
    if media_provenance is not None and not masking_df.empty:
        for _, row in masking_df.iterrows():
            g = str(row["genome"])
            if g in media_provenance:
                exp_h = media_provenance[g]["medium_hash"]
                act_h = str(row.get("medium_hash", ""))
                if act_h != exp_h:
                    med_hash_ok = False
                    break

    non_eval = masking_df[masking_df["is_evaluable"] == False]
    if not non_eval.empty:
        non_eval_clean = bool(
            non_eval["n_newly_induced_gapfill_rxns"].isna().all()
            and non_eval["n_true_masked_rxns_recovered_in_gapfill"].isna().all()
            and non_eval["n_alternative_candidate_rxns_induced"].isna().all()
            and non_eval["masked_rxn_recovery_rate_under_medium"].isna().all()
            and non_eval["alternative_candidate_fraction_under_medium"].isna().all()
        )
    else:
        non_eval_clean = True

    zero_removed = masking_df[masking_df["n_removed_gpr_reactions"] == 0]
    if not zero_removed.empty:
        zero_removed_ok = bool(zero_removed["masked_rxn_recovery_rate_under_medium"].isna().all())
    else:
        zero_removed_ok = True

    eval_rows = masking_df[masking_df["is_evaluable"] == True]
    if not eval_rows.empty:
        rec_alt_sum = eval_rows["n_true_masked_rxns_recovered_in_gapfill"] + eval_rows["n_alternative_candidate_rxns_induced"]
        eval_consistent = bool(np.isclose(eval_rows["n_newly_induced_gapfill_rxns"], rec_alt_sum).all())
    else:
        eval_consistent = True

    ledger_ok = True
    if ledger_df is not None:
        if ledger_df.empty:
            ledger_ok = False
        else:
            acc_l = ledger_df[ledger_df["validity_disposition"] == "accepted_feasible"]
            if not acc_l.empty:
                res_ok = bool(((acc_l["stoichiometric_max_residual"].notna()) & (acc_l["stoichiometric_max_residual"] <= 1e-5)).all())
                bnd_ok = bool(((acc_l["max_bound_violation"].notna()) & (acc_l["max_bound_violation"] <= 1e-5)).all())
                gr_ok = bool(((acc_l["verified_growth"].notna()) & (acc_l["verified_growth"] >= 0.00999)).all())
                spur_ok = bool((acc_l["spurious_atp_cycle"] == False).all())
                ledger_ok = bool(res_ok and bnd_ok and gr_ok and spur_ok)

    checks = {
        "expected_7_hq_genomes_matched": genomes_ok,
        "exact_63_declared_tuples_matched": tuples_exact_match,
        "stable_seeds_verified_against_tuple_derivation": seeds_derived_ok,
        "same_s4_medium_and_b2_biomass_verified": same_medium_biomass_ok,
        "medium_hash_provenance_verified": med_hash_ok,
        "non_evaluable_handled_strictly_as_na": non_eval_clean,
        "zero_removed_recovery_rate_strictly_nan": zero_removed_ok,
        "evaluable_rows_arithmetically_consistent": eval_consistent,
        "per_solution_ledger_integrity_verified": ledger_ok,
    }
    all_passed = all(checks.values())

    return {
        "validation_passed": all_passed,
        "n_genomes_evaluated": len(actual_genomes),
        "n_total_mask_rows": len(masking_df),
        "n_evaluable_rows": int(len(eval_rows)),
        "n_non_evaluable_rows": int(len(non_eval)),
        "mean_recovery_rate": float(eval_rows["masked_rxn_recovery_rate_under_medium"].dropna().mean()) if not eval_rows.empty else math.nan,
        "mean_alternative_fraction": float(eval_rows["alternative_candidate_fraction_under_medium"].dropna().mean()) if not eval_rows.empty else math.nan,
        "scientific_qualification": (
            "Masking recovery evaluates candidate recovery under declared pool and S4 medium; "
            "candidate weights 100 for removed reactions reflect withheld sequence evidence. "
            "Non-recovered reactions reflect lack of essentiality under this nutrient condition."
        ),
        "checks": checks,
    }


def atomic_write_tsv(df: pd.DataFrame, dest_path: Path) -> str:
    """Atomically write DataFrame to TSV via .tmp file and return SHA-256."""
    tmp_path = dest_path.with_name(f"{dest_path.name}.tmp")
    df.to_csv(tmp_path, sep="\t", index=False)
    os.replace(tmp_path, dest_path)
    return compute_file_sha256(dest_path)


def atomic_write_json(data: Any, dest_path: Path) -> str:
    """Atomically write JSON data via .tmp file and return SHA-256."""
    tmp_path = dest_path.with_name(f"{dest_path.name}.tmp")
    with open(tmp_path, "w", encoding="utf-8") as fh:
        json.dump(data, fh, indent=2, sort_keys=True)
    os.replace(tmp_path, dest_path)
    return compute_file_sha256(dest_path)


def main() -> None:
    ap = argparse.ArgumentParser(description="Run Resumable Same-Medium Masking Controls on 7 HQ External Genomes (v4).")
    ap.add_argument("--manifest-tsv", type=Path, required=True)
    ap.add_argument("--pool-pkl", type=Path, required=True)
    ap.add_argument("--gapseq-out-dir", type=Path, required=True)
    ap.add_argument("--supp-models-dir", type=Path, required=True)
    ap.add_argument("--out-tsv", type=Path, default=None)
    ap.add_argument("--baselines-tsv", type=Path, default=None)
    ap.add_argument("--ledger-tsv", type=Path, default=None)
    ap.add_argument("--media-json", type=Path, default=None)
    ap.add_argument("--checkpoint-json", type=Path, default=None)
    ap.add_argument("--report-json", type=Path, default=None)
    ap.add_argument("--workers", type=int, default=32)
    ap.add_argument("--milp-time-limit-s", type=float, default=15.0)
    args = ap.parse_args()

    t_start = time.time()
    spawn_ctx = mp.get_context("spawn")

    # Authoritative fresh v4 namespace defaults
    results_dir = args.out_tsv.parent if args.out_tsv else Path("results")
    out_tsv = args.out_tsv or (results_dir / "exp1_masking_controls_7hq_v4.tsv")
    out_tsv.parent.mkdir(parents=True, exist_ok=True)

    baselines_tsv = args.baselines_tsv or (out_tsv.parent / "exp1_masking_controls_7hq_v4_baselines.tsv")
    ledger_tsv = args.ledger_tsv or (out_tsv.parent / "exp1_masking_controls_7hq_v4_per_solution_ledger.tsv")
    media_json = args.media_json or (out_tsv.parent / "exp1_masking_controls_v4_media_provenance.json")
    checkpoint_json = args.checkpoint_json or (out_tsv.parent / "exp1_masking_controls_v4_checkpoint.json")
    report_json = args.report_json or (out_tsv.parent / "exp1_masking_controls_7hq_v4_report.json")

    manifest_df = pd.read_csv(args.manifest_tsv, sep="\t")
    hq_df = manifest_df[manifest_df["genome"].isin(EXACT_7_HQ_EXTERNAL_IDS)].drop_duplicates("genome").copy()
    if len(hq_df) != 7:
        raise ValueError(f"Expected exactly 7 HQ external genomes, found {len(hq_df)} in manifest")

    print(f"[{time.strftime('%X')}] [v4 Resumable / Spawn] Loaded {len(hq_df)} HQ external genomes for masking controls.")
    print(f"[{time.strftime('%X')}] Namespace outputs:")
    print(f"       out_tsv:         {out_tsv}")
    print(f"       baselines_tsv:   {baselines_tsv}")
    print(f"       ledger_tsv:      {ledger_tsv}")
    print(f"       media_json:      {media_json}")
    print(f"       checkpoint_json: {checkpoint_json}")
    print(f"       report_json:     {report_json}")

    pool_sha256 = compute_file_sha256(args.pool_pkl)
    print(f"[{time.strftime('%X')}] Candidate pool {args.pool_pkl.name} SHA-256: {pool_sha256[:16]}...")

    with open(args.pool_pkl, "rb") as fh:
        pool_data = pickle.load(fh)
    universal_rxns = pool_data["reactions"]

    # Compute engine and solver configuration hashes for strict resume validation
    engine_file = CODE_DIR / "exp1_sbml_highs_engine_v4.py"
    if not engine_file.exists():
        engine_file = CODE_DIR / "exp1_sbml_highs_engine.py"
    engine_sha256 = compute_file_sha256(engine_file)
    solver_config = {
        "milp_time_limit_s": args.milp_time_limit_s,
        "mip_rel_gap": 0.05,
        "mip_feasibility_tolerance": 1e-9,
        "primal_feasibility_tolerance": 1e-9,
        "threads": 1,
        "biomass_recipe_id": BIO_B2_OMIT_7_KEEP_F420,
        "scenario_id": SCEN_S4_OXIC_DESERT,
    }
    solver_config_str = json.dumps(solver_config, sort_keys=True)
    solver_config_hash = hashlib.sha256(solver_config_str.encode("utf-8")).hexdigest()

    # Step 0: Checkpoint validation for resume
    committed_checkpoint: dict[str, Any] = {}
    valid_resume = False
    if checkpoint_json.exists():
        try:
            with open(checkpoint_json, "r", encoding="utf-8") as fh:
                chk = json.load(fh)
            # Verify engine and solver configuration match
            if (
                chk.get("engine_sha256") == engine_sha256
                and chk.get("solver_config_hash") == solver_config_hash
                and baselines_tsv.exists()
                and compute_file_sha256(baselines_tsv) == chk.get("baselines_tsv_sha256")
                and ledger_tsv.exists()
                and compute_file_sha256(ledger_tsv) == chk.get("ledger_tsv_sha256")
            ):
                committed_checkpoint = chk
                valid_resume = True
                print(f"[{time.strftime('%X')}] Valid commit marker verified. Resuming from prior v4 checkpoint.")
            else:
                print(f"[{time.strftime('%X')}] Notice: commit marker invalid or code/config changed; recomputing fresh.")
        except Exception as e:
            print(f"[{time.strftime('%X')}] Notice: could not load checkpoint marker: {e}")

    # Build genome payloads and media provenance
    media_provenance: dict[str, Any] = {}
    genome_payloads: dict[str, dict[str, Any]] = {}
    for _, row in hq_df.iterrows():
        g = str(row["genome"])
        tar_path = args.gapseq_out_dir / f"{g}.tar.gz"
        archive_sha = compute_file_sha256(tar_path)
        payload = load_genome_archive_payload(g, args.gapseq_out_dir, args.supp_models_dir)
        draft_mod = payload["draft_model"]
        cand_weights = build_genome_candidate_weights(universal_rxns, payload["rxn_evi"], payload["trans_evi"])

        d_med = payload.get("desert_med")
        if d_med:
            s4_cpds = {**d_med, "cpd00007": ("O2", 10.0)}
        else:
            s4_cpds = dict(NUTRIENT_SCENARIOS[SCEN_S4_OXIC_DESERT]["compounds"])

        ser_json, med_hash = serialize_and_hash_medium(s4_cpds)
        media_provenance[g] = {
            "genome": g,
            "archive_sha256": archive_sha,
            "medium_hash": med_hash,
            "compounds_json": ser_json,
        }

        inp_hash = compute_exact_input_provenance_hash(
            archive_sha, draft_mod.model_id, med_hash, pool_sha256, BIO_B2_OMIT_7_KEEP_F420
        )
        base_case_hash = compute_case_config_hash(inp_hash, solver_config_hash, engine_sha256, f"{g}_base")

        genome_payloads[g] = {
            "genome": g,
            "cohort": str(row["cohort"]),
            "canonical_order": str(row.get("canonical_order", "")),
            "canonical_family": str(row.get("canonical_family", "")),
            "completeness": row["completeness"],
            "contamination": row["contamination"],
            "draft_model": draft_mod,
            "s4_compounds": s4_cpds,
            "medium_hash": med_hash,
            "archive_sha256": archive_sha,
            "candidate_weights": cand_weights,
            "input_hash": inp_hash,
            "case_config_hash": base_case_hash,
        }

    # Load pre-existing valid data if resume is valid
    existing_baselines: dict[str, dict[str, Any]] = {}
    existing_ledger_rows: list[dict[str, Any]] = []
    existing_mask_rows: list[dict[str, Any]] = []
    existing_mask_tuples: set[tuple[str, float, int]] = set()

    if valid_resume:
        chk_bases = committed_checkpoint.get("baseline_records", {})
        if baselines_tsv.exists():
            b_df = pd.read_csv(baselines_tsv, sep="\t")
            for _, r in b_df.iterrows():
                g = str(r["genome"])
                # Only accept if case_config_hash matches
                if chk_bases.get(g) == genome_payloads[g]["case_config_hash"]:
                    rxns_str = str(r.get("gapfill_reaction_ids", ""))
                    gap_rxns = [x.strip() for x in rxns_str.split(",") if x.strip()] if pd.notna(rxns_str) and rxns_str != "nan" else []
                    existing_baselines[g] = {
                        "summary_row": r.to_dict(),
                        "status": str(r["solver_status"]),
                        "n_valid_solutions": int(r["n_valid_solutions"]),
                        "gapfill_rxns": gap_rxns,
                    }
        if ledger_tsv.exists():
            l_df = pd.read_csv(ledger_tsv, sep="\t")
            existing_ledger_rows = l_df.to_dict("records")
        if out_tsv.exists() and compute_file_sha256(out_tsv) == committed_checkpoint.get("summary_tsv_sha256"):
            m_df = pd.read_csv(out_tsv, sep="\t")
            chk_masks = committed_checkpoint.get("mask_records", {})
            for _, r in m_df.iterrows():
                tup = (str(r["genome"]), round(float(r["mask_fraction"]), 2), int(r["replicate"]))
                cid = f"{tup[0]}_{tup[1]:.2f}_{tup[2]}"
                if chk_masks.get(cid) == r.get("case_config_hash"):
                    existing_mask_tuples.add(tup)
                    existing_mask_rows.append(r.to_dict())

    # Step 1: Execute remaining baselines IN PARALLEL with spawn
    baseline_tasks: list[dict[str, Any]] = []
    for g, payload in genome_payloads.items():
        if g not in existing_baselines:
            baseline_tasks.append({
                **payload,
                "candidate_reactions": universal_rxns,
                "milp_time_limit_s": args.milp_time_limit_s,
            })

    all_baseline_rows: list[dict[str, Any]] = [v["summary_row"] for v in existing_baselines.values()]
    all_ledger_rows: list[dict[str, Any]] = list(existing_ledger_rows)
    committed_base_records: dict[str, str] = {g: existing_baselines[g]["summary_row"].get("case_config_hash", "") for g in existing_baselines}
    committed_mask_records: dict[str, str] = {f"{r['genome']}_{float(r['mask_fraction']):.2f}_{int(r['replicate'])}": r.get("case_config_hash", "") for r in existing_mask_rows}

    if baseline_tasks:
        n_workers_base = min(len(baseline_tasks), args.workers, 7)
        print(f"\n[{time.strftime('%X')}] Step 1: Solving {len(baseline_tasks)} baselines in parallel across {n_workers_base} spawn workers...")
        with ProcessPoolExecutor(max_workers=n_workers_base, mp_context=spawn_ctx) as ex:
            futs = {ex.submit(evaluate_single_baseline_case, t): t["genome"] for t in baseline_tasks}
            for fut in as_completed(futs):
                g = futs[fut]
                res = fut.result()
                existing_baselines[g] = {
                    "summary_row": res["summary_row"],
                    "status": res["status"],
                    "n_valid_solutions": res["n_valid_solutions"],
                    "gapfill_rxns": res["gapfill_rxns"],
                }
                all_baseline_rows.append(res["summary_row"])
                all_ledger_rows.extend(res["ledger_rows"])
                committed_base_records[g] = res["summary_row"]["case_config_hash"]
                print(f"  [BASELINE DONE] {g}: status={res['status']}, n_valid_solutions={res['n_valid_solutions']}, n_gapfill_rxns={len(res['gapfill_rxns'])} in {res['runtime_seconds']}s", flush=True)

        # Atomic commit of all 7 baselines
        base_sha = atomic_write_tsv(pd.DataFrame(all_baseline_rows).sort_values("genome"), baselines_tsv)
        med_sha = atomic_write_json(media_provenance, media_json)
        led_sha = atomic_write_tsv(pd.DataFrame(all_ledger_rows).sort_values(["condition_type", "genome", "mask_fraction", "replicate", "solution_index"]), ledger_tsv)
        chk_init = {
            "committed_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
            "engine_sha256": engine_sha256,
            "solver_config_hash": solver_config_hash,
            "n_baselines_committed": len(all_baseline_rows),
            "n_mask_cases_committed": len(existing_mask_rows),
            "n_ledger_rows_committed": len(all_ledger_rows),
            "baselines_tsv_sha256": base_sha,
            "summary_tsv_sha256": compute_file_sha256(out_tsv),
            "ledger_tsv_sha256": led_sha,
            "baseline_records": committed_base_records,
            "mask_records": committed_mask_records,
        }
        atomic_write_json(chk_init, checkpoint_json)
        print(f"[{time.strftime('%X')}] All 7 baselines complete and committed atomically to {baselines_tsv.name}.")
    else:
        print(f"[{time.strftime('%X')}] All 7 baselines verified and loaded from checkpoint.")

    # Step 2: Prepare 63 masking tasks
    mask_tasks: list[dict[str, Any]] = []
    for g, payload in genome_payloads.items():
        base_info = existing_baselines[g]
        for frac in (0.10, 0.20, 0.30):
            for rep in (0, 1, 2):
                tup = (g, round(frac, 2), rep)
                if tup in existing_mask_tuples:
                    continue
                cid = f"{g}_{frac:.2f}_{rep}"
                mask_case_hash = compute_case_config_hash(payload["input_hash"], solver_config_hash, engine_sha256, cid)
                mask_tasks.append({
                    "genome": g,
                    "cohort": payload["cohort"],
                    "canonical_order": payload["canonical_order"],
                    "canonical_family": payload["canonical_family"],
                    "completeness": payload["completeness"],
                    "contamination": payload["contamination"],
                    "mask_fraction": frac,
                    "replicate": rep,
                    "draft_model": payload["draft_model"],
                    "s4_compounds": payload["s4_compounds"],
                    "medium_hash": payload["medium_hash"],
                    "candidate_reactions": universal_rxns,
                    "candidate_weights": payload["candidate_weights"],
                    "milp_time_limit_s": args.milp_time_limit_s,
                    "input_hash": payload["input_hash"],
                    "case_config_hash": mask_case_hash,
                    "baseline_status": base_info["status"],
                    "baseline_n_valid_solutions": base_info["n_valid_solutions"],
                    "baseline_gapfill_rxns": base_info["gapfill_rxns"],
                })

    masking_summary_rows: list[dict[str, Any]] = list(existing_mask_rows)

    if mask_tasks:
        n_workers_mask = min(args.workers, len(mask_tasks))
        print(f"\n[{time.strftime('%X')}] Step 2: Executing {len(mask_tasks)} masking control runs across {n_workers_mask} spawn workers...")
        with ProcessPoolExecutor(max_workers=n_workers_mask, mp_context=spawn_ctx) as ex:
            futs = {ex.submit(evaluate_single_mask_case, t): (t["genome"], t["mask_fraction"], t["replicate"], t["case_config_hash"]) for t in mask_tasks}
            count = 0
            for fut in as_completed(futs):
                g, frac, rep, c_hash = futs[fut]
                res = fut.result()
                masking_summary_rows.append(res["summary_row"])
                all_ledger_rows.extend(res["ledger_rows"])
                committed_mask_records[f"{g}_{frac:.2f}_{rep}"] = c_hash
                srow = res["summary_row"]
                count += 1
                rec_val = srow['masked_rxn_recovery_rate_under_medium']
                rec_disp = f"{rec_val:.3f}" if (rec_val is not None and not (isinstance(rec_val, float) and math.isnan(rec_val))) else "NaN"
                alt_val = srow['alternative_candidate_fraction_under_medium']
                alt_disp = f"{alt_val:.3f}" if (alt_val is not None and not (isinstance(alt_val, float) and math.isnan(alt_val))) else "NaN"
                print(f"  [{count}/{len(mask_tasks)}] {g} frac={frac:.2f} rep={rep}: evaluable={srow['is_evaluable']}, rec={rec_disp}, alt={alt_disp} ({srow['runtime_seconds']}s)", flush=True)

                # Periodic atomic checkpointing every 5 runs or on final run
                if count % 5 == 0 or count == len(mask_tasks):
                    out_sha = atomic_write_tsv(pd.DataFrame(masking_summary_rows).sort_values(["genome", "mask_fraction", "replicate"]), out_tsv)
                    led_sha = atomic_write_tsv(pd.DataFrame(all_ledger_rows).sort_values(["condition_type", "genome", "mask_fraction", "replicate", "solution_index"]), ledger_tsv)
                    chk_upd = {
                        "committed_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
                        "engine_sha256": engine_sha256,
                        "solver_config_hash": solver_config_hash,
                        "n_baselines_committed": len(all_baseline_rows),
                        "n_mask_cases_committed": len(masking_summary_rows),
                        "n_ledger_rows_committed": len(all_ledger_rows),
                        "baselines_tsv_sha256": compute_file_sha256(baselines_tsv),
                        "summary_tsv_sha256": out_sha,
                        "ledger_tsv_sha256": led_sha,
                        "baseline_records": committed_base_records,
                        "mask_records": committed_mask_records,
                    }
                    atomic_write_json(chk_upd, checkpoint_json)
    else:
        print(f"[{time.strftime('%X')}] All 63 masking runs already complete and verified!")

    out_df = pd.DataFrame(masking_summary_rows).sort_values(["genome", "mask_fraction", "replicate"]).reset_index(drop=True)
    atomic_write_tsv(out_df, out_tsv)

    ledger_df = pd.DataFrame(all_ledger_rows).sort_values(["condition_type", "genome", "mask_fraction", "replicate", "solution_index"]).reset_index(drop=True)
    atomic_write_tsv(ledger_df, ledger_tsv)

    report = validate_masking_controls_results(
        masking_df=out_df,
        expected_genomes=set(EXACT_7_HQ_EXTERNAL_IDS),
        ledger_df=ledger_df,
        media_provenance=media_provenance,
    )
    report["elapsed_seconds"] = round(time.time() - t_start, 2)
    report["candidate_pool_sha256"] = pool_sha256
    report["engine_sha256"] = engine_sha256
    report["solver_config_hash"] = solver_config_hash

    atomic_write_json(report, report_json)

    print(f"\n[{time.strftime('%X')}] Masking Controls Report:")
    print(json.dumps(report, indent=2))

    if not report["validation_passed"]:
        print(f"[ERROR] Masking controls validation FAILED checks: {report['checks']}", file=sys.stderr)
        sys.exit(2)

    print(f"\n[SUCCESS] Completed {len(out_df)} masking control runs across 7 HQ genomes in {report['elapsed_seconds']}s.")


if __name__ == "__main__":
    main()
