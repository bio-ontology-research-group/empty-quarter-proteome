#!/usr/bin/env python3
"""Core Engine and Post-Hoc Correction Logic for Experiment 1 Counterfactual Closure
(`EXP1-COUNTERFACTUAL-CLOSURE-IMPLEMENTATION.txt`, `EXP1-COUNTERFACTUAL-ROOT-SERIALIZATION-FINDING.txt`,
`EXP1-COUNTERFACTUAL-PRELAUNCH-REVIEW.txt`, and `EXP1-FVA-THREAD-ORDER-DIAGNOSIS.txt`).

Owned exclusively inside `reproducibility/mechanism_20261109/exp1_metabolic/counterfactual_review/`.
Does not modify any existing `exp1_metabolic/code/` or collector files.
"""
from __future__ import annotations

import hashlib
import json
import math
import multiprocessing as mp
import re
import sys
import warnings
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.optimize import OptimizeWarning, linprog

REVALIDATED_FEASIBLE_STATUSES = {
    "optimal",
    "optimal_after_cycle_qc_filtering",
    "optimal_draft_no_gapfill_needed",
    "heuristic_incumbent_bounded",
    "post_pruning_cost_changed_not_proven_optimal",
}

WKF_SCHEDULED_SCENARIOS = (
    "S3_anoxic_desert_oa_aa",
    "S4_oxic_desert_oa_aa",
    "S5_anoxic_defined_aa_precursors",
    "S6_oxic_defined_aa_precursors",
)

BIO_B0_ID = "B0_methanosarcina_template_baseline"
BIO_B2_ID = "B2_perturbation_omit_7_non_f420_template_cofactors_keep_f420"

WKF_SCHEDULED_BIOMASS = (
    BIO_B0_ID,
    BIO_B2_ID,
)

EXPECTED_CODE_V2_SHA256 = {
    "exp1_sbml_highs_engine.py": "c4ef47caa6460d46db7915f42128bc0c9903bd42ae60fcdf75f8e0b0c3713547",
    "stage2_pilot24_ensemble.py": "bdf6e6eb1d75a1ae774c18edc40523d189602b550d40bb32aaf7dc0f718b56ec",
    "stage3_full484_worker.py": "273bd9355bec84a45806f1aef890074079d54fe22936df3a8dd4650be4d3fee7",
}


def compute_file_sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        while chunk := fh.read(1024 * 1024):
            h.update(chunk)
    return h.hexdigest()


def verify_code_v2_integrity(code_v2_dir: Path) -> dict[str, str]:
    """Verify that the frozen `code_v2` directory matches the exact executed SHA-256 digests."""
    actual_hashes: dict[str, str] = {}
    for fname, exp_sha in EXPECTED_CODE_V2_SHA256.items():
        fpath = code_v2_dir / fname
        if not fpath.is_file():
            raise FileNotFoundError(f"Missing frozen code_v2 file: {fpath}")
        act_sha = compute_file_sha256(fpath)
        if act_sha != exp_sha:
            raise RuntimeError(
                f"SHA-256 mismatch on frozen {fpath}: expected {exp_sha}, got {act_sha}"
            )
        actual_hashes[fname] = act_sha
    return actual_hashes


def initialize_highs_single_thread_scheduler() -> dict[str, Any]:
    """Initialize the process-global HiGHS thread scheduler to `threads=1` BEFORE any default-option
    HiGHS call can occur in the current process (`EXP1-FVA-THREAD-ORDER-DIAGNOSIS.txt` Item 1).

    In SciPy's HiGHS wrapper, if the first `linprog(..., method='highs')` call in a process uses default
    options (hardware threads) and a subsequent call passes `options={'threads': 1}`, HiGHS fails with
    `status=4 ('HiGHS Status 0: Not Set')`. Conversely, when the first HiGHS call in the process sets
    `options={'threads': 1}`, all subsequent default-option (`solve_fba_highs`, `detect_closed_uptake_egc`)
    and `threads=1` FVA calls succeed with `status=0`.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=OptimizeWarning)
        res = linprog(
            c=np.array([1.0], dtype=np.float64),
            bounds=[(1.0, 2.0)],
            method="highs",
            options={
                "threads": 1,
                "presolve": True,
                "primal_feasibility_tolerance": 1e-9,
            },
        )
    raw_st = int(getattr(res, "status", -1))
    msg = str(getattr(res, "message", ""))
    x0 = float(res.x[0]) if (res.success and res.x is not None and len(res.x) == 1) else math.nan
    if not res.success or raw_st != 0 or not math.isfinite(x0) or abs(x0 - 1.0) > 1e-12:
        raise RuntimeError(
            f"Failed to initialize HiGHS single-thread scheduler at worker startup: "
            f"status={raw_st}, success={res.success}, x0={x0}, message={msg!r}"
        )
    return {
        "initialized_threads": 1,
        "status": raw_st,
        "success": bool(res.success),
        "x0": x0,
        "message": msg,
    }


def _child_probe_default_first_then_threads1(q: Any) -> None:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=OptimizeWarning)
        r1 = linprog(c=[1.0], bounds=[(1.0, 2.0)], method="highs")
        r2 = linprog(c=[1.0], bounds=[(1.0, 2.0)], method="highs", options={"threads": 1})
        r3 = linprog(c=[1.0], bounds=[(1.0, 2.0)], method="highs")
    q.put(
        {
            "step1_default": {"status": int(r1.status), "success": bool(r1.success), "message": str(r1.message)},
            "step2_threads1": {"status": int(r2.status), "success": bool(r2.success), "message": str(r2.message)},
            "step3_default": {"status": int(r3.status), "success": bool(r3.success), "message": str(r3.message)},
        }
    )


def _child_probe_threads1_first_then_default_then_threads1(q: Any) -> None:
    init_res = initialize_highs_single_thread_scheduler()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=OptimizeWarning)
        r2 = linprog(c=[1.0], bounds=[(1.0, 2.0)], method="highs")
        r3 = linprog(c=[1.0], bounds=[(1.0, 2.0)], method="highs", options={"threads": 1})
    q.put(
        {
            "step1_threads1_init": init_res,
            "step2_default": {"status": int(r2.status), "success": bool(r2.success), "message": str(r2.message)},
            "step3_threads1": {"status": int(r3.status), "success": bool(r3.success), "message": str(r3.message)},
        }
    )


def reproduce_highs_thread_order_defect_in_spawned_subprocesses() -> dict[str, Any]:
    """Reproduce the exact HiGHS global thread scheduler call-order defect and verify the fix in fresh
    spawned subprocesses (`EXP1-FVA-THREAD-ORDER-DIAGNOSIS.txt` Item 1).
    """
    ctx = mp.get_context("spawn")

    q1 = ctx.Queue()
    p1 = ctx.Process(target=_child_probe_default_first_then_threads1, args=(q1,))
    p1.start()
    uninit_res = q1.get(timeout=30)
    p1.join(timeout=10)

    q2 = ctx.Queue()
    p2 = ctx.Process(target=_child_probe_threads1_first_then_default_then_threads1, args=(q2,))
    p2.start()
    fixed_res = q2.get(timeout=30)
    p2.join(timeout=10)

    defect_reproduced = bool(
        uninit_res["step1_default"]["status"] == 0
        and uninit_res["step2_threads1"]["status"] == 4
        and uninit_res["step3_default"]["status"] == 0
    )
    fix_verified = bool(
        fixed_res["step1_threads1_init"]["status"] == 0
        and fixed_res["step2_default"]["status"] == 0
        and fixed_res["step3_threads1"]["status"] == 0
    )
    if not fix_verified:
        raise RuntimeError(f"HiGHS thread-1 initialization fix failed verification: {fixed_res}")

    return {
        "default_first_then_threads1_probe": uninit_res,
        "threads1_first_then_default_then_threads1_probe": fixed_res,
        "defect_reproduced_status4_when_default_first": defect_reproduced,
        "fix_verified_all_optimal_when_threads1_first": fix_verified,
    }


def classify_withdrawal_counterfactual_v2(
    baseline_growth: float,
    withdrawn_growth_primary: float,
    withdrawn_growth_union: float,
    addback_rescue_growth: float,
    min_growth: float = 0.01,
) -> dict[str, Any]:
    """Classify a `nutrient_withdrawal_addback` row without threshold gaps or rescue conflation
    (`EXP1-COUNTERFACTUAL-CLOSURE-IMPLEMENTATION.txt` Item 2 & `EXP1-COUNTERFACTUAL-PRELAUNCH-REVIEW.txt` Item 2).
    """
    if (
        not math.isfinite(baseline_growth)
        or baseline_growth < 1e-6
        or not math.isfinite(withdrawn_growth_primary)
        or not math.isfinite(withdrawn_growth_union)
    ):
        return {
            "primary_growth_lost_lt_1e-4": None,
            "addback_rescues_growth_v2": None,
            "addback_rescue_disposition": "unknown_solver_or_invalid_baseline",
            "addback_to_baseline_growth_ratio": math.nan,
            "corrected_classification": "unknown_solver_or_invalid_baseline",
            "union_growth_regime": "unknown_solver_or_invalid_baseline",
        }

    addback_ratio = (
        float(addback_rescue_growth / baseline_growth)
        if math.isfinite(addback_rescue_growth)
        else math.nan
    )

    if withdrawn_growth_union < 1e-4:
        u_regime = "lt_1e-4_abolished_in_union"
    elif withdrawn_growth_union < 0.5 * min_growth:
        u_regime = "1e-4_to_0.005_partial_subthreshold_in_union"
    elif withdrawn_growth_union < min_growth:
        u_regime = "0.005_to_0.01_below_model_min_growth_in_union"
    elif withdrawn_growth_union < 0.95 * baseline_growth - 1e-6:
        u_regime = "0.01_to_95pct_g0_above_model_min_growth_reduced_in_union"
    elif withdrawn_growth_union <= baseline_growth + 1e-6:
        u_regime = "95pct_g0_to_g0_near_baseline_in_union"
    else:
        u_regime = "gt_g0_union_exceeds_primary_baseline"

    if withdrawn_growth_primary < 1e-4:
        if not math.isfinite(addback_rescue_growth):
            rescue_flag: bool | None = None
            rescue_disp = "unknown_addback_solver_failure"
        elif addback_rescue_growth >= min_growth * 0.99:
            rescue_flag = True
            rescue_disp = "restored_primary_growth_above_declared_threshold_0.99_min_growth"
        else:
            rescue_flag = False
            rescue_disp = "failed_addback_threshold_restoration"

        if withdrawn_growth_union < 1e-4:
            corr_cls = "required_in_primary_and_sampled_ensemble_union_lt_1e-4"
        elif withdrawn_growth_union < 0.5 * min_growth:
            corr_cls = "primary_growth_lost_lt_1e-4__partial_subthreshold_union_growth_1e-4_to_0.005"
        elif withdrawn_growth_union < min_growth:
            corr_cls = "primary_growth_lost_lt_1e-4__below_model_min_growth_union_0.005_to_0.01"
        else:
            corr_cls = "primary_growth_lost_lt_1e-4__above_model_min_growth_union_bypass_ge_0.01_relaxed_union"

        return {
            "primary_growth_lost_lt_1e-4": True,
            "addback_rescues_growth_v2": rescue_flag,
            "addback_rescue_disposition": rescue_disp,
            "addback_to_baseline_growth_ratio": addback_ratio,
            "corrected_classification": corr_cls,
            "union_growth_regime": u_regime,
        }

    if not math.isfinite(addback_rescue_growth):
        rescue_disp = "unknown_addback_solver_failure"
    else:
        rescue_disp = "not_applicable_primary_growth_retained_ge_1e-4"

    if withdrawn_growth_primary < min_growth:
        corr_cls = "primary_growth_retained_ge_1e-4_below_model_min_growth_lt_0.01"
    else:
        corr_cls = "primary_growth_retained_above_model_min_growth_ge_0.01_non_essential_or_stimulatory"

    return {
        "primary_growth_lost_lt_1e-4": False,
        "addback_rescues_growth_v2": None,
        "addback_rescue_disposition": rescue_disp,
        "addback_to_baseline_growth_ratio": addback_ratio,
        "corrected_classification": corr_cls,
        "union_growth_regime": u_regime,
    }


def classify_reaction_knockout_counterfactual_v2(
    baseline_growth: float,
    ko_growth_union: float,
    min_growth: float = 0.01,
) -> dict[str, Any]:
    """Classify a `reaction_knockout_and_fva` row while marking unperformed primary-KO and
    addback fields as `NaN`/`None` with explicit `"not_performed_in_executed_code_v2"` status.
    """
    if not math.isfinite(baseline_growth) or baseline_growth < 1e-6 or not math.isfinite(ko_growth_union):
        u_regime = "unknown_solver_or_invalid_baseline"
        corr_cls = "unknown_solver_or_invalid_baseline"
    elif ko_growth_union < 1e-4:
        u_regime = "lt_1e-4_essential_in_union"
        corr_cls = "essential_in_sampled_ensemble_union_lt_1e-4"
    elif ko_growth_union < min_growth:
        u_regime = "1e-4_to_0.01_below_model_min_growth_in_union"
        corr_cls = "subthreshold_in_sampled_ensemble_union_1e-4_to_0.01_below_min_growth"
    elif ko_growth_union < 0.95 * baseline_growth - 1e-6:
        u_regime = "0.01_to_95pct_g0_above_model_min_growth_reduced_in_union"
        corr_cls = "bypassable_in_sampled_ensemble_union_reduced_0.01_to_95pct_g0_relaxed_union"
    elif ko_growth_union <= baseline_growth + 1e-6:
        u_regime = "95pct_g0_to_g0_near_baseline_in_union"
        corr_cls = "bypassable_in_sampled_ensemble_union_near_baseline_95pct_g0_to_g0_relaxed_union"
    else:
        u_regime = "gt_g0_union_exceeds_primary_baseline"
        corr_cls = "bypassable_in_sampled_ensemble_union_exceeds_primary_baseline_gt_g0_relaxed_union"

    return {
        "perturbed_growth_primary_v2": math.nan,
        "primary_perturbation_measurement_status": "not_performed_in_executed_code_v2",
        "addback_rescue_growth_v2": math.nan,
        "addback_measurement_status": "not_performed_in_executed_code_v2",
        "addback_rescues_growth_v2": None,
        "addback_rescue_disposition": "not_performed_in_executed_code_v2",
        "addback_to_baseline_growth_ratio": math.nan,
        "corrected_classification": corr_cls,
        "union_growth_regime": u_regime,
    }


def classify_union_network_qc(
    uni_fba_ok: bool,
    uni_fba_status: str,
    uni_egc: dict[str, Any],
    uni_cyc_ok: bool,
) -> str:
    """Classify the mathematical union network's numerical and closed-uptake EGC status
    (`EXP1-COUNTERFACTUAL-PRELAUNCH-REVIEW.txt` Item 3).
    Never claims the union is physically established; distinguishes unknown EGC from positive cycle failure.
    """
    if not uni_fba_ok:
        return f"mathematical_union_failed_fba_or_numerical_checks__{uni_fba_status}"
    spur = uni_egc.get("has_spurious_atp_cycle", None)
    egc_st = str(uni_egc.get("egc_check_status", "unknown"))
    if spur is False and uni_cyc_ok:
        return "mathematical_union_passed_declared_numerical_and_closed_uptake_egc_checks"
    if spur is True:
        return "mathematical_union_failed_closed_uptake_atp_cycle_check"
    return f"mathematical_union_unknown_closed_uptake_egc_check__{egc_st}"


def verify_primary_reconstruction_gate(
    sol0_status: str,
    prim_cyc_ok: bool,
    g0: float,
    ledger_g0: float,
    wkf_baseline_g0: float,
    max_residual: float,
    max_bound_viol: float,
    primary_identity_matches_wkf: bool = True,
    g0_tol: float = 1e-6,
    residual_tol: float = 1e-5,
) -> tuple[bool, str]:
    """Fail-closed gate for primary (`solution_index == 0`) model reconstruction
    (`EXP1-COUNTERFACTUAL-PRELAUNCH-REVIEW.txt` Item 1).
    """
    if not primary_identity_matches_wkf:
        return False, "unverified_primary_reaction_set_mismatch_vs_wkf"
    if sol0_status != "optimal" or not math.isfinite(g0):
        return False, f"unverified_primary_fba_solver_{sol0_status}"
    if not prim_cyc_ok:
        return False, "unverified_primary_failed_closed_uptake_egc_check"
    if g0 < 0.00999:
        return False, "unverified_primary_growth_below_min_growth"
    if not math.isfinite(max_residual) or max_residual > residual_tol:
        return False, "unverified_primary_residual_exceeds_1e-5"
    if not math.isfinite(max_bound_viol) or max_bound_viol > residual_tol:
        return False, "unverified_primary_bound_violation_exceeds_1e-5"
    if not math.isfinite(ledger_g0) or abs(g0 - ledger_g0) > g0_tol:
        return False, "unverified_primary_g0_mismatch_vs_solution_ledger"
    if not math.isfinite(wkf_baseline_g0) or abs(g0 - wkf_baseline_g0) > g0_tol:
        return False, "unverified_primary_g0_mismatch_vs_wkf_baseline"
    return True, "verified_primary_reconstruction"


def classify_fva_interval_v2(
    fva_min: float,
    fva_max: float,
    min_status: str = "optimal",
    max_status: str = "optimal",
    ordering_tol: float = 1e-9,
) -> dict[str, Any]:
    """Classify an FVA `[fva_min, fva_max]` interval with strict finite, status, and ordering checks.
    Never coerces failed/unknown LPs to 0.0 or False (`EXP1-COUNTERFACTUAL-CLOSURE-IMPLEMENTATION.txt` Item 5).
    """
    if (
        min_status != "optimal"
        or max_status != "optimal"
        or not math.isfinite(fva_min)
        or not math.isfinite(fva_max)
    ):
        return {
            "fva_both_optimal_and_finite": False,
            "fva_interval_ordered": None,
            "fva_span": math.nan,
            "fva_strictly_nonzero_1e-6": None,
            "fva_strictly_nonzero_1e-12": None,
            "fva_interval_class": f"unknown_fva_solver__min_{min_status}__max_{max_status}",
        }

    if fva_min > fva_max + ordering_tol:
        return {
            "fva_both_optimal_and_finite": True,
            "fva_interval_ordered": False,
            "fva_span": math.nan,
            "fva_strictly_nonzero_1e-6": None,
            "fva_strictly_nonzero_1e-12": None,
            "fva_interval_class": "invalid_inverted_fva_interval",
        }

    span = float(max(0.0, fva_max - fva_min))
    strictly_1e6 = bool(fva_min > 1e-6 or fva_max < -1e-6)
    strictly_1e12 = bool(fva_min > 1e-12 or fva_max < -1e-12)

    if fva_min > 1e-6:
        icls = "strictly_positive_gt_1e-6"
    elif fva_min > 1e-12:
        icls = "strictly_positive_trace_1e-12_to_1e-6"
    elif fva_max < -1e-6:
        icls = "strictly_negative_lt_-1e-6"
    elif fva_max < -1e-12:
        icls = "strictly_negative_trace_-1e-6_to_-1e-12"
    elif abs(fva_min) <= 1e-12 and abs(fva_max) <= 1e-12:
        icls = "locked_zero_within_1e-12"
    else:
        icls = "spans_zero_contains_zero_flux"

    return {
        "fva_both_optimal_and_finite": True,
        "fva_interval_ordered": True,
        "fva_span": span,
        "fva_strictly_nonzero_1e-6": strictly_1e6,
        "fva_strictly_nonzero_1e-12": strictly_1e12,
        "fva_interval_class": icls,
    }


def make_unverified_primary_fva_result(reason: str) -> dict[str, Any]:
    """Return an explicit fail-closed unknown FVA result when the primary baseline reconstruction
    is failed or unverified (`EXP1-COUNTERFACTUAL-PRELAUNCH-REVIEW.txt` Item 1).
    """
    return {
        "status": f"unknown_{reason}",
        "raw_status": -1,
        "message": f"FVA not executed because primary baseline failed verification: {reason}",
        "target_flux": math.nan,
        "achieved_biomass_growth": math.nan,
        "max_residual": math.nan,
        "max_bound_viol": math.nan,
        "n_nonzero_fluxes": 0,
        "flux_witness_sha256": None,
        "sparse_fluxes": {},
    }


def solve_verified_linear_objective_highs(
    sp_ids: list[str],
    rxn_ids: list[str],
    S: sparse.csc_matrix,
    lb: np.ndarray,
    ub: np.ndarray,
    obj_j: int,
    target_j: int,
    sense: str,
    min_Required_growth: float,
    time_limit_s: float = 10.0,
    feasibility_tol: float = 1e-9,
    residual_tol: float = 1e-5,
) -> dict[str, Any]:
    """Solve a single minimization (`sense == 'min'`) or maximization (`sense == 'max'`) LP on `rxn_ids[target_j]`
    using HiGHS (`scipy.optimize.linprog`) with `threads=1` and the exact independent post-solve residual,
    bound-violation, finite-vector, and biomass-floor checks as `solve_fba_highs`.
    """
    if sense not in {"min", "max"}:
        raise ValueError(f"sense must be 'min' or 'max', got {sense!r}")

    if (
        not math.isfinite(min_Required_growth)
        or min_Required_growth <= 0.0
        or not math.isfinite(float(lb[obj_j]))
        or float(lb[obj_j]) < min_Required_growth - 1e-12
    ):
        return make_unverified_primary_fva_result("unverified_primary_growth_floor")

    c = np.zeros(len(rxn_ids), dtype=np.float64)
    c[target_j] = 1.0 if sense == "min" else -1.0
    b_eq = np.zeros(len(sp_ids), dtype=np.float64)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=OptimizeWarning)
        res = linprog(
            c=c,
            A_eq=S,
            b_eq=b_eq,
            bounds=list(zip(lb, ub)),
            method="highs",
            options={
                "time_limit": float(time_limit_s),
                "presolve": True,
                "primal_feasibility_tolerance": float(feasibility_tol),
                "threads": 1,
            },
        )
    raw_st = int(getattr(res, "status", -1))
    if not res.success or res.x is None:
        status_map = {
            1: "timeout",
            2: "infeasible",
            3: "unbounded",
            4: "numerical_error",
        }
        return {
            "status": status_map.get(raw_st, "solver_failed"),
            "raw_status": raw_st,
            "message": str(getattr(res, "message", "")),
            "target_flux": math.nan,
            "achieved_biomass_growth": math.nan,
            "max_residual": math.nan,
            "max_bound_viol": math.nan,
            "n_nonzero_fluxes": 0,
            "flux_witness_sha256": None,
            "sparse_fluxes": {},
        }

    x_arr = np.asarray(res.x, dtype=np.float64)
    if not np.all(np.isfinite(x_arr)):
        return {
            "status": "numerical_error",
            "raw_status": raw_st,
            "message": "Solution vector contains non-finite values (NaN/Inf)",
            "target_flux": math.nan,
            "achieved_biomass_growth": math.nan,
            "max_residual": math.nan,
            "max_bound_viol": math.nan,
            "n_nonzero_fluxes": 0,
            "flux_witness_sha256": None,
            "sparse_fluxes": {},
        }

    target_val = float(x_arr[target_j])
    bio_val = float(x_arr[obj_j])
    if not math.isfinite(target_val) or not math.isfinite(bio_val):
        return {
            "status": "numerical_error",
            "raw_status": raw_st,
            "message": "Target flux or biomass growth is non-finite",
            "target_flux": math.nan,
            "achieved_biomass_growth": math.nan,
            "max_residual": math.nan,
            "max_bound_viol": math.nan,
            "n_nonzero_fluxes": 0,
            "flux_witness_sha256": None,
            "sparse_fluxes": {},
        }

    residuals = S.dot(x_arr) - b_eq
    max_res = float(np.max(np.abs(residuals))) if len(residuals) > 0 else 0.0
    lb_viol = np.maximum(0.0, lb - x_arr)
    ub_viol = np.maximum(0.0, x_arr - ub)
    max_bnd_viol = float(
        max(
            np.max(lb_viol) if len(lb_viol) > 0 else 0.0,
            np.max(ub_viol) if len(ub_viol) > 0 else 0.0,
        )
    )

    if max_res > residual_tol or max_bnd_viol > residual_tol or bio_val < min_Required_growth - residual_tol:
        return {
            "status": "numerical_coupling_violation",
            "raw_status": raw_st,
            "message": (
                f"Independent residual/bound/growth validation failed: max|Sx|={max_res:.2e}, "
                f"max_bound_viol={max_bnd_viol:.2e}, bio={bio_val:.6g} vs floor={min_Required_growth:.6g}"
            ),
            "target_flux": math.nan,
            "achieved_biomass_growth": bio_val,
            "max_residual": max_res,
            "max_bound_viol": max_bnd_viol,
            "n_nonzero_fluxes": 0,
            "flux_witness_sha256": None,
            "sparse_fluxes": {},
        }

    sparse_fluxes = {
        rid: float(x_arr[j])
        for j, rid in enumerate(rxn_ids)
        if x_arr[j] != 0.0
    }
    witness_json = json.dumps(sparse_fluxes, sort_keys=True, separators=(",", ":"))
    witness_sha = hashlib.sha256(witness_json.encode("utf-8")).hexdigest()

    return {
        "status": "optimal",
        "raw_status": raw_st,
        "message": str(getattr(res, "message", "")),
        "target_flux": target_val,
        "achieved_biomass_growth": bio_val,
        "max_residual": max_res,
        "max_bound_viol": max_bnd_viol,
        "n_nonzero_fluxes": len(sparse_fluxes),
        "flux_witness_sha256": witness_sha,
        "sparse_fluxes": sparse_fluxes,
    }
