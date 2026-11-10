#!/usr/bin/env python3
"""Targeted Candidate FVA Closure & Post-Hoc Counterfactual Table Correction
(`EXP1-COUNTERFACTUAL-CLOSURE-IMPLEMENTATION.txt`, `EXP1-COUNTERFACTUAL-ROOT-SERIALIZATION-FINDING.txt`,
and `EXP1-COUNTERFACTUAL-PRELAUNCH-REVIEW.txt`).

Modes:
  1. `--mode pilot`:
     - Verifies SHA-256 of frozen `code_v2` and input artifacts.
     - Runs deterministic synthetic/toy LP verification checks (min/max interval, infeasible/unbounded
       handling without 0.0 coercion, fail-closed primary g0=NaN and finite-mismatch g0 guards,
       union unknown-EGC vs positive-ATP-cycle distinction, threshold-restoration withdrawal labels,
       and exact float64 sparse flux witness round-trip).
     - Runs 1 real `(genome, scenario_id, B2)` case from the 306 candidate deficiency cases end-to-end.
     - Writes `pilot_verification_receipt.json` and `pilot_case_instances.tsv`.
  2. `--mode shard --task-id <0..7> --num-tasks 8`:
     - Evaluates the assigned slice of the 306 `(genome, scenario_id, B2)` candidate deficiency cases
       (grouped by genome so each genome archive is loaded and hashed once).
     - Uses `mp.get_context("spawn")` worker pool with `threads=1` per HiGHS LP call.
     - Writes `shards/fva_closure_shard_0{0..7}_cases.tsv`, `shards/fva_closure_shard_0{0..7}_instances.tsv.gz`,
       `shards/fva_closure_shard_0{0..7}_witnesses.jsonl.gz`, and `shards/fva_closure_shard_0{0..7}_meta.json`.
  3. `--mode collect`:
     - Collects all 8 shards (`306` cases, `13,750` candidate deficiency instances).
     - Verifies all `306` cases and `13,750` target instances are accounted for (`100%` target accounting,
       without forcing an artificial 95%/100% LP success gate).
     - Writes:
       * `exp1_targeted_candidate_fva_306_cases_union_qc.tsv`
       * `exp1_targeted_candidate_fva_13750_instances.tsv.gz`
       * `exp1_targeted_candidate_fva_63_reactions_summary.tsv`
       * `exp1_targeted_candidate_fva_244_strata_summary.tsv`
       * `exp1_counterfactual_complete_3872_case_ledger.tsv`
       * `exp1_full484_counterfactual_wkf_corrected_v2.tsv.gz`
       * `exp1_counterfactual_closure_summary.json`
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
import multiprocessing as mp
import os
import pickle
import re
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy import sparse

THIS_DIR = Path(__file__).resolve().parent
if str(THIS_DIR) not in sys.path:
    sys.path.insert(0, str(THIS_DIR))

from counterfactual_closure_engine import (
    BIO_B0_ID,
    BIO_B2_ID,
    REVALIDATED_FEASIBLE_STATUSES,
    WKF_SCHEDULED_BIOMASS,
    WKF_SCHEDULED_SCENARIOS,
    classify_fva_interval_v2,
    classify_reaction_knockout_counterfactual_v2,
    classify_union_network_qc,
    classify_withdrawal_counterfactual_v2,
    compute_file_sha256,
    make_unverified_primary_fva_result,
    solve_verified_linear_objective_highs,
    verify_code_v2_integrity,
    verify_primary_reconstruction_gate,
)

# Globals initialized once per worker process
_WORKER_CTX: dict[str, Any] = {}


def _init_worker(code_v2_dir_str: str, pool_pkl_str: str, gapseq_out_dir_str: str, supp_models_dir_str: str) -> None:
    global _WORKER_CTX
    code_v2_dir = Path(code_v2_dir_str)
    if str(code_v2_dir) not in sys.path:
        sys.path.insert(0, str(code_v2_dir))

    import exp1_sbml_highs_engine as eng
    import stage2_pilot24_ensemble as st2

    with open(pool_pkl_str, "rb") as fh:
        pool_data = pickle.load(fh)

    _WORKER_CTX = {
        "eng": eng,
        "st2": st2,
        "universal_rxns": pool_data["reactions"],
        "universal_species": pool_data["species"],
        "gapseq_out_dir": Path(gapseq_out_dir_str),
        "supp_models_dir": Path(supp_models_dir_str),
    }


def _resolve_genome_archive_path_and_sha256(
    genome: str, gapseq_out_dir: Path, supp_models_dir: Path
) -> tuple[str, str]:
    tar_path = gapseq_out_dir / f"{genome}.tar.gz"
    if tar_path.is_file():
        return str(tar_path), compute_file_sha256(tar_path)
    supp_draft = supp_models_dir / genome / f"{genome}-draft.xml"
    if supp_draft.is_file():
        return str(supp_draft), compute_file_sha256(supp_draft)
    raise FileNotFoundError(f"Neither {tar_path} nor {supp_draft} found for genome {genome}")


def evaluate_genome_b2_cases_worker(task: dict[str, Any]) -> dict[str, Any]:
    """Evaluate all `(genome, scenario_id, B2)` cases for a single genome.
    Loads and hashes the genome archive ONCE per genome.
    """
    t0 = time.time()
    eng = _WORKER_CTX["eng"]
    st2 = _WORKER_CTX["st2"]
    universal_rxns = _WORKER_CTX["universal_rxns"]
    universal_species = _WORKER_CTX["universal_species"]
    gapseq_out_dir = _WORKER_CTX["gapseq_out_dir"]
    supp_models_dir = _WORKER_CTX["supp_models_dir"]

    genome: str = task["genome"]
    archive_path, archive_sha256 = _resolve_genome_archive_path_and_sha256(
        genome, gapseq_out_dir, supp_models_dir
    )

    payload = st2.load_genome_archive_payload(genome, gapseq_out_dir, supp_models_dir)
    draft_mod = payload["draft_model"]
    for sid, sp in draft_mod.species.items():
        if sid in universal_species:
            usp = universal_species[sid]
            if sp.formula is None and usp.formula is not None:
                sp.formula = usp.formula
            if sp.charge is None and usp.charge is not None:
                sp.charge = usp.charge

    scenarios_map: dict[str, dict[str, Any]] = dict(eng.NUTRIENT_SCENARIOS)
    if payload.get("desert_med"):
        d_med = payload["desert_med"]
        scenarios_map[eng.SCEN_S3_ANOXIC_DESERT] = {
            "description": "Per-genome gapseq desert medium (anoxic)",
            "redox_regime": "anoxic",
            "compounds": {cid: (cname, u) for cid, (cname, u) in d_med.items() if cid != "cpd00007"},
        }
        scenarios_map[eng.SCEN_S4_OXIC_DESERT] = {
            "description": "Per-genome gapseq desert medium (oxic with O2)",
            "redox_regime": "oxic",
            "compounds": {**d_med, "cpd00007": ("O2", 10.0)},
        }

    case_qc_rows: list[dict[str, Any]] = []
    instance_rows: list[dict[str, Any]] = []
    witness_records: list[dict[str, Any]] = []

    for sc_task in task["scenario_tasks"]:
        scen_id: str = sc_task["scenario_id"]
        bio_id: str = sc_task["biomass_recipe_id"]
        primary_added: list[str] = sc_task["primary_added"]
        union_added: list[str] = sc_task["union_added"]
        primary_identity_matches_wkf: bool = bool(sc_task.get("primary_identity_matches_wkf", True))
        n_accepted_solutions: int = int(sc_task["n_accepted_solutions"])
        ledger_g0: float = float(sc_task["ledger_sol0_verified_growth"])
        ledger_sol0_cost: float = float(sc_task["ledger_sol0_total_cost"])
        wkf_baseline_g0: float = float(sc_task["wkf_baseline_growth"])
        wkf_ko_map: dict[str, dict[str, Any]] = sc_task["wkf_ko_map"]
        cand_instances: list[dict[str, Any]] = sc_task["candidate_instances"]

        scen_cpds = scenarios_map[scen_id]["compounds"]
        redox_regime = scenarios_map[scen_id]["redox_regime"]

        mod_scenario = eng.apply_biomass_perturbation(draft_mod, bio_id)
        mod_scenario = eng.apply_nutrient_scenario_bounds(mod_scenario, scen_cpds)
        fixed_medium_base, _mod_super_unused, _cand_unused = eng.build_fixed_medium_supernetwork_and_base(
            mod_base=mod_scenario,
            candidate_reactions=universal_rxns,
            scenario_compounds=scen_cpds,
        )

        def _build_augmented(added_list: list[str]):
            m = fixed_medium_base.copy()
            for rid in added_list:
                crxn = universal_rxns[rid]
                for sp_id in crxn.stoichiometry:
                    if sp_id not in m.species:
                        m_sp = re.match(r"^M_(cpd\d+)_([cep]0)$", sp_id)
                        m.species[sp_id] = eng.SpeciesRecord(
                            id=sp_id,
                            cpd_id=m_sp.group(1) if m_sp else sp_id,
                            name=sp_id,
                            compartment=m_sp.group(2) if m_sp else "c0",
                            charge=None,
                            formula=None,
                        )
                m.reactions[rid] = eng.ReactionRecord(**vars(crxn))
            return m

        mod_primary = _build_augmented(primary_added)
        mod_union = _build_augmented(union_added)

        # 1. Verify mod_primary FBA and closed-uptake EGC with fail-closed gate (Prelaunch Review Item 1)
        prim_egc = eng.detect_closed_uptake_egc(mod_primary, time_limit_s=5.0)
        prim_cyc_ok, prim_cyc_st = eng._evaluate_cycle_qc_eligibility(prim_egc)
        sol0 = eng.solve_fba_highs(mod_primary, time_limit_s=5.0)
        g0 = float(sol0["growth_rate"]) if sol0["status"] == "optimal" else math.nan
        g0_abs_diff_ledger = abs(g0 - ledger_g0) if math.isfinite(g0) and math.isfinite(ledger_g0) else math.nan
        g0_abs_diff_wkf = abs(g0 - wkf_baseline_g0) if math.isfinite(g0) and math.isfinite(wkf_baseline_g0) else math.nan

        prim_verified_ok, prim_verification_status = verify_primary_reconstruction_gate(
            sol0_status=str(sol0["status"]),
            prim_cyc_ok=prim_cyc_ok,
            g0=g0,
            ledger_g0=ledger_g0,
            wkf_baseline_g0=wkf_baseline_g0,
            max_residual=float(sol0.get("max_residual", math.nan)),
            max_bound_viol=float(sol0.get("max_bound_viol", math.nan)),
            primary_identity_matches_wkf=primary_identity_matches_wkf,
        )

        # 2. Evaluate mod_union FBA and closed-uptake EGC (Prelaunch Review Item 3)
        uni_egc = eng.detect_closed_uptake_egc(mod_union, time_limit_s=5.0)
        uni_cyc_ok, uni_cyc_st = eng._evaluate_cycle_qc_eligibility(uni_egc)
        sol_uni = eng.solve_fba_highs(mod_union, time_limit_s=5.0)
        g_uni = float(sol_uni["growth_rate"]) if sol_uni["status"] == "optimal" else math.nan
        uni_fba_ok = bool(
            sol_uni["status"] == "optimal"
            and math.isfinite(g_uni)
            and g_uni >= 0.00999
            and float(sol_uni.get("max_residual", 1.0)) <= 1e-5
            and float(sol_uni.get("max_bound_viol", 1.0)) <= 1e-5
        )
        uni_qc_class = classify_union_network_qc(
            uni_fba_ok=uni_fba_ok,
            uni_fba_status=str(sol_uni["status"]),
            uni_egc=uni_egc,
            uni_cyc_ok=uni_cyc_ok,
        )

        # 3. Build stoichiometric matrix of mod_union and run verified FVA across all candidate instances
        #    FAIL-CLOSED GUARD: If primary reconstruction is NOT verified, do NOT compute max(lb, 0.95*g0)
        #    and do NOT solve FVA with a zero or unverified floor; mark every target unknown!
        sp_ids, rxn_ids, S, lb, ub, obj_j = eng.build_sparse_stoichiometric_matrix(mod_union)
        rxn_idx = {r: j for j, r in enumerate(rxn_ids)}
        if prim_verified_ok:
            min_req_growth = float(max(lb[obj_j], 0.95 * g0))
            lb_fva = lb.copy()
            lb_fva[obj_j] = min_req_growth
        else:
            min_req_growth = math.nan
            lb_fva = lb.copy()

        prim_rank_map = {r: idx for idx, r in enumerate(primary_added)}
        n_cand_case = len(cand_instances)
        n_both_opt = 0
        n_strict_1e6 = 0
        n_strict_1e12 = 0
        n_top30_checked = 0
        max_top30_abs_diff = 0.0

        for c_inst in cand_instances:
            rid = str(c_inst["rxn_id"])
            alpha_rank = prim_rank_map.get(rid, -1)
            was_top30 = bool(0 <= alpha_rank < 30)

            if not prim_verified_ok or rid not in rxn_idx:
                fail_reason = prim_verification_status if not prim_verified_ok else "target_rxn_missing_in_union"
                res_min = make_unverified_primary_fva_result(fail_reason)
                res_max = make_unverified_primary_fva_result(fail_reason)
            else:
                j_target = rxn_idx[rid]
                res_min = solve_verified_linear_objective_highs(
                    sp_ids=sp_ids,
                    rxn_ids=rxn_ids,
                    S=S,
                    lb=lb_fva,
                    ub=ub,
                    obj_j=obj_j,
                    target_j=j_target,
                    sense="min",
                    min_Required_growth=min_req_growth,
                    time_limit_s=5.0,
                )
                res_max = solve_verified_linear_objective_highs(
                    sp_ids=sp_ids,
                    rxn_ids=rxn_ids,
                    S=S,
                    lb=lb_fva,
                    ub=ub,
                    obj_j=obj_j,
                    target_j=j_target,
                    sense="max",
                    min_Required_growth=min_req_growth,
                    time_limit_s=5.0,
                )

            v_min = float(res_min["target_flux"])
            v_max = float(res_max["target_flux"])
            fva_cls = classify_fva_interval_v2(
                fva_min=v_min,
                fva_max=v_max,
                min_status=str(res_min["status"]),
                max_status=str(res_max["status"]),
            )
            if fva_cls["fva_both_optimal_and_finite"] and fva_cls["fva_interval_ordered"]:
                n_both_opt += 1
            if fva_cls["fva_strictly_nonzero_1e-6"] is True:
                n_strict_1e6 += 1
            if fva_cls["fva_strictly_nonzero_1e-12"] is True:
                n_strict_1e12 += 1

            wkf_orig = wkf_ko_map.get(rid, {})
            orig_fmin = float(wkf_orig.get("fva_min_flux", math.nan))
            orig_fmax = float(wkf_orig.get("fva_max_flux", math.nan))
            orig_ko_growth = float(wkf_orig.get("perturbed_growth_ensemble_union", math.nan))
            orig_ko_class = str(wkf_orig.get("classification", ""))

            if was_top30 and math.isfinite(orig_fmin) and math.isfinite(orig_fmax) and math.isfinite(v_min) and math.isfinite(v_max):
                n_top30_checked += 1
                diff_min = abs(v_min - orig_fmin)
                diff_max = abs(v_max - orig_fmax)
                max_top30_abs_diff = max(max_top30_abs_diff, diff_min, diff_max)
            else:
                diff_min = math.nan
                diff_max = math.nan

            witness_records.append(
                {
                    "genome": genome,
                    "scenario_id": scen_id,
                    "biomass_recipe_id": bio_id,
                    "rxn_id": rid,
                    "primary_verification_status": prim_verification_status,
                    "min_Required_growth_95pct_g0": min_req_growth,
                    "witness_representation": "exact_float64_nonzero_entries",
                    "fva_min_status": str(res_min["status"]),
                    "fva_min_flux": v_min,
                    "fva_min_witness_sha256": res_min["flux_witness_sha256"],
                    "fva_min_sparse_fluxes": res_min["sparse_fluxes"],
                    "fva_max_status": str(res_max["status"]),
                    "fva_max_flux": v_max,
                    "fva_max_witness_sha256": res_max["flux_witness_sha256"],
                    "fva_max_sparse_fluxes": res_max["sparse_fluxes"],
                }
            )

            instance_rows.append(
                {
                    "genome": genome,
                    "cohort": str(c_inst["cohort"]),
                    "canonical_order": str(c_inst["canonical_order"]),
                    "canonical_family": str(c_inst["canonical_family"]),
                    "canonical_genus": str(c_inst["canonical_genus"]),
                    "completeness": float(c_inst["completeness"]),
                    "contamination": float(c_inst["contamination"]),
                    "scenario_id": scen_id,
                    "redox_regime": redox_regime,
                    "biomass_recipe_id": bio_id,
                    "rxn_id": rid,
                    "base_id": str(c_inst["base_id"]),
                    "rxn_name": str(c_inst["rxn_name"]),
                    "is_transport": bool(c_inst["is_transport"]),
                    "candidate_weight": float(c_inst["candidate_weight"]),
                    "blast_status_in_genome": str(c_inst["blast_status_in_genome"]),
                    "blast_bitscore_in_genome": float(c_inst["blast_bitscore_in_genome"]),
                    "freq_B0_template_baseline": float(c_inst["freq_B0_template_baseline"]),
                    "freq_B2_omit_7_keep_f420": float(c_inst["freq_B2_omit_7_keep_f420"]),
                    "supernetwork_ko_necessity_B2": str(c_inst["supernetwork_ko_necessity_B2"]),
                    "biomass_sensitivity_classification": str(c_inst["biomass_sensitivity_classification"]),
                    "followup_selection_tier": str(c_inst["followup_selection_tier"]),
                    "alphabetical_rank_in_primary_added": alpha_rank,
                    "n_primary_added_rxns": len(primary_added),
                    "n_union_added_rxns": len(union_added),
                    "n_accepted_solutions_in_ensemble": n_accepted_solutions,
                    "was_in_original_top30_fva": was_top30,
                    "primary_reconstruction_verified": prim_verified_ok,
                    "primary_verification_status": prim_verification_status,
                    "primary_baseline_growth_g0": g0,
                    "fva_biomass_floor_95pct_g0": min_req_growth,
                    "union_baseline_growth": g_uni,
                    "union_egc_check_status": str(uni_egc["egc_check_status"]),
                    "union_cycle_qc_status": uni_cyc_st,
                    "union_has_spurious_atp_cycle": uni_egc["has_spurious_atp_cycle"],
                    "union_closed_uptake_atp_dissipation_flux": float(uni_egc["closed_uptake_atp_dissipation_flux"]),
                    "union_network_qc_class": uni_qc_class,
                    "wkf_union_knockout_growth": orig_ko_growth,
                    "wkf_union_knockout_class": orig_ko_class,
                    "orig_wkf_fva_min_flux": orig_fmin,
                    "orig_wkf_fva_max_flux": orig_fmax,
                    "fva_min_status": str(res_min["status"]),
                    "fva_min_raw_status": int(res_min["raw_status"]),
                    "fva_min_flux_at_95pct_opt": v_min,
                    "fva_min_achieved_biomass_growth": float(res_min["achieved_biomass_growth"]),
                    "fva_min_max_residual": float(res_min["max_residual"]),
                    "fva_min_max_bound_viol": float(res_min["max_bound_viol"]),
                    "fva_min_n_nonzero_fluxes": int(res_min["n_nonzero_fluxes"]),
                    "fva_min_flux_witness_sha256": res_min["flux_witness_sha256"],
                    "fva_max_status": str(res_max["status"]),
                    "fva_max_raw_status": int(res_max["raw_status"]),
                    "fva_max_flux_at_95pct_opt": v_max,
                    "fva_max_achieved_biomass_growth": float(res_max["achieved_biomass_growth"]),
                    "fva_max_max_residual": float(res_max["max_residual"]),
                    "fva_max_max_bound_viol": float(res_max["max_bound_viol"]),
                    "fva_max_n_nonzero_fluxes": int(res_max["n_nonzero_fluxes"]),
                    "fva_max_flux_witness_sha256": res_max["flux_witness_sha256"],
                    "fva_both_optimal_and_finite": fva_cls["fva_both_optimal_and_finite"],
                    "fva_interval_ordered": fva_cls["fva_interval_ordered"],
                    "fva_span": fva_cls["fva_span"],
                    "fva_strictly_nonzero_1e-6": fva_cls["fva_strictly_nonzero_1e-6"],
                    "fva_strictly_nonzero_1e-12": fva_cls["fva_strictly_nonzero_1e-12"],
                    "fva_interval_class": fva_cls["fva_interval_class"],
                    "top30_abs_diff_fva_min": diff_min,
                    "top30_abs_diff_fva_max": diff_max,
                    "genome_archive_sha256": archive_sha256,
                }
            )

        case_qc_rows.append(
            {
                "genome": genome,
                "scenario_id": scen_id,
                "redox_regime": redox_regime,
                "biomass_recipe_id": bio_id,
                "n_accepted_solutions": n_accepted_solutions,
                "n_primary_added_rxns": len(primary_added),
                "n_union_added_rxns": len(union_added),
                "n_extra_union_rxns_beyond_primary": len(union_added) - len(primary_added),
                "primary_identity_matches_wkf": primary_identity_matches_wkf,
                "ledger_sol0_total_cost": ledger_sol0_cost,
                "ledger_sol0_verified_growth": ledger_g0,
                "wkf_baseline_growth": wkf_baseline_g0,
                "reconstructed_primary_fba_status": str(sol0["status"]),
                "reconstructed_primary_g0": g0,
                "primary_g0_abs_diff_vs_ledger": g0_abs_diff_ledger,
                "primary_g0_abs_diff_vs_wkf": g0_abs_diff_wkf,
                "primary_fba_max_residual": float(sol0.get("max_residual", math.nan)),
                "primary_fba_max_bound_viol": float(sol0.get("max_bound_viol", math.nan)),
                "primary_egc_check_status": str(prim_egc["egc_check_status"]),
                "primary_has_spurious_atp_cycle": prim_egc["has_spurious_atp_cycle"],
                "primary_closed_uptake_atp_dissipation_flux": float(prim_egc["closed_uptake_atp_dissipation_flux"]),
                "primary_cycle_qc_status": prim_cyc_st,
                "primary_reconstruction_verified": prim_verified_ok,
                "primary_verification_status": prim_verification_status,
                "union_fba_status": str(sol_uni["status"]),
                "union_fba_growth": g_uni,
                "union_fba_max_residual": float(sol_uni.get("max_residual", math.nan)),
                "union_fba_max_bound_viol": float(sol_uni.get("max_bound_viol", math.nan)),
                "union_egc_check_status": str(uni_egc["egc_check_status"]),
                "union_has_spurious_atp_cycle": uni_egc["has_spurious_atp_cycle"],
                "union_closed_uptake_atp_dissipation_flux": float(uni_egc["closed_uptake_atp_dissipation_flux"]),
                "union_cycle_qc_status": uni_cyc_st,
                "union_cycle_qc_passed": uni_cyc_ok,
                "union_network_qc_class": uni_qc_class,
                "n_candidate_deficiency_instances": n_cand_case,
                "n_fva_both_optimal_and_ordered": n_both_opt,
                "n_fva_strictly_nonzero_1e-6": n_strict_1e6,
                "n_fva_strictly_nonzero_1e-12": n_strict_1e12,
                "n_top30_checked_against_orig_wkf": n_top30_checked,
                "max_top30_abs_diff_vs_orig_wkf": max_top30_abs_diff,
                "genome_archive_path": archive_path,
                "genome_archive_sha256": archive_sha256,
                "elapsed_s": round(time.time() - t0, 3),
            }
        )

    return {
        "genome": genome,
        "case_qc_rows": case_qc_rows,
        "instance_rows": instance_rows,
        "witness_records": witness_records,
    }


def build_306_case_tasks_by_genome(
    results_dir: Path,
) -> tuple[list[dict[str, Any]], dict[str, str]]:
    """Build deterministic per-genome task payloads for all 306 `(genome, scenario_id, B2)` cases
    and `13,750` candidate deficiency instances, and compute input SHA-256 digests ONCE.
    """
    cand_inst_path = results_dir / "exp1_stage4_supported_candidate_deficiencies_instances.tsv"
    ledger_path = results_dir / "exp1_full484_per_solution_ledger.tsv"
    wkf_path = results_dir / "exp1_full484_withdrawal_knockout_fva.tsv.gz"
    pool_pkl_path = results_dir / "exp1_candidate_reaction_pool.pkl"

    input_hashes = {
        "exp1_stage4_supported_candidate_deficiencies_instances.tsv": compute_file_sha256(cand_inst_path),
        "exp1_full484_per_solution_ledger.tsv": compute_file_sha256(ledger_path),
        "exp1_full484_withdrawal_knockout_fva.tsv.gz": compute_file_sha256(wkf_path),
        "exp1_candidate_reaction_pool.pkl": compute_file_sha256(pool_pkl_path),
    }

    cand_df = pd.read_csv(cand_inst_path, sep="\t")
    if len(cand_df) != 13750:
        raise RuntimeError(f"Expected 13,750 candidate deficiency instances, got {len(cand_df)}")

    target_genomes = sorted(cand_df["genome"].unique())
    if len(target_genomes) != 81:
        raise RuntimeError(f"Expected 81 unique HQ EQ genomes in candidate instances, got {len(target_genomes)}")

    # Filter ledger to accepted_feasible B2 rows for those 81 genomes
    ledger_df = pd.read_csv(ledger_path, sep="\t")
    ledger_b2 = ledger_df[
        (ledger_df["genome"].isin(target_genomes))
        & (ledger_df["biomass_recipe_id"] == BIO_B2_ID)
        & (ledger_df["validity_disposition"] == "accepted_feasible")
    ].copy()

    # Filter wkf to B2 reaction_knockout_and_fva rows for those 81 genomes
    wkf_df = pd.read_csv(wkf_path, sep="\t", low_memory=False)
    wkf_b2_ko = wkf_df[
        (wkf_df["genome"].isin(target_genomes))
        & (wkf_df["biomass_recipe_id"] == BIO_B2_ID)
        & (wkf_df["record_type"] == "reaction_knockout_and_fva")
    ].copy()

    # Index ledger by (genome, scenario_id)
    ledger_by_case: dict[tuple[str, str], pd.DataFrame] = {
        (str(g), str(s)): grp.sort_values("solution_index")
        for (g, s), grp in ledger_b2.groupby(["genome", "scenario_id"], sort=False)
    }
    wkf_by_case: dict[tuple[str, str], pd.DataFrame] = {
        (str(g), str(s)): grp
        for (g, s), grp in wkf_b2_ko.groupby(["genome", "scenario_id"], sort=False)
    }
    cand_by_case: dict[tuple[str, str], pd.DataFrame] = {
        (str(g), str(s)): grp.sort_values("rxn_id")
        for (g, s), grp in cand_df.groupby(["genome", "scenario_id"], sort=False)
    }

    if len(cand_by_case) != 306:
        raise RuntimeError(f"Expected 306 unique (genome, scenario_id) cases, got {len(cand_by_case)}")

    tasks_by_genome: list[dict[str, Any]] = []
    for g in target_genomes:
        scen_tasks: list[dict[str, Any]] = []
        for s in WKF_SCHEDULED_SCENARIOS:
            if (g, s) not in cand_by_case:
                continue
            c_grp = cand_by_case[(g, s)]
            l_grp = ledger_by_case.get((g, s), pd.DataFrame())
            w_grp = wkf_by_case.get((g, s), pd.DataFrame())

            sol0_rows = l_grp[l_grp["solution_index"] == 0] if not l_grp.empty else pd.DataFrame()
            if len(sol0_rows) == 1:
                sol0 = sol0_rows.iloc[0]
                primary_added = sorted(
                    r for r in str(sol0["selected_reaction_ids"]).split(";") if r
                )
                ledger_g0 = float(sol0["verified_growth"])
                ledger_cost = float(sol0["total_cost"])
            else:
                primary_added = []
                ledger_g0 = math.nan
                ledger_cost = math.nan

            union_set: set[str] = set()
            if not l_grp.empty:
                for r_str in l_grp["selected_reaction_ids"].dropna():
                    union_set.update(r for r in str(r_str).split(";") if r)
            union_added = sorted(union_set)

            wkf_rids = sorted(w_grp["target_id"].astype(str).tolist()) if not w_grp.empty else []
            primary_identity_matches_wkf = bool(primary_added == wkf_rids and len(primary_added) > 0)

            wkf_ko_map = {
                str(r["target_id"]): {
                    "perturbed_growth_ensemble_union": float(r["perturbed_growth_ensemble_union"]),
                    "classification": str(r["classification"]),
                    "fva_min_flux": float(r["fva_min_flux"]) if pd.notna(r["fva_min_flux"]) else math.nan,
                    "fva_max_flux": float(r["fva_max_flux"]) if pd.notna(r["fva_max_flux"]) else math.nan,
                }
                for _, r in w_grp.iterrows()
            }

            scen_tasks.append(
                {
                    "scenario_id": s,
                    "biomass_recipe_id": BIO_B2_ID,
                    "primary_added": primary_added,
                    "union_added": union_added,
                    "primary_identity_matches_wkf": primary_identity_matches_wkf,
                    "n_accepted_solutions": len(l_grp),
                    "ledger_sol0_verified_growth": ledger_g0,
                    "ledger_sol0_total_cost": ledger_cost,
                    "wkf_baseline_growth": float(w_grp["baseline_growth"].iloc[0]) if not w_grp.empty else math.nan,
                    "wkf_ko_map": wkf_ko_map,
                    "candidate_instances": c_grp.to_dict(orient="records"),
                }
            )
        tasks_by_genome.append({"genome": g, "scenario_tasks": scen_tasks})

    return tasks_by_genome, input_hashes


def run_synthetic_toy_checks() -> dict[str, Any]:
    """Run deterministic synthetic LP and classification checks before any real case evaluation
    (including Prelaunch Review Items 1, 2, 3, 4).
    """
    # Toy LP: S v = 0 with 3 reactions: R_in (-> A), R_mid (A -> B), R_bio (B ->)
    sp_ids = ["M_A", "M_B"]
    rxn_ids = ["R_in", "R_mid", "R_bio"]
    S = sparse.csc_matrix(
        np.array(
            [
                [1.0, -1.0, 0.0],
                [0.0, 1.0, -1.0],
            ],
            dtype=np.float64,
        )
    )
    lb = np.array([0.0, 0.0, 0.0], dtype=np.float64)
    ub = np.array([10.0, 10.0, 10.0], dtype=np.float64)
    obj_j = 2
    min_req = 0.95 * 10.0
    lb_fva = lb.copy()
    lb_fva[obj_j] = min_req

    r_min = solve_verified_linear_objective_highs(
        sp_ids, rxn_ids, S, lb_fva, ub, obj_j=2, target_j=1, sense="min", min_Required_growth=min_req
    )
    r_max = solve_verified_linear_objective_highs(
        sp_ids, rxn_ids, S, lb_fva, ub, obj_j=2, target_j=1, sense="max", min_Required_growth=min_req
    )
    assert r_min["status"] == "optimal" and abs(r_min["target_flux"] - 9.5) <= 1e-7
    assert r_max["status"] == "optimal" and abs(r_max["target_flux"] - 10.0) <= 1e-7
    f_cls = classify_fva_interval_v2(r_min["target_flux"], r_max["target_flux"], r_min["status"], r_max["status"])
    assert f_cls["fva_both_optimal_and_finite"] is True
    assert f_cls["fva_strictly_nonzero_1e-6"] is True

    # Prelaunch Item 1a: Failed primary solve (g0 = NaN) MUST fail closed and never use lb=0
    ok_nan, st_nan = verify_primary_reconstruction_gate(
        sol0_status="infeasible",
        prim_cyc_ok=True,
        g0=math.nan,
        ledger_g0=0.05,
        wkf_baseline_g0=0.05,
        max_residual=0.0,
        max_bound_viol=0.0,
    )
    assert ok_nan is False and st_nan == "unverified_primary_fba_solver_infeasible"
    # Even if someone passed min_Required_growth=NaN and lb[obj_j]=0.0 into solve_verified_linear_objective_highs:
    r_nan_floor = solve_verified_linear_objective_highs(
        sp_ids, rxn_ids, S, lb, ub, obj_j=2, target_j=1, sense="min", min_Required_growth=math.nan
    )
    assert r_nan_floor["status"].startswith("unknown_") and math.isnan(r_nan_floor["target_flux"])

    # Prelaunch Item 1b: Finite g0 that mismatches ledger_g0 MUST fail closed
    ok_mis, st_mis = verify_primary_reconstruction_gate(
        sol0_status="optimal",
        prim_cyc_ok=True,
        g0=0.050000,
        ledger_g0=0.052000,
        wkf_baseline_g0=0.052000,
        max_residual=0.0,
        max_bound_viol=0.0,
    )
    assert ok_mis is False and st_mis == "unverified_primary_g0_mismatch_vs_solution_ledger"
    res_unv = make_unverified_primary_fva_result(st_mis)
    cls_unv = classify_fva_interval_v2(
        res_unv["target_flux"], res_unv["target_flux"], res_unv["status"], res_unv["status"]
    )
    assert cls_unv["fva_both_optimal_and_finite"] is False
    assert cls_unv["fva_strictly_nonzero_1e-6"] is None

    # Infeasible toy LP check -> must return infeasible/NaN/None, NEVER 0.0 or False
    lb_infeas = np.array([0.0, 0.0, 20.0], dtype=np.float64)
    r_inf = solve_verified_linear_objective_highs(
        sp_ids, rxn_ids, S, lb_infeas, ub, obj_j=2, target_j=1, sense="min", min_Required_growth=20.0
    )
    assert r_inf["status"] == "infeasible" and math.isnan(r_inf["target_flux"])
    f_inf_cls = classify_fva_interval_v2(r_inf["target_flux"], r_max["target_flux"], r_inf["status"], r_max["status"])
    assert f_inf_cls["fva_both_optimal_and_finite"] is False
    assert f_inf_cls["fva_strictly_nonzero_1e-6"] is None

    # Prelaunch Item 2: Withdrawal classification edge cases (threshold restoration + continuous ratio)
    w_gap = classify_withdrawal_counterfactual_v2(0.05, 0.0, 0.002, 0.01)
    assert w_gap["primary_growth_lost_lt_1e-4"] is True
    assert w_gap["corrected_classification"] == "primary_growth_lost_lt_1e-4__partial_subthreshold_union_growth_1e-4_to_0.005"
    assert w_gap["addback_rescues_growth_v2"] is True
    assert w_gap["addback_rescue_disposition"] == "restored_primary_growth_above_declared_threshold_0.99_min_growth"
    assert abs(w_gap["addback_to_baseline_growth_ratio"] - 0.2) <= 1e-12

    w_ret = classify_withdrawal_counterfactual_v2(0.05, 0.04, 0.05, 0.05)
    assert w_ret["primary_growth_lost_lt_1e-4"] is False
    assert w_ret["addback_rescues_growth_v2"] is None
    assert w_ret["addback_rescue_disposition"] == "not_applicable_primary_growth_retained_ge_1e-4"

    # Prelaunch Item 3: Union QC class distinguishes pass vs positive ATP cycle vs unknown EGC
    u_pass = classify_union_network_qc(
        True, "optimal", {"has_spurious_atp_cycle": False, "egc_check_status": "evaluated_optimal"}, True
    )
    u_fail = classify_union_network_qc(
        True, "optimal", {"has_spurious_atp_cycle": True, "egc_check_status": "evaluated_optimal"}, False
    )
    u_unk = classify_union_network_qc(
        True, "optimal", {"has_spurious_atp_cycle": None, "egc_check_status": "unknown_timeout"}, False
    )
    assert u_pass == "mathematical_union_passed_declared_numerical_and_closed_uptake_egc_checks"
    assert u_fail == "mathematical_union_failed_closed_uptake_atp_cycle_check"
    assert u_unk == "mathematical_union_unknown_closed_uptake_egc_check__unknown_timeout"

    # Knockout classification -> unperformed primary/addback fields are NaN/None
    ko_c = classify_reaction_knockout_counterfactual_v2(0.05, 0.0)
    assert math.isnan(ko_c["perturbed_growth_primary_v2"])
    assert math.isnan(ko_c["addback_rescue_growth_v2"])
    assert ko_c["addback_rescues_growth_v2"] is None
    assert ko_c["primary_perturbation_measurement_status"] == "not_performed_in_executed_code_v2"

    return {
        "toy_lp_coupled_interval": [r_min["target_flux"], r_max["target_flux"]],
        "toy_infeasible_lp_status": r_inf["status"],
        "toy_nan_primary_floor_status": r_nan_floor["status"],
        "toy_mismatch_primary_status": st_mis,
        "union_qc_classes": {"pass": u_pass, "cycle_fail": u_fail, "unknown_egc": u_unk},
        "withdrawal_partial_union_class": w_gap["corrected_classification"],
        "withdrawal_threshold_restoration_disposition": w_gap["addback_rescue_disposition"],
        "withdrawal_retained_addback_disposition": w_ret["addback_rescue_disposition"],
        "ko_unperformed_status": ko_c["primary_perturbation_measurement_status"],
        "all_toy_assertions_passed": True,
    }


def run_pilot_mode(args: argparse.Namespace) -> None:
    t0 = time.time()
    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    code_hashes = verify_code_v2_integrity(args.code_v2_dir)
    toy_res = run_synthetic_toy_checks()
    tasks_by_genome, input_hashes = build_306_case_tasks_by_genome(args.results_dir)

    # Pick the first genome and its first scenario case for the 1-case pilot
    first_g_task = tasks_by_genome[0]
    expected_pilot_targets = len(first_g_task["scenario_tasks"][0]["candidate_instances"])
    pilot_task = {
        "genome": first_g_task["genome"],
        "scenario_tasks": [first_g_task["scenario_tasks"][0]],
    }

    _init_worker(
        str(args.code_v2_dir),
        str(args.results_dir / "exp1_candidate_reaction_pool.pkl"),
        str(args.gapseq_out_dir),
        str(args.supp_models_dir),
    )
    res = evaluate_genome_b2_cases_worker(pilot_task)
    c_row = res["case_qc_rows"][0]
    inst_df = pd.DataFrame(res["instance_rows"])

    # Require 100% target accounting (no silent drops)
    if len(inst_df) != expected_pilot_targets:
        raise RuntimeError(f"Pilot target accounting mismatch: expected {expected_pilot_targets}, got {len(inst_df)}")

    inst_path = out_dir / "pilot_case_instances.tsv"
    inst_df.to_csv(inst_path, sep="\t", index=False)

    receipt = {
        "mode": "pilot",
        "elapsed_s": round(time.time() - t0, 3),
        "code_v2_sha256": code_hashes,
        "input_sha256": input_hashes,
        "synthetic_toy_checks": toy_res,
        "pilot_case_qc": c_row,
        "pilot_expected_targets": expected_pilot_targets,
        "pilot_accounted_targets": len(inst_df),
        "pilot_fva_interval_class_counts": {
            str(k): int(v) for k, v in inst_df["fva_interval_class"].value_counts().items()
        },
        "pilot_passed_all_gates": True,
    }
    receipt_path = out_dir / "pilot_verification_receipt.json"
    receipt_path.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print(f"[Pilot Complete in {receipt['elapsed_s']}s] Wrote {receipt_path}", flush=True)


def run_shard_mode(args: argparse.Namespace) -> None:
    t0 = time.time()
    shard_dir = args.out_dir / "shards"
    shard_dir.mkdir(parents=True, exist_ok=True)

    code_hashes = verify_code_v2_integrity(args.code_v2_dir)
    tasks_by_genome, input_hashes = build_306_case_tasks_by_genome(args.results_dir)
    assigned_tasks = tasks_by_genome[args.task_id :: args.num_tasks]

    n_cases = sum(len(t["scenario_tasks"]) for t in assigned_tasks)
    n_insts = sum(
        len(st["candidate_instances"])
        for t in assigned_tasks
        for st in t["scenario_tasks"]
    )
    print(
        f"[FVA Closure Shard {args.task_id}/{args.num_tasks}] Assigned {len(assigned_tasks)} genomes, "
        f"{n_cases} (genome, scenario, B2) cases, {n_insts} candidate FVA instances",
        flush=True,
    )

    pool_pkl_str = str(args.results_dir / "exp1_candidate_reaction_pool.pkl")
    max_workers = min(args.workers, max(1, len(assigned_tasks)), os.cpu_count() or 16)
    ctx = mp.get_context("spawn")

    all_case_rows: list[dict[str, Any]] = []
    all_inst_rows: list[dict[str, Any]] = []
    all_witnesses: list[dict[str, Any]] = []

    with ProcessPoolExecutor(
        max_workers=max_workers,
        mp_context=ctx,
        initializer=_init_worker,
        initargs=(
            str(args.code_v2_dir),
            pool_pkl_str,
            str(args.gapseq_out_dir),
            str(args.supp_models_dir),
        ),
    ) as ex:
        futs = {ex.submit(evaluate_genome_b2_cases_worker, t): t["genome"] for t in assigned_tasks}
        for idx, fut in enumerate(as_completed(futs), 1):
            g_name = futs[fut]
            res = fut.result()
            all_case_rows.extend(res["case_qc_rows"])
            all_inst_rows.extend(res["instance_rows"])
            all_witnesses.extend(res["witness_records"])
            print(
                f"[FVA Closure Shard {args.task_id}: {idx}/{len(assigned_tasks)}] Completed {g_name} "
                f"({len(res['case_qc_rows'])} cases, {len(res['instance_rows'])} instances)",
                flush=True,
            )

    tag = f"fva_closure_shard_{args.task_id:02d}"
    cases_df = pd.DataFrame(all_case_rows).sort_values(["genome", "scenario_id"]).reset_index(drop=True)
    inst_df = pd.DataFrame(all_inst_rows).sort_values(["genome", "scenario_id", "rxn_id"]).reset_index(drop=True)

    if len(cases_df) != n_cases or len(inst_df) != n_insts:
        raise RuntimeError(
            f"Shard {args.task_id} target accounting mismatch: cases {len(cases_df)}/{n_cases}, "
            f"instances {len(inst_df)}/{n_insts}"
        )

    cases_path = shard_dir / f"{tag}_cases.tsv"
    inst_path = shard_dir / f"{tag}_instances.tsv.gz"
    wit_path = shard_dir / f"{tag}_witnesses.jsonl.gz"
    meta_path = shard_dir / f"{tag}_meta.json"

    cases_df.to_csv(cases_path, sep="\t", index=False)
    inst_df.to_csv(inst_path, sep="\t", index=False, compression="gzip")

    all_witnesses.sort(key=lambda r: (r["genome"], r["scenario_id"], r["rxn_id"]))
    with gzip.open(wit_path, "wt", encoding="utf-8") as gz:
        for w in all_witnesses:
            gz.write(json.dumps(w, sort_keys=True, separators=(",", ":")) + "\n")

    meta = {
        "task_id": args.task_id,
        "num_tasks": args.num_tasks,
        "n_genomes": len(assigned_tasks),
        "n_cases_expected": n_cases,
        "n_cases_accounted": len(cases_df),
        "n_instances_expected": n_insts,
        "n_instances_accounted": len(inst_df),
        "n_primary_verified": int(cases_df["primary_reconstruction_verified"].sum()),
        "n_primary_unverified": int((~cases_df["primary_reconstruction_verified"]).sum()),
        "n_fva_both_optimal_and_ordered": int(cases_df["n_fva_both_optimal_and_ordered"].sum()),
        "max_top30_abs_diff_vs_orig_wkf": float(cases_df["max_top30_abs_diff_vs_orig_wkf"].max()),
        "cases_tsv_sha256": compute_file_sha256(cases_path),
        "instances_tsv_gz_sha256": compute_file_sha256(inst_path),
        "witnesses_jsonl_gz_sha256": compute_file_sha256(wit_path),
        "code_v2_sha256": code_hashes,
        "input_sha256": input_hashes,
        "elapsed_s": round(time.time() - t0, 3),
    }
    meta_path.write_text(json.dumps(meta, indent=2) + "\n", encoding="utf-8")
    print(f"[FVA Closure Shard {args.task_id} Complete in {meta['elapsed_s']}s] Wrote {meta_path}", flush=True)


def build_complete_3872_case_ledger(
    results_dir: Path,
    wkf_df: pd.DataFrame,
    out_path: Path,
) -> dict[str, Any]:
    """Build `exp1_counterfactual_complete_3872_case_ledger.tsv` covering all 3,872 scheduled
    `(genome, scenario_id, biomass_recipe_id)` cases (`484 x 4 x 2`), including the 9 non-evaluable cases.
    """
    scen_sum_path = results_dir / "exp1_full484_ensemble_summary.tsv.gz"
    if not scen_sum_path.is_file():
        scen_sum_path = results_dir / "exp1_full484_scenario_summary.tsv"
    scen_sum_df = pd.read_csv(scen_sum_path, sep="\t")
    sched_df = scen_sum_df[
        scen_sum_df["scenario_id"].isin(WKF_SCHEDULED_SCENARIOS)
        & scen_sum_df["biomass_recipe_id"].isin(WKF_SCHEDULED_BIOMASS)
    ].copy()
    if len(sched_df) != 3872:
        raise RuntimeError(f"Expected 3,872 scheduled scenario_summary rows, got {len(sched_df)}")

    w_sub = wkf_df[wkf_df["record_type"] == "nutrient_withdrawal_addback"]
    k_sub = wkf_df[wkf_df["record_type"] == "reaction_knockout_and_fva"]

    w_cnt = w_sub.groupby(["genome", "scenario_id", "biomass_recipe_id"]).size().to_dict()
    k_cnt = k_sub.groupby(["genome", "scenario_id", "biomass_recipe_id"]).size().to_dict()
    f_cnt = (
        k_sub[k_sub["fva_min_flux"].notna() & k_sub["fva_max_flux"].notna()]
        .groupby(["genome", "scenario_id", "biomass_recipe_id"])
        .size()
        .to_dict()
    )
    g0_map = (
        wkf_df.groupby(["genome", "scenario_id", "biomass_recipe_id"])["baseline_growth"]
        .first()
        .to_dict()
    )

    rows: list[dict[str, Any]] = []
    for _, r in sched_df.sort_values(["genome", "scenario_id", "biomass_recipe_id"]).iterrows():
        key = (str(r["genome"]), str(r["scenario_id"]), str(r["biomass_recipe_id"]))
        solver_st = str(r["solver_status"])
        is_accepted = solver_st in REVALIDATED_FEASIBLE_STATUSES
        n_w = int(w_cnt.get(key, 0))
        n_k = int(k_cnt.get(key, 0))
        n_f = int(f_cnt.get(key, 0))
        in_wkf = (n_w + n_k) > 0

        if is_accepted and in_wkf:
            case_disp = "evaluated_feasible_baseline"
        elif not is_accepted and not in_wkf:
            case_disp = f"non_evaluable_unresolved_{solver_st}"
        else:
            case_disp = f"unexpected_case_state_accepted_{is_accepted}_in_wkf_{in_wkf}"

        rows.append(
            {
                "genome": key[0],
                "cohort": str(r["cohort"]),
                "canonical_order": str(r["canonical_order"]),
                "canonical_family": str(r["canonical_family"]),
                "canonical_genus": str(r["canonical_genus"]),
                "completeness": float(r["completeness"]) if pd.notna(r["completeness"]) else math.nan,
                "contamination": float(r["contamination"]) if pd.notna(r["contamination"]) else math.nan,
                "scenario_id": key[1],
                "redox_regime": str(r["redox_regime"]),
                "biomass_recipe_id": key[2],
                "solver_status": solver_st,
                "enumeration_completeness": str(r["enumeration_completeness"]),
                "pool_exhaustion_type": str(r["pool_exhaustion_type"]),
                "n_solutions_attempted": int(r["n_solutions_attempted"]),
                "n_optimization_feasible_solutions": int(r["n_optimization_feasible_solutions"]),
                "n_cycle_qc_rejected_solutions": int(r["n_cycle_qc_rejected_solutions"]),
                "n_distinct_valid_solutions": int(r["n_distinct_valid_solutions"]),
                "fewer_than_target_reason": str(r["fewer_than_target_reason"]),
                "first_solve_mip_gap": float(r["first_solve_mip_gap"]) if pd.notna(r["first_solve_mip_gap"]) else math.nan,
                "first_solve_dual_bound": float(r["first_solve_dual_bound"]) if pd.notna(r["first_solve_dual_bound"]) else math.nan,
                "in_executed_wkf_table": in_wkf,
                "wkf_case_disposition": case_disp,
                "wkf_baseline_growth": float(g0_map.get(key, math.nan)),
                "n_wkf_total_rows": n_w + n_k,
                "n_wkf_withdrawal_rows": n_w,
                "n_wkf_knockout_rows": n_k,
                "n_wkf_original_finite_fva_rows": n_f,
                "n_wkf_original_unassayed_fva_rows": n_k - n_f,
            }
        )

    out_df = pd.DataFrame(rows)
    out_df.to_csv(out_path, sep="\t", index=False)
    return {
        "n_scheduled_cases_total": len(out_df),
        "n_evaluated_feasible_cases": int((out_df["in_executed_wkf_table"] == True).sum()),
        "n_non_evaluable_cases": int((out_df["in_executed_wkf_table"] == False).sum()),
        "disposition_counts": {
            str(k): int(v) for k, v in out_df["wkf_case_disposition"].value_counts().items()
        },
        "non_evaluable_records": out_df[out_df["in_executed_wkf_table"] == False][
            [
                "genome",
                "cohort",
                "canonical_family",
                "scenario_id",
                "biomass_recipe_id",
                "solver_status",
                "n_solutions_attempted",
                "n_optimization_feasible_solutions",
                "n_cycle_qc_rejected_solutions",
                "n_distinct_valid_solutions",
                "fewer_than_target_reason",
                "wkf_case_disposition",
            ]
        ].to_dict(orient="records"),
    }


def build_corrected_wkf_v2_table(
    wkf_df: pd.DataFrame,
    cand_inst_df: pd.DataFrame,
    out_path: Path,
) -> dict[str, Any]:
    """Create `exp1_full484_counterfactual_wkf_corrected_v2.tsv.gz` (`692,922` rows)."""
    df = wkf_df.copy()
    if len(df) != 692922:
        raise RuntimeError(f"Expected 692,922 rows in wkf_df, got {len(df)}")

    df["raw_perturbed_growth_primary"] = df["perturbed_growth_primary"]
    df["raw_addback_rescue_growth"] = df["addback_rescue_growth"]
    df["raw_addback_rescues_growth"] = df["addback_rescues_growth"]
    df["raw_classification"] = df["classification"]
    df["raw_fva_min_flux"] = df["fva_min_flux"]
    df["raw_fva_max_flux"] = df["fva_max_flux"]

    is_w = df["record_type"] == "nutrient_withdrawal_addback"
    is_k = df["record_type"] == "reaction_knockout_and_fva"

    g0 = df["baseline_growth"].to_numpy(dtype=np.float64)
    gp = df["raw_perturbed_growth_primary"].to_numpy(dtype=np.float64)
    gu = df["perturbed_growth_ensemble_union"].to_numpy(dtype=np.float64)
    ga = df["raw_addback_rescue_growth"].to_numpy(dtype=np.float64)

    n = len(df)
    pert_prim_v2 = np.full(n, np.nan, dtype=np.float64)
    addback_gr_v2 = np.full(n, np.nan, dtype=np.float64)
    addback_ratio_v2 = np.full(n, np.nan, dtype=np.float64)
    addback_flag_v2 = np.full(n, None, dtype=object)
    prim_lost_v2 = np.full(n, None, dtype=object)
    prim_meas_st = np.full(n, "not_performed_in_executed_code_v2", dtype=object)
    addback_meas_st = np.full(n, "not_performed_in_executed_code_v2", dtype=object)
    addback_disp = np.full(n, "not_performed_in_executed_code_v2", dtype=object)
    corr_cls = np.full(n, "", dtype=object)
    uni_regime = np.full(n, "", dtype=object)

    w_indices = np.flatnonzero(is_w.to_numpy())
    for i in w_indices:
        pert_prim_v2[i] = gp[i]
        addback_gr_v2[i] = ga[i]
        prim_meas_st[i] = "measured_fba_on_primary_sol0_network"
        addback_meas_st[i] = "measured_fba_single_metabolite_addback_on_primary_sol0_network"
        c_res = classify_withdrawal_counterfactual_v2(
            baseline_growth=float(g0[i]),
            withdrawn_growth_primary=float(gp[i]),
            withdrawn_growth_union=float(gu[i]),
            addback_rescue_growth=float(ga[i]),
        )
        prim_lost_v2[i] = c_res["primary_growth_lost_lt_1e-4"]
        addback_flag_v2[i] = c_res["addback_rescues_growth_v2"]
        addback_disp[i] = c_res["addback_rescue_disposition"]
        addback_ratio_v2[i] = float(c_res["addback_to_baseline_growth_ratio"])
        corr_cls[i] = c_res["corrected_classification"]
        uni_regime[i] = c_res["union_growth_regime"]

    k_indices = np.flatnonzero(is_k.to_numpy())
    for i in k_indices:
        c_res = classify_reaction_knockout_counterfactual_v2(
            baseline_growth=float(g0[i]),
            ko_growth_union=float(gu[i]),
        )
        corr_cls[i] = c_res["corrected_classification"]
        uni_regime[i] = c_res["union_growth_regime"]

    df["perturbed_growth_primary"] = pert_prim_v2
    df["perturbed_growth_primary_v2"] = pert_prim_v2
    df["primary_perturbation_measurement_status"] = prim_meas_st
    df["primary_growth_lost_lt_1e-4"] = prim_lost_v2
    df["addback_rescue_growth"] = addback_gr_v2
    df["addback_rescue_growth_v2"] = addback_gr_v2
    df["addback_to_baseline_growth_ratio"] = addback_ratio_v2
    df["addback_measurement_status"] = addback_meas_st
    df["addback_rescues_growth"] = addback_flag_v2
    df["addback_rescues_growth_v2"] = addback_flag_v2
    df["addback_rescue_disposition"] = addback_disp
    df["classification"] = corr_cls
    df["corrected_classification"] = corr_cls
    df["union_growth_regime"] = uni_regime

    cand_fva_map: dict[tuple[str, str, str, str], tuple[float, float, str]] = {}
    for _, r in cand_inst_df.iterrows():
        key = (
            str(r["genome"]),
            str(r["scenario_id"]),
            str(r["biomass_recipe_id"]),
            str(r["rxn_id"]),
        )
        cand_fva_map[key] = (
            float(r["fva_min_flux_at_95pct_opt"]) if pd.notna(r["fva_min_flux_at_95pct_opt"]) else math.nan,
            float(r["fva_max_flux_at_95pct_opt"]) if pd.notna(r["fva_max_flux_at_95pct_opt"]) else math.nan,
            str(r["fva_interval_class"]),
        )

    fmin_v2 = df["raw_fva_min_flux"].to_numpy(dtype=np.float64).copy()
    fmax_v2 = df["raw_fva_max_flux"].to_numpy(dtype=np.float64).copy()
    fva_src = np.full(n, "not_applicable_withdrawal_row", dtype=object)
    fva_icls = np.full(n, "not_applicable_withdrawal_row", dtype=object)

    genomes_arr = df["genome"].astype(str).to_numpy()
    scens_arr = df["scenario_id"].astype(str).to_numpy()
    bios_arr = df["biomass_recipe_id"].astype(str).to_numpy()
    tids_arr = df["target_id"].astype(str).to_numpy()

    n_cand_matched = 0
    for i in k_indices:
        orig_has_fva = math.isfinite(fmin_v2[i]) and math.isfinite(fmax_v2[i])
        key = (genomes_arr[i], scens_arr[i], bios_arr[i], tids_arr[i])
        if key in cand_fva_map:
            n_cand_matched += 1
            c_fmin, c_fmax, c_icls = cand_fva_map[key]
            fmin_v2[i] = c_fmin
            fmax_v2[i] = c_fmax
            fva_icls[i] = c_icls
            fva_src[i] = (
                "original_top30_and_verified_targeted_candidate_fva_closure"
                if orig_has_fva
                else "targeted_candidate_fva_closure_beyond_top30"
            )
        elif orig_has_fva:
            c_res = classify_fva_interval_v2(float(fmin_v2[i]), float(fmax_v2[i]))
            fva_icls[i] = c_res["fva_interval_class"]
            fva_src[i] = "original_executed_top30_fva"
        else:
            fva_icls[i] = "not_assayed_outside_top30_non_candidate_rxn"
            fva_src[i] = "not_assayed_outside_top30_non_candidate_rxn"

    if n_cand_matched != 13750:
        raise RuntimeError(f"Expected 13,750 candidate deficiency rows matched in wkf_df, got {n_cand_matched}")

    df["fva_min_flux"] = fmin_v2
    df["fva_max_flux"] = fmax_v2
    df["fva_min_flux_v2"] = fmin_v2
    df["fva_max_flux_v2"] = fmax_v2
    df["fva_measurement_source_v2"] = fva_src
    df["fva_interval_class_v2"] = fva_icls

    df.to_csv(out_path, sep="\t", index=False, compression="gzip")

    w_df = df[is_w]
    k_df = df[is_k]
    return {
        "total_rows": len(df),
        "withdrawal_rows": len(w_df),
        "knockout_rows": len(k_df),
        "withdrawal_corrected_classification_counts": {
            str(k): int(v) for k, v in w_df["corrected_classification"].value_counts().items()
        },
        "withdrawal_addback_disposition_counts": {
            str(k): int(v) for k, v in w_df["addback_rescue_disposition"].value_counts().items()
        },
        "knockout_corrected_classification_counts": {
            str(k): int(v) for k, v in k_df["corrected_classification"].value_counts().items()
        },
        "knockout_fva_measurement_source_counts": {
            str(k): int(v) for k, v in k_df["fva_measurement_source_v2"].value_counts().items()
        },
        "knockout_fva_interval_class_counts": {
            str(k): int(v) for k, v in k_df["fva_interval_class_v2"].value_counts().items()
        },
        "knockout_non_nan_perturbed_growth_primary_count": int(k_df["perturbed_growth_primary"].notna().sum()),
        "knockout_non_nan_addback_rescue_growth_count": int(k_df["addback_rescue_growth"].notna().sum()),
        "knockout_non_null_addback_rescues_growth_count": int(k_df["addback_rescues_growth"].notna().sum()),
    }


def build_candidate_fva_summaries(
    cand_inst_df: pd.DataFrame,
    results_dir: Path,
    rxn_63_out_path: Path,
    strata_244_out_path: Path,
) -> dict[str, Any]:
    """Build both:
    - `exp1_targeted_candidate_fva_63_reactions_summary.tsv` (63 distinct candidate reactions)
    - `exp1_targeted_candidate_fva_244_strata_summary.tsv` (244 `(scenario_id, rxn_id)` summary strata)
    """
    rxn_rows: list[dict[str, Any]] = []
    for (rid, bid, rname, is_tr, tier), grp in cand_inst_df.groupby(
        ["rxn_id", "base_id", "rxn_name", "is_transport", "followup_selection_tier"], sort=True
    ):
        n_tot = len(grp)
        n_orig = int((grp["was_in_original_top30_fva"] == True).sum())
        n_new = int((grp["was_in_original_top30_fva"] == False).sum())
        n_opt = int((grp["fva_both_optimal_and_finite"] == True).sum())
        n_unk = int((grp["fva_both_optimal_and_finite"] == False).sum())
        n_str_1e6 = int((grp["fva_strictly_nonzero_1e-6"] == True).sum())
        n_str_1e12 = int((grp["fva_strictly_nonzero_1e-12"] == True).sum())
        n_trace = int((grp["fva_interval_class"] == "strictly_positive_trace_1e-12_to_1e-6").sum())
        n_spans0 = int((grp["fva_interval_class"] == "spans_zero_contains_zero_flux").sum())
        n_uni_cyc_pass = int(
            (grp["union_network_qc_class"] == "mathematical_union_passed_declared_numerical_and_closed_uptake_egc_checks").sum()
        )
        n_uni_cyc_fail = int(
            (grp["union_network_qc_class"] == "mathematical_union_failed_closed_uptake_atp_cycle_check").sum()
        )

        rxn_rows.append(
            {
                "rxn_id": str(rid),
                "base_id": str(bid),
                "rxn_name": str(rname),
                "is_transport": bool(is_tr),
                "followup_selection_tier": str(tier),
                "n_candidate_instances_total": n_tot,
                "n_scenarios_present": int(grp["scenario_id"].nunique()),
                "n_hq_eq_genomes_present": int(grp["genome"].nunique()),
                "min_alphabetical_rank_in_primary_added": int(grp["alphabetical_rank_in_primary_added"].min()),
                "median_alphabetical_rank_in_primary_added": float(grp["alphabetical_rank_in_primary_added"].median()),
                "max_alphabetical_rank_in_primary_added": int(grp["alphabetical_rank_in_primary_added"].max()),
                "n_assayed_in_original_top30": n_orig,
                "n_newly_evaluated_in_closure": n_new,
                "n_fva_both_optimal_and_ordered": n_opt,
                "n_fva_unknown_or_failed": n_unk,
                "pct_fva_optimal_after_closure": round(100.0 * n_opt / n_tot, 4),
                "n_fva_strictly_nonzero_1e-6": n_str_1e6,
                "n_fva_strictly_nonzero_1e-12": n_str_1e12,
                "n_fva_trace_positive_1e-12_to_1e-6": n_trace,
                "n_fva_spans_zero": n_spans0,
                "n_instances_in_egc_passing_mathematical_union": n_uni_cyc_pass,
                "n_instances_in_egc_failing_mathematical_union": n_uni_cyc_fail,
                "n_instances_in_egc_unknown_mathematical_union": n_tot - n_uni_cyc_pass - n_uni_cyc_fail,
                "median_fva_min_flux": float(grp["fva_min_flux_at_95pct_opt"].median()),
                "median_fva_max_flux": float(grp["fva_max_flux_at_95pct_opt"].median()),
                "median_fva_span": float(grp["fva_span"].median()),
            }
        )
    rxn_63_df = pd.DataFrame(rxn_rows).sort_values("rxn_id").reset_index(drop=True)
    rxn_63_df.to_csv(rxn_63_out_path, sep="\t", index=False)

    orig_sum_df = pd.read_csv(
        results_dir / "exp1_stage4_supported_candidate_deficiencies_summary.tsv", sep="\t"
    )
    strata_rows: list[dict[str, Any]] = []
    inst_by_stratum = {
        (str(s), str(r)): grp
        for (s, r), grp in cand_inst_df.groupby(["scenario_id", "rxn_id"], sort=False)
    }
    for _, srow in orig_sum_df.iterrows():
        key = (str(srow["scenario_id"]), str(srow["rxn_id"]))
        grp = inst_by_stratum[key]
        n_tot = len(grp)
        n_orig = int((grp["was_in_original_top30_fva"] == True).sum())
        n_opt = int((grp["fva_both_optimal_and_finite"] == True).sum())
        n_unk = int((grp["fva_both_optimal_and_finite"] == False).sum())
        n_str_1e6 = int((grp["fva_strictly_nonzero_1e-6"] == True).sum())
        n_str_1e12 = int((grp["fva_strictly_nonzero_1e-12"] == True).sum())
        n_trace = int((grp["fva_interval_class"] == "strictly_positive_trace_1e-12_to_1e-6").sum())
        n_spans0 = int((grp["fva_interval_class"] == "spans_zero_contains_zero_flux").sum())
        n_uni_cyc_pass = int(
            (grp["union_network_qc_class"] == "mathematical_union_passed_declared_numerical_and_closed_uptake_egc_checks").sum()
        )
        n_uni_cyc_fail = int(
            (grp["union_network_qc_class"] == "mathematical_union_failed_closed_uptake_atp_cycle_check").sum()
        )

        d = srow.to_dict()
        d.update(
            {
                "n_instances_in_stratum": n_tot,
                "n_assayed_in_original_top30": n_orig,
                "n_newly_evaluated_in_closure": n_tot - n_orig,
                "n_fva_both_optimal_and_ordered": n_opt,
                "n_fva_unknown_or_failed": n_unk,
                "n_fva_strictly_nonzero_1e-6": n_str_1e6,
                "n_fva_strictly_nonzero_1e-12": n_str_1e12,
                "n_fva_trace_positive_1e-12_to_1e-6": n_trace,
                "n_fva_spans_zero": n_spans0,
                "n_instances_in_egc_passing_mathematical_union": n_uni_cyc_pass,
                "n_instances_in_egc_failing_mathematical_union": n_uni_cyc_fail,
                "n_instances_in_egc_unknown_mathematical_union": n_tot - n_uni_cyc_pass - n_uni_cyc_fail,
                "median_fva_min_flux": float(grp["fva_min_flux_at_95pct_opt"].median()),
                "median_fva_max_flux": float(grp["fva_max_flux_at_95pct_opt"].median()),
                "median_fva_span": float(grp["fva_span"].median()),
            }
        )
        strata_rows.append(d)

    strata_244_df = pd.DataFrame(strata_rows)
    strata_244_df.to_csv(strata_244_out_path, sep="\t", index=False)

    return {
        "n_distinct_reactions": len(rxn_63_df),
        "n_reactions_all_targets_accounted": int((rxn_63_df["n_candidate_instances_total"] > 0).sum()),
        "n_reactions_100pct_optimal_fva": int((rxn_63_df["pct_fva_optimal_after_closure"] == 100.0).sum()),
        "n_strata": len(strata_244_df),
    }


def run_collect_mode(args: argparse.Namespace) -> None:
    t0 = time.time()
    out_dir = args.out_dir
    shard_dir = out_dir / "shards"
    out_dir.mkdir(parents=True, exist_ok=True)

    code_hashes = verify_code_v2_integrity(args.code_v2_dir)

    case_dfs: list[pd.DataFrame] = []
    inst_dfs: list[pd.DataFrame] = []
    shard_metas: list[dict[str, Any]] = []

    for tid in range(args.num_tasks):
        tag = f"fva_closure_shard_{tid:02d}"
        c_path = shard_dir / f"{tag}_cases.tsv"
        i_path = shard_dir / f"{tag}_instances.tsv.gz"
        m_path = shard_dir / f"{tag}_meta.json"
        if not c_path.is_file() or not i_path.is_file() or not m_path.is_file():
            raise FileNotFoundError(f"Missing shard output for task {tid}: {c_path}, {i_path}, or {m_path}")
        case_dfs.append(pd.read_csv(c_path, sep="\t"))
        inst_dfs.append(pd.read_csv(i_path, sep="\t"))
        shard_metas.append(json.loads(m_path.read_text(encoding="utf-8")))

    cases_306_df = pd.concat(case_dfs, ignore_index=True).sort_values(["genome", "scenario_id"]).reset_index(drop=True)
    inst_13750_df = (
        pd.concat(inst_dfs, ignore_index=True)
        .sort_values(["genome", "scenario_id", "rxn_id"])
        .reset_index(drop=True)
    )

    # Prelaunch Review Item 1 & 5: require 100% target accounting (306 cases, 13,750 instances),
    # do NOT crash collector or force an artificial 95%/100% LP success gate.
    if len(cases_306_df) != 306:
        raise RuntimeError(f"Expected 306 accounted (genome, scenario_id, B2) cases, got {len(cases_306_df)}")
    if len(inst_13750_df) != 13750:
        raise RuntimeError(f"Expected 13,750 accounted candidate deficiency instances, got {len(inst_13750_df)}")

    cases_out_path = out_dir / "exp1_targeted_candidate_fva_306_cases_union_qc.tsv"
    inst_out_path = out_dir / "exp1_targeted_candidate_fva_13750_instances.tsv.gz"
    rxn_63_out_path = out_dir / "exp1_targeted_candidate_fva_63_reactions_summary.tsv"
    strata_244_out_path = out_dir / "exp1_targeted_candidate_fva_244_strata_summary.tsv"
    ledger_3872_out_path = out_dir / "exp1_counterfactual_complete_3872_case_ledger.tsv"
    wkf_v2_out_path = out_dir / "exp1_full484_counterfactual_wkf_corrected_v2.tsv.gz"

    cases_306_df.to_csv(cases_out_path, sep="\t", index=False)
    inst_13750_df.to_csv(inst_out_path, sep="\t", index=False, compression="gzip")

    cand_sum_stats = build_candidate_fva_summaries(
        cand_inst_df=inst_13750_df,
        results_dir=args.results_dir,
        rxn_63_out_path=rxn_63_out_path,
        strata_244_out_path=strata_244_out_path,
    )

    wkf_raw_path = args.results_dir / "exp1_full484_withdrawal_knockout_fva.tsv.gz"
    wkf_raw_sha256 = compute_file_sha256(wkf_raw_path)
    wkf_df = pd.read_csv(wkf_raw_path, sep="\t", low_memory=False)

    ledger_stats = build_complete_3872_case_ledger(
        results_dir=args.results_dir,
        wkf_df=wkf_df,
        out_path=ledger_3872_out_path,
    )

    wkf_v2_stats = build_corrected_wkf_v2_table(
        wkf_df=wkf_df,
        cand_inst_df=inst_13750_df,
        out_path=wkf_v2_out_path,
    )

    ct_fva_by_union_qc = (
        pd.crosstab(
            inst_13750_df["union_network_qc_class"],
            inst_13750_df["fva_interval_class"],
        )
        .to_dict(orient="index")
    )
    ct_fva_by_tier = (
        pd.crosstab(
            inst_13750_df["followup_selection_tier"],
            inst_13750_df["fva_interval_class"],
        )
        .to_dict(orient="index")
    )

    summary_payload = {
        "closure_timestamp_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "elapsed_s": round(time.time() - t0, 3),
        "code_v2_sha256": code_hashes,
        "immutable_raw_wkf_sha256": wkf_raw_sha256,
        "complete_3872_case_ledger": ledger_stats,
        "corrected_wkf_v2_table": wkf_v2_stats,
        "targeted_candidate_fva_closure": {
            "n_unique_b2_cases_accounted": len(cases_306_df),
            "n_primary_reconstructions_verified": int(cases_306_df["primary_reconstruction_verified"].sum()),
            "n_primary_reconstructions_unverified": int((~cases_306_df["primary_reconstruction_verified"]).sum()),
            "primary_verification_status_counts": {
                str(k): int(v) for k, v in cases_306_df["primary_verification_status"].value_counts().items()
            },
            "max_primary_g0_abs_diff_vs_ledger": float(cases_306_df["primary_g0_abs_diff_vs_ledger"].max()),
            "max_primary_g0_abs_diff_vs_wkf": float(cases_306_df["primary_g0_abs_diff_vs_wkf"].max()),
            "union_network_qc_class_counts_306_cases": {
                str(k): int(v) for k, v in cases_306_df["union_network_qc_class"].value_counts().items()
            },
            "union_egc_status_counts_306_cases": {
                str(k): int(v) for k, v in cases_306_df["union_egc_check_status"].value_counts().items()
            },
            "n_candidate_instances_accounted": len(inst_13750_df),
            "n_previously_assayed_in_top30": int((inst_13750_df["was_in_original_top30_fva"] == True).sum()),
            "n_newly_evaluated_beyond_top30": int((inst_13750_df["was_in_original_top30_fva"] == False).sum()),
            "n_fva_both_optimal_and_ordered": int((inst_13750_df["fva_both_optimal_and_finite"] == True).sum()),
            "n_fva_failed_or_unknown": int((inst_13750_df["fva_both_optimal_and_finite"] == False).sum()),
            "max_fva_min_residual": float(inst_13750_df["fva_min_max_residual"].max()),
            "max_fva_max_residual": float(inst_13750_df["fva_max_max_residual"].max()),
            "max_fva_min_bound_viol": float(inst_13750_df["fva_min_max_bound_viol"].max()),
            "max_fva_max_bound_viol": float(inst_13750_df["fva_max_max_bound_viol"].max()),
            "max_top30_abs_diff_vs_orig_wkf": float(cases_306_df["max_top30_abs_diff_vs_orig_wkf"].max()),
            "fva_interval_class_counts_13750": {
                str(k): int(v) for k, v in inst_13750_df["fva_interval_class"].value_counts().items()
            },
            "fva_strictly_nonzero_1e-6_counts": {
                str(k): int(v) for k, v in inst_13750_df["fva_strictly_nonzero_1e-6"].value_counts(dropna=False).items()
            },
            "fva_strictly_nonzero_1e-12_counts": {
                str(k): int(v) for k, v in inst_13750_df["fva_strictly_nonzero_1e-12"].value_counts(dropna=False).items()
            },
            "crosstab_fva_interval_by_union_qc_class": ct_fva_by_union_qc,
            "crosstab_fva_interval_by_followup_tier": ct_fva_by_tier,
            **cand_sum_stats,
        },
        "output_artifacts_sha256": {
            cases_out_path.name: compute_file_sha256(cases_out_path),
            inst_out_path.name: compute_file_sha256(inst_out_path),
            rxn_63_out_path.name: compute_file_sha256(rxn_63_out_path),
            strata_244_out_path.name: compute_file_sha256(strata_244_out_path),
            ledger_3872_out_path.name: compute_file_sha256(ledger_3872_out_path),
            wkf_v2_out_path.name: compute_file_sha256(wkf_v2_out_path),
        },
        "shard_metadata": shard_metas,
    }

    summary_json_path = out_dir / "exp1_counterfactual_closure_summary.json"
    summary_json_path.write_text(json.dumps(summary_payload, indent=2) + "\n", encoding="utf-8")
    print(f"[Collect Complete in {summary_payload['elapsed_s']}s] Wrote {summary_json_path}", flush=True)


def main() -> None:
    ap = argparse.ArgumentParser(description="Experiment 1 Counterfactual Closure & Targeted Candidate FVA.")
    ap.add_argument("--mode", choices=("pilot", "shard", "collect"), required=True)
    ap.add_argument(
        "--results-dir",
        type=Path,
        default=Path("/scratch/project/k10018/eq-mechanism-20261109/exp1_metabolic/results"),
    )
    ap.add_argument(
        "--code-v2-dir",
        type=Path,
        default=Path("/scratch/project/k10018/eq-mechanism-20261109/exp1_metabolic/code_v2"),
    )
    ap.add_argument(
        "--gapseq-out-dir",
        type=Path,
        default=Path("/project/k10018/work/eq-annot-20261030/results/c11_gapseq/out"),
    )
    ap.add_argument(
        "--supp-models-dir",
        type=Path,
        default=Path("/scratch/project/k10018/eq-archaeal-function-20261106/results/supplementary_models"),
    )
    ap.add_argument(
        "--out-dir",
        type=Path,
        default=Path("/scratch/project/k10018/eq-mechanism-20261109/exp1_metabolic/results/counterfactual_closure"),
    )
    ap.add_argument("--task-id", type=int, default=0)
    ap.add_argument("--num-tasks", type=int, default=8)
    ap.add_argument("--workers", type=int, default=16)
    args = ap.parse_args()

    if args.mode == "pilot":
        run_pilot_mode(args)
    elif args.mode == "shard":
        run_shard_mode(args)
    elif args.mode == "collect":
        run_collect_mode(args)


if __name__ == "__main__":
    main()
