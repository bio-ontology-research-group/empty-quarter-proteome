#!/usr/bin/env python3
"""Stage 2 of `exp1_metabolic`: Pilot24 exact binary MILP alternative-gap ensemble (`>=20` solutions target),
oxic/anoxic nutrient scenarios (`S1..S6` + native/desert), 5 Methanosarcina biomass perturbations (`B0..B4`),
withdrawal/add-back/KO/FVA, GPR-aware artificial genome-incompleteness masking controls (10%, 20%, 30%),
and strict Pilot24 Scientific QC Gate.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
import pickle
import sys
import tarfile
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

CODE_DIR = Path(__file__).resolve().parent
if str(CODE_DIR) not in sys.path:
    sys.path.insert(0, str(CODE_DIR))

from exp1_sbml_highs_engine import (
    BIO_B0_BASELINE,
    BIO_B1_OMIT_8,
    BIO_B2_OMIT_7_KEEP_F420,
    BIO_B3_OMIT_LIPIDS,
    BIO_B4_CORE_SKELETON,
    BIOMASS_PERTURBATION_RECIPES,
    DEFINED_20AA_AND_PRECURSORS,
    NUTRIENT_SCENARIOS,
    OBSERVED_COHORT_538_POOL_SCOPE,
    SCEN_S0_NATIVE,
    SCEN_S1_ANOXIC_MINIMAL,
    SCEN_S3_ANOXIC_DESERT,
    SCEN_S4_OXIC_DESERT,
    SCEN_S5_ANOXIC_DEFINED,
    SCEN_S6_OXIC_DEFINED,
    ReactionRecord,
    SBMLModelData,
    apply_artificial_gene_masking,
    classify_biomass_requirement_provenance,
    enumerate_alternative_gapfill_ensemble,
    get_planned_protocol_scenario_biomass_tuples,
    parse_gapseq_medium_csv,
    parse_gapseq_sbml_bytes,
    parse_genome_reaction_and_transporter_evidence,
    run_withdrawal_knockout_fva,
    solve_nutrient_supplement_rescue_milp,
)


def load_genome_archive_payload(
    genome: str,
    gapseq_out_dir: Path,
    supp_models_dir: Path,
) -> dict[str, Any]:
    """Load `<g>-draft.xml`, `<g>-all-Reactions.tbl`, `<g>-Transporter.tbl`, and medium CSVs for `genome`."""
    tar_path = gapseq_out_dir / f"{genome}.tar.gz"
    supp_dir = supp_models_dir / genome

    draft_bytes: bytes | None = None
    rxn_tbl_text = ""
    trans_tbl_text = ""
    desert_med_text = ""
    native_med_text = ""

    if tar_path.exists():
        with tarfile.open(tar_path, "r:gz") as tf:
            names = set(tf.getnames())
            d_key = f"{genome}/{genome}-draft.xml"
            r_key = f"{genome}/{genome}-all-Reactions.tbl"
            t_key = f"{genome}/{genome}-Transporter.tbl"
            dm_key = f"{genome}/{genome}-desert_medium.csv"
            nm_key = f"{genome}/{genome}-native_medium.csv"
            if d_key in names:
                draft_bytes = tf.extractfile(d_key).read()  # type: ignore[union-attr]
            if r_key in names:
                rxn_tbl_text = tf.extractfile(r_key).read().decode("utf-8", errors="replace")  # type: ignore[union-attr]
            if t_key in names:
                trans_tbl_text = tf.extractfile(t_key).read().decode("utf-8", errors="replace")  # type: ignore[union-attr]
            if dm_key in names:
                desert_med_text = tf.extractfile(dm_key).read().decode("utf-8", errors="replace")  # type: ignore[union-attr]
            if nm_key in names:
                native_med_text = tf.extractfile(nm_key).read().decode("utf-8", errors="replace")  # type: ignore[union-attr]
    elif supp_dir.exists():
        d_file = supp_dir / f"{genome}-draft.xml"
        r_file = supp_dir / f"{genome}-all-Reactions.tbl"
        t_file = supp_dir / f"{genome}-Transporter.tbl"
        dm_file = supp_dir / f"{genome}-desert_medium.csv"
        nm_file = supp_dir / f"{genome}-native_medium.csv"
        if d_file.exists():
            draft_bytes = d_file.read_bytes()
        if r_file.exists():
            rxn_tbl_text = r_file.read_text(encoding="utf-8", errors="replace")
        if t_file.exists():
            trans_tbl_text = t_file.read_text(encoding="utf-8", errors="replace")
        if dm_file.exists():
            desert_med_text = dm_file.read_text(encoding="utf-8", errors="replace")
        if nm_file.exists():
            native_med_text = nm_file.read_text(encoding="utf-8", errors="replace")

    if draft_bytes is None:
        raise FileNotFoundError(f"Draft SBML not found for {genome}")

    draft_mod = parse_gapseq_sbml_bytes(draft_bytes, provenance_default="draft_sequence_supported")
    rxn_evi, trans_evi = parse_genome_reaction_and_transporter_evidence(rxn_tbl_text, trans_tbl_text)
    desert_med = parse_gapseq_medium_csv(desert_med_text) if desert_med_text.strip() else {}
    native_med = parse_gapseq_medium_csv(native_med_text) if native_med_text.strip() else {}

    return {
        "genome": genome,
        "draft_model": draft_mod,
        "rxn_evi": rxn_evi,
        "trans_evi": trans_evi,
        "desert_med": desert_med,
        "native_med": native_med,
    }


def build_genome_candidate_weights(
    universal_rxns: dict[str, ReactionRecord],
    rxn_evi: dict[str, dict[str, Any]],
    trans_evi: dict[str, dict[str, Any]],
) -> dict[str, float]:
    """Assign gapseq-compatible reaction addition weights (`[0.005, 100.0]`) to each reaction in `universal_rxns`
    using the genome's own sequence alignment evidence (`rxn_evi` and `trans_evi`).
    """
    weights: dict[str, float] = {}
    for rid, rrec in universal_rxns.items():
        bid = rrec.base_id
        m = bid.split("_")[0] if bid.startswith("rxn") else bid
        if m in rxn_evi:
            weights[rid] = float(rxn_evi[m]["weight"])
        elif m in trans_evi:
            weights[rid] = float(trans_evi[m]["weight"])
        elif rrec.is_transport:
            weights[rid] = 25.0
        else:
            weights[rid] = 100.0
    return weights


def evaluate_single_genome_full_protocol(task: dict[str, Any]) -> dict[str, Any]:
    """Run the complete Experiment 1 MILP ensemble, biomass sensitivity, withdrawal/KO/FVA, and supplement
    rescue protocol on a single genome. Shared between Stage 2 (Pilot24) and Stage 3 (Full 484).
    """
    g = task["genome"]
    cohort = task["cohort"]
    canon_order = task["canonical_order"]
    canon_family = task["canonical_family"]
    canon_genus = task["canonical_genus"]
    completeness = float(task["completeness"]) if pd.notna(task["completeness"]) else math.nan
    contamination = float(task["contamination"]) if pd.notna(task["contamination"]) else math.nan

    gapseq_out_dir = Path(task["gapseq_out_dir"])
    supp_models_dir = Path(task["supp_models_dir"])
    pool_pkl_path = Path(task["pool_pkl_path"])
    n_target_solutions = int(task.get("n_target_solutions", 20))
    milp_time_limit_s = float(task.get("milp_time_limit_s", 6.0))

    with open(pool_pkl_path, "rb") as fh:
        pool_data = pickle.load(fh)
    universal_rxns: dict[str, ReactionRecord] = pool_data["reactions"]
    universal_species = pool_data["species"]

    payload = load_genome_archive_payload(g, gapseq_out_dir, supp_models_dir)
    draft_mod: SBMLModelData = payload["draft_model"]
    # Enrich species formula/charge metadata from universal_species where missing
    for sid, sp in draft_mod.species.items():
        if sid in universal_species:
            usp = universal_species[sid]
            if sp.formula is None and usp.formula is not None:
                sp.formula = usp.formula
            if sp.charge is None and usp.charge is not None:
                sp.charge = usp.charge

    rxn_evi = payload["rxn_evi"]
    trans_evi = payload["trans_evi"]
    cand_weights = build_genome_candidate_weights(universal_rxns, rxn_evi, trans_evi)

    summary_rows: list[dict[str, Any]] = []
    solution_ledger_rows: list[dict[str, Any]] = []
    rxn_freq_rows: list[dict[str, Any]] = []
    wkf_rows: list[dict[str, Any]] = []
    supplement_rows: list[dict[str, Any]] = []

    # 1. Run all 5 biomass perturbations (B0..B4) on S3_anoxic_desert_oa_aa and S4_oxic_desert_oa_aa
    #    plus S1..S6 under B0 and B2
    scenarios_to_run: dict[str, dict[str, Any]] = dict(NUTRIENT_SCENARIOS)
    if payload.get("desert_med"):
        d_med = payload["desert_med"]
        # Anoxic desert medium: gapseq desert medium without O2
        scenarios_to_run[SCEN_S3_ANOXIC_DESERT] = {
            "description": "Per-genome gapseq desert medium (anoxic)",
            "redox_regime": "anoxic",
            "compounds": {cid: (cname, u) for cid, (cname, u) in d_med.items() if cid != "cpd00007"},
        }
        # Oxic desert medium: gapseq desert medium with O2
        scenarios_to_run[SCEN_S4_OXIC_DESERT] = {
            "description": "Per-genome gapseq desert medium (oxic with O2)",
            "redox_regime": "oxic",
            "compounds": {**d_med, "cpd00007": ("O2", 10.0)},
        }
    if payload.get("native_med"):
        scenarios_to_run[SCEN_S0_NATIVE] = {
            "description": "Per-genome gapseq predicted native medium",
            "redox_regime": "oxic" if "cpd00007" in payload["native_med"] else "anoxic",
            "compounds": payload["native_med"],
        }

    ens_cache: dict[tuple[str, str], dict[str, Any]] = {}

    for scen_id, scen_info in sorted(scenarios_to_run.items()):
        scen_cpds = scen_info["compounds"]
        redox = scen_info["redox_regime"]
        if scen_id in {SCEN_S3_ANOXIC_DESERT, SCEN_S4_OXIC_DESERT, SCEN_S5_ANOXIC_DEFINED, SCEN_S6_OXIC_DEFINED}:
            bio_list = list(BIOMASS_PERTURBATION_RECIPES.keys())
        else:
            bio_list = [
                BIO_B0_BASELINE,
                BIO_B2_OMIT_7_KEEP_F420,
            ]

        for bio_id in bio_list:
            n_sol_cap = n_target_solutions if bio_id in {
                BIO_B0_BASELINE,
                BIO_B2_OMIT_7_KEEP_F420,
            } else max(5, n_target_solutions // 2)

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
            ens_cache[(scen_id, bio_id)] = ens

            sols = ens["solutions"]
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

            summary_rows.append(
                {
                    "genome": g,
                    "cohort": cohort,
                    "canonical_order": canon_order,
                    "canonical_family": canon_family,
                    "canonical_genus": canon_genus,
                    "completeness": completeness,
                    "contamination": contamination,
                    "scenario_id": scen_id,
                    "redox_regime": redox,
                    "biomass_recipe_id": bio_id,
                    "solver_status": ens["solver_status"],
                    "optimizer_status": ens["optimizer_status"],
                    "raw_solver_status": ens["raw_solver_status"],
                    "enumeration_completeness": ens["enumeration_completeness"],
                    "pool_exhaustion_type": ens["pool_exhaustion_type"],
                    "cost_window_anchor_type": ens["cost_window_anchor_type"],
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
                }
            )

            # Build per-solution ledger rows for this scenario evaluation (per EXP1-FULL-RELEASE-DECISION.txt)
            sol_list = ens.get("solutions", [])
            opt_sols = ens.get("optimization_feasible_solutions", [])
            solver_st = ens["solver_status"]
            is_reval_feasible = solver_st in {
                "optimal",
                "optimal_after_cycle_qc_filtering",
                "optimal_draft_no_gapfill_needed",
                "heuristic_incumbent_bounded",
                "post_pruning_cost_changed_not_proven_optimal",
            }
            if sol_list:
                for s in sol_list:
                    r_ids = sorted(s.get("added_rxns", []))
                    r_str = ";".join(r_ids)
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
                    solution_ledger_rows.append(
                        {
                            "genome": g,
                            "cohort": cohort,
                            "canonical_order": canon_order,
                            "canonical_family": canon_family,
                            "scenario_id": scen_id,
                            "biomass_recipe_id": bio_id,
                            "solution_index": int(s.get("solution_index", 0)),
                            "validity_disposition": disposition,
                            "solver_status": solver_st,
                            "selected_reaction_ids": r_str,
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
                        }
                    )
            elif opt_sols:
                for s in opt_sols:
                    r_ids = sorted(s.get("added_rxns", []))
                    res_val = float(s.get("fba_check_max_residual", math.nan))
                    bnd_val = float(s.get("fba_check_max_bound_viol", math.nan))
                    gr_val = float(s.get("verified_fba_growth", math.nan))
                    cyc_st = str(s.get("cycle_qc_status", "rejected"))
                    spur_atp = s.get("post_pruning_egc_has_spurious_atp_cycle", True)
                    solution_ledger_rows.append(
                        {
                            "genome": g,
                            "cohort": cohort,
                            "canonical_order": canon_order,
                            "canonical_family": canon_family,
                            "scenario_id": scen_id,
                            "biomass_recipe_id": bio_id,
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
                        }
                    )
            else:
                solution_ledger_rows.append(
                    {
                        "genome": g,
                        "cohort": cohort,
                        "canonical_order": canon_order,
                        "canonical_family": canon_family,
                        "scenario_id": scen_id,
                        "biomass_recipe_id": bio_id,
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
                    }
                )

            # Run withdrawal / add-back rescue / ensemble-union KO / FVA on key desert & rich scenarios
            if scen_id in {SCEN_S3_ANOXIC_DESERT, SCEN_S4_OXIC_DESERT, SCEN_S5_ANOXIC_DEFINED, SCEN_S6_OXIC_DEFINED} and bio_id in {
                BIO_B0_BASELINE,
                BIO_B2_OMIT_7_KEEP_F420,
            }:
                wkf = run_withdrawal_knockout_fva(
                    base_model=draft_mod,
                    ensemble_result=ens,
                    candidate_reactions=universal_rxns,
                    scenario_compounds=scen_cpds,
                    biomass_recipe_id=bio_id,
                    min_growth=0.01,
                    time_limit_s=5.0,
                )
                for w in wkf["withdrawals"]:
                    wkf_rows.append(
                        {
                            "genome": g,
                            "cohort": cohort,
                            "canonical_family": canon_family,
                            "scenario_id": scen_id,
                            "biomass_recipe_id": bio_id,
                            "evaluation_status": wkf["evaluation_status"],
                            "record_type": "nutrient_withdrawal_addback",
                            "target_id": w["compound_id"],
                            "target_name": w["compound_name"],
                            "baseline_growth": w["baseline_growth"],
                            "perturbed_growth_primary": w["withdrawn_growth_primary"],
                            "perturbed_growth_ensemble_union": w["withdrawn_growth_ensemble_union"],
                            "addback_rescue_growth": w["addback_rescue_growth"],
                            "addback_rescues_growth": w["addback_rescues_growth"],
                            "classification": w["dependency_classification"],
                            "fva_min_flux": math.nan,
                            "fva_max_flux": math.nan,
                        }
                    )
                fva_map = {f["rxn_id"]: f for f in wkf["fva"]}
                for ko in wkf["knockouts"]:
                    f_entry = fva_map.get(ko["rxn_id"], {})
                    wkf_rows.append(
                        {
                            "genome": g,
                            "cohort": cohort,
                            "canonical_family": canon_family,
                            "scenario_id": scen_id,
                            "biomass_recipe_id": bio_id,
                            "evaluation_status": wkf["evaluation_status"],
                            "record_type": "reaction_knockout_and_fva",
                            "target_id": ko["rxn_id"],
                            "target_name": ko["base_id"],
                            "baseline_growth": ko["baseline_growth"],
                            "perturbed_growth_primary": 0.0,
                            "perturbed_growth_ensemble_union": ko["knockout_growth_in_ensemble_union"],
                            "addback_rescue_growth": ko["baseline_growth"],
                            "addback_rescues_growth": True,
                            "classification": ko["knockout_class"],
                            "fva_min_flux": f_entry.get("fva_min_flux_at_95pct_opt", math.nan),
                            "fva_max_flux": f_entry.get("fva_max_flux_at_95pct_opt", math.nan),
                        }
                    )

            # Separate exploratory nutrient-supplement rescue if fixed minimal/desert scenario was infeasible
            if ens["solver_status"] == "infeasible_under_fixed_medium" and scen_id in {SCEN_S1_ANOXIC_MINIMAL, SCEN_S3_ANOXIC_DESERT, SCEN_S4_OXIC_DESERT} and bio_id == BIO_B2_OMIT_7_KEEP_F420:
                supp_cands = {
                    cid: (cname, float(max_u), 50.0)
                    for cid, (cname, max_u) in DEFINED_20AA_AND_PRECURSORS.items()
                    if cid not in scen_cpds
                }
                supp_res = solve_nutrient_supplement_rescue_milp(
                    draft_model=draft_mod,
                    candidate_reactions=universal_rxns,
                    candidate_weights=cand_weights,
                    fixed_scenario_compounds=scen_cpds,
                    supplement_candidate_compounds=supp_cands,
                    biomass_recipe_id=bio_id,
                    redox_regime=redox,
                    min_growth=0.01,
                    time_limit_s=milp_time_limit_s,
                )
                supplement_rows.append(
                    {
                        "genome": g,
                        "cohort": cohort,
                        "scenario_id": scen_id,
                        "redox_regime": redox,
                        "biomass_recipe_id": bio_id,
                        "experiment_type": supp_res["experiment_type"],
                        "rescued": supp_res["rescued"],
                        "solver_status": supp_res["solver_status"],
                        "primary_supplements_required_csv": ",".join(supp_res["primary_supplements_required"]),
                    }
                )

    # 2. Cross-reference B0..B4 biomass perturbations on scenarios with B0..B4 evaluated to classify every reaction requirement's biomass provenance
    focus_scenarios = [
        s for s in (SCEN_S3_ANOXIC_DESERT, SCEN_S4_OXIC_DESERT, SCEN_S5_ANOXIC_DEFINED, SCEN_S6_OXIC_DEFINED)
        if any(bool(ens_cache.get((s, b), {}).get("reaction_frequencies")) for b in BIOMASS_PERTURBATION_RECIPES)
    ]
    if not focus_scenarios:
        focus_scenarios = [SCEN_S3_ANOXIC_DESERT, SCEN_S4_OXIC_DESERT]

    for focus_scen in focus_scenarios:
        b0_freqs = ens_cache.get((focus_scen, BIO_B0_BASELINE), {}).get("reaction_frequencies", {})
        b1_freqs = ens_cache.get((focus_scen, BIO_B1_OMIT_8), {}).get("reaction_frequencies", {})
        b2_freqs = ens_cache.get((focus_scen, BIO_B2_OMIT_7_KEEP_F420), {}).get("reaction_frequencies", {})
        b3_freqs = ens_cache.get((focus_scen, BIO_B3_OMIT_LIPIDS), {}).get("reaction_frequencies", {})
        b4_freqs = ens_cache.get((focus_scen, BIO_B4_CORE_SKELETON), {}).get("reaction_frequencies", {})

        all_rids = sorted(set(b0_freqs) | set(b1_freqs) | set(b2_freqs) | set(b3_freqs) | set(b4_freqs))
        for rid in all_rids:
            in_b0 = rid in b0_freqs
            in_b1 = rid in b1_freqs
            in_b2 = rid in b2_freqs
            in_b3 = rid in b3_freqs
            in_b4 = rid in b4_freqs
            bio_prov_class = classify_biomass_requirement_provenance(
                req_in_b0_baseline=in_b0,
                req_in_b1_omit_all8=in_b1,
                req_in_b2_omit_7_keep_f420=in_b2,
                req_in_b3_omit_surface_lipids=in_b3,
                req_in_b4_core_skeleton=in_b4,
            )
            ref_entry = b0_freqs.get(rid) or b2_freqs.get(rid) or b1_freqs.get(rid) or b3_freqs.get(rid) or b4_freqs[rid]
            bid_clean = ref_entry["base_id"].split("_")[0] if ref_entry["base_id"].startswith("rxn") else ref_entry["base_id"]
            ev = rxn_evi.get(bid_clean, {})
            rxn_freq_rows.append(
                {
                    "genome": g,
                    "cohort": cohort,
                    "canonical_order": canon_order,
                    "canonical_family": canon_family,
                    "canonical_genus": canon_genus,
                    "completeness": completeness,
                    "contamination": contamination,
                    "scenario_id": focus_scen,
                    "rxn_id": rid,
                    "base_id": ref_entry["base_id"],
                    "rxn_name": ref_entry["name"],
                    "is_transport": ref_entry["is_transport"],
                    "pool_source_tier": universal_rxns[rid].provenance if rid in universal_rxns else "unknown",
                    "candidate_pool_scope": OBSERVED_COHORT_538_POOL_SCOPE,
                    "candidate_weight": ref_entry["candidate_weight"],
                    "blast_status_in_genome": ev.get("status", "no_blast"),
                    "blast_bitscore_in_genome": ev.get("bitscore", 0.0),
                    "best_gene_in_genome": ev.get("best_gene", ""),
                    "pathway_id": ev.get("pathway", ""),
                    "freq_B0_template_baseline": b0_freqs[rid]["ensemble_frequency"] if in_b0 else 0.0,
                    "freq_B1_omit_all8_cofactors": b1_freqs[rid]["ensemble_frequency"] if in_b1 else 0.0,
                    "freq_B2_omit_7_keep_f420": b2_freqs[rid]["ensemble_frequency"] if in_b2 else 0.0,
                    "freq_B3_omit_surface_lipids": b3_freqs[rid]["ensemble_frequency"] if in_b3 else 0.0,
                    "freq_B4_core_skeleton": b4_freqs[rid]["ensemble_frequency"] if in_b4 else 0.0,
                    "stability_tier_B0": b0_freqs[rid]["stability_tier"] if in_b0 else "absent",
                    "stability_tier_B2": b2_freqs[rid]["stability_tier"] if in_b2 else "absent",
                    "supernetwork_ko_necessity_B0": b0_freqs[rid]["supernetwork_ko_necessity"] if in_b0 else "not_in_B0_ensemble",
                    "supernetwork_ko_necessity_B2": b2_freqs[rid]["supernetwork_ko_necessity"] if in_b2 else "not_in_B2_ensemble",
                    "biomass_sensitivity_classification": bio_prov_class,
                    "enumeration_completeness_B0": b0_freqs[rid]["enumeration_completeness"] if in_b0 else "absent",
                    "enumeration_completeness_B2": b2_freqs[rid]["enumeration_completeness"] if in_b2 else "absent",
                }
            )

    # 3. Run artificial genome-incompleteness masking controls (10%, 20%, 30%) if requested for this HQ genome
    masking_rows: list[dict[str, Any]] = []
    if task.get("run_masking_control", False):
        base_ens_s4 = ens_cache[(SCEN_S4_OXIC_DESERT, BIO_B2_OMIT_7_KEEP_F420)]
        base_rxns_s4 = set(base_ens_s4["reaction_frequencies"].keys())
        s4_cpds = NUTRIENT_SCENARIOS[SCEN_S4_OXIC_DESERT]["compounds"]

        for mask_frac in (0.10, 0.20, 0.30):
            for rep in range(3):
                seed_val = abs(hash((g, int(mask_frac * 100), rep))) % (2**31 - 1)
                rng = np.random.default_rng(seed_val)
                masked_draft, mask_meta = apply_artificial_gene_masking(draft_mod, mask_fraction=mask_frac, rng=rng)
                masked_weights = dict(cand_weights)
                for rem_rid in mask_meta["removed_rxn_ids"]:
                    masked_weights[rem_rid] = 100.0

                m_ens = enumerate_alternative_gapfill_ensemble(
                    draft_model=masked_draft,
                    candidate_reactions=universal_rxns,
                    candidate_weights=masked_weights,
                    scenario_compounds=s4_cpds,
                    biomass_recipe_id=BIO_B2_OMIT_7_KEEP_F420,
                    min_growth=0.01,
                    n_Target_solutions=8,
                    cost_slack_frac=0.25,
                    candidate_pool_scope="artificial_gene_masking_control",
                    time_limit_per_milp_s=milp_time_limit_s,
                )
                m_rxns = set(m_ens["reaction_frequencies"].keys())
                newly_induced = m_rxns - base_rxns_s4
                removed_set = set(mask_meta["removed_rxn_ids"])
                recovered_true_masked = newly_induced & removed_set
                false_bypass_rxns = newly_induced - removed_set

                masking_rows.append(
                    {
                        "genome": g,
                        "cohort": cohort,
                        "completeness": completeness,
                        "contamination": contamination,
                        "mask_fraction": mask_frac,
                        "replicate": rep,
                        "n_total_genes": mask_meta["n_total_genes"],
                        "n_masked_genes": mask_meta["n_masked_genes"],
                        "n_removed_gpr_reactions": mask_meta["n_removed_gpr_reactions"],
                        "masked_solver_status": m_ens["solver_status"],
                        "n_baseline_gapfill_rxns": len(base_rxns_s4),
                        "n_masked_gapfill_rxns": len(m_rxns),
                        "n_newly_induced_gapfill_rxns": len(newly_induced),
                        "n_true_masked_rxns_recovered": len(recovered_true_masked),
                        "n_alternative_bypass_rxns_induced": len(false_bypass_rxns),
                    }
                )

    return {
        "genome": g,
        "summary_rows": summary_rows,
        "solution_ledger_rows": solution_ledger_rows,
        "rxn_freq_rows": rxn_freq_rows,
        "wkf_rows": wkf_rows,
        "supplement_rows": supplement_rows,
        "masking_rows": masking_rows,
    }


REQUIRED_VERIFICATION_COLUMNS: tuple[str, ...] = (
    "n_distinct_valid_solutions",
    "min_added_rxns_count",
    "min_total_cost",
    "min_solution_verified_growth",
    "max_solution_residual",
    "max_solution_bound_viol",
)


def validate_pilot24_scientific_gate(
    summary_df: pd.DataFrame,
    rxn_freq_df: pd.DataFrame,
    wkf_df: pd.DataFrame,
    masking_df: pd.DataFrame,
    expected_n_genomes: int,
    expected_genome_ids: list[str] | set[str] | None = None,
    expected_masking_genome_ids: list[str] | set[str] | None = None,
    solution_ledger_df: pd.DataFrame | None = None,
) -> dict[str, Any]:
    """Enforce the Pilot24 Scientific QC Gate (`EXP1-ACCEPTANCE-AND-COLLECTOR-REVIEW.txt` Item 1,
    `EXP1-FOLLOWUP-AND-GATE-COUNTEREXAMPLES.txt` Item A, and `EXP1-FULL-RELEASE-DECISION.txt`).
    Replaces arbitrary 95% aggregate completion gate with exact execution coverage and correctness:
    - Exact planned genome IDs evaluated (no missing, no unexpected).
    - Exact planned scenario/biomass tuples per genome (exact set comparison against declared manifest protocol).
    - No duplicate evaluations across (genome, scenario_id, biomass_recipe_id).
    - Statuses exhaustive and honest (all belong to declared status sets).
    - Independently validated accepted solutions:
      Every accepted feasible solution MUST have finite verified residual <= 1e-5, bound violation <= 1e-5,
      verified growth >= 0.00999, non-negative cost, valid added reaction count, and distinct solutions > 0.
      Missing columns or NA values MUST fail acceptance, NEVER fillna(0.0).
      Non-feasible/unknown outcomes must have NA added reaction counts and NA costs.
    - Per-solution ledger consistency verified if provided.
    - Masking controls executed and audited on exact eligible HQ reference genomes (7 genomes x 3 fractions x 3 replicates).
    - Differentiates feasible, proven infeasible, timeouts, numerical failures, and cycle failures.
    """
    n_genomes_evaluated = int(summary_df["genome"].nunique()) if not summary_df.empty else 0
    n_summary_rows = len(summary_df)

    # Status classification per EXP1-QC-DESIGN-CORRECTION.txt & EXP1-ACCEPTANCE-AND-COLLECTOR-REVIEW.txt
    revalidated_feasible_statuses = {
        "optimal",
        "optimal_after_cycle_qc_filtering",
        "optimal_draft_no_gapfill_needed",
        "heuristic_incumbent_bounded",
        "post_pruning_cost_changed_not_proven_optimal",
    }
    proven_infeasible_statuses = {
        "infeasible",
        "infeasible_under_fixed_medium",
    }
    timeout_statuses = {
        "time_limit_no_incumbent",
        "time_limit_heuristic_incumbent",
    }
    numerical_failure_statuses = {
        "numerical_coupling_violation",
        "fba_verification_failed",
        "numerical_error",
    }
    known_statuses = revalidated_feasible_statuses | proven_infeasible_statuses | timeout_statuses | numerical_failure_statuses

    status_counts = summary_df["solver_status"].value_counts().to_dict() if not summary_df.empty else {}
    n_revalidated_feasible = sum(status_counts.get(s, 0) for s in revalidated_feasible_statuses)
    n_proven_infeasible = sum(status_counts.get(s, 0) for s in proven_infeasible_statuses)
    n_timeouts = sum(status_counts.get(s, 0) for s in timeout_statuses)
    n_numerical_failures = sum(status_counts.get(s, 0) for s in numerical_failure_statuses)
    n_cycle_failures = int(summary_df["n_cycle_qc_rejected_solutions"].sum()) if not summary_df.empty and "n_cycle_qc_rejected_solutions" in summary_df.columns else 0

    handled_statuses = revalidated_feasible_statuses | proven_infeasible_statuses | timeout_statuses
    n_handled = sum(status_counts.get(s, 0) for s in handled_statuses)
    frac_handled = float(n_handled) / float(max(1, n_summary_rows))
    frac_revalidated_feasible = float(n_revalidated_feasible) / float(max(1, n_summary_rows))
    frac_proven_infeasible = float(n_proven_infeasible) / float(max(1, n_summary_rows))
    frac_timeouts = float(n_timeouts) / float(max(1, n_summary_rows))
    frac_numerical_failures = float(n_numerical_failures) / float(max(1, n_summary_rows))

    n_b0_vs_b2_sensitive = int(
        rxn_freq_df["biomass_sensitivity_classification"].str.startswith("biomass_assumption_sensitive").sum()
    ) if not rxn_freq_df.empty else 0
    n_robust_rxns = int(
        (rxn_freq_df["biomass_sensitivity_classification"] == "robust_across_all_tested_biomass_perturbations").sum()
    ) if not rxn_freq_df.empty else 0

    # Counts by genome and scenario
    counts_by_genome = summary_df.groupby("genome")["solver_status"].value_counts().unstack(fill_value=0).to_dict(orient="index") if not summary_df.empty else {}
    counts_by_scenario = summary_df.groupby("scenario_id")["solver_status"].value_counts().unstack(fill_value=0).to_dict(orient="index") if not summary_df.empty else {}

    # 1. Exact planned genome IDs
    actual_genome_set = set(summary_df["genome"].unique()) if not summary_df.empty else set()
    if expected_genome_ids is not None:
        exact_genome_ids_matched = bool(actual_genome_set == set(expected_genome_ids))
    else:
        exact_genome_ids_matched = bool(n_genomes_evaluated == expected_n_genomes and expected_n_genomes > 0)

    # 2. No duplicate evaluations
    no_duplicate_evaluations = bool(not summary_df.duplicated(subset=["genome", "scenario_id", "biomass_recipe_id"]).any() if not summary_df.empty else False)

    # 3. Exact planned scenario/biomass combinations per genome (Root Item A: exact set comparison against declared manifest protocol)
    target_genomes = actual_genome_set if expected_genome_ids is None else set(expected_genome_ids)
    expected_tuples = get_planned_protocol_scenario_biomass_tuples(target_genomes)
    actual_tuples = set(zip(summary_df["genome"].astype(str), summary_df["scenario_id"].astype(str), summary_df["biomass_recipe_id"].astype(str))) if not summary_df.empty else set()
    exact_tuples_per_genome = bool(actual_tuples == expected_tuples) if expected_tuples else False

    # 4. Statuses exhaustive and honest
    actual_statuses = set(summary_df["solver_status"].unique()) if not summary_df.empty else set()
    statuses_exhaustive_and_honest = bool(actual_statuses.issubset(known_statuses) and len(actual_statuses) > 0)

    # 5. Correctness of accepted solutions vs non-feasible outcomes
    # Fixed acceptance hole per EXP1-FULL-RELEASE-DECISION.txt:
    # All required verification columns MUST exist, and accepted rows MUST have finite verified evidence (no fillna(0.0)).
    required_verification_columns = (
        "n_distinct_valid_solutions",
        "min_added_rxns_count",
        "min_total_cost",
        "min_solution_verified_growth",
        "max_solution_residual",
        "max_solution_bound_viol",
    )
    missing_req_cols = [c for c in required_verification_columns if c not in summary_df.columns]
    has_all_req_cols = len(missing_req_cols) == 0

    accepted_mask = summary_df["solver_status"].isin(revalidated_feasible_statuses)
    accepted_df = summary_df[accepted_mask]
    non_accepted_df = summary_df[~accepted_mask]

    accepted_solutions_valid = True
    if not accepted_df.empty:
        if not has_all_req_cols:
            accepted_solutions_valid = False
        else:
            distinct_ok = (accepted_df["n_distinct_valid_solutions"].notna()) & (accepted_df["n_distinct_valid_solutions"] > 0)
            added_ok = (accepted_df["min_added_rxns_count"].notna()) & (accepted_df["min_added_rxns_count"] >= 0)
            cost_ok = (accepted_df["min_total_cost"].notna()) & (accepted_df["min_total_cost"] >= 0.0)
            growth_ok = (accepted_df["min_solution_verified_growth"].notna()) & np.isfinite(accepted_df["min_solution_verified_growth"]) & (accepted_df["min_solution_verified_growth"] >= 0.00999)
            residual_ok = (accepted_df["max_solution_residual"].notna()) & np.isfinite(accepted_df["max_solution_residual"]) & (accepted_df["max_solution_residual"] >= 0.0) & (accepted_df["max_solution_residual"] <= 1e-5)
            bound_ok = (accepted_df["max_solution_bound_viol"].notna()) & np.isfinite(accepted_df["max_solution_bound_viol"]) & (accepted_df["max_solution_bound_viol"] >= 0.0) & (accepted_df["max_solution_bound_viol"] <= 1e-5)

            if "n_cycle_qc_rejected_solutions" in accepted_df.columns and "n_optimization_feasible_solutions" in accepted_df.columns:
                cycle_ok = ~((accepted_df["n_cycle_qc_rejected_solutions"] > 0) & (accepted_df["n_cycle_qc_rejected_solutions"] == accepted_df["n_optimization_feasible_solutions"]) & (accepted_df["n_distinct_valid_solutions"] == 0))
            else:
                cycle_ok = pd.Series(True, index=accepted_df.index)

            accepted_solutions_valid = bool(
                distinct_ok.all() and added_ok.all() and cost_ok.all() and growth_ok.all() and residual_ok.all() and bound_ok.all() and cycle_ok.all()
            )

    non_accepted_clean = True
    if not non_accepted_df.empty:
        if "n_distinct_valid_solutions" in non_accepted_df.columns:
            non_accepted_clean = non_accepted_clean and bool((non_accepted_df["n_distinct_valid_solutions"].fillna(0) == 0).all())
        if "min_added_rxns_count" in non_accepted_df.columns:
            non_accepted_clean = non_accepted_clean and bool(non_accepted_df["min_added_rxns_count"].isna().all())
        if "min_total_cost" in non_accepted_df.columns:
            non_accepted_clean = non_accepted_clean and bool(non_accepted_df["min_total_cost"].isna().all())

    ledger_valid = True
    if solution_ledger_df is not None:
        if solution_ledger_df.empty:
            ledger_valid = (n_summary_rows == 0)
        else:
            # 1. Exact correspondence of keys (genome, scenario_id, biomass_recipe_id)
            sum_keys = set(
                zip(
                    summary_df["genome"].astype(str),
                    summary_df["scenario_id"].astype(str),
                    summary_df["biomass_recipe_id"].astype(str),
                )
            )
            led_keys = set(
                zip(
                    solution_ledger_df["genome"].astype(str),
                    solution_ledger_df["scenario_id"].astype(str),
                    solution_ledger_df["biomass_recipe_id"].astype(str),
                )
            )
            keys_match = bool(sum_keys == led_keys)

            # 2. Join exact keys against summary counts
            acc_ledger = solution_ledger_df[solution_ledger_df["validity_disposition"] == "accepted_feasible"]
            acc_counts = (
                acc_ledger.groupby(["genome", "scenario_id", "biomass_recipe_id"]).size().to_dict()
                if not acc_ledger.empty
                else {}
            )

            counts_match = True
            for _, srow in summary_df.iterrows():
                k = (str(srow["genome"]), str(srow["scenario_id"]), str(srow["biomass_recipe_id"]))
                exp_count = (
                    int(srow["n_distinct_valid_solutions"])
                    if srow["solver_status"] in revalidated_feasible_statuses and pd.notna(srow["n_distinct_valid_solutions"])
                    else 0
                )
                act_count = acc_counts.get(k, 0)
                if exp_count != act_count:
                    counts_match = False
                    break

            # 3. Strict verification of all accepted rows in ledger
            acc_rows_valid = True
            if not acc_ledger.empty:
                l_res_ok = (
                    acc_ledger["stoichiometric_max_residual"].notna()
                    & np.isfinite(acc_ledger["stoichiometric_max_residual"])
                    & (acc_ledger["stoichiometric_max_residual"] >= 0.0)
                    & (acc_ledger["stoichiometric_max_residual"] <= 1e-5)
                )
                l_bnd_ok = (
                    acc_ledger["max_bound_violation"].notna()
                    & np.isfinite(acc_ledger["max_bound_violation"])
                    & (acc_ledger["max_bound_violation"] >= 0.0)
                    & (acc_ledger["max_bound_violation"] <= 1e-5)
                )
                l_gr_ok = (
                    acc_ledger["verified_growth"].notna()
                    & np.isfinite(acc_ledger["verified_growth"])
                    & (acc_ledger["verified_growth"] >= 0.00999)
                )
                l_ver_ok = acc_ledger["residual_bounds_growth_verified"] == True
                l_cyc_ok = acc_ledger["cycle_qc_status"].isin(
                    {"passed_closed_uptake_egc_check", "passed_no_free_atp_flux", "evaluated_no_atp_species", "not_required_no_atp_species"}
                )
                l_spur_ok = acc_ledger["spurious_atp_cycle"] != True
                acc_rows_valid = bool(
                    l_res_ok.all()
                    and l_bnd_ok.all()
                    and l_gr_ok.all()
                    and l_ver_ok.all()
                    and l_cyc_ok.all()
                    and l_spur_ok.all()
                )
            elif n_revalidated_feasible > 0:
                acc_rows_valid = False

            # 4. Strict check on non-accepted rows in ledger
            non_acc_ledger = solution_ledger_df[solution_ledger_df["validity_disposition"] != "accepted_feasible"]
            non_acc_rows_valid = True
            if not non_acc_ledger.empty:
                non_acc_rows_valid = bool((non_acc_ledger["residual_bounds_growth_verified"] == False).all())

            ledger_valid = bool(keys_match and counts_match and acc_rows_valid and non_acc_rows_valid)

    per_solution_correctness_verified = bool(
        accepted_solutions_valid and non_accepted_clean and ledger_valid and (has_all_req_cols or accepted_df.empty)
    )

    # 6. Masking controls executed and audited on exact eligible HQ reference genomes
    if expected_n_genomes > 1:
        if expected_masking_genome_ids is not None:
            expected_mask_genomes = set(expected_masking_genome_ids)
        else:
            expected_mask_genomes = {g for g in actual_genome_set if str(g).startswith("MGYG")}  # 7 external HQ reference genomes
        actual_mask_genomes = set(masking_df["genome"].unique()) if not masking_df.empty else set()
        
        mask_genomes_ok = bool(actual_mask_genomes == expected_mask_genomes) if expected_mask_genomes else (len(actual_mask_genomes) >= 7)
        mask_statuses_ok = bool(set(masking_df["masked_solver_status"].unique()).issubset(known_statuses)) if not masking_df.empty and "masked_solver_status" in masking_df.columns else True
        mask_no_dups = bool(not masking_df.duplicated(subset=["genome", "mask_fraction", "replicate"]).any()) if not masking_df.empty else False
        expected_min_mask_rows = len(expected_mask_genomes) * 9 if expected_mask_genomes else 63
        masking_controls_audited = bool(mask_genomes_ok and mask_statuses_ok and mask_no_dups and len(masking_df) >= expected_min_mask_rows)
    else:
        masking_controls_audited = not masking_df.empty

    checks = {
        "all_manifest_genomes_evaluated": exact_genome_ids_matched,
        "exact_planned_genome_ids_matched": exact_genome_ids_matched,
        "exact_planned_scenario_biomass_tuples_per_genome": exact_tuples_per_genome,
        "no_duplicate_evaluations": no_duplicate_evaluations,
        "statuses_exhaustive_and_honest": statuses_exhaustive_and_honest,
        "per_solution_correctness_verified": per_solution_correctness_verified,
        "masking_controls_executed_and_audited": masking_controls_audited,
        "solver_clean_completion_rate_ge_95pct": bool(frac_handled >= 0.95),
        "biomass_perturbations_b0_to_b4_populated": int(summary_df["biomass_recipe_id"].nunique()) == 5 if not summary_df.empty else False,
        "nutrient_scenarios_s1_to_s6_populated": int(summary_df["scenario_id"].nunique()) >= 6 if not summary_df.empty else False,
        "reaction_frequencies_and_biomass_sensitivity_non_empty": len(rxn_freq_df) > 0,
        "withdrawal_knockout_fva_non_empty": len(wkf_df) > 0,
        "masking_controls_non_empty": len(masking_df) > 0,
    }
    # Per Root EXP1-ACCEPTANCE-AND-COLLECTOR-REVIEW.txt Item 1: Gate release is determined by exact coverage,
    # exhaustive honest accounting, per-solution correctness verification, and masking audits,
    # NOT by an arbitrary 95% aggregate fraction.
    gating_keys = [
        "exact_planned_genome_ids_matched",
        "exact_planned_scenario_biomass_tuples_per_genome",
        "no_duplicate_evaluations",
        "statuses_exhaustive_and_honest",
        "per_solution_correctness_verified",
        "masking_controls_executed_and_audited",
    ]
    gate_passed = all(checks[k] for k in gating_keys)
    return {
        "gate_passed": gate_passed,
        "n_genomes_evaluated": n_genomes_evaluated,
        "expected_n_genomes": expected_n_genomes,
        "n_summary_rows": n_summary_rows,
        "solver_clean_completion_rate": frac_handled,
        "rate_handled_outcomes": frac_handled,
        "rate_feasible_independently_revalidated": frac_revalidated_feasible,
        "rate_proven_infeasible": frac_proven_infeasible,
        "rate_unknown_timeouts": frac_timeouts,
        "rate_numerical_failures": frac_numerical_failures,
        "n_cycle_qc_rejected_solutions": n_cycle_failures,
        "counts_by_status": status_counts,
        "counts_by_genome": counts_by_genome,
        "counts_by_scenario": counts_by_scenario,
        "n_reaction_stability_rows": len(rxn_freq_df),
        "n_biomass_assumption_sensitive_rxn_instances": n_b0_vs_b2_sensitive,
        "n_robust_rxn_instances": n_robust_rxns,
        "n_withdrawal_knockout_fva_rows": len(wkf_df),
        "n_masking_control_rows": len(masking_df),
        "checks": checks,
    }


def main() -> None:
    ap = argparse.ArgumentParser(description="Stage 2: Pilot24 + HQ8 exact binary MILP ensembles, masking controls, and QC gate.")
    ap.add_argument("--manifest-tsv", type=Path, required=True)
    ap.add_argument("--pool-pkl", type=Path, required=True)
    ap.add_argument("--gapseq-out-dir", type=Path, default=Path("/project/k10018/work/eq-annot-20261030/results/c11_gapseq/out"))
    ap.add_argument("--supp-models-dir", type=Path, default=Path("/scratch/project/k10018/eq-archaeal-function-20261106/results/supplementary_models"))
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--n-target-solutions", type=int, default=20)
    ap.add_argument("--milp-time-limit-s", type=float, default=12.0)
    ap.add_argument("--workers", type=int, default=32)
    args = ap.parse_args()

    t0 = time.time()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    manifest_df = pd.read_csv(args.manifest_tsv, sep="\t")

    tasks: list[dict[str, Any]] = []
    for _, row in manifest_df.iterrows():
        comp_v = float(row["completeness"]) if pd.notna(row["completeness"]) else 0.0
        cont_v = float(row["contamination"]) if pd.notna(row["contamination"]) else 100.0
        is_hq_mask = (str(row.get("pilot_role", "")) == "external_HQ_masking_reference")
        tasks.append(
            {
                "genome": str(row["genome"]),
                "cohort": str(row["cohort"]),
                "canonical_order": str(row["canonical_order"]),
                "canonical_family": str(row["canonical_family"]),
                "canonical_genus": str(row["canonical_genus"]),
                "completeness": row["completeness"],
                "contamination": row["contamination"],
                "gapseq_out_dir": str(args.gapseq_out_dir),
                "supp_models_dir": str(args.supp_models_dir),
                "pool_pkl_path": str(args.pool_pkl),
                "n_target_solutions": args.n_target_solutions,
                "milp_time_limit_s": args.milp_time_limit_s,
                "run_masking_control": is_hq_mask,
            }
        )

    all_summary: list[dict[str, Any]] = []
    all_solution_ledger: list[dict[str, Any]] = []
    all_rxn_freq: list[dict[str, Any]] = []
    all_wkf: list[dict[str, Any]] = []
    all_supp: list[dict[str, Any]] = []
    all_masking: list[dict[str, Any]] = []

    with ProcessPoolExecutor(max_workers=min(args.workers, max(1, len(tasks)), os.cpu_count() or 32)) as ex:
        futs = {ex.submit(evaluate_single_genome_full_protocol, t): t["genome"] for t in tasks}
        for idx, fut in enumerate(as_completed(futs), 1):
            g_name = futs[fut]
            res = fut.result()
            all_summary.extend(res["summary_rows"])
            all_solution_ledger.extend(res.get("solution_ledger_rows", []))
            all_rxn_freq.extend(res["rxn_freq_rows"])
            all_wkf.extend(res["wkf_rows"])
            all_supp.extend(res["supplement_rows"])
            all_masking.extend(res["masking_rows"])
            print(f"[Pilot24 Progress {idx}/{len(tasks)}] Completed {g_name}", flush=True)

    summary_df = pd.DataFrame(all_summary)
    if not summary_df.empty:
        summary_df = summary_df.sort_values(["cohort", "genome", "scenario_id", "biomass_recipe_id"])
    solution_ledger_df = pd.DataFrame(all_solution_ledger)
    if not solution_ledger_df.empty:
        solution_ledger_df = solution_ledger_df.sort_values(["cohort", "genome", "scenario_id", "biomass_recipe_id", "solution_index"])
    rxn_freq_df = pd.DataFrame(all_rxn_freq)
    if not rxn_freq_df.empty:
        rxn_freq_df = rxn_freq_df.sort_values(["cohort", "genome", "scenario_id", "rxn_id"])
    wkf_df = pd.DataFrame(all_wkf)
    if not wkf_df.empty:
        wkf_df = wkf_df.sort_values(["cohort", "genome", "scenario_id", "biomass_recipe_id", "record_type", "target_id"])
    supp_df = pd.DataFrame(all_supp)
    masking_df = pd.DataFrame(all_masking)
    if not masking_df.empty:
        masking_df = masking_df.sort_values(["cohort", "genome", "mask_fraction", "replicate"])

    summary_df.to_csv(args.out_dir / "exp1_pilot24_ensemble_summary.tsv", sep="\t", index=False)
    solution_ledger_df.to_csv(args.out_dir / "exp1_pilot24_per_solution_ledger.tsv", sep="\t", index=False)
    rxn_freq_df.to_csv(args.out_dir / "exp1_pilot24_reaction_stability_and_biomass_sensitivity.tsv", sep="\t", index=False)
    wkf_df.to_csv(args.out_dir / "exp1_pilot24_withdrawal_knockout_fva.tsv", sep="\t", index=False)
    supp_df.to_csv(args.out_dir / "exp1_pilot24_supplement_rescue.tsv", sep="\t", index=False)
    masking_df.to_csv(args.out_dir / "exp1_pilot24_masking_controls.tsv", sep="\t", index=False)

    gate_report = validate_pilot24_scientific_gate(
        summary_df=summary_df,
        rxn_freq_df=rxn_freq_df,
        wkf_df=wkf_df,
        masking_df=masking_df,
        expected_n_genomes=len(tasks),
        expected_genome_ids=[t["genome"] for t in tasks],
        expected_masking_genome_ids=[t["genome"] for t in tasks if t.get("run_masking_control", False)],
        solution_ledger_df=solution_ledger_df,
    )
    gate_report["elapsed_seconds"] = round(time.time() - t0, 2)
    with open(args.out_dir / "exp1_pilot24_qc_gate.json", "w", encoding="utf-8") as fh:
        json.dump(gate_report, fh, indent=2, sort_keys=True)

    print(f"[Stage 2 Pilot24 Gate] {json.dumps(gate_report)}", flush=True)
    if not gate_report["gate_passed"]:
        sys.exit(2)


if __name__ == "__main__":
    main()
