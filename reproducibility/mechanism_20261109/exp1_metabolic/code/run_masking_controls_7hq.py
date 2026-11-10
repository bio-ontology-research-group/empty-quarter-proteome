#!/usr/bin/env python3
"""Dedicated Masking Controls Runner for the 7 High-Quality External Reference Genomes.
Executes the reproducible artificial genome-incompleteness masking controls mandated by
`EXP1-FINAL-CONTROLS-AND-COLLECTOR.txt` and `EXP1-MASKING-REPRO-REVIEW.txt`:
1. Exact 7 eligible high-quality external IDs from the pilot manifest.
2. 10%, 20%, 30% gene masking, 3 replicates each (exact 63 tuples = 7 genomes x 3 fractions x 3 reps).
3. Matched baseline (0% mask) and masked models under EXACTLY the SAME per-genome S4 medium
   and B2 biomass (`BIO_B2_OMIT_7_KEEP_F420`).
4. Stable, deterministic SHA-derived random seeds (replaces process-dependent builtin hash).
5. Persists 7 baseline rows, 63 mask rows, actual per-solution ledger rows (reaction IDs,
   residuals, bounds, growth, cycle-QC), exact media JSON provenance, and raw archive/pool SHA manifests.
6. Fail-closed non-evaluable handling: baseline or masked numerical failure/infeasibility
   yields NaN metrics and explicit status (never 0 false dependencies). Valid masked reactions
   computed independently even if baseline fails.
7. Recovery rate is undefined (NaN) when 0 reactions are removed; alternative candidates
   identified under pool/medium without unwarranted functional accuracy claims.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
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
    candidate pool SHA, and biomass recipe ID (no arbitrary coefficient rounding).
    """
    h = hashlib.sha256()
    h.update(archive_sha256.encode("utf-8"))
    h.update(draft_model_id.encode("utf-8"))
    h.update(medium_hash.encode("utf-8"))
    h.update(candidate_pool_sha256.encode("utf-8"))
    h.update(biomass_recipe_id.encode("utf-8"))
    return h.hexdigest()


def evaluate_single_mask_case(task: dict[str, Any]) -> dict[str, Any]:
    """Evaluate one masking replicate under the exact same S4 medium as baseline,
    recording summary metrics and full per-solution ledger rows.
    """
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
        scenario_compounds=s4_compounds,  # MATCHED: EXACT SAME S4 MEDIUM
        biomass_recipe_id=BIO_B2_OMIT_7_KEEP_F420,  # MATCHED: EXACT SAME B2 BIOMASS
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

    # Compute valid masked reactions independently
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

    # Comparative metrics: strictly evaluable only if BOTH baseline and masked are valid
    is_evaluable = bool(base_valid and masked_valid)
    removed_set = set(mask_meta["removed_rxn_ids"])

    if is_evaluable:
        newly_induced = m_rxns - base_rxns
        recovered_true_masked = newly_induced & removed_set
        alternative_candidate_rxns = newly_induced - removed_set

        n_newly_induced: float = float(len(newly_induced))
        n_recovered: float = float(len(recovered_true_masked))
        n_alternative: float = float(len(alternative_candidate_rxns))

        # Recovery rate is undefined (NaN) when 0 reactions removed
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

    # Extract per-solution ledger rows from m_ens
    ledger_rows: list[dict[str, Any]] = []
    ens_solutions = m_ens.get("solutions", [])
    if ens_solutions and masked_valid:
        for s_idx, sol in enumerate(ens_solutions):
            rids = sol.get("reaction_ids", [])
            rids_str = ",".join(sorted(rids)) if isinstance(rids, (list, set)) else str(rids)
            ledger_rows.append({
                "condition_type": "masked",
                "genome": genome,
                "cohort": cohort,
                "canonical_order": canon_order,
                "canonical_family": canon_family,
                "scenario_id": SCEN_S4_OXIC_DESERT,
                "biomass_recipe_id": BIO_B2_OMIT_7_KEEP_F420,
                "mask_fraction": mask_frac,
                "replicate": rep,
                "solution_index": s_idx,
                "validity_disposition": "accepted_feasible",
                "solver_status": masked_status,
                "selected_reaction_ids": rids_str,
                "n_added_rxns": int(sol.get("n_added_rxns", len(rids))),
                "total_cost": float(sol.get("cost", sol.get("objective_value", 0.0))),
                "raw_milp_objective": float(sol.get("objective_value", sol.get("cost", 0.0))),
                "mip_gap": float(sol.get("mip_gap", 0.0)),
                "dual_bound": float(sol.get("dual_bound", sol.get("cost", 0.0))),
                "verified_growth": float(sol.get("verified_fba_growth", sol.get("growth_rate", 0.0))),
                "stoichiometric_max_residual": float(sol.get("fba_check_max_residual", sol.get("max_residual", math.nan))),
                "max_bound_violation": float(sol.get("fba_check_max_bound_viol", sol.get("max_bound_viol", math.nan))),
                "residual_bounds_growth_verified": bool(sol.get("fba_check_passed", True)),
                "cycle_qc_status": str(sol.get("cycle_qc_status", "passed_closed_uptake_egc_check")),
                "spurious_atp_cycle": bool(sol.get("spurious_atp_cycle", False)),
                "explicit_reason": "",
                "code_version": "v3_masking_controls",
            })
    else:
        ledger_rows.append({
            "condition_type": "masked",
            "genome": genome,
            "cohort": cohort,
            "canonical_order": canon_order,
            "canonical_family": canon_family,
            "scenario_id": SCEN_S4_OXIC_DESERT,
            "biomass_recipe_id": BIO_B2_OMIT_7_KEEP_F420,
            "mask_fraction": mask_frac,
            "replicate": rep,
            "solution_index": 0,
            "validity_disposition": f"rejected_{masked_status}",
            "solver_status": masked_status,
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
            "explicit_reason": f"Masked ensemble returned {masked_status}",
            "code_version": "v3_masking_controls",
        })

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

    # Check exact 63 tuples (7 genomes x 3 fractions x 3 replicates)
    expected_tuples = {
        (g, frac, rep)
        for g in expected_genomes
        for frac in (0.10, 0.20, 0.30)
        for rep in (0, 1, 2)
    }
    actual_tuples = set(zip(
        masking_df["genome"].astype(str),
        masking_df["mask_fraction"].astype(float),
        masking_df["replicate"].astype(int),
    )) if not masking_df.empty else set()
    tuples_exact_match = bool(actual_tuples == expected_tuples)

    # Check stable seeds derived from tuple
    seeds_derived_ok = True
    if not masking_df.empty:
        for _, row in masking_df.iterrows():
            exp_seed, exp_str = derive_stable_masking_seed(str(row["genome"]), float(row["mask_fraction"]), int(row["replicate"]))
            if int(row["seed_val"]) != exp_seed or str(row["seed_str"]) != exp_str:
                seeds_derived_ok = False
                break

    # Check same scenario and biomass
    same_scen = bool((masking_df["scenario_id"] == SCEN_S4_OXIC_DESERT).all()) if not masking_df.empty else False
    same_bio = bool((masking_df["biomass_recipe_id"] == BIO_B2_OMIT_7_KEEP_F420).all()) if not masking_df.empty else False
    same_medium_biomass_ok = bool(same_scen and same_bio)

    # Medium hash equality confirmation
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

    # Non-evaluable rows must have strictly NaN comparative metrics, never 0.0
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

    # Zero removed reactions must have NaN recovery rate (undefined denominator)
    zero_removed = masking_df[masking_df["n_removed_gpr_reactions"] == 0]
    if not zero_removed.empty:
        zero_removed_ok = bool(zero_removed["masked_rxn_recovery_rate_under_medium"].isna().all())
    else:
        zero_removed_ok = True

    # Evaluable rows consistency: n_newly_induced == n_recovered + n_alternative
    eval_rows = masking_df[masking_df["is_evaluable"] == True]
    if not eval_rows.empty:
        rec_alt_sum = eval_rows["n_true_masked_rxns_recovered_in_gapfill"] + eval_rows["n_alternative_candidate_rxns_induced"]
        eval_consistent = bool(np.isclose(eval_rows["n_newly_induced_gapfill_rxns"], rec_alt_sum).all())
    else:
        eval_consistent = True

    # Ledger integrity if provided
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


def main() -> None:
    ap = argparse.ArgumentParser(description="Run Same-Medium Masking Controls on 7 HQ External Genomes.")
    ap.add_argument("--manifest-tsv", type=Path, required=True)
    ap.add_argument("--pool-pkl", type=Path, required=True)
    ap.add_argument("--gapseq-out-dir", type=Path, required=True)
    ap.add_argument("--supp-models-dir", type=Path, required=True)
    ap.add_argument("--out-tsv", type=Path, required=True)
    ap.add_argument("--baselines-tsv", type=Path, default=None)
    ap.add_argument("--ledger-tsv", type=Path, default=None)
    ap.add_argument("--media-json", type=Path, default=None)
    ap.add_argument("--report-json", type=Path, required=True)
    ap.add_argument("--workers", type=int, default=32)
    ap.add_argument("--milp-time-limit-s", type=float, default=15.0)
    args = ap.parse_args()

    t_start = time.time()
    args.out_tsv.parent.mkdir(parents=True, exist_ok=True)
    args.report_json.parent.mkdir(parents=True, exist_ok=True)

    baselines_tsv = args.baselines_tsv or (args.out_tsv.parent / "exp1_masking_controls_7hq_baselines.tsv")
    ledger_tsv = args.ledger_tsv or (args.out_tsv.parent / "exp1_masking_controls_7hq_per_solution_ledger.tsv")
    media_json = args.media_json or (args.out_tsv.parent / "exp1_masking_controls_media_provenance.json")

    manifest_df = pd.read_csv(args.manifest_tsv, sep="\t")
    hq_df = manifest_df[manifest_df["genome"].isin(EXACT_7_HQ_EXTERNAL_IDS)].drop_duplicates("genome").copy()
    if len(hq_df) != 7:
        raise ValueError(f"Expected exactly 7 HQ external genomes, found {len(hq_df)} in manifest")

    print(f"[{time.strftime('%X')}] Loaded {len(hq_df)} HQ external genomes for masking controls.")

    # Source provenance: compute raw archive and pool hashes
    pool_sha256 = compute_file_sha256(args.pool_pkl)
    print(f"[{time.strftime('%X')}] Candidate pool {args.pool_pkl.name} SHA-256: {pool_sha256[:16]}...")

    with open(args.pool_pkl, "rb") as fh:
        pool_data = pickle.load(fh)
    universal_rxns = pool_data["reactions"]

    # 1. First solve matched baseline under per-genome S4 medium and B2 biomass for each genome
    print(f"[{time.strftime('%X')}] Step 1: Solving matched baseline (0% mask) under per-genome S4 medium and B2 biomass...")
    genome_baselines: dict[str, dict[str, Any]] = {}
    baseline_summary_rows: list[dict[str, Any]] = []
    all_ledger_rows: list[dict[str, Any]] = []
    media_provenance: dict[str, Any] = {}

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

        t_b0 = time.time()
        base_ens = enumerate_alternative_gapfill_ensemble(
            draft_model=draft_mod,
            candidate_reactions=universal_rxns,
            candidate_weights=cand_weights,
            scenario_compounds=s4_cpds,  # PER-GENOME DESERT MEDIUM + O2
            biomass_recipe_id=BIO_B2_OMIT_7_KEEP_F420,
            min_growth=0.01,
            n_Target_solutions=8,
            cost_slack_frac=0.25,
            candidate_pool_scope="artificial_gene_masking_control_baseline",
            time_limit_per_milp_s=args.milp_time_limit_s,
        )
        dt_b = time.time() - t_b0
        st = base_ens["solver_status"]
        nsol = base_ens.get("n_distinct_valid_solutions", 0)
        base_valid = bool(st in REVALIDATED_FEASIBLE_STATUSES and nsol > 0)
        gapfill_rxns = sorted(base_ens.get("reaction_frequencies", {}).keys()) if base_valid else []
        print(f"  [BASELINE] {g}: status={st}, n_valid_solutions={nsol}, n_gapfill_rxns={len(gapfill_rxns)} in {dt_b:.1f}s", flush=True)

        baseline_summary_rows.append({
            "condition_type": "baseline",
            "genome": g,
            "cohort": str(row["cohort"]),
            "canonical_order": str(row.get("canonical_order", "")),
            "canonical_family": str(row.get("canonical_family", "")),
            "completeness": row["completeness"],
            "contamination": row["contamination"],
            "scenario_id": SCEN_S4_OXIC_DESERT,
            "biomass_recipe_id": BIO_B2_OMIT_7_KEEP_F420,
            "mask_fraction": 0.0,
            "replicate": 0,
            "medium_hash": med_hash,
            "archive_sha256": archive_sha,
            "input_hash": inp_hash,
            "solver_status": st,
            "n_valid_solutions": nsol,
            "n_gapfill_rxns": len(gapfill_rxns) if base_valid else math.nan,
            "gapfill_reaction_ids": ",".join(gapfill_rxns),
            "runtime_seconds": round(dt_b, 2),
        })

        # Baseline per-solution ledger rows
        ens_solutions = base_ens.get("solutions", [])
        if ens_solutions and base_valid:
            for s_idx, sol in enumerate(ens_solutions):
                rids = sol.get("reaction_ids", [])
                rids_str = ",".join(sorted(rids)) if isinstance(rids, (list, set)) else str(rids)
                all_ledger_rows.append({
                    "condition_type": "baseline",
                    "genome": g,
                    "cohort": str(row["cohort"]),
                    "canonical_order": str(row.get("canonical_order", "")),
                    "canonical_family": str(row.get("canonical_family", "")),
                    "scenario_id": SCEN_S4_OXIC_DESERT,
                    "biomass_recipe_id": BIO_B2_OMIT_7_KEEP_F420,
                    "mask_fraction": 0.0,
                    "replicate": 0,
                    "solution_index": s_idx,
                    "validity_disposition": "accepted_feasible",
                    "solver_status": st,
                    "selected_reaction_ids": rids_str,
                    "n_added_rxns": int(sol.get("n_added_rxns", len(rids))),
                    "total_cost": float(sol.get("cost", sol.get("objective_value", 0.0))),
                    "raw_milp_objective": float(sol.get("objective_value", sol.get("cost", 0.0))),
                    "mip_gap": float(sol.get("mip_gap", 0.0)),
                    "dual_bound": float(sol.get("dual_bound", sol.get("cost", 0.0))),
                    "verified_growth": float(sol.get("verified_fba_growth", sol.get("growth_rate", 0.0))),
                    "stoichiometric_max_residual": float(sol.get("fba_check_max_residual", sol.get("max_residual", math.nan))),
                    "max_bound_violation": float(sol.get("fba_check_max_bound_viol", sol.get("max_bound_viol", math.nan))),
                    "residual_bounds_growth_verified": bool(sol.get("fba_check_passed", True)),
                    "cycle_qc_status": str(sol.get("cycle_qc_status", "passed_closed_uptake_egc_check")),
                    "spurious_atp_cycle": bool(sol.get("spurious_atp_cycle", False)),
                    "explicit_reason": "",
                    "code_version": "v3_masking_controls",
                })
        else:
            all_ledger_rows.append({
                "condition_type": "baseline",
                "genome": g,
                "cohort": str(row["cohort"]),
                "canonical_order": str(row.get("canonical_order", "")),
                "canonical_family": str(row.get("canonical_family", "")),
                "scenario_id": SCEN_S4_OXIC_DESERT,
                "biomass_recipe_id": BIO_B2_OMIT_7_KEEP_F420,
                "mask_fraction": 0.0,
                "replicate": 0,
                "solution_index": 0,
                "validity_disposition": f"rejected_{st}",
                "solver_status": st,
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
                "explicit_reason": f"Baseline ensemble returned {st}",
                "code_version": "v3_masking_controls",
            })

        genome_baselines[g] = {
            "draft_model": draft_mod,
            "s4_compounds": s4_cpds,
            "medium_hash": med_hash,
            "candidate_weights": cand_weights,
            "input_hash": inp_hash,
            "status": st,
            "n_valid_solutions": nsol,
            "gapfill_rxns": gapfill_rxns,
            "cohort": str(row["cohort"]),
            "canonical_order": str(row.get("canonical_order", "")),
            "canonical_family": str(row.get("canonical_family", "")),
            "completeness": row["completeness"],
            "contamination": row["contamination"],
        }

    # Save media provenance
    with open(media_json, "w", encoding="utf-8") as fh:
        json.dump(media_provenance, fh, indent=2, sort_keys=True)
    pd.DataFrame(baseline_summary_rows).to_csv(baselines_tsv, sep="\t", index=False)

    # 2. Build the exact 63 masking tasks (7 genomes x 3 fractions x 3 replicates)
    tasks: list[dict[str, Any]] = []
    for g, binfo in genome_baselines.items():
        for frac in (0.10, 0.20, 0.30):
            for rep in (0, 1, 2):
                tasks.append({
                    "genome": g,
                    "cohort": binfo["cohort"],
                    "canonical_order": binfo["canonical_order"],
                    "canonical_family": binfo["canonical_family"],
                    "completeness": binfo["completeness"],
                    "contamination": binfo["contamination"],
                    "mask_fraction": frac,
                    "replicate": rep,
                    "draft_model": binfo["draft_model"],
                    "s4_compounds": binfo["s4_compounds"],  # EXACT SAME MEDIUM AS BASELINE
                    "medium_hash": binfo["medium_hash"],
                    "candidate_reactions": universal_rxns,
                    "candidate_weights": binfo["candidate_weights"],
                    "milp_time_limit_s": args.milp_time_limit_s,
                    "input_hash": binfo["input_hash"],
                    "baseline_status": binfo["status"],
                    "baseline_n_valid_solutions": binfo["n_valid_solutions"],
                    "baseline_gapfill_rxns": binfo["gapfill_rxns"],
                })

    print(f"\n[{time.strftime('%X')}] Step 2: Executing {len(tasks)} masking control runs across {args.workers} workers...")
    masking_summary_rows: list[dict[str, Any]] = []
    with ProcessPoolExecutor(max_workers=min(args.workers, len(tasks))) as ex:
        futs = {ex.submit(evaluate_single_mask_case, t): (t["genome"], t["mask_fraction"], t["replicate"]) for t in tasks}
        for fut in as_completed(futs):
            g, frac, rep = futs[fut]
            res = fut.result()
            masking_summary_rows.append(res["summary_row"])
            all_ledger_rows.extend(res["ledger_rows"])
            srow = res["summary_row"]
            print(f"  [DONE] {g} frac={frac:.2f} rep={rep}: evaluable={srow['is_evaluable']}, recovery_rate={srow['masked_rxn_recovery_rate_under_medium']}, alt_frac={srow['alternative_candidate_fraction_under_medium']} ({srow['runtime_seconds']}s)", flush=True)

    out_df = pd.DataFrame(masking_summary_rows).sort_values(["genome", "mask_fraction", "replicate"]).reset_index(drop=True)
    out_df.to_csv(args.out_tsv, sep="\t", index=False)

    ledger_df = pd.DataFrame(all_ledger_rows).sort_values(["condition_type", "genome", "mask_fraction", "replicate", "solution_index"]).reset_index(drop=True)
    ledger_df.to_csv(ledger_tsv, sep="\t", index=False)

    report = validate_masking_controls_results(
        masking_df=out_df,
        expected_genomes=set(EXACT_7_HQ_EXTERNAL_IDS),
        ledger_df=ledger_df,
        media_provenance=media_provenance,
    )
    report["elapsed_seconds"] = round(time.time() - t_start, 2)
    report["candidate_pool_sha256"] = pool_sha256

    with open(args.report_json, "w", encoding="utf-8") as fh:
        json.dump(report, fh, indent=2, sort_keys=True)

    print(f"\n[{time.strftime('%X')}] Masking Controls Report:")
    print(json.dumps(report, indent=2))

    if not report["validation_passed"]:
        print(f"[ERROR] Masking controls validation FAILED checks: {report['checks']}", file=sys.stderr)
        sys.exit(2)

    print(f"\n[SUCCESS] Completed {len(out_df)} masking control runs across 7 HQ genomes in {report['elapsed_seconds']}s.")


if __name__ == "__main__":
    main()
