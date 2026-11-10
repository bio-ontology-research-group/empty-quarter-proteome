#!/usr/bin/env python3
"""Standalone Frozen Output & Reproducibility Gate Verifier for Experiment 1 Metabolic Modelling.
Validates the frozen outputs, manifests, receipts, and scientific boundaries staged in
`reproducibility/mechanism_20261109/exp1_metabolic/` per the project acceptance contract
and `EXP1-FINAL-PACKAGE-ACCEPTANCE-REVIEW.txt`:
1. Exact 484-genome cohort universe & ID uniqueness: 237 Rub' al Khali MAGs (`EQ_desert`) +
   247 external reference genomes (`external_comparator`: 218 GTDB r232 + 29 MGnify) with heterogeneous
   quality (7 external genomes meeting numeric >=90% completeness / <=5% contamination for masking controls).
2. Candidate reaction pool integrity: 2,271 reactions audited, required conflict files present and empty
   (0 reaction ID conflicts, 0 species conflicts), and explicit unknown formula/charge status accounting.
3. Biomass template audit: consistent archaeal (`Methanosarcina` iAF692) cofactor assignments across all 484 genomes.
4. Full 484-genome production case & per-solution ledger verification:
   - Accounts for all 12,580 planned cases across 484 genomes (10,578 feasible ensembles, 1,936 fixed-medium
     infeasible, 56 FBA verification-failed, 10 solver timeouts = 66 explicit unknown/failed cases retained as NaN).
   - Verifies all 93,808 per-solution ledger rows and 91,806 `accepted_feasible` solutions (finite verified growth
     >= 0.01*0.999, stoichiometric residual <= 1e-5, bound violation <= 1e-5, closed-uptake cycle check passed,
     no duplicate solution keys, exact per-case accepted count correspondence, and exact 8-shard source hash
     reconstruction). Explicitly notes that raw flux vectors are not persisted (`raw_flux_replay_not_claimed`)
     and preserves `exp1_pilot24_qc_gate.json` as an immutable historical Stage 2 record (not an aggregate >=95%
     or zero-timeout final acceptance gate).
5. Four changed witness Slurm validation: exact model IDs and status revalidated.
6. Masking controls baseline, provenance & v4 post-hoc validation: exact 7 eligible external IDs, matched S4
   medium hash, 63 v4 cases validated.
7. Slurm execution receipts, Stage 3 collector v4 receipt, Stage 4 v5 search receipt, and bounded post-hoc
   3,713-unique-pair scientific disposition & 47-sample assembly input provenance audit.
8. Whole-package scoped SHA-256 manifest verification: verifies paths, byte sizes, and SHA-256 hashes for all
   owned core package files while excluding mutable caches (`.pytest_cache`, `__pycache__`, `*.pyc`) and
   independently owned subpackages (`covariates/`, `interpretation/`, `counterfactual_review/`).
"""
from __future__ import annotations

import gzip
import hashlib
import json
import math
import sys
from pathlib import Path
from typing import Any

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent
EVIDENCE_DIR = REPO_ROOT / "compact_evidence"
MANIFESTS_DIR = REPO_ROOT / "manifests"

EXPECTED_484_GENOME_COUNT = 484
EXPECTED_EQ_COUNT = 237
EXPECTED_EXT_COUNT = 247
EXPECTED_EXT_GTDB_COUNT = 218
EXPECTED_EXT_MGNIFY_COUNT = 29

EXPECTED_FULL484_CASES = 12580
EXPECTED_FULL484_LEDGER_ROWS = 93808
EXPECTED_FULL484_ACCEPTED_SOLUTIONS = 91806
EXPECTED_FULL484_STATUS_COUNTS = {
    "heuristic_incumbent_bounded": 10578,
    "infeasible_under_fixed_medium": 1936,
    "fba_verification_failed": 56,
    "timeout_no_incumbent": 10,
}
OPTIONAL_S0_ABSENT_GENOMES = {
    "GCA_016180755.1",
    "S22PRr1__self_SemiBin_29",
}

REVALIDATED_FEASIBLE_STATUSES = {
    "optimal",
    "optimal_after_cycle_qc_filtering",
    "optimal_draft_no_gapfill_needed",
    "heuristic_incumbent_bounded",
    "post_pruning_cost_changed_not_proven_optimal",
}

EXPECTED_24_CORE_SCENARIO_BIOMASS_COMBINATIONS = {
    ("S1_anoxic_minimal_h2_co2_nh4", "B0_methanosarcina_template_baseline"),
    ("S1_anoxic_minimal_h2_co2_nh4", "B2_perturbation_omit_7_non_f420_template_cofactors_keep_f420"),
    ("S2_oxic_minimal_h2_co2_nh4", "B0_methanosarcina_template_baseline"),
    ("S2_oxic_minimal_h2_co2_nh4", "B2_perturbation_omit_7_non_f420_template_cofactors_keep_f420"),
    ("S3_anoxic_desert_oa_aa", "B0_methanosarcina_template_baseline"),
    ("S3_anoxic_desert_oa_aa", "B1_perturbation_omit_8_template_cofactors"),
    ("S3_anoxic_desert_oa_aa", "B2_perturbation_omit_7_non_f420_template_cofactors_keep_f420"),
    ("S3_anoxic_desert_oa_aa", "B3_perturbation_omit_template_cofactors_and_surface_glycolipids"),
    ("S3_anoxic_desert_oa_aa", "B4_perturbation_core_amino_nucleotide_energy_skeleton"),
    ("S4_oxic_desert_oa_aa", "B0_methanosarcina_template_baseline"),
    ("S4_oxic_desert_oa_aa", "B1_perturbation_omit_8_template_cofactors"),
    ("S4_oxic_desert_oa_aa", "B2_perturbation_omit_7_non_f420_template_cofactors_keep_f420"),
    ("S4_oxic_desert_oa_aa", "B3_perturbation_omit_template_cofactors_and_surface_glycolipids"),
    ("S4_oxic_desert_oa_aa", "B4_perturbation_core_amino_nucleotide_energy_skeleton"),
    ("S5_anoxic_defined_aa_precursors", "B0_methanosarcina_template_baseline"),
    ("S5_anoxic_defined_aa_precursors", "B1_perturbation_omit_8_template_cofactors"),
    ("S5_anoxic_defined_aa_precursors", "B2_perturbation_omit_7_non_f420_template_cofactors_keep_f420"),
    ("S5_anoxic_defined_aa_precursors", "B3_perturbation_omit_template_cofactors_and_surface_glycolipids"),
    ("S5_anoxic_defined_aa_precursors", "B4_perturbation_core_amino_nucleotide_energy_skeleton"),
    ("S6_oxic_defined_aa_precursors", "B0_methanosarcina_template_baseline"),
    ("S6_oxic_defined_aa_precursors", "B1_perturbation_omit_8_template_cofactors"),
    ("S6_oxic_defined_aa_precursors", "B2_perturbation_omit_7_non_f420_template_cofactors_keep_f420"),
    ("S6_oxic_defined_aa_precursors", "B3_perturbation_omit_template_cofactors_and_surface_glycolipids"),
    ("S6_oxic_defined_aa_precursors", "B4_perturbation_core_amino_nucleotide_energy_skeleton"),
}
OPTIONAL_S0_COMBINATIONS = {
    ("S0_gapseq_native_predicted", "B0_methanosarcina_template_baseline"),
    ("S0_gapseq_native_predicted", "B2_perturbation_omit_7_non_f420_template_cofactors_keep_f420"),
}

EXACT_7_HQ_EXTERNAL_IDS = {
    "MGYG000519815",
    "MGYG000520046",
    "MGYG000519748",
    "MGYG000510876",
    "MGYG000495096",
    "MGYG000526805",
    "MGYG000523080",
}

EXCLUDED_SUBPACKAGE_PREFIXES = (
    "covariates/",
    "interpretation/",
    "counterfactual_review/",
)
EXCLUDED_CACHE_PARTS = {
    ".pytest_cache",
    "__pycache__",
}


def compute_sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        while chunk := fh.read(1024 * 1024):
            h.update(chunk)
    return h.hexdigest()


def verify_gate1_genome_inventory(inv_file: Path | None = None) -> tuple[bool, str]:
    """Gate 1: Exact 484-genome cohort universe & ID uniqueness (237 EQ_desert + 247 external_comparator: 218 GTDB + 29 MGnify)."""
    inv_file = inv_file or (EVIDENCE_DIR / "exp1_model_inventory_and_qc_484.tsv")
    if not inv_file.exists():
        return False, f"Missing inventory file: {inv_file}"
    df = pd.read_csv(inv_file, sep="\t")
    if len(df) != EXPECTED_484_GENOME_COUNT:
        return False, f"Expected {EXPECTED_484_GENOME_COUNT} rows, found {len(df)}"
    if df["genome"].isna().any() or (df["genome"].astype(str).str.strip() == "").any():
        return False, "Inventory contains empty or NaN genome IDs"
    if df["genome"].astype(str).duplicated().any():
        return False, "Inventory contains duplicate genome IDs"

    eq_count = int((df["cohort"] == "EQ_desert").sum())
    ext_count = int((df["cohort"] == "external_comparator").sum())
    if eq_count != EXPECTED_EQ_COUNT:
        return False, f"Expected {EXPECTED_EQ_COUNT} EQ MAGs, found {eq_count}"
    if ext_count != EXPECTED_EXT_COUNT:
        return False, f"Expected {EXPECTED_EXT_COUNT} external reference genomes, found {ext_count}"

    if "source" in df.columns:
        src_counts = df["source"].value_counts().to_dict()
        if src_counts.get("EQ_MAG", 0) != EXPECTED_EQ_COUNT:
            return False, f"Expected {EXPECTED_EQ_COUNT} EQ_MAG source rows, found {src_counts.get('EQ_MAG', 0)}"
        if src_counts.get("GTDB_r232_rep", 0) != EXPECTED_EXT_GTDB_COUNT:
            return False, f"Expected {EXPECTED_EXT_GTDB_COUNT} GTDB_r232_rep source rows, found {src_counts.get('GTDB_r232_rep', 0)}"
        if src_counts.get("MGnify", 0) != EXPECTED_EXT_MGNIFY_COUNT:
            return False, f"Expected {EXPECTED_EXT_MGNIFY_COUNT} MGnify source rows, found {src_counts.get('MGnify', 0)}"

    ext_hq_ids = set(
        df[
            (df["cohort"] == "external_comparator")
            & (df["completeness"] >= 90.0)
            & (df["contamination"] <= 5.0)
        ]["genome"].astype(str)
    )
    if ext_hq_ids != EXACT_7_HQ_EXTERNAL_IDS:
        return False, f"External numeric >=90%/<=5% set mismatch: {ext_hq_ids}"

    return True, (
        f"484 unique genomes verified ({eq_count} EQ_MAG + {ext_count} external_comparator "
        f"[{EXPECTED_EXT_GTDB_COUNT} GTDB_r232_rep + {EXPECTED_EXT_MGNIFY_COUNT} MGnify]; "
        f"7 external references meet numeric >=90%/<=5% masking filter)"
    )


def audit_candidate_pool_files(
    audit_file: Path,
    rxn_conflicts_file: Path,
    sp_conflicts_file: Path,
) -> tuple[bool, str]:
    """Fail-closed candidate pool and conflict file verifier (`EXP1-FINAL-PACKAGE-ACCEPTANCE-REVIEW.txt` Item 6)."""
    for req_file, label in (
        (audit_file, "candidate pool audit"),
        (rxn_conflicts_file, "reaction conflicts"),
        (sp_conflicts_file, "species conflicts"),
    ):
        if not req_file.exists():
            return False, f"Missing required {label} file: {req_file}"

    audit_df = pd.read_csv(audit_file, sep="\t")
    if len(audit_df) != 2271:
        return False, f"Expected 2,271 candidate pool reactions, found {len(audit_df)}"
    if audit_df["rxn_id"].isna().any() or audit_df["rxn_id"].astype(str).duplicated().any():
        return False, "Candidate pool audit contains missing or duplicate rxn_id entries"
    if (audit_df["has_stoichiometry_or_compartment_conflict"] != False).any():
        return False, "Candidate pool audit contains stoichiometry or compartment conflicts"

    bal_counts = audit_df["balance_status"].value_counts(dropna=False).to_dict()
    expected_bal = {
        "strict_atom_and_charge_balanced": 1836,
        "unknown_formula_or_charge": 228,
        "boundary_or_biomass_exempt": 206,
        "elemental_imbalance": 1,
    }
    if bal_counts != expected_bal:
        return False, f"Candidate pool balance_status mismatch: {bal_counts} != {expected_bal}"

    elig_counts = audit_df["eligible_for_mechanism_gapfill"].value_counts(dropna=False).to_dict()
    if elig_counts != {True: 2270, False: 1}:
        return False, f"Candidate pool eligibility counts mismatch: {elig_counts}"

    rxn_conf_df = pd.read_csv(rxn_conflicts_file, sep="\t")
    if not rxn_conf_df.empty:
        return False, f"Found {len(rxn_conf_df)} reaction conflicts in candidate pool"

    sp_conf_df = pd.read_csv(sp_conflicts_file, sep="\t")
    if not sp_conf_df.empty:
        return False, f"Found {len(sp_conf_df)} species conflicts in candidate pool"

    return True, (
        "Candidate pool clean: 2,271 reactions audited (2,270 eligible, 1 elemental imbalance excluded; "
        "1,836 strict balanced, 228 unknown formula/charge, 206 boundary/biomass exempt), "
        "both conflict files present with 0 conflicts"
    )


def verify_gate2_candidate_pool() -> tuple[bool, str]:
    """Gate 2: Candidate reaction pool integrity (2,271 rows, required conflict files present & empty, balance statuses preserved)."""
    return audit_candidate_pool_files(
        EVIDENCE_DIR / "exp1_candidate_reaction_pool_audit.tsv",
        EVIDENCE_DIR / "exp1_candidate_reaction_pool_conflicts.tsv",
        EVIDENCE_DIR / "exp1_candidate_species_conflicts.tsv",
    )


def verify_gate3_biomass_template_audit() -> tuple[bool, str]:
    """Gate 3: Biomass template audit across all 484 unique genomes."""
    bio_file = EVIDENCE_DIR / "exp1_biomass_template_audit_484.tsv"
    inv_file = EVIDENCE_DIR / "exp1_model_inventory_and_qc_484.tsv"
    if not bio_file.exists():
        return False, f"Missing biomass template audit file: {bio_file}"
    if not inv_file.exists():
        return False, f"Missing inventory file: {inv_file}"
    bio_df = pd.read_csv(bio_file, sep="\t")
    inv_df = pd.read_csv(inv_file, sep="\t")
    if len(bio_df) != EXPECTED_484_GENOME_COUNT or bio_df["genome"].astype(str).duplicated().any():
        return False, f"Expected {EXPECTED_484_GENOME_COUNT} unique biomass rows, found {len(bio_df)}"
    if set(bio_df["genome"].astype(str)) != set(inv_df["genome"].astype(str)):
        return False, "Biomass template audit genome IDs do not match 484 inventory genome IDs"
    if (bio_df["n_all8_template_trace_cofactors_demanded"] != 8).any():
        return False, "Not all 484 draft models demand the 8 Methanosarcina template trace cofactors in B0"
    if (bio_df["n_non_f420_methanosarcina_cofactors_demanded"] != 7).any():
        return False, "Not all 484 draft models demand the 7 non-F420 Methanosarcina template cofactors in B0"
    if (bio_df["f420_demanded_in_template"] != True).any():
        return False, "Not all 484 draft models demand F420 (cpd00649) in B0"
    return True, "Biomass template audit verified across all 484 unique genomes (8 Methanosarcina template trace cofactors in B0)"


def audit_full484_summary_and_solution_ledger(
    inv_df: pd.DataFrame,
    summary_df: pd.DataFrame,
    ledger_df: pd.DataFrame,
) -> tuple[bool, str, dict[str, Any]]:
    """Perform full per-case and per-solution verification on a summary and solution ledger table.
    Enforces the scientific completion rule (`EXP1-FINAL-PACKAGE-ACCEPTANCE-REVIEW.txt` Items 1 & 2):
    - No arbitrary >=95% aggregate completion or zero-timeout requirement.
    - All planned cases must be present without duplicates; infeasible and unknown/timeout/failed cases
      must be explicitly retained with NaN summary metrics (never false zeros).
    - Every `accepted_feasible` solution row in the ledger must have finite verified numerics,
      `verified_growth >= 0.01 * 0.999`, `stoichiometric_max_residual <= 1e-5`,
      `max_bound_violation <= 1e-5`, `residual_bounds_growth_verified == True`,
      `cycle_qc_status == 'passed_closed_uptake_egc_check'`, `spurious_atp_cycle == False`,
      exact reaction count agreement, no duplicate solution keys, and exact per-case count agreement
      with `n_distinct_valid_solutions` in `summary_df`.
    """
    if inv_df.empty or inv_df["genome"].isna().any() or inv_df["genome"].astype(str).duplicated().any():
        return False, "Inventory must contain unique non-null genome IDs", {}

    expected_genomes = set(inv_df["genome"].astype(str))
    summary_genomes = set(summary_df["genome"].astype(str)) if not summary_df.empty else set()
    if summary_genomes != expected_genomes:
        return False, f"Summary genome set mismatch: expected {len(expected_genomes)}, got {len(summary_genomes)}", {}

    if summary_df.duplicated(subset=["genome", "scenario_id", "biomass_recipe_id"]).any():
        return False, "Duplicate (genome, scenario_id, biomass_recipe_id) rows in summary", {}

    if ledger_df.duplicated(subset=["genome", "scenario_id", "biomass_recipe_id", "solution_index"]).any():
        return False, "Duplicate (genome, scenario_id, biomass_recipe_id, solution_index) rows in solution ledger", {}

    # Check 24 core combinations for every genome and optional S0 combinations
    for g, grp in summary_df.groupby(summary_df["genome"].astype(str)):
        pairs = set(zip(grp["scenario_id"].astype(str), grp["biomass_recipe_id"].astype(str)))
        missing_core = EXPECTED_24_CORE_SCENARIO_BIOMASS_COMBINATIONS - pairs
        if missing_core:
            return False, f"Genome {g} missing core (scenario, biomass) combinations: {sorted(missing_core)[:3]}", {}
        extra = pairs - EXPECTED_24_CORE_SCENARIO_BIOMASS_COMBINATIONS
        if not extra.issubset(OPTIONAL_S0_COMBINATIONS):
            return False, f"Genome {g} has unexpected (scenario, biomass) combinations: {sorted(extra)}", {}
        if extra and extra != OPTIONAL_S0_COMBINATIONS:
            return False, f"Genome {g} has partial S0 combinations: {sorted(extra)}", {}

    # Check case universe equality between summary and ledger
    summary_cases = set(
        zip(
            summary_df["genome"].astype(str),
            summary_df["scenario_id"].astype(str),
            summary_df["biomass_recipe_id"].astype(str),
        )
    )
    ledger_cases = set(
        zip(
            ledger_df["genome"].astype(str),
            ledger_df["scenario_id"].astype(str),
            ledger_df["biomass_recipe_id"].astype(str),
        )
    )
    if summary_cases != ledger_cases:
        return False, "Case universe mismatch between summary_df and ledger_df", {}

    # Check summary statuses and strict NaN retention on non-feasible / unknown cases
    allowed_summary_statuses = REVALIDATED_FEASIBLE_STATUSES | {
        "infeasible_under_fixed_medium",
        "fba_verification_failed",
        "timeout_no_incumbent",
    }
    if not set(summary_df["solver_status"].astype(str)).issubset(allowed_summary_statuses):
        bad_st = set(summary_df["solver_status"].astype(str)) - allowed_summary_statuses
        return False, f"Unexpected summary solver_status values: {sorted(bad_st)}", {}

    feas_mask = summary_df["solver_status"].astype(str).isin(REVALIDATED_FEASIBLE_STATUSES)
    non_feas_df = summary_df[~feas_mask]
    if (non_feas_df["n_distinct_valid_solutions"].fillna(0).astype(int) != 0).any():
        return False, "Non-feasible/unknown summary row has n_distinct_valid_solutions != 0", {}

    na_metric_cols = [
        "min_added_rxns_count",
        "max_added_rxns_count",
        "n_rxns_observed_in_all_enumerated_models",
        "n_rxns_essential_in_declared_supernetwork",
        "min_total_cost",
        "max_ensemble_cost",
    ]
    for col in na_metric_cols:
        if col in non_feas_df.columns and non_feas_df[col].notna().any():
            return False, f"Non-feasible/unknown summary rows must have NaN for {col}", {}

    feas_df = summary_df[feas_mask]
    if (feas_df["n_distinct_valid_solutions"].fillna(0).astype(int) <= 0).any():
        return False, "Feasible summary row has n_distinct_valid_solutions <= 0", {}
    for col in na_metric_cols:
        if col in feas_df.columns and feas_df[col].isna().any():
            return False, f"Feasible summary row has NaN for {col}", {}

    # Check per-solution ledger rows
    acc_mask = ledger_df["validity_disposition"].astype(str) == "accepted_feasible"
    acc_df = ledger_df[acc_mask]
    unres_df = ledger_df[~acc_mask]

    allowed_unres_dispositions = {
        "unresolved_infeasible_under_fixed_medium",
        "unresolved_fba_verification_failed",
        "unresolved_timeout_no_incumbent",
        "rejected_cycle_qc",
        "rejected_fba_verification",
    }
    if not set(unres_df["validity_disposition"].astype(str)).issubset(allowed_unres_dispositions):
        bad_disp = set(unres_df["validity_disposition"].astype(str)) - allowed_unres_dispositions
        return False, f"Unexpected non-accepted ledger validity_disposition: {sorted(bad_disp)}", {}

    if not set(acc_df["solver_status"].astype(str)).issubset(REVALIDATED_FEASIBLE_STATUSES):
        return False, "Accepted ledger row has non-feasible solver_status", {}

    for num_col in ("verified_growth", "stoichiometric_max_residual", "max_bound_violation", "total_cost", "n_added_rxns"):
        vals = pd.to_numeric(acc_df[num_col], errors="coerce")
        if vals.isna().any() or not vals.map(math.isfinite).all():
            return False, f"Accepted ledger row has NaN or non-finite {num_col}", {}

    gr = pd.to_numeric(acc_df["verified_growth"], errors="coerce")
    if (gr < 0.01 * 0.999).any():
        return False, f"Accepted ledger row fails verified_growth >= 0.00999 (min={gr.min()})", {}

    res = pd.to_numeric(acc_df["stoichiometric_max_residual"], errors="coerce")
    if ((res < 0.0) | (res > 1e-5)).any():
        return False, f"Accepted ledger row fails stoichiometric_max_residual <= 1e-5 (max={res.max()})", {}

    bnd = pd.to_numeric(acc_df["max_bound_violation"], errors="coerce")
    if ((bnd < 0.0) | (bnd > 1e-5)).any():
        return False, f"Accepted ledger row fails max_bound_violation <= 1e-5 (max={bnd.max()})", {}

    if (acc_df["residual_bounds_growth_verified"] != True).any():
        return False, "Accepted ledger row has residual_bounds_growth_verified != True", {}

    if (acc_df["cycle_qc_status"].astype(str) != "passed_closed_uptake_egc_check").any():
        return False, "Accepted ledger row failed closed_uptake_egc_check", {}

    if (acc_df["spurious_atp_cycle"] != False).any():
        return False, "Accepted ledger row has spurious_atp_cycle != False", {}

    # Check reaction count matches selected_reaction_ids
    sel_series = acc_df["selected_reaction_ids"].fillna("").astype(str)
    rxn_counts = sel_series.map(lambda s: 0 if not s.strip() else len([x for x in s.split(";") if x.strip()]))
    recorded_counts = pd.to_numeric(acc_df["n_added_rxns"], errors="coerce").astype(int)
    if (rxn_counts != recorded_counts).any():
        return False, "Accepted ledger row n_added_rxns does not match selected_reaction_ids length", {}

    # Check exact accepted count correspondence between summary and ledger
    acc_counts_by_case = (
        acc_df.groupby(
            [
                acc_df["genome"].astype(str),
                acc_df["scenario_id"].astype(str),
                acc_df["biomass_recipe_id"].astype(str),
            ]
        )
        .size()
        .to_dict()
    )
    for _, srow in summary_df.iterrows():
        key = (str(srow["genome"]), str(srow["scenario_id"]), str(srow["biomass_recipe_id"]))
        exp_n = int(srow["n_distinct_valid_solutions"]) if str(srow["solver_status"]) in REVALIDATED_FEASIBLE_STATUSES else 0
        act_n = int(acc_counts_by_case.get(key, 0))
        if exp_n != act_n:
            return False, f"Case {key} accepted solution count mismatch: summary={exp_n}, ledger={act_n}", {}

    status_counts = summary_df["solver_status"].astype(str).value_counts().to_dict()
    stats_out = {
        "n_genomes": len(expected_genomes),
        "n_cases": len(summary_df),
        "n_ledger_rows": len(ledger_df),
        "n_accepted_solutions": len(acc_df),
        "solver_status_counts": status_counts,
        "raw_flux_replay_claimed": False,
    }
    return True, "Full case and per-solution ledger audit passed", stats_out


def verify_gate4_full484_case_and_solution_ledgers() -> tuple[bool, str]:
    """Gate 4: Full 484-genome production case & per-solution ledger audit + historical Pilot24 receipt preservation."""
    inv_file = EVIDENCE_DIR / "exp1_model_inventory_and_qc_484.tsv"
    summary_file = EVIDENCE_DIR / "exp1_full484_ensemble_summary.tsv.gz"
    ledger_file = EVIDENCE_DIR / "exp1_full484_per_solution_ledger.tsv.gz"
    root_audit_file = EVIDENCE_DIR / "provenance" / "exp1_full484_independent_table_audit.json"
    pilot_gate_file = EVIDENCE_DIR / "exp1_pilot24_qc_gate.json"

    for req in (inv_file, summary_file, ledger_file, root_audit_file, pilot_gate_file):
        if not req.exists():
            return False, f"Missing required Gate 4 file: {req}"

    inv_df = pd.read_csv(inv_file, sep="\t")
    summary_df = pd.read_csv(summary_file, sep="\t")
    ledger_df = pd.read_csv(ledger_file, sep="\t")

    ok, msg, stats_out = audit_full484_summary_and_solution_ledger(inv_df, summary_df, ledger_df)
    if not ok:
        return False, msg

    if stats_out["n_genomes"] != EXPECTED_484_GENOME_COUNT:
        return False, f"Expected {EXPECTED_484_GENOME_COUNT} genomes, got {stats_out['n_genomes']}"
    if stats_out["n_cases"] != EXPECTED_FULL484_CASES:
        return False, f"Expected {EXPECTED_FULL484_CASES} planned cases, got {stats_out['n_cases']}"
    if stats_out["n_ledger_rows"] != EXPECTED_FULL484_LEDGER_ROWS:
        return False, f"Expected {EXPECTED_FULL484_LEDGER_ROWS} ledger rows, got {stats_out['n_ledger_rows']}"
    if stats_out["n_accepted_solutions"] != EXPECTED_FULL484_ACCEPTED_SOLUTIONS:
        return False, f"Expected {EXPECTED_FULL484_ACCEPTED_SOLUTIONS} accepted solutions, got {stats_out['n_accepted_solutions']}"
    if stats_out["solver_status_counts"] != EXPECTED_FULL484_STATUS_COUNTS:
        return False, f"Full484 status counts mismatch: {stats_out['solver_status_counts']} != {EXPECTED_FULL484_STATUS_COUNTS}"

    # Verify optional S0 absence is restricted to the exact 2 empty-native-medium genomes
    s0_genomes = set(
        summary_df[summary_df["scenario_id"] == "S0_gapseq_native_predicted"]["genome"].astype(str)
    )
    absent_s0 = set(inv_df["genome"].astype(str)) - s0_genomes
    if absent_s0 != OPTIONAL_S0_ABSENT_GENOMES:
        return False, f"Optional S0 absent genomes mismatch: {absent_s0} != {OPTIONAL_S0_ABSENT_GENOMES}"

    # Verify exact reconstruction of the 8 shard solution ledger SHA-256 hashes recorded in Root's independent audit
    with open(root_audit_file, "r", encoding="utf-8") as fh:
        root_audit = json.load(fh)
    if (
        root_audit.get("genomes") != EXPECTED_484_GENOME_COUNT
        or root_audit.get("cases") != EXPECTED_FULL484_CASES
        or root_audit.get("ledger_rows") != EXPECTED_FULL484_LEDGER_ROWS
        or root_audit.get("accepted_solutions") != EXPECTED_FULL484_ACCEPTED_SOLUTIONS
    ):
        return False, "Root independent full484 table audit JSON header counts mismatch"

    sorted_inv = inv_df.sort_values("genome").reset_index(drop=True)
    genome_to_shard: dict[str, int] = {}
    for tid in range(8):
        for g in sorted_inv.iloc[tid::8]["genome"].astype(str):
            genome_to_shard[g] = tid

    with gzip.open(ledger_file, "rb") as gz:
        raw_lines = gz.read().splitlines(keepends=True)
    header_line = raw_lines[0]
    shard_lines: dict[int, list[bytes]] = {tid: [header_line] for tid in range(8)}
    for line in raw_lines[1:]:
        g = line.split(b"\t", 1)[0].decode("utf-8")
        shard_lines[genome_to_shard[g]].append(line)

    for sh_entry in root_audit.get("shards", []):
        tid = int(sh_entry["shard"])
        exp_led_sha = sh_entry["source_sha256"][f"exp1_shard_{tid:02d}_solution_ledger.tsv"]
        act_led_sha = hashlib.sha256(b"".join(shard_lines[tid])).hexdigest()
        if act_led_sha != exp_led_sha:
            return False, f"Reconstructed shard {tid:02d} solution ledger SHA-256 mismatch: {act_led_sha} != {exp_led_sha}"

    # Verify historical Pilot24 QC gate JSON is preserved intact (as historical Stage 2 record only, NOT final acceptance rule)
    with open(pilot_gate_file, "r", encoding="utf-8") as fh:
        pilot_data = json.load(fh)
    if pilot_data.get("n_genomes_evaluated") != 31 or pilot_data.get("n_summary_rows") != 806:
        return False, "Historical Stage 2 Pilot24 QC gate JSON record corrupted"

    return True, (
        "Full484 case & solution ledger verified: 484 genomes, 12,580 cases "
        "(10,578 feasible ensembles, 1,936 fixed-medium infeasible, 56 verification-failed, 10 timeouts = "
        "66 explicit unknown/failed cases retained as NaN), 93,808 ledger rows, 91,806 accepted solutions "
        "(all residual/bounds/growth/cycle checks passed, 8 shard source hashes matched; raw-flux replay not claimed; "
        "historical Pilot24 JSON preserved)"
    )


def verify_gate5_four_witness_validation() -> tuple[bool, str]:
    """Gate 5: Four changed witness Slurm validation."""
    wit_json = EVIDENCE_DIR / "focused_release_validation" / "exp1_four_changed_witnesses_verification.json"
    if not wit_json.exists():
        return False, f"Missing four witness validation file: {wit_json}"
    with open(wit_json, "r", encoding="utf-8") as fh:
        data = json.load(fh)

    if not data.get("all_witnesses_verified", False):
        return False, "Four witness validation failed"

    n_wit = data.get("n_witnesses_evaluated", 0)
    if n_wit != 4:
        return False, f"Expected 4 witnesses, found {n_wit}"

    return True, f"Four witnesses revalidated on Slurm with zero drift across {n_wit} cases"


def verify_gate6_masking_controls_provenance() -> tuple[bool, str]:
    """Gate 6: Masking controls baseline, media provenance, and v4 post-hoc validation receipt."""
    base_file = EVIDENCE_DIR / "exp1_masking_controls_7hq_baselines.tsv"
    med_file = EVIDENCE_DIR / "exp1_masking_controls_media_provenance.json"
    v4_receipt_file = EVIDENCE_DIR / "exp1_masking_controls_7hq_v4_validation_receipt.json"
    v4_ctrl_file = EVIDENCE_DIR / "exp1_masking_controls_7hq_v4_validated.tsv"

    if not base_file.exists():
        return False, f"Missing masking baselines file: {base_file}"
    if not med_file.exists():
        return False, f"Missing masking media provenance file: {med_file}"
    if not v4_receipt_file.exists():
        return False, f"Missing v4 masking validation receipt: {v4_receipt_file}"
    if not v4_ctrl_file.exists():
        return False, f"Missing v4 validated masking controls table: {v4_ctrl_file}"

    base_df = pd.read_csv(base_file, sep="\t")
    actual_genomes = set(base_df["genome"].astype(str).unique())
    if actual_genomes != EXACT_7_HQ_EXTERNAL_IDS:
        return False, f"Masking baseline genomes mismatch: {actual_genomes} != {EXACT_7_HQ_EXTERNAL_IDS}"

    with open(med_file, "r", encoding="utf-8") as fh:
        med_data = json.load(fh)
    if set(med_data.keys()) != EXACT_7_HQ_EXTERNAL_IDS:
        return False, f"Masking media provenance genomes mismatch: {set(med_data.keys())}"

    with open(v4_receipt_file, "r", encoding="utf-8") as fh:
        v4_rec = json.load(fh)
    if v4_rec.get("computational_receipt_integrity") != "PASS" or not v4_rec.get("validation_passed", False):
        return False, "v4 masking controls computational receipt integrity did not pass"
    v4_df = pd.read_csv(v4_ctrl_file, sep="\t")
    if len(v4_df) != 63:
        return False, f"Expected 63 v4 masking control rows, found {len(v4_df)}"

    return True, (
        f"Masking controls verified: 7 external reference genomes (numeric >=90%/<=5%, formal MIMAG unknown), 63 v4 cases "
        f"(integrity={v4_rec.get('computational_receipt_integrity')}, "
        f"calibration={v4_rec.get('control_calibration_status')}, "
        f"raw_flux_status={v4_rec.get('raw_flux_vector_verification_status', 'UNKNOWN_NOT_PERSISTED_IN_V4_TABLES')})"
    )


def verify_gate7_slurm_and_stage3_stage4_receipts() -> tuple[bool, str]:
    """Gate 7: Slurm execution receipts, Stage 3 collector v4 receipt, Stage 4 v5 search receipt, and post-hoc 3,713-pair scientific disposition audit."""
    receipts_file = MANIFESTS_DIR / "slurm_execution_receipts.tsv"
    fetch_file = MANIFESTS_DIR / "institutional_fetch_instructions.txt"
    st3_receipt_file = EVIDENCE_DIR / "exp1_stage3_shard_collection_receipts.json"
    st4_receipt_file = EVIDENCE_DIR / "provenance" / "stage4_v5_execution_receipts.json"
    st4_crosswalk_file = EVIDENCE_DIR / "provenance" / "stage4_v5_reference_namespace_crosswalk.tsv"
    st4_summary_file = EVIDENCE_DIR / "exp1_stage4_v5_deficiency_verification_summary.tsv"
    st4_posthoc_receipt_file = EVIDENCE_DIR / "provenance" / "stage4_v5_posthoc_scientific_disposition_receipt.json"
    st4_posthoc_pairs_file = EVIDENCE_DIR / "exp1_stage4_v5_posthoc_unique_pair_dispositions.tsv"
    st4_posthoc_summary_file = EVIDENCE_DIR / "exp1_stage4_v5_posthoc_scientific_disposition_summary.tsv"
    st4_posthoc_asm_manifest_file = EVIDENCE_DIR / "provenance" / "stage4_v5_posthoc_47_sample_assembly_input_manifest.tsv"

    for req in (
        receipts_file,
        fetch_file,
        st3_receipt_file,
        st4_receipt_file,
        st4_crosswalk_file,
        st4_summary_file,
        st4_posthoc_receipt_file,
        st4_posthoc_pairs_file,
        st4_posthoc_summary_file,
        st4_posthoc_asm_manifest_file,
    ):
        if not req.exists():
            return False, f"Missing required receipt/evidence file: {req}"

    df = pd.read_csv(receipts_file, sep="\t")
    if len(df) < 8:
        return False, f"Expected at least 8 Slurm execution receipts, found {len(df)}"

    with open(st3_receipt_file, "r", encoding="utf-8") as fh:
        st3 = json.load(fh)
    if st3.get("n_total_collected_genomes") != EXPECTED_484_GENOME_COUNT or not st3.get("exact_genome_ids_matched"):
        return False, "Stage 3 collector v4 receipt failed 484-genome verification"

    with open(st4_receipt_file, "r", encoding="utf-8") as fh:
        st4 = json.load(fh)
    if st4.get("status") != "completed_verified_search" or st4.get("output_version") != "v5":
        return False, f"Stage 4 v5 receipt invalid status/version: {st4.get('status')} / {st4.get('output_version')}"
    ref_acc = st4.get("reference_namespace_resolution_accounting", {})
    n_searched_rxns = ref_acc.get("n_reactions_with_reference_sequences_searched", 0)
    n_unsearched_rxns = ref_acc.get("n_reactions_without_reference_sequences_unsearched", 0)
    n_searched_inst = ref_acc.get("n_instances_searched", 0)
    n_unsearched_inst = ref_acc.get("n_instances_unsearched_missing_reference", 0)
    if n_searched_rxns + n_unsearched_rxns != 63 or n_searched_inst + n_unsearched_inst != 13750:
        return False, "Stage 4 v5 reaction/instance accounting does not sum to 63 reactions and 13,750 instances"

    with open(st4_posthoc_receipt_file, "r", encoding="utf-8") as fh:
        posthoc = json.load(fh)
    denom = posthoc.get("denominator_policy", {})
    if denom.get("n_unique_genome_reaction_pairs") != 3713 or denom.get("n_within_pair_disposition_conflicts") != 0:
        return False, "Stage 4 v5 post-hoc unique pair denominator audit failed (expected 3,713 unique pairs, 0 conflicts)"
    nohit = posthoc.get("no_hit_unique_pair_audit", {})
    if nohit.get("n_no_hit_unique_pairs_total") != 59 or set(nohit.get("reactions", {}).keys()) != {"rxn10457", "rxn02322"}:
        return False, "Stage 4 v5 post-hoc 59 no-hit pair audit failed (expected rxn10457=47 and rxn02322=12)"
    asm_addendum = posthoc.get("sample_assembly_input_provenance_addendum", {})
    if asm_addendum.get("n_sample_assemblies_verified") != 47 or asm_addendum.get("original_execution_receipt_modified") is not False:
        return False, "Stage 4 v5 post-hoc 47-sample assembly input manifest verification failed"

    return True, (
        f"Slurm ({len(df)} jobs), Stage 3 (484 genomes, 8 shards), & Stage 4 v5 "
        f"({n_searched_rxns}/63 reactions, 3,713 unique pairs / 13,750 scenario instances, 47 assembly input hashes, 59 no-hit pairs qualified) verified"
    )


def is_excluded_path(rel_path_str: str) -> bool:
    """Return True if `rel_path_str` is an independently owned subpackage or mutable cache/bytecode file."""
    if rel_path_str == "manifests/sha256_manifest.tsv":
        return True
    if rel_path_str.startswith(EXCLUDED_SUBPACKAGE_PREFIXES):
        return True
    parts = Path(rel_path_str).parts
    if any(p in EXCLUDED_CACHE_PARTS for p in parts):
        return True
    if rel_path_str.endswith((".pyc", ".pyo")):
        return True
    return False


def audit_scoped_package_manifest(repo_root: Path, manifest_file: Path) -> tuple[bool, str]:
    """Verify that `manifests/sha256_manifest.tsv` covers all owned package files, excludes mutable
    caches and independently owned subpackages (`covariates/`, `interpretation/`, `counterfactual_review/`),
    and matches exact byte sizes and SHA-256 hashes on disk (`EXP1-FINAL-PACKAGE-ACCEPTANCE-REVIEW.txt` Item 6).
    """
    if not manifest_file.exists():
        return False, f"Missing package SHA-256 manifest: {manifest_file}"

    df = pd.read_csv(manifest_file, sep="\t")
    required_cols = {"sha256", "size_bytes", "file_path"}
    if not required_cols.issubset(df.columns):
        return False, f"Manifest {manifest_file} missing required columns {sorted(required_cols)}: found {list(df.columns)}"
    if df.empty or df["file_path"].isna().any() or df["file_path"].astype(str).duplicated().any():
        return False, "Manifest is empty or contains duplicate/null file_path entries"

    manifest_paths: set[str] = set()
    for _, row in df.iterrows():
        rel_str = str(row["file_path"])
        if is_excluded_path(rel_str):
            return False, f"Manifest includes forbidden cache or independently owned subpackage path: {rel_str}"
        abs_p = repo_root / rel_str
        if not abs_p.is_file():
            return False, f"Manifest file missing on disk: {rel_str}"
        actual_size = int(abs_p.stat().st_size)
        expected_size = int(row["size_bytes"])
        if actual_size != expected_size:
            return False, f"Byte size mismatch for {rel_str}: disk={actual_size}, manifest={expected_size}"
        actual_sha = compute_sha256(abs_p)
        expected_sha = str(row["sha256"]).strip()
        if actual_sha != expected_sha:
            return False, f"SHA-256 mismatch for {rel_str}: disk={actual_sha}, manifest={expected_sha}"
        manifest_paths.add(rel_str)

    # Discover all owned files on disk and ensure none are omitted from the manifest
    disk_owned_paths: set[str] = set()
    for p in sorted(repo_root.rglob("*")):
        if not p.is_file():
            continue
        rel_str = p.relative_to(repo_root).as_posix()
        if is_excluded_path(rel_str):
            continue
        disk_owned_paths.add(rel_str)

    untracked = sorted(disk_owned_paths - manifest_paths)
    if untracked:
        return False, f"Owned files on disk missing from sha256_manifest.tsv ({len(untracked)}): {untracked[:5]}"

    return True, (
        f"Scoped package SHA-256 manifest verified: {len(manifest_paths)} owned files matched "
        f"by path, byte size, and SHA-256 (mutable caches and independent subpackages excluded)"
    )


def verify_gate8_scoped_package_manifest() -> tuple[bool, str]:
    """Gate 8: Whole-package scoped SHA-256 manifest verification."""
    return audit_scoped_package_manifest(REPO_ROOT, MANIFESTS_DIR / "sha256_manifest.tsv")


def main() -> None:
    print("=" * 80)
    print("Experiment 1 Metabolic Modelling — Frozen Outputs & Reproducibility Gate Verifier")
    print(f"Package directory: {REPO_ROOT}")
    print("=" * 80)

    gates = [
        ("Gate 1: 484-Genome Cohort Universe & ID Uniqueness", verify_gate1_genome_inventory),
        ("Gate 2: Candidate Reaction Pool (2,271 Rows) & Fail-Closed Conflict Files", verify_gate2_candidate_pool),
        ("Gate 3: Biomass Template Audit Consistency", verify_gate3_biomass_template_audit),
        ("Gate 4: Full 484-Genome Case (12,580) & Per-Solution Ledger (93,808) Verification", verify_gate4_full484_case_and_solution_ledgers),
        ("Gate 5: Four Changed Witness Slurm Revalidation", verify_gate5_four_witness_validation),
        ("Gate 6: Masking Controls Baselines, Provenance & v4 Validation", verify_gate6_masking_controls_provenance),
        ("Gate 7: Slurm, Stage 3 Collector v4 & Stage 4 v5 Search Receipts", verify_gate7_slurm_and_stage3_stage4_receipts),
        ("Gate 8: Whole-Package Scoped SHA-256 & Byte-Size Manifest", verify_gate8_scoped_package_manifest),
    ]

    all_passed = True
    for name, func in gates:
        try:
            ok, msg = func()
            if ok:
                print(f"[PASS] {name}")
                print(f"       -> {msg}")
            else:
                print(f"[FAIL] {name}")
                print(f"       -> {msg}")
                all_passed = False
        except Exception as e:
            print(f"[ERROR] {name}")
            print(f"        -> Exception: {e}")
            all_passed = False

    print("=" * 80)
    if all_passed:
        print("[SUCCESS] All 8 frozen output and reproducibility gates PASSED cleanly.")
        sys.exit(0)
    else:
        print("[FAILURE] One or more gates failed verification.")
        sys.exit(1)


if __name__ == "__main__":
    main()
