from __future__ import annotations

import hashlib
import json
import math
import tarfile
from argparse import Namespace
from pathlib import Path

import pandas as pd
import pytest

from validate_masking_controls_v4_outputs import (
    EXPECTED_GENOMES, PASSED_CYCLE, ValidationError, compute_case_hash,
    compute_input_hash, sha256_file, sha256_text, stable_seed, validate_bundle,
    validate_ledger_row,
)
from stage3_collect_and_summarize_v4 import load_validated_masking_bundle
from run_masking_controls_7hq_v5_resumable import (
    build_empty_solver_ledger_row as build_v5_empty_ledger,
    build_solution_ledger_row as build_v5_solution_ledger,
    compute_case_config_hash as compute_v5_case_hash,
)

BIOMASS = "B2_perturbation_omit_7_non_f420_template_cofactors_keep_f420"
SCENARIO = "S4_oxic_desert_oa_aa"


def write_json(path: Path, value: dict) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True), encoding="utf-8")


def solution(condition: str, genome: str, fraction: float, replicate: int, status: str) -> dict:
    return {
        "condition_type": condition, "genome": genome, "cohort": "external_comparator",
        "canonical_order": "o__Test", "canonical_family": "f__Test", "scenario_id": SCENARIO,
        "biomass_recipe_id": BIOMASS, "mask_fraction": fraction, "replicate": replicate,
        "solution_index": 0, "validity_disposition": "accepted_feasible", "solver_status": status,
        "selected_reaction_ids": "R2" if condition == "masked" else "R1", "n_added_rxns": 1,
        "total_cost": 1.0, "raw_milp_objective": 1.0, "mip_gap": 0.0, "dual_bound": 1.0,
        "verified_growth": 0.02, "stoichiometric_max_residual": 1e-10,
        "max_bound_violation": 0.0, "residual_bounds_growth_verified": True,
        "cycle_qc_status": PASSED_CYCLE, "spurious_atp_cycle": False, "explicit_reason": "",
        "code_version": "v4_masking_controls_resumable",
    }


def placeholder(condition: str, genome: str, fraction: float, replicate: int, status: str) -> dict:
    row = solution(condition, genome, fraction, replicate, status)
    row.update({
        "validity_disposition": f"rejected_{status}", "selected_reaction_ids": "", "n_added_rxns": 0,
        "total_cost": math.nan, "raw_milp_objective": math.nan, "mip_gap": math.nan,
        "dual_bound": math.nan, "verified_growth": math.nan, "stoichiometric_max_residual": math.nan,
        "max_bound_violation": math.nan, "residual_bounds_growth_verified": False,
        "cycle_qc_status": "not_applicable_solver_failed", "spurious_atp_cycle": math.nan,
        "explicit_reason": f"Solver returned {status}",
    })
    return row


def complete_fixture(root: Path, valid_genomes: int = 7, valid_masks: bool = True) -> Namespace:
    results = root / "results"
    results.mkdir(parents=True)
    archives = root / "archives"
    archives.mkdir()
    code = Path(__file__).resolve().parents[1] / "code"
    engine = code / "exp1_sbml_highs_engine_v4.py"
    runner = code / "run_masking_controls_7hq_v4_resumable.py"
    helper = code / "stage2_pilot24_ensemble.py"
    pool = results / "pool.pkl"
    pool.write_bytes(b"fixture candidate pool")
    pool_sha = sha256_file(pool)
    engine_sha = sha256_file(engine)
    solver_config = {
        "milp_time_limit_s": 15.0, "mip_rel_gap": 0.05, "mip_feasibility_tolerance": 1e-9,
        "primal_feasibility_tolerance": 1e-9, "threads": 1, "biomass_recipe_id": BIOMASS,
        "scenario_id": SCENARIO,
    }
    solver_sha = sha256_text(json.dumps(solver_config, sort_keys=True))
    media: dict[str, dict] = {}
    input_hashes: dict[str, str] = {}
    for genome in EXPECTED_GENOMES:
        archive = archives / f"{genome}.tar.gz"
        xml = f'<sbml><model id="MODEL_{genome}"/></sbml>'.encode()
        medium_csv = b"compounds,name,maxflux\ncpd00001,H2O,10\n"
        with tarfile.open(archive, "w:gz") as tf:
            for name, payload in ((f"{genome}/{genome}-draft.xml", xml), (f"{genome}/{genome}-desert_medium.csv", medium_csv)):
                info = tarfile.TarInfo(name)
                info.size = len(payload)
                tf.addfile(info, __import__("io").BytesIO(payload))
        archive_sha = sha256_file(archive)
        compounds = {"cpd00001": ["H2O", 10.0], "cpd00007": ["O2", 10.0]}
        compounds_json = json.dumps(compounds, sort_keys=True)
        medium_hash = sha256_text(compounds_json)
        media[genome] = {"genome": genome, "archive_sha256": archive_sha, "medium_hash": medium_hash, "compounds_json": compounds_json}
        input_hashes[genome] = compute_input_hash(archive_sha, f"MODEL_{genome}", medium_hash, pool_sha, BIOMASS)

    baselines: list[dict] = []
    controls: list[dict] = []
    ledger: list[dict] = []
    base_records: dict[str, str] = {}
    mask_records: dict[str, str] = {}
    valid_set = set(EXPECTED_GENOMES[:valid_genomes])
    for genome in EXPECTED_GENOMES:
        base_valid = genome in valid_set
        status = "optimal" if base_valid else "solver_failed_status_4"
        case_hash = compute_case_hash(input_hashes[genome], solver_sha, engine_sha, f"{genome}_base")
        base_records[genome] = case_hash
        baselines.append({
            "condition_type": "baseline", "genome": genome, "cohort": "external_comparator",
            "canonical_order": "o__Test", "canonical_family": "f__Test", "completeness": 95.0,
            "contamination": 1.0, "scenario_id": SCENARIO, "biomass_recipe_id": BIOMASS,
            "mask_fraction": 0.0, "replicate": 0, "medium_hash": media[genome]["medium_hash"],
            "archive_sha256": media[genome]["archive_sha256"], "input_hash": input_hashes[genome],
            "case_config_hash": case_hash, "solver_status": status, "n_valid_solutions": int(base_valid),
            "n_gapfill_rxns": 1 if base_valid else math.nan, "gapfill_reaction_ids": "R1" if base_valid else "",
            "runtime_seconds": 1.0,
        })
        ledger.append(solution("baseline", genome, 0.0, 0, status) if base_valid else placeholder("baseline", genome, 0.0, 0, status))
        for fraction in (0.1, 0.2, 0.3):
            for replicate in (0, 1, 2):
                mask_valid = base_valid and valid_masks
                mask_status = "optimal" if mask_valid else "solver_failed_status_4"
                cid = f"{genome}_{fraction:.2f}_{replicate}"
                mask_hash = compute_case_hash(input_hashes[genome], solver_sha, engine_sha, cid)
                mask_records[cid] = mask_hash
                seed, seed_string = stable_seed(genome, fraction, replicate)
                usable = base_valid and mask_valid
                controls.append({
                    "condition_type": "masked", "genome": genome, "cohort": "external_comparator",
                    "canonical_order": "o__Test", "canonical_family": "f__Test", "completeness": 95.0,
                    "contamination": 1.0, "scenario_id": SCENARIO, "biomass_recipe_id": BIOMASS,
                    "mask_fraction": fraction, "replicate": replicate, "seed_val": seed, "seed_str": seed_string,
                    "medium_hash": media[genome]["medium_hash"], "input_hash": input_hashes[genome],
                    "case_config_hash": mask_hash, "n_total_genes": 10, "n_masked_genes": 1,
                    "masked_gene_ids": "G1", "n_removed_gpr_reactions": 1, "removed_rxn_ids": "R2",
                    "baseline_solver_status": status, "baseline_n_valid_solutions": int(base_valid),
                    "n_baseline_gapfill_rxns": 1 if base_valid else math.nan,
                    "masked_solver_status": mask_status, "masked_n_valid_solutions": int(mask_valid),
                    "n_masked_gapfill_rxns": 1 if mask_valid else math.nan, "is_evaluable": usable,
                    "evaluation_status": "evaluable_valid_ensembles" if usable else "non_evaluable_invalid_masked",
                    "n_newly_induced_gapfill_rxns": 1 if usable else math.nan,
                    "n_true_masked_rxns_recovered_in_gapfill": 1 if usable else math.nan,
                    "n_alternative_candidate_rxns_induced": 0 if usable else math.nan,
                    "masked_rxn_recovery_rate_under_medium": 1.0 if usable else math.nan,
                    "alternative_candidate_fraction_under_medium": 0.0 if usable else math.nan,
                    "runtime_seconds": 1.0,
                })
                ledger.append(solution("masked", genome, fraction, replicate, mask_status) if mask_valid else placeholder("masked", genome, fraction, replicate, mask_status))

    paths = {
        "controls_tsv": results / "controls.tsv", "baselines_tsv": results / "baselines.tsv",
        "ledger_tsv": results / "ledger.tsv", "media_json": results / "media.json",
        "checkpoint_json": results / "checkpoint.json", "original_report_json": results / "report.json",
    }
    pd.DataFrame(baselines).to_csv(paths["baselines_tsv"], sep="\t", index=False)
    pd.DataFrame(controls).to_csv(paths["controls_tsv"], sep="\t", index=False)
    pd.DataFrame(ledger).to_csv(paths["ledger_tsv"], sep="\t", index=False)
    write_json(paths["media_json"], media)
    write_json(paths["checkpoint_json"], {
        "baseline_records": base_records, "mask_records": mask_records, "engine_sha256": engine_sha,
        "solver_config_hash": solver_sha, "baselines_tsv_sha256": sha256_file(paths["baselines_tsv"]),
        "summary_tsv_sha256": sha256_file(paths["controls_tsv"]), "ledger_tsv_sha256": sha256_file(paths["ledger_tsv"]),
        "n_baselines_committed": 7, "n_mask_cases_committed": 63, "n_ledger_rows_committed": len(ledger),
    })
    write_json(paths["original_report_json"], {
        "validation_passed": True, "engine_sha256": engine_sha, "candidate_pool_sha256": pool_sha,
        "solver_config_hash": solver_sha,
    })
    return Namespace(
        **paths, engine_py=engine, runner_py=runner, helper_py=helper, pool_pkl=pool,
        gapseq_out_dir=archives, out_receipt=results / "validation.json", out_dir=results,
        milp_time_limit_s=15.0,
    )


def refresh_checkpoint(args: Namespace) -> None:
    checkpoint = json.loads(args.checkpoint_json.read_text())
    ledger = pd.read_csv(args.ledger_tsv, sep="\t")
    checkpoint.update({
        "baselines_tsv_sha256": sha256_file(args.baselines_tsv),
        "summary_tsv_sha256": sha256_file(args.controls_tsv),
        "ledger_tsv_sha256": sha256_file(args.ledger_tsv),
        "n_baselines_committed": len(pd.read_csv(args.baselines_tsv, sep="\t")),
        "n_mask_cases_committed": len(pd.read_csv(args.controls_tsv, sep="\t")),
        "n_ledger_rows_committed": len(ledger),
    })
    write_json(args.checkpoint_json, checkpoint)


def test_complete_and_partial_bundles_do_not_require_100_percent(tmp_path: Path) -> None:
    complete = complete_fixture(tmp_path / "complete")
    complete.out_receipt.parent.mkdir(exist_ok=True)
    receipt = validate_bundle(complete)
    assert receipt["computational_receipt_integrity"] == "PASS"
    assert receipt["n_biologically_usable_mask_cases"] == 63
    summary = load_validated_masking_bundle(complete.out_receipt, tmp_path / "collected")
    assert summary["control_calibration_status"] == "RESOLVED_ALL_CASES_USABLE"

    partial = complete_fixture(tmp_path / "partial", valid_genomes=2)
    partial.out_receipt.parent.mkdir(exist_ok=True)
    receipt = validate_bundle(partial)
    assert receipt["n_valid_matched_baselines"] == 2
    assert receipt["n_biologically_usable_mask_cases"] == 18
    assert receipt["control_calibration_status"] == "RESOLVED_PARTIAL_USABLE_CONTROLS"


def test_all_non_evaluable_is_complete_but_calibration_unresolved(tmp_path: Path) -> None:
    args = complete_fixture(tmp_path / "all_unknown", valid_genomes=0)
    args.out_receipt.parent.mkdir(exist_ok=True)
    receipt = validate_bundle(args)
    assert receipt["computational_receipt_integrity"] == "PASS"
    assert receipt["control_calibration_status"] == "UNRESOLVED_NO_USABLE_CONTROLS"
    assert receipt["statistical_claim_eligibility"] == "NOT_ELIGIBLE_NO_USABLE_CONTROLS"
    summary = load_validated_masking_bundle(args.out_receipt, tmp_path / "collected")
    assert summary["bounded_conditional_analysis_allowed"] is False
    assert summary["mean_recovery_rate"] is None
    corrected = pd.read_csv(receipt["derived_files"]["ledger"]["path"], sep="\t")
    assert corrected["n_added_rxns"].isna().all()


def test_partial_duplicate_orphan_and_hash_mutation_fail_closed(tmp_path: Path) -> None:
    args = complete_fixture(tmp_path / "bad")
    args.out_receipt.parent.mkdir(exist_ok=True)
    controls = pd.read_csv(args.controls_tsv, sep="\t")
    controls.head(1).to_csv(args.controls_tsv, sep="\t", index=False)
    refresh_checkpoint(args)
    with pytest.raises(ValidationError, match="Exact 63"):
        validate_bundle(args)

    args = complete_fixture(tmp_path / "duplicate")
    args.out_receipt.parent.mkdir(exist_ok=True)
    controls = pd.read_csv(args.controls_tsv, sep="\t")
    pd.concat([controls, controls.head(1)]).to_csv(args.controls_tsv, sep="\t", index=False)
    refresh_checkpoint(args)
    with pytest.raises(ValidationError, match="Exact 63"):
        validate_bundle(args)

    args = complete_fixture(tmp_path / "orphan")
    args.out_receipt.parent.mkdir(exist_ok=True)
    ledger = pd.read_csv(args.ledger_tsv, sep="\t")
    orphan = ledger.iloc[[0]].copy()
    orphan["genome"] = "ORPHAN"
    pd.concat([ledger, orphan]).to_csv(args.ledger_tsv, sep="\t", index=False)
    refresh_checkpoint(args)
    with pytest.raises(ValidationError, match="Orphaned ledger"):
        validate_bundle(args)

    args = complete_fixture(tmp_path / "mutation")
    args.out_receipt.parent.mkdir(exist_ok=True)
    with args.controls_tsv.open("a", encoding="utf-8") as fh:
        fh.write("# mutation\n")
    with pytest.raises(ValidationError, match="final commit marker"):
        validate_bundle(args)


def test_exact_boolean_types_and_actual_engine_solution_fields() -> None:
    authentic = solution("masked", "G1", 0.1, 0, "optimal")
    accepted, corrections = validate_ledger_row(pd.Series(authentic), "optimal")
    assert accepted and not corrections
    for field in ("spurious_atp_cycle", "residual_bounds_growth_verified"):
        for bad in (None, 0, 1, "False", float("nan")):
            corrupted = dict(authentic)
            corrupted[field] = bad
            with pytest.raises(ValidationError):
                validate_ledger_row(pd.Series(corrupted), "optimal")
    unknown = dict(authentic, science_qc_eligible=True, post_pruning_egc_has_spurious_atp_cycle=None)
    assert unknown["post_pruning_egc_has_spurious_atp_cycle"] is None


def test_future_v5_unknowns_actual_hash_and_schema() -> None:
    empty = build_v5_empty_ledger(
        "masked", "G1", "cohort", "order", "family", 0.1, 0, "solver_failed_status_4",
        input_hash="input", case_config_hash="case", engine_sha256="engine",
        runner_sha256="runner", helper_sha256="helper", solver_config_hash="solver",
    )
    assert math.isnan(empty["n_added_rxns"])
    assert math.isnan(empty["residual_bounds_growth_verified"])
    assert empty["case_config_hash"] == "case"

    engine_solution = {
        "added_rxns": ["R1"], "total_cost": 1.0, "raw_milp_objective": 1.0,
        "verified_fba_growth": 0.02, "fba_check_max_residual": 0.0,
        "fba_check_max_bound_viol": 0.0, "mip_gap": 0.0, "dual_bound": 1.0,
        "post_pruning_egc_has_spurious_atp_cycle": None,
        "cycle_qc_status": "unknown_cycle_check_timeout", "science_qc_eligible": True,
    }
    unknown = build_v5_solution_ledger(
        engine_solution, "masked", "G1", "cohort", "order", "family", 0.1, 0, 0, "optimal",
    )
    assert unknown["validity_disposition"] == "rejected_cycle_qc_ineligible"
    assert math.isnan(unknown["spurious_atp_cycle"])
    for field, bad in (("science_qc_eligible", 1), ("post_pruning_egc_has_spurious_atp_cycle", "False")):
        corrupted = dict(engine_solution)
        corrupted[field] = bad
        with pytest.raises(TypeError):
            build_v5_solution_ledger(corrupted, "masked", "G1", "cohort", "order", "family", 0.1, 0, 0, "optimal")

    first = compute_v5_case_hash("input", "solver", "engine-runner-helper-A", "case")
    second = compute_v5_case_hash("input", "solver", "engine-runner-helper-B", "case")
    assert first != second
