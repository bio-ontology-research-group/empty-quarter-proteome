#!/usr/bin/env python3
"""Independent, post-hoc validator for immutable exp1 masking-control v4 outputs.

This validates computational receipt integrity separately from biological control
usability. It never requires seven valid baselines, 63 valid masks, or any
recovery fraction. Missing/inconsistent infrastructure fails nonzero; a complete
bundle with zero usable cases exits zero with control calibration UNRESOLVED.
Original v4 files are never modified. Derived tables correct only explicit v4
unknown-as-zero/unknown-cycle placeholder defects and retain source hashes.
"""
from __future__ import annotations

import argparse
import ast
import hashlib
import importlib.util
import json
import math
import os
import sys
import tarfile
import time
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

SCHEMA_VERSION = "exp1_masking_v4_posthoc_validation_v1"
EXPECTED_GENOMES = (
    "MGYG000519815", "MGYG000520046", "MGYG000519748", "MGYG000510876",
    "MGYG000495096", "MGYG000526805", "MGYG000523080",
)
EXPECTED_FRACTIONS = (0.10, 0.20, 0.30)
EXPECTED_REPLICATES = (0, 1, 2)
VALID_SOLVER_STATUSES = {
    "optimal", "optimal_after_cycle_qc_filtering", "optimal_draft_no_gapfill_needed",
    "heuristic_incumbent_bounded", "post_pruning_cost_changed_not_proven_optimal",
}
PASSED_CYCLE = "passed_closed_uptake_egc_check"
FAILED_CYCLE = "failed_closed_uptake_atp_cycle"
NO_SOLUTION_CYCLE = "not_applicable_solver_failed"


class ValidationError(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValidationError(message)


def sha256_file(path: Path) -> str:
    require(path.is_file(), f"Required file missing: {path}")
    h = hashlib.sha256()
    with path.open("rb") as fh:
        while chunk := fh.read(1024 * 1024):
            h.update(chunk)
    return h.hexdigest()


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def atomic_json(data: Any, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    with tmp.open("w", encoding="utf-8") as fh:
        json.dump(data, fh, indent=2, sort_keys=True, allow_nan=False)
    os.replace(tmp, path)


def atomic_tsv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    df.to_csv(tmp, sep="\t", index=False)
    os.replace(tmp, path)


def load_json(path: Path) -> dict[str, Any]:
    require(path.is_file(), f"Required JSON missing: {path}")
    with path.open("r", encoding="utf-8") as fh:
        value = json.load(fh)
    require(isinstance(value, dict), f"JSON root must be an object: {path}")
    return value


def exact_bool(value: Any, field: str) -> bool:
    require(isinstance(value, (bool, np.bool_)), f"{field} must be boolean, got {type(value).__name__}: {value!r}")
    return bool(value)


def finite_float(value: Any, field: str, allow_missing: bool = False) -> float:
    if allow_missing and (value is None or (isinstance(value, float) and math.isnan(value)) or pd.isna(value)):
        return math.nan
    require(not isinstance(value, (bool, np.bool_)), f"{field} must be numeric, not boolean")
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise ValidationError(f"{field} is not numeric: {value!r}") from exc
    require(math.isfinite(result), f"{field} must be finite: {value!r}")
    return result


def integer_value(value: Any, field: str, minimum: int = 0) -> int:
    result = finite_float(value, field)
    require(result.is_integer() and result >= minimum, f"{field} must be an integer >= {minimum}: {value!r}")
    return int(result)


def missing(value: Any) -> bool:
    if value is None:
        return True
    try:
        return bool(pd.isna(value))
    except (TypeError, ValueError):
        return False


def split_ids(value: Any) -> list[str]:
    if missing(value):
        return []
    text = str(value)
    return [item for item in text.split(",") if item]


def mask_key(genome: str, fraction: float, replicate: int) -> str:
    return f"{genome}_{fraction:.2f}_{replicate}"


def ledger_key(row: pd.Series) -> tuple[str, float, int]:
    return str(row["genome"]), round(finite_float(row["mask_fraction"], "ledger mask_fraction"), 2), integer_value(row["replicate"], "ledger replicate")


def stable_seed(genome: str, fraction: float, replicate: int) -> tuple[int, str]:
    seed_string = f"exp1_masking_{genome}_{int(round(fraction * 100))}_{replicate}"
    seed = int(hashlib.sha256(seed_string.encode("utf-8")).hexdigest()[:8], 16)
    return seed, seed_string


def compute_input_hash(archive_sha: str, model_id: str, medium_hash: str, pool_sha: str, biomass: str) -> str:
    h = hashlib.sha256()
    for value in (archive_sha, model_id, medium_hash, pool_sha, biomass):
        h.update(value.encode("utf-8"))
    return h.hexdigest()


def compute_case_hash(input_hash: str, solver_hash: str, engine_sha: str, condition: str) -> str:
    h = hashlib.sha256()
    for value in (input_hash, solver_hash, engine_sha, condition):
        h.update(value.encode("utf-8"))
    return h.hexdigest()


def function_default(path: Path, function_name: str, argument: str) -> Any:
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    for node in ast.walk(tree):
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and node.name == function_name:
            args = node.args.args
            defaults = node.args.defaults
            mapping = {arg.arg: default for arg, default in zip(args[-len(defaults):], defaults)}
            require(argument in mapping, f"{path}: {function_name} lacks default {argument}")
            return ast.literal_eval(mapping[argument])
    raise ValidationError(f"{path}: function not found: {function_name}")


def verify_execution_correction(engine_py: Path, runner_py: Path, helper_py: Path, checkpoint: dict[str, Any], report: dict[str, Any], milp_time_limit_s: float) -> dict[str, Any]:
    """Verify the executed default statically; never rewrite the launch-time hash."""
    actual_gap = finite_float(function_default(engine_py, "enumerate_alternative_gapfill_ensemble", "mip_rel_gap_tol"), "executed mip_rel_gap_tol")
    tree = ast.parse(runner_py.read_text(encoding="utf-8"), filename=str(runner_py))
    ensemble_calls = [node for node in ast.walk(tree) if isinstance(node, ast.Call) and isinstance(node.func, ast.Name) and node.func.id == "enumerate_alternative_gapfill_ensemble"]
    require(len(ensemble_calls) == 2, f"Expected two ensemble calls in immutable v4 runner, found {len(ensemble_calls)}")
    for call in ensemble_calls:
        require(all(kw.arg != "mip_rel_gap_tol" for kw in call.keywords), "v4 runner unexpectedly overrides mip_rel_gap_tol")
        keywords = {kw.arg: kw.value for kw in call.keywords if kw.arg}
        for name, expected in (("n_Target_solutions", 8), ("cost_slack_frac", 0.25), ("min_growth", 0.01)):
            require(name in keywords and ast.literal_eval(keywords[name]) == expected, f"Unexpected executed {name}")

    recorded_config = {
        "milp_time_limit_s": float(milp_time_limit_s),
        "mip_rel_gap": 0.05,
        "mip_feasibility_tolerance": 1e-9,
        "primal_feasibility_tolerance": 1e-9,
        "threads": 1,
        "biomass_recipe_id": "B2_perturbation_omit_7_non_f420_template_cofactors_keep_f420",
        "scenario_id": "S4_oxic_desert_oa_aa",
    }
    recorded_hash = sha256_text(json.dumps(recorded_config, sort_keys=True))
    require(recorded_hash == checkpoint.get("solver_config_hash"), "Recorded v4 solver configuration does not reproduce checkpoint hash")
    require(report.get("solver_config_hash") == recorded_hash, "Original report and checkpoint solver hashes differ")
    actual_config = dict(recorded_config)
    actual_config["mip_rel_gap"] = actual_gap
    actual_config.update({"n_Target_solutions": 8, "cost_slack_frac": 0.25, "min_growth": 0.01})
    return {
        "schema_version": "exp1_v4_execution_correction_v1",
        "supersession": "POSTHOC_EXECUTION_CORRECTION_NOT_A_RETROACTIVE_LAUNCH_HASH",
        "original_report_sha256": None,
        "original_checkpoint_solver_config_sha256": recorded_hash,
        "original_recorded_mip_rel_gap": 0.05,
        "independently_verified_executed_mip_rel_gap_tol": actual_gap,
        "verification_method": "AST engine default plus both un-overridden runner call sites, bound to immutable source SHA-256",
        "executed_source_sha256": {
            "engine": sha256_file(engine_py), "runner": sha256_file(runner_py), "archive_weight_helper": sha256_file(helper_py),
        },
        "executed_engine_sha256_matches_checkpoint": sha256_file(engine_py) == checkpoint.get("engine_sha256"),
        "actual_execution_config": actual_config,
        "posthoc_actual_execution_config_sha256": sha256_text(json.dumps(actual_config, sort_keys=True)),
        "case_hash_interpretation": "Original per-case hashes intentionally use original checkpoint solver hash; corrected hash is descriptive only.",
    }


def archive_inputs(archive: Path, genome: str, engine_module: Any) -> tuple[str, str, dict[str, tuple[str, float]]]:
    require(archive.is_file(), f"Current source archive missing for {genome}: {archive}")
    with tarfile.open(archive, "r:gz") as tf:
        names = set(tf.getnames())
        draft_name = f"{genome}/{genome}-draft.xml"
        medium_name = f"{genome}/{genome}-desert_medium.csv"
        require(draft_name in names, f"Missing {draft_name} in {archive}")
        draft_bytes = tf.extractfile(draft_name).read()  # type: ignore[union-attr]
        medium_text = ""
        if medium_name in names:
            medium_text = tf.extractfile(medium_name).read().decode("utf-8", errors="replace")  # type: ignore[union-attr]
    root = ET.fromstring(draft_bytes)
    model = next((child for child in root if child.tag.rsplit("}", 1)[-1] == "model"), None)
    require(model is not None, f"No SBML model element for {genome}")
    model_id = model.attrib.get("id", "unknown_model")
    desert = engine_module.parse_gapseq_medium_csv(medium_text) if medium_text.strip() else {}
    if desert:
        medium = {**desert, "cpd00007": ("O2", 10.0)}
    else:
        medium = dict(engine_module.NUTRIENT_SCENARIOS[engine_module.SCEN_S4_OXIC_DESERT]["compounds"])
    return model_id, medium_text, medium


def canonical_medium(medium: dict[str, tuple[str, float]]) -> tuple[str, str]:
    normalized = {cid: [str(name), float(limit)] for cid, (name, limit) in sorted(medium.items())}
    text = json.dumps(normalized, sort_keys=True)
    return text, sha256_text(text)


def validate_ledger_row(row: pd.Series, expected_solver_status: str) -> tuple[bool, dict[str, Any]]:
    """Return (accepted, corrections); reject malformed/unknown states fail-closed."""
    corrections: dict[str, Any] = {}
    disposition = str(row.get("validity_disposition", ""))
    cycle = str(row.get("cycle_qc_status", ""))
    require(str(row.get("solver_status")) == expected_solver_status, "Ledger/summary solver status mismatch")
    accepted = disposition == "accepted_feasible"
    if accepted:
        require(cycle == PASSED_CYCLE, f"Accepted row lacks known cycle-pass label: {cycle}")
        require(exact_bool(row.get("spurious_atp_cycle"), "accepted spurious_atp_cycle") is False, "Accepted row has unknown/spurious cycle")
        require(exact_bool(row.get("residual_bounds_growth_verified"), "accepted residual_bounds_growth_verified") is True, "Accepted row lacks verified residual flag")
    elif cycle == NO_SOLUTION_CYCLE:
        require(not split_ids(row.get("selected_reaction_ids")), "No-solution placeholder unexpectedly has selected reactions")
        for column in ("total_cost", "raw_milp_objective", "verified_growth", "stoichiometric_max_residual", "max_bound_violation"):
            require(missing(row.get(column)), f"No-solution placeholder must have NA {column}")
        require(missing(row.get("spurious_atp_cycle")), "No-solution cycle state must be NA")
        # Executed v4 wrote zero here without a solution. Correct derived output to NA.
        if not missing(row.get("n_added_rxns")):
            require(finite_float(row.get("n_added_rxns"), "placeholder n_added_rxns") == 0.0, "Unexpected placeholder reaction count")
            corrections["n_added_rxns"] = math.nan
        if not missing(row.get("residual_bounds_growth_verified")):
            require(exact_bool(row.get("residual_bounds_growth_verified"), "placeholder residual flag") is False, "Unexpected placeholder residual flag")
            corrections["residual_bounds_growth_verified"] = math.nan
        return False, corrections
    else:
        require(disposition.startswith("rejected_"), f"Unknown ledger validity disposition: {disposition}")
        require(cycle in {PASSED_CYCLE, FAILED_CYCLE} or cycle.startswith("unknown_cycle_check_"), f"Unknown cycle QC label: {cycle}")
        if cycle.startswith("unknown_cycle_check_"):
            require(not accepted, "Unknown cycle cannot be accepted")
            # v4 coerced None to False. The unknown label takes precedence in derived output.
            if not missing(row.get("spurious_atp_cycle")):
                corrections["spurious_atp_cycle"] = math.nan
        else:
            exact_bool(row.get("spurious_atp_cycle"), "rejected spurious_atp_cycle")
        if disposition == "rejected_spurious_atp_cycle":
            require(cycle == FAILED_CYCLE and exact_bool(row.get("spurious_atp_cycle"), "spurious flag") is True, "Inconsistent spurious-cycle rejection")

    n_added = integer_value(row.get("n_added_rxns"), "n_added_rxns")
    require(len(split_ids(row.get("selected_reaction_ids"))) == n_added, "selected_reaction_ids/n_added_rxns mismatch")
    finite_float(row.get("total_cost"), "total_cost")
    finite_float(row.get("raw_milp_objective"), "raw_milp_objective")
    finite_float(row.get("mip_gap"), "mip_gap", allow_missing=True)
    finite_float(row.get("dual_bound"), "dual_bound", allow_missing=True)
    growth = finite_float(row.get("verified_growth"), "verified_growth")
    residual = finite_float(row.get("stoichiometric_max_residual"), "stoichiometric_max_residual")
    bound = finite_float(row.get("max_bound_violation"), "max_bound_violation")
    require(residual >= 0.0 and bound >= 0.0, "Negative residual/bound diagnostic")
    if accepted:
        require(growth >= 0.00999 and residual <= 1e-5 and bound <= 1e-5, "Accepted solution fails independent finite growth/residual/bound checks")
    return accepted, corrections


def validate_bundle(args: argparse.Namespace) -> dict[str, Any]:
    paths = {
        "controls": args.controls_tsv, "baselines": args.baselines_tsv, "ledger": args.ledger_tsv,
        "media": args.media_json, "checkpoint": args.checkpoint_json, "original_report": args.original_report_json,
        "engine": args.engine_py, "runner": args.runner_py, "helper": args.helper_py, "pool": args.pool_pkl,
    }
    source_sha = {name: sha256_file(path) for name, path in paths.items()}
    checkpoint = load_json(args.checkpoint_json)
    original_report = load_json(args.original_report_json)
    media = load_json(args.media_json)
    require(exact_bool(original_report.get("validation_passed"), "original report validation_passed") is True, "Original computational report did not pass")
    require(original_report.get("engine_sha256") == checkpoint.get("engine_sha256") == source_sha["engine"], "Executed engine hash mismatch")
    require(original_report.get("candidate_pool_sha256") == source_sha["pool"], "Current candidate pool differs from original report")
    require(checkpoint.get("baselines_tsv_sha256") == source_sha["baselines"], "Baseline file does not match final commit marker")
    require(checkpoint.get("summary_tsv_sha256") == source_sha["controls"], "Control file does not match final commit marker")
    require(checkpoint.get("ledger_tsv_sha256") == source_sha["ledger"], "Ledger file does not match final commit marker")

    correction = verify_execution_correction(args.engine_py, args.runner_py, args.helper_py, checkpoint, original_report, args.milp_time_limit_s)
    correction["original_report_sha256"] = source_sha["original_report"]
    require(correction["executed_engine_sha256_matches_checkpoint"] is True, "Execution correction engine mismatch")

    spec = importlib.util.spec_from_file_location("exp1_executed_v4_validation_source", args.engine_py)
    require(spec is not None and spec.loader is not None, "Cannot load executed engine for medium parsing")
    engine = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = engine
    spec.loader.exec_module(engine)
    biomass = engine.BIO_B2_OMIT_7_KEEP_F420
    scenario = engine.SCEN_S4_OXIC_DESERT

    baselines = pd.read_csv(args.baselines_tsv, sep="\t")
    controls = pd.read_csv(args.controls_tsv, sep="\t")
    ledger = pd.read_csv(args.ledger_tsv, sep="\t")
    require(len(baselines) == 7 and set(baselines["genome"].astype(str)) == set(EXPECTED_GENOMES), "Exact seven baseline rows/genomes required")
    require(not baselines["genome"].duplicated().any(), "Duplicate baseline genome")
    expected_masks = {(g, f, r) for g in EXPECTED_GENOMES for f in EXPECTED_FRACTIONS for r in EXPECTED_REPLICATES}
    actual_masks = {(str(row["genome"]), round(float(row["mask_fraction"]), 2), int(row["replicate"])) for _, row in controls.iterrows()}
    require(len(controls) == 63 and actual_masks == expected_masks and not controls.duplicated(["genome", "mask_fraction", "replicate"]).any(), "Exact 63 unique mask cases required")
    require(set(media) == set(EXPECTED_GENOMES), "Media provenance must contain exactly seven genomes")

    current: dict[str, dict[str, Any]] = {}
    for genome in EXPECTED_GENOMES:
        archive = args.gapseq_out_dir / f"{genome}.tar.gz"
        archive_sha = sha256_file(archive)
        model_id, _medium_text, medium = archive_inputs(archive, genome, engine)
        medium_json, medium_hash = canonical_medium(medium)
        provenance = media[genome]
        require(provenance.get("archive_sha256") == archive_sha, f"Current archive hash differs for {genome}")
        require(provenance.get("medium_hash") == medium_hash, f"Current archive-derived medium differs for {genome}")
        require(provenance.get("compounds_json") == medium_json, f"Persisted medium JSON differs from current archive for {genome}")
        input_hash = compute_input_hash(archive_sha, model_id, medium_hash, source_sha["pool"], biomass)
        current[genome] = {"archive_sha": archive_sha, "model_id": model_id, "medium_hash": medium_hash, "input_hash": input_hash}

    base_by_genome: dict[str, pd.Series] = {}
    expected_ledger_keys: set[tuple[str, float, int]] = set()
    for _, row in baselines.iterrows():
        genome = str(row["genome"])
        require(str(row["condition_type"]) == "baseline" and float(row["mask_fraction"]) == 0.0 and int(row["replicate"]) == 0, f"Bad baseline identity: {genome}")
        require(str(row["scenario_id"]) == scenario and str(row["biomass_recipe_id"]) == biomass, f"Baseline scenario/biomass mismatch: {genome}")
        require(str(row["archive_sha256"]) == current[genome]["archive_sha"] and str(row["medium_hash"]) == current[genome]["medium_hash"], f"Baseline input provenance mismatch: {genome}")
        require(str(row["input_hash"]) == current[genome]["input_hash"], f"Baseline input hash mismatch: {genome}")
        condition = f"{genome}_base"
        expected_case = compute_case_hash(current[genome]["input_hash"], str(checkpoint["solver_config_hash"]), source_sha["engine"], condition)
        require(str(row["case_config_hash"]) == expected_case and checkpoint.get("baseline_records", {}).get(genome) == expected_case, f"Baseline case hash mismatch: {genome}")
        base_by_genome[genome] = row
        expected_ledger_keys.add((genome, 0.0, 0))

    mask_by_key: dict[tuple[str, float, int], pd.Series] = {}
    for _, row in controls.iterrows():
        genome = str(row["genome"])
        fraction = round(finite_float(row["mask_fraction"], "mask_fraction"), 2)
        replicate = integer_value(row["replicate"], "replicate")
        key = (genome, fraction, replicate)
        require(str(row["condition_type"]) == "masked" and str(row["scenario_id"]) == scenario and str(row["biomass_recipe_id"]) == biomass, f"Mask identity mismatch: {key}")
        seed, seed_string = stable_seed(genome, fraction, replicate)
        require(integer_value(row["seed_val"], "seed_val") == seed and str(row["seed_str"]) == seed_string, f"Seed mismatch: {key}")
        require(str(row["medium_hash"]) == current[genome]["medium_hash"] and str(row["input_hash"]) == current[genome]["input_hash"], f"Mask input provenance mismatch: {key}")
        cid = mask_key(genome, fraction, replicate)
        expected_case = compute_case_hash(current[genome]["input_hash"], str(checkpoint["solver_config_hash"]), source_sha["engine"], cid)
        require(str(row["case_config_hash"]) == expected_case and checkpoint.get("mask_records", {}).get(cid) == expected_case, f"Mask case hash mismatch: {key}")
        mask_by_key[key] = row
        expected_ledger_keys.add(key)

    require(integer_value(checkpoint.get("n_baselines_committed"), "checkpoint baseline count") == 7, "Checkpoint baseline count mismatch")
    require(integer_value(checkpoint.get("n_mask_cases_committed"), "checkpoint mask count") == 63, "Checkpoint mask count mismatch")
    require(integer_value(checkpoint.get("n_ledger_rows_committed"), "checkpoint ledger count") == len(ledger), "Checkpoint ledger count mismatch")

    groups: dict[tuple[str, float, int], list[tuple[int, pd.Series]]] = {}
    for index, row in ledger.iterrows():
        key = ledger_key(row)
        require(key in expected_ledger_keys, f"Orphaned ledger row {index}: {key}")
        groups.setdefault(key, []).append((index, row))
    require(set(groups) == expected_ledger_keys, "Missing per-case ledger group")

    corrected_ledger = ledger.copy(deep=True)
    for column in ("n_added_rxns", "residual_bounds_growth_verified", "spurious_atp_cycle"):
        corrected_ledger[column] = corrected_ledger[column].astype(object)
    accepted_by_case: dict[tuple[str, float, int], list[pd.Series]] = {}
    unpersisted_rejected_cases: list[str] = []
    for key, members in groups.items():
        indices = [integer_value(row["solution_index"], "solution_index") for _, row in members]
        require(len(indices) == len(set(indices)), f"Duplicate solution_index in ledger case {key}")
        summary = base_by_genome[key[0]] if key[1] == 0.0 else mask_by_key[key]
        solver_status = str(summary["solver_status"] if key[1] == 0.0 else summary["masked_solver_status"])
        accepted_rows: list[pd.Series] = []
        for index, row in members:
            accepted, corrections = validate_ledger_row(row, solver_status)
            for column, value in corrections.items():
                corrected_ledger.at[index, column] = value
            if accepted:
                accepted_rows.append(row)
        accepted_by_case[key] = accepted_rows
        if not accepted_rows and str(summary.get("solver_status", summary.get("masked_solver_status", ""))).startswith("all_optimization_feasible_solutions_rejected_by_cycle_qc"):
            unpersisted_rejected_cases.append(mask_key(*key))

    baseline_valid: dict[str, bool] = {}
    for genome, row in base_by_genome.items():
        key = (genome, 0.0, 0)
        accepted = accepted_by_case[key]
        declared = integer_value(row["n_valid_solutions"], "baseline n_valid_solutions")
        require(declared == len(accepted), f"Baseline accepted ledger count mismatch: {genome}")
        valid = str(row["solver_status"]) in VALID_SOLVER_STATUSES and len(accepted) > 0
        require(valid == (str(row["solver_status"]) in VALID_SOLVER_STATUSES and declared > 0), f"Baseline status/count inconsistency: {genome}")
        union = {rid for solution in accepted for rid in split_ids(solution["selected_reaction_ids"])}
        require(set(split_ids(row["gapfill_reaction_ids"])) == union, f"Baseline reaction union mismatch: {genome}")
        if valid:
            require(integer_value(row["n_gapfill_rxns"], "baseline n_gapfill_rxns") == len(union), f"Baseline reaction count mismatch: {genome}")
        else:
            require(missing(row["n_gapfill_rxns"]), f"Invalid baseline reaction count must be NA: {genome}")
        baseline_valid[genome] = valid

    corrected_controls = controls.copy(deep=True)
    usable_keys: list[tuple[str, float, int]] = []
    for key, row in mask_by_key.items():
        genome = key[0]
        accepted = accepted_by_case[key]
        declared = integer_value(row["masked_n_valid_solutions"], "masked_n_valid_solutions")
        require(declared == len(accepted), f"Masked accepted ledger count mismatch: {key}")
        masked_valid = str(row["masked_solver_status"]) in VALID_SOLVER_STATUSES and len(accepted) > 0
        require(masked_valid == (str(row["masked_solver_status"]) in VALID_SOLVER_STATUSES and declared > 0), f"Masked status/count inconsistency: {key}")
        require(str(row["baseline_solver_status"]) == str(base_by_genome[genome]["solver_status"]), f"Baseline status linkage mismatch: {key}")
        require(integer_value(row["baseline_n_valid_solutions"], "baseline_n_valid_solutions") == len(accepted_by_case[(genome, 0.0, 0)]), f"Baseline count linkage mismatch: {key}")
        union = {rid for solution in accepted for rid in split_ids(solution["selected_reaction_ids"])}
        if masked_valid:
            require(integer_value(row["n_masked_gapfill_rxns"], "n_masked_gapfill_rxns") == len(union), f"Masked reaction union/count mismatch: {key}")
        else:
            require(missing(row["n_masked_gapfill_rxns"]), f"Invalid masked reaction count must be NA: {key}")
        usable = baseline_valid[genome] and masked_valid
        original_evaluable = exact_bool(row["is_evaluable"], f"is_evaluable {key}")
        require(original_evaluable == usable, f"Summary usability contradicts validated solution evidence: {key}")
        metric_columns = ("n_newly_induced_gapfill_rxns", "n_true_masked_rxns_recovered_in_gapfill", "n_alternative_candidate_rxns_induced", "masked_rxn_recovery_rate_under_medium", "alternative_candidate_fraction_under_medium")
        if usable:
            newly = integer_value(row["n_newly_induced_gapfill_rxns"], "n_newly_induced")
            recovered = integer_value(row["n_true_masked_rxns_recovered_in_gapfill"], "n_recovered")
            alternative = integer_value(row["n_alternative_candidate_rxns_induced"], "n_alternative")
            require(newly == recovered + alternative, f"Recovery arithmetic mismatch: {key}")
            removed = split_ids(row["removed_rxn_ids"])
            require(integer_value(row["n_removed_gpr_reactions"], "n_removed_gpr_reactions") == len(removed), f"Removed reaction count mismatch: {key}")
            if removed:
                require(math.isclose(finite_float(row["masked_rxn_recovery_rate_under_medium"], "recovery"), recovered / len(removed), rel_tol=1e-10, abs_tol=1e-12), f"Recovery fraction mismatch: {key}")
            else:
                require(missing(row["masked_rxn_recovery_rate_under_medium"]), f"Zero-denominator recovery must be NA: {key}")
            expected_alt = alternative / newly if newly else 0.0
            require(math.isclose(finite_float(row["alternative_candidate_fraction_under_medium"], "alternative fraction"), expected_alt, rel_tol=1e-10, abs_tol=1e-12), f"Alternative fraction mismatch: {key}")
            require(str(row["evaluation_status"]) == "evaluable_valid_ensembles", f"Usable case has wrong evaluation status: {key}")
            usable_keys.append(key)
        else:
            for column in metric_columns:
                require(missing(row[column]), f"Non-usable case metric must be NA: {key} {column}")
            require(str(row["evaluation_status"]).startswith("non_evaluable_"), f"Non-usable case has wrong evaluation status: {key}")

    args.out_dir.mkdir(parents=True, exist_ok=True)
    derived_controls_path = args.out_dir / "exp1_masking_controls_7hq_v4_validated.tsv"
    derived_baselines_path = args.out_dir / "exp1_masking_controls_7hq_v4_validated_baselines.tsv"
    derived_ledger_path = args.out_dir / "exp1_masking_controls_7hq_v4_validated_per_solution_ledger.tsv"
    correction_path = args.out_dir / "exp1_masking_controls_7hq_v4_execution_correction.json"
    # Preserve values; the independent validation above establishes usability. Corrections are only explicit unknown placeholders.
    atomic_tsv(corrected_controls, derived_controls_path)
    atomic_tsv(baselines, derived_baselines_path)
    atomic_tsv(corrected_ledger, derived_ledger_path)
    atomic_json(correction, correction_path)

    n_usable = len(usable_keys)
    n_valid_baselines = sum(baseline_valid.values())
    if n_usable:
        calibration = "RESOLVED_ALL_CASES_USABLE" if n_usable == 63 else "RESOLVED_PARTIAL_USABLE_CONTROLS"
        claim_eligibility = "ELIGIBLE_BOUNDED_CONDITIONAL_WITH_REPORTED_DENOMINATOR"
    else:
        calibration = "UNRESOLVED_NO_USABLE_CONTROLS"
        claim_eligibility = "NOT_ELIGIBLE_NO_USABLE_CONTROLS"
    usable_set = set(usable_keys)
    evaluable = controls.apply(lambda row: (str(row["genome"]), round(float(row["mask_fraction"]), 2), int(row["replicate"])) in usable_set, axis=1)
    recovery = controls.loc[evaluable, "masked_rxn_recovery_rate_under_medium"].dropna()
    alternative = controls.loc[evaluable, "alternative_candidate_fraction_under_medium"].dropna()
    derived_files = {
        "controls": derived_controls_path, "baselines": derived_baselines_path,
        "ledger": derived_ledger_path, "execution_correction": correction_path,
    }
    receipt = {
        "schema_version": SCHEMA_VERSION,
        "validated_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "computational_receipt_integrity": "PASS",
        "validation_passed": True,
        "control_calibration_status": calibration,
        "statistical_claim_eligibility": claim_eligibility,
        "n_declared_baselines": 7,
        "n_valid_matched_baselines": n_valid_baselines,
        "n_declared_mask_cases": 63,
        "n_biologically_usable_mask_cases": n_usable,
        "n_unusable_or_unknown_mask_cases": 63 - n_usable,
        "usable_case_denominator": 63,
        "mean_recovery_over_usable_non_NA": float(recovery.mean()) if len(recovery) else None,
        "n_recovery_values": int(len(recovery)),
        "mean_alternative_fraction_over_usable_non_NA": float(alternative.mean()) if len(alternative) else None,
        "n_alternative_fraction_values": int(len(alternative)),
        "zero_recovery_is_not_required_for_usability": True,
        "no_aggregate_percent_success_gate": True,
        "raw_flux_vectors_persisted": False,
        "raw_flux_vector_status": "UNKNOWN_NOT_PERSISTED_IN_V4_TABLES",
        "science_qc_eligible_source_field_persisted": False,
        "accepted_solution_rule": "accepted_disposition + exact booleans + known cycle-pass label + finite independently checked growth/residual/bounds",
        "unpersisted_rejected_solution_measurements_status": "UNKNOWN_NOT_PERSISTED_BY_EXECUTED_V4" if unpersisted_rejected_cases else "NOT_APPLICABLE_NO_SUCH_CASE_IDENTIFIED",
        "unpersisted_rejected_solution_cases": sorted(unpersisted_rejected_cases),
        "source_files": {name: {"path": str(path), "sha256": source_sha[name]} for name, path in paths.items()},
        "current_input_provenance": current,
        "original_solver_config_sha256": checkpoint["solver_config_hash"],
        "posthoc_execution_correction_sha256": sha256_file(correction_path),
        "derived_files": {name: {"path": str(path), "sha256": sha256_file(path)} for name, path in derived_files.items()},
        "original_files_modified": False,
        "interpretation_boundary": "Complete receipt integrity is not equivalent to biological control usability. Partial usable controls support only bounded conditional analysis with denominators; zero usable controls leaves calibration unresolved.",
    }
    atomic_json(receipt, args.out_receipt)
    return receipt


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("controls-tsv", "baselines-tsv", "ledger-tsv", "media-json", "checkpoint-json", "original-report-json", "engine-py", "runner-py", "helper-py", "pool-pkl"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--gapseq-out-dir", type=Path, required=True)
    parser.add_argument("--out-receipt", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--milp-time-limit-s", type=float, default=15.0)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        receipt = validate_bundle(args)
    except Exception as exc:
        failure = {
            "schema_version": SCHEMA_VERSION,
            "validated_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
            "computational_receipt_integrity": "FAIL",
            "validation_passed": False,
            "control_calibration_status": "UNKNOWN_INFRASTRUCTURE_FAILURE",
            "statistical_claim_eligibility": "NOT_ELIGIBLE_INFRASTRUCTURE_FAILURE",
            "error_type": type(exc).__name__,
            "error": str(exc),
            "original_files_modified": False,
        }
        try:
            atomic_json(failure, args.out_receipt)
        finally:
            print(f"[V4 POSTHOC VALIDATION FAILED] {type(exc).__name__}: {exc}", file=sys.stderr)
        return 2
    print(json.dumps({
        "computational_receipt_integrity": receipt["computational_receipt_integrity"],
        "control_calibration_status": receipt["control_calibration_status"],
        "n_valid_matched_baselines": receipt["n_valid_matched_baselines"],
        "n_biologically_usable_mask_cases": receipt["n_biologically_usable_mask_cases"],
        "receipt": str(args.out_receipt),
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
