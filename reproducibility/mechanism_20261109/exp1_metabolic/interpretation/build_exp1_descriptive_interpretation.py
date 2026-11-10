#!/usr/bin/env python3
"""Reproducible final descriptive interpretation of the completed 484 metabolic models.

Implements the fixed scientific design in:
- mechanism-20261109/EXP1-DESCRIPTIVE-INTERPRETATION-IMPLEMENTATION.txt
- mechanism-20261109/EXP1-FULL-COHORT-INTERPRETATION-DESIGN.txt
- mechanism-20261109/EXP1-INTERPRETATION-FIGURE-REVIEW.txt

Scope & invariants:
- Owns ONLY exp1_metabolic/interpretation/ and its review-20261109/interpretation/ mirror.
- Does not rerun metabolic models or modify existing engine/collector/deficiency files.
- Never uses /tmp; all temporary/matplotlib cache directories live inside the output tree and are cleaned up.
- Portable relative default paths anchored to the package directory; explicit CLI paths fail closed if missing.
- Unknown / timeout / verification-failed cases are separated from fixed-medium infeasibility; neither is coerced to zero.
- MGnify unknown-method quality metrics are never pooled into CheckM2 method-family strata.
- The 7 masking control genomes are described as high-completeness, low-contamination reference genomes
  (numeric >=90% completeness, <=5% contamination; formal MIMAG status unknown for MGnify).
- 63 masking evaluations are summarized across the 7 reference genomes (9 conditions/replicates per genome:
  3 mask fractions x 3 random seeds), never as 63 independent genomes.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import shutil
import statistics
from collections import Counter
from pathlib import Path
from typing import Any, Iterable

import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
PACKAGE_DIR = SCRIPT_DIR.parent

# Portable relative default paths within the reproducibility package / review mirror
DEFAULT_COMBINED_SUMMARY = SCRIPT_DIR / "exp1_full484_combined_summary_12580.tsv"
DEFAULT_REVIEW_SHARDS_CANDIDATES = (
    SCRIPT_DIR.parent / "full484" / "shards",
    Path("/home/leechuck/Documents/papers/empty-quarter/mechanism-20261109/review-20261109/full484/shards"),
)
DEFAULT_COVARIATE_SUPPLEMENT = (
    PACKAGE_DIR / "covariates" / "exp1_484_covariate_metadata_supplement.tsv"
)
DEFAULT_INVENTORY = (
    PACKAGE_DIR / "compact_evidence" / "exp1_model_inventory_and_qc_484.tsv"
    if (PACKAGE_DIR / "compact_evidence" / "exp1_model_inventory_and_qc_484.tsv").is_file()
    else PACKAGE_DIR / "exp1_model_inventory_and_qc_484.tsv"
)
DEFAULT_MASKING_CASES = SCRIPT_DIR / "inputs" / "exp1_masking_controls_7hq_v4.tsv"
DEFAULT_MASKING_BASELINES = SCRIPT_DIR / "inputs" / "exp1_masking_controls_7hq_v4_baselines.tsv"
DEFAULT_ROOT_PAIRED_AUDIT_JSON = (
    SCRIPT_DIR / "inputs" / "exp1_full484_paired_biomass_sensitivity.json"
)
DEFAULT_OUTDIR = SCRIPT_DIR

EXPECTED_SOURCE_SHA256 = {
    "exp1_shard_00_summary.tsv": "5182219833aac046f33beccfc7483e355e6073c8b6cd4c815cfc0cd5e3ac1b3e",
    "exp1_shard_01_summary.tsv": "b8bc464b2c9c72d83c8a5ecaa7b735614a7e15cada379c738aa845d442000740",
    "exp1_shard_02_summary.tsv": "15886a3eed4ddd7790626a54b0b08b626c9c2ed59bc5e36f03ee678faacb00c3",
    "exp1_shard_03_summary.tsv": "3b28adffb53c6ca35380178f86e873425f420a3ef9f7747b4bfe18e33032904e",
    "exp1_shard_04_summary.tsv": "94372f5d9a8899557cfe06af5ee251e4a1f1337c6beba08fa19a624342c9d1bd",
    "exp1_shard_05_summary.tsv": "820a9572bea37e748aa837407199ee302c755891a73d25fc528916823aa21c8e",
    "exp1_shard_06_summary.tsv": "3ee7456574b6fe5dc9f0d413f29f175e1607aef7d32542d92408c987d768cc21",
    "exp1_shard_07_summary.tsv": "0e30d8d6b2503b4c944913c2ff17cad6bf202929c3e0cd16a5e5be8e727e438b",
    "exp1_484_covariate_metadata_supplement.tsv": "c3492b38be9978d3b81fd269ae81ab7a0baa41e7d8554a2772d312e5083a016c",
    "exp1_model_inventory_and_qc_484.tsv": "044014df6aa1f7bf9b5d4fcaea3a2daea0ef62f3dfcb3c79ef08e54b0347eb2b",
    "exp1_masking_controls_7hq_v4.tsv": "8268f1f5372c0f60190e5b21b6cb25b8b72e0f342fe65bf52f70b8ce5013a53f",
    "exp1_masking_controls_7hq_v4_baselines.tsv": "80fbbe6761b2c1c4bd0d5a2093cb70c0770c55b84f803c9fca2f5b8e132ab528",
}

B0 = "B0_methanosarcina_template_baseline"
B1 = "B1_perturbation_omit_8_template_cofactors"
B2 = "B2_perturbation_omit_7_non_f420_template_cofactors_keep_f420"
B3 = "B3_perturbation_omit_template_cofactors_and_surface_glycolipids"
B4 = "B4_perturbation_core_amino_nucleotide_energy_skeleton"

S0 = "S0_gapseq_native_predicted"
S1 = "S1_anoxic_minimal_h2_co2_nh4"
S2 = "S2_oxic_minimal_h2_co2_nh4"
S3 = "S3_anoxic_desert_oa_aa"
S4 = "S4_oxic_desert_oa_aa"
S5 = "S5_anoxic_defined_aa_precursors"
S6 = "S6_oxic_defined_aa_precursors"

EXPECTED_24_CORE_COMBINATIONS: list[tuple[str, str]] = [
    (S1, B0),
    (S1, B2),
    (S2, B0),
    (S2, B2),
    (S3, B0),
    (S3, B1),
    (S3, B2),
    (S3, B3),
    (S3, B4),
    (S4, B0),
    (S4, B1),
    (S4, B2),
    (S4, B3),
    (S4, B4),
    (S5, B0),
    (S5, B1),
    (S5, B2),
    (S5, B3),
    (S5, B4),
    (S6, B0),
    (S6, B1),
    (S6, B2),
    (S6, B3),
    (S6, B4),
]

EXPECTED_OPTIONAL_S0_COMBINATIONS: list[tuple[str, str]] = [
    (S0, B0),
    (S0, B2),
]

VERIFIED_EMPTY_S0_GENOMES = frozenset({"GCA_016180755.1", "S22PRr1__self_SemiBin_29"})
SHARED_ORDERS = ("o__JACQPN01", "o__JACPAO01", "o__SW-10-69-26")

COMPLETENESS_BINS = ("[50,70)", "[70,90)", "[90,101]", "missing", "out_of_range")
CONTAMINATION_BINS = ("[0,5)", "[5,10)", "boundary_10", "missing", "out_of_range")

VALID_SOLVER_STATUSES = frozenset({"heuristic_incumbent_bounded", "optimal"})
INFEASIBLE_SOLVER_STATUSES = frozenset({"infeasible_under_fixed_medium"})
NUMERICAL_UNKNOWN_STATUSES = frozenset({"fba_verification_failed", "timeout_no_incumbent"})
ALL_KNOWN_SOLVER_STATUSES = (
    VALID_SOLVER_STATUSES | INFEASIBLE_SOLVER_STATUSES | NUMERICAL_UNKNOWN_STATUSES
)

RESIDUAL_TOL = 1e-6
BOUND_VIOL_TOL = 1e-6
MIN_VERIFIED_GROWTH = 0.01 - 1e-9


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        raise FileNotFoundError(f"Required TSV file does not exist: {path}")
    with path.open("r", encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path: Path, fieldnames: list[str], rows: Iterable[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: ("" if row.get(k) is None else row.get(k)) for k in fieldnames})


def classify_completeness_bin(val_str: str | float | None) -> str:
    """Classify completeness into fixed bins [50,70), [70,90), [90,101], missing, or out_of_range."""
    if val_str is None or str(val_str).strip() == "":
        return "missing"
    try:
        v = float(val_str)
    except ValueError:
        return "missing"
    if not math.isfinite(v):
        return "missing"
    if 50.0 <= v < 70.0:
        return "[50,70)"
    if 70.0 <= v < 90.0:
        return "[70,90)"
    if 90.0 <= v <= 101.0:
        return "[90,101]"
    return "out_of_range"


def classify_contamination_bin(val_str: str | float | None) -> str:
    """Classify contamination into [0,5), [5,10), boundary_10, missing, or out_of_range."""
    if val_str is None or str(val_str).strip() == "":
        return "missing"
    try:
        v = float(val_str)
    except ValueError:
        return "missing"
    if not math.isfinite(v):
        return "missing"
    if 0.0 <= v < 5.0:
        return "[0,5)"
    if 5.0 <= v < 10.0:
        return "[5,10)"
    if v == 10.0:
        return "boundary_10"
    return "out_of_range"


def parse_finite_float(val_str: str | float | None) -> float | None:
    if val_str is None or str(val_str).strip() == "":
        return None
    try:
        v = float(val_str)
    except ValueError:
        return None
    return v if math.isfinite(v) else None


def is_valid_summary_endpoint(row: dict[str, str], *, strict_diagnostics: bool = True) -> bool:
    """Verify whether a summary row represents a valid feasible model ensemble endpoint."""
    status = (row.get("solver_status") or "").strip()
    if status not in ALL_KNOWN_SOLVER_STATUSES:
        raise ValueError(f"Unknown solver_status={status!r} for genome={row.get('genome')}")

    n_valid_f = parse_finite_float(row.get("n_distinct_valid_solutions"))
    n_valid = int(n_valid_f) if n_valid_f is not None else 0

    if status in INFEASIBLE_SOLVER_STATUSES | NUMERICAL_UNKNOWN_STATUSES:
        if n_valid > 0:
            raise ValueError(
                f"Inconsistent summary row: solver_status={status} but n_distinct_valid_solutions={n_valid} "
                f"for ({row.get('genome')}, {row.get('scenario_id')}, {row.get('biomass_recipe_id')})"
            )
        return False

    if n_valid <= 0:
        if strict_diagnostics:
            raise ValueError(
                f"Inconsistent summary row: solver_status={status} has n_distinct_valid_solutions={n_valid} "
                f"for ({row.get('genome')}, {row.get('scenario_id')}, {row.get('biomass_recipe_id')})"
            )
        return False

    min_added = parse_finite_float(row.get("min_added_rxns_count"))
    max_res = parse_finite_float(row.get("max_solution_residual"))
    max_bv = parse_finite_float(row.get("max_solution_bound_viol"))
    min_gr = parse_finite_float(row.get("min_solution_verified_growth"))

    diagnostics_ok = (
        min_added is not None
        and min_added >= 0.0
        and max_res is not None
        and max_res <= RESIDUAL_TOL
        and max_bv is not None
        and max_bv <= BOUND_VIOL_TOL
        and min_gr is not None
        and min_gr >= MIN_VERIFIED_GROWTH
    )
    if not diagnostics_ok:
        if strict_diagnostics:
            raise ValueError(
                f"Failed verified diagnostics on feasible summary row "
                f"({row.get('genome')}, {row.get('scenario_id')}, {row.get('biomass_recipe_id')}): "
                f"min_added={min_added}, max_res={max_res}, max_bv={max_bv}, min_gr={min_gr}"
            )
        return False
    return True


def quantile_summary(values: list[float]) -> dict[str, float | None]:
    if not values:
        return {
            "median": None,
            "q25": None,
            "q75": None,
            "iqr": None,
            "min": None,
            "max": None,
            "mean": None,
        }
    arr = np.asarray(values, dtype=float)
    q25, med, q75 = [float(x) for x in np.percentile(arr, [25, 50, 75])]
    return {
        "median": med,
        "q25": q25,
        "q75": q75,
        "iqr": float(q75 - q25),
        "min": float(np.min(arr)),
        "max": float(np.max(arr)),
        "mean": float(np.mean(arr)),
    }


def fmt_num(val: float | int | None, digits: int = 4) -> str:
    if val is None:
        return ""
    if isinstance(val, int):
        return str(val)
    if not math.isfinite(val):
        return ""
    if float(val).is_integer() and digits <= 2:
        return f"{val:.1f}"
    return f"{val:.{digits}f}".rstrip("0").rstrip(".") if digits > 2 else f"{val:.{digits}f}"


def load_summary_rows(
    *,
    shards_dir: Path | None,
    combined_summary: Path | None,
    strict_source_sha256: bool,
) -> tuple[list[dict[str, str]], dict[str, str]]:
    """Load summary rows from either explicit combined TSV or 8 shard summary TSVs.

    Never falls back if an explicit path is missing or invalid.
    """
    input_hashes: dict[str, str] = {}
    if combined_summary is not None:
        if not combined_summary.is_file():
            raise FileNotFoundError(f"Explicit --combined-summary file not found: {combined_summary}")
        input_hashes[combined_summary.name] = sha256_file(combined_summary)
        rows = read_tsv(combined_summary)
        if not rows:
            raise ValueError(f"Combined summary TSV is empty: {combined_summary}")
        return rows, input_hashes

    if shards_dir is None:
        raise ValueError("Either --shards-dir or --combined-summary must be provided")
    if not shards_dir.is_dir():
        raise FileNotFoundError(f"Explicit --shards-dir directory not found: {shards_dir}")

    shard_paths = [shards_dir / f"exp1_shard_0{i}_summary.tsv" for i in range(8)]
    rows: list[dict[str, str]] = []
    for sp in shard_paths:
        if not sp.is_file():
            raise FileNotFoundError(f"Missing required shard summary file: {sp}")
        sh = sha256_file(sp)
        if strict_source_sha256 and sp.name in EXPECTED_SOURCE_SHA256:
            exp_sh = EXPECTED_SOURCE_SHA256[sp.name]
            if sh != exp_sh:
                raise ValueError(f"SHA-256 mismatch for {sp.name}: got {sh}, expected {exp_sh}")
        input_hashes[sp.name] = sh
        rows.extend(read_tsv(sp))
    return rows, input_hashes


def validate_universe_and_join_covariates(
    summary_rows: list[dict[str, str]],
    cov_rows: list[dict[str, str]],
    inv_rows: list[dict[str, str]] | None = None,
    *,
    expected_Empty_s0_genomes: frozenset[str] = VERIFIED_EMPTY_S0_GENOMES,
) -> tuple[dict[tuple[str, str, str], dict[str, str]], dict[str, dict[str, str]]]:
    """Validate exact case universe, no duplicates, no missing core cases, and exact covariate join."""
    if not cov_rows:
        raise ValueError("Covariate supplement is empty")

    cov_by_genome: dict[str, dict[str, str]] = {}
    for r in cov_rows:
        g = (r.get("genome") or "").strip()
        if not g:
            raise ValueError("Empty genome ID in covariate supplement")
        if g in cov_by_genome:
            raise ValueError(f"Duplicate genome in covariate supplement: {g}")
        q_method = (r.get("reported_quality_method") or "").strip()
        if r.get("source") == "MGnify" or "unknown" in q_method.lower():
            if (r.get("checkm2_method_family_compatible") or "").strip() == "True":
                raise ValueError(
                    f"Unknown quality method genome {g} ({q_method}) illegally marked checkm2_method_family_compatible=True"
                )
            if (r.get("method_family_checkm2_completeness") or "").strip() != "" or (
                r.get("method_family_checkm2_contamination") or ""
            ).strip() != "":
                raise ValueError(
                    f"Unknown quality method genome {g} has non-empty method_family_checkm2 quality columns"
                )
        cov_by_genome[g] = r

    if inv_rows is not None:
        inv_by_genome = {r["genome"].strip(): r for r in inv_rows}
        if set(inv_by_genome.keys()) != set(cov_by_genome.keys()):
            raise ValueError("Inventory genome set does not match covariate supplement genome set")
        for g, c_row in cov_by_genome.items():
            i_row = inv_by_genome[g]
            for col in ("cohort", "source", "canonical_order", "canonical_family", "canonical_genus"):
                if (c_row.get(col) or "").strip() != (i_row.get(col) or "").strip():
                    raise ValueError(f"Covariate vs inventory mismatch on {col} for {g}")

    summary_by_key: dict[tuple[str, str, str], dict[str, str]] = {}
    cases_by_genome: dict[str, set[tuple[str, str]]] = {g: set() for g in cov_by_genome}

    for r in summary_rows:
        g = (r.get("genome") or "").strip()
        scen = (r.get("scenario_id") or "").strip()
        bio = (r.get("biomass_recipe_id") or "").strip()
        if g not in cov_by_genome:
            raise ValueError(f"Summary row genome {g!r} not present in covariate supplement")
        c_row = cov_by_genome[g]
        for col in ("cohort", "canonical_order", "canonical_family", "canonical_genus"):
            if (r.get(col) or "").strip() != (c_row.get(col) or "").strip():
                raise ValueError(
                    f"Summary vs covariate metadata mismatch on {col} for {g}: "
                    f"{r.get(col)!r} vs {c_row.get(col)!r}"
                )
        key = (g, scen, bio)
        if key in summary_by_key:
            raise ValueError(f"Duplicate summary case key: {key}")
        is_valid_summary_endpoint(r, strict_diagnostics=True)
        summary_by_key[key] = r
        cases_by_genome[g].add((scen, bio))

    expected_core_set = set(EXPECTED_24_CORE_COMBINATIONS)
    expected_opt_set = set(EXPECTED_OPTIONAL_S0_COMBINATIONS)
    expected_all_allowed = expected_core_set | expected_opt_set

    for g, pairs in cases_by_genome.items():
        unknown_pairs = pairs - expected_all_allowed
        if unknown_pairs:
            raise ValueError(f"Unexpected scenario/biomass pairs for {g}: {sorted(unknown_pairs)}")
        missing_core = expected_core_set - pairs
        if missing_core:
            raise ValueError(f"Missing core scenario/biomass pairs for {g}: {sorted(missing_core)}")
        opt_present = pairs & expected_opt_set
        if g in expected_Empty_s0_genomes:
            if opt_present:
                raise ValueError(
                    f"Genome {g} expected to have empty optional S0 medium, but found {sorted(opt_present)}"
                )
        else:
            if opt_present != expected_opt_set:
                raise ValueError(
                    f"Genome {g} missing optional S0 combinations: expected {sorted(expected_opt_set)}, got {sorted(opt_present)}"
                )

    return summary_by_key, cov_by_genome


def build_table_a_case_universe(
    summary_rows: list[dict[str, str]],
    cov_by_genome: dict[str, dict[str, str]],
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Build Part A table: exact case universe and solver status by cohort/source/order/scenario/biomass."""
    rows_out: list[dict[str, Any]] = []

    def summarize_subset(
        stratum_type: str,
        stratum_id: str,
        cohort_label: str,
        source_label: str,
        order_label: str,
        scenario_label: str,
        biomass_label: str,
        genomes_in_scope: set[str],
        subset_rows: list[dict[str, str]],
        core_per_genome: int,
        opt_per_genome: int,
    ) -> dict[str, Any]:
        n_genomes = len(genomes_in_scope)
        n_planned_core = n_genomes * core_per_genome
        n_empty_s0_genomes = len(genomes_in_scope & VERIFIED_EMPTY_S0_GENOMES)
        n_opt_absent = n_empty_s0_genomes * opt_per_genome
        n_opt_present = (n_genomes - n_empty_s0_genomes) * opt_per_genome
        n_eval = len(subset_rows)
        st_counts = Counter(r["solver_status"] for r in subset_rows)
        n_valid = sum(1 for r in subset_rows if is_valid_summary_endpoint(r))
        n_infeas = st_counts.get("infeasible_under_fixed_medium", 0)
        n_fba_fail = st_counts.get("fba_verification_failed", 0)
        n_timeout = st_counts.get("timeout_no_incumbent", 0)
        n_num_unk = n_fba_fail + n_timeout
        assert n_valid + n_infeas + n_num_unk == n_eval

        sol_counts = [
            int(float(r["n_distinct_valid_solutions"]))
            for r in subset_rows
            if is_valid_summary_endpoint(r)
        ]
        min_sol = min(sol_counts) if sol_counts else None
        med_sol = float(statistics.median(sol_counts)) if sol_counts else None
        max_sol = max(sol_counts) if sol_counts else None

        return {
            "stratum_type": stratum_type,
            "stratum_id": stratum_id,
            "cohort": cohort_label,
            "source": source_label,
            "canonical_order": order_label,
            "scenario_id": scenario_label,
            "biomass_recipe_id": biomass_label,
            "n_genomes": n_genomes,
            "n_planned_core_cases": n_planned_core,
            "n_optional_s0_present": n_opt_present,
            "n_optional_s0_absent_empty_medium": n_opt_absent,
            "n_total_evaluated_cases": n_eval,
            "n_valid_feasible_cases": n_valid,
            "n_fixed_medium_infeasible_cases": n_infeas,
            "n_fba_verification_failed_unknown": n_fba_fail,
            "n_timeout_no_incumbent_unknown": n_timeout,
            "n_total_numerical_unknown_cases": n_num_unk,
            "prop_valid_feasible_of_evaluated": fmt_num(n_valid / n_eval if n_eval else None, 6),
            "prop_fixed_medium_infeasible_of_evaluated": fmt_num(n_infeas / n_eval if n_eval else None, 6),
            "prop_numerical_unknown_of_evaluated": fmt_num(n_num_unk / n_eval if n_eval else None, 6),
            "min_valid_solutions_among_feasible": fmt_num(min_sol, 0),
            "median_valid_solutions_among_feasible": fmt_num(med_sol, 2),
            "max_valid_solutions_among_feasible": fmt_num(max_sol, 0),
        }

    all_genomes = set(cov_by_genome.keys())
    rows_out.append(
        summarize_subset(
            "overall_all_484",
            "all_484",
            "all_484",
            "all_sources",
            "all_orders",
            "all_scenarios",
            "all_biomass_recipes",
            all_genomes,
            summary_rows,
            24,
            2,
        )
    )

    for coh in ("EQ_desert", "external_comparator"):
        gset = {g for g, c in cov_by_genome.items() if c["cohort"] == coh}
        sub = [r for r in summary_rows if r["genome"] in gset]
        rows_out.append(
            summarize_subset(
                "by_cohort",
                coh,
                coh,
                "all_cohort_sources",
                "all_orders",
                "all_scenarios",
                "all_biomass_recipes",
                gset,
                sub,
                24,
                2,
            )
        )

    for src in ("EQ_MAG", "GTDB_r232_rep", "MGnify"):
        gset = {g for g, c in cov_by_genome.items() if c["source"] == src}
        if not gset:
            continue
        coh = next(iter(cov_by_genome[g]["cohort"] for g in gset))
        sub = [r for r in summary_rows if r["genome"] in gset]
        rows_out.append(
            summarize_subset(
                "by_source",
                src,
                coh,
                src,
                "all_orders",
                "all_scenarios",
                "all_biomass_recipes",
                gset,
                sub,
                24,
                2,
            )
        )

    all_combos = EXPECTED_OPTIONAL_S0_COMBINATIONS + EXPECTED_24_CORE_COMBINATIONS
    for scen, bio in all_combos:
        is_opt = scen == S0
        sub_all = [r for r in summary_rows if r["scenario_id"] == scen and r["biomass_recipe_id"] == bio]
        rows_out.append(
            summarize_subset(
                "by_scenario_biomass_all",
                f"all_484|{scen}|{bio}",
                "all_484",
                "all_sources",
                "all_orders",
                scen,
                bio,
                all_genomes,
                sub_all,
                0 if is_opt else 1,
                1 if is_opt else 0,
            )
        )
        for coh in ("EQ_desert", "external_comparator"):
            gset = {g for g, c in cov_by_genome.items() if c["cohort"] == coh}
            sub_coh = [r for r in sub_all if r["genome"] in gset]
            rows_out.append(
                summarize_subset(
                    "by_cohort_scenario_biomass",
                    f"{coh}|{scen}|{bio}",
                    coh,
                    "all_cohort_sources",
                    "all_orders",
                    scen,
                    bio,
                    gset,
                    sub_coh,
                    0 if is_opt else 1,
                    1 if is_opt else 0,
                )
            )

    orders_sorted = sorted({c["canonical_order"] for c in cov_by_genome.values()})
    for coh in ("EQ_desert", "external_comparator"):
        for ord_name in orders_sorted:
            gset = {
                g
                for g, c in cov_by_genome.items()
                if c["cohort"] == coh and c["canonical_order"] == ord_name
            }
            if not gset:
                continue
            sub_ord = [r for r in summary_rows if r["genome"] in gset]
            rows_out.append(
                summarize_subset(
                    "by_cohort_and_order",
                    f"{coh}|{ord_name}",
                    coh,
                    "all_cohort_sources",
                    ord_name,
                    "all_scenarios",
                    "all_biomass_recipes",
                    gset,
                    sub_ord,
                    24,
                    2,
                )
            )

    overall_row = rows_out[0]
    part_a_summary = {
        "n_genomes": overall_row["n_genomes"],
        "n_planned_core_cases": overall_row["n_planned_core_cases"],
        "n_optional_s0_present": overall_row["n_optional_s0_present"],
        "n_optional_s0_absent_empty_medium": overall_row["n_optional_s0_absent_empty_medium"],
        "optional_s0_absent_genomes": sorted(all_genomes & VERIFIED_EMPTY_S0_GENOMES),
        "n_total_evaluated_cases": overall_row["n_total_evaluated_cases"],
        "solver_status_counts": {
            "heuristic_incumbent_bounded": overall_row["n_valid_feasible_cases"],
            "infeasible_under_fixed_medium": overall_row["n_fixed_medium_infeasible_cases"],
            "fba_verification_failed": overall_row["n_fba_verification_failed_unknown"],
            "timeout_no_incumbent": overall_row["n_timeout_no_incumbent_unknown"],
        },
        "n_total_numerical_unknown_cases": overall_row["n_total_numerical_unknown_cases"],
    }
    return rows_out, part_a_summary


def build_table_b_paired_interventions(
    summary_by_key: dict[tuple[str, str, str], dict[str, str]],
    cov_by_genome: dict[str, dict[str, str]],
    root_paired_audit_path: Path | None,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Build Part B table: within-genome paired biomass and redox sensitivity.

    Separates non-estimable pairs explicitly into:
    - optionalmedium_unavailable_pairs (S0 empty native/desert medium)
    - fixedmedium_infeasible_pairs (infeasible_under_fixed_medium; never labelled numerically unresolved)
    - numerical_unknown_pairs (fba_verification_failed or timeout_no_incumbent)
    """
    biomass_pairs = [
        ("biomass_B0_minus_B2", B0, B2, [S0, S1, S2, S3, S4, S5, S6]),
        ("biomass_B1_minus_B2_f420_contrast", B1, B2, [S3, S4, S5, S6]),
        ("biomass_B0_minus_B1", B0, B1, [S3, S4, S5, S6]),
        ("biomass_B0_minus_B3", B0, B3, [S3, S4, S5, S6]),
        ("biomass_B0_minus_B4", B0, B4, [S3, S4, S5, S6]),
    ]
    redox_pairs = [
        ("redox_S3_anoxic_minus_S4_oxic", S3, S4, [B0, B1, B2, B3, B4]),
        ("redox_S5_anoxic_minus_S6_oxic", S5, S6, [B0, B1, B2, B3, B4]),
    ]

    cohort_scopes: list[tuple[str, set[str]]] = [
        ("EQ_desert", {g for g, c in cov_by_genome.items() if c["cohort"] == "EQ_desert"}),
        (
            "external_comparator",
            {g for g, c in cov_by_genome.items() if c["cohort"] == "external_comparator"},
        ),
        ("all_484", set(cov_by_genome.keys())),
        (
            "shared_orders_checkm2_293_EQ_desert",
            {
                g
                for g, c in cov_by_genome.items()
                if c["cohort"] == "EQ_desert"
                and c["shared_order_stratum"] == "True"
                and c["checkm2_method_family_compatible"] == "True"
            },
        ),
        (
            "shared_orders_checkm2_293_external_gtdb",
            {
                g
                for g, c in cov_by_genome.items()
                if c["cohort"] == "external_comparator"
                and c["shared_order_stratum"] == "True"
                and c["checkm2_method_family_compatible"] == "True"
            },
        ),
    ]

    rows_out: list[dict[str, Any]] = []
    b0_vs_b2_report_dict: dict[str, dict[str, Any]] = {}
    b0_vs_b2_crosscheck_dict: dict[str, dict[str, Any]] = {}

    def evaluate_pair_group(
        intervention_type: str,
        contrast_id: str,
        cohort_scope: str,
        genomes_in_scope: set[str],
        fixed_context_id: str,
        ref_id: str,
        alt_id: str,
        is_biomass_intervention: bool,
    ) -> dict[str, Any]:
        n_scope = len(genomes_in_scope)
        n_present = 0
        n_opt_unavail = 0
        n_infeas = 0
        n_num_unk = 0
        ref_vals: list[float] = []
        alt_vals: list[float] = []
        diffs: list[float] = []

        for g in sorted(genomes_in_scope):
            if is_biomass_intervention:
                k_ref = (g, fixed_context_id, ref_id)
                k_alt = (g, fixed_context_id, alt_id)
            else:
                k_ref = (g, ref_id, fixed_context_id)
                k_alt = (g, alt_id, fixed_context_id)
            r_ref = summary_by_key.get(k_ref)
            r_alt = summary_by_key.get(k_alt)
            if r_ref is None and r_alt is None:
                n_opt_unavail += 1
                continue
            if r_ref is None or r_alt is None:
                raise ValueError(
                    f"Unpaired single endpoint present for genome={g}: ref={k_ref in summary_by_key}, alt={k_alt in summary_by_key}"
                )
            n_present += 1
            v_ref = is_valid_summary_endpoint(r_ref)
            v_alt = is_valid_summary_endpoint(r_alt)
            if v_ref and v_alt:
                rv = float(r_ref["min_added_rxns_count"])
                av = float(r_alt["min_added_rxns_count"])
                ref_vals.append(rv)
                alt_vals.append(av)
                diffs.append(rv - av)
            else:
                st_ref = r_ref["solver_status"]
                st_alt = r_alt["solver_status"]
                if st_ref in INFEASIBLE_SOLVER_STATUSES or st_alt in INFEASIBLE_SOLVER_STATUSES:
                    n_infeas += 1
                else:
                    n_num_unk += 1

        n_usable = len(diffs)
        n_nonestimable_evaluated = n_present - n_usable
        assert n_nonestimable_evaluated == n_infeas + n_num_unk

        ref_q = quantile_summary(ref_vals)
        alt_q = quantile_summary(alt_vals)
        diff_q = quantile_summary(diffs)

        n_pos = sum(1 for d in diffs if d > 0.0)
        n_tie = sum(1 for d in diffs if d == 0.0)
        n_neg = sum(1 for d in diffs if d < 0.0)

        return {
            "intervention_type": intervention_type,
            "contrast_id": contrast_id,
            "cohort_scope": cohort_scope,
            "fixed_context_id": fixed_context_id,
            "ref_endpoint_id": ref_id,
            "alt_endpoint_id": alt_id,
            "n_cohort_genomes": n_scope,
            "n_present_pairs": n_present,
            "n_optionalmedium_unavailable_pairs": n_opt_unavail,
            "n_paired_usable": n_usable,
            "n_nonestimable_evaluated_pairs": n_nonestimable_evaluated,
            "n_fixedmedium_infeasible_pairs": n_infeas,
            "n_numerical_unknown_pairs": n_num_unk,
            "ref_median_min_added_rxns": fmt_num(ref_q["median"], 2),
            "ref_q25": fmt_num(ref_q["q25"], 2),
            "ref_q75": fmt_num(ref_q["q75"], 2),
            "ref_iqr": fmt_num(ref_q["iqr"], 2),
            "alt_median_min_added_rxns": fmt_num(alt_q["median"], 2),
            "alt_q25": fmt_num(alt_q["q25"], 2),
            "alt_q75": fmt_num(alt_q["q75"], 2),
            "alt_iqr": fmt_num(alt_q["iqr"], 2),
            "diff_median_ref_minus_alt": fmt_num(diff_q["median"], 2),
            "diff_q25": fmt_num(diff_q["q25"], 2),
            "diff_q75": fmt_num(diff_q["q75"], 2),
            "diff_iqr": fmt_num(diff_q["iqr"], 2),
            "diff_min": fmt_num(diff_q["min"], 2),
            "diff_max": fmt_num(diff_q["max"], 2),
            "n_diff_gt_0_reduction": n_pos,
            "n_diff_eq_0_tie": n_tie,
            "n_diff_lt_0_increase": n_neg,
            "prop_diff_gt_0": fmt_num(n_pos / n_usable if n_usable else None, 4),
            "prop_diff_eq_0": fmt_num(n_tie / n_usable if n_usable else None, 4),
            "prop_diff_lt_0": fmt_num(n_neg / n_usable if n_usable else None, 4),
            "metric_scope_note": "min_added_rxns_among_found_valid_solutions_not_proven_global_minimum",
            "_report_dict": {
                "cohort_genomes": n_scope,
                "present_pairs": n_present,
                "optionalmedium_unavailable_pairs": n_opt_unavail,
                "paired_usable": n_usable,
                "nonestimable_pairs": n_nonestimable_evaluated,
                "fixedmedium_infeasible_pairs": n_infeas,
                "numerical_unknown_pairs": n_num_unk,
                "median_added_reaction_count_reduction": diff_q["median"],
                "q25": diff_q["q25"],
                "q75": diff_q["q75"],
                "n_reduction": n_pos,
                "n_tie": n_tie,
                "n_increase": n_neg,
            },
            "_crosscheck_dict": {
                "genomes": n_present,
                "paired_usable": n_usable,
                "unresolved_pairs": n_nonestimable_evaluated,
                "median_added_reaction_count_reduction": diff_q["median"],
                "q25": diff_q["q25"],
                "q75": diff_q["q75"],
                "n_reduction": n_pos,
                "n_tie": n_tie,
                "n_increase": n_neg,
            },
        }

    for contrast_id, b_ref, b_alt, scens in biomass_pairs:
        for coh_scope, gset in cohort_scopes:
            if not gset:
                continue
            for scen in scens:
                res = evaluate_pair_group(
                    "within_genome_same_medium_biomass_perturbation",
                    contrast_id,
                    coh_scope,
                    gset,
                    scen,
                    b_ref,
                    b_alt,
                    True,
                )
                if contrast_id == "biomass_B0_minus_B2" and coh_scope in ("EQ_desert", "external_comparator"):
                    k = f"{coh_scope}|{scen}"
                    b0_vs_b2_report_dict[k] = res["_report_dict"]
                    b0_vs_b2_crosscheck_dict[k] = res["_crosscheck_dict"]
                rows_out.append(
                    {k: v for k, v in res.items() if not k.startswith("_")}
                )

    for contrast_id, s_ref, s_alt, bios in redox_pairs:
        for coh_scope, gset in cohort_scopes:
            if not gset:
                continue
            for bio in bios:
                res = evaluate_pair_group(
                    "within_genome_same_biomass_redox_contrast",
                    contrast_id,
                    coh_scope,
                    gset,
                    bio,
                    s_ref,
                    s_alt,
                    False,
                )
                rows_out.append(
                    {k: v for k, v in res.items() if not k.startswith("_")}
                )

    crosscheck_passed = False
    if root_paired_audit_path is not None:
        if not root_paired_audit_path.is_file():
            raise FileNotFoundError(
                f"Explicit --root-paired-audit-json file not found: {root_paired_audit_path}"
            )
        with root_paired_audit_path.open("r", encoding="utf-8") as fh:
            root_audit = json.load(fh)["results"]
        if set(root_audit.keys()) != set(b0_vs_b2_crosscheck_dict.keys()):
            raise ValueError("Keys in root paired biomass audit JSON do not match computed B0 vs B2 keys")
        for k, exp_dict in root_audit.items():
            got_dict = b0_vs_b2_crosscheck_dict[k]
            if got_dict != exp_dict:
                raise ValueError(f"Root paired biomass sensitivity mismatch on {k}: got {got_dict} vs {exp_dict}")
        crosscheck_passed = True

    part_b_summary = {
        "root_paired_biomass_sensitivity_crosscheck_passed": crosscheck_passed,
        "b0_minus_b2_by_cohort_and_scenario": b0_vs_b2_report_dict,
    }
    return rows_out, part_b_summary


def build_tables_c_cohort_descriptive(
    summary_by_key: dict[tuple[str, str, str], dict[str, str]],
    cov_by_genome: dict[str, dict[str, str]],
) -> tuple[
    list[dict[str, Any]],
    list[dict[str, Any]],
    list[dict[str, Any]],
    list[dict[str, Any]],
    dict[str, Any],
]:
    """Build Part C descriptive cohort tables."""
    genome_bins: dict[str, dict[str, str]] = {}
    for g, c in cov_by_genome.items():
        is_cm2 = (c.get("checkm2_method_family_compatible") or "").strip() == "True"
        if is_cm2:
            cb = classify_completeness_bin(c.get("method_family_checkm2_completeness"))
            kb = classify_contamination_bin(c.get("method_family_checkm2_contamination"))
        else:
            cb = "missing"
            kb = "missing"
        rcb = classify_completeness_bin(c.get("reported_quality_completeness"))
        rkb = classify_contamination_bin(c.get("reported_quality_contamination"))
        genome_bins[g] = {
            "checkm2_comp_bin": cb,
            "checkm2_cont_bin": kb,
            "reported_comp_bin": rcb,
            "reported_cont_bin": rkb,
        }

    orders_all = sorted(
        {c["canonical_order"] for c in cov_by_genome.values()},
        key=lambda o: (0 if o in SHARED_ORDERS else 1, o),
    )
    comp_rows: list[dict[str, Any]] = []
    for ord_name in orders_all:
        is_shared = ord_name in SHARED_ORDERS
        for cb in COMPLETENESS_BINS:
            for kb in CONTAMINATION_BINS:
                eq_gs = [
                    g
                    for g, c in cov_by_genome.items()
                    if c["canonical_order"] == ord_name
                    and c["cohort"] == "EQ_desert"
                    and genome_bins[g]["checkm2_comp_bin"] == cb
                    and genome_bins[g]["checkm2_cont_bin"] == kb
                ]
                gtdb_gs = [
                    g
                    for g, c in cov_by_genome.items()
                    if c["canonical_order"] == ord_name
                    and c["source"] == "GTDB_r232_rep"
                    and genome_bins[g]["checkm2_comp_bin"] == cb
                    and genome_bins[g]["checkm2_cont_bin"] == kb
                ]
                mgnify_gs = [
                    g
                    for g, c in cov_by_genome.items()
                    if c["canonical_order"] == ord_name
                    and c["source"] == "MGnify"
                    and genome_bins[g]["checkm2_comp_bin"] == cb
                    and genome_bins[g]["checkm2_cont_bin"] == kb
                ]
                if not eq_gs and not gtdb_gs and not mgnify_gs:
                    continue
                is_cm2_cell = cb not in ("missing", "out_of_range") and kb not in ("missing", "out_of_range")
                both_cm2 = is_cm2_cell and len(eq_gs) > 0 and len(gtdb_gs) > 0
                sparse_cm2 = both_cm2 and (len(eq_gs) < 5 or len(gtdb_gs) < 5)

                eq_comp_vals = [
                    float(cov_by_genome[g]["method_family_checkm2_completeness"])
                    for g in eq_gs
                    if cov_by_genome[g]["method_family_checkm2_completeness"] != ""
                ]
                gtdb_comp_vals = [
                    float(cov_by_genome[g]["method_family_checkm2_completeness"])
                    for g in gtdb_gs
                    if cov_by_genome[g]["method_family_checkm2_completeness"] != ""
                ]

                comp_rows.append(
                    {
                        "canonical_order": ord_name,
                        "shared_order_stratum": str(is_shared),
                        "checkm2_completeness_bin": cb,
                        "checkm2_contamination_bin": kb,
                        "n_eq_desert_checkm2": len(eq_gs),
                        "n_external_gtdb_checkm2": len(gtdb_gs),
                        "n_external_mgnify_unknown_method": len(mgnify_gs),
                        "n_total_genomes_in_cell": len(eq_gs) + len(gtdb_gs) + len(mgnify_gs),
                        "both_cohorts_present_in_checkm2_cell": str(both_cm2),
                        "is_sparse_checkm2_overlap_n_lt_5": str(sparse_cm2),
                        "eq_checkm2_completeness_median": fmt_num(
                            statistics.median(eq_comp_vals) if eq_comp_vals else None, 2
                        ),
                        "gtdb_checkm2_completeness_median": fmt_num(
                            statistics.median(gtdb_comp_vals) if gtdb_comp_vals else None, 2
                        ),
                        "cell_provenance_note": (
                            "mgnify_unknown_quality_method_excluded_from_checkm2_strata"
                            if mgnify_gs
                            else "checkm2_method_family_subset_version_and_model_may_differ"
                        ),
                    }
                )

    strat_rows: list[dict[str, Any]] = []
    target_scens = [S3, S4, S5, S6]
    target_bios = [B0, B2]

    for ord_name in orders_all:
        is_shared = ord_name in SHARED_ORDERS
        for cb in COMPLETENESS_BINS:
            for kb in CONTAMINATION_BINS:
                if cb == "missing" or kb == "missing":
                    continue
                eq_gs = sorted(
                    g
                    for g, c in cov_by_genome.items()
                    if c["canonical_order"] == ord_name
                    and c["cohort"] == "EQ_desert"
                    and c["checkm2_method_family_compatible"] == "True"
                    and genome_bins[g]["checkm2_comp_bin"] == cb
                    and genome_bins[g]["checkm2_cont_bin"] == kb
                )
                ext_gs = sorted(
                    g
                    for g, c in cov_by_genome.items()
                    if c["canonical_order"] == ord_name
                    and c["cohort"] == "external_comparator"
                    and c["checkm2_method_family_compatible"] == "True"
                    and genome_bins[g]["checkm2_comp_bin"] == cb
                    and genome_bins[g]["checkm2_cont_bin"] == kb
                )
                if not eq_gs and not ext_gs:
                    continue

                for scen in target_scens:
                    for bio in target_bios:
                        eq_valid_vals: list[float] = []
                        eq_infeas = 0
                        eq_unk = 0
                        for g in eq_gs:
                            r = summary_by_key[(g, scen, bio)]
                            if is_valid_summary_endpoint(r):
                                eq_valid_vals.append(float(r["min_added_rxns_count"]))
                            elif r["solver_status"] in INFEASIBLE_SOLVER_STATUSES:
                                eq_infeas += 1
                            else:
                                eq_unk += 1

                        ext_valid_vals: list[float] = []
                        ext_infeas = 0
                        ext_unk = 0
                        for g in ext_gs:
                            r = summary_by_key[(g, scen, bio)]
                            if is_valid_summary_endpoint(r):
                                ext_valid_vals.append(float(r["min_added_rxns_count"]))
                            elif r["solver_status"] in INFEASIBLE_SOLVER_STATUSES:
                                ext_infeas += 1
                            else:
                                ext_unk += 1

                        eq_q = quantile_summary(eq_valid_vals)
                        ext_q = quantile_summary(ext_valid_vals)
                        both_valid = len(eq_valid_vals) > 0 and len(ext_valid_vals) > 0
                        is_sparse = both_valid and (len(eq_valid_vals) < 5 or len(ext_valid_vals) < 5)
                        if not both_valid:
                            contrast_status = "no_contrast_missing_cohort"
                            diff_med = None
                        elif is_sparse:
                            contrast_status = "descriptive_sparse_n_lt_5_flagged"
                            diff_med = float(eq_q["median"] - ext_q["median"])  # type: ignore[operator]
                        else:
                            contrast_status = "descriptive_non_sparse_n_ge_5"
                            diff_med = float(eq_q["median"] - ext_q["median"])  # type: ignore[operator]

                        strat_rows.append(
                            {
                                "canonical_order": ord_name,
                                "shared_order_stratum": str(is_shared),
                                "checkm2_completeness_bin": cb,
                                "checkm2_contamination_bin": kb,
                                "scenario_id": scen,
                                "biomass_recipe_id": bio,
                                "metric": "min_added_rxns_count",
                                "n_eq_total_in_stratum": len(eq_gs),
                                "n_eq_valid": len(eq_valid_vals),
                                "n_eq_fixed_medium_infeasible": eq_infeas,
                                "n_eq_numerical_unknown": eq_unk,
                                "eq_median": fmt_num(eq_q["median"], 2),
                                "eq_q25": fmt_num(eq_q["q25"], 2),
                                "eq_q75": fmt_num(eq_q["q75"], 2),
                                "eq_iqr": fmt_num(eq_q["iqr"], 2),
                                "n_ext_checkm2_total_in_stratum": len(ext_gs),
                                "n_ext_checkm2_valid": len(ext_valid_vals),
                                "n_ext_fixed_medium_infeasible": ext_infeas,
                                "n_ext_numerical_unknown": ext_unk,
                                "ext_checkm2_median": fmt_num(ext_q["median"], 2),
                                "ext_checkm2_q25": fmt_num(ext_q["q25"], 2),
                                "ext_checkm2_q75": fmt_num(ext_q["q75"], 2),
                                "ext_checkm2_iqr": fmt_num(ext_q["iqr"], 2),
                                "both_cohorts_valid_present": str(both_valid),
                                "is_sparse_n_lt_5": str(is_sparse),
                                "contrast_status": contrast_status,
                                "descriptive_eq_minus_ext_median_diff": fmt_num(diff_med, 2),
                                "interpretation_note": (
                                    "descriptive_only_no_adaptation_or_phylogeny_controlled_claim;"
                                    "checkm2_version_and_training_model_may_differ"
                                ),
                            }
                        )

    # Table C3: MGnify 29 unknown-method reported-only sensitivity (never pooled into CheckM2)
    mgnify_rows: list[dict[str, Any]] = []
    mgnify_all_gs = sorted(g for g, c in cov_by_genome.items() if c["source"] == "MGnify")
    for ord_scope_label, gs_filter in [
        ("all_mgnify_29", mgnify_all_gs),
        ("shared_orders_mgnify_8", [g for g in mgnify_all_gs if cov_by_genome[g]["shared_order_stratum"] == "True"]),
        ("non_shared_orders_mgnify_21", [g for g in mgnify_all_gs if cov_by_genome[g]["shared_order_stratum"] == "False"]),
    ]:
        for rcb in ("all_reported_bins", "[50,70)", "[70,90)", "[90,101]"):
            for rkb in ("all_reported_bins", "[0,5)", "[5,10)"):
                if (rcb == "all_reported_bins") != (rkb == "all_reported_bins"):
                    continue
                if rcb == "all_reported_bins":
                    sub_gs = gs_filter
                else:
                    sub_gs = [
                        g
                        for g in gs_filter
                        if genome_bins[g]["reported_comp_bin"] == rcb
                        and genome_bins[g]["reported_cont_bin"] == rkb
                    ]
                if not sub_gs:
                    continue
                for scen in target_scens:
                    for bio in target_bios:
                        vals: list[float] = []
                        n_infeas = 0
                        n_num_unk = 0
                        for g in sub_gs:
                            r = summary_by_key[(g, scen, bio)]
                            if is_valid_summary_endpoint(r):
                                vals.append(float(r["min_added_rxns_count"]))
                            elif r["solver_status"] in INFEASIBLE_SOLVER_STATUSES:
                                n_infeas += 1
                            else:
                                n_num_unk += 1
                        q = quantile_summary(vals)
                        mgnify_rows.append(
                            {
                                "mgnify_stratum": ord_scope_label,
                                "reported_quality_method": "unknown_mgnify_catalogue_pipeline",
                                "reported_completeness_bin": rcb,
                                "reported_contamination_bin": rkb,
                                "scenario_id": scen,
                                "biomass_recipe_id": bio,
                                "n_mgnify_total": len(sub_gs),
                                "n_mgnify_valid": len(vals),
                                "n_mgnify_fixedmedium_infeasible": n_infeas,
                                "n_mgnify_numerical_unknown": n_num_unk,
                                "mgnify_min_added_rxns_median": fmt_num(q["median"], 2),
                                "mgnify_min_added_rxns_q25": fmt_num(q["q25"], 2),
                                "mgnify_min_added_rxns_q75": fmt_num(q["q75"], 2),
                                "mgnify_min_added_rxns_iqr": fmt_num(q["iqr"], 2),
                                "separation_note": "reported_quality_method_unknown_never_pooled_into_checkm2_strata",
                            }
                        )

    # Table C4: Unadjusted all-genome descriptive table (no Mann-Whitney p-values)
    unadj_rows: list[dict[str, Any]] = []
    eq_all_gs = sorted(g for g, c in cov_by_genome.items() if c["cohort"] == "EQ_desert")
    ext_all_gs = sorted(g for g, c in cov_by_genome.items() if c["cohort"] == "external_comparator")
    gtdb_all_gs = sorted(g for g, c in cov_by_genome.items() if c["source"] == "GTDB_r232_rep")

    for scen, bio in EXPECTED_OPTIONAL_S0_COMBINATIONS + EXPECTED_24_CORE_COMBINATIONS:
        for metric_col in (
            "min_added_rxns_count",
            "min_total_cost",
            "n_rxns_observed_in_all_enumerated_models",
            "n_rxns_essential_in_declared_supernetwork",
        ):
            def _group_stats(gs: list[str]) -> tuple[int, int, int, int, dict[str, float | None]]:
                n_pres = 0
                n_inf = 0
                n_unk = 0
                vlist: list[float] = []
                for g in gs:
                    r = summary_by_key.get((g, scen, bio))
                    if r is None:
                        continue
                    n_pres += 1
                    if is_valid_summary_endpoint(r):
                        vlist.append(float(r[metric_col]))
                    elif r["solver_status"] in INFEASIBLE_SOLVER_STATUSES:
                        n_inf += 1
                    else:
                        n_unk += 1
                return n_pres, len(vlist), n_inf, n_unk, quantile_summary(vlist)

            eq_pres, eq_val, eq_inf, eq_unk, eq_q = _group_stats(eq_all_gs)
            ext_pres, ext_val, ext_inf, ext_unk, ext_q = _group_stats(ext_all_gs)
            gtdb_pres, gtdb_val, _, _, gtdb_q = _group_stats(gtdb_all_gs)
            mgn_pres, mgn_val, _, _, mgn_q = _group_stats(mgnify_all_gs)

            both_val = eq_val > 0 and ext_val > 0
            diff_unadj = (
                float(eq_q["median"] - ext_q["median"])  # type: ignore[operator]
                if both_val
                else None
            )

            unadj_rows.append(
                {
                    "scenario_id": scen,
                    "biomass_recipe_id": bio,
                    "metric": metric_col,
                    "n_eq_present": eq_pres,
                    "n_eq_valid": eq_val,
                    "n_eq_fixedmedium_infeasible": eq_inf,
                    "n_eq_numerical_unknown": eq_unk,
                    "eq_median_unadjusted": fmt_num(eq_q["median"], 2),
                    "eq_q25_unadjusted": fmt_num(eq_q["q25"], 2),
                    "eq_q75_unadjusted": fmt_num(eq_q["q75"], 2),
                    "eq_iqr_unadjusted": fmt_num(eq_q["iqr"], 2),
                    "n_ext_all_present": ext_pres,
                    "n_ext_all_valid": ext_val,
                    "n_ext_all_fixedmedium_infeasible": ext_inf,
                    "n_ext_all_numerical_unknown": ext_unk,
                    "ext_all_median_unadjusted": fmt_num(ext_q["median"], 2),
                    "ext_all_q25_unadjusted": fmt_num(ext_q["q25"], 2),
                    "ext_all_q75_unadjusted": fmt_num(ext_q["q75"], 2),
                    "ext_all_iqr_unadjusted": fmt_num(ext_q["iqr"], 2),
                    "n_ext_gtdb218_valid": gtdb_val,
                    "ext_gtdb218_median": fmt_num(gtdb_q["median"], 2),
                    "n_ext_mgnify29_valid": mgn_val,
                    "ext_mgnify29_median": fmt_num(mgn_q["median"], 2),
                    "unadjusted_eq_minus_ext_median_diff": fmt_num(diff_unadj, 2),
                    "confounding_warning": (
                        "descriptive_unadjusted_only_confounded_by_deep_lineage_composition_and_completeness"
                    ),
                }
            )

    non_sparse_s3_b0 = [
        r
        for r in strat_rows
        if r["scenario_id"] == S3
        and r["biomass_recipe_id"] == B0
        and r["contrast_status"] == "descriptive_non_sparse_n_ge_5"
    ]
    non_sparse_s3_b2 = [
        r
        for r in strat_rows
        if r["scenario_id"] == S3
        and r["biomass_recipe_id"] == B2
        and r["contrast_status"] == "descriptive_non_sparse_n_ge_5"
    ]

    part_c_summary = {
        "checkm2_method_family_compatible_total": sum(
            1 for c in cov_by_genome.values() if c["checkm2_method_family_compatible"] == "True"
        ),
        "checkm2_shared_orders_total": sum(
            1
            for c in cov_by_genome.values()
            if c["checkm2_method_family_compatible"] == "True" and c["shared_order_stratum"] == "True"
        ),
        "mgnify_unknown_method_separate_total": len(mgnify_all_gs),
        "n_order_quality_cells_total": len(comp_rows),
        "n_checkm2_cells_with_both_cohorts": sum(
            1 for r in comp_rows if r["both_cohorts_present_in_checkm2_cell"] == "True"
        ),
        "n_checkm2_cells_non_sparse_n_ge_5": sum(
            1
            for r in comp_rows
            if r["both_cohorts_present_in_checkm2_cell"] == "True"
            and r["is_sparse_checkm2_overlap_n_lt_5"] == "False"
        ),
        "non_sparse_s3_b0_cells": [
            {
                "order": r["canonical_order"],
                "comp_bin": r["checkm2_completeness_bin"],
                "cont_bin": r["checkm2_contamination_bin"],
                "n_eq": r["n_eq_valid"],
                "n_ext": r["n_ext_checkm2_valid"],
                "eq_median": float(r["eq_median"]),
                "ext_median": float(r["ext_checkm2_median"]),
                "eq_minus_ext": float(r["descriptive_eq_minus_ext_median_diff"]),
            }
            for r in non_sparse_s3_b0
        ],
        "non_sparse_s3_b2_cells": [
            {
                "order": r["canonical_order"],
                "comp_bin": r["checkm2_completeness_bin"],
                "cont_bin": r["checkm2_contamination_bin"],
                "n_eq": r["n_eq_valid"],
                "n_ext": r["n_ext_checkm2_valid"],
                "eq_median": float(r["eq_median"]),
                "ext_median": float(r["ext_checkm2_median"]),
                "eq_minus_ext": float(r["descriptive_eq_minus_ext_median_diff"]),
            }
            for r in non_sparse_s3_b2
        ],
    }
    return comp_rows, strat_rows, mgnify_rows, unadj_rows, part_c_summary


def build_table_d_masking_7genomes(
    mask_rows: list[dict[str, str]],
    base_rows: list[dict[str, str]],
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Summarize the 63 masking evaluations across the 7 high-completeness, low-contamination reference genomes."""
    if len(base_rows) != 7:
        raise ValueError(f"Expected 7 masking baseline genomes, got {len(base_rows)}")
    if len(mask_rows) != 63:
        raise ValueError(f"Expected 63 masking control rows, got {len(mask_rows)}")

    base_by_genome = {r["genome"].strip(): r for r in base_rows}
    if len(base_by_genome) != 7:
        raise ValueError("Duplicate genome IDs in masking baselines")

    seen_tuples: set[tuple[str, str, str]] = set()
    by_genome: dict[str, list[dict[str, str]]] = {g: [] for g in base_by_genome}
    for r in mask_rows:
        g = r["genome"].strip()
        if g not in base_by_genome:
            raise ValueError(f"Masked case genome {g} missing from baselines")
        mf = r["mask_fraction"].strip()
        rep = r["replicate"].strip()
        key = (g, mf, rep)
        if key in seen_tuples:
            raise ValueError(f"Duplicate masking control tuple: {key}")
        seen_tuples.add(key)
        by_genome[g].append(r)

    rows_out: list[dict[str, Any]] = []
    genome_mean_recoveries: list[float] = []
    genome_median_recoveries: list[float] = []
    genome_mean_alt_fracs: list[float] = []

    for g in sorted(base_by_genome.keys()):
        b = base_by_genome[g]
        sub = by_genome[g]
        if len(sub) != 9:
            raise ValueError(f"Expected 9 masked replicates for {g}, got {len(sub)}")
        evaluable = [
            r
            for r in sub
            if r.get("is_evaluable") == "True"
            and int(float(r.get("baseline_n_valid_solutions") or 0)) > 0
            and int(float(r.get("masked_n_valid_solutions") or 0)) > 0
            and int(float(r.get("n_removed_gpr_reactions") or 0)) > 0
        ]
        recs = [float(r["masked_rxn_recovery_rate_under_medium"]) for r in evaluable]
        alts = [float(r["alternative_candidate_fraction_under_medium"]) for r in evaluable]
        rq = quantile_summary(recs)
        aq = quantile_summary(alts)
        if rq["mean"] is not None:
            genome_mean_recoveries.append(float(rq["mean"]))
        if rq["median"] is not None:
            genome_median_recoveries.append(float(rq["median"]))
        if aq["mean"] is not None:
            genome_mean_alt_fracs.append(float(aq["mean"]))

        rows_out.append(
            {
                "summary_level": "per_genome_across_9_nested_replicates",
                "genome": g,
                "canonical_order": b["canonical_order"],
                "canonical_family": b["canonical_family"],
                "scenario_id": b["scenario_id"],
                "biomass_recipe_id": b["biomass_recipe_id"],
                "mask_fraction": "all_3_fractions_0.1_0.2_0.3",
                "n_independent_genomes_in_row": 1,
                "n_nested_replicates_total": len(sub),
                "n_nested_replicates_evaluable": len(evaluable),
                "baseline_n_valid_solutions": int(float(b["n_valid_solutions"])),
                "baseline_n_gapfill_rxns": int(float(b["n_gapfill_rxns"])),
                "recovery_rate_mean": fmt_num(rq["mean"], 6),
                "recovery_rate_median": fmt_num(rq["median"], 6),
                "recovery_rate_q25": fmt_num(rq["q25"], 6),
                "recovery_rate_q75": fmt_num(rq["q75"], 6),
                "recovery_rate_min": fmt_num(rq["min"], 6),
                "recovery_rate_max": fmt_num(rq["max"], 6),
                "alternative_fraction_mean": fmt_num(aq["mean"], 6),
                "alternative_fraction_median": fmt_num(aq["median"], 6),
                "alternative_fraction_min": fmt_num(aq["min"], 6),
                "alternative_fraction_max": fmt_num(aq["max"], 6),
                "unit_of_analysis_note": "9_replicates_nested_within_1_reference_genome_3_fractions_x_3_seeds",
            }
        )

        for mf in ("0.1", "0.2", "0.3"):
            sub_mf = [r for r in evaluable if r["mask_fraction"] == mf]
            recs_mf = [float(r["masked_rxn_recovery_rate_under_medium"]) for r in sub_mf]
            alts_mf = [float(r["alternative_candidate_fraction_under_medium"]) for r in sub_mf]
            rq_mf = quantile_summary(recs_mf)
            aq_mf = quantile_summary(alts_mf)
            rows_out.append(
                {
                    "summary_level": "per_genome_by_mask_fraction_3_replicates",
                    "genome": g,
                    "canonical_order": b["canonical_order"],
                    "canonical_family": b["canonical_family"],
                    "scenario_id": b["scenario_id"],
                    "biomass_recipe_id": b["biomass_recipe_id"],
                    "mask_fraction": mf,
                    "n_independent_genomes_in_row": 1,
                    "n_nested_replicates_total": 3,
                    "n_nested_replicates_evaluable": len(sub_mf),
                    "baseline_n_valid_solutions": int(float(b["n_valid_solutions"])),
                    "baseline_n_gapfill_rxns": int(float(b["n_gapfill_rxns"])),
                    "recovery_rate_mean": fmt_num(rq_mf["mean"], 6),
                    "recovery_rate_median": fmt_num(rq_mf["median"], 6),
                    "recovery_rate_q25": fmt_num(rq_mf["q25"], 6),
                    "recovery_rate_q75": fmt_num(rq_mf["q75"], 6),
                    "recovery_rate_min": fmt_num(rq_mf["min"], 6),
                    "recovery_rate_max": fmt_num(rq_mf["max"], 6),
                    "alternative_fraction_mean": fmt_num(aq_mf["mean"], 6),
                    "alternative_fraction_median": fmt_num(aq_mf["median"], 6),
                    "alternative_fraction_min": fmt_num(aq_mf["min"], 6),
                    "alternative_fraction_max": fmt_num(aq_mf["max"], 6),
                    "unit_of_analysis_note": "3_random_seed_replicates_nested_within_1_genome_and_mask_fraction",
                }
            )

    gq_rec = quantile_summary(genome_mean_recoveries)
    gq_alt = quantile_summary(genome_mean_alt_fracs)
    rows_out.insert(
        0,
        {
            "summary_level": "cohort_7_reference_genomes_summary_of_genome_means",
            "genome": "ALL_7_REFERENCE_GENOMES",
            "canonical_order": "o__DHVEG-1,o__UBA202,o__JACQPN01",
            "canonical_family": "multiple",
            "scenario_id": S4,
            "biomass_recipe_id": B2,
            "mask_fraction": "all_3_fractions_0.1_0.2_0.3",
            "n_independent_genomes_in_row": len(genome_mean_recoveries),
            "n_nested_replicates_total": 63,
            "n_nested_replicates_evaluable": sum(
                1 for r in mask_rows if r.get("is_evaluable") == "True"
            ),
            "baseline_n_valid_solutions": "",
            "baseline_n_gapfill_rxns": "",
            "recovery_rate_mean": fmt_num(gq_rec["mean"], 6),
            "recovery_rate_median": fmt_num(gq_rec["median"], 6),
            "recovery_rate_q25": fmt_num(gq_rec["q25"], 6),
            "recovery_rate_q75": fmt_num(gq_rec["q75"], 6),
            "recovery_rate_min": fmt_num(gq_rec["min"], 6),
            "recovery_rate_max": fmt_num(gq_rec["max"], 6),
            "alternative_fraction_mean": fmt_num(gq_alt["mean"], 6),
            "alternative_fraction_median": fmt_num(gq_alt["median"], 6),
            "alternative_fraction_min": fmt_num(gq_alt["min"], 6),
            "alternative_fraction_max": fmt_num(gq_alt["max"], 6),
            "unit_of_analysis_note": "primary_unit_is_7_reference_genomes_with_9_nested_replicates_each",
        },
    )

    part_d_summary = {
        "control_genome_selection_note": (
            "7 high-completeness, low-contamination MGnify reference genomes "
            "(reported numeric completeness >=90%, contamination <=5%; formal MIMAG status unknown "
            "because rRNA/tRNA criteria are not recorded in the metabolic inventory)."
        ),
        "n_independent_control_genomes": 7,
        "n_nested_replicates_per_genome": 9,
        "n_total_masked_evaluations": 63,
        "mean_of_7_genome_mean_recovery_rates": gq_rec["mean"],
        "median_of_7_genome_mean_recovery_rates": gq_rec["median"],
        "min_of_7_genome_mean_recovery_rates": gq_rec["min"],
        "max_of_7_genome_mean_recovery_rates": gq_rec["max"],
        "mean_of_7_genome_mean_alternative_fractions": gq_alt["mean"],
        "median_of_7_genome_mean_alternative_fractions": gq_alt["median"],
        "per_genome_mean_recovery_rate": {
            r["genome"]: float(r["recovery_rate_mean"])
            for r in rows_out
            if r["summary_level"] == "per_genome_across_9_nested_replicates"
        },
    }
    return rows_out, part_d_summary


def generate_review_figures(
    outdir: Path,
    summary_by_key: dict[tuple[str, str, str], dict[str, str]],
    cov_by_genome: dict[str, dict[str, str]],
    comp_rows: list[dict[str, Any]],
    strat_rows: list[dict[str, Any]],
    mask_rows: list[dict[str, str]],
    base_rows: list[dict[str, str]],
) -> list[Path]:
    """Generate the 3 required review figures in PDF (and PNG)."""
    mpl_dir = outdir / ".mpl_cache"
    mpl_dir.mkdir(parents=True, exist_ok=True)
    os.environ["MPLCONFIGDIR"] = str(mpl_dir)
    os.environ["TMPDIR"] = str(mpl_dir)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update(
        {
            "font.size": 8.5,
            "axes.titlesize": 9.5,
            "axes.labelsize": 8.8,
            "xtick.labelsize": 7.8,
            "ytick.labelsize": 8.0,
            "legend.fontsize": 7.8,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )

    generated: list[Path] = []
    try:
        # ---------------------------------------------------------------------
        # Figure 1: Paired B0 vs B2 distributions & recipe sensitivity display
        # ---------------------------------------------------------------------
        fig1, (ax1a, ax1b) = plt.subplots(
            1, 2, figsize=(13.2, 5.1), gridspec_kw={"width_ratios": [1.38, 1.0]}
        )

        scens_feasible = [S0, S3, S4, S5, S6]
        scen_short = {
            S0: "S0 native (opt.)",
            S3: "S3 desert anoxic",
            S4: "S4 desert oxic",
            S5: "S5 defined anoxic",
            S6: "S6 defined oxic",
        }

        positions = []
        data_boxes = []
        colors_boxes = []
        xtick_pos = []
        xtick_labels = []

        eq_gs = sorted(g for g, c in cov_by_genome.items() if c["cohort"] == "EQ_desert")
        ext_gs = sorted(g for g, c in cov_by_genome.items() if c["cohort"] == "external_comparator")

        for idx, scen in enumerate(scens_feasible):
            base_x = idx * 3.1
            xtick_pos.append(base_x + 0.42)
            for offset, coh_label, gs, col in [
                (0.0, "EQ_desert", eq_gs, "#d95f02"),
                (0.84, "external", ext_gs, "#1b9e77"),
            ]:
                diffs = []
                n_pres = 0
                for g in gs:
                    r0 = summary_by_key.get((g, scen, B0))
                    r2 = summary_by_key.get((g, scen, B2))
                    if r0 is None and r2 is None:
                        continue
                    n_pres += 1
                    if (
                        r0 is not None
                        and r2 is not None
                        and is_valid_summary_endpoint(r0)
                        and is_valid_summary_endpoint(r2)
                    ):
                        diffs.append(float(r0["min_added_rxns_count"]) - float(r2["min_added_rxns_count"]))
                positions.append(base_x + offset)
                data_boxes.append(diffs)
                colors_boxes.append(col)
                if offset == 0.0:
                    eq_n_str = f"EQ: {len(diffs)}/{n_pres}"
                else:
                    ext_n_str = f"Ext: {len(diffs)}/{n_pres}"
            xtick_labels.append(f"{scen_short[scen]}\n{eq_n_str}\n{ext_n_str}")

        bp = ax1a.boxplot(
            data_boxes,
            positions=positions,
            widths=0.64,
            patch_artist=True,
            showfliers=True,
            flierprops={"marker": "o", "markersize": 2.4, "alpha": 0.48},
            medianprops={"color": "black", "linewidth": 1.5},
        )
        for patch, col in zip(bp["boxes"], colors_boxes):
            patch.set_facecolor(col)
            patch.set_alpha(0.75)

        ax1a.axhline(0.0, color="#444444", linestyle="--", linewidth=1.0, label="No change (B0 - B2 = 0)")
        ax1a.set_xticks(xtick_pos)
        ax1a.set_xticklabels(xtick_labels)
        ax1a.set_ylabel("Paired difference in min added reactions (B0 minus B2)\namong found valid solutions")
        ax1a.set_title(
            "A. Within-genome paired B0 vs B2 sensitivity by cohort & medium\n"
            "(S1/S2 minimal H2+CO2+NH4: 0/237 EQ & 0/247 Ext feasible under fixed medium)"
        )
        ax1a.plot([], [], color="#d95f02", linewidth=6, alpha=0.75, label="EQ_desert (237 genomes)")
        ax1a.plot([], [], color="#1b9e77", linewidth=6, alpha=0.75, label="external_comparator (247 genomes)")
        # Place legend in lower-right under S5/S6 where all points are between +58 and +92 (zero overlap with outliers)
        ax1a.legend(loc="lower right", frameon=True)
        ax1a.grid(axis="y", linestyle=":", alpha=0.5)

        # Panel 1B: Discrete biomass recipe sensitivity summary in S3 and S5
        recipes_order = [
            (B0, "B0\n(template\nbaseline)"),
            (B2, "B2\n(omit 7 cof,\nkeep F420)"),
            (B1, "B1\n(omit 8 cof\nincl. F420)"),
            (B3, "B3\n(omit cof &\nglycolipids)"),
            (B4, "B4\n(core AA/nt\nskeleton)"),
        ]
        x_rec = np.arange(len(recipes_order))
        for scen, ls, marker, s_lbl in [(S3, "-", "o", "S3 desert"), (S5, "--", "s", "S5 defined")]:
            for coh_lbl, gs, col in [("EQ_desert", eq_gs, "#d95f02"), ("external", ext_gs, "#1b9e77")]:
                meds = []
                q25s = []
                q75s = []
                ns = []
                for b_id, _ in recipes_order:
                    vals = [
                        float(summary_by_key[(g, scen, b_id)]["min_added_rxns_count"])
                        for g in gs
                        if is_valid_summary_endpoint(summary_by_key[(g, scen, b_id)])
                    ]
                    ns.append(len(vals))
                    q = quantile_summary(vals)
                    meds.append(float(q["median"]))  # type: ignore[arg-type]
                    q25s.append(float(q["q25"]))  # type: ignore[arg-type]
                    q75s.append(float(q["q75"]))  # type: ignore[arg-type]
                n_min, n_max = min(ns), max(ns)
                n_tag = f"n={n_min}/{len(gs)}" if n_min == n_max else f"n={n_min}-{n_max}/{len(gs)}"
                ax1b.plot(
                    x_rec,
                    meds,
                    linestyle=ls,
                    marker=marker,
                    color=col,
                    label=f"{coh_lbl} {s_lbl} ({n_tag})",
                )
                ax1b.fill_between(x_rec, q25s, q75s, color=col, alpha=0.12)

        ax1b.set_xticks(x_rec)
        ax1b.set_xticklabels([lbl for _, lbl in recipes_order])
        ax1b.set_ylabel("Median (IQR) min added reactions among available valid solutions")
        ax1b.set_title(
            "B. Biomass recipe sensitivity summary (S3 & S5 available valid cases)\n"
            "(Similar median sensitivity in both cohorts)"
        )
        ax1b.legend(loc="upper right", frameon=True)
        ax1b.grid(axis="y", linestyle=":", alpha=0.5)

        fig1.tight_layout()
        fig1_pdf = outdir / "fig1_paired_b0_vs_b2_distributions.pdf"
        fig1_png = outdir / "fig1_paired_b0_vs_b2_distributions.png"
        fig1.savefig(fig1_pdf, bbox_inches="tight")
        fig1.savefig(fig1_png, dpi=200, bbox_inches="tight")
        plt.close(fig1)
        generated.extend([fig1_pdf, fig1_png])

        # ---------------------------------------------------------------------
        # Figure 2: Quality/order composition and overlap & stratified S3 medians
        # ---------------------------------------------------------------------
        fig2, (ax2a, ax2b, ax2c) = plt.subplots(1, 3, figsize=(15.2, 5.0))

        orders_plot = sorted(
            {c["canonical_order"] for c in cov_by_genome.values()},
            key=lambda o: (
                0 if o in SHARED_ORDERS else 1,
                -sum(1 for c in cov_by_genome.values() if c["canonical_order"] == o),
                o,
            ),
        )
        y_ord = np.arange(len(orders_plot))
        eq_ord_counts = [
            sum(1 for c in cov_by_genome.values() if c["canonical_order"] == o and c["source"] == "EQ_MAG")
            for o in orders_plot
        ]
        gtdb_ord_counts = [
            sum(1 for c in cov_by_genome.values() if c["canonical_order"] == o and c["source"] == "GTDB_r232_rep")
            for o in orders_plot
        ]
        mgn_ord_counts = [
            sum(1 for c in cov_by_genome.values() if c["canonical_order"] == o and c["source"] == "MGnify")
            for o in orders_plot
        ]

        ax2a.barh(y_ord, eq_ord_counts, color="#d95f02", alpha=0.85, label="EQ_MAG (CheckM2, n=237)")
        ax2a.barh(
            y_ord,
            gtdb_ord_counts,
            left=eq_ord_counts,
            color="#1b9e77",
            alpha=0.85,
            label="GTDB_r232_rep (CheckM2, n=218)",
        )
        ax2a.barh(
            y_ord,
            mgn_ord_counts,
            left=np.array(eq_ord_counts) + np.array(gtdb_ord_counts),
            color="#7570b3",
            alpha=0.85,
            label="MGnify (unknown method, n=29)",
        )
        ax2a.set_yticks(y_ord)
        ax2a.set_yticklabels([f"{o} *" if o in SHARED_ORDERS else o for o in orders_plot])
        ax2a.invert_yaxis()
        ax2a.set_xlabel("Number of genomes (total = 484)")
        ax2a.set_title("A. Order & quality-method composition\n(* = 3 shared orders: 234 EQ, 59 GTDB, 8 MGnify)")
        ax2a.legend(loc="lower right", frameon=True)
        ax2a.grid(axis="x", linestyle=":", alpha=0.5)

        # Panel 2B: CheckM2 completeness vs contamination in shared orders (n=293)
        for src, col, marker, lbl in [
            ("EQ_MAG", "#d95f02", "o", "EQ_MAG shared orders (n=234)"),
            ("GTDB_r232_rep", "#1b9e77", "^", "GTDB_r232_rep shared orders (n=59)"),
        ]:
            xs = [
                float(c["method_family_checkm2_completeness"])
                for c in cov_by_genome.values()
                if c["source"] == src and c["shared_order_stratum"] == "True"
            ]
            ys = [
                float(c["method_family_checkm2_contamination"])
                for c in cov_by_genome.values()
                if c["source"] == src and c["shared_order_stratum"] == "True"
            ]
            ax2b.scatter(xs, ys, c=col, marker=marker, s=22, alpha=0.65, label=lbl)

        for xv in (50.0, 70.0, 90.0):
            ax2b.axvline(xv, color="#555555", linestyle="--", linewidth=0.9)
        for yv in (5.0, 10.0):
            ax2b.axhline(yv, color="#555555", linestyle="--", linewidth=0.9)
        ax2b.set_xlim(48, 102)
        ax2b.set_ylim(-0.3, 10.5)
        ax2b.set_xlabel("CheckM2 completeness (%) [bins: [50,70), [70,90), [90,101]]")
        ax2b.set_ylabel("CheckM2 contamination (%) [bins: [0,5), [5,10)]")
        ax2b.set_title("B. Shared-order CheckM2 quality overlap (n=293)\n(MGnify n=29 unknown method kept separate)")
        ax2b.legend(loc="upper right", frameon=True)
        ax2b.grid(linestyle=":", alpha=0.4)

        # Panel 2C: Stratified S3 min_added_rxns_count (B0) in overlapping CheckM2 cells
        overlap_cells_s3_b0 = [
            r
            for r in strat_rows
            if r["scenario_id"] == S3
            and r["biomass_recipe_id"] == B0
            and r["both_cohorts_valid_present"] == "True"
        ]
        y_cells = np.arange(len(overlap_cells_s3_b0))
        cell_labels = []
        for i, r in enumerate(overlap_cells_s3_b0):
            sparse_tag = " [n<5 sparse]" if r["is_sparse_n_lt_5"] == "True" else ""
            ord_short = r["canonical_order"].replace("o__", "")
            cell_labels.append(
                f"{ord_short} {r['checkm2_completeness_bin']} {r['checkm2_contamination_bin']}\n"
                f"(EQ={r['n_eq_valid']}, GTDB={r['n_ext_checkm2_valid']}){sparse_tag}"
            )
            eq_m = float(r["eq_median"])
            ext_m = float(r["ext_checkm2_median"])
            ax2c.plot([eq_m, ext_m], [i, i], color="#666666", linewidth=1.2, zorder=1)
            ax2c.scatter([eq_m], [i], color="#d95f02", s=45, zorder=3, label="EQ median (S3, B0)" if i == 0 else "")
            ax2c.scatter(
                [ext_m],
                [i],
                color="#1b9e77",
                marker="^",
                s=48,
                zorder=3,
                label="GTDB median (S3, B0)" if i == 0 else "",
            )

        ax2c.set_yticks(y_cells)
        ax2c.set_yticklabels(cell_labels)
        ax2c.invert_yaxis()
        ax2c.set_xlabel("Median min added reactions in S3 (B0 baseline)")
        ax2c.set_title(
            "C. Within-stratum cohort medians (CheckM2 overlap cells)\n"
            "(Descriptive differences vary across order/quality strata)"
        )
        ax2c.legend(loc="lower right", frameon=True)
        ax2c.grid(axis="x", linestyle=":", alpha=0.5)

        fig2.tight_layout()
        fig2_pdf = outdir / "fig2_quality_order_composition_and_overlap.pdf"
        fig2_png = outdir / "fig2_quality_order_composition_and_overlap.png"
        fig2.savefig(fig2_pdf, bbox_inches="tight")
        fig2.savefig(fig2_png, dpi=200, bbox_inches="tight")
        plt.close(fig2)
        generated.extend([fig2_pdf, fig2_png])

        # ---------------------------------------------------------------------
        # Figure 3: Masking recovery & alternative fraction across 7 reference genomes
        # ---------------------------------------------------------------------
        fig3, (ax3a, ax3b) = plt.subplots(1, 2, figsize=(12.6, 4.7), sharey=True)

        genomes_7 = sorted(r["genome"].strip() for r in base_rows)
        base_map = {r["genome"].strip(): r for r in base_rows}
        y_g = np.arange(len(genomes_7))
        frac_colors = {"0.1": "#1f78b4", "0.2": "#33a02c", "0.3": "#e31a1c"}
        frac_offsets = {"0.1": -0.18, "0.2": 0.0, "0.3": 0.18}

        for mf in ("0.1", "0.2", "0.3"):
            lbl = f"Mask {int(float(mf)*100)}% (3 random seeds/genome)"
            ax3a.scatter([], [], color=frac_colors[mf], s=28, alpha=0.8, label=lbl)
            ax3b.scatter([], [], color=frac_colors[mf], s=28, alpha=0.8, label=lbl)

        g_labels = []
        for idx, g in enumerate(genomes_7):
            b = base_map[g]
            g_labels.append(f"{g}\n({b['canonical_order']}, comp={b['completeness']}%, cont={b['contamination']}%)")
            sub = [r for r in mask_rows if r["genome"].strip() == g]
            recs_all = [float(r["masked_rxn_recovery_rate_under_medium"]) for r in sub]
            alts_all = [float(r["alternative_candidate_fraction_under_medium"]) for r in sub]

            for r in sub:
                mf = r["mask_fraction"].strip()
                y_pos = idx + frac_offsets[mf]
                rec_v = float(r["masked_rxn_recovery_rate_under_medium"])
                alt_v = float(r["alternative_candidate_fraction_under_medium"])
                ax3a.scatter([rec_v], [y_pos], color=frac_colors[mf], s=26, alpha=0.78, zorder=2)
                ax3b.scatter([alt_v], [y_pos], color=frac_colors[mf], s=26, alpha=0.78, zorder=2)

            ax3a.scatter(
                [float(np.mean(recs_all))],
                [idx],
                color="black",
                marker="D",
                s=46,
                zorder=4,
                label="Per-genome mean (9 conditions)" if idx == 0 else "",
            )
            ax3a.scatter(
                [float(np.median(recs_all))],
                [idx],
                color="#ff7f00",
                marker="|",
                s=110,
                linewidths=2.2,
                zorder=5,
                label="Per-genome median" if idx == 0 else "",
            )

            ax3b.scatter(
                [float(np.mean(alts_all))],
                [idx],
                color="black",
                marker="D",
                s=46,
                zorder=4,
                label="Per-genome mean (9 conditions)" if idx == 0 else "",
            )
            ax3b.scatter(
                [float(np.median(alts_all))],
                [idx],
                color="#ff7f00",
                marker="|",
                s=110,
                linewidths=2.2,
                zorder=5,
                label="Per-genome median" if idx == 0 else "",
            )

        overall_rec_mean = float(np.mean([float(r["masked_rxn_recovery_rate_under_medium"]) for r in mask_rows]))
        overall_alt_mean = float(np.mean([float(r["alternative_candidate_fraction_under_medium"]) for r in mask_rows]))
        ax3a.axvline(
            overall_rec_mean, color="#333333", linestyle="--", linewidth=1.1, label=f"7-genome mean = {overall_rec_mean:.3f}"
        )
        ax3b.axvline(
            overall_alt_mean, color="#333333", linestyle="--", linewidth=1.1, label=f"7-genome mean = {overall_alt_mean:.3f}"
        )

        ax3a.set_yticks(y_g)
        ax3a.set_yticklabels(g_labels)
        ax3a.invert_yaxis()
        ax3a.set_xlim(0.0, 0.52)
        ax3a.set_xlabel("Masked GPR reaction recovery rate under matched S4 + B2 medium")
        ax3a.set_title(
            "A. Masked-reaction recovery across 7 reference genomes\n"
            "(7 genomes; 9 masking conditions/replicates per genome: 3 fractions × 3 seeds)"
        )
        ax3a.legend(loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=2, frameon=True)
        ax3a.grid(axis="x", linestyle=":", alpha=0.5)

        ax3b.set_xlim(0.40, 1.02)
        ax3b.set_xlabel("Alternative candidate fraction among newly induced gap-fill reactions")
        ax3b.set_title(
            "B. Alternative non-masked pathway fraction (7 reference genomes)\n"
            "(Alternative reactions account for 76.7% of new additions)"
        )
        ax3b.legend(loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=2, frameon=True)
        ax3b.grid(axis="x", linestyle=":", alpha=0.5)

        fig3.tight_layout()
        fig3.subplots_adjust(wspace=0.10)
        fig3_pdf = outdir / "fig3_masking_recovery_by_7genomes.pdf"
        fig3_png = outdir / "fig3_masking_recovery_by_7genomes.png"
        fig3.savefig(fig3_pdf, bbox_inches="tight")
        fig3.savefig(fig3_png, dpi=200, bbox_inches="tight")
        plt.close(fig3)
        generated.extend([fig3_pdf, fig3_png])
    finally:
        shutil.rmtree(mpl_dir, ignore_errors=True)

    return generated


def run_interpretation_pipeline(
    *,
    shards_dir: Path | None,
    combined_summary: Path | None,
    covariate_supplement: Path,
    inventory: Path | None,
    masking_cases: Path,
    masking_baselines: Path,
    root_paired_audit_json: Path | None,
    outdir: Path,
    strict_source_sha256: bool = True,
) -> dict[str, Any]:
    """Execute the full descriptive interpretation pipeline and write all artifacts to outdir."""
    if str(outdir).startswith("/tmp"):
        raise ValueError(f"Output directory must never be inside /tmp: {outdir}")
    outdir.mkdir(parents=True, exist_ok=True)

    for label, p in [
        ("--covariate-supplement", covariate_supplement),
        ("--masking-cases", masking_cases),
        ("--masking-baselines", masking_baselines),
    ]:
        if not p.is_file():
            raise FileNotFoundError(f"Explicit {label} file not found: {p}")

    if inventory is not None and not inventory.is_file():
        raise FileNotFoundError(f"Explicit --inventory file not found: {inventory}")
    if root_paired_audit_json is not None and not root_paired_audit_json.is_file():
        raise FileNotFoundError(f"Explicit --root-paired-audit-json file not found: {root_paired_audit_json}")

    summary_rows, input_hashes = load_summary_rows(
        shards_dir=shards_dir,
        combined_summary=combined_summary,
        strict_source_sha256=strict_source_sha256,
    )

    for p in [covariate_supplement, masking_cases, masking_baselines] + ([inventory] if inventory else []):
        sh = sha256_file(p)
        if strict_source_sha256 and p.name in EXPECTED_SOURCE_SHA256:
            exp_sh = EXPECTED_SOURCE_SHA256[p.name]
            if sh != exp_sh:
                raise ValueError(f"SHA-256 mismatch for {p.name}: got {sh}, expected {exp_sh}")
        input_hashes[p.name] = sh

    if root_paired_audit_json is not None:
        input_hashes[root_paired_audit_json.name] = sha256_file(root_paired_audit_json)

    cov_rows = read_tsv(covariate_supplement)
    inv_rows = read_tsv(inventory) if inventory is not None else None
    mask_rows = read_tsv(masking_cases)
    base_rows = read_tsv(masking_baselines)

    summary_by_key, cov_by_genome = validate_universe_and_join_covariates(
        summary_rows,
        cov_rows,
        inv_rows,
    )

    combined_out_path = outdir / "exp1_full484_combined_summary_12580.tsv"
    if summary_rows:
        summary_fields = list(summary_rows[0].keys())
        sorted_summary_rows = sorted(
            summary_rows,
            key=lambda r: (r["genome"], r["scenario_id"], r["biomass_recipe_id"]),
        )
        write_tsv(combined_out_path, summary_fields, sorted_summary_rows)

    rows_a, summary_a = build_table_a_case_universe(summary_rows, cov_by_genome)
    path_a = outdir / "exp1_case_universe_and_solver_status_summary.tsv"
    write_tsv(path_a, list(rows_a[0].keys()), rows_a)

    rows_b, summary_b = build_table_b_paired_interventions(
        summary_by_key,
        cov_by_genome,
        root_paired_audit_json,
    )
    path_b = outdir / "exp1_paired_biomass_and_redox_sensitivity.tsv"
    write_tsv(path_b, list(rows_b[0].keys()), rows_b)

    comp_rows, strat_rows, mgnify_rows, unadj_rows, summary_c = build_tables_c_cohort_descriptive(
        summary_by_key,
        cov_by_genome,
    )
    path_c1 = outdir / "exp1_order_and_quality_composition_overlap.tsv"
    write_tsv(path_c1, list(comp_rows[0].keys()), comp_rows)

    path_c2 = outdir / "exp1_checkm2_stratified_cohort_comparisons.tsv"
    write_tsv(path_c2, list(strat_rows[0].keys()), strat_rows)

    path_c3 = outdir / "exp1_mgnify29_reported_quality_sensitivity.tsv"
    write_tsv(path_c3, list(mgnify_rows[0].keys()), mgnify_rows)

    path_c4 = outdir / "exp1_unadjusted_all484_descriptive_summary.tsv"
    write_tsv(path_c4, list(unadj_rows[0].keys()), unadj_rows)

    rows_d, summary_d = build_table_d_masking_7genomes(mask_rows, base_rows)
    path_d = outdir / "exp1_masking_7genomes_descriptive_summary.tsv"
    write_tsv(path_d, list(rows_d[0].keys()), rows_d)

    fig_paths = generate_review_figures(
        outdir,
        summary_by_key,
        cov_by_genome,
        comp_rows,
        strat_rows,
        mask_rows,
        base_rows,
    )

    report_path = outdir / "exp1_descriptive_interpretation_report.json"
    report_payload = {
        "module": "exp1_metabolic/interpretation",
        "design_authority": [
            "mechanism-20261109/EXP1-DESCRIPTIVE-INTERPRETATION-IMPLEMENTATION.txt",
            "mechanism-20261109/EXP1-FULL-COHORT-INTERPRETATION-DESIGN.txt",
            "mechanism-20261109/EXP1-INTERPRETATION-FIGURE-REVIEW.txt",
        ],
        "scientific_scope_and_limitations": {
            "primary_axis": "Within-genome sensitivity to declared biomass template and medium recipes.",
            "optimality_caveat": (
                "All min_added_rxns_count values are minima among found valid time-limited MILP solutions "
                "(heuristic_incumbent_bounded), not proven global minima. Weighted reaction costs allow "
                "counter-directional reaction-count differences under recipe relaxation."
            ),
            "infeasibility_vs_unknown": (
                "1,936 fixed-medium infeasible cases (all 484 genomes in S1/S2 minimal H2+CO2+NH4 under B0/B2) "
                "are separated from 66 numerical unknown cases (56 FBA verification failures, 10 timeouts with no incumbent) "
                "and 4 optional S0 unavailable cases (2 genomes with verified empty native/desert media). "
                "Fixed-medium infeasibility is never described as numerically unresolved, and neither category is "
                "coerced to zero added reactions."
            ),
            "cohort_comparison_caveat": (
                "Between-cohort comparisons are descriptive only. Even within the 455-genome CheckM2 method-family "
                "subset and 293-genome shared-order subset (o__JACQPN01, o__JACPAO01, o__SW-10-69-26), CheckM2 software "
                "versions/models may differ, deep-lineage structure within orders is not fully phylogeny-controlled, "
                "and no environmental adaptation claim is made. The 29 MGnify genomes have unknown quality assessment "
                "method and are kept strictly separate."
            ),
            "masking_control_caveat": (
                "Across the 7 selected high-completeness, low-contamination MGnify reference genomes "
                "(numeric >=90% completeness, <=5% contamination; formal MIMAG status unknown due to unrecorded rRNA/tRNA criteria; "
                "9 nested masking conditions/replicates per genome = 63 evaluations), mean true masked GPR reaction recovery "
                "is 16.82% while alternative reactions account for 76.73% of newly induced additions. Gap-filling alone "
                "does not establish specific enzyme functions for dark proteins; these observations apply to the 7 tested "
                "reference genomes and are not generalized to all archaea."
            ),
        },
        "input_sha256": input_hashes,
        "part_a_case_universe": summary_a,
        "part_b_paired_interventions": summary_b,
        "part_c_cohort_descriptive": summary_c,
        "part_d_masking_7genomes": summary_d,
    }
    with report_path.open("w", encoding="utf-8") as fh:
        json.dump(report_payload, fh, indent=2, sort_keys=False)
        fh.write("\n")

    output_files = [
        combined_out_path,
        path_a,
        path_b,
        path_c1,
        path_c2,
        path_c3,
        path_c4,
        path_d,
        *fig_paths,
        report_path,
    ]
    manifest_rows = [
        {
            "artifact_role": "output",
            "filename": p.name,
            "bytes": p.stat().st_size,
            "sha256": sha256_file(p),
        }
        for p in output_files
    ]
    for code_name in (
        "build_exp1_descriptive_interpretation.py",
        "test_descriptive_interpretation.py",
        "README.md",
    ):
        cp = SCRIPT_DIR / code_name
        if cp.is_file():
            manifest_rows.append(
                {
                    "artifact_role": "code_or_doc",
                    "filename": cp.name,
                    "bytes": cp.stat().st_size,
                    "sha256": sha256_file(cp),
                }
            )
    for name, sh in sorted(input_hashes.items()):
        manifest_rows.append(
            {
                "artifact_role": "input_verified",
                "filename": name,
                "bytes": "",
                "sha256": sh,
            }
        )
    manifest_path = outdir / "exp1_interpretation_sha256_manifest.tsv"
    write_tsv(manifest_path, ["artifact_role", "filename", "bytes", "sha256"], manifest_rows)

    return report_payload


def _resolve_default_summary_source() -> tuple[Path | None, Path | None]:
    """Resolve portable default summary input when neither --shards-dir nor --combined-summary is explicitly given."""
    for cand in DEFAULT_REVIEW_SHARDS_CANDIDATES:
        if cand.is_dir() and all((cand / f"exp1_shard_0{i}_summary.tsv").is_file() for i in range(8)):
            return cand, None
    if DEFAULT_COMBINED_SUMMARY.is_file():
        return None, DEFAULT_COMBINED_SUMMARY
    return DEFAULT_REVIEW_SHARDS_CANDIDATES[0], None


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build reproducible descriptive interpretation tables, figures, and report for Exp1 484 models."
    )
    parser.add_argument(
        "--shards-dir",
        type=Path,
        default=None,
        help="Explicit directory containing exp1_shard_00..07_summary.tsv (fails closed if missing).",
    )
    parser.add_argument(
        "--combined-summary",
        type=Path,
        default=None,
        help="Explicit combined 12,580-row summary TSV path (fails closed if missing).",
    )
    parser.add_argument(
        "--covariate-supplement",
        type=Path,
        default=DEFAULT_COVARIATE_SUPPLEMENT,
        help="Path to corrected exp1_484_covariate_metadata_supplement.tsv.",
    )
    parser.add_argument(
        "--inventory",
        type=Path,
        default=DEFAULT_INVENTORY,
        help="Path to immutable exp1_model_inventory_and_qc_484.tsv.",
    )
    parser.add_argument(
        "--masking-cases",
        type=Path,
        default=DEFAULT_MASKING_CASES,
        help="Path to exp1_masking_controls_7hq_v4.tsv (63 rows).",
    )
    parser.add_argument(
        "--masking-baselines",
        type=Path,
        default=DEFAULT_MASKING_BASELINES,
        help="Path to exp1_masking_controls_7hq_v4_baselines.tsv (7 rows).",
    )
    parser.add_argument(
        "--root-paired-audit-json",
        type=Path,
        default=DEFAULT_ROOT_PAIRED_AUDIT_JSON,
        help="Path to root exp1_full484_paired_biomass_sensitivity.json for programmatic crosscheck.",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=DEFAULT_OUTDIR,
        help="Output directory for descriptive interpretation artifacts.",
    )
    parser.add_argument(
        "--no-strict-source-sha256",
        action="store_true",
        help="Disable production SHA-256 hash enforcement (for synthetic unit tests).",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    shards_dir = args.shards_dir
    combined_summary = args.combined_summary
    if shards_dir is None and combined_summary is None:
        shards_dir, combined_summary = _resolve_default_summary_source()

    run_interpretation_pipeline(
        shards_dir=shards_dir,
        combined_summary=combined_summary,
        covariate_supplement=args.covariate_supplement,
        inventory=args.inventory,
        masking_cases=args.masking_cases,
        masking_baselines=args.masking_baselines,
        root_paired_audit_json=args.root_paired_audit_json,
        outdir=args.outdir,
        strict_source_sha256=not args.no_strict_source_sha256,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
