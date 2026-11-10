#!/usr/bin/env python3
"""Unit, mutation, and acceptance tests for exp1_metabolic/interpretation."""

from __future__ import annotations

import copy
import json
from pathlib import Path

import pytest

from build_exp1_descriptive_interpretation import (
    B0,
    B1,
    B2,
    B3,
    B4,
    DEFAULT_COMBINED_SUMMARY,
    DEFAULT_COVARIATE_SUPPLEMENT,
    DEFAULT_INVENTORY,
    DEFAULT_MASKING_BASELINES,
    DEFAULT_MASKING_CASES,
    DEFAULT_OUTDIR,
    DEFAULT_ROOT_PAIRED_AUDIT_JSON,
    EXPECTED_24_CORE_COMBINATIONS,
    EXPECTED_OPTIONAL_S0_COMBINATIONS,
    S0,
    S1,
    S2,
    S3,
    S4,
    S5,
    S6,
    VERIFIED_EMPTY_S0_GENOMES,
    build_table_a_case_universe,
    build_table_b_paired_interventions,
    build_table_d_masking_7genomes,
    build_tables_c_cohort_descriptive,
    classify_completeness_bin,
    classify_contamination_bin,
    is_valid_summary_endpoint,
    load_summary_rows,
    read_tsv,
    run_interpretation_pipeline,
    validate_universe_and_join_covariates,
)


def _make_synthetic_dataset() -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    """Build a small synthetic 4-genome dataset (including 1 empty-S0 genome and 1 MGnify genome)."""
    cov_rows = [
        {
            "genome": "EQ_G1",
            "cohort": "EQ_desert",
            "source": "EQ_MAG",
            "canonical_order": "o__JACQPN01",
            "canonical_family": "f__JAIVGT01",
            "canonical_genus": "g__JAIVGT01",
            "shared_order_stratum": "True",
            "reported_quality_completeness": "92.0",
            "reported_quality_contamination": "2.0",
            "reported_quality_method": "CheckM2_v1.1.0",
            "method_family_checkm2_completeness": "92.0",
            "method_family_checkm2_contamination": "2.0",
            "checkm2_method_family_compatible": "True",
            "mimag_high_quality": "unknown",
        },
        {
            "genome": "S22PRr1__self_SemiBin_29",
            "cohort": "EQ_desert",
            "source": "EQ_MAG",
            "canonical_order": "o__JACQPN01",
            "canonical_family": "f__JAIVGT01",
            "canonical_genus": "g__JAIVGT01",
            "shared_order_stratum": "True",
            "reported_quality_completeness": "75.0",
            "reported_quality_contamination": "4.0",
            "reported_quality_method": "CheckM2_v1.1.0",
            "method_family_checkm2_completeness": "75.0",
            "method_family_checkm2_contamination": "4.0",
            "checkm2_method_family_compatible": "True",
            "mimag_high_quality": "unknown",
        },
        {
            "genome": "GTDB_G1",
            "cohort": "external_comparator",
            "source": "GTDB_r232_rep",
            "canonical_order": "o__JACQPN01",
            "canonical_family": "f__JAIVGT01",
            "canonical_genus": "g__JAIVGT01",
            "shared_order_stratum": "True",
            "reported_quality_completeness": "94.0",
            "reported_quality_contamination": "1.5",
            "reported_quality_method": "CheckM_v1_and_CheckM2_GTDB_r232",
            "method_family_checkm2_completeness": "94.0",
            "method_family_checkm2_contamination": "1.5",
            "checkm2_method_family_compatible": "True",
            "mimag_high_quality": "True",
        },
        {
            "genome": "MGYG_G1",
            "cohort": "external_comparator",
            "source": "MGnify",
            "canonical_order": "o__JACQPN01",
            "canonical_family": "f__JAIVGT01",
            "canonical_genus": "g__JAIVGT01",
            "shared_order_stratum": "True",
            "reported_quality_completeness": "95.0",
            "reported_quality_contamination": "1.0",
            "reported_quality_method": "unknown_mgnify_catalogue_pipeline",
            "method_family_checkm2_completeness": "",
            "method_family_checkm2_contamination": "",
            "checkm2_method_family_compatible": "False",
            "mimag_high_quality": "unknown",
        },
    ]

    summary_rows: list[dict[str, str]] = []
    for c in cov_rows:
        g = c["genome"]
        combos = list(EXPECTED_24_CORE_COMBINATIONS)
        if g not in VERIFIED_EMPTY_S0_GENOMES:
            combos = list(EXPECTED_OPTIONAL_S0_COMBINATIONS) + combos
        for scen, bio in combos:
            if scen in (S1, S2):
                status = "infeasible_under_fixed_medium"
                n_val = "0"
                min_add = ""
                res = ""
                bv = ""
                gr = ""
            else:
                status = "heuristic_incumbent_bounded"
                n_val = "5"
                base_rxn = 170.0 if c["cohort"] == "EQ_desert" else 180.0
                offset = 0.0 if bio == B0 else (-73.0 if bio == B2 else -84.0)
                min_add = str(base_rxn + offset)
                res = "1e-12"
                bv = "0.0"
                gr = "0.015"
            summary_rows.append(
                {
                    "genome": g,
                    "cohort": c["cohort"],
                    "canonical_order": c["canonical_order"],
                    "canonical_family": c["canonical_family"],
                    "canonical_genus": c["canonical_genus"],
                    "scenario_id": scen,
                    "biomass_recipe_id": bio,
                    "solver_status": status,
                    "n_distinct_valid_solutions": n_val,
                    "min_added_rxns_count": min_add,
                    "min_total_cost": "5000.0" if n_val != "0" else "",
                    "n_rxns_observed_in_all_enumerated_models": "80.0" if n_val != "0" else "",
                    "n_rxns_essential_in_declared_supernetwork": "45.0" if n_val != "0" else "",
                    "max_solution_residual": res,
                    "max_solution_bound_viol": bv,
                    "min_solution_verified_growth": gr,
                }
            )
    return summary_rows, cov_rows


def test_1_exact_case_universe_and_fail_closed_mutations() -> None:
    """Verify exact case universe rules: 24 core cases per genome, optional S0 only absent for verified IDs."""
    summary_rows, cov_rows = _make_synthetic_dataset()
    summary_by_key, cov_by_genome = validate_universe_and_join_covariates(summary_rows, cov_rows)
    assert len(cov_by_genome) == 4
    assert len(summary_by_key) == 3 * 26 + 1 * 24

    # Mutation 1: Remove one core case -> must fail closed
    bad_missing_core = [
        r for r in summary_rows if not (r["genome"] == "EQ_G1" and r["scenario_id"] == S3 and r["biomass_recipe_id"] == B0)
    ]
    with pytest.raises(ValueError, match="Missing core scenario/biomass pairs"):
        validate_universe_and_join_covariates(bad_missing_core, cov_rows)

    # Mutation 2: Duplicate a summary row -> must fail closed
    bad_dup = summary_rows + [copy.deepcopy(summary_rows[0])]
    with pytest.raises(ValueError, match="Duplicate summary case key"):
        validate_universe_and_join_covariates(bad_dup, cov_rows)

    # Mutation 3: Remove S0 from a genome that is NOT in VERIFIED_EMPTY_S0_GENOMES -> must fail closed
    bad_s0_missing = [
        r for r in summary_rows if not (r["genome"] == "EQ_G1" and r["scenario_id"] == S0)
    ]
    with pytest.raises(ValueError, match="missing optional S0 combinations"):
        validate_universe_and_join_covariates(bad_s0_missing, cov_rows)

    # Mutation 4: Add S0 to S22PRr1__self_SemiBin_29 (verified empty medium) -> must fail closed
    extra_s0 = copy.deepcopy(summary_rows[0])
    extra_s0["genome"] = "S22PRr1__self_SemiBin_29"
    extra_s0["scenario_id"] = S0
    extra_s0["biomass_recipe_id"] = B0
    with pytest.raises(ValueError, match="expected to have empty optional S0 medium"):
        validate_universe_and_join_covariates(summary_rows + [extra_s0], cov_rows)


def test_2_unknown_and_infeasible_never_treated_as_zero() -> None:
    """Verify unknown/timeout/verification-failed and fixed-medium-infeasible cases are never coerced to 0."""
    summary_rows, cov_rows = _make_synthetic_dataset()

    # Mutate EQ_G1 S3 B2 to be fba_verification_failed with a bogus "0.0" in min_added_rxns_count
    for r in summary_rows:
        if r["genome"] == "EQ_G1" and r["scenario_id"] == S3 and r["biomass_recipe_id"] == B2:
            r["solver_status"] = "fba_verification_failed"
            r["n_distinct_valid_solutions"] = "0"
            r["min_added_rxns_count"] = "0.0"

    assert not is_valid_summary_endpoint(
        next(r for r in summary_rows if r["genome"] == "EQ_G1" and r["scenario_id"] == S3 and r["biomass_recipe_id"] == B2)
    )

    summary_by_key, cov_by_genome = validate_universe_and_join_covariates(summary_rows, cov_rows)
    rows_b, _ = build_table_b_paired_interventions(summary_by_key, cov_by_genome, None)

    eq_s3_b0_b2 = next(
        r
        for r in rows_b
        if r["contrast_id"] == "biomass_B0_minus_B2"
        and r["cohort_scope"] == "EQ_desert"
        and r["fixed_context_id"] == S3
    )
    # Only S22PRr1__self_SemiBin_29 is usable; EQ_G1 is counted as numerical_unknown_pairs, NOT usable with 0!
    assert eq_s3_b0_b2["n_present_pairs"] == 2
    assert eq_s3_b0_b2["n_paired_usable"] == 1
    assert eq_s3_b0_b2["n_nonestimable_evaluated_pairs"] == 1
    assert eq_s3_b0_b2["n_numerical_unknown_pairs"] == 1
    assert eq_s3_b0_b2["n_fixedmedium_infeasible_pairs"] == 0
    assert float(eq_s3_b0_b2["diff_median_ref_minus_alt"]) == 73.0

    # Check S1 minimal medium: must be classified as fixedmedium_infeasible_pairs, never numerical_unknown_pairs
    eq_s1_b0_b2 = next(
        r
        for r in rows_b
        if r["contrast_id"] == "biomass_B0_minus_B2"
        and r["cohort_scope"] == "EQ_desert"
        and r["fixed_context_id"] == S1
    )
    assert eq_s1_b0_b2["n_paired_usable"] == 0
    assert eq_s1_b0_b2["n_fixedmedium_infeasible_pairs"] == 2
    assert eq_s1_b0_b2["n_numerical_unknown_pairs"] == 0

    # Check diagnostic violation on a supposedly valid row -> must fail closed
    bad_diag_row = copy.deepcopy(summary_rows[0])
    bad_diag_row["max_solution_residual"] = "1e-2"
    with pytest.raises(ValueError, match="Failed verified diagnostics"):
        is_valid_summary_endpoint(bad_diag_row, strict_diagnostics=True)


def test_3_paired_key_alignment_under_shuffled_rows() -> None:
    """Verify paired biomass and redox calculations align strictly by (genome, scenario_id, biomass_recipe_id)."""
    summary_rows, cov_rows = _make_synthetic_dataset()
    # Give EQ_G1 and S22PRr1 distinct B0 and B2 values so any misalignment changes the paired differences
    for r in summary_rows:
        if r["genome"] == "EQ_G1" and r["scenario_id"] == S3:
            if r["biomass_recipe_id"] == B0:
                r["min_added_rxns_count"] = "200.0"
            elif r["biomass_recipe_id"] == B2:
                r["min_added_rxns_count"] = "110.0"  # diff = 90.0
        if r["genome"] == "S22PRr1__self_SemiBin_29" and r["scenario_id"] == S3:
            if r["biomass_recipe_id"] == B0:
                r["min_added_rxns_count"] = "150.0"
            elif r["biomass_recipe_id"] == B2:
                r["min_added_rxns_count"] = "100.0"  # diff = 50.0

    summary_by_key_1, cov_by_genome_1 = validate_universe_and_join_covariates(summary_rows, cov_rows)
    rows_b_1, _ = build_table_b_paired_interventions(summary_by_key_1, cov_by_genome_1, None)

    reversed_rows = list(reversed(summary_rows))
    summary_by_key_2, cov_by_genome_2 = validate_universe_and_join_covariates(reversed_rows, list(reversed(cov_rows)))
    rows_b_2, _ = build_table_b_paired_interventions(summary_by_key_2, cov_by_genome_2, None)

    assert rows_b_1 == rows_b_2
    eq_s3 = next(
        r
        for r in rows_b_1
        if r["contrast_id"] == "biomass_B0_minus_B2"
        and r["cohort_scope"] == "EQ_desert"
        and r["fixed_context_id"] == S3
    )
    assert float(eq_s3["diff_median_ref_minus_alt"]) == 70.0  # median of 90.0 and 50.0
    assert float(eq_s3["diff_min"]) == 50.0
    assert float(eq_s3["diff_max"]) == 90.0


def test_4_unknown_quality_method_never_pooled_into_checkm2_strata() -> None:
    """Verify MGnify unknown quality method is never pooled into CheckM2 strata and illegal relabeling fails closed."""
    summary_rows, cov_rows = _make_synthetic_dataset()
    summary_by_key, cov_by_genome = validate_universe_and_join_covariates(summary_rows, cov_rows)
    comp_rows, strat_rows, mgnify_rows, _, summary_c = build_tables_c_cohort_descriptive(
        summary_by_key, cov_by_genome
    )

    assert summary_c["checkm2_method_family_compatible_total"] == 3  # 2 EQ + 1 GTDB; MGYG_G1 excluded
    assert summary_c["mgnify_unknown_method_separate_total"] == 1

    # In [90,101] x [0,5) CheckM2 stratum for o__JACQPN01, only GTDB_G1 (n=1) is external, NOT MGYG_G1!
    cell_90 = next(
        r
        for r in strat_rows
        if r["canonical_order"] == "o__JACQPN01"
        and r["checkm2_completeness_bin"] == "[90,101]"
        and r["checkm2_contamination_bin"] == "[0,5)"
        and r["scenario_id"] == S3
        and r["biomass_recipe_id"] == B0
    )
    assert cell_90["n_eq_valid"] == 1
    assert cell_90["n_ext_checkm2_valid"] == 1
    assert cell_90["is_sparse_n_lt_5"] == "True"

    # In [70,90) x [0,5) CheckM2 stratum, EQ has S22PRr1 (n=1) and external has 0 -> no contrast!
    cell_70 = next(
        r
        for r in strat_rows
        if r["canonical_order"] == "o__JACQPN01"
        and r["checkm2_completeness_bin"] == "[70,90)"
        and r["checkm2_contamination_bin"] == "[0,5)"
        and r["scenario_id"] == S3
        and r["biomass_recipe_id"] == B0
    )
    assert cell_70["n_eq_valid"] == 1
    assert cell_70["n_ext_checkm2_valid"] == 0
    assert cell_70["both_cohorts_valid_present"] == "False"
    assert cell_70["contrast_status"] == "no_contrast_missing_cohort"
    assert cell_70["descriptive_eq_minus_ext_median_diff"] == ""

    # Mutation: Attempt to mark MGnify unknown-method genome as checkm2_method_family_compatible=True -> fail closed
    bad_cov = copy.deepcopy(cov_rows)
    bad_cov[3]["checkm2_method_family_compatible"] = "True"
    with pytest.raises(ValueError, match="Unknown quality method genome"):
        validate_universe_and_join_covariates(summary_rows, bad_cov)


def test_5_completeness_and_contamination_bin_boundaries() -> None:
    """Verify exact completeness [50,70), [70,90), [90,101] and contamination [0,5), [5,10), boundary_10 rules."""
    assert classify_completeness_bin("") == "missing"
    assert classify_completeness_bin(None) == "missing"
    assert classify_completeness_bin("49.999") == "out_of_range"
    assert classify_completeness_bin("50.0") == "[50,70)"
    assert classify_completeness_bin("69.999") == "[50,70)"
    assert classify_completeness_bin("70.0") == "[70,90)"
    assert classify_completeness_bin("89.999") == "[70,90)"
    assert classify_completeness_bin("90.0") == "[90,101]"
    assert classify_completeness_bin("100.0") == "[90,101]"
    assert classify_completeness_bin("101.0") == "[90,101]"
    assert classify_completeness_bin("101.001") == "out_of_range"

    assert classify_contamination_bin("") == "missing"
    assert classify_contamination_bin(None) == "missing"
    assert classify_contamination_bin("-0.001") == "out_of_range"
    assert classify_contamination_bin("0.0") == "[0,5)"
    assert classify_contamination_bin("4.999") == "[0,5)"
    assert classify_contamination_bin("5.0") == "[5,10)"
    assert classify_contamination_bin("9.999") == "[5,10)"
    assert classify_contamination_bin("10.0") == "boundary_10"
    assert classify_contamination_bin("10.001") == "out_of_range"


def test_6_explicit_bad_paths_fail_closed_no_fallback(tmp_path: Path) -> None:
    """Verify that bad explicit CLI paths immediately raise FileNotFoundError without silent fallback."""
    nonexistent = tmp_path / "does_not_exist.tsv"
    with pytest.raises(FileNotFoundError, match="Explicit --combined-summary"):
        load_summary_rows(shards_dir=None, combined_summary=nonexistent, strict_source_sha256=False)

    with pytest.raises(FileNotFoundError, match="Explicit --shards-dir"):
        load_summary_rows(shards_dir=tmp_path / "missing_shards", combined_summary=None, strict_source_sha256=False)

    with pytest.raises(FileNotFoundError, match="Explicit --covariate-supplement"):
        run_interpretation_pipeline(
            shards_dir=None,
            combined_summary=DEFAULT_COMBINED_SUMMARY,
            covariate_supplement=nonexistent,
            inventory=DEFAULT_INVENTORY,
            masking_cases=DEFAULT_MASKING_CASES,
            masking_baselines=DEFAULT_MASKING_BASELINES,
            root_paired_audit_json=DEFAULT_ROOT_PAIRED_AUDIT_JSON,
            outdir=tmp_path / "out",
            strict_source_sha256=False,
        )


def test_7_production_artifacts_and_portable_combined_summary_parity(tmp_path: Path) -> None:
    """Verify production report numbers, 7-genome masking hierarchy, and portable combined-summary parity."""
    report_path = DEFAULT_OUTDIR / "exp1_descriptive_interpretation_report.json"
    assert report_path.is_file()
    with report_path.open("r", encoding="utf-8") as fh:
        rep = json.load(fh)

    # Part A checks
    pa = rep["part_a_case_universe"]
    assert pa["n_genomes"] == 484
    assert pa["n_planned_core_cases"] == 11616
    assert pa["n_optional_s0_present"] == 964
    assert pa["n_optional_s0_absent_empty_medium"] == 4
    assert pa["n_total_evaluated_cases"] == 12580
    assert pa["solver_status_counts"] == {
        "heuristic_incumbent_bounded": 10578,
        "infeasible_under_fixed_medium": 1936,
        "fba_verification_failed": 56,
        "timeout_no_incumbent": 10,
    }
    assert pa["n_total_numerical_unknown_cases"] == 66

    # Part B checks: nonestimable_pairs split into fixedmedium_infeasible_pairs vs numerical_unknown_pairs
    pb = rep["part_b_paired_interventions"]
    assert pb["root_paired_biomass_sensitivity_crosscheck_passed"] is True
    s1_eq = pb["b0_minus_b2_by_cohort_and_scenario"]["EQ_desert|S1_anoxic_minimal_h2_co2_nh4"]
    assert s1_eq["nonestimable_pairs"] == 237
    assert s1_eq["fixedmedium_infeasible_pairs"] == 237
    assert s1_eq["numerical_unknown_pairs"] == 0
    assert "unresolved_pairs" not in s1_eq

    s0_eq = pb["b0_minus_b2_by_cohort_and_scenario"]["EQ_desert|S0_gapseq_native_predicted"]
    assert s0_eq["cohort_genomes"] == 237
    assert s0_eq["present_pairs"] == 236
    assert s0_eq["optionalmedium_unavailable_pairs"] == 1
    assert s0_eq["paired_usable"] == 235
    assert s0_eq["nonestimable_pairs"] == 1
    assert s0_eq["fixedmedium_infeasible_pairs"] == 0
    assert s0_eq["numerical_unknown_pairs"] == 1

    # Part C checks
    pc = rep["part_c_cohort_descriptive"]
    assert pc["checkm2_method_family_compatible_total"] == 455
    assert pc["checkm2_shared_orders_total"] == 293
    assert pc["mgnify_unknown_method_separate_total"] == 29
    assert pc["n_checkm2_cells_with_both_cohorts"] == 9
    assert pc["n_checkm2_cells_non_sparse_n_ge_5"] == 5

    # Part D checks: 7 reference genomes, 9 nested replicates per genome = 63
    pd = rep["part_d_masking_7genomes"]
    assert pd["n_independent_control_genomes"] == 7
    assert pd["n_nested_replicates_per_genome"] == 9
    assert pd["n_total_masked_evaluations"] == 63
    assert abs(pd["mean_of_7_genome_mean_recovery_rates"] - 0.1682162622416794) < 1e-12
    assert abs(pd["mean_of_7_genome_mean_alternative_fractions"] - 0.7673003790562458) < 1e-12
    assert "formal MIMAG status unknown" in pd["control_genome_selection_note"]

    # Verify portable execution from --combined-summary produces identical TSV tables
    portable_out = tmp_path / "portable_run"
    rep_portable = run_interpretation_pipeline(
        shards_dir=None,
        combined_summary=DEFAULT_COMBINED_SUMMARY,
        covariate_supplement=DEFAULT_COVARIATE_SUPPLEMENT,
        inventory=DEFAULT_INVENTORY,
        masking_cases=DEFAULT_MASKING_CASES,
        masking_baselines=DEFAULT_MASKING_BASELINES,
        root_paired_audit_json=DEFAULT_ROOT_PAIRED_AUDIT_JSON,
        outdir=portable_out,
        strict_source_sha256=True,
    )
    assert rep_portable["part_a_case_universe"] == rep["part_a_case_universe"]
    assert rep_portable["part_b_paired_interventions"] == rep["part_b_paired_interventions"]
    assert rep_portable["part_c_cohort_descriptive"] == rep["part_c_cohort_descriptive"]
    assert rep_portable["part_d_masking_7genomes"] == rep["part_d_masking_7genomes"]

    for tsv_name in (
        "exp1_case_universe_and_solver_status_summary.tsv",
        "exp1_paired_biomass_and_redox_sensitivity.tsv",
        "exp1_order_and_quality_composition_overlap.tsv",
        "exp1_checkm2_stratified_cohort_comparisons.tsv",
        "exp1_mgnify29_reported_quality_sensitivity.tsv",
        "exp1_unadjusted_all484_descriptive_summary.tsv",
        "exp1_masking_7genomes_descriptive_summary.tsv",
    ):
        assert read_tsv(portable_out / tsv_name) == read_tsv(DEFAULT_OUTDIR / tsv_name)
