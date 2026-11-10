#!/usr/bin/env python3
"""8-Gate Standalone Reproducibility & Independent Replay Verifier for Stage 4 v4 (`exp3_viruses/stats_v4`)."""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "code"))

from stage4_statistical_audit_v4 import (  # noqa: E402
    DESIGN_LEAD_FAMILIES,
    bh_adjust,
    build_design_matrix_audit,
    build_family_assay_maps,
    compute_abundance_metrics,
    compute_sha256,
    conditional_family_tests,
    evaluate_environment_matched_associations,
    find_matched_alternative_controls,
    join_candidate_taxonomy,
    load_and_orient_abundance_matrices,
    load_field_map,
    load_metadata,
    load_provisional_candidate,
    load_root_unresolved_targets,
    partition_detection_opportunity,
    resolve_graded_candidate_pairs,
    select_primary_environment_design,
    spearman_permutation_result,
    stratified_permutation_indices,
    validate_assay_provenance,
)


def assert_true(cond: bool, msg: str) -> None:
    if not cond:
        raise AssertionError(msg)


def verify_all(root: Path = ROOT) -> dict:
    res_dir = root / "results"
    in_dir = root / "inputs"
    ind_dir = root / "independent_audit"
    man_dir = root / "manifests"

    # [1/8] SHA256 Manifest & Remote Output Manifest Integrity
    sha_manifest = pd.read_csv(man_dir / "sha256_manifest.tsv", sep="\t")
    assert_true(len(sha_manifest) >= 40, f"Expected >=40 entries in sha256_manifest.tsv, got {len(sha_manifest)}")
    for row in sha_manifest.itertuples(index=False):
        target_path = root / row.relative_path
        assert_true(target_path.is_file(), f"Missing packaged file: {row.relative_path}")
        assert_true(target_path.stat().st_size == int(row.bytes), f"Byte size mismatch for {row.relative_path}")
        assert_true(compute_sha256(target_path) == row.sha256, f"SHA256 mismatch for {row.relative_path}")

    out_manifest = pd.read_csv(res_dir / "stage4_v4_output_manifest.tsv", sep="\t")
    assert_true(len(out_manifest) == 20, f"Expected 20 files in stage4_v4_output_manifest.tsv, got {len(out_manifest)}")
    for row in out_manifest.itertuples(index=False):
        p = res_dir / row.file
        assert_true(p.is_file(), f"Missing output file: {row.file}")
        assert_true(p.stat().st_size == int(row.bytes), f"Output size mismatch: {row.file}")
        assert_true(compute_sha256(p) == row.sha256, f"Output SHA256 mismatch: {row.file}")

    code_sha = compute_sha256(root / "code/stage4_statistical_audit_v4.py")
    assert_true(
        code_sha == "f2f17db56a68b195326b492951df00ff4fd0a466b621ac4a792b28d2f083e63a",
        f"Unexpected stage4_statistical_audit_v4.py SHA256: {code_sha}",
    )
    print("[1/8] PASS: SHA256 package manifest, 20-file remote output manifest, and executed code hash verified.")

    # [2/8] Candidate Pair Ledger Derivation & Individual Root Dispositions
    pairs_raw = pd.read_csv(res_dir / "stage3_archaeal_virus_host_pairs_all_tiers_v2.tsv", sep="\t")
    taxonomy = pd.read_csv(in_dir / "taxonomy_reconciliation_all5438.tsv.gz", sep="\t")
    unresolved_targets = load_root_unresolved_targets(res_dir / "stage4_v4_root_disposition_41Sr3.json")
    provisional_key, provisional_payload = load_provisional_candidate(
        res_dir / "stage4_v4_provisional_disposition_second_thermo.json"
    )
    ledger, primary, provisional = resolve_graded_candidate_pairs(pairs_raw, unresolved_targets, {provisional_key})
    primary = join_candidate_taxonomy(primary, taxonomy)
    provisional = join_candidate_taxonomy(provisional, taxonomy)

    assert_true(len(ledger) == 2121, f"Expected 2121 raw ledger rows, got {len(ledger)}")
    assert_true(int(ledger["evidence_tier"].eq("Tier1_verified_direct").sum()) == 8, "Expected 8 Tier1 rows")
    assert_true(int(ledger["is_alias_duplicate"].sum()) == 1, "Expected 1 alias duplicate row")
    assert_true(int(ledger["is_unresolved_root_disposition"].sum()) == 1, "Expected 1 unresolved root disposition row")
    assert_true(int(ledger["is_boundary_sensitive_provisional"].sum()) == 1, "Expected 1 boundary-sensitive provisional row")
    assert_true(len(primary) == 5 and int(primary["canonical_591_host_mag"].nunique()) == 4, "Expected 5 primary pairs across 4 hosts")
    assert_true(len(provisional) == 6 and int(provisional["canonical_591_host_mag"].nunique()) == 5, "Expected 6 provisional pairs across 5 hosts")
    pd.testing.assert_frame_equal(ledger, pd.read_csv(res_dir / "stage4_v4_candidate_pairs_graded_ledger.tsv", sep="\t"), check_dtype=False)
    print("[2/8] PASS: Graded candidate pair ledger (2,121 rows; 5 primary pairs/4 hosts; 6 provisional pairs/5 hosts) verified.")

    # [3/8] Detection Opportunity Partition (447 zero-spacer untested vs 144 eligible)
    screening = pd.read_csv(res_dir / "stage1_mag_screening_status.tsv", sep="\t")
    evidence = pd.read_csv(res_dir / "stage3_all591_archaeal_mag_host_evidence_v2.tsv", sep="\t")
    partition = partition_detection_opportunity(screening, evidence, taxonomy, primary, provisional)
    assert_true(len(partition) == 591, "Expected 591 canonical archaeal MAGs")
    assert_true(int(partition["n_valid_spacers"].eq(0).sum()) == 447, "Expected 447 zero-spacer untested MAGs")
    assert_true(int(partition["is_screen_positive_eligible"].sum()) == 144, "Expected 144 screen-positive eligible MAGs")
    assert_true(partition.loc[partition["n_valid_spacers"].eq(0), "primary_DETECTEDLINK"].isna().all(), "Zero-spacer outcomes must be NA")
    assert_true(int(partition["is_primary_detected_link_host"].sum()) == 4, "Expected 4 primary detected hosts")
    assert_true(int(partition["is_provisional_detected_link_host"].sum()) == 5, "Expected 5 provisional detected hosts")
    print("[3/8] PASS: Detection opportunity partition (591 = 447 zero-spacer untested + 144 eligible) verified.")

    # [4/8] Complete Protein-Family Assignment Scope (237 Thermoplasmatota, 44/47 families, 3,184 gate-false retained)
    members = pd.read_csv(in_dir / "exp4_47_member_evidence.tsv", sep="\t")
    full_loci = pd.read_csv(in_dir / "task_a_stage1b_all_candidate_member_loci_verified.tsv.gz", sep="\t")
    thermo237 = pd.read_csv(in_dir / "taxonomy_reconciliation_thermo237.tsv", sep="\t")
    assert_true(set(full_loci["genome"].astype(str)) == set(thermo237["genome"].astype(str)), "Full-locus genome set must equal exact 237 Thermoplasmatota set")
    assay_prov = validate_assay_provenance(full_loci, in_dir / "stage1_summary.json", in_dir / "stage1b_summary.json")
    assert_true(assay_prov["n_screened_genomes"] == 237, "Expected 237 screened genomes")
    assert_true(assay_prov["n_original_proteins_in_stage1_scope"] == 777177, "Expected 777,177 proteins in Stage 1 scope")
    assert_true(assay_prov["n_candidate_member_locus_rows_complete_membership"] == 41494, "Expected 41,494 candidate member locus rows")
    assert_true(assay_prov["n_final_promotion_gate_false_rows_retained_as_membership_evidence"] == 3184, "Expected 3,184 gate-false rows retained")

    fam_res, fam_cov = conditional_family_tests(partition, members, full_loci, assay_prov)
    prim_fam = fam_res.loc[fam_res["outcome_set"].eq("primary")]
    assert_true(len(prim_fam) == 47, "Expected 47 primary family rows")
    assert_true(int(prim_fam["fisher_exact_pvalue"].notna().sum()) == 44, "Expected 44 estimable families")
    assert_true(int(prim_fam["fisher_exact_pvalue"].isna().sum()) == 3, "Expected 3 dossier-only non-estimable families")
    assert_true(prim_fam.loc[prim_fam["fisher_exact_pvalue"].notna(), "n_assay_covered_eligible"].eq(58).all(), "Expected 58 eligible covered genomes")
    assert_true(prim_fam.loc[prim_fam["fisher_exact_pvalue"].notna(), "n_assay_covered_detected_links"].eq(1).all(), "Expected 1 primary event in 58 covered genomes")
    print("[4/8] PASS: Complete 237-genome protein-family assignment scope (44/47 estimable, 3 NA, 3,184 gate-false retained) verified.")

    # [5/8] Design Matrix Rank, Separation LP Witness, & Zero Imputation
    elig_raw, elig_audit = build_design_matrix_audit(
        partition, members, full_loci, "primary_DETECTEDLINK", "screen_positive_eligible_primary", assay_prov
    )
    assert_true(elig_audit["raw_matrix_rows"] == 144, "Expected 144 raw eligible matrix rows")
    assert_true(elig_audit["complete_case_rows_no_imputation"] == 58, "Expected 58 complete-case rows")
    assert_true(elig_audit["complete_case_detected_link_events"] == 1, "Expected 1 complete-case primary event")
    assert_true(elig_audit["design_matrix_rank"] == 5 and elig_audit["n_parameters_including_intercept"] == 9, "Expected rank 5 of 9")
    assert_true(elig_audit["rank_deficient"] is True, "Expected rank_deficient=True")
    assert_true(elig_audit["separation_LP"]["quasi_complete"] is True, "Expected LP quasi-complete separation witness")
    assert_true(elig_audit["model_fit_performed"] is False, "Overparameterized model must not be fit")
    print("[5/8] PASS: Design matrix diagnostics (144 raw -> 58 complete cases, 1 event, rank 5/9, LP witness, no model fit) verified.")

    # [6/8] Matched Alternative MAG Control Selection Replay from Compact 591-MAG Subset
    metadata = load_metadata(in_dir / "metadata.tsv")
    mag_sub_df, vir_sub_df = load_and_orient_abundance_matrices(
        in_dir / "mag_abundance_arch591_441samples.tsv.gz",
        in_dir / "viral_tpm_candidate_targets_441samples.tsv.gz",
        set(metadata["sample_id"]),
    )
    shared_samples = [s for s in metadata["sample_id"] if s in mag_sub_df.columns and s in vir_sub_df.columns]
    assert_true(len(shared_samples) == 441, f"Expected 441 shared samples, got {len(shared_samples)}")

    metrics = compute_abundance_metrics(mag_sub_df, partition["genome"], shared_samples)
    prod_metrics = pd.read_csv(res_dir / "stage4_v4_host_abundance_assay_metrics.tsv", sep="\t")
    pd.testing.assert_frame_equal(metrics, prod_metrics, check_exact=False, rtol=1e-12, atol=1e-15)
    assert_true(int((metrics["abundance_assay_status"] == "observed_matrix_row").sum()) == 580, "Expected 580 observed MAG rows")
    assert_true(int((metrics["abundance_assay_status"] == "missing_matrix_row_not_zero").sum()) == 11, "Expected 11 missing MAG rows retained as NA")

    unresolved_hosts = set(ledger.loc[ledger["is_unresolved_root_disposition"].eq(1), "canonical_591_host_mag"])
    controls = find_matched_alternative_controls(
        set(provisional["canonical_591_host_mag"]),
        set(primary["canonical_591_host_mag"]),
        set(provisional["canonical_591_host_mag"]),
        unresolved_hosts,
        partition,
        metrics,
    )
    prod_controls = pd.read_csv(res_dir / "stage4_v4_matched_alternative_MAG_controls.tsv", sep="\t")
    pd.testing.assert_frame_equal(controls, prod_controls, check_exact=False, rtol=1e-12, atol=1e-15)
    print("[6/8] PASS: Compact 591-MAG abundance metrics (580 observed, 11 NA) and 5 matched alternative MAG controls verified.")

    # [7/8] Environment Selection (trip5 x PlantRoot, 47 sites) & Independent 47-Site Permutation Replay
    field_map = load_field_map(in_dir / "matched_profiles.tsv")
    group_counts, selected, selection = select_primary_environment_design(metadata, field_map, set(shared_samples))
    assert_true(selection["selected_trip"] == "trip5" and selection["selected_compartment"] == "PlantRoot", "Expected trip5 x PlantRoot")
    assert_true(selection["n_selected_sites"] == 47 and selection["n_selected_samples"] == 47, "Expected 47 sites/samples")
    assert_true(selection["outcome_values_used_for_selection"] is False, "Selection must be outcome-independent")

    abund, diff, env_design = evaluate_environment_matched_associations(
        primary, provisional, controls, mag_sub_df, vir_sub_df, selected, 9999, 20261109
    )
    prod_abund = pd.read_csv(res_dir / "stage4_v4_environment_matched_abundance.tsv", sep="\t")
    prod_diff = pd.read_csv(res_dir / "stage4_v4_candidate_control_correlation_differences.tsv", sep="\t")
    pd.testing.assert_frame_equal(abund, prod_abund, check_exact=False, rtol=1e-12, atol=1e-15)
    pd.testing.assert_frame_equal(diff, prod_diff, check_exact=False, rtol=1e-12, atol=1e-15)
    assert_true("difference_permutation_pvalue" not in diff.columns, "Correlation difference table must be descriptive only (no p-values)")

    # Direct replay from independent_audit/stage4_v4_independent_audit_selected47_raw_site_observations.tsv
    raw47 = pd.read_csv(ind_dir / "stage4_v4_independent_audit_selected47_raw_site_observations.tsv", sep="\t")
    assert_true(len(raw47) == 282, f"Expected 282 rows (6 pairs x 47 sites), got {len(raw47)}")
    focal = raw47.loc[
        raw47["candidate_host_mag"].eq("41PRr2__single_SemiBin_236")
        & raw47["target_id"].eq("V60PRr2::V60PRr2_715976")
    ].sort_values("site")
    assert_true(len(focal) == 47, "Expected 47 sites for focal V60 pair")
    perms = stratified_permutation_indices(focal["dune_field"].to_numpy(dtype=str), 9999, 20261109)
    focal_res = spearman_permutation_result(
        focal["candidate_host_site_mean_abundance"].to_numpy(dtype=float),
        focal["viral_votu_site_mean_tpm"].to_numpy(dtype=float),
        perms,
    )
    assert_true(np.isclose(focal_res["observed_rho"], 0.25415577869613154, atol=1e-12), "Focal rho mismatch")
    assert_true(np.isclose(focal_res["pvalue"], 0.1692, atol=1e-12), "Focal p-value mismatch")
    assert_true(focal_res["n_extreme_including_ties"] == 1691, "Focal n_extreme mismatch")
    print("[7/8] PASS: Environment design (trip5 x PlantRoot, 47 sites) and independent Spearman + within-field permutation replay verified.")

    # [8/8] Slurm Execution Receipts & Extraction Provenance
    receipts = pd.read_csv(man_dir / "slurm_execution_receipts.tsv", sep="\t")
    job_ids = set(receipts["job_id"].astype(str))
    assert_true({"15894451", "15894468", "15894666"} <= job_ids, f"Missing expected Slurm job IDs: {job_ids}")
    ext_receipt = json.loads((ind_dir / "stage4_v4_independent_audit_extraction_receipt.json").read_text())
    assert_true(ext_receipt["slurm_job_id"] == "15894666", "Expected extraction job 15894666")
    assert_true(ext_receipt["production_v4_job_id_preserved"] == "15894468", "Expected preserved production job 15894468")
    print("[8/8] PASS: Slurm execution receipts (15894451, 15894468, 15894666) and source extraction provenance verified.")
    print("ALL 8 STAGE 4 v4 VERIFICATION GATES PASSED.")
    return {"status": "PASS", "gates_passed": 8}


if __name__ == "__main__":
    verify_all()
