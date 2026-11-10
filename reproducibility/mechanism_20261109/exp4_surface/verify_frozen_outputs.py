#!/usr/bin/env python3
"""Standalone Frozen-Output Reproducibility Verifier for Experiment 4 (`exp4_surface`).

Lightweight verification suite executing in ~2-3 seconds locally without cluster compute or GPU.
Validates:
1. Exact SHA256 checksums of all package artifacts against manifests/sha256_manifest.tsv.
2. Global model availability partition: 565 exact, 802 missing, 93 mismatch (sum = 1,460).
3. Full-universe genus joint BH-FDR correction: 8,760 tests, min q = 1.0000, 0 significant at q <= 0.10.
4. Candidate dossiers (47 families): honest deferred counts (0 aligned), no biochemical assignments.
5. Member-level evidence (428 members): 326 verified interior, 32 exact models, 27 verified interior with exact model.
6. Flank-level evidence (2,883 loci): annotated guild support semantics, unknown confidence markers.
7. Coordinate audit: 24/24 lead context pairs immediately adjacent (offset +/-1), cooriented, tight intergenic gaps [-4, 37] bp.
8. Lead scientific dispositions: bounded narratives, no overstatements, classified as uncertain_insufficient_evidence.
9. Independent replay receipts: 100% agreement across 428 members and byte-exact MSA sequence identity.
10. Slurm execution receipts and raw scheduler accounting: owner validation (`hohndor`), allowed job names, exact v2 array tasks 0..16, and required stage correspondence.
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import pandas as pd

PKG_ROOT = Path(__file__).resolve().parent
MANIFEST_P = PKG_ROOT / "manifests" / "sha256_manifest.tsv"
COMPACT_DIR = PKG_ROOT / "compact_evidence"


def sha256_file(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        while chunk := fh.read(65536):
            h.update(chunk)
    return h.hexdigest()


def test_sha256_checksums() -> None:
    print("[1/10] Checking SHA256 manifest integrity...")
    assert MANIFEST_P.exists(), f"Missing manifest: {MANIFEST_P}"
    manifest_df = pd.read_csv(MANIFEST_P, sep="\t")
    assert len(manifest_df) >= 30, f"Expected at least 30 checksum entries, got {len(manifest_df)}"

    for _, row in manifest_df.iterrows():
        expected_sha = row["sha256"]
        rel_p = row["relative_path"]
        target_p = PKG_ROOT / rel_p
        assert target_p.exists(), f"Manifest file missing on disk: {target_p}"
        actual_sha = sha256_file(target_p)
        assert actual_sha == expected_sha, (
            f"Checksum mismatch for {rel_p}:\n  expected: {expected_sha}\n  actual:   {actual_sha}"
        )
    print(f"  PASS: All {len(manifest_df)} package artifacts verified byte-exact.")


def test_global_summary_and_fdr() -> None:
    print("[2/10] Recomputing full-universe FDR and model availability from raw evaluation table...")
    eval_p = COMPACT_DIR / "exp4_surface_family_evaluation_v2_audited.tsv"
    frozen_p = COMPACT_DIR / "frozen_family_universe.tsv"
    summary_p = COMPACT_DIR / "exp4_surface_final_summary_v2_audited.json"
    audit_proof_p = COMPACT_DIR / "exp4_independent_polymer_and_msa_replay_audit.json"

    assert eval_p.exists(), f"Evaluation table missing: {eval_p}"
    assert frozen_p.exists(), f"Frozen universe table missing: {frozen_p}"
    df = pd.read_csv(eval_p, sep="\t")
    frozen_df = pd.read_csv(frozen_p, sep="\t")

    assert len(df) == 1460, f"Expected 1460 families in evaluation table, got {len(df)}"
    assert df["family_id"].nunique() == 1460, "Duplicate family IDs found in evaluation table"
    assert len(frozen_df) == 1460, f"Expected 1460 families in frozen universe, got {len(frozen_df)}"
    assert frozen_df["family_id"].nunique() == 1460, "Duplicate family IDs found in frozen universe"

    # Exact ID set match (ensures any wrong ID substitution fails immediately)
    eval_fids = set(df["family_id"])
    frozen_fids = set(frozen_df["family_id"])
    assert eval_fids == frozen_fids, (
        f"Family ID set mismatch between evaluation table and frozen universe:\n"
        f"  Missing from evaluation: {len(frozen_fids - eval_fids)}\n"
        f"  Unexpected in evaluation: {len(eval_fids - frozen_fids)}"
    )

    # Exact cohort group match per family between evaluation and frozen universe
    merged_cohorts = df[["family_id", "cohort_group"]].merge(
        frozen_df[["family_id", "cohort_group"]],
        on="family_id",
        suffixes=("_eval", "_frozen"),
    )
    assert len(merged_cohorts) == 1460
    assert (merged_cohorts["cohort_group_eval"] == merged_cohorts["cohort_group_frozen"]).all(), (
        "Cohort group mismatch for one or more families between evaluation and frozen universe"
    )

    null_guilds = [
        "nutrient_acquisition_peptidase",
        "substrate_binding_or_uptake",
        "surface_assembly_pilin_archaellin",
        "cell_coat_or_envelope",
        "host_defense_antiviral",
        "viral_or_mobile_element_origin",
    ]

    # Recompute genus-dereplicated joint BH-FDR across all 8,760 tests from raw p-values
    import numpy as np
    raw_genus_p: list[float] = []
    raw_primary_p: list[float] = []
    test_keys: list[tuple[str, str]] = []

    for _, r in df.iterrows():
        fid = str(r["family_id"])
        for g in null_guilds:
            # Genus p: missing/NaN treated conservatively as 1.0
            g_val = r.get(f"genus_derep_p_exact_{g}")
            gp = float(g_val) if pd.notna(g_val) else 1.0
            raw_genus_p.append(max(0.0, min(1.0, gp)))

            # Primary p
            p_val = r.get(f"null_p_exact_{g}")
            pp = float(p_val) if pd.notna(p_val) else 1.0
            raw_primary_p.append(max(0.0, min(1.0, pp)))
            test_keys.append((fid, g))

    assert len(raw_genus_p) == 8760, f"Expected 8760 genus tests, got {len(raw_genus_p)}"
    assert len(raw_primary_p) == 8760, f"Expected 8760 primary tests, got {len(raw_primary_p)}"

    # Direct BH-FDR implementation (no dependencies)
    def compute_bh(p_list: list[float]) -> np.ndarray:
        m = len(p_list)
        p_arr = np.array(p_list, dtype=float)
        order = np.argsort(p_arr)
        p_sorted = p_arr[order]
        q_sorted = np.zeros(m)
        min_q = 1.0
        for i in range(m - 1, -1, -1):
            rank = i + 1
            val = (p_sorted[i] * m) / rank
            min_q = min(min_q, val)
            q_sorted[i] = min(1.0, min_q)
        q_orig = np.empty(m)
        q_orig[order] = q_sorted
        return q_orig

    recomp_genus_q = compute_bh(raw_genus_p)
    recomp_primary_q = compute_bh(raw_primary_p)

    # Recomputed Genus Sensitivity: min q == 1.0000, 0 significant
    min_genus_q = float(np.min(recomp_genus_q))
    n_genus_sig_10 = int(np.sum(recomp_genus_q <= 0.10))
    n_genus_sig_05 = int(np.sum(recomp_genus_q <= 0.05))
    assert min_genus_q >= 1.0 - 1e-9, f"Recomputed min genus q must be >= 1.0, got {min_genus_q}"
    assert n_genus_sig_10 == 0, f"Expected 0 genus-significant tests at q<=0.10, got {n_genus_sig_10}"
    assert n_genus_sig_05 == 0

    # Recomputed Primary Joint FDR: exactly 2 leads significant at q <= 0.05
    primary_sig_05 = [test_keys[i] for i, q in enumerate(recomp_primary_q) if q <= 0.05]
    assert len(primary_sig_05) == 2, f"Expected exactly 2 primary leads at q<=0.05, got {len(primary_sig_05)}: {primary_sig_05}"
    assert ("22PRr1::22PRr1_2531955_3", "surface_assembly_pilin_archaellin") in primary_sig_05
    assert ("V60PRr2::V60PRr2_1314979_4", "cell_coat_or_envelope") in primary_sig_05

    # Recompute model availability directly from rows
    exact_recomp = int((df["model_polymer_hash_status"].isin(["exact_match", "exact_sequence_hash_match"])).sum())
    missing_recomp = int((df["model_polymer_hash_status"] == "missing_prediction_not_negative").sum())
    mismatch_recomp = int(
        (
            df["model_polymer_hash_status"].str.startswith("invalid_sequence_mismatch")
            | df["model_polymer_hash_status"].isin(["mismatch", "sequence_mismatch"])
        ).sum()
    )
    assert exact_recomp == 565, f"Expected 565 exact models, got {exact_recomp}"
    assert missing_recomp == 802, f"Expected 802 missing models, got {missing_recomp}"
    assert mismatch_recomp == 93, f"Expected 93 mismatch models, got {mismatch_recomp}"
    assert exact_recomp + missing_recomp + mismatch_recomp == 1460

    # Validate against summary JSON
    with open(summary_p) as fh:
        summary = json.load(fh)
    assert summary["n_families_total"] == 1460
    assert summary["n_tests_total"] == 8760
    fdr_s = summary["full_universe_genus_joint_fdr"]
    assert fdr_s["min_genus_joint_q"] == min_genus_q
    assert fdr_s["n_genus_joint_significant_q_le_0_10"] == n_genus_sig_10
    part_s = summary["model_availability_partition_global"]
    assert part_s["exact_sequence_hash_match"] == exact_recomp
    assert part_s["missing_prediction_not_negative"] == missing_recomp
    assert part_s["invalid_sequence_mismatch"] == mismatch_recomp

    print(f"  PASS: Recomputed 8,760 genus FDR tests (min q={min_genus_q:.4f}, 0 sig), 2 primary leads, and 565/802/93 partition directly from raw evaluation rows.")


def test_cohort_architecture_comparison() -> None:
    print("[3/10] Checking cohort comparison table...")
    cohort_p = COMPACT_DIR / "exp4_cohort_architecture_and_role_comparison_v2_audited.tsv"
    df = pd.read_csv(cohort_p, sep="\t")

    assert len(df) == 5
    assert df["n_families"].sum() == 1460
    assert df["n_models_exact_match"].sum() == 565
    assert df["n_models_missing"].sum() == 802
    assert df["n_models_mismatch"].sum() == 93

    for _, row in df.iterrows():
        total = row["n_families"]
        s = row["n_models_exact_match"] + row["n_models_missing"] + row["n_models_mismatch"]
        assert s == total, f"Cohort {row['cohort_group']} sum mismatch: {s} != {total}"
    print("  PASS: All 5 cohort groups partition exactly to 1,460 families.")


def test_candidate_dossiers() -> None:
    print("[4/10] Checking 47 candidate dossiers...")
    dossiers_p = COMPACT_DIR / "exp4_47_candidate_dossiers.tsv"
    df = pd.read_csv(dossiers_p, sep="\t")

    assert len(df) == 47
    # 29 completed, 18 deferred
    assert (df["msa_status"] == "mafft_completed").sum() == 29
    deferred = df[df["msa_status"] == "deferred_non_focal_no_role_hit"]
    assert len(deferred) == 18
    # Honest sequence count: deferred must have 0 aligned sequences
    assert (deferred["n_aligned_sequences"] == 0).all()
    assert (deferred["n_member_sequences"] > 0).all()

    # No biochemical assignments
    assert (df["is_biochemical_assignment"] == False).all()
    assert (df["candidate_computational_role_hypothesis"].isin(["none", "surface_assembly_pilin_archaellin", "cell_coat_or_envelope", "nutrient_acquisition_peptidase", "substrate_binding_or_uptake", "host_defense_antiviral", "viral_or_mobile_element_origin", "exploratory_fold_match"])).all()
    print("  PASS: 47 dossiers verified (29 completed MSAs, 18 honest deferred with 0 aligned, 0 biochemical assignments).")


def test_member_level_evidence() -> None:
    print("[5/10] Checking member-level evidence table...")
    mem_p = COMPACT_DIR / "exp4_47_member_evidence.tsv"
    df = pd.read_csv(mem_p, sep="\t")

    assert len(df) == 428
    assert int(df["verified_complete_interior_member"].sum()) == 326
    assert int(df["is_exact_model_match"].sum()) == 32
    assert int(df["is_verified_interior_with_exact_model"].sum()) == 27

    # Reject identity-by-length: exact match must require exact sequence hash status
    for _, row in df.iterrows():
        if row["is_exact_model_match"]:
            assert row["best_model_verification_status"] == "exact_sequence_hash_match"
        if row["is_verified_interior_with_exact_model"]:
            assert bool(row["verified_complete_interior_member"]) is True
            assert bool(row["is_exact_model_match"]) is True

    # Coordinates
    coord_joined = df[df["coordinate_provenance"].str.startswith("verified_loci_input_join")]
    assert len(coord_joined) == 310
    print("  PASS: 428 members verified (326 verified interior, 32 exact models, 27 verified interior with exact model).")


def test_flank_level_evidence() -> None:
    print("[6/10] Checking flank-level evidence table...")
    flank_p = COMPACT_DIR / "exp4_47_flank_evidence.tsv"
    df = pd.read_csv(flank_p, sep="\t")

    assert len(df) == 2883
    assert "annotated_guild_support" in df.columns
    assert (df["machinery_experimental_verification_status"] == "unverified_computational_annotation_confidence_unknown").all()
    assert (df["annotation_confidence_or_evalue"] == "unknown_not_in_shard_bundle").all()

    coord_joined = df[df["coordinate_provenance"].str.startswith("verified_loci_input_join")]
    assert len(coord_joined) == 346
    missing_coords = df[df["coordinate_provenance"] == "unknown_not_in_verified_loci_input"]
    assert len(missing_coords) == 2883 - 346
    assert (missing_coords["intergenic_distance_bp"] == "unknown_not_in_verified_loci_input").all()
    print("  PASS: 2,883 flanking loci verified with honest annotation provenance and confidence markers.")


def test_lead_coordinate_audit() -> None:
    print("[7/10] Checking Lead 1 and Lead 2 coordinate evidence...")
    lead_coord_p = COMPACT_DIR / "exp4_lead_flank_coordinates_independent.tsv"
    flank_p = COMPACT_DIR / "exp4_47_flank_evidence.tsv"
    ref_df = pd.read_csv(lead_coord_p, sep="\t")
    flank_df = pd.read_csv(flank_p, sep="\t")

    assert len(ref_df) == 24
    l1_ref = ref_df[ref_df["family_id"] == "22PRr1::22PRr1_2531955_3"]
    l2_ref = ref_df[ref_df["family_id"] == "V60PRr2::V60PRr2_1314979_4"]
    assert len(l1_ref) == 13
    assert len(l2_ref) == 11

    # Verify all 24 pairs are immediately adjacent (offset +/-1), cooriented, tight gaps
    assert (ref_df["same_strand"] == True).all()
    assert (ref_df["gene_offset"].abs() == 1).all()
    assert l1_ref["gap_bp_negative_overlap"].min() >= -4 and l1_ref["gap_bp_negative_overlap"].max() <= 37
    assert l2_ref["gap_bp_negative_overlap"].min() >= -4 and l2_ref["gap_bp_negative_overlap"].max() <= 6

    # Verify exact match in flank evidence table
    l1_flank = flank_df[(flank_df["family_id"] == "22PRr1::22PRr1_2531955_3") & (flank_df["flanking_specific_guilds"].str.contains("surface_assembly_pilin_archaellin"))]
    l2_flank = flank_df[(flank_df["family_id"] == "V60PRr2::V60PRr2_1314979_4") & (flank_df["flanking_specific_guilds"].str.contains("cell_coat_or_envelope"))]
    assert len(l1_flank) == 13
    assert len(l2_flank) == 11
    assert (l1_flank["same_strand"].astype(str) == "True").all()
    assert (l2_flank["same_strand"].astype(str) == "True").all()
    assert (l1_flank["relative_gene_offset"].abs() == 1).all()
    assert (l2_flank["relative_gene_offset"].abs() == 1).all()
    print("  PASS: All 24 lead context pairs verified immediately adjacent, cooriented, gaps [-4, 37] bp.")


def test_lead_dispositions_and_boundaries() -> None:
    print("[8/10] Checking bounded scientific dispositions for Lead 1 and Lead 2...")
    dossiers_p = COMPACT_DIR / "exp4_47_candidate_dossiers.tsv"
    df = pd.read_csv(dossiers_p, sep="\t")

    # Lead 1
    r1 = df[df["family_id"] == "22PRr1::22PRr1_2531955_3"].iloc[0]
    assert bool(r1["primary_significant_at_0_10"]) is True
    assert bool(r1["genus_significant_at_0_10"]) is False
    assert bool(r1["seed_is_verified_interior"]) is False  # contig-edge
    assert bool(r1["rep_is_verified_interior"]) is True
    assert r1["rep_model_verification_status"] == "no_model_available"
    assert int(r1["n_observed_flanking_loci"]) == 13
    assert int(r1["n_context_genus_blocks"]) == 6
    disp1 = str(r1["scientific_disposition_summary"])
    assert "CRITICAL FALSIFICATION" not in disp1
    assert "absence of structural homology" not in disp1
    assert "completely lacks" not in disp1
    assert "experimental AlphaFold" not in disp1
    assert "genus-restricted sensitivity does not" in disp1
    assert "F37Dr2" in disp1
    assert "S20PRr1" in disp1
    assert "PilJ" in disp1
    assert "Co-transcription" in disp1
    assert "archaeal pilin identity" in disp1
    assert r1["supported_evidence_route"] == "uncertain_insufficient_evidence"

    # Lead 2
    r2 = df[df["family_id"] == "V60PRr2::V60PRr2_1314979_4"].iloc[0]
    assert bool(r2["primary_significant_at_0_10"]) is True
    assert bool(r2["genus_significant_at_0_10"]) is False
    assert bool(r2["seed_is_verified_interior"]) is False  # contig-edge
    assert bool(r2["rep_is_verified_interior"]) is True
    assert r2["rep_model_verification_status"] == "no_model_available"
    assert int(r2["n_observed_flanking_loci"]) == 11
    assert int(r2["n_context_genus_blocks"]) == 5
    disp2 = str(r2["scientific_disposition_summary"])
    assert "CRITICAL FALSIFICATION" not in disp2
    assert "absence of structural homology" not in disp2
    assert "completely lacks" not in disp2
    assert "experimental AlphaFold" not in disp2
    assert "genus-restricted sensitivity does not" in disp2
    assert "must not be interpreted as absence of homology" in disp2
    assert "Co-transcription" in disp2
    assert "S-layer membership" in disp2
    assert r2["supported_evidence_route"] == "uncertain_insufficient_evidence"
    print("  PASS: Lead 1 and Lead 2 scientific dispositions adhere strictly to evidence boundaries.")


def test_independent_replay_and_msa_identity() -> None:
    print("[9/10] Checking independent replay and MSA byte-level sequence identity...")
    replay_json_p = COMPACT_DIR / "exp4_independent_polymer_and_msa_replay_audit.json"
    with open(replay_json_p) as fh:
        rep_data = json.load(fh)

    assert rep_data["all_17_shards_present"] is True
    assert rep_data["n_members_replayed"] == 428
    assert rep_data["n_verified_interior_members"] == 326
    assert rep_data["n_members_exact_sequence_hash_match"] == 32
    assert rep_data["n_verified_interior_members_exact_match"] == 27
    assert rep_data["member_evidence_agreement_pass"] is True
    assert rep_data["all_msa_taxa_valid"] is True
    assert rep_data["all_ungapped_sequences_byte_match"] is True

    # Validate Lead 1 and Lead 2 treefiles exist and are non-empty
    t1 = COMPACT_DIR / "22PRr1__22PRr1_2531955_3.iqtree.treefile"
    t2 = COMPACT_DIR / "V60PRr2__V60PRr2_1314979_4.iqtree.treefile"
    assert t1.exists() and t1.stat().st_size > 0
    assert t2.exists() and t2.stat().st_size > 0
    print("  PASS: Independent replay confirms 100% polymer and MSA sequence identity.")


def test_slurm_execution_receipts_and_owner_validation() -> None:
    print("[10/10] Checking Slurm execution receipts and scheduler owner validation...")
    receipts_p = PKG_ROOT / "manifests" / "slurm_execution_receipts.tsv"
    raw_records_p = PKG_ROOT / "manifests" / "raw_slurm_execution_records.tsv"

    assert receipts_p.exists(), f"Execution receipts table missing: {receipts_p}"
    assert raw_records_p.exists(), f"Raw scheduler records table missing: {raw_records_p}"

    receipts_df = pd.read_csv(receipts_p, sep="\t")
    raw_df = pd.read_csv(raw_records_p, sep="\t")

    assert len(receipts_df) >= 9, f"Expected at least 9 receipt rows, got {len(receipts_df)}"
    assert len(raw_df) == 42, f"Expected exactly 42 raw scheduler records for hohndor protein jobs, got {len(raw_df)}"

    # 1. Owner validation: ALL records must belong to expected owner 'hohndor'
    expected_owner = "hohndor"
    assert (receipts_df["user"] == expected_owner).all(), (
        f"Unauthorized users found in receipts: {receipts_df['user'].unique()}"
    )
    assert (raw_df["user"] == expected_owner).all(), (
        f"Unauthorized users found in raw scheduler rows: {raw_df['user'].unique()}"
    )

    # 2. Allowed protein job names validation
    allowed_job_names = {"eq_exp4_gspa", "eq_exp4_pilot", "eq_exp4_arr", "eq_exp4_col", "eq_exp4_faudit"}
    assert set(receipts_df["job_name"]).issubset(allowed_job_names), (
        f"Unexpected job names in receipts: {set(receipts_df['job_name']) - allowed_job_names}"
    )
    assert set(raw_df["job_name"]).issubset(allowed_job_names), (
        f"Unexpected job names in raw scheduler records: {set(raw_df['job_name']) - allowed_job_names}"
    )

    # 3. Exact required v2 array tasks 0..16 and collector ID
    raw_job_ids = set(raw_df["job_id"].astype(str))
    expected_v2_tasks = {f"15890701_{i}" for i in range(17)}
    assert expected_v2_tasks.issubset(raw_job_ids), (
        f"Missing required v2 array task IDs: {expected_v2_tasks - raw_job_ids}"
    )
    assert "15890718" in raw_job_ids, "Missing required v2 collector job ID 15890718"

    # All v2 tasks and collector must be COMPLETED with exit_code 0:0
    v2_records = raw_df[raw_df["job_id"].astype(str).isin(expected_v2_tasks | {"15890718"})]
    assert (v2_records["state"] == "COMPLETED").all(), (
        f"Non-completed v2 tasks: {v2_records[v2_records['state'] != 'COMPLETED']}"
    )
    assert (v2_records["exit_code"] == "0:0").all(), (
        f"Non-zero exit codes in v2 tasks: {v2_records[v2_records['exit_code'] != '0:0']}"
    )

    # 4. Proper primary-stage correspondence
    required_stage_mappings = {
        "stage2_gspa": ("15872818", "eq_exp4_gspa"),
        "stage3_pilot_v1": ("15873659", "eq_exp4_pilot"),
        "stage3_production_v1": ("15873660", "eq_exp4_arr"),
        "stage4_collect_v1": ("15873661", "eq_exp4_col"),
        "stage3_production_v2": ("15890701", "eq_exp4_arr"),
        "stage4_collect_v2": ("15890718", "eq_exp4_col"),
        "stage5_targeted_followup_lead_phylo": ("15892410", "eq_exp4_faudit"),
        "stage5_independent_replay": ("15892822", "eq_exp4_faudit"),
        "stage5_coordinate_audit": ("15892964", "eq_exp4_faudit"),
    }

    stages_present = set(receipts_df["pipeline_stage"])
    for stage, (exp_jid, exp_jname) in required_stage_mappings.items():
        assert stage in stages_present, f"Missing required pipeline stage in receipts: {stage}"
        row = receipts_df[receipts_df["pipeline_stage"] == stage].iloc[0]
        assert str(row["job_id"]) == exp_jid, f"Stage {stage} expected job_id {exp_jid}, got {row['job_id']}"
        assert str(row["job_name"]) == exp_jname, f"Stage {stage} expected job_name {exp_jname}, got {row['job_name']}"
        assert str(row["user"]) == expected_owner, f"Stage {stage} expected user {expected_owner}, got {row['user']}"
        assert str(row["state"]) == "COMPLETED", f"Stage {stage} not COMPLETED: {row['state']}"
        assert str(row["exit_code"]) == "0:0", f"Stage {stage} exit code non-zero: {row['exit_code']}"

    # 5. GSPA calibration limit honest reporting check
    gspa_row = receipts_df[receipts_df["pipeline_stage"] == "stage2_gspa"].iloc[0]
    gspa_ctx = str(gspa_row["submission_context"]).lower()
    assert "uncalibrated" in gspa_ctx or "exploratory" in gspa_ctx, (
        "GSPA receipt must explicitly document uncalibrated/exploratory status"
    )
    assert "blocked" in gspa_ctx, "GSPA receipt must document heldout calibration blocker"

    # 6. Targeted followup node and elapsed verification
    phylo_row = receipts_df[receipts_df["pipeline_stage"] == "stage5_targeted_followup_lead_phylo"].iloc[0]
    assert phylo_row["nodes"] == "nid00082", f"Targeted audit phylo expected nid00082, got {phylo_row['nodes']}"
    assert phylo_row["elapsed_hms"] == "00:00:14", f"Targeted audit phylo expected 00:00:14, got {phylo_row['elapsed_hms']}"

    # 7. Remote log manifest integrity check (lightweight owned small logs inventory)
    log_manifest_p = PKG_ROOT / "manifests" / "remote_slurm_log_manifest.tsv"
    assert log_manifest_p.exists(), f"Remote log manifest missing: {log_manifest_p}"
    log_df = pd.read_csv(log_manifest_p, sep="\t")
    assert len(log_df) >= 84, f"Expected at least 84 log records, got {len(log_df)}"
    assert (log_df["remote_path"].str.startswith("/scratch/project/k10018/eq-mechanism-20261109/exp4_surface/logs/")).all()
    assert (log_df["sha256"].str.len() == 64).all()
    assert (log_df["observed_at_utc"].str.len() > 10).all()
    assert set(log_df["status"]).issubset({"production_execution_clean", "failed_historical_iteration"})

    print(
        f"  PASS: All {len(receipts_df)} receipts, {len(raw_df)} raw scheduler rows, and {len(log_df)} owned log records "
        f"verified for owner '{expected_owner}', allowed job names, exact v2 tasks 0..16, and required stage correspondence."
    )


def main() -> None:
    print("================================================================================")
    print("Experiment 4 (`exp4_surface`): Frozen-Output Reproducibility Verification")
    print(f"Package Root: {PKG_ROOT}")
    print("================================================================================")
    test_sha256_checksums()
    test_global_summary_and_fdr()
    test_cohort_architecture_comparison()
    test_candidate_dossiers()
    test_member_level_evidence()
    test_flank_level_evidence()
    test_lead_coordinate_audit()
    test_lead_dispositions_and_boundaries()
    test_independent_replay_and_msa_identity()
    test_slurm_execution_receipts_and_owner_validation()
    print("================================================================================")
    print("ALL 10 VERIFICATION GATES PASSED CLEANLY (Lightweight Replay Successful).")
    print("================================================================================")


if __name__ == "__main__":
    main()
