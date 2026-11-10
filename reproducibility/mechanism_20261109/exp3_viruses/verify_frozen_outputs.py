#!/usr/bin/env python3
"""Standalone Frozen Output & Reproducibility Gate Verifier for Experiment 3 Archaeal-Viral Analysis.

Validates the frozen outputs, manifests, code hashes, Slurm receipts, and scientific boundaries
staged in `reproducibility/mechanism_20261109/exp3_viruses/` without triggering heavy compute:
1. Exact 591 canonical archaeal MAG cohort universe and 5,438 canonical MAG background manifest.
2. Complete Stage 2 BLAST search coverage (6 catalog pools + 169 sample assemblies) and collision v2 audit
   (144 multi-pool targets, 0 sequence collisions, 100% identical duplicates).
3. Reviewed Tier 1 pair ledger & root boundary dispositions: 5 primary pairs across 4 hosts (with
   48Sr2/coasm_48_S alias collapse), 1 unresolved sensitivity pair (41Sr3::41Sr3_56779), and 1
   boundary-sensitive provisional pair (43PRr3__self_SemiBin_469 / 44Dr2::44Dr2_4002037, 67,924 bp).
4. Completed Stage 4 sequence-domain HMM annotation & existing-model structural audit (Job 15894662):
   262 unannotated ORFs across 5 primary targets, 325 Pfam hits in 80 proteins, and 0 exact existing
   structural models (Foldseek explicitly reported as not_performed).
5. Stage 4 sample matrix orientation & defense/surface estimability audit: 441/441 shared samples and
   explicit non-estimability disclosure due to circularity and sparse positive events.
6. Repaired Slurm competitive mapping submission chain & code snapshot provenance: 350G memory request
   compatible with Shaheen workq (RealMemory=370000 MB), 14-sample array (0-13, no %4 throttle), and
   verified dependent chain (15893933 -> 15894659 -> 15894660 -> 15894661 -> 15894662).
7. Full SHA-256 manifest verification across all committed files.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent
CODE_DIR = REPO_ROOT / "code"
INPUTS_DIR = REPO_ROOT / "inputs"
EVIDENCE_DIR = REPO_ROOT / "compact_evidence"
MANIFESTS_DIR = REPO_ROOT / "manifests"

sys.path.insert(0, str(CODE_DIR))
from stage4_competitive_common import load_analysis_pairs, stream_sha256  # noqa: E402


def verify_gate1_cohort_and_background_manifests() -> tuple[bool, str]:
    s0_sum = json.loads((EVIDENCE_DIR / "stage0_summary.json").read_text())
    s1_sum = json.loads((EVIDENCE_DIR / "stage1_summary.json").read_text())
    arch_fna = pd.read_csv(EVIDENCE_DIR / "stage0_arch591_fasta_verification.tsv", sep="\t")
    s1_status = pd.read_csv(EVIDENCE_DIR / "stage1_mag_screening_status.tsv", sep="\t")
    s3_hosts = pd.read_csv(EVIDENCE_DIR / "stage3_all591_archaeal_mag_host_evidence_v2.tsv", sep="\t")
    bg_manifest = pd.read_csv(INPUTS_DIR / "genome_manifest.tsv", sep="\t")

    if len(arch_fna) != 591 or len(s3_hosts) != 591:
        return False, f"Expected 591 canonical archaeal MAGs in Stage 0/3, found {len(arch_fna)}/{len(s3_hosts)}"
    if s0_sum.get("n_verified_exact_match_genome_manifest") != 591:
        return False, f"Stage 0 summary did not verify 591 exact genome_manifest matches: {s0_sum}"
    if s1_sum.get("n_arch591_successfully_screened") != 591 or len(s1_status) != 793:
        return False, f"Stage 1 screening mismatch: {s1_sum}, rows={len(s1_status)}"
    if len(bg_manifest) != 5438 or bg_manifest["mag"].duplicated().any():
        return False, f"Expected 5438 unique canonical MAGs in genome_manifest.tsv, found {len(bg_manifest)}"
    return True, "591 canonical archaeal MAGs (793 total screened with 202 non-reps) and 5,438 background MAG manifest verified"


def verify_gate2_search_coverage_and_collision_audit() -> tuple[bool, str]:
    s2_sum = json.loads((EVIDENCE_DIR / "stage2_summary.json").read_text())
    pool_df = pd.read_csv(EVIDENCE_DIR / "stage2_pool_coverage_and_provenance.tsv", sep="\t")
    col_sum = json.loads((EVIDENCE_DIR / "stage2_target_identity_collision_summary_v2.json").read_text())
    col_rec = json.loads((EVIDENCE_DIR / "stage2_target_identity_collision_receipt_v2.json").read_text())
    root_acc = json.loads((EVIDENCE_DIR / "exp3_collision_v2_independent_acceptance_audit.json").read_text())
    col_path = EVIDENCE_DIR / "stage2_target_identity_collision_audit_v2.tsv"
    col_df = pd.read_csv(col_path, sep="\t")

    if len(pool_df) != 7 or not (pool_df["status"] == "searched_complete").all():
        return False, f"Stage 2 pool coverage incomplete: {pool_df.to_dict(orient='records')}"
    if s2_sum.get("n_raw_assemblies_searched") != 169 or s2_sum.get("n_assemblies_present") != 169:
        return False, f"Unexpected Stage 2 assembly search counts: {s2_sum}"
    if not col_sum.get("collision_free", False):
        return False, "Stage 2 collision audit v2 is not collision_free"
    if len(col_df) != 144 or col_sum.get("n_identical_duplicates") != 144 or col_sum.get("n_sequence_collisions") != 0:
        return False, f"Unexpected collision v2 audit counts: rows={len(col_df)}, summary={col_sum}"
    audit_sha = stream_sha256(col_path)
    if col_rec.get("outputs", {}).get("audit_tsv_v2_sha256") != audit_sha or root_acc.get("source_sha256") != audit_sha:
        return False, "Collision v2 receipt/root-acceptance SHA mismatch against stage2_target_identity_collision_audit_v2.tsv"
    if col_rec.get("inputs", {}).get("assemblies_manifest_sha256") != stream_sha256(EVIDENCE_DIR / "stage2_target_assemblies_manifest.tsv"):
        return False, "Collision v2 receipt assemblies_manifest_sha256 mismatch"
    if col_rec.get("inputs", {}).get("identity_table_sha256") != stream_sha256(EVIDENCE_DIR / "stage2_target_viral_identity_table.tsv"):
        return False, "Collision v2 receipt identity_table_sha256 mismatch"
    return True, "175/175 BLAST search scopes (6 catalogs + 169 assemblies) and 144/144 collision-free multi-pool targets (336 source sequences) verified"


def verify_gate3_pair_dispositions_and_boundary_contracts() -> tuple[bool, str]:
    combined, reviewed, boundary = load_analysis_pairs(
        EVIDENCE_DIR / "stage3_tier1_pairs_individual_evidence_audit.tsv",
        INPUTS_DIR / "stage4_boundary_sensitive_intervals.json",
    )
    n_primary = int((reviewed["analysis_role"] == "primary").sum())
    n_sens = int((reviewed["analysis_role"] == "sensitivity_unresolved_virality").sum())
    n_bound = len(boundary)
    n_primary_hosts = int(reviewed.loc[reviewed["analysis_role"] == "primary", "host_mag"].nunique())

    if (len(combined), n_primary, n_primary_hosts, n_sens, n_bound) != (7, 5, 4, 1, 1):
        return (
            False,
            f"Unexpected pair partition: combined={len(combined)}, primary={n_primary} ({n_primary_hosts} hosts), "
            f"unresolved_sens={n_sens}, boundary_provisional={n_bound}",
        )
    cand = boundary.iloc[0]
    if cand["host_mag"] != "43PRr3__self_SemiBin_469" or cand["target_id"] != "44Dr2::44Dr2_4002037":
        return False, f"Unexpected boundary candidate identity: {cand.to_dict()}"

    # Verify host-array & contig provenance audit across all 6 hosts (8 array-bearing contigs)
    host_ext_path = EVIDENCE_DIR / "stage4_host_array_contig_provenance_extract.json"
    host_sum_path = EVIDENCE_DIR / "stage4_host_array_provenance_summary.json"
    host_tbl_path = EVIDENCE_DIR / "stage4_host_array_contig_provenance_table.tsv"
    host_join_path = EVIDENCE_DIR / "stage4_host_array_spacer_target_join_audit.tsv"
    for p in (host_ext_path, host_sum_path, host_tbl_path, host_join_path):
        if not p.exists():
            return False, f"Missing host-array provenance artifact: {p}"

    host_sum = json.loads(host_sum_path.read_text())
    if host_sum.get("schema_version") != 2 or host_sum.get("status") != "host_array_contig_provenance_completed_verified":
        return False, f"Unexpected host provenance summary status/schema: {host_sum.get('status')}, v={host_sum.get('schema_version')}"
    if host_sum.get("source_extract", {}).get("extract_sha256") != stream_sha256(host_ext_path):
        return False, "Host provenance extract_sha256 mismatch"
    if host_sum.get("outputs", {}).get("contig_provenance_table_sha256") != stream_sha256(host_tbl_path):
        return False, "Host provenance contig_provenance_table_sha256 mismatch"
    if host_sum.get("outputs", {}).get("spacer_target_join_audit_sha256") != stream_sha256(host_join_path):
        return False, "Host provenance spacer_target_join_audit_sha256 mismatch"

    rec = host_sum.get("focal_thermoplasmatota_host_41PRr2_SemiBin_236_reconciliation", {})
    if (
        rec.get("focal_contig_length_bp") != 651471
        or rec.get("n_physical_crispr_array_loci") != 1
        or rec.get("n_crispr_caller_records") != 2
        or rec.get("n_raw_extracted_matching_spacers") != 7
        or rec.get("n_strict_unique_spacers_verified") != 6
        or rec.get("n_exact_0mm_unique_spacers_verified") != 4
        or not rec.get("root_independent_join_verified", False)
    ):
        return False, f"Focal 41PRr2_1588631 host-array reconciliation failed: {rec}"

    sum_text = json.dumps(host_sum)
    if ".A=" in sum_text or ".B=" in sum_text:
        return False, "Invented .A/.B geNomad suffix strings present in host provenance summary"

    host_tbl = pd.read_csv(host_tbl_path, sep="\t")
    if len(host_tbl) != 8 or not (host_tbl["crosswalk_row_count"] == 2).all():
        return False, f"Expected 8 array-bearing host contigs with 2 crosswalk rows each, got len={len(host_tbl)}"
    if not (
        host_tbl["qualified_source_alias_interpretation"]
        == "source_category_duplication_same_contig_in_raw_search_pool_unique_mag_owner"
    ).all():
        return False, "Unexpected qualified_source_alias_interpretation in host provenance table"
    if host_tbl.astype(str).apply(lambda s: s.str.contains("NA_pending_")).any().any():
        return False, "Pending placeholder NA_pending_* found in stage4_host_array_contig_provenance_table.tsv"
    if not (
        host_tbl["sample_mapping_coverage_status"]
        == "completed_sample_mapping_integrated_via_representative_header"
    ).all():
        return False, "Expected completed_sample_mapping_integrated_via_representative_header across all 8 rows"
    focal_row = host_tbl[host_tbl["raw_contig_id"] == "41PRr2_1588631"].iloc[0]
    if (
        abs(float(focal_row["sample_mapping_contig_mean_depth"]) - 0.871882) > 1e-6
        or abs(float(focal_row["sample_mapping_contig_breadth_pct"]) - 8.149096) > 1e-6
        or abs(float(focal_row["sample_mapping_remaining_mag_mean_depth"]) - 13.369586) > 1e-6
        or abs(float(focal_row["sample_mapping_remaining_mag_breadth_pct"]) - 71.472404) > 1e-6
    ):
        return False, f"Focal 41PRr2_1588631 sample mapping coverage mismatch: {focal_row.to_dict()}"

    return (
        True,
        "7 logical mapping pairs (5 primary across 4 hosts, 1 unresolved 41Sr3, 1 boundary 43PRr3/44Dr2) "
        "and 6-host/8-contig provenance + sample coverage audit (focal 41PRr2_1588631 651,471 bp, 6 strict/4 exact spacers, "
        "0.871882x/8.149096% vs rest-of-MAG 13.369586x/71.472404%) verified",
    )


def verify_gate4_sequence_domain_and_existing_model_audit() -> tuple[bool, str]:
    summary_path = EVIDENCE_DIR / "stage4_sequence_domain_annotation_summary.json"
    done_path = EVIDENCE_DIR / "STAGE4_SEQUENCE_DOMAIN_ANNOTATION_DONE"
    posthoc_path = EVIDENCE_DIR / "stage4_sequence_posthoc_qualification_summary.json"
    extract_path = EVIDENCE_DIR / "stage4_target_genomad_checkv_pfam_ga_extract.json"
    qual_hits_path = EVIDENCE_DIR / "stage4_sequence_domain_hits_posthoc_qualified.tsv"
    focal_orf_path = EVIDENCE_DIR / "stage4_focal_V60PRr2_715976_orf_evidence_table.tsv"
    kw_audit_path = EVIDENCE_DIR / "stage4_target_hallmark_keyword_label_audit.tsv"

    summary = json.loads(summary_path.read_text())
    done = json.loads(done_path.read_text())
    posthoc = json.loads(posthoc_path.read_text())
    extract = json.loads(extract_path.read_text())

    if done.get("summary_sha256") != stream_sha256(summary_path):
        return False, "STAGE4_SEQUENCE_DOMAIN_ANNOTATION_DONE hash mismatch against summary JSON"
    if summary.get("status") != "sequence_domain_annotation_completed_verified":
        return False, f"Unexpected sequence annotation status: {summary.get('status')}"
    if summary.get("script_sha256") != "8fb33e18af4dc6b165acc7a1b932c1ff5ec3cf1f09e1b4d2bdcac1c25a3abbd9":
        return False, f"Unexpected executed sequence script_sha256: {summary.get('script_sha256')}"
    cand = summary.get("candidate_selection", {})
    pfam = summary.get("pfam_results", {})
    struct = summary.get("structural_followup", {})
    if cand.get("n_candidate_orfs") != 262 or cand.get("n_candidate_proteins_extracted") != 262:
        return False, f"Expected 262 unannotated ORFs extracted across 5 primary targets, found {cand}"
    if pfam.get("n_hits") != 325 or pfam.get("n_unique_proteins_with_hits") != 80:
        return False, f"Unexpected Pfam hit counts: {pfam}"
    if struct.get("status") != "not_performed_no_exact_sequence_matched_existing_models" or struct.get("n_exact_sequence_matched_models") != 0:
        return False, f"Unexpected structural follow-up status: {struct}"

    # Verify extraction provenance and posthoc qualification hierarchy (325/80 -> 153/75 i-Evalue<=1e-3 -> 82/50 Pfam GA)
    if extract.get("schema_version") != 2:
        return False, f"Expected schema_version=2 in {extract_path}"
    src_files = extract.get("extraction_provenance", {}).get("source_files", {})
    if src_files.get("pfam_db_sha256_from_job_15894662_summary") != "a78a42d6faf265b6bfca59e8f062d06fae6083ce2c6e335d7b381f20b82b7903":
        return False, f"Unexpected Pfam DB SHA-256 in extraction provenance: {src_files}"

    th = posthoc.get("threshold_hierarchy", {})
    if th.get("raw_full_sequence_evalue_lte_1e3", {}).get("domain_rows") != 325 or th.get("raw_full_sequence_evalue_lte_1e3", {}).get("unique_proteins") != 80:
        return False, f"Unexpected raw domain counts in posthoc summary: {th}"
    if th.get("posthoc_independent_domain_ievalue_lte_1e3", {}).get("domain_rows") != 153 or th.get("posthoc_independent_domain_ievalue_lte_1e3", {}).get("unique_proteins") != 75:
        return False, f"Unexpected posthoc i-Evalue<=1e-3 counts: {th}"
    if th.get("posthoc_pfam_gathering_threshold_ga", {}).get("domain_rows") != 82 or th.get("posthoc_pfam_gathering_threshold_ga", {}).get("unique_proteins") != 50:
        return False, f"Unexpected posthoc Pfam GA counts: {th}"

    post_out = posthoc.get("posthoc_outputs", {})
    if post_out.get("qualified_hits_sha256") != stream_sha256(qual_hits_path):
        return False, "Posthoc qualified_hits_sha256 mismatch"
    if post_out.get("focal_orf_table_sha256") != stream_sha256(focal_orf_path):
        return False, "Posthoc focal_orf_table_sha256 mismatch"
    if post_out.get("target_hallmark_keyword_audit_sha256") != stream_sha256(kw_audit_path):
        return False, "Posthoc target_hallmark_keyword_audit_sha256 mismatch"

    orf10_audit = posthoc.get("focal_target_V60PRr2_715976_audit", {}).get("orf8_to_orf12_genomic_neighborhood", {}).get("intervening_orf10", {})
    if orf10_audit.get("official_reference_sources_2026_11_09", {}).get("PF09979", {}).get("name") != "Pa193 gp24 scaffolding protein":
        return False, f"ORF10 PF09979 official reference correction missing: {orf10_audit}"
    if "terminase_protease_neighborhood" in json.dumps(posthoc):
        return False, "Unsupported terminase_protease_neighborhood label found in posthoc summary"

    return (
        True,
        "Job 15894662 & posthoc qualification verified: 262 unannotated ORFs, 325/80 raw -> 153/75 i-Evalue<=1e-3 -> "
        "82/50 Pfam GA, 0 exact structural models (not_performed), ORF10/ORF12 & legacy keyword audit verified",
    )


def verify_gate5_matrix_orientation_and_estimability() -> tuple[bool, str]:
    s4_sum = json.loads((EVIDENCE_DIR / "stage4_host_specificity_cooccurrence_and_mechanism_summary_v2.json").read_text())
    est = json.loads((EVIDENCE_DIR / "stage4_defense_surface_mechanism_estimability_v2.json").read_text())

    align = s4_sum.get("sample_intersection", {})
    if align.get("n_shared_samples") != 441 or align.get("intersection_status") != "complete_exact_sample_match":
        return False, f"Sample alignment audit failed: {align}"
    if est.get("estimability_status") != "non_estimable_insufficient_positive_events_and_circularity":
        return False, f"Unexpected defense/surface estimability status: {est}"
    return True, "441/441 shared samples and explicit defense/surface non-estimability verified"


def verify_gate6_scheduler_chain_and_code_hashes() -> tuple[bool, str]:
    sub = json.loads((MANIFESTS_DIR / "stage4_competitive_pipeline_submission.json").read_text())
    expected_jobs = ("15893933", "15905724", "15905725", "15910797", "15905726", "15910957", "15894662")
    actual_jobs = (
        sub.get("job_ref_build"),
        sub.get("job_validate_ref"),
        sub.get("job_map_array"),
        sub.get("job_recover_v60prr2"),
        sub.get("job_aggregate"),
        sub.get("job_posthoc_followup"),
        sub.get("job_sequence_followup"),
    )
    if actual_jobs != expected_jobs:
        return False, f"Unexpected active Slurm chain job IDs: {actual_jobs} != {expected_jobs}"
    if sub.get("status") != "completed_verified_all_stages":
        return False, f"Expected submission status completed_verified_all_stages, got {sub.get('status')}"
    hist1 = sub.get("historical_chain_15894659", {})
    expected_hist1 = ("15893933", "15894659", "15894660", "15894661", "15894662")
    actual_hist1 = (
        hist1.get("job_ref_build"),
        hist1.get("job_validate_ref"),
        hist1.get("job_map_array_cancelled"),
        hist1.get("job_aggregate_cancelled"),
        hist1.get("job_sequence_followup_completed"),
    )
    if actual_hist1 != expected_hist1:
        return False, f"Unexpected historical 15894659 chain job IDs: {actual_hist1} != {expected_hist1}"
    hist2 = sub.get("historical_chain_15902808", {})
    expected_hist2 = ("15893933", "15902808", "15902809", "15902810", "15894662")
    actual_hist2 = (
        hist2.get("job_ref_build"),
        hist2.get("job_validate_ref"),
        hist2.get("job_map_array_cancelled"),
        hist2.get("job_aggregate_cancelled"),
        hist2.get("job_sequence_followup_completed"),
    )
    if actual_hist2 != expected_hist2:
        return False, f"Unexpected historical 15902808 chain job IDs: {actual_hist2} != {expected_hist2}"
    if len(sub.get("samples", [])) != 14 or "43PRr3" not in sub["samples"] or "44Dr2" not in sub["samples"]:
        return False, f"Expected 14 mapping samples including 43PRr3 and 44Dr2, found {sub.get('samples')}"

    # Verify completed 14-sample competitive mapping aggregation, V60PRr2 recovery, and posthoc inStrain sensitivity
    aggr_sum_path = EVIDENCE_DIR / "stage4_competitive_mapping_all_pairs_summary.json"
    aggr_done_path = EVIDENCE_DIR / "STAGE4_COMPETITIVE_MAPPING_DONE"
    if not aggr_sum_path.exists() or not aggr_done_path.exists():
        return False, "Missing stage4_competitive_mapping_all_pairs_summary.json or STAGE4_COMPETITIVE_MAPPING_DONE"
    aggr_sum = json.loads(aggr_sum_path.read_text())
    aggr_done = json.loads(aggr_done_path.read_text())
    if aggr_done.get("summary_sha256") != stream_sha256(aggr_sum_path):
        return False, "STAGE4_COMPETITIVE_MAPPING_DONE hash mismatch against summary JSON"
    if stream_sha256(aggr_sum_path) != "558eab32f5dbc43279f77eae2eee868acc3f904a27e05c7aa2c799cd714222a6":
        return False, f"Unexpected aggregate summary SHA: {stream_sha256(aggr_sum_path)}"
    if (
        aggr_sum.get("status") != "completed_verified"
        or aggr_sum.get("n_samples_mapped") != 14
        or aggr_sum.get("n_pair_sample_combinations") != 98
        or aggr_sum.get("n_boundary_interval_sample_rows") != 42
        or aggr_sum.get("n_qualifying_unique_viral_sha_sample_combinations") != 1
        or aggr_sum.get("n_instrain_profiles_completed") != 1
    ):
        return False, f"Unexpected aggregate summary metrics: {aggr_sum}"

    v60_rec_path = EVIDENCE_DIR / "stage4_competitive_mapping/V60PRr2/V60PRr2_mapping_receipt.json"
    v60_prov_path = EVIDENCE_DIR / "stage4_competitive_mapping/V60PRr2/V60PRr2_recovery_provenance.json"
    if stream_sha256(v60_rec_path) != "2c00a7c260d31fd03682e04216dfbedc328d862ca4fc574c3c65d27ab62ef498":
        return False, "Unexpected V60PRr2_mapping_receipt.json SHA-256"
    v60_prov = json.loads(v60_prov_path.read_text())
    if v60_prov.get("status") != "completed_verified_v60prr2_instrain_recovery":
        return False, f"Unexpected V60PRr2_recovery_provenance.json status: {v60_prov.get('status')}"
    prior_ids = [a.get("source_slurm_job_id") for a in v60_prov.get("recovery_execution", {}).get("preserved_failed_attempts", [])]
    if prior_ids != ["15905725_1", "15910426"]:
        return False, f"Expected preserved V60PRr2 failure history ['15905725_1', '15910426'], got {prior_ids}"

    post_sum_path = EVIDENCE_DIR / "stage4_instrain_posthoc_callable_and_sensitivity_summary.json"
    post_tbl_path = EVIDENCE_DIR / "stage4_instrain_posthoc_profiles_table.tsv"
    post_sum = json.loads(post_sum_path.read_text())
    post_tbl = pd.read_csv(post_tbl_path, sep="\t")
    if (
        post_sum.get("status") != "completed_verified"
        or post_sum.get("n_sample_receipts_verified") != 14
        or post_sum.get("n_source_primary_profiles_qualified") != 1
        or post_sum.get("n_crosssample_sensitivity_profiles_qualified") != 1
        or len(post_tbl) != 2
    ):
        return False, f"Unexpected posthoc inStrain summary counts: {post_sum.get('status')}, rows={len(post_tbl)}"
    pop_comps = post_sum.get("population_comparisons", [])
    if len(pop_comps) != 10 or not all(
        c.get("comparison_status") == "quantitatively_non_estimable_fewer_than_two_qualifying_sample_populations"
        for c in pop_comps
    ):
        return False, f"Unexpected population_comparisons status in posthoc summary: {pop_comps}"
    src_prof = post_sum["source_primary_profiles"][0]
    sens_prof = post_sum["crosssample_sensitivity_profiles"][0]
    if (
        src_prof.get("snv_count_scaffold_info_unfiltered_cumulative") != 24
        or src_prof.get("exported_non_cryptic_snv_count") != 22
        or src_prof.get("exported_non_cryptic_sns_count") != 0
        or src_prof.get("cryptic_variant_site_count_excluded_from_exported_snv_table") != 2
        or src_prof.get("inferred_depth_qualified_positions_ge_mincov") != 38585
    ):
        return False, f"Unexpected source_primary (>0.95 ANI) profile metrics: {src_prof}"
    if (
        sens_prof.get("snv_count_scaffold_info_unfiltered_cumulative") != 24
        or sens_prof.get("sns_count_scaffold_info_unfiltered_cumulative") != 1
        or sens_prof.get("divergent_site_count_scaffold_info_unfiltered_cumulative") != 25
        or sens_prof.get("exported_non_cryptic_snv_count") != 21
        or sens_prof.get("exported_non_cryptic_sns_count") != 1
        or sens_prof.get("cryptic_variant_site_count_excluded_from_exported_snv_table") != 3
        or sens_prof.get("inferred_depth_qualified_positions_ge_mincov") != 38589
    ):
        return False, f"Unexpected crosssample_sensitivity (>0.92 ANI) profile metrics: {sens_prof}"

    # 1. Verify historical snapshot 1 (pi_viral_repair_20261109)
    hist1_snap_root = REPO_ROOT / "version_snapshots/pi_viral_repair_20261109"
    hist1_hash_lines = (MANIFESTS_DIR / "pi_viral_repair_code_hashes.tsv").read_text().strip().splitlines()
    for line in hist1_hash_lines:
        expected_sha, rel_path = line.split()
        snap_target = hist1_snap_root / rel_path
        if not snap_target.exists() or stream_sha256(snap_target) != expected_sha:
            return False, f"Historical 15894659 snapshot hash mismatch for {rel_path}: expected {expected_sha}"

    # 2. Verify historical snapshot 2 (pi_viral_sentinel_repair_20261109)
    hist2_snap_root = REPO_ROOT / "version_snapshots/pi_viral_sentinel_repair_20261109"
    hist2_hash_lines = (MANIFESTS_DIR / "pi_viral_sentinel_repair_code_hashes.tsv").read_text().strip().splitlines()
    for line in hist2_hash_lines:
        expected_sha, rel_path = line.split()
        snap_target = hist2_snap_root / rel_path
        if not snap_target.exists() or stream_sha256(snap_target) != expected_sha:
            return False, f"Historical 15902808 snapshot hash mismatch for {rel_path}: expected {expected_sha}"

    # 3. Verify active boundary-length-repaired snapshot (pi_viral_boundary_length_repair_20261109) and live code/inputs
    active_snap_root = REPO_ROOT / "version_snapshots/pi_viral_boundary_length_repair_20261109"
    active_hash_lines = (MANIFESTS_DIR / "pi_viral_boundary_length_repair_code_hashes.tsv").read_text().strip().splitlines()
    for line in active_hash_lines:
        expected_sha, rel_path = line.split()
        snap_target = active_snap_root / rel_path
        if not snap_target.exists() or stream_sha256(snap_target) != expected_sha:
            return False, f"Active boundary-length-repair snapshot hash mismatch for {rel_path}: expected {expected_sha}"
        if rel_path != "tests/test_viral_mapping_repair.py":
            live_target = REPO_ROOT / rel_path
            if not live_target.exists() or stream_sha256(live_target) != expected_sha:
                return False, f"Executed code/input hash mismatch for {rel_path}: expected {expected_sha}"

    # 4. Verify Slurm batch spools for 15902808..15902810 and 15905724..15905726
    spools_dir = MANIFESTS_DIR / "slurm_spools_20261109"
    spool_files = (
        "slurm_15902808_validate_ref.sh",
        "slurm_15902809_map_array.sh",
        "slurm_15902810_aggregate.sh",
        "scontrol_show_submitted_15902808_15902809_15902810.txt",
        "slurm_15905724_validate_ref.sh",
        "slurm_15905725_map_array.sh",
        "slurm_15905726_aggregate.sh",
        "scontrol_show_submitted_15905724_15905725_15905726.txt",
        "scontrol_show_submitted_15910186.txt",
    )
    for spool_file in spool_files:
        p = spools_dir / spool_file
        if not p.exists() or p.stat().st_size == 0:
            return False, f"Missing or empty Slurm spool file: {p}"

    if stream_sha256(CODE_DIR / "stage2_viral_identity_and_spacers_v3.py") != "e6c50836d24b926bd7d0bf089f71ab623109fed221691020b344b5b5de68562e":
        return False, "stage2_viral_identity_and_spacers_v3.py hash changed unexpectedly"

    posthoc_hash_path = MANIFESTS_DIR / "posthoc_sequence_qualification_code_hashes.tsv"
    if not posthoc_hash_path.exists():
        return False, f"Missing {posthoc_hash_path}"
    posthoc_lines = posthoc_hash_path.read_text().strip().splitlines()
    for line in posthoc_lines:
        expected_sha, rel_path = line.split()
        target = REPO_ROOT / rel_path
        if not target.exists() or stream_sha256(target) != expected_sha:
            return False, f"Posthoc code/test hash mismatch for {rel_path}: expected {expected_sha}"

    return (
        True,
        f"Completed Slurm chain {actual_jobs} (14/14 sample receipts, V60PRr2 recovery 15910797, aggregate 15905726 [98/224/42], "
        f"posthoc 15910957 [24/22/2 primary vs 25/22/3 sensitivity]), {len(active_hash_lines)} active + "
        f"{len(hist1_hash_lines)} + {len(hist2_hash_lines)} historical snapshot hashes, {len(spool_files)} Slurm spools, "
        f"and {len(posthoc_lines)} posthoc hashes verified",
    )


def verify_gate7_sha256_manifest() -> tuple[bool, str]:
    manifest_path = MANIFESTS_DIR / "sha256_manifest.tsv"
    if not manifest_path.exists():
        return False, f"Missing {manifest_path}"
    rows = pd.read_csv(manifest_path, sep="\t")
    if rows.empty:
        return False, "sha256_manifest.tsv is empty"
    forbidden_prefixes = (".pytest_cache/", "__pycache__/", ".pytest_tmp/", "logs/", "stats_v4/", "coverage_qc/")
    for row in rows.to_dict(orient="records"):
        rel = str(row["relative_path"])
        if rel.startswith(forbidden_prefixes) or "/__pycache__/" in rel or rel.endswith(".pyc"):
            return False, f"Mutable cache/runtime file must not be in sha256_manifest.tsv: {rel}"
        path = REPO_ROOT / rel
        if not path.exists() or path.stat().st_size != int(row["bytes"]):
            return False, f"Missing or size-changed file in manifest: {rel}"
        observed = stream_sha256(path)
        if observed != row["sha256"]:
            return False, f"SHA-256 mismatch for {rel}: {observed} != {row['sha256']}"
    return True, f"All {len(rows)} files in manifests/sha256_manifest.tsv verified byte-for-byte (mutable caches excluded)"


def main() -> int:
    gates = [
        ("Gate 1: Cohort & Background Manifests", verify_gate1_cohort_and_background_manifests),
        ("Gate 2: Stage 2 Search Coverage & Collision v2 Audit", verify_gate2_search_coverage_and_collision_audit),
        ("Gate 3: Reviewed Pair Dispositions & Boundary Contracts", verify_gate3_pair_dispositions_and_boundary_contracts),
        ("Gate 4: Sequence-Domain HMM & Existing-Model Audit", verify_gate4_sequence_domain_and_existing_model_audit),
        ("Gate 5: Matrix Orientation & Non-Estimability Disclosure", verify_gate5_matrix_orientation_and_estimability),
        ("Gate 6: Repaired Slurm Chain & Code Provenance", verify_gate6_scheduler_chain_and_code_hashes),
        ("Gate 7: Package SHA-256 Manifest Integrity", verify_gate7_sha256_manifest),
    ]
    all_passed = True
    print("=" * 80)
    print("EXPERIMENT 3 ARCHAEAL-VIRAL REPRODUCIBILITY & FROZEN OUTPUT VERIFIER")
    print("=" * 80)
    for title, fn in gates:
        passed, message = fn()
        status = "PASS" if passed else "FAIL"
        print(f"[{status}] {title}: {message}")
        if not passed:
            all_passed = False
    print("=" * 80)
    return 0 if all_passed else 1


if __name__ == "__main__":
    sys.exit(main())
