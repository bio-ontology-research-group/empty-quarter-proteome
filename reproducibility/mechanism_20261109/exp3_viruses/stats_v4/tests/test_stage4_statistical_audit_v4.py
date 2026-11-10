"""Regression tests for the scoped Stage 4 v4 statistical repair."""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

BASE = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(BASE / "code"))
from stage4_statistical_audit_v4 import (  # noqa: E402
    DESIGN_LEAD_FAMILIES,
    aggregate_selected_sites,
    bh_adjust,
    build_family_assay_maps,
    conditional_family_tests,
    compute_abundance_metrics,
    diagnose_binary_separation,
    find_matched_alternative_controls,
    join_candidate_taxonomy,
    load_root_unresolved_targets,
    load_provisional_candidate,
    partition_detection_opportunity,
    resolve_graded_candidate_pairs,
    select_primary_environment_design,
    separation_lp_diagnostic,
    stratified_permutation_indices,
    validate_assay_provenance,
)

PROVISIONAL = ("43PRr3__self_SemiBin_469", "44Dr2::44Dr2_4002037", "6293eefbf49a4cf8cb58cb6aa2c6dd877e3e81450917bd018483580bf3028bb1")
UNRESOLVED = "41Sr3::41Sr3_56779"


def sha(label: str) -> str:
    return hashlib.sha256(label.encode()).hexdigest()


def pair(host: str, target: str, sequence: str, tier: str = "Tier1_verified_direct", host_genes: int = 0) -> dict:
    return {
        "host_mag_raw": host, "canonical_591_host_mag": host, "host_reconciliation_status": "canonical_591_direct",
        "target_id": target, "target_seq_sha256": sha(sequence), "votu_id": target,
        "evidence_tier": tier, "checkv_host_genes": host_genes, "checkv_viral_genes": 2,
        "is_primary_sequence_representative": 1,
    }


def taxonomy(genomes: list[str]) -> pd.DataFrame:
    return pd.DataFrame({
        "genome": genomes,
        "canonical_phylum": ["p__Thermoplasmatota"] * len(genomes),
        "canonical_class": ["c__x"] * len(genomes),
        "canonical_order": ["o__JACQPN01"] * len(genomes),
        "canonical_family": ["f__x"] * len(genomes),
        "canonical_genus": ["g__"] * len(genomes),
        "canonical_species": ["s__"] * len(genomes),
        "canonical_completeness": [91.62] * len(genomes),
        "canonical_contamination": [4.85] * len(genomes),
    })


def test_schema_and_exact_hash_are_strict() -> None:
    rows = pd.DataFrame([pair("h1", "t1", "a")])
    with pytest.raises(ValueError, match="missing required columns"):
        resolve_graded_candidate_pairs(rows.drop(columns=["votu_id"]), set())
    bad = rows.copy()
    bad.loc[0, "target_seq_sha256"] = "not-a-sha"
    with pytest.raises(ValueError, match="invalid target_seq_sha256"):
        resolve_graded_candidate_pairs(bad, set())
    missing = rows.copy()
    missing.loc[0, "target_seq_sha256"] = pd.NA
    with pytest.raises(ValueError, match="missing/blank target_seq_sha256"):
        resolve_graded_candidate_pairs(missing, set())


def test_root_disposition_must_be_individual_unresolved(tmp_path: Path) -> None:
    path = tmp_path / "root.json"
    path.write_text(json.dumps({"target": "t1", "root_review_disposition": "unresolved_viral_identity"}))
    assert load_root_unresolved_targets(path) == {"t1"}
    path.write_text(json.dumps({"target": "t1", "root_review_disposition": "accepted"}))
    with pytest.raises(ValueError, match="unresolved"):
        load_root_unresolved_targets(path)


def test_boundary_sensitive_provisional_root_witness(tmp_path: Path) -> None:
    payload = {
        "host": "h3", "target": "t3", "target_sha256": sha("c"),
        "disposition": "boundary_sensitive_provisional_not_primary_host_link",
        "canonical_spacer_sequence": "A" * 32, "exact_sequence_witness_verified": True,
        "spacer_target_interval_1based_inclusive": [37670, 37701],
        "catalog_viral_interval": [1, 34828], "fresh_checkv_viral_interval": [1, 43323],
    }
    path = tmp_path / "provisional.json"
    path.write_text(json.dumps(payload))
    key, loaded = load_provisional_candidate(path)
    assert key == ("h3", "t3", sha("c"))
    assert loaded["disposition"].endswith("not_primary_host_link")
    payload["catalog_viral_interval"] = [1, 40000]
    path.write_text(json.dumps(payload))
    with pytest.raises(ValueError, match="outside catalog core"):
        load_provisional_candidate(path)


def test_no_global_checkv_threshold_dedup_and_provisional_derivation() -> None:
    rows = pd.DataFrame([
        pair("h1", "high_cellular_but_not_root_excluded", "a", host_genes=500),
        pair("h2", "unresolved_target", "b", host_genes=70),
        pair("h3", "pending_target", "c"),
        pair("h4", "alias_a", "same"),
        pair("h4", "alias_b", "same"),
        pair("h5", "other_tier", "d", tier="Tier2_corroborated_single_spacer_or_provirus"),
    ])
    ledger, primary, provisional = resolve_graded_candidate_pairs(rows, {"unresolved_target"}, {("h3", "pending_target", sha("c"))})
    assert len(ledger) == len(rows)
    high = ledger.loc[ledger["target_id"].eq("high_cellular_but_not_root_excluded")].iloc[0]
    assert high["is_primary_detected_link"] == 1
    assert ledger["is_alias_duplicate"].sum() == 1
    assert ledger["is_unresolved_root_disposition"].sum() == 1
    assert set(primary["canonical_591_host_mag"]) == {"h1", "h4"}
    assert set(provisional["canonical_591_host_mag"]) == {"h1", "h3", "h4"}
    unresolved = ledger.loc[ledger["target_id"].eq("unresolved_target")].iloc[0]
    assert "negative_control" in unresolved["review_disposition_reason"]
    assert ledger.loc[ledger["target_id"].eq("other_tier"), "is_primary_detected_link"].eq(0).all()


def make_partition() -> pd.DataFrame:
    return pd.DataFrame({
        "genome": ["positive", "negative", "zero", "unassayed_negative"],
        "n_valid_spacers": [2, 3, 0, 4],
        "is_screen_positive_eligible": [1, 1, 0, 1],
        "primary_DETECTEDLINK": [1.0, 0.0, np.nan, 0.0],
        "provisional_DETECTEDLINK": [1.0, 0.0, np.nan, 0.0],
        "canonical_phylum": ["p__x"] * 4,
        "canonical_completeness": [90.0] * 4,
        "canonical_order": ["o__x"] * 4,
        "canonical_family": ["f__x"] * 4,
        "canonical_genus": ["g__"] * 4,
        "detection_opportunity_category": ["screen_positive_conditional_detection_cohort", "screen_positive_conditional_detection_cohort", "zero_spacers_untested_not_negative", "screen_positive_conditional_detection_cohort"],
    })


def test_zero_spacer_untested_and_family_absence_denominator() -> None:
    partition = make_partition()
    members = pd.DataFrame({"family_id": ["full_family", "dossier_family"], "genome": ["positive", "positive"]})
    full = pd.DataFrame({"family_id": ["full_family", "assay_universe", "assay_universe", "assay_universe"],
                         "genome": ["positive", "negative", "zero", "other"],
                         "final_promotion_gate_pass": [True] * 4})
    result, coverage = conditional_family_tests(partition, members, full)
    full_row = result.loc[result["family_id"].eq("full_family") & result["outcome_set"].eq("primary")].iloc[0]
    assert full_row["n_assay_covered_eligible"] == 2  # positive + negative; zero and other are excluded
    assert full_row["present_positive"] == 1 and full_row["absent_negative"] == 1
    dossier = result.loc[result["family_id"].eq("dossier_family") & result["outcome_set"].eq("primary")].iloc[0]
    assert pd.isna(dossier["fisher_exact_pvalue"]) and pd.isna(dossier["bh_fdr_qvalue_fixed47"])
    assert dossier["test_status"] == "non_estimable_no_verified_absence_denominator"
    assert result["fixed_BH_family_size"].eq(2).all()
    assert coverage["assay_coverage_status"].str.contains("dossier_only").sum() == 1


def test_fixed_full_bh_family_counts_nonestimable_hypotheses() -> None:
    q = bh_adjust(np.array([0.001, 0.02, np.nan]), family_size=47)
    assert np.isclose(q[0], 0.047)
    assert np.isclose(q[1], 0.47)
    assert np.isnan(q[2])
    with pytest.raises(ValueError, match="family_size"):
        bh_adjust(np.array([0.1, 0.2]), family_size=1)


def test_zero_variance_is_not_separation_and_lp_witness_is_real() -> None:
    y = np.array([0, 0, 1, 1])
    constant = diagnose_binary_separation(y, np.zeros(4))
    assert constant["status"] == "constant_predictor_no_contrast_not_separation"
    assert not constant["quasi_complete"]
    all_zero_y = diagnose_binary_separation(np.zeros(4), np.array([0, 1, 0, 1]))
    assert all_zero_y["status"] == "constant_outcome_not_predictor_separation"
    complete = diagnose_binary_separation(y, y)
    assert complete["status"] == "complete_separation"
    X = np.column_stack([np.ones(4), y])
    lp = separation_lp_diagnostic(y, X)
    assert lp["complete"] and lp["minimum_signed_margin"] >= 1
    overlap = separation_lp_diagnostic(y, np.column_stack([np.ones(4), np.array([0, 1, 1, 0])]))
    assert not overlap["complete"] and not overlap["quasi_complete"]


def control_inputs() -> tuple[pd.DataFrame, pd.DataFrame]:
    partition = pd.DataFrame({
        "genome": ["candidate", "unknown_related", "same_genus", "missing_prevalence"],
        "canonical_phylum": ["p__x"] * 4,
        "canonical_order": ["o__1", "o__1", "o__1", "o__1"],
        "canonical_family": ["f__1"] * 4,
        "canonical_genus": ["g__known", "g__", "g__known", "g__other"],
        "canonical_completeness": [90.0, 91.0, 91.0, 91.0],
        "detection_opportunity_category": ["screen_positive_conditional_detection_cohort"] * 4,
    })
    # The same-known-genus alternative must be excluded; blank genus remains unknown relatedness.
    metrics = pd.DataFrame({
        "genome": partition["genome"], "sample_prevalence": [10, 11, 11, np.nan],
        "mean_sample_abundance": [0.1] * 4, "abundance_assay_status": ["observed_matrix_row"] * 4,
    })
    return partition, metrics


def test_unknown_genus_is_not_different_genus_and_no_imputation() -> None:
    partition, metrics = control_inputs()
    matched = find_matched_alternative_controls(["candidate"], {"candidate"}, {"candidate"}, set(), partition, metrics)
    row = matched.iloc[0]
    assert row["matched_alternative_mag"] == "unknown_related"
    assert row["relatedness_status"] == "unknown_genus_relatedness_not_established"
    assert "not_unrelated_taxon" in row["matching_status"]
    missing = partition.copy()
    missing.loc[missing["genome"].eq("candidate"), "canonical_completeness"] = np.nan
    unmatched = find_matched_alternative_controls(["candidate"], {"candidate"}, {"candidate"}, set(), missing, metrics)
    assert pd.isna(unmatched.iloc[0]["matched_alternative_mag"])
    assert "NA_candidate_missing" in unmatched.iloc[0]["matching_status"]


def test_prevalence_missing_is_not_zero(tmp_path: Path) -> None:
    # Matrix row absence is NA, never a zero-prevalence pseudo-observation.
    mag = pd.DataFrame([[0.0, 1.0]], index=["observed"], columns=["s1", "s2"])
    metrics = compute_abundance_metrics(mag, ["observed", "missing"], ["s1", "s2"])
    missing = metrics.loc[metrics["genome"].eq("missing")].iloc[0]
    assert pd.isna(missing["sample_prevalence"])
    assert missing["abundance_assay_status"] == "missing_matrix_row_not_zero"


def metadata_and_fields() -> tuple[pd.DataFrame, pd.DataFrame]:
    metadata = pd.DataFrame({
        "sample_id": ["s1", "s2", "s3", "s4", "s5"], "site": ["1", "1", "2", "3", "4"],
        "compartment": ["Deep", "Deep", "Deep", "Surface", "Surface"], "replicate": [1, 2, 1, 1, 1],
        "trip": ["trip1"] * 5, "lat": [0.0] * 5, "lon": [0.0] * 5,
    })
    fields = pd.DataFrame({
        "sample_id": metadata["sample_id"], "dune_field": ["A", "A", "A", "B", "B"],
        "site_mag": ["1", "1", "2", "3", "4"], "trip_mag": ["trip1"] * 5,
        "compartment_mag": metadata["compartment"],
    })
    return metadata, fields


def test_metadata_only_selection_is_order_and_outcome_invariant() -> None:
    metadata, fields = metadata_and_fields()
    groups1, selected1, rule1 = select_primary_environment_design(metadata, fields, set(metadata["sample_id"]))
    # Arbitrary outcome permutations cannot enter this function; row-order shuffle must not matter.
    groups2, selected2, rule2 = select_primary_environment_design(metadata.sample(frac=1, random_state=4), fields.sample(frac=1, random_state=5), set(metadata["sample_id"]))
    pd.testing.assert_frame_equal(groups1.sort_values(["trip", "compartment"]).reset_index(drop=True), groups2.sort_values(["trip", "compartment"]).reset_index(drop=True))
    assert rule1["selected_sample_ids"] == rule2["selected_sample_ids"]
    assert rule1["selected_compartment"] == "Deep" and not rule1["outcome_values_used_for_selection"]


def test_current_rund_compartment_synonym_is_audited_not_guessed() -> None:
    metadata, fields = metadata_and_fields()
    metadata["compartment"] = "PlantRoot"
    fields["compartment_mag"] = "RootAdjacent"
    _, _, selection = select_primary_environment_design(metadata, fields, set(metadata["sample_id"]))
    assert "RootAdjacent" in selection["compartment_crosswalk"]
    fields["compartment_mag"] = "Deep"
    with pytest.raises(ValueError, match="disagree"):
        select_primary_environment_design(metadata, fields, set(metadata["sample_id"]))


def test_field_permutations_never_cross_strata_and_singletons_fixed() -> None:
    fields = np.array(["A", "A", "A", "B", "B", "singleton"])
    perms = stratified_permutation_indices(fields, n_perms=200, seed=17)
    assert perms.shape == (200, 6)
    for row in perms:
        assert np.array_equal(np.sort(row[:3]), np.arange(3))
        assert np.array_equal(np.sort(row[3:5]), np.arange(3, 5))
        assert row[5] == 5
    assert (np.unique(perms[:, 5]) == [5]).all()


def test_site_aggregation_uses_same_trip_compartment_samples() -> None:
    metadata, _ = metadata_and_fields()
    selected = metadata.loc[metadata["compartment"].eq("Deep")].copy()
    selected["dune_field"] = ["A", "A", "A"]
    agg = aggregate_selected_sites(selected, np.array([1.0, 3.0, 5.0]), np.array([2.0, 4.0, 8.0]))
    assert len(agg) == 2
    assert agg.loc[agg["site"].eq("1"), "host"].iloc[0] == 2.0
    assert agg.loc[agg["site"].eq("1"), "virus"].iloc[0] == 3.0


def test_actual_tier1_counts_provisional_status_and_exact_taxonomy_join() -> None:
    pairs_path = BASE / "results/stage3_archaeal_virus_host_pairs_all_tiers_v2.tsv"
    tax_path = BASE / "inputs/taxonomy_reconciliation_all5438.tsv.gz"
    if not pairs_path.exists() or not tax_path.exists():
        pytest.skip("compact actual inputs unavailable")
    pairs = pd.read_csv(pairs_path, sep="\t")
    ledger, primary, provisional = resolve_graded_candidate_pairs(pairs, {UNRESOLVED}, {PROVISIONAL})
    assert len(ledger) == len(pairs)
    assert ledger["evidence_tier"].eq("Tier1_verified_direct").sum() == 8
    assert len(primary) == 5 and primary["canonical_591_host_mag"].nunique() == 4
    assert len(provisional) == 6 and provisional["canonical_591_host_mag"].nunique() == 5
    new = ledger.loc[ledger["target_id"].eq(PROVISIONAL[1])].iloc[0]
    assert new["is_boundary_sensitive_provisional"] == 1 and new["is_primary_detected_link"] == 0
    assert "outside_catalog_core" in new["review_disposition_reason"]
    tax = pd.read_csv(tax_path, sep="\t")
    joined = join_candidate_taxonomy(provisional, tax)
    new_tax = joined.loc[joined["canonical_591_host_mag"].eq(PROVISIONAL[0])].iloc[0]
    assert new_tax["canonical_phylum"] == "p__Thermoplasmatota"
    assert new_tax["canonical_order"] == "o__JACQPN01"
    assert new_tax["canonical_genus"] == "g__"
    assert new_tax["canonical_completeness"] == 91.62 and new_tax["canonical_contamination"] == 4.85
    old = joined.loc[joined["canonical_591_host_mag"].eq("41PRr2__single_SemiBin_236")].iloc[0]
    assert old["canonical_order"] == "o__JACPAO01"


def test_actual_detection_and_frozen_family_coverage() -> None:
    required = [BASE / "results/stage1_mag_screening_status.tsv", BASE / "results/stage3_all591_archaeal_mag_host_evidence_v2.tsv",
                BASE / "results/stage3_archaeal_virus_host_pairs_all_tiers_v2.tsv", BASE / "inputs/exp4_47_member_evidence.tsv",
                BASE / "inputs/task_a_stage1b_all_candidate_member_loci_verified.tsv.gz", BASE / "inputs/taxonomy_reconciliation_all5438.tsv.gz"]
    if not all(p.exists() for p in required):
        pytest.skip("compact actual inputs unavailable")
    pairs = pd.read_csv(required[2], sep="\t")
    _, primary, provisional = resolve_graded_candidate_pairs(pairs, {UNRESOLVED}, {PROVISIONAL})
    partition = partition_detection_opportunity(pd.read_csv(required[0], sep="\t"), pd.read_csv(required[1], sep="\t"),
                                                 pd.read_csv(required[5], sep="\t"), primary, provisional)
    assert len(partition) == 591
    assert partition["n_valid_spacers"].eq(0).sum() == 447
    assert partition["is_screen_positive_eligible"].sum() == 144
    assert partition.loc[partition["n_valid_spacers"].eq(0), "primary_DETECTEDLINK"].isna().all()
    members, full = pd.read_csv(required[3], sep="\t"), pd.read_csv(required[4], sep="\t")
    families, _, universes, coverage = build_family_assay_maps(members, full)
    assert len(families) == 47
    assert sum(x is not None for x in universes.values()) == 44
    assert full["genome"].nunique() == 237
    result, _ = conditional_family_tests(partition, members, full)
    primary_result = result.loc[result["outcome_set"].eq("primary")]
    assert len(primary_result) == 47 and primary_result["fisher_exact_pvalue"].notna().sum() == 44
    assert primary_result["bh_fdr_qvalue_fixed47"].notna().sum() == 44
    assert primary_result.loc[primary_result["assay_coverage_status"].str.contains("complete_assignment_exact_family"), "n_assay_covered_eligible"].eq(58).all()
    assert primary_result["n_eligible_detected_links_total"].eq(4).all()
    assert set(DESIGN_LEAD_FAMILIES) <= set(primary_result["family_id"])


def test_actual_complete_assignment_receipts_and_failed_gate_retention() -> None:
    full_path = BASE / "inputs/task_a_stage1b_all_candidate_member_loci_verified.tsv.gz"
    stage1 = BASE.parent.parent / "science-20261108/canonical-results/stage1_summary.json"
    stage1b = BASE.parent.parent / "science-20261108/canonical-results/stage1b_summary.json"
    if not full_path.exists() or not stage1.exists() or not stage1b.exists():
        pytest.skip("full-locus table or source completion receipts unavailable")
    full = pd.read_csv(full_path, sep="\t")
    provenance = validate_assay_provenance(full, stage1, stage1b)
    assert provenance["n_screened_genomes"] == 237
    assert provenance["n_original_proteins_in_stage1_scope"] == 777177
    assert provenance["n_candidate_member_locus_rows_complete_membership"] == 41494
    assert provenance["n_exact_family_ids"] == 8776
    assert provenance["n_final_promotion_gate_false_rows_retained_as_membership_evidence"] == 3184
    failed_members = full.loc[~full["final_promotion_gate_pass"].astype(bool), ["family_id", "genome"]]
    assert len(failed_members) == 3184
    members = pd.read_csv(BASE / "inputs/exp4_47_member_evidence.tsv", sep="\t")
    _, member_maps, _, _ = build_family_assay_maps(members, full, provenance)
    frozen_failed = failed_members.loc[failed_members["family_id"].isin(member_maps)]
    assert len(frozen_failed) > 0
    # Presence uses complete original assignment; no failed-promotion row is dropped.
    assert all(row.genome in member_maps[row.family_id] for row in frozen_failed.itertuples(index=False))
