#!/usr/bin/env python3
"""Pytest suite for Cross-Paper Shared Foundation.

Verifies:
1. Frozen foundation data passes all 6 empirical verification gates cleanly.
2. Adversarial / mutation tests prove verifier catches:
   - wrong ID preserving count (manifest or taxonomy substitution)
   - duplicate ID (manifest, taxonomy, or inventory)
   - invalid Thermoplasmatota subset (dropped member or substituted non-Thermo archaeon)
   - full taxonomy string or rank column mismatch between thermo237 and all5438
   - invalid GTDB rank prefixes (e.g. not d__, p__, c__, o__, f__, g__, s__)
   - wrong rank counts (e.g. 6 ranks instead of 7)
   - inconsistency between individual canonical rank columns and classification string
   - unexpected cohort values (must strictly be 'EQ_desert' or 'external_comparator')
   - altered cohort membership (EQ <-> external flip, or ID overlap with EQ MAGs)
   - invalid external comparator source labels (must be 'GTDB_r232_rep' or 'MGnify')
   - corrupted SHA256 checksums in manifest.
"""
from __future__ import annotations

import copy
from pathlib import Path

import pytest

from verify_shared_foundation import (
    PKG_ROOT,
    INPUTS_DIR,
    MANIFEST_P,
    load_tsv_rows,
    verify_all_gates,
    verify_archaea_and_thermo_subsets_and_concordance,
    verify_canonical_taxonomy_syntax_and_ranks,
    verify_external_source_labels_and_harmonization_boundary,
    verify_manifest_and_taxonomy_ids,
    verify_metabolic_cohort_join,
    verify_sha256_manifest,
)


@pytest.fixture(scope="module")
def foundation_data():
    """Load all four frozen foundation tables once for the test module."""
    manifest_rows = load_tsv_rows(INPUTS_DIR / "genome_manifest.tsv")
    tax_rows = load_tsv_rows(INPUTS_DIR / "taxonomy_reconciliation_all5438.tsv.gz", is_gzip=True)
    thermo_rows = load_tsv_rows(INPUTS_DIR / "taxonomy_reconciliation_thermo237.tsv")
    exp1_rows = load_tsv_rows(INPUTS_DIR / "exp1_model_inventory_and_qc_484.tsv")
    return {
        "manifest_rows": manifest_rows,
        "tax_rows": tax_rows,
        "thermo_rows": thermo_rows,
        "exp1_rows": exp1_rows,
    }


def test_real_foundation_all_gates_pass():
    """Gates 1-6: Full verification passes cleanly on actual frozen files."""
    res = verify_all_gates(PKG_ROOT)
    assert res["status"] == "PASS"
    assert res["gates_passed"] == 6
    assert res["total_mags"] == 5438
    assert res["archaeal_mags"] == 591
    assert res["thermoplasmatota_mags"] == 237
    assert res["other_archaeal_mags"] == 354
    assert res["metabolic_inventory_total"] == 484
    assert res["metabolic_eq_cohort"] == 237
    assert res["metabolic_external_cohort"] == 247
    assert res["external_source_counts"]["GTDB_r232_rep"] == 218
    assert res["external_source_counts"]["MGnify"] == 29


def test_catches_wrong_id_preserving_count_in_manifest(foundation_data):
    """Verifier must catch substituted wrong ID preserving exact count 5,438."""
    manifest_mutated = copy.deepcopy(foundation_data["manifest_rows"])
    tax_rows = foundation_data["tax_rows"]
    manifest_mutated[0]["mag"] = "WRONG_SUBSTITUTED_MAG_9999"
    with pytest.raises(AssertionError, match="ID mismatch between manifest and taxonomy reconciliation"):
        verify_manifest_and_taxonomy_ids(manifest_mutated, tax_rows)


def test_catches_wrong_id_preserving_count_in_taxonomy(foundation_data):
    """Verifier must catch substituted wrong ID preserving exact count 5,438."""
    manifest_rows = foundation_data["manifest_rows"]
    tax_mutated = copy.deepcopy(foundation_data["tax_rows"])
    tax_mutated[0]["genome"] = "WRONG_SUBSTITUTED_MAG_9999"
    with pytest.raises(AssertionError, match="ID mismatch between manifest and taxonomy reconciliation"):
        verify_manifest_and_taxonomy_ids(manifest_rows, tax_mutated)


def test_catches_duplicate_id_in_manifest(foundation_data):
    """Verifier must catch duplicate MAG ID preserving total count 5,438."""
    manifest_mutated = copy.deepcopy(foundation_data["manifest_rows"])
    tax_rows = foundation_data["tax_rows"]
    manifest_mutated[1]["mag"] = manifest_mutated[0]["mag"]
    with pytest.raises(AssertionError, match="Duplicate MAG IDs detected in manifest"):
        verify_manifest_and_taxonomy_ids(manifest_mutated, tax_rows)


def test_catches_duplicate_id_in_taxonomy(foundation_data):
    """Verifier must catch duplicate genome ID in taxonomy table."""
    manifest_rows = foundation_data["manifest_rows"]
    tax_mutated = copy.deepcopy(foundation_data["tax_rows"])
    tax_mutated[1]["genome"] = tax_mutated[0]["genome"]
    with pytest.raises(AssertionError, match="Duplicate genome IDs detected in taxonomy"):
        verify_manifest_and_taxonomy_ids(manifest_rows, tax_mutated)


def test_catches_duplicate_id_in_metabolic_inventory(foundation_data):
    """Verifier must catch duplicate genome ID in metabolic 484 inventory."""
    exp1_mutated = copy.deepcopy(foundation_data["exp1_rows"])
    exp1_mutated[1]["genome"] = exp1_mutated[0]["genome"]
    thermo_ids = set(r["genome"] for r in foundation_data["thermo_rows"])
    all_ids = set(r["mag"] for r in foundation_data["manifest_rows"])
    with pytest.raises(AssertionError, match="Duplicate genome IDs found in exp1 inventory"):
        verify_metabolic_cohort_join(exp1_mutated, thermo_ids, all_ids)


def test_catches_invalid_thermo_subset_dropped_member(foundation_data):
    """Verifier must catch dropped Thermoplasmatota genome (count 236 != 237)."""
    tax_rows = foundation_data["tax_rows"]
    thermo_mutated = copy.deepcopy(foundation_data["thermo_rows"][:-1])
    with pytest.raises(AssertionError, match="Expected 237 rows in thermo237 table"):
        verify_archaea_and_thermo_subsets_and_concordance(tax_rows, thermo_mutated)


def test_catches_invalid_thermo_subset_substituted_non_thermo(foundation_data):
    """Verifier must catch substituting a non-Thermo archaeon into thermo237 while keeping count 237."""
    tax_rows = foundation_data["tax_rows"]
    thermo_mutated = copy.deepcopy(foundation_data["thermo_rows"])
    non_thermo = next(
        r["genome"] for r in tax_rows
        if r["canonical_domain"] == "d__Archaea" and r["canonical_phylum"] != "p__Thermoplasmatota"
    )
    thermo_mutated[0]["genome"] = non_thermo
    with pytest.raises(AssertionError, match="Mismatch between Thermoplasmatota in all5438 and thermo237 table"):
        verify_archaea_and_thermo_subsets_and_concordance(tax_rows, thermo_mutated)


def test_catches_full_taxonomy_mismatch_between_thermo237_and_all5438(foundation_data):
    """Verifier must catch mismatched full taxonomy string for a genome between thermo237 and all5438."""
    tax_rows = foundation_data["tax_rows"]
    thermo_mutated = copy.deepcopy(foundation_data["thermo_rows"])
    # Corrupt the genus in the classification string of the first genome in thermo237
    orig = thermo_mutated[0]["canonical_gtdb_r232_classification"]
    thermo_mutated[0]["canonical_gtdb_r232_classification"] = orig.replace(";g__", ";g__MUTATED_")
    with pytest.raises(AssertionError, match="Full taxonomy string mismatch for genome"):
        verify_archaea_and_thermo_subsets_and_concordance(tax_rows, thermo_mutated)


def test_catches_rank_column_mismatch_between_thermo237_and_all5438(foundation_data):
    """Verifier must catch mismatched rank column for a genome between thermo237 and all5438."""
    tax_rows = foundation_data["tax_rows"]
    thermo_mutated = copy.deepcopy(foundation_data["thermo_rows"])
    # Corrupt the canonical_genus column value in thermo237 while leaving classification string identical
    thermo_mutated[0]["canonical_genus"] = "g__CORRUPTED_GENUS"
    with pytest.raises(AssertionError, match="Rank column 'canonical_genus' mismatch for genome"):
        verify_archaea_and_thermo_subsets_and_concordance(tax_rows, thermo_mutated)


def test_catches_invalid_rank_prefix_in_canonical_taxonomy(foundation_data):
    """Verifier must catch invalid GTDB rank prefix (e.g. 'phylum__' instead of 'p__')."""
    tax_mutated = copy.deepcopy(foundation_data["tax_rows"])
    orig = tax_mutated[0]["canonical_gtdb_r232_classification"]
    # Replace 'p__' with 'phylum__'
    tax_mutated[0]["canonical_gtdb_r232_classification"] = orig.replace(";p__", ";phylum__")
    tax_mutated[0]["canonical_phylum"] = "phylum__" + tax_mutated[0]["canonical_phylum"][3:]
    with pytest.raises(AssertionError, match="phylum rank does not start with prefix 'p__'"):
        verify_canonical_taxonomy_syntax_and_ranks(tax_mutated)


def test_catches_invalid_rank_count_in_canonical_taxonomy(foundation_data):
    """Verifier must catch lineage missing ranks (e.g. 6 ranks instead of 7)."""
    tax_mutated = copy.deepcopy(foundation_data["tax_rows"])
    parts = tax_mutated[0]["canonical_gtdb_r232_classification"].split(";")
    tax_mutated[0]["canonical_gtdb_r232_classification"] = ";".join(parts[:6])  # drop species
    with pytest.raises(AssertionError, match="has 6 taxonomy ranks; expected exactly 7"):
        verify_canonical_taxonomy_syntax_and_ranks(tax_mutated)


def test_catches_rank_column_inconsistency_with_classification(foundation_data):
    """Verifier must catch inconsistency between individual rank column and classification string."""
    tax_mutated = copy.deepcopy(foundation_data["tax_rows"])
    tax_mutated[0]["canonical_genus"] = "g__INCONSISTENT_GENUS"
    with pytest.raises(AssertionError, match="rank column 'canonical_genus' .* inconsistent with lineage rank"):
        verify_canonical_taxonomy_syntax_and_ranks(tax_mutated)


def test_catches_empty_taxonomy_string(foundation_data):
    """Verifier must catch empty canonical GTDB r232 classification."""
    tax_mutated = copy.deepcopy(foundation_data["tax_rows"])
    tax_mutated[0]["canonical_gtdb_r232_classification"] = ""
    with pytest.raises(AssertionError, match="has empty canonical taxonomy string"):
        verify_canonical_taxonomy_syntax_and_ranks(tax_mutated)


def test_catches_unexpected_cohort_value_in_inventory(foundation_data):
    """Verifier must catch cohort value that is neither 'EQ_desert' nor 'external_comparator'."""
    exp1_mutated = copy.deepcopy(foundation_data["exp1_rows"])
    thermo_ids = set(r["genome"] for r in foundation_data["thermo_rows"])
    all_ids = set(r["mag"] for r in foundation_data["manifest_rows"])
    ext_idx = next(i for i, r in enumerate(exp1_mutated) if r["cohort"] == "external_comparator")
    exp1_mutated[ext_idx]["cohort"] = "external_unverified"
    with pytest.raises(AssertionError, match="has unexpected cohort value 'external_unverified'"):
        verify_metabolic_cohort_join(exp1_mutated, thermo_ids, all_ids)


def test_catches_altered_cohort_membership_external_flipped_to_eq(foundation_data):
    """Verifier must catch external genome reassigned to EQ_desert cohort."""
    exp1_mutated = copy.deepcopy(foundation_data["exp1_rows"])
    thermo_ids = set(r["genome"] for r in foundation_data["thermo_rows"])
    all_ids = set(r["mag"] for r in foundation_data["manifest_rows"])
    ext_idx = next(i for i, r in enumerate(exp1_mutated) if r["cohort"] == "external_comparator")
    exp1_mutated[ext_idx]["cohort"] = "EQ_desert"
    with pytest.raises(AssertionError, match="Expected exactly 237 EQ_desert rows"):
        verify_metabolic_cohort_join(exp1_mutated, thermo_ids, all_ids)


def test_catches_altered_cohort_membership_eq_flipped_to_external(foundation_data):
    """Verifier must catch EQ genome reassigned to external cohort."""
    exp1_mutated = copy.deepcopy(foundation_data["exp1_rows"])
    thermo_ids = set(r["genome"] for r in foundation_data["thermo_rows"])
    all_ids = set(r["mag"] for r in foundation_data["manifest_rows"])
    eq_idx = next(i for i, r in enumerate(exp1_mutated) if r["cohort"] == "EQ_desert")
    exp1_mutated[eq_idx]["cohort"] = "external_comparator"
    with pytest.raises(AssertionError, match="Expected exactly 237 EQ_desert rows"):
        verify_metabolic_cohort_join(exp1_mutated, thermo_ids, all_ids)


def test_catches_altered_cohort_membership_overlap(foundation_data):
    """Verifier must catch external genome ID colliding with an EQ MAG ID."""
    exp1_mutated = copy.deepcopy(foundation_data["exp1_rows"])
    thermo_ids = set(r["genome"] for r in foundation_data["thermo_rows"])
    all_ids = set(r["mag"] for r in foundation_data["manifest_rows"])
    ext_idx = next(i for i, r in enumerate(exp1_mutated) if r["cohort"] == "external_comparator")
    some_bacterial_mag = next(r["mag"] for r in foundation_data["manifest_rows"] if r["mag"] not in thermo_ids)
    exp1_mutated[ext_idx]["genome"] = some_bacterial_mag
    with pytest.raises(AssertionError, match="External comparator IDs overlap with Empty Quarter MAG catalogue"):
        verify_metabolic_cohort_join(exp1_mutated, thermo_ids, all_ids)


def test_catches_invalid_external_source_label(foundation_data):
    """Verifier must catch unexpected source label for an external comparator."""
    exp1_mutated = copy.deepcopy(foundation_data["exp1_rows"])
    ext_idx = next(i for i, r in enumerate(exp1_mutated) if r["cohort"] == "external_comparator")
    exp1_mutated[ext_idx]["source"] = "NCBI_RefSeq"
    with pytest.raises(AssertionError, match="Expected external source labels"):
        verify_external_source_labels_and_harmonization_boundary(exp1_mutated)


def test_catches_corrupted_manifest_sha(tmp_path):
    """Verifier must catch altered SHA256 checksum in manifest."""
    corrupted_manifest = tmp_path / "corrupted_manifest.tsv"
    content = MANIFEST_P.read_text().replace(
        "835ded5ac9aee7e50be1405ebf5e9794a4382157b07d30a157d593e94b004bb5",
        "0000000000000000000000000000000000000000000000000000000000000000",
    )
    corrupted_manifest.write_text(content)
    with pytest.raises(AssertionError, match="Checksum mismatch"):
        verify_sha256_manifest(corrupted_manifest, PKG_ROOT)
