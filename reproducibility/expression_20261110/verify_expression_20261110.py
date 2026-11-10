#!/usr/bin/env python3
"""Standalone verifier for reproducibility/expression_20261110 (Karaoz et al. 2026 external metaproteome reanalysis).

Verifies:
1. Stage 1 v2 homology & search-database eligibility audit counts (2,321 derep MAGs, 113 Archaea,
   38 Thermoplasmatota, 274,006 ORFs, 122,370 valid EMSL50964 proteins, 333 archaeal-linked proteins,
   78 Thermoplasmatota-linked proteins, 1,514 query families).
2. Stage 2 expression validation audit (30 .mzTab files, 0 discordant PSM Q rows, 0/15 split composite
   search FASTAs deposited -> full_background_verified=False).
3. Non-estimability of primary Lead 1 (22PRr1::22PRr1_2531955_3 -> nonestimable_search_database_inclusion_unverified),
   Lead 2 (V60PRr2::V60PRr2_1314979_4 -> nonestimable_no_homolog_in_study_genomes), and primary family/clade warming
   (nonestimable_full_search_background_unavailable).
4. Conservative vs dataset-reported observed-subset Thermoplasmatota protein detection breakdown
   (2 confirmed >=2 peptides, 15 single peptide, 4 shared peptide, 34 rejected spectral count without Q<=0.01
   unambiguous peptide, 23 undetected).
5. Byte-exact SHA-256 integrity of all packaged files against manifests/sha256_manifest.tsv.
"""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
import sys
import unittest

PACKAGE_DIR = Path(__file__).resolve().parent


def verify_expression_package(quiet: bool = False) -> dict:
    s1_path = PACKAGE_DIR / "results" / "stage1_v2_homology_eligibility" / "stage1_v2_homology_eligibility_summary.json"
    s2_path = PACKAGE_DIR / "results" / "stage2_expression_validation" / "stage2_expression_validation_summary.json"
    bg_path = PACKAGE_DIR / "results" / "stage2_expression_validation" / "karaoz2026_search_database_background_completeness_audit.json"
    manifest_path = PACKAGE_DIR / "manifests" / "sha256_manifest.tsv"

    for p in (s1_path, s2_path, bg_path, manifest_path):
        if not p.is_file():
            raise FileNotFoundError(f"Missing required expression_20261110 file: {p}")

    s1 = json.loads(s1_path.read_text(encoding="utf-8"))
    s2 = json.loads(s2_path.read_text(encoding="utf-8"))
    bg = json.loads(bg_path.read_text(encoding="utf-8"))

    assert s1["counts"]["total_derep_mags"] == 2321
    assert s1["counts"]["derep_archaeal_mags"] == 113
    assert s1["counts"]["derep_thermoplasmatota_mags"] == 38
    assert s1["counts"]["predicted_orfs_in_113_archaeal_mags_v2"] == 274006
    assert s1["counts"]["emsl50964_total_valid_detected_proteins"] == 122370
    assert s1["counts"]["emsl50964_proteins_linked_to_archaeal_mags"] == 333
    assert s1["counts"]["emsl50964_proteins_linked_to_thermoplasmatota_mags"] == 78
    assert s1["counts"]["total_verified_query_families"] == 1514

    assert s2["runtime_verified_audit_counts"]["n_mztab_files_audited"] == 30
    assert s2["runtime_verified_audit_counts"]["total_discordant_psm_q_rows"] == 0
    assert s2["runtime_verified_audit_counts"]["full_15split_composite_search_fasta_deposited"] is False
    assert bg["total_deposited_files_in_massive"] == 604
    assert s2["primary_clade_warming_estimability_status"] == "nonestimable_full_search_background_unavailable"

    roles_cons = s2["eq_family_conservative_detection_summary_by_cohort_role"]
    assert roles_cons["primary_lead_1_pilj_archaellum_context"] == {
        "nonestimable_search_database_inclusion_unverified": 1
    }
    assert roles_cons["primary_lead_2_slayer_coat_context"] == {
        "nonestimable_no_homolog_in_study_genomes": 1
    }

    thermo_obs = s2["thermoplasmatota_proteins_in_emsl50964_dataset_reported_observed_subset_summary"]
    assert thermo_obs == {
        "undetected_in_soluble_metaproteome_not_proof_of_absence": 23,
        "rejected_spectral_count_without_qualifying_q01_unambiguous_peptide": 34,
        "provisional_single_unambiguous_peptide_q01": 15,
        "confirmed_detected_ge2_unambiguous_peptides_q01": 2,
        "provisional_clade_or_family_single_shared_peptide": 4,
    }

    with manifest_path.open("r", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert len(rows) >= 20
    for row in rows:
        target = PACKAGE_DIR / row["relative_path"]
        assert target.is_file(), f"Missing file in sha256_manifest.tsv: {row['relative_path']}"
        assert target.stat().st_size == int(row["size_bytes"]), f"Size mismatch: {row['relative_path']}"
        digest = hashlib.sha256(target.read_bytes()).hexdigest()
        assert digest == row["sha256"], f"SHA-256 mismatch: {row['relative_path']}"

    # Run the unit test suite programmatically
    loader = unittest.TestLoader()
    suite = loader.discover(str(PACKAGE_DIR / "tests"), pattern="test_*.py")
    runner = unittest.TextTestRunner(stream=sys.stderr if not quiet else open("/dev/null", "w"), verbosity=1)
    result = runner.run(suite)
    if not result.wasSuccessful():
        raise RuntimeError(f"Unit test suite failed: {len(result.failures)} failures, {len(result.errors)} errors")

    summary = {
        "module": "expression_20261110",
        "assertion_id": "expression_karaoz2026",
        "status": "PASS",
        "files_verified": len(rows),
        "unit_tests_run": result.testsRun,
        "derep_archaeal_mags": 113,
        "derep_thermoplasmatota_mags": 38,
        "mztab_files_audited": 30,
        "full_background_verified": False,
        "lead1_status": "nonestimable_search_database_inclusion_unverified",
        "lead2_status": "nonestimable_no_homolog_in_study_genomes",
        "primary_warming_status": "nonestimable_full_search_background_unavailable",
    }
    if not quiet:
        print(json.dumps(summary, indent=2))
    return summary


if __name__ == "__main__":
    verify_expression_package(quiet="--quiet" in sys.argv)
