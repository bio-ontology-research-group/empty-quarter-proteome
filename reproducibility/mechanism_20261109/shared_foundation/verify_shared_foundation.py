#!/usr/bin/env python3
"""Standalone Frozen-Output Reproducibility Verifier for Cross-Paper Shared Foundation.

Lightweight verification suite executing in <1-2 seconds locally without cluster compute or GPU.
Validates 6 empirical verification gates:
1. Exact SHA256 checksums of all frozen foundation assets against manifests/sha256_manifest.tsv.
2. Genome manifest: exact 5,438 MAG IDs, all unique (no duplicates), required columns present.
3. Full taxonomy reconciliation: exact 5,438 entries, exact set and order match with genome manifest.
   Canonical GTDB r232 classification syntax: exact seven rank prefixes (d__, p__, c__, o__, f__, g__, s__)
   and strict consistency of canonical rank columns with the classification string.
4. Archaea and Thermoplasmatota subsets & cross-table concordance: exact 591 canonical archaea,
   exact 237 Thermoplasmatota, exact 354 other archaea across verified GTDB phyla.
   Thermoplasmatota 237 subset matches taxonomy_reconciliation_thermo237.tsv 1:1 in ID set and full
   canonical GTDB r232 taxonomy string.
5. Metabolic 484 model inventory exact cohort join: exact 484 unique IDs; exact 237 EQ desert IDs
   (cohort == "EQ_desert") matching Thermoplasmatota 237 1:1; exact 247 external comparator IDs
   (cohort == "external_comparator" explicitly verified); exact disjointness between EQ and external
   cohorts, and disjointness between external IDs and all 5,438 EQ MAGs.
6. External comparator source label verification: exact source labels "GTDB_r232_rep" (218) and
   "MGnify" (29); documented uncertainty regarding cross-release genus/order harmonization.

Informational boundary (not a verification gate):
- Preserves boundary note: frozen inputs form the cross-paper foundation, not proof that ongoing
  metabolic and viral analyses finished. Root owns scientific interpretation.
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Set, Tuple

PKG_ROOT = Path(__file__).resolve().parent
INPUTS_DIR = PKG_ROOT / "inputs"
MANIFESTS_DIR = PKG_ROOT / "manifests"
PROVENANCE_DIR = PKG_ROOT / "provenance"
MANIFEST_P = MANIFESTS_DIR / "sha256_manifest.tsv"

GTDB_RANK_PREFIXES = ("d__", "p__", "c__", "o__", "f__", "g__", "s__")
GTDB_RANK_COLUMNS = (
    "canonical_domain",
    "canonical_phylum",
    "canonical_class",
    "canonical_order",
    "canonical_family",
    "canonical_genus",
    "canonical_species",
)


def sha256_file(p: Path) -> str:
    """Compute SHA256 hex digest of a file."""
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        while chunk := fh.read(65536):
            h.update(chunk)
    return h.hexdigest()


def verify_sha256_manifest(manifest_path: Path, root_dir: Path) -> int:
    """Gate 1: Verify SHA256 integrity of all files listed in manifest."""
    assert manifest_path.exists(), f"Missing manifest: {manifest_path}"
    with open(manifest_path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        rows = list(reader)
    assert len(rows) >= 6, f"Expected at least 6 manifest rows, got {len(rows)}"

    for row in rows:
        expected_sha = row["sha256"]
        rel_p = row["relative_path"]
        expected_size = int(row["file_size_bytes"])
        target_p = root_dir / rel_p
        assert target_p.exists(), f"Manifest file missing on disk: {target_p}"
        actual_size = target_p.stat().st_size
        assert actual_size == expected_size, (
            f"File size mismatch for {rel_p}:\n"
            f"  expected: {expected_size}\n"
            f"  actual:   {actual_size}"
        )
        actual_sha = sha256_file(target_p)
        assert actual_sha == expected_sha, (
            f"Checksum mismatch for {rel_p}:\n"
            f"  expected: {expected_sha}\n"
            f"  actual:   {actual_sha}"
        )
    return len(rows)


def load_tsv_rows(p: Path, is_gzip: bool = False) -> List[Dict[str, str]]:
    """Load TSV file into list of dicts."""
    opener = gzip.open if is_gzip else open
    with opener(p, "rt", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        return list(reader)


def verify_manifest_and_taxonomy_ids(
    manifest_rows: List[Dict[str, str]],
    tax_rows: List[Dict[str, str]],
) -> Tuple[List[str], List[str]]:
    """Gate 2 & 3a: Prove exact 5,438 ID count, zero duplicates, and exact 1:1 ID match."""
    # Gate 2: Genome manifest
    manifest_ids = [r["mag"] for r in manifest_rows]
    assert len(manifest_ids) == 5438, f"Expected exactly 5438 MAGs in manifest, got {len(manifest_ids)}"
    manifest_unique = set(manifest_ids)
    assert len(manifest_unique) == 5438, (
        f"Duplicate MAG IDs detected in manifest: {5438 - len(manifest_unique)} duplicates"
    )

    # Gate 3a: Taxonomy reconciliation
    tax_ids = [r["genome"] for r in tax_rows]
    assert len(tax_ids) == 5438, f"Expected exactly 5438 genomes in taxonomy reconciliation, got {len(tax_ids)}"
    tax_unique = set(tax_ids)
    assert len(tax_unique) == 5438, (
        f"Duplicate genome IDs detected in taxonomy reconciliation: {5438 - len(tax_unique)} duplicates"
    )

    # Exact set equality
    missing_in_tax = manifest_unique - tax_unique
    unexpected_in_tax = tax_unique - manifest_unique
    assert manifest_unique == tax_unique, (
        f"ID mismatch between manifest and taxonomy reconciliation:\n"
        f"  Missing in taxonomy: {len(missing_in_tax)}\n"
        f"  Unexpected in taxonomy: {len(unexpected_in_tax)}"
    )

    # Exact row order correspondence
    assert manifest_ids == tax_ids, "Row order between manifest and taxonomy table differs"

    return manifest_ids, tax_ids


def verify_canonical_taxonomy_syntax_and_ranks(tax_rows: List[Dict[str, str]]) -> None:
    """Gate 3b: Validate exact seven rank prefixes and consistency of canonical rank columns with classification."""
    for col in GTDB_RANK_COLUMNS:
        assert col in tax_rows[0], f"Required taxonomy rank column '{col}' missing from taxonomy reconciliation"
    assert "canonical_gtdb_r232_classification" in tax_rows[0], (
        "Column 'canonical_gtdb_r232_classification' missing from taxonomy reconciliation"
    )

    for idx, r in enumerate(tax_rows):
        gid = r.get("genome", f"row_{idx}")
        tax_str = r.get("canonical_gtdb_r232_classification", "")
        assert tax_str and tax_str.strip(), f"Genome {gid} (row {idx}) has empty canonical taxonomy string"

        parts = tax_str.split(";")
        assert len(parts) == 7, (
            f"Genome {gid} (row {idx}) has {len(parts)} taxonomy ranks; expected exactly 7 (d,p,c,o,f,g,s): '{tax_str}'"
        )

        for pfx, rank_name, rank_col, part in zip(GTDB_RANK_PREFIXES, ("domain", "phylum", "class", "order", "family", "genus", "species"), GTDB_RANK_COLUMNS, parts):
            assert part.startswith(pfx), (
                f"Genome {gid} (row {idx}) {rank_name} rank does not start with prefix '{pfx}': '{part}' in '{tax_str}'"
            )
            col_val = r.get(rank_col, "")
            assert col_val == part, (
                f"Genome {gid} (row {idx}) rank column '{rank_col}' ('{col_val}') inconsistent with lineage rank '{part}'"
            )


def verify_archaea_and_thermo_subsets_and_concordance(
    tax_rows: List[Dict[str, str]],
    thermo_rows: List[Dict[str, str]],
) -> Tuple[Set[str], Set[str], Set[str]]:
    """Gate 4: Verify exact 591 Archaea, exact 237 Thermoplasmatota, 354 other Archaea, and cross-table full taxonomy concordance."""
    archaeal_rows = [r for r in tax_rows if r["canonical_domain"] == "d__Archaea"]
    thermo_in_all = [r for r in tax_rows if r["canonical_phylum"] == "p__Thermoplasmatota"]
    other_archaea = [r for r in archaeal_rows if r["canonical_phylum"] != "p__Thermoplasmatota"]

    assert len(archaeal_rows) == 591, f"Expected exactly 591 archaea, got {len(archaeal_rows)}"
    assert len(thermo_in_all) == 237, f"Expected exactly 237 Thermoplasmatota in all5438, got {len(thermo_in_all)}"
    assert len(other_archaea) == 354, f"Expected exactly 354 other archaea, got {len(other_archaea)}"

    # Check that all Thermoplasmatota are classified as Archaea
    for r in thermo_in_all:
        assert r["canonical_domain"] == "d__Archaea", (
            f"Thermoplasmatota genome {r['genome']} not classified as Archaea"
        )

    # Check taxonomy_reconciliation_thermo237.tsv row count and uniqueness
    thermo_ids_table = [r["genome"] for r in thermo_rows]
    assert len(thermo_ids_table) == 237, f"Expected 237 rows in thermo237 table, got {len(thermo_ids_table)}"
    thermo_unique = set(thermo_ids_table)
    assert len(thermo_unique) == 237, "Duplicate genome IDs found in thermo237 table"

    # Exact ID set match
    thermo_in_all_map = {r["genome"]: r for r in thermo_in_all}
    thermo_in_all_ids = set(thermo_in_all_map.keys())
    assert thermo_in_all_ids == thermo_unique, (
        f"Mismatch between Thermoplasmatota in all5438 and thermo237 table:\n"
        f"  In all5438 only: {len(thermo_in_all_ids - thermo_unique)}\n"
        f"  In thermo237 only: {len(thermo_unique - thermo_in_all_ids)}"
    )

    # Cross-table full taxonomy string and rank-by-rank concordance for every single Thermoplasmatota genome
    for r237 in thermo_rows:
        gid = r237["genome"]
        rall = thermo_in_all_map[gid]
        assert r237["canonical_gtdb_r232_classification"] == rall["canonical_gtdb_r232_classification"], (
            f"Full taxonomy string mismatch for genome {gid} between thermo237 and all5438:\n"
            f"  thermo237: {r237['canonical_gtdb_r232_classification']}\n"
            f"  all5438:   {rall['canonical_gtdb_r232_classification']}"
        )
        for col in GTDB_RANK_COLUMNS:
            assert r237[col] == rall[col], (
                f"Rank column '{col}' mismatch for genome {gid} between thermo237 ('{r237[col]}') and all5438 ('{rall[col]}')"
            )

    return set(r["genome"] for r in archaeal_rows), thermo_unique, set(r["genome"] for r in other_archaea)


def verify_metabolic_cohort_join(
    exp1_rows: List[Dict[str, str]],
    thermo_ids: Set[str],
    all_mag_ids: Set[str],
) -> Tuple[Set[str], Set[str]]:
    """Gate 5: Verify exact 484 model inventory, explicit cohort values (EQ_desert vs external_comparator), and disjointness."""
    assert len(exp1_rows) == 484, f"Expected exactly 484 rows in exp1 inventory, got {len(exp1_rows)}"
    exp1_ids = [r["genome"] for r in exp1_rows]
    assert len(set(exp1_ids)) == 484, f"Duplicate genome IDs found in exp1 inventory: {484 - len(set(exp1_ids))}"

    # Explicit cohort validation: every row must be either "EQ_desert" or "external_comparator"
    for r in exp1_rows:
        cohort = r.get("cohort", "")
        assert cohort in ("EQ_desert", "external_comparator"), (
            f"Genome {r.get('genome')} has unexpected cohort value '{cohort}'; must be 'EQ_desert' or 'external_comparator'"
        )

    eq_rows = [r for r in exp1_rows if r["cohort"] == "EQ_desert"]
    ext_rows = [r for r in exp1_rows if r["cohort"] == "external_comparator"]

    assert len(eq_rows) == 237, f"Expected exactly 237 EQ_desert rows in exp1 inventory, got {len(eq_rows)}"
    assert len(ext_rows) == 247, f"Expected exactly 247 external_comparator rows in exp1 inventory, got {len(ext_rows)}"
    assert len(eq_rows) + len(ext_rows) == 484, "Sum of EQ_desert and external_comparator does not equal 484"

    eq_ids = set(r["genome"] for r in eq_rows)
    ext_ids = set(r["genome"] for r in ext_rows)

    # 1:1 join between EQ_desert cohort and Thermoplasmatota 237
    assert eq_ids == thermo_ids, (
        f"EQ_desert cohort in exp1 inventory does not match canonical Thermoplasmatota 237:\n"
        f"  In exp1 EQ only: {len(eq_ids - thermo_ids)}\n"
        f"  In thermo237 only: {len(thermo_ids - eq_ids)}"
    )

    # Disjointness between EQ and external
    overlap_eq_ext = eq_ids & ext_ids
    assert len(overlap_eq_ext) == 0, f"Overlap between EQ and external cohorts in exp1: {overlap_eq_ext}"

    # Disjointness between external IDs and all 5,438 Empty Quarter MAGs
    overlap_ext_all = ext_ids & all_mag_ids
    assert len(overlap_ext_all) == 0, (
        f"External comparator IDs overlap with Empty Quarter MAG catalogue: {overlap_ext_all}"
    )

    # Validate that all 484 models have draft SBML presence flag True
    for r in exp1_rows:
        assert r.get("has_draft_sbml") == "True", f"Genome {r['genome']} missing draft SBML model flag"

    return eq_ids, ext_ids


def verify_external_source_labels_and_harmonization_boundary(exp1_rows: List[Dict[str, str]]) -> Dict[str, int]:
    """Gate 6: Verify declared external comparator source labels and record documented harmonization uncertainty."""
    ext_rows = [r for r in exp1_rows if r["cohort"] == "external_comparator"]
    source_counts: Dict[str, int] = {}
    for r in ext_rows:
        s = r.get("source", "")
        source_counts[s] = source_counts.get(s, 0) + 1

    expected_sources = {"GTDB_r232_rep", "MGnify"}
    actual_sources = set(source_counts.keys())
    assert actual_sources == expected_sources, (
        f"Expected external source labels {expected_sources}, got {actual_sources}"
    )
    assert source_counts["GTDB_r232_rep"] == 218, (
        f"Expected 218 GTDB_r232_rep sources, got {source_counts['GTDB_r232_rep']}"
    )
    assert source_counts["MGnify"] == 29, (
        f"Expected 29 MGnify sources, got {source_counts['MGnify']}"
    )

    # Confirm all external comparators are classified under Thermoplasmatota in inventory
    for r in ext_rows:
        assert r["canonical_classification"].startswith("d__Archaea;p__Thermoplasmatota"), (
            f"External comparator {r['genome']} not classified as Thermoplasmatota: {r['canonical_classification']}"
        )

    return source_counts


def verify_all_gates(pkg_root: Path = PKG_ROOT) -> Dict[str, Any]:
    """Execute all 6 empirical verification gates on the shared foundation package."""
    inputs_dir = pkg_root / "inputs"
    manifest_p = pkg_root / "manifests" / "sha256_manifest.tsv"

    print("=" * 80)
    print("Cross-Paper Shared Foundation: Exact ID / Taxonomy / Cohort Verifier")
    print(f"Package Root: {pkg_root}")
    print("=" * 80)

    # Gate 1: Checksum integrity
    print("[1/6] Checking SHA256 manifest integrity...")
    n_files = verify_sha256_manifest(manifest_p, pkg_root)
    print(f"  PASS: All {n_files} foundation artifacts verified byte-exact.")

    # Load frozen inputs
    manifest_p_tsv = inputs_dir / "genome_manifest.tsv"
    tax_p_gz = inputs_dir / "taxonomy_reconciliation_all5438.tsv.gz"
    thermo_p_tsv = inputs_dir / "taxonomy_reconciliation_thermo237.tsv"
    exp1_p_tsv = inputs_dir / "exp1_model_inventory_and_qc_484.tsv"

    manifest_rows = load_tsv_rows(manifest_p_tsv)
    tax_rows = load_tsv_rows(tax_p_gz, is_gzip=True)
    thermo_rows = load_tsv_rows(thermo_p_tsv)
    exp1_rows = load_tsv_rows(exp1_p_tsv)

    # Gate 2: Genome manifest IDs
    print("[2/6] Verifying 5,438 genome manifest IDs and uniqueness...")
    manifest_ids, tax_ids = verify_manifest_and_taxonomy_ids(manifest_rows, tax_rows)
    print("  PASS: Exactly 5,438 unique MAG IDs verified with zero duplicates.")

    # Gate 3: Taxonomy reconciliation syntax and rank consistency
    print("[3/6] Verifying canonical GTDB r232 taxonomy syntax and rank consistency...")
    verify_canonical_taxonomy_syntax_and_ranks(tax_rows)
    print(
        "  PASS: All 5,438 genomes have valid GTDB seven-rank prefixes (d__, p__, c__, o__, f__, g__, s__)\n"
        "        and strict consistency between canonical rank columns and classification strings."
    )

    # Gate 4: Archaea and Thermoplasmatota subsets and cross-table concordance
    print("[4/6] Verifying Archaea (591) / Thermoplasmatota (237) subsets & cross-table concordance...")
    archaea_ids, thermo_ids, other_archaea_ids = verify_archaea_and_thermo_subsets_and_concordance(
        tax_rows, thermo_rows
    )
    print(
        f"  PASS: Exact 591 Archaea, exact 237 Thermoplasmatota, and exact 354 other Archaea verified;\n"
        f"        Thermoplasmatota 237 subset matches taxonomy_reconciliation_thermo237.tsv 1:1 in ID set\n"
        f"        and full canonical GTDB r232 taxonomy strings."
    )

    # Gate 5: Metabolic 484 model inventory exact cohort join
    print("[5/6] Verifying metabolic 484 model inventory and exact cohort join...")
    all_mag_ids = set(manifest_ids)
    eq_ids, ext_ids = verify_metabolic_cohort_join(exp1_rows, thermo_ids, all_mag_ids)
    print(
        f"  PASS: Exactly 484 unique models verified:\n"
        f"        - Exact 237 EQ desert IDs (cohort == 'EQ_desert') equal canonical Thermoplasmatota 237 (1:1 join);\n"
        f"        - Exact 247 external comparator IDs (cohort == 'external_comparator') disjoint from EQ IDs and all 5,438 MAGs."
    )

    # Gate 6: External comparator source label verification & harmonization boundary
    print("[6/6] Verifying external comparator source labels and harmonization boundary...")
    source_counts = verify_external_source_labels_and_harmonization_boundary(exp1_rows)
    print(
        f"  PASS: Declared source labels verified: {source_counts['GTDB_r232_rep']} GTDB_r232_rep, {source_counts['MGnify']} MGnify.\n"
        f"        Documented limitation: source labels are checked as declared in inventory, but underlying\n"
        f"        cross-release genus/order harmonization between external comparator genomes and EQ MAGs is not\n"
        f"        independently audited by this join."
    )

    print("=" * 80)
    print("ALL 6 VERIFICATION GATES PASSED CLEANLY (Shared Foundation Verified).")
    print("=" * 80)

    # Informational notice (not a verification gate)
    print("[INFO] Scientific Interpretation Boundary & Pending Status:")
    print("       - Root owns scientific interpretation.")
    print("       - Frozen inputs form the common cross-paper foundation, not proof that downstream analyses finished.")
    print("       - Ongoing metabolic modeling (Exp 1) and viral searches (Exp 3) remain pending independent review.")
    print("=" * 80)

    return {
        "status": "PASS",
        "gates_passed": 6,
        "total_mags": len(manifest_ids),
        "archaeal_mags": len(archaea_ids),
        "thermoplasmatota_mags": len(thermo_ids),
        "other_archaeal_mags": len(other_archaea_ids),
        "metabolic_inventory_total": len(exp1_rows),
        "metabolic_eq_cohort": len(eq_ids),
        "metabolic_external_cohort": len(ext_ids),
        "external_source_counts": source_counts,
    }


def main() -> None:
    try:
        verify_all_gates(PKG_ROOT)
    except AssertionError as e:
        print(f"\nVERIFICATION FAILED: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
