"""Unit, acceptance, and negative mutation/fixture tests for Exp1 484-genome covariate metadata supplement.

Verifies:
  1. Exact 484-genome universe and 1:1 join against frozen model inventory.
  2. Exact 218 GTDB release-232 versioned join with 0 taxonomy mismatches and explicit taxonomy_comparison_type.
  3. Unknown quality method stays empty/NA (29 MGnify genomes keep checkm_* and checkm2_* empty; reported values in reported_quality_* only) and exact usable counts (455 CheckM2 method-family-compatible, 293 in shared orders).
  4. Formal MIMAG flags preserved from official GTDB r232 release only (24 HQ, 194 MQ, 0 LQ) and recorded as 'unknown' (never False or inferred from >=90%/<5%) for EQ and MGnify.
  5. Exact even-n median via statistics.median (90.63% CheckM v1, 88.05% CheckM2 on 218 GTDB reps).
  6. Fixture tests for wrong accession version, unversioned candidate, multi-version stem preservation without overwrite, and unverified assembly identity.
  7. Fixture tests for missing/NaN/inf/out-of-range quality rejection (never defaulting to 0.0) and inventory-vs-source numeric discrepancy rejection.
  8. Fixture tests for duplicate IDs, duplicate MGnify NCBI accessions, malformed TSV columns, invalid explicit CLI paths (no silent fallback), and canonical checksum verification.
"""

from __future__ import annotations

import csv
import gzip
import importlib.util
import json
from pathlib import Path
from typing import Any, Dict, List
import pytest

COV_DIR = Path(__file__).resolve().parent
TSV_PATH = COV_DIR / "exp1_484_covariate_metadata_supplement.tsv"
JSON_PATH = COV_DIR / "exp1_484_covariate_metadata_report.json"
SCRIPT_PATH = COV_DIR / "build_exp1_covariate_metadata_supplement.py"

INV_CANDIDATES = [
    COV_DIR.parent / "compact_evidence" / "exp1_model_inventory_and_qc_484.tsv",
    COV_DIR.parent / "exp1_model_inventory_and_qc_484.tsv",
    Path("mechanism-20261109/review-20261109/exp1_model_inventory_and_qc_484.tsv"),
    Path("/home/leechuck/Documents/papers/empty-quarter/mechanism-20261109/review-20261109/exp1_model_inventory_and_qc_484.tsv"),
]
INV_PATH = next((p for p in INV_CANDIDATES if p.is_file()), None)


def _load_builder_module():
    spec = importlib.util.spec_from_file_location("build_exp1_covariate_metadata_supplement", SCRIPT_PATH)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


builder = _load_builder_module()


@pytest.fixture(scope="module")
def supplement_rows() -> List[Dict[str, str]]:
    assert TSV_PATH.is_file(), f"Missing supplement TSV at {TSV_PATH}"
    with open(TSV_PATH, "r", encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


@pytest.fixture(scope="module")
def report_data() -> Dict[str, Any]:
    assert JSON_PATH.is_file(), f"Missing report JSON at {JSON_PATH}"
    with open(JSON_PATH, "r", encoding="utf-8") as f:
        return json.load(f)


@pytest.fixture(scope="module")
def inventory_rows() -> List[Dict[str, str]]:
    assert INV_PATH is not None and INV_PATH.is_file(), "Missing inventory TSV in expected locations"
    with open(INV_PATH, "r", encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def _write_tsv(path: Path, fieldnames: List[str], rows: List[Dict[str, Any]], gzipped: bool = False) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if gzipped:
        with gzip.open(path, "wt", encoding="utf-8", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
            w.writeheader()
            w.writerows(rows)
    else:
        with open(path, "w", encoding="utf-8", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
            w.writeheader()
            w.writerows(rows)


def _make_valid_mini_fixture(tmp_path: Path) -> Dict[str, Path]:
    """Create a valid 3-genome fixture (1 EQ_MAG, 1 GTDB_r232_rep, 1 MGnify) for mutation testing."""
    inv_path = tmp_path / "inventory.tsv"
    gtdb_path = tmp_path / "gtdb.tsv.gz"
    cm2_path = tmp_path / "checkm2.tsv"
    mgn_comp_path = tmp_path / "mgn_comp.tsv"
    mgn_tk_path = tmp_path / "mgn_tk.tsv"
    out_tsv = tmp_path / "out.tsv"
    out_json = tmp_path / "out.json"

    inv_cols = [
        "genome", "cohort", "source", "canonical_classification",
        "canonical_order", "canonical_family", "canonical_genus",
        "completeness", "contamination", "quality_tier",
    ]
    inv_rows = [
        {
            "genome": "EQ_TEST_1",
            "cohort": "EQ_desert",
            "source": "EQ_MAG",
            "canonical_classification": "d__Archaea;p__Thermoplasmatota;c__SW-10-69-26;o__JACQPN01;f__JAIVGT01;g__JAIVGT01;s__",
            "canonical_order": "o__JACQPN01",
            "canonical_family": "f__JAIVGT01",
            "canonical_genus": "g__JAIVGT01",
            "completeness": "95.5",
            "contamination": "1.2",
            "quality_tier": "HQ",
        },
        {
            "genome": "GCA_001595915.1",
            "cohort": "external_comparator",
            "source": "GTDB_r232_rep",
            "canonical_classification": "d__Archaea;p__Thermoplasmatota;c__E2;o__DHVEG-1;f__DHVEG-1;g__SG8-52-3;s__SG8-52-3 sp001595915",
            "canonical_order": "o__DHVEG-1",
            "canonical_family": "f__DHVEG-1",
            "canonical_genus": "g__SG8-52-3",
            "completeness": "",
            "contamination": "",
            "quality_tier": "MQ_or_unknown",
        },
        {
            "genome": "MGYG000479240",
            "cohort": "external_comparator",
            "source": "MGnify",
            "canonical_classification": "d__Archaea;p__Thermoplasmatota;c__E2;o__UBA9212;f__UBA9212;g__UBA9212;s__UBA9212 sp9212u",
            "canonical_order": "o__UBA9212",
            "canonical_family": "f__UBA9212",
            "canonical_genus": "g__UBA9212",
            "completeness": "76.49",
            "contamination": "1.12",
            "quality_tier": "MQ_or_unknown",
        },
    ]
    _write_tsv(inv_path, inv_cols, inv_rows)

    gtdb_cols = [
        "accession", "gtdb_taxonomy", "checkm_completeness", "checkm_contamination",
        "checkm_marker_count", "checkm_marker_lineage", "checkm_marker_set_count",
        "checkm_strain_heterogeneity", "checkm2_completeness", "checkm2_contamination",
        "checkm2_model", "mimag_high_quality", "mimag_medium_quality",
        "mimag_low_quality", "ncbi_assembly_level",
    ]
    gtdb_rows = [
        {
            "accession": "GB_GCA_001595915.1",
            "gtdb_taxonomy": "d__Archaea;p__Thermoplasmatota;c__E2;o__DHVEG-1;f__DHVEG-1;g__SG8-52-3;s__SG8-52-3 sp001595915",
            "checkm_completeness": "85.2",
            "checkm_contamination": "2.4",
            "checkm_marker_count": "188",
            "checkm_marker_lineage": "p__Euryarchaeota (UID3)",
            "checkm_marker_set_count": "125",
            "checkm_strain_heterogeneity": "0.0",
            "checkm2_completeness": "78.98",
            "checkm2_contamination": "1.06",
            "checkm2_model": "General",
            "mimag_high_quality": "f",
            "mimag_medium_quality": "t",
            "mimag_low_quality": "f",
            "ncbi_assembly_level": "Contig",
        },
        {
            "accession": "GB_GCA_013329505.1",
            "gtdb_taxonomy": "d__Archaea;p__Thermoplasmatota;c__E2;o__UBA9212;f__UBA9212;g__UBA9212;s__UBA9212 sp9212u",
            "checkm_completeness": "76.01",
            "checkm_contamination": "3.6",
            "checkm_marker_count": "188",
            "checkm_marker_lineage": "p__Euryarchaeota (UID3)",
            "checkm_marker_set_count": "125",
            "checkm_strain_heterogeneity": "0.0",
            "checkm2_completeness": "76.49",
            "checkm2_contamination": "1.12",
            "checkm2_model": "General",
            "mimag_high_quality": "f",
            "mimag_medium_quality": "t",
            "mimag_low_quality": "f",
            "ncbi_assembly_level": "Scaffold",
        },
    ]
    _write_tsv(gtdb_path, gtdb_cols, gtdb_rows, gzipped=True)

    cm2_cols = [
        "Name", "Completeness", "Contamination",
        "Completeness_Model_Used", "Genome_Size", "Total_Contigs",
    ]
    cm2_rows = [
        {
            "Name": "EQ_TEST_1",
            "Completeness": "95.5",
            "Contamination": "1.2",
            "Completeness_Model_Used": "Gradient Boost (General Model)",
            "Genome_Size": "2500000",
            "Total_Contigs": "120",
        }
    ]
    _write_tsv(cm2_path, cm2_cols, cm2_rows)

    mgn_comp_cols = [
        "genome_id", "source", "catalogue", "biome", "lineage",
        "ncbi_accession", "completeness", "contamination", "mgnify_type",
        "same_assembly_as_gtdb_rep", "status", "md5", "md5_ok",
    ]
    mgn_comp_rows = [
        {
            "genome_id": "MGYG000479240",
            "source": "MGnify",
            "catalogue": "marine-v2-0",
            "biome": "root:Environmental:Aquatic:Marine",
            "lineage": "d__Archaea;p__Thermoplasmatota;c__E2;o__UBA9212;f__UBA9212;g__UBA9212;s__UBA9212 sp9212u",
            "ncbi_accession": "GCA_013329505",
            "completeness": "76.49",
            "contamination": "1.12",
            "mgnify_type": "MAG",
            "same_assembly_as_gtdb_rep": "no",
            "status": "ok",
            "md5": "362c749a001c80a0ec6cd69a2d3a30ba",
            "md5_ok": "no_published_md5",
        }
    ]
    _write_tsv(mgn_comp_path, mgn_comp_cols, mgn_comp_rows)

    mgn_tk_cols = ["user_genome", "classification", "classification_method"]
    mgn_tk_rows = [
        {
            "user_genome": "MGYG000479240",
            "classification": "d__Archaea;p__Thermoplasmatota;c__E2;o__UBA9212;f__UBA9212;g__UBA9212;s__UBA9212 sp9212u",
            "classification_method": "ani_screen",
        }
    ]
    _write_tsv(mgn_tk_path, mgn_tk_cols, mgn_tk_rows)

    return {
        "inventory": inv_path,
        "gtdb": gtdb_path,
        "cm2": cm2_path,
        "mgn_comp": mgn_comp_path,
        "mgn_tk": mgn_tk_path,
        "out_tsv": out_tsv,
        "out_json": out_json,
    }


def test_exact_484_universe_and_uniqueness(supplement_rows, inventory_rows):
    """Test that supplement table contains exactly 484 unique genomes matching inventory 1:1."""
    assert len(supplement_rows) == 484
    assert len(inventory_rows) == 484

    sup_ids = [r["genome"] for r in supplement_rows]
    inv_ids = [r["genome"] for r in inventory_rows]

    assert len(set(sup_ids)) == 484, "Duplicate genome IDs found in supplement"
    assert sup_ids == inv_ids, "Genome ID sequence or membership does not match inventory 1:1"

    sources = [r["source"] for r in supplement_rows]
    assert sources.count("EQ_MAG") == 237
    assert sources.count("GTDB_r232_rep") == 218
    assert sources.count("MGnify") == 29


def test_gtdb_exact_join_and_taxonomy_comparison_types(supplement_rows, report_data):
    """Test 218 GTDB versioned join and explicit taxonomy comparison types across all 3 sources."""
    gtdb_rows = [r for r in supplement_rows if r["source"] == "GTDB_r232_rep"]
    eq_rows = [r for r in supplement_rows if r["source"] == "EQ_MAG"]
    mgn_rows = [r for r in supplement_rows if r["source"] == "MGnify"]

    assert len(gtdb_rows) == 218
    for r in gtdb_rows:
        assert r["taxonomy_comparison_type"] == "official_gtdb_r232_release_metadata"
        assert r["order_match_inventory_vs_independent_r232"] == "True"
        assert r["full_taxonomy_match_inventory_vs_independent_r232"] == "True"
        assert r["canonical_classification"] == r["independent_r232_classification"]

    # EQ copied inventory is not treated as an independent taxonomy rerun
    assert len(eq_rows) == 237
    for r in eq_rows:
        assert r["taxonomy_comparison_type"] == "inventory_only_no_independent_source"
        assert r["independent_r232_classification"] == ""
        assert r["order_match_inventory_vs_independent_r232"] == "unknown"
        assert r["full_taxonomy_match_inventory_vs_independent_r232"] == "unknown"

    # MGnify has native GTDB-Tk v2.7.2 r232 rerun: 29/29 order match, 21 full match, 8 genus/species diffs
    assert len(mgn_rows) == 29
    for r in mgn_rows:
        assert r["taxonomy_comparison_type"] == "gtdbtk_v2.7.2_r232_classify_wf_rerun"
        assert r["order_match_inventory_vs_independent_r232"] == "True"
        assert r["independent_r232_classification"].startswith("d__Archaea;p__Thermoplasmatota;")
    assert sum(r["full_taxonomy_match_inventory_vs_independent_r232"] == "True" for r in mgn_rows) == 21
    assert sum(r["full_taxonomy_match_inventory_vs_independent_r232"] == "False" for r in mgn_rows) == 8

    assert report_data["gtdb_join"]["target_count"] == 218
    assert report_data["gtdb_join"]["matched_count"] == 218
    assert report_data["gtdb_join"]["taxonomy_mismatches"] == 0


def test_unknown_method_stays_na_and_exact_usable_counts(supplement_rows, report_data):
    """Verify unknown quality method stays empty in method-specific columns and check usable counts."""
    eq_rows = [r for r in supplement_rows if r["source"] == "EQ_MAG"]
    gtdb_rows = [r for r in supplement_rows if r["source"] == "GTDB_r232_rep"]
    mgn_rows = [r for r in supplement_rows if r["source"] == "MGnify"]

    for r in eq_rows:
        assert r["checkm_completeness"] == ""
        assert r["checkm_contamination"] == ""
        assert r["checkm_marker_count"] == ""
        assert r["checkm_completeness_contamination_tier"] == "unknown"
        assert r["checkm2_completeness"] != ""
        assert r["checkm2_version"] == "1.1.0"
        assert r["method_family_checkm2_completeness"] == r["checkm2_completeness"]
        assert r["checkm2_method_family_compatible"] == "True"

    for r in gtdb_rows:
        assert r["reported_quality_completeness"] == ""
        assert r["reported_quality_contamination"] == ""
        assert r["reported_quality_method"] == "none_in_frozen_inventory"
        assert r["checkm_completeness"] != ""
        assert r["checkm2_completeness"] != ""
        assert r["checkm2_version"] == "unknown_gtdb_r232_release"
        assert r["method_family_checkm2_completeness"] == r["checkm2_completeness"]
        assert r["checkm2_method_family_compatible"] == "True"

    # Crucial: MGnify has unknown quality method and only unversioned/missing NCBI accessions;
    # method-specific CheckM and CheckM2 columns MUST remain empty!
    for r in mgn_rows:
        assert r["reported_quality_completeness"] != ""
        assert r["reported_quality_contamination"] != ""
        assert r["reported_quality_method"] == "unknown_mgnify_catalogue_pipeline"
        assert r["checkm_completeness"] == ""
        assert r["checkm_contamination"] == ""
        assert r["checkm2_completeness"] == ""
        assert r["checkm2_contamination"] == ""
        assert r["checkm2_model"] == ""
        assert r["checkm2_version"] == ""
        assert r["checkm2_completeness_contamination_tier"] == "unknown"
        assert r["method_family_checkm2_completeness"] == ""
        assert r["method_family_checkm2_contamination"] == ""
        assert r["checkm2_method_family_compatible"] == "False"

    uc = report_data["usable_counts"]
    assert uc["checkm_v1_usable_total"] == 218
    assert uc["checkm2_method_family_compatible_total"] == 455
    assert uc["checkm2_method_family_compatible_by_source"] == {"EQ_MAG": 237, "GTDB_r232_rep": 218, "MGnify": 0}
    assert uc["checkm2_shared_order_stratum_usable_total"] == 293
    assert uc["checkm2_shared_order_stratum_usable_by_cohort"] == {"EQ_desert": 234, "external_comparator": 59}
    assert uc["reported_quality_present_total"] == 266
    assert uc["formal_mimag_known_total"] == 218


def test_formal_mimag_not_inferred_and_missing_is_unknown(supplement_rows, report_data):
    """Verify formal MIMAG is attributed only to GTDB release and is 'unknown' (not False) for EQ/MGnify."""
    gtdb_rows = [r for r in supplement_rows if r["source"] == "GTDB_r232_rep"]
    non_gtdb_rows = [r for r in supplement_rows if r["source"] != "GTDB_r232_rep"]

    for r in non_gtdb_rows:
        assert r["mimag_high_quality"] == "unknown"
        assert r["mimag_medium_quality"] == "unknown"
        assert r["mimag_low_quality"] == "unknown"
        assert r["mimag_assessment_source"].startswith("unknown_rrna_trna_")

    assert sum(r["mimag_high_quality"] == "True" for r in gtdb_rows) == 24
    assert sum(r["mimag_medium_quality"] == "True" for r in gtdb_rows) == 194
    assert sum(r["mimag_low_quality"] == "True" for r in gtdb_rows) == 0

    # Numeric >=90/<5 counts are tracked separately from formal MIMAG
    assert report_data["gtdb_join"]["official_mimag_counts"] == {"high": 24, "medium": 194, "low": 0}
    assert report_data["gtdb_join"]["numeric_comp90_cont5_counts"] == {"checkm_v1": 108, "checkm2": 82}


def test_checkm_vs_checkm2_separation_and_even_median(supplement_rows, report_data):
    """Test statistics.median on even-n arrays and verify exact 218-genome GTDB medians."""
    # Direct unit test of calc_stats with even n
    even_stats = builder.calc_stats([10.0, 20.0, 30.0, 40.0])
    assert even_stats["n"] == 4
    assert even_stats["median"] == 25.0  # (20 + 30) / 2, not upper-middle 30.0

    # On the 218 GTDB rows (even n=218):
    gtdb_metrics = report_data["quality_metrics"]["GTDB_r232_rep"]
    assert gtdb_metrics["checkm_completeness"]["n"] == 218
    assert gtdb_metrics["checkm_completeness"]["median"] == 90.63
    assert gtdb_metrics["checkm2_completeness"]["n"] == 218
    assert gtdb_metrics["checkm2_completeness"]["median"] == 88.05


def test_fixture_wrong_accession_version_and_unversioned_candidates(tmp_path):
    """Mutation tests for wrong accession version, unversioned candidate, multi-version stems, and unverified identity."""
    # 1. Base mini-fixture: MGnify has unversioned 'GCA_013329505' -> unversioned_candidate_only, CheckM2 stays empty
    paths = _make_valid_mini_fixture(tmp_path / "case1")
    rows, rep = builder.run_pipeline(paths, verify_checksums=False, enforce_484_counts=False)
    mgn_row = next(r for r in rows if r["source"] == "MGnify")
    assert mgn_row["mgnify_ncbi_accession_match_status"] == "unversioned_candidate_only"
    assert mgn_row["candidate_gtdb_r232_accession"] == "GB_GCA_013329505.1"
    assert mgn_row["checkm_completeness"] == ""
    assert mgn_row["checkm2_completeness"] == ""
    assert mgn_row["method_family_checkm2_completeness"] == ""
    assert mgn_row["checkm2_method_family_compatible"] == "False"
    # Even though EQ_TEST_1 has 95.5% completeness, 1.2% contamination, and quality_tier='HQ', formal MIMAG is 'unknown'
    eq_row = next(r for r in rows if r["source"] == "EQ_MAG")
    assert eq_row["mimag_high_quality"] == "unknown"
    assert eq_row["checkm2_completeness_contamination_tier"] == "comp90_cont5_high"

    # 2. Multiple GTDB versions for the same stem ('GB_GCA_013329505.1' and 'GB_GCA_013329505.2') -> neither overwrites
    paths2 = _make_valid_mini_fixture(tmp_path / "case2")
    with gzip.open(paths2["gtdb"], "rt", encoding="utf-8") as f:
        gtdb_rows = list(csv.DictReader(f, delimiter="\t"))
        gtdb_cols = list(gtdb_rows[0].keys())
    second_ver = dict(gtdb_rows[1])
    second_ver["accession"] = "GB_GCA_013329505.2"
    second_ver["checkm2_completeness"] = "99.99"
    gtdb_rows.append(second_ver)
    _write_tsv(paths2["gtdb"], gtdb_cols, gtdb_rows, gzipped=True)

    rows2, _ = builder.run_pipeline(paths2, verify_checksums=False, enforce_484_counts=False)
    mgn_row2 = next(r for r in rows2 if r["source"] == "MGnify")
    assert mgn_row2["mgnify_ncbi_accession_match_status"] == "unversioned_candidate_only"
    assert mgn_row2["candidate_gtdb_r232_accession"] == "GB_GCA_013329505.1;GB_GCA_013329505.2"
    assert mgn_row2["checkm2_completeness"] == ""

    # 3. MGnify comparator has versioned accession with wrong version ('GCA_013329505.2' vs '.1' in GTDB)
    paths3 = _make_valid_mini_fixture(tmp_path / "case3")
    with open(paths3["mgn_comp"], "r", encoding="utf-8") as f:
        mc_rows = list(csv.DictReader(f, delimiter="\t"))
        mc_cols = list(mc_rows[0].keys())
    mc_rows[0]["ncbi_accession"] = "GCA_013329505.2"
    _write_tsv(paths3["mgn_comp"], mc_cols, mc_rows)

    rows3, _ = builder.run_pipeline(paths3, verify_checksums=False, enforce_484_counts=False)
    mgn_row3 = next(r for r in rows3 if r["source"] == "MGnify")
    assert mgn_row3["mgnify_ncbi_accession_match_status"] == "version_mismatch_candidate_only"
    assert mgn_row3["candidate_gtdb_r232_accession"] == "GB_GCA_013329505.1"
    assert mgn_row3["checkm2_completeness"] == ""
    assert mgn_row3["method_family_checkm2_completeness"] == ""

    # 4. MGnify comparator has matching version ('GCA_013329505.1') but unverified assembly identity (md5_ok != 'yes')
    paths4 = _make_valid_mini_fixture(tmp_path / "case4")
    mc_rows[0]["ncbi_accession"] = "GCA_013329505.1"
    mc_rows[0]["md5_ok"] = "no_published_md5"
    _write_tsv(paths4["mgn_comp"], mc_cols, mc_rows)

    rows4, _ = builder.run_pipeline(paths4, verify_checksums=False, enforce_484_counts=False)
    mgn_row4 = next(r for r in rows4 if r["source"] == "MGnify")
    assert mgn_row4["mgnify_ncbi_accession_match_status"] == "version_matched_unverified_assembly_identity_candidate_only"
    assert mgn_row4["checkm2_completeness"] == ""
    assert mgn_row4["checkm2_method_family_compatible"] == "False"

    # 5. GTDB_r232_rep in inventory has wrong version ('GCA_001595915.2' vs 'GB_GCA_001595915.1' in GTDB) -> fails closed
    paths5 = _make_valid_mini_fixture(tmp_path / "case5")
    with open(paths5["inventory"], "r", encoding="utf-8") as f:
        inv_r = list(csv.DictReader(f, delimiter="\t"))
        inv_c = list(inv_r[0].keys())
    inv_r[1]["genome"] = "GCA_001595915.2"
    _write_tsv(paths5["inventory"], inv_c, inv_r)
    with pytest.raises(ValueError, match="Failed to match all GTDB versioned accessions"):
        builder.run_pipeline(paths5, verify_checksums=False, enforce_484_counts=False)

    # 6. Unversioned GTDB_r232_rep in inventory ('GCA_001595915') -> fails closed
    inv_r[1]["genome"] = "GCA_001595915"
    _write_tsv(paths5["inventory"], inv_c, inv_r)
    with pytest.raises(ValueError, match="must be a versioned assembly accession"):
        builder.run_pipeline(paths5, verify_checksums=False, enforce_484_counts=False)


def test_fixture_missing_nonfinite_and_value_mismatch_rejection(tmp_path):
    """Mutation tests verifying missing, NaN, inf, out-of-range, or mismatched quality values fail closed."""
    for bad_val in ("", "NA", "NaN", "inf", "-5.0", "105.0"):
        paths = _make_valid_mini_fixture(tmp_path / f"bad_cm2_{bad_val or 'empty'}")
        with open(paths["cm2"], "r", encoding="utf-8") as f:
            cm2_rows = list(csv.DictReader(f, delimiter="\t"))
            cm2_cols = list(cm2_rows[0].keys())
        cm2_rows[0]["Completeness"] = bad_val
        _write_tsv(paths["cm2"], cm2_cols, cm2_rows)
        with pytest.raises(ValueError, match="Completeness"):
            builder.run_pipeline(paths, verify_checksums=False, enforce_484_counts=False)

    # Discrepancy between inventory completeness and CheckM2 report Completeness -> fails closed
    paths_mis = _make_valid_mini_fixture(tmp_path / "mismatch_cm2")
    with open(paths_mis["cm2"], "r", encoding="utf-8") as f:
        cm2_rows = list(csv.DictReader(f, delimiter="\t"))
        cm2_cols = list(cm2_rows[0].keys())
    cm2_rows[0]["Completeness"] = "91.0"  # Inventory has 95.5
    _write_tsv(paths_mis["cm2"], cm2_cols, cm2_rows)
    with pytest.raises(ValueError, match="Completeness mismatch for EQ_MAG"):
        builder.run_pipeline(paths_mis, verify_checksums=False, enforce_484_counts=False)


def test_fixture_duplicate_ids_malformed_columns_and_bad_explicit_path(tmp_path):
    """Mutation tests for duplicate IDs, duplicate MGnify NCBI accessions, malformed columns, and bad paths."""
    # 1. Duplicate genome ID in inventory
    p_dup_inv = _make_valid_mini_fixture(tmp_path / "dup_inv")
    with open(p_dup_inv["inventory"], "r", encoding="utf-8") as f:
        inv_r = list(csv.DictReader(f, delimiter="\t"))
        inv_c = list(inv_r[0].keys())
    inv_r.append(dict(inv_r[0]))
    _write_tsv(p_dup_inv["inventory"], inv_c, inv_r)
    with pytest.raises(ValueError, match="Duplicate genome ID in inventory"):
        builder.run_pipeline(p_dup_inv, verify_checksums=False, enforce_484_counts=False)

    # 2. Duplicate Name in CheckM2 report
    p_dup_cm2 = _make_valid_mini_fixture(tmp_path / "dup_cm2")
    with open(p_dup_cm2["cm2"], "r", encoding="utf-8") as f:
        cm2_r = list(csv.DictReader(f, delimiter="\t"))
        cm2_c = list(cm2_r[0].keys())
    cm2_r.append(dict(cm2_r[0]))
    _write_tsv(p_dup_cm2["cm2"], cm2_c, cm2_r)
    with pytest.raises(ValueError, match="Duplicate Name in CheckM2 report"):
        builder.run_pipeline(p_dup_cm2, verify_checksums=False, enforce_484_counts=False)

    # 3. Duplicate NCBI accession across two MGnify comparators
    p_dup_mgn = _make_valid_mini_fixture(tmp_path / "dup_mgn")
    with open(p_dup_mgn["inventory"], "r", encoding="utf-8") as f:
        inv_r = list(csv.DictReader(f, delimiter="\t"))
        inv_c = list(inv_r[0].keys())
    inv2 = dict(inv_r[2])
    inv2["genome"] = "MGYG000999999"
    inv_r.append(inv2)
    _write_tsv(p_dup_mgn["inventory"], inv_c, inv_r)

    with open(p_dup_mgn["mgn_comp"], "r", encoding="utf-8") as f:
        mc_r = list(csv.DictReader(f, delimiter="\t"))
        mc_c = list(mc_r[0].keys())
    mc2 = dict(mc_r[0])
    mc2["genome_id"] = "MGYG000999999"  # Same ncbi_accession GCA_013329505
    mc_r.append(mc2)
    _write_tsv(p_dup_mgn["mgn_comp"], mc_c, mc_r)
    with pytest.raises(ValueError, match="Duplicate NCBI accession"):
        builder.run_pipeline(p_dup_mgn, verify_checksums=False, enforce_484_counts=False)

    # 4. Malformed columns (missing required column and ragged row)
    p_bad_col = _make_valid_mini_fixture(tmp_path / "bad_col")
    with open(p_bad_col["cm2"], "w", encoding="utf-8") as f:
        f.write("Name\tCompleteness\nEQ_TEST_1\t95.5\n")
    with pytest.raises(ValueError, match="missing required columns"):
        builder.run_pipeline(p_bad_col, verify_checksums=False, enforce_484_counts=False)

    p_ragged = _make_valid_mini_fixture(tmp_path / "ragged_row")
    with open(p_ragged["cm2"], "a", encoding="utf-8") as f:
        f.write("EQ_TEST_2\t90.0\n")  # Ragged row with too few columns
    with pytest.raises(ValueError, match="malformed column count"):
        builder.run_pipeline(p_ragged, verify_checksums=False, enforce_484_counts=False)

    # 4b. Duplicate stripped versioned accession in GTDB metadata
    p_dup_gtdb = _make_valid_mini_fixture(tmp_path / "dup_gtdb")
    with gzip.open(p_dup_gtdb["gtdb"], "rt", encoding="utf-8") as f:
        g_r = list(csv.DictReader(f, delimiter="\t"))
        g_c = list(g_r[0].keys())
    dup_g = dict(g_r[0])
    dup_g["accession"] = "RS_GCA_001595915.1"  # Strips to same clean_acc GCA_001595915.1
    g_r.append(dup_g)
    _write_tsv(p_dup_gtdb["gtdb"], g_c, g_r, gzipped=True)
    with pytest.raises(ValueError, match="Duplicate accession in GTDB metadata"):
        builder.run_pipeline(p_dup_gtdb, verify_checksums=False, enforce_484_counts=False)

    # 5. Explicitly supplied non-existent path must fail closed without falling back
    bad_args = builder.parse_args(["--inventory-tsv", str(tmp_path / "does_not_exist.tsv")])
    with pytest.raises(FileNotFoundError, match="Explicitly supplied inventory TSV path does not exist"):
        builder.resolve_default_paths(bad_args)

    # 6. Checksum verification rejects mutated file when verify_checksums=True
    p_hash = _make_valid_mini_fixture(tmp_path / "bad_hash")
    with pytest.raises(ValueError, match="Source checksum mismatch"):
        builder.run_pipeline(p_hash, verify_checksums=True, enforce_484_counts=False)


def test_lineage_stratification_and_report_receipts(supplement_rows, report_data):
    """Test shared order counts and cryptographic file receipts in report JSON."""
    shared_eq = [r for r in supplement_rows if r["cohort"] == "EQ_desert" and r["shared_order_stratum"] == "True"]
    shared_ext = [r for r in supplement_rows if r["cohort"] == "external_comparator" and r["shared_order_stratum"] == "True"]

    assert len(shared_eq) == 234
    assert len(shared_ext) == 67
    assert report_data["shared_order_counts"] == {"EQ_desert": 234, "external_comparator": 67, "total": 301}

    receipts = report_data["file_receipts"]
    assert receipts["gtdb_r232"]["md5"] == "0f51fe6f01b4ccbf0d9bab6c07df87af"
    for k, expected_sha in [
        ("inventory", builder.EXPECTED_CANONICAL_SHA256["inventory"]),
        ("gtdb_r232", builder.EXPECTED_CANONICAL_SHA256["gtdb"]),
        ("checkm2_refined", builder.EXPECTED_CANONICAL_SHA256["cm2"]),
        ("mgnify_comparators", builder.EXPECTED_CANONICAL_SHA256["mgn_comp"]),
        ("mgnify_gtdbtk", builder.EXPECTED_CANONICAL_SHA256["mgn_tk"]),
    ]:
        assert receipts[k]["sha256"] == expected_sha
