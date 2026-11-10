#!/usr/bin/env python3
"""Build post-hoc covariate metadata supplement for the Exp1 484-genome metabolic cohort.

This script recovers and documents covariate metadata across all 484 genomes in the
metabolic cohort without mutating the frozen exp1_model_inventory_and_qc_484.tsv inventory:
  1. Exact 218 GTDB r232 species representative join with official ar53_metadata_r232.tsv.gz
     by versioned assembly accession (preserves CheckM v1 and CheckM2 estimates separately;
     0 taxonomy string mismatches; preserves official GTDB MIMAG flags attributed to release).
  2. Provenance mapping for 237 EQ desert MAGs to canonical CheckM2 v1.1.0 quality estimates
     (checkm2_refined_quality_report.tsv), verifying exact numeric agreement with the inventory
     and failing closed on missing, non-finite, or out-of-range values.
  3. Provenance mapping for 29 MGnify external comparators: preserves original MGnify reported
     estimates in separate reported_quality_* columns with unknown method provenance, keeps
     method-specific CheckM/CheckM2 columns empty unless exact versioned accession and assembly
     identity provenance are verified, and records candidate unversioned/version-mismatched
     NCBI-to-GTDB accession audits explicitly.
  4. Separation of numeric completeness/contamination threshold tiers from formal MIMAG
     classifications (which require rRNA and tRNA criteria not present in EQ/MGnify metabolic
     inputs; missing MIMAG flags are recorded as 'unknown', never inferred or set to False).
  5. Explicit lineage stratification and method-family compatibility flags (retaining software
     version and model differences between EQ CheckM2 v1.1.0 and GTDB r232 CheckM2).
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import logging
import math
import re
import statistics
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Set, Tuple

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger("exp1_covariates")

EXPECTED_COHORT_COUNTS = {
    "EQ_MAG": 237,
    "GTDB_r232_rep": 218,
    "MGnify": 29,
    "TOTAL": 484,
}

SHARED_ORDERS = {
    "o__JACQPN01",
    "o__JACPAO01",
    "o__SW-10-69-26",
}

EXPECTED_CANONICAL_SHA256 = {
    "inventory": "044014df6aa1f7bf9b5d4fcaea3a2daea0ef62f3dfcb3c79ef08e54b0347eb2b",
    "gtdb": "d6ddc0f95d269a6e31af8329e0c3c2af0d90d5dad48f5b3e9d820f340951d8c3",
    "cm2": "9cbea621317ed8cd0d780d569b9c8fa2d560ddc9957c2659ba3b13bb2dfdcfe0",
    "mgn_comp": "fa9eacde99aba982756a16903818125b371e2b5dc94afb9ad358c5ab3a999bb6",
    "mgn_tk": "d411d7b1205736bc4a5e3c524b244fc1fe642e99fc628b450a134ad1f93cfe51",
}
EXPECTED_GTDB_MD5 = "0f51fe6f01b4ccbf0d9bab6c07df87af"

REQUIRED_INVENTORY_COLS = {
    "genome",
    "cohort",
    "source",
    "canonical_classification",
    "canonical_order",
    "canonical_family",
    "canonical_genus",
    "completeness",
    "contamination",
    "quality_tier",
}

REQUIRED_GTDB_COLS = {
    "accession",
    "gtdb_taxonomy",
    "checkm_completeness",
    "checkm_contamination",
    "checkm_marker_count",
    "checkm_marker_lineage",
    "checkm_marker_set_count",
    "checkm_strain_heterogeneity",
    "checkm2_completeness",
    "checkm2_contamination",
    "checkm2_model",
    "mimag_high_quality",
    "mimag_medium_quality",
    "mimag_low_quality",
    "ncbi_assembly_level",
}

REQUIRED_CHECKM2_COLS = {
    "Name",
    "Completeness",
    "Contamination",
    "Completeness_Model_Used",
    "Genome_Size",
    "Total_Contigs",
}

REQUIRED_MGNIFY_COMP_COLS = {
    "genome_id",
    "source",
    "catalogue",
    "biome",
    "lineage",
    "ncbi_accession",
    "completeness",
    "contamination",
    "mgnify_type",
    "same_assembly_as_gtdb_rep",
    "status",
    "md5",
    "md5_ok",
}

REQUIRED_MGNIFY_TK_COLS = {
    "user_genome",
    "classification",
    "classification_method",
}

VERSIONED_NCBI_ACC_RE = re.compile(r"^GC[AF]_\d+\.\d+$")
UNVERSIONED_NCBI_ACC_RE = re.compile(r"^GC[AF]_\d+$")


def sha256_file(path: Path) -> str:
    """Compute hex SHA-256 digest of a file."""
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while chunk := fh.read(65536):
            h.update(chunk)
    return h.hexdigest()


def md5_file(path: Path) -> str:
    """Compute hex MD5 digest of a file."""
    h = hashlib.md5()
    with open(path, "rb") as fh:
        while chunk := fh.read(65536):
            h.update(chunk)
    return h.hexdigest()


def clean_gtdb_accession(acc: str) -> str:
    """Strip GB_ or RS_ prefix from GTDB accession while preserving version."""
    if acc.startswith("GB_") or acc.startswith("RS_"):
        return acc[3:]
    return acc


def classify_ncbi_accession_format(acc: str, context_id: str) -> str:
    """Classify an NCBI assembly accession string as 'none', 'versioned', or 'unversioned'."""
    val = (acc or "").strip()
    if not val:
        return "none"
    if VERSIONED_NCBI_ACC_RE.match(val):
        return "versioned"
    if UNVERSIONED_NCBI_ACC_RE.match(val):
        return "unversioned"
    raise ValueError(f"Malformed NCBI assembly accession for {context_id}: {acc!r}")


def parse_required_finite_float(
    raw_val: Optional[str],
    field_name: str,
    record_id: str,
    min_val: float = 0.0,
    max_val: Optional[float] = 100.0,
) -> float:
    """Parse a required numeric field strictly, rejecting missing, NaN, infinite, or out-of-range values."""
    if raw_val is None:
        raise ValueError(f"Missing required numeric field '{field_name}' for {record_id}")
    s = str(raw_val).strip()
    if not s or s.lower() in {"na", "n/a", "nan", "null", "none", "inf", "-inf", "+inf"}:
        raise ValueError(f"Empty or non-finite required numeric field '{field_name}' for {record_id}: {raw_val!r}")
    try:
        val = float(s)
    except ValueError as exc:
        raise ValueError(f"Invalid float in '{field_name}' for {record_id}: {raw_val!r}") from exc
    if not math.isfinite(val):
        raise ValueError(f"Non-finite float in '{field_name}' for {record_id}: {raw_val!r}")
    if val < min_val:
        raise ValueError(f"Value in '{field_name}' for {record_id} below minimum {min_val}: {val}")
    if max_val is not None and val > max_val:
        raise ValueError(f"Value in '{field_name}' for {record_id} above maximum {max_val}: {val}")
    return val


def parse_required_positive_int(raw_val: Optional[str], field_name: str, record_id: str, min_val: int = 0) -> int:
    """Parse a required integer field strictly."""
    if raw_val is None:
        raise ValueError(f"Missing required integer field '{field_name}' for {record_id}")
    s = str(raw_val).strip()
    if not s:
        raise ValueError(f"Empty required integer field '{field_name}' for {record_id}")
    try:
        val = int(s)
    except ValueError as exc:
        raise ValueError(f"Invalid integer in '{field_name}' for {record_id}: {raw_val!r}") from exc
    if val < min_val:
        raise ValueError(f"Integer in '{field_name}' for {record_id} below minimum {min_val}: {val}")
    return val


def parse_gtdb_bool_flag(raw_val: Optional[str], field_name: str, record_id: str) -> str:
    """Parse official GTDB boolean flag ('t'/'f' or 'true'/'false') into 'True'/'False'."""
    if raw_val is None:
        raise ValueError(f"Missing GTDB boolean flag '{field_name}' for {record_id}")
    s = str(raw_val).strip().lower()
    if s in {"t", "true", "1"}:
        return "True"
    if s in {"f", "false", "0"}:
        return "False"
    raise ValueError(f"Invalid GTDB boolean flag '{field_name}' for {record_id}: {raw_val!r}")


def numeric_completeness_contamination_tier(
    comp: Optional[float], cont: Optional[float]
) -> str:
    """Classify numeric completeness/contamination thresholds without inferring formal MIMAG."""
    if comp is None or cont is None:
        return "unknown"
    if comp >= 90.0 and cont < 5.0:
        return "comp90_cont5_high"
    if comp >= 50.0 and cont < 10.0:
        return "comp50_cont10_medium"
    return "comp_below50_or_cont10_low"


def validate_reader_schema_and_row(
    reader: csv.DictReader,
    required_cols: Set[str],
    source_label: str,
) -> None:
    """Verify required columns exist in TSV header."""
    if reader.fieldnames is None:
        raise ValueError(f"{source_label} has no header row")
    missing = required_cols - set(reader.fieldnames)
    if missing:
        raise ValueError(f"{source_label} missing required columns: {sorted(missing)}")


def check_row_well_formed(row: Dict[Any, Any], source_label: str, row_num: int) -> None:
    """Reject ragged TSV rows with extra or missing columns."""
    if None in row or any(v is None for v in row.values()):
        raise ValueError(f"{source_label} has malformed column count at data row {row_num}")


def extract_rank_token(tax: str, rank_prefix: str) -> str:
    """Extract a specific rank token (e.g. 'o__...') from a semicolon-delimited taxonomy string."""
    for part in tax.split(";"):
        part_s = part.strip()
        if part_s.startswith(rank_prefix):
            return part_s
    return ""


def calc_stats(vals: Sequence[float]) -> Dict[str, Any]:
    """Compute descriptive statistics using statistics.median for exact even-n medians."""
    if not vals:
        return {}
    s = sorted(float(x) for x in vals)
    n = len(s)
    mean_val = statistics.mean(s)
    std_val = statistics.pstdev(s) if n > 1 else 0.0
    med_val = statistics.median(s)
    return {
        "n": n,
        "mean": round(mean_val, 2),
        "std": round(std_val, 2),
        "min": round(s[0], 2),
        "median": round(med_val, 2),
        "max": round(s[-1], 2),
    }


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build post-hoc covariate metadata supplement for Exp1 484 metabolic models."
    )
    parser.add_argument(
        "--inventory-tsv",
        type=Path,
        default=None,
        help="Path to exp1_model_inventory_and_qc_484.tsv",
    )
    parser.add_argument(
        "--gtdb-r232-gz",
        type=Path,
        default=None,
        help="Path to official ar53_metadata_r232.tsv.gz",
    )
    parser.add_argument(
        "--checkm2-refined-tsv",
        type=Path,
        default=None,
        help="Path to checkm2_refined_quality_report.tsv",
    )
    parser.add_argument(
        "--mgnify-comparators-tsv",
        type=Path,
        default=None,
        help="Path to v2_thermo_relatives_comparators.tsv",
    )
    parser.add_argument(
        "--mgnify-gtdbtk-tsv",
        type=Path,
        default=None,
        help="Path to v2_thermo_relatives_gtdbtk_mgnify.tsv",
    )
    parser.add_argument(
        "--output-tsv",
        type=Path,
        default=None,
        help="Output TSV path for covariate metadata supplement",
    )
    parser.add_argument(
        "--output-report-json",
        type=Path,
        default=None,
        help="Output JSON path for summary report",
    )
    parser.add_argument(
        "--skip-canonical-checksum-verification",
        action="store_true",
        help="Skip canonical SHA-256/MD5 verification (for synthetic test fixtures only).",
    )
    return parser.parse_args(argv)


def _resolve_single_path(
    explicit_path: Optional[Path],
    default_candidates: Sequence[Path],
    label: str,
) -> Path:
    """Resolve an input path: fail closed if an explicit path is supplied and missing."""
    if explicit_path is not None:
        if not explicit_path.is_file():
            raise FileNotFoundError(f"Explicitly supplied {label} path does not exist: {explicit_path}")
        return explicit_path

    found = next((p for p in default_candidates if p.is_file()), None)
    if found is None:
        raise FileNotFoundError(f"Could not locate {label} in default locations: {[str(p) for p in default_candidates]}")
    return found


def resolve_default_paths(args: argparse.Namespace) -> Dict[str, Path]:
    """Resolve file paths, rejecting any non-existent explicit path without silent fallback."""
    script_dir = Path(__file__).resolve().parent
    repo_root = script_dir.parents[2]  # reproducibility/

    inv_path = _resolve_single_path(
        args.inventory_tsv,
        [
            script_dir.parent / "compact_evidence" / "exp1_model_inventory_and_qc_484.tsv",
            script_dir.parent / "exp1_model_inventory_and_qc_484.tsv",
            Path("mechanism-20261109/review-20261109/exp1_model_inventory_and_qc_484.tsv"),
            Path("/home/leechuck/Documents/papers/empty-quarter/mechanism-20261109/review-20261109/exp1_model_inventory_and_qc_484.tsv"),
        ],
        "inventory TSV",
    )

    gtdb_path = _resolve_single_path(
        args.gtdb_r232_gz,
        [
            Path("mechanism-20261109/review-20261109/gtdb-r232-source/ar53_metadata_r232.tsv.gz"),
            Path("/home/leechuck/Documents/papers/empty-quarter/mechanism-20261109/review-20261109/gtdb-r232-source/ar53_metadata_r232.tsv.gz"),
            script_dir / "ar53_metadata_r232.tsv.gz",
        ],
        "GTDB r232 metadata GZ",
    )

    cm2_path = _resolve_single_path(
        args.checkm2_refined_tsv,
        [
            Path("/home/leechuck/Public/software/empty-quarter-paper-repos/data-paper/evidence/mag-catalogue-20261024/source/checkm2_refined_quality_report.tsv"),
        ],
        "CheckM2 refined quality report TSV",
    )

    mgn_comp_path = _resolve_single_path(
        args.mgnify_comparators_tsv,
        [
            Path("science-20261108/inputs/v2_thermo_relatives_comparators.tsv"),
            Path("/home/leechuck/Documents/papers/empty-quarter/science-20261108/inputs/v2_thermo_relatives_comparators.tsv"),
            repo_root / "reframed_20261108" / "science_20261108_provisional" / "inputs" / "v2_thermo_relatives_comparators.tsv",
        ],
        "MGnify comparators TSV",
    )

    mgn_tk_path = _resolve_single_path(
        args.mgnify_gtdbtk_tsv,
        [
            Path("science-20261108/inputs/v2_thermo_relatives_gtdbtk_mgnify.tsv"),
            Path("/home/leechuck/Documents/papers/empty-quarter/science-20261108/inputs/v2_thermo_relatives_gtdbtk_mgnify.tsv"),
            repo_root / "reframed_20261108" / "science_20261108_provisional" / "inputs" / "v2_thermo_relatives_gtdbtk_mgnify.tsv",
        ],
        "MGnify GTDB-Tk TSV",
    )

    out_tsv = args.output_tsv or (script_dir / "exp1_484_covariate_metadata_supplement.tsv")
    out_json = args.output_report_json or (script_dir / "exp1_484_covariate_metadata_report.json")

    return {
        "inventory": inv_path,
        "gtdb": gtdb_path,
        "cm2": cm2_path,
        "mgn_comp": mgn_comp_path,
        "mgn_tk": mgn_tk_path,
        "out_tsv": out_tsv,
        "out_json": out_json,
    }


def verify_source_checksums(paths: Dict[str, Path]) -> Dict[str, Dict[str, str]]:
    """Compute and verify source file checksums against canonical expectations."""
    receipts: Dict[str, Dict[str, str]] = {}
    for key in ("inventory", "gtdb", "cm2", "mgn_comp", "mgn_tk"):
        p = paths[key]
        actual_sha = sha256_file(p)
        expected_sha = EXPECTED_CANONICAL_SHA256[key]
        if actual_sha != expected_sha:
            raise ValueError(
                f"Source checksum mismatch for {key} ({p}): got sha256={actual_sha}, expected {expected_sha}"
            )
        entry: Dict[str, str] = {"path": str(p), "sha256": actual_sha}
        if key == "gtdb":
            actual_md5 = md5_file(p)
            if actual_md5 != EXPECTED_GTDB_MD5:
                raise ValueError(
                    f"GTDB MD5 mismatch for {p}: got md5={actual_md5}, expected {EXPECTED_GTDB_MD5}"
                )
            entry["md5"] = actual_md5
        receipts[key] = entry
    return receipts


def compute_source_receipts(paths: Dict[str, Path]) -> Dict[str, Dict[str, str]]:
    """Compute source file checksums without enforcing canonical hashes (for custom fixtures)."""
    receipts: Dict[str, Dict[str, str]] = {}
    for key, rkey in (
        ("inventory", "inventory"),
        ("gtdb", "gtdb_r232"),
        ("cm2", "checkm2_refined"),
        ("mgn_comp", "mgnify_comparators"),
        ("mgn_tk", "mgnify_gtdbtk"),
    ):
        p = paths[key]
        entry: Dict[str, str] = {"path": str(p), "sha256": sha256_file(p)}
        if key == "gtdb":
            entry["md5"] = md5_file(p)
        receipts[rkey] = entry
    return receipts


def load_inventory(path: Path, enforce_484_counts: bool = True) -> List[Dict[str, str]]:
    """Load and validate the model inventory, rejecting duplicate IDs or malformed rows."""
    rows: List[Dict[str, str]] = []
    seen_ids: Set[str] = set()
    with open(path, "r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        validate_reader_schema_and_row(reader, REQUIRED_INVENTORY_COLS, f"Inventory ({path})")
        for idx, r in enumerate(reader, start=1):
            check_row_well_formed(r, f"Inventory ({path})", idx)
            gid = r["genome"].strip()
            if not gid:
                raise ValueError(f"Empty genome ID at row {idx} in {path}")
            if gid in seen_ids:
                raise ValueError(f"Duplicate genome ID in inventory {path}: {gid}")
            seen_ids.add(gid)
            rows.append(r)

    if enforce_484_counts:
        if len(rows) != EXPECTED_COHORT_COUNTS["TOTAL"]:
            raise ValueError(
                f"Inventory row count mismatch: got {len(rows)}, expected {EXPECTED_COHORT_COUNTS['TOTAL']}"
            )
        counts: Dict[str, int] = {}
        for r in rows:
            counts[r["source"]] = counts.get(r["source"], 0) + 1
        for src, expected in EXPECTED_COHORT_COUNTS.items():
            if src == "TOTAL":
                continue
            if counts.get(src, 0) != expected:
                raise ValueError(f"Cohort count mismatch for {src}: got {counts.get(src, 0)}, expected {expected}")
        logger.info("Successfully validated 484 inventory: %s", counts)

    return rows


def load_mgnify_comparators(comp_path: Path, mgn_ids: Set[str]) -> Dict[str, Dict[str, Any]]:
    """Load MGnify comparators, rejecting duplicate IDs, duplicate NCBI accessions, or missing quality."""
    comp_data: Dict[str, Dict[str, Any]] = {}
    seen_all_ids: Set[str] = set()
    seen_mgn_ncbi: Dict[str, str] = {}

    with open(comp_path, "r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        validate_reader_schema_and_row(reader, REQUIRED_MGNIFY_COMP_COLS, f"MGnify comparators ({comp_path})")
        for idx, r in enumerate(reader, start=1):
            check_row_well_formed(r, f"MGnify comparators ({comp_path})", idx)
            gid = r["genome_id"].strip()
            if not gid:
                raise ValueError(f"Empty genome_id at row {idx} in {comp_path}")
            if gid in seen_all_ids:
                raise ValueError(f"Duplicate genome_id in {comp_path}: {gid}")
            seen_all_ids.add(gid)

            if gid in mgn_ids:
                if r["source"].strip() != "MGnify":
                    raise ValueError(f"Expected source='MGnify' for {gid} in {comp_path}, got {r['source']!r}")
                ncbi_acc = r.get("ncbi_accession", "").strip()
                acc_format = classify_ncbi_accession_format(ncbi_acc, gid)
                if ncbi_acc:
                    clean_ncbi = clean_gtdb_accession(ncbi_acc)
                    if clean_ncbi in seen_mgn_ncbi:
                        raise ValueError(
                            f"Duplicate NCBI accession {ncbi_acc!r} across MGnify comparators "
                            f"{seen_mgn_ncbi[clean_ncbi]} and {gid}"
                        )
                    seen_mgn_ncbi[clean_ncbi] = gid

                rep_comp = parse_required_finite_float(r.get("completeness"), "completeness", gid, 0.0, 100.0)
                rep_cont = parse_required_finite_float(r.get("contamination"), "contamination", gid, 0.0, 100.0)

                comp_data[gid] = {
                    "catalogue": r.get("catalogue", "").strip(),
                    "biome": r.get("biome", "").strip(),
                    "catalogue_lineage": r.get("lineage", "").strip(),
                    "ncbi_accession": ncbi_acc,
                    "ncbi_accession_format": acc_format,
                    "geographic_origin": r.get("geographic_origin", "").strip(),
                    "reported_completeness": rep_comp,
                    "reported_contamination": rep_cont,
                    "mgnify_type": r.get("mgnify_type", "").strip(),
                    "same_assembly_as_gtdb_rep": r.get("same_assembly_as_gtdb_rep", "").strip(),
                    "status": r.get("status", "").strip(),
                    "md5": r.get("md5", "").strip(),
                    "md5_ok": r.get("md5_ok", "").strip(),
                }

    if len(comp_data) != len(mgn_ids):
        missing = sorted(mgn_ids - set(comp_data.keys()))
        raise ValueError(f"Missing {len(missing)} MGnify genomes in {comp_path}: {missing[:5]}")

    return comp_data


def load_gtdb_metadata(
    path: Path,
    gtdb_target_ids: Set[str],
    mgn_comp: Dict[str, Dict[str, Any]],
) -> Tuple[Dict[str, Dict[str, Any]], Dict[str, Dict[str, Any]]]:
    """Stream official GTDB r232 metadata for 218 exact versioned reps and audit MGnify NCBI candidates.

    Rules enforced:
      - Every GTDB_r232_rep genome must have a versioned accession (GC[AF]_*.version) matching
        ar53_metadata_r232.tsv.gz after stripping only GB_ or RS_.
      - Duplicate raw or stripped versioned accessions in GTDB metadata fail closed.
      - MGnify NCBI accessions are audited without overwriting multiple stem matches:
        * Unversioned accessions (e.g. GCA_013329505) remain 'unversioned_candidate_only'.
        * Versioned accessions with a mismatched version in GTDB (e.g. .2 vs .1) are recorded
          as 'version_mismatch_candidate_only'.
        * Versioned accessions matching a GTDB version without verified assembly sequence identity
          (md5_ok == 'yes' and assembly_identity_verified) remain
          'version_matched_unverified_assembly_identity_candidate_only'.
        * Only a versioned accession with verified exact assembly identity can be promoted to
          'verified_exact_version_and_identity'.
    """
    for gid in gtdb_target_ids:
        if classify_ncbi_accession_format(gid, gid) != "versioned":
            raise ValueError(f"GTDB_r232_rep genome ID must be a versioned assembly accession: {gid!r}")

    # Collect stems of interest from MGnify comparators
    stems_of_interest: Set[str] = set()
    for m_info in mgn_comp.values():
        ncbi_acc = m_info["ncbi_accession"]
        if ncbi_acc:
            stems_of_interest.add(clean_gtdb_accession(ncbi_acc).split(".")[0])

    matched_gtdb: Dict[str, Dict[str, Any]] = {}
    gtdb_by_stem: Dict[str, List[Dict[str, Any]]] = {}
    seen_raw_acc: Set[str] = set()
    seen_clean_acc: Set[str] = set()

    with gzip.open(path, "rt", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        validate_reader_schema_and_row(reader, REQUIRED_GTDB_COLS, f"GTDB metadata ({path})")
        for idx, row in enumerate(reader, start=1):
            check_row_well_formed(row, f"GTDB metadata ({path})", idx)
            raw_acc = row["accession"].strip()
            if not raw_acc:
                raise ValueError(f"Empty accession at row {idx} in {path}")
            if not (raw_acc.startswith("GB_") or raw_acc.startswith("RS_")):
                raise ValueError(f"GTDB accession missing GB_/RS_ prefix at row {idx}: {raw_acc!r}")
            clean_acc = clean_gtdb_accession(raw_acc)
            if classify_ncbi_accession_format(clean_acc, raw_acc) != "versioned":
                raise ValueError(f"GTDB stripped accession is not versioned at row {idx}: {raw_acc!r}")

            if raw_acc in seen_raw_acc or clean_acc in seen_clean_acc:
                raise ValueError(f"Duplicate accession in GTDB metadata {path}: {raw_acc} ({clean_acc})")
            seen_raw_acc.add(raw_acc)
            seen_clean_acc.add(clean_acc)

            stem = clean_acc.split(".")[0]
            if clean_acc in gtdb_target_ids or stem in stems_of_interest:
                rec = {
                    "raw_accession": raw_acc,
                    "clean_accession": clean_acc,
                    "stem": stem,
                    "gtdb_taxonomy": row["gtdb_taxonomy"].strip(),
                    "checkm_completeness": parse_required_finite_float(
                        row.get("checkm_completeness"), "checkm_completeness", raw_acc, 0.0, 100.0
                    ),
                    "checkm_contamination": parse_required_finite_float(
                        row.get("checkm_contamination"), "checkm_contamination", raw_acc, 0.0, 100.0
                    ),
                    "checkm_marker_count": parse_required_positive_int(
                        row.get("checkm_marker_count"), "checkm_marker_count", raw_acc, 0
                    ),
                    "checkm_marker_lineage": row.get("checkm_marker_lineage", "").strip(),
                    "checkm_marker_set_count": parse_required_positive_int(
                        row.get("checkm_marker_set_count"), "checkm_marker_set_count", raw_acc, 0
                    ),
                    "checkm_strain_heterogeneity": parse_required_finite_float(
                        row.get("checkm_strain_heterogeneity"), "checkm_strain_heterogeneity", raw_acc, 0.0, 100.0
                    ),
                    "checkm2_completeness": parse_required_finite_float(
                        row.get("checkm2_completeness"), "checkm2_completeness", raw_acc, 0.0, 100.0
                    ),
                    "checkm2_contamination": parse_required_finite_float(
                        row.get("checkm2_contamination"), "checkm2_contamination", raw_acc, 0.0, 100.0
                    ),
                    "checkm2_model": row.get("checkm2_model", "").strip(),
                    "mimag_high_quality": parse_gtdb_bool_flag(
                        row.get("mimag_high_quality"), "mimag_high_quality", raw_acc
                    ),
                    "mimag_medium_quality": parse_gtdb_bool_flag(
                        row.get("mimag_medium_quality"), "mimag_medium_quality", raw_acc
                    ),
                    "mimag_low_quality": parse_gtdb_bool_flag(
                        row.get("mimag_low_quality"), "mimag_low_quality", raw_acc
                    ),
                    "ncbi_assembly_level": row.get("ncbi_assembly_level", "").strip(),
                }
                if not rec["checkm2_model"]:
                    raise ValueError(f"Empty checkm2_model for {raw_acc} in {path}")
                if not rec["ncbi_assembly_level"]:
                    raise ValueError(f"Empty ncbi_assembly_level for {raw_acc} in {path}")

                if clean_acc in gtdb_target_ids:
                    matched_gtdb[clean_acc] = rec
                if stem in stems_of_interest:
                    gtdb_by_stem.setdefault(stem, []).append(rec)

    if len(matched_gtdb) != len(gtdb_target_ids):
        missing = sorted(gtdb_target_ids - set(matched_gtdb.keys()))
        raise ValueError(f"Failed to match all GTDB versioned accessions: {len(missing)} missing ({missing[:5]})")

    # Resolve MGnify NCBI audit per MGnify genome without overwriting
    mgn_ncbi_audit: Dict[str, Dict[str, Any]] = {}
    for gid, m_info in mgn_comp.items():
        ncbi_acc = m_info["ncbi_accession"]
        acc_fmt = m_info["ncbi_accession_format"]
        if acc_fmt == "none":
            mgn_ncbi_audit[gid] = {
                "match_status": "no_ncbi_accession",
                "candidate_gtdb_accessions": [],
                "verified_gtdb_record": None,
            }
            continue

        clean_ncbi = clean_gtdb_accession(ncbi_acc)
        stem = clean_ncbi.split(".")[0]
        stem_matches = sorted(gtdb_by_stem.get(stem, []), key=lambda x: x["raw_accession"])
        candidate_raw_accs = [r["raw_accession"] for r in stem_matches]

        if acc_fmt == "unversioned":
            if stem_matches:
                mgn_ncbi_audit[gid] = {
                    "match_status": "unversioned_candidate_only",
                    "candidate_gtdb_accessions": candidate_raw_accs,
                    "verified_gtdb_record": None,
                }
            else:
                mgn_ncbi_audit[gid] = {
                    "match_status": "unversioned_no_gtdb_r232_match",
                    "candidate_gtdb_accessions": [],
                    "verified_gtdb_record": None,
                }
        elif acc_fmt == "versioned":
            exact_matches = [r for r in stem_matches if r["clean_accession"] == clean_ncbi]
            if not exact_matches:
                status = (
                    "version_mismatch_candidate_only"
                    if stem_matches
                    else "versioned_no_gtdb_r232_match"
                )
                mgn_ncbi_audit[gid] = {
                    "match_status": status,
                    "candidate_gtdb_accessions": candidate_raw_accs,
                    "verified_gtdb_record": None,
                }
            else:
                # Even when version matches, require explicit verified assembly checksum/identity
                exact_rec = exact_matches[0]
                identity_verified = (
                    m_info.get("md5_ok") == "yes"
                    and m_info.get("same_assembly_as_gtdb_rep") == "yes"
                )
                if identity_verified:
                    mgn_ncbi_audit[gid] = {
                        "match_status": "verified_exact_version_and_identity",
                        "candidate_gtdb_accessions": [exact_rec["raw_accession"]],
                        "verified_gtdb_record": exact_rec,
                    }
                else:
                    mgn_ncbi_audit[gid] = {
                        "match_status": "version_matched_unverified_assembly_identity_candidate_only",
                        "candidate_gtdb_accessions": candidate_raw_accs,
                        "verified_gtdb_record": None,
                    }

    logger.info("Matched all %d GTDB r232 representatives from official metadata", len(matched_gtdb))
    return matched_gtdb, mgn_ncbi_audit


def load_checkm2_refined(path: Path, eq_ids: Set[str]) -> Dict[str, Dict[str, Any]]:
    """Load canonical CheckM2 quality report for EQ MAGs, rejecting missing/non-finite/duplicate rows."""
    data: Dict[str, Dict[str, Any]] = {}
    seen_names: Set[str] = set()

    with open(path, "r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        validate_reader_schema_and_row(reader, REQUIRED_CHECKM2_COLS, f"CheckM2 report ({path})")
        for idx, r in enumerate(reader, start=1):
            check_row_well_formed(r, f"CheckM2 report ({path})", idx)
            name = r["Name"].strip()
            if not name:
                raise ValueError(f"Empty Name at row {idx} in {path}")
            if name in seen_names:
                raise ValueError(f"Duplicate Name in CheckM2 report {path}: {name}")
            seen_names.add(name)

            if name in eq_ids:
                comp = parse_required_finite_float(r.get("Completeness"), "Completeness", name, 0.0, 100.0)
                cont = parse_required_finite_float(r.get("Contamination"), "Contamination", name, 0.0, 100.0)
                model = r.get("Completeness_Model_Used", "").strip()
                if not model or model.lower() in {"none", "na", "null"}:
                    raise ValueError(f"Missing Completeness_Model_Used for {name} in {path}")
                gsize = parse_required_positive_int(r.get("Genome_Size"), "Genome_Size", name, 1)
                tcontigs = parse_required_positive_int(r.get("Total_Contigs"), "Total_Contigs", name, 1)

                data[name] = {
                    "checkm2_completeness": comp,
                    "checkm2_contamination": cont,
                    "checkm2_model": model,
                    "genome_size": gsize,
                    "total_contigs": tcontigs,
                }

    if len(data) != len(eq_ids):
        missing = sorted(eq_ids - set(data.keys()))
        raise ValueError(f"Missing CheckM2 records for {len(missing)} EQ MAGs: {missing[:5]}")

    logger.info("Loaded canonical CheckM2 records for all %d EQ MAGs", len(data))
    return data


def load_mgnify_gtdbtk(tk_path: Path, mgn_ids: Set[str]) -> Dict[str, Dict[str, Any]]:
    """Load native GTDB-Tk r232 classifications for MGnify genomes, rejecting duplicates/missing."""
    tk_data: Dict[str, Dict[str, Any]] = {}
    seen_ug: Set[str] = set()

    with open(tk_path, "r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        validate_reader_schema_and_row(reader, REQUIRED_MGNIFY_TK_COLS, f"MGnify GTDB-Tk ({tk_path})")
        for idx, r in enumerate(reader, start=1):
            check_row_well_formed(r, f"MGnify GTDB-Tk ({tk_path})", idx)
            ug = r["user_genome"].strip()
            if not ug:
                raise ValueError(f"Empty user_genome at row {idx} in {tk_path}")
            if ug in seen_ug:
                raise ValueError(f"Duplicate user_genome in {tk_path}: {ug}")
            seen_ug.add(ug)

            if ug in mgn_ids:
                classification = r["classification"].strip()
                if not classification or not classification.startswith("d__"):
                    raise ValueError(f"Invalid GTDB-Tk classification for {ug} in {tk_path}: {classification!r}")
                tk_data[ug] = {
                    "gtdbtk_r232_classification": classification,
                    "classification_method": r.get("classification_method", "").strip(),
                    "closest_reference": r.get("closest_genome_reference", "").strip(),
                    "closest_ani": r.get("closest_genome_ani", "").strip(),
                }

    if len(tk_data) != len(mgn_ids):
        missing = sorted(mgn_ids - set(tk_data.keys()))
        raise ValueError(f"Missing GTDB-Tk classification for {len(missing)} MGnify genomes: {missing[:5]}")
    return tk_data


def build_covariate_table(
    inventory: List[Dict[str, str]],
    gtdb_meta: Dict[str, Dict[str, Any]],
    eq_meta: Dict[str, Dict[str, Any]],
    mgn_comp: Dict[str, Dict[str, Any]],
    mgn_tk: Dict[str, Dict[str, Any]],
    mgn_ncbi_audit: Dict[str, Dict[str, Any]],
) -> Tuple[List[Dict[str, Any]], Dict[str, Any]]:
    """Assemble the covariate metadata supplement table and audit report."""
    output_rows: List[Dict[str, Any]] = []

    mgn_status_counts: Dict[str, int] = {}
    mgn_candidates_detail: List[Dict[str, Any]] = []

    report: Dict[str, Any] = {
        "total_genomes": len(inventory),
        "cohort_counts": {},
        "usable_counts": {
            "checkm_v1_usable_total": 0,
            "checkm_v1_usable_by_source": {"EQ_MAG": 0, "GTDB_r232_rep": 0, "MGnify": 0},
            "checkm2_method_family_compatible_total": 0,
            "checkm2_method_family_compatible_by_source": {"EQ_MAG": 0, "GTDB_r232_rep": 0, "MGnify": 0},
            "checkm2_shared_order_stratum_usable_total": 0,
            "checkm2_shared_order_stratum_usable_by_cohort": {"EQ_desert": 0, "external_comparator": 0},
            "reported_quality_present_total": 0,
            "formal_mimag_known_total": 0,
        },
        "gtdb_join": {
            "target_count": len(gtdb_meta),
            "matched_count": len(gtdb_meta),
            "taxonomy_mismatches": 0,
            "official_mimag_counts": {"high": 0, "medium": 0, "low": 0},
            "numeric_comp90_cont5_counts": {"checkm_v1": 0, "checkm2": 0},
        },
        "taxonomy_verification_summary": {
            "EQ_MAG": {
                "comparison_type": "inventory_only_no_independent_source",
                "independently_verified": 0,
                "unverified_copied_inventory": 0,
            },
            "GTDB_r232_rep": {
                "comparison_type": "official_gtdb_r232_release_metadata",
                "exact_full_string_matches": 0,
                "order_matches": 0,
                "mismatches": 0,
            },
            "MGnify": {
                "comparison_type": "gtdbtk_v2.7.2_r232_classify_wf_rerun",
                "exact_full_string_matches": 0,
                "order_matches": 0,
                "genus_or_species_differences": 0,
                "order_mismatches": 0,
            },
        },
        "mgnify_provenance": {
            "total_mgnify": len(mgn_comp),
            "verified_exact_version_and_identity": 0,
            "ncbi_accession_match_status_counts": mgn_status_counts,
            "candidate_matches": mgn_candidates_detail,
        },
        "order_breakdown": {"EQ_desert": {}, "external_comparator": {}},
        "shared_order_counts": {"EQ_desert": 0, "external_comparator": 0, "total": 0},
        "quality_metrics": {
            "EQ_desert": {},
            "GTDB_r232_rep": {},
            "MGnify": {},
        },
    }

    for inv in inventory:
        gid = inv["genome"].strip()
        src = inv["source"].strip()
        cohort = inv["cohort"].strip()
        can_tax = inv["canonical_classification"].strip()
        can_order = inv["canonical_order"].strip()

        if extract_rank_token(can_tax, "o__") != can_order:
            raise ValueError(
                f"Inventory canonical_order {can_order!r} does not match canonical_classification {can_tax!r} for {gid}"
            )

        report["cohort_counts"][src] = report["cohort_counts"].get(src, 0) + 1
        report["order_breakdown"][cohort][can_order] = (
            report["order_breakdown"][cohort].get(can_order, 0) + 1
        )

        is_shared_order = can_order in SHARED_ORDERS
        if is_shared_order:
            report["shared_order_counts"][cohort] += 1
            report["shared_order_counts"]["total"] += 1

        row: Dict[str, Any] = {
            "genome": gid,
            "cohort": cohort,
            "source": src,
            "canonical_order": can_order,
            "canonical_family": inv["canonical_family"].strip(),
            "canonical_genus": inv["canonical_genus"].strip(),
            "canonical_classification": can_tax,
            "lineage_stratification_order": can_order,
            "shared_order_stratum": "True" if is_shared_order else "False",
        }

        if src == "EQ_MAG":
            cm2_info = eq_meta[gid]
            inv_comp = parse_required_finite_float(inv.get("completeness"), "inventory.completeness", gid, 0.0, 100.0)
            inv_cont = parse_required_finite_float(inv.get("contamination"), "inventory.contamination", gid, 0.0, 100.0)
            if abs(inv_comp - cm2_info["checkm2_completeness"]) > 1e-4:
                raise ValueError(
                    f"Completeness mismatch for EQ_MAG {gid}: inventory={inv_comp} vs checkm2={cm2_info['checkm2_completeness']}"
                )
            if abs(inv_cont - cm2_info["checkm2_contamination"]) > 1e-4:
                raise ValueError(
                    f"Contamination mismatch for EQ_MAG {gid}: inventory={inv_cont} vs checkm2={cm2_info['checkm2_contamination']}"
                )

            comp_val = cm2_info["checkm2_completeness"]
            cont_val = cm2_info["checkm2_contamination"]
            num_tier = numeric_completeness_contamination_tier(comp_val, cont_val)

            report["taxonomy_verification_summary"]["EQ_MAG"]["unverified_copied_inventory"] += 1
            report["usable_counts"]["checkm2_method_family_compatible_total"] += 1
            report["usable_counts"]["checkm2_method_family_compatible_by_source"]["EQ_MAG"] += 1
            report["usable_counts"]["reported_quality_present_total"] += 1
            if is_shared_order:
                report["usable_counts"]["checkm2_shared_order_stratum_usable_total"] += 1
                report["usable_counts"]["checkm2_shared_order_stratum_usable_by_cohort"]["EQ_desert"] += 1

            row.update({
                "taxonomy_comparison_type": "inventory_only_no_independent_source",
                "independent_r232_classification": "",
                "order_match_inventory_vs_independent_r232": "unknown",
                "full_taxonomy_match_inventory_vs_independent_r232": "unknown",
                "reported_quality_completeness": comp_val,
                "reported_quality_contamination": cont_val,
                "reported_quality_method": "CheckM2_v1.1.0",
                "reported_completeness_contamination_tier": num_tier,
                "checkm_completeness": "",
                "checkm_contamination": "",
                "checkm_marker_count": "",
                "checkm_marker_lineage": "",
                "checkm_marker_set_count": "",
                "checkm_strain_heterogeneity": "",
                "checkm_completeness_contamination_tier": "unknown",
                "checkm2_completeness": comp_val,
                "checkm2_contamination": cont_val,
                "checkm2_model": cm2_info["checkm2_model"],
                "checkm2_version": "1.1.0",
                "checkm2_completeness_contamination_tier": num_tier,
                "method_family_checkm2_completeness": comp_val,
                "method_family_checkm2_contamination": cont_val,
                "checkm2_method_family_compatible": "True",
                "quality_method_sensitivity_flag": "checkm2_v1.1.0_eq_catalogue",
                "mimag_high_quality": "unknown",
                "mimag_medium_quality": "unknown",
                "mimag_low_quality": "unknown",
                "mimag_assessment_source": "unknown_rrna_trna_not_in_metabolic_inventory",
                "ncbi_assembly_level": "unknown",
                "mgnify_catalogue": "",
                "mgnify_biome": "",
                "mgnify_reported_ncbi_accession": "",
                "mgnify_ncbi_accession_match_status": "not_applicable",
                "candidate_gtdb_r232_accession": "",
            })

        elif src == "GTDB_r232_rep":
            if inv.get("completeness", "").strip() != "" or inv.get("contamination", "").strip() != "":
                raise ValueError(f"Expected empty completeness/contamination in frozen inventory for GTDB rep {gid}")

            g_info = gtdb_meta[gid]
            rel_tax = g_info["gtdb_taxonomy"]
            rel_order = extract_rank_token(rel_tax, "o__")
            order_match = can_order == rel_order
            full_tax_match = can_tax == rel_tax
            if not full_tax_match:
                report["gtdb_join"]["taxonomy_mismatches"] += 1
                report["taxonomy_verification_summary"]["GTDB_r232_rep"]["mismatches"] += 1
                raise ValueError(f"Taxonomy mismatch for GTDB rep {gid}: inv={can_tax!r} vs gtdb={rel_tax!r}")
            else:
                report["taxonomy_verification_summary"]["GTDB_r232_rep"]["exact_full_string_matches"] += 1
                report["taxonomy_verification_summary"]["GTDB_r232_rep"]["order_matches"] += 1

            cm1_comp = g_info["checkm_completeness"]
            cm1_cont = g_info["checkm_contamination"]
            cm2_comp = g_info["checkm2_completeness"]
            cm2_cont = g_info["checkm2_contamination"]

            cm1_tier = numeric_completeness_contamination_tier(cm1_comp, cm1_cont)
            cm2_tier = numeric_completeness_contamination_tier(cm2_comp, cm2_cont)

            if cm1_tier == "comp90_cont5_high":
                report["gtdb_join"]["numeric_comp90_cont5_counts"]["checkm_v1"] += 1
            if cm2_tier == "comp90_cont5_high":
                report["gtdb_join"]["numeric_comp90_cont5_counts"]["checkm2"] += 1

            if g_info["mimag_high_quality"] == "True":
                report["gtdb_join"]["official_mimag_counts"]["high"] += 1
            if g_info["mimag_medium_quality"] == "True":
                report["gtdb_join"]["official_mimag_counts"]["medium"] += 1
            if g_info["mimag_low_quality"] == "True":
                report["gtdb_join"]["official_mimag_counts"]["low"] += 1

            report["usable_counts"]["checkm_v1_usable_total"] += 1
            report["usable_counts"]["checkm_v1_usable_by_source"]["GTDB_r232_rep"] += 1
            report["usable_counts"]["checkm2_method_family_compatible_total"] += 1
            report["usable_counts"]["checkm2_method_family_compatible_by_source"]["GTDB_r232_rep"] += 1
            report["usable_counts"]["formal_mimag_known_total"] += 1
            if is_shared_order:
                report["usable_counts"]["checkm2_shared_order_stratum_usable_total"] += 1
                report["usable_counts"]["checkm2_shared_order_stratum_usable_by_cohort"]["external_comparator"] += 1

            row.update({
                "taxonomy_comparison_type": "official_gtdb_r232_release_metadata",
                "independent_r232_classification": rel_tax,
                "order_match_inventory_vs_independent_r232": "True" if order_match else "False",
                "full_taxonomy_match_inventory_vs_independent_r232": "True" if full_tax_match else "False",
                "reported_quality_completeness": "",
                "reported_quality_contamination": "",
                "reported_quality_method": "none_in_frozen_inventory",
                "reported_completeness_contamination_tier": "unknown",
                "checkm_completeness": cm1_comp,
                "checkm_contamination": cm1_cont,
                "checkm_marker_count": g_info["checkm_marker_count"],
                "checkm_marker_lineage": g_info["checkm_marker_lineage"],
                "checkm_marker_set_count": g_info["checkm_marker_set_count"],
                "checkm_strain_heterogeneity": g_info["checkm_strain_heterogeneity"],
                "checkm_completeness_contamination_tier": cm1_tier,
                "checkm2_completeness": cm2_comp,
                "checkm2_contamination": cm2_cont,
                "checkm2_model": g_info["checkm2_model"],
                "checkm2_version": "unknown_gtdb_r232_release",
                "checkm2_completeness_contamination_tier": cm2_tier,
                "method_family_checkm2_completeness": cm2_comp,
                "method_family_checkm2_contamination": cm2_cont,
                "checkm2_method_family_compatible": "True",
                "quality_method_sensitivity_flag": "gtdb_r232_official_checkm2_version_unrecorded",
                "mimag_high_quality": g_info["mimag_high_quality"],
                "mimag_medium_quality": g_info["mimag_medium_quality"],
                "mimag_low_quality": g_info["mimag_low_quality"],
                "mimag_assessment_source": "GTDB_r232_official_release",
                "ncbi_assembly_level": g_info["ncbi_assembly_level"],
                "mgnify_catalogue": "",
                "mgnify_biome": "",
                "mgnify_reported_ncbi_accession": "",
                "mgnify_ncbi_accession_match_status": "not_applicable",
                "candidate_gtdb_r232_accession": "",
            })

        elif src == "MGnify":
            m_info = mgn_comp[gid]
            t_info = mgn_tk[gid]
            audit_info = mgn_ncbi_audit[gid]

            inv_comp = parse_required_finite_float(inv.get("completeness"), "inventory.completeness", gid, 0.0, 100.0)
            inv_cont = parse_required_finite_float(inv.get("contamination"), "inventory.contamination", gid, 0.0, 100.0)
            if abs(inv_comp - m_info["reported_completeness"]) > 1e-4:
                raise ValueError(
                    f"Completeness mismatch for MGnify {gid}: inventory={inv_comp} vs comparators={m_info['reported_completeness']}"
                )
            if abs(inv_cont - m_info["reported_contamination"]) > 1e-4:
                raise ValueError(
                    f"Contamination mismatch for MGnify {gid}: inventory={inv_cont} vs comparators={m_info['reported_contamination']}"
                )

            tk_tax = t_info["gtdbtk_r232_classification"]
            tk_order = extract_rank_token(tk_tax, "o__")
            order_match = can_order == tk_order
            full_tax_match = can_tax == tk_tax

            if order_match:
                report["taxonomy_verification_summary"]["MGnify"]["order_matches"] += 1
            else:
                report["taxonomy_verification_summary"]["MGnify"]["order_mismatches"] += 1
                raise ValueError(f"Order mismatch for MGnify {gid}: inv={can_order!r} vs gtdbtk={tk_order!r}")

            if full_tax_match:
                report["taxonomy_verification_summary"]["MGnify"]["exact_full_string_matches"] += 1
            else:
                report["taxonomy_verification_summary"]["MGnify"]["genus_or_species_differences"] += 1

            match_status = audit_info["match_status"]
            mgn_status_counts[match_status] = mgn_status_counts.get(match_status, 0) + 1
            cand_accs = audit_info["candidate_gtdb_accessions"]
            cand_str = ";".join(cand_accs)
            if cand_accs:
                mgn_candidates_detail.append({
                    "genome": gid,
                    "reported_ncbi_accession": m_info["ncbi_accession"],
                    "match_status": match_status,
                    "candidate_gtdb_r232_accessions": cand_accs,
                    "mgnify_reported_completeness": m_info["reported_completeness"],
                    "mgnify_reported_contamination": m_info["reported_contamination"],
                })

            rep_comp = m_info["reported_completeness"]
            rep_cont = m_info["reported_contamination"]
            rep_tier = numeric_completeness_contamination_tier(rep_comp, rep_cont)
            report["usable_counts"]["reported_quality_present_total"] += 1

            verified_rec = audit_info["verified_gtdb_record"]
            if verified_rec is not None:
                report["mgnify_provenance"]["verified_exact_version_and_identity"] += 1
                cm1_comp_v: Any = verified_rec["checkm_completeness"]
                cm1_cont_v: Any = verified_rec["checkm_contamination"]
                cm1_mc_v: Any = verified_rec["checkm_marker_count"]
                cm1_ml_v = verified_rec["checkm_marker_lineage"]
                cm1_msc_v: Any = verified_rec["checkm_marker_set_count"]
                cm1_sh_v: Any = verified_rec["checkm_strain_heterogeneity"]
                cm1_tier_v = numeric_completeness_contamination_tier(cm1_comp_v, cm1_cont_v)

                cm2_comp_v: Any = verified_rec["checkm2_completeness"]
                cm2_cont_v: Any = verified_rec["checkm2_contamination"]
                cm2_mod_v = verified_rec["checkm2_model"]
                cm2_ver_v = "unknown_gtdb_r232_release"
                cm2_tier_v = numeric_completeness_contamination_tier(cm2_comp_v, cm2_cont_v)
                mf_comp_v: Any = cm2_comp_v
                mf_cont_v: Any = cm2_cont_v
                mf_ok_v = "True"
                sens_flag = "verified_exact_version_and_identity_gtdb_r232_checkm2"

                mimag_hq_v = verified_rec["mimag_high_quality"]
                mimag_mq_v = verified_rec["mimag_medium_quality"]
                mimag_lq_v = verified_rec["mimag_low_quality"]
                mimag_src_v = "GTDB_r232_official_release"
                asm_lvl_v = verified_rec["ncbi_assembly_level"]

                report["usable_counts"]["checkm_v1_usable_total"] += 1
                report["usable_counts"]["checkm_v1_usable_by_source"]["MGnify"] += 1
                report["usable_counts"]["checkm2_method_family_compatible_total"] += 1
                report["usable_counts"]["checkm2_method_family_compatible_by_source"]["MGnify"] += 1
                report["usable_counts"]["formal_mimag_known_total"] += 1
                if is_shared_order:
                    report["usable_counts"]["checkm2_shared_order_stratum_usable_total"] += 1
                    report["usable_counts"]["checkm2_shared_order_stratum_usable_by_cohort"]["external_comparator"] += 1
            else:
                # Unknown method in MGnify and no verified exact assembly join: leave method-specific columns empty
                cm1_comp_v = ""
                cm1_cont_v = ""
                cm1_mc_v = ""
                cm1_ml_v = ""
                cm1_msc_v = ""
                cm1_sh_v = ""
                cm1_tier_v = "unknown"

                cm2_comp_v = ""
                cm2_cont_v = ""
                cm2_mod_v = ""
                cm2_ver_v = ""
                cm2_tier_v = "unknown"
                mf_comp_v = ""
                mf_cont_v = ""
                mf_ok_v = "False"
                sens_flag = f"mgnify_reported_method_unknown_{match_status}"

                mimag_hq_v = "unknown"
                mimag_mq_v = "unknown"
                mimag_lq_v = "unknown"
                mimag_src_v = "unknown_rrna_trna_not_in_mgnify_comparators"
                asm_lvl_v = "unknown"

            row.update({
                "taxonomy_comparison_type": "gtdbtk_v2.7.2_r232_classify_wf_rerun",
                "independent_r232_classification": tk_tax,
                "order_match_inventory_vs_independent_r232": "True" if order_match else "False",
                "full_taxonomy_match_inventory_vs_independent_r232": "True" if full_tax_match else "False",
                "reported_quality_completeness": rep_comp,
                "reported_quality_contamination": rep_cont,
                "reported_quality_method": "unknown_mgnify_catalogue_pipeline",
                "reported_completeness_contamination_tier": rep_tier,
                "checkm_completeness": cm1_comp_v,
                "checkm_contamination": cm1_cont_v,
                "checkm_marker_count": cm1_mc_v,
                "checkm_marker_lineage": cm1_ml_v,
                "checkm_marker_set_count": cm1_msc_v,
                "checkm_strain_heterogeneity": cm1_sh_v,
                "checkm_completeness_contamination_tier": cm1_tier_v,
                "checkm2_completeness": cm2_comp_v,
                "checkm2_contamination": cm2_cont_v,
                "checkm2_model": cm2_mod_v,
                "checkm2_version": cm2_ver_v,
                "checkm2_completeness_contamination_tier": cm2_tier_v,
                "method_family_checkm2_completeness": mf_comp_v,
                "method_family_checkm2_contamination": mf_cont_v,
                "checkm2_method_family_compatible": mf_ok_v,
                "quality_method_sensitivity_flag": sens_flag,
                "mimag_high_quality": mimag_hq_v,
                "mimag_medium_quality": mimag_mq_v,
                "mimag_low_quality": mimag_lq_v,
                "mimag_assessment_source": mimag_src_v,
                "ncbi_assembly_level": asm_lvl_v,
                "mgnify_catalogue": m_info["catalogue"],
                "mgnify_biome": m_info["biome"],
                "mgnify_reported_ncbi_accession": m_info["ncbi_accession"],
                "mgnify_ncbi_accession_match_status": match_status,
                "candidate_gtdb_r232_accession": cand_str,
            })

        else:
            raise ValueError(f"Unexpected source {src!r} for genome {gid}")

        output_rows.append(row)

    # Summary statistics using statistics.median
    eq_comp = [
        float(r["checkm2_completeness"])
        for r in output_rows
        if r["source"] == "EQ_MAG" and r["checkm2_completeness"] != ""
    ]
    eq_cont = [
        float(r["checkm2_contamination"])
        for r in output_rows
        if r["source"] == "EQ_MAG" and r["checkm2_contamination"] != ""
    ]
    report["quality_metrics"]["EQ_desert"]["checkm2_completeness"] = calc_stats(eq_comp)
    report["quality_metrics"]["EQ_desert"]["checkm2_contamination"] = calc_stats(eq_cont)

    gtdb_cm1_comp = [
        float(r["checkm_completeness"])
        for r in output_rows
        if r["source"] == "GTDB_r232_rep" and r["checkm_completeness"] != ""
    ]
    gtdb_cm1_cont = [
        float(r["checkm_contamination"])
        for r in output_rows
        if r["source"] == "GTDB_r232_rep" and r["checkm_contamination"] != ""
    ]
    gtdb_cm2_comp = [
        float(r["checkm2_completeness"])
        for r in output_rows
        if r["source"] == "GTDB_r232_rep" and r["checkm2_completeness"] != ""
    ]
    gtdb_cm2_cont = [
        float(r["checkm2_contamination"])
        for r in output_rows
        if r["source"] == "GTDB_r232_rep" and r["checkm2_contamination"] != ""
    ]
    report["quality_metrics"]["GTDB_r232_rep"]["checkm_completeness"] = calc_stats(gtdb_cm1_comp)
    report["quality_metrics"]["GTDB_r232_rep"]["checkm_contamination"] = calc_stats(gtdb_cm1_cont)
    report["quality_metrics"]["GTDB_r232_rep"]["checkm2_completeness"] = calc_stats(gtdb_cm2_comp)
    report["quality_metrics"]["GTDB_r232_rep"]["checkm2_contamination"] = calc_stats(gtdb_cm2_cont)

    mgn_rep_comp = [
        float(r["reported_quality_completeness"])
        for r in output_rows
        if r["source"] == "MGnify" and r["reported_quality_completeness"] != ""
    ]
    mgn_rep_cont = [
        float(r["reported_quality_contamination"])
        for r in output_rows
        if r["source"] == "MGnify" and r["reported_quality_contamination"] != ""
    ]
    report["quality_metrics"]["MGnify"]["reported_quality_completeness"] = calc_stats(mgn_rep_comp)
    report["quality_metrics"]["MGnify"]["reported_quality_contamination"] = calc_stats(mgn_rep_cont)

    return output_rows, report


def run_pipeline(
    paths: Dict[str, Path],
    verify_checksums: bool = True,
    enforce_484_counts: bool = True,
) -> Tuple[List[Dict[str, Any]], Dict[str, Any]]:
    """Execute full validation and build pipeline, returning rows and report dictionary."""
    if verify_checksums:
        verified = verify_source_checksums(paths)
        receipts = {
            "inventory": verified["inventory"],
            "gtdb_r232": verified["gtdb"],
            "checkm2_refined": verified["cm2"],
            "mgnify_comparators": verified["mgn_comp"],
            "mgnify_gtdbtk": verified["mgn_tk"],
        }
    else:
        receipts = compute_source_receipts(paths)

    inventory = load_inventory(paths["inventory"], enforce_484_counts=enforce_484_counts)

    gtdb_ids = {r["genome"].strip() for r in inventory if r["source"].strip() == "GTDB_r232_rep"}
    eq_ids = {r["genome"].strip() for r in inventory if r["source"].strip() == "EQ_MAG"}
    mgn_ids = {r["genome"].strip() for r in inventory if r["source"].strip() == "MGnify"}

    mgn_comp = load_mgnify_comparators(paths["mgn_comp"], mgn_ids) if mgn_ids else {}
    gtdb_meta, mgn_ncbi_audit = load_gtdb_metadata(paths["gtdb"], gtdb_ids, mgn_comp)
    eq_meta = load_checkm2_refined(paths["cm2"], eq_ids) if eq_ids else {}
    mgn_tk = load_mgnify_gtdbtk(paths["mgn_tk"], mgn_ids) if mgn_ids else {}

    supplement_rows, report = build_covariate_table(
        inventory, gtdb_meta, eq_meta, mgn_comp, mgn_tk, mgn_ncbi_audit
    )
    report["file_receipts"] = receipts

    paths["out_tsv"].parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(supplement_rows[0].keys())
    with open(paths["out_tsv"], "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(supplement_rows)
    logger.info("Wrote %d rows to %s", len(supplement_rows), paths["out_tsv"])

    paths["out_json"].parent.mkdir(parents=True, exist_ok=True)
    with open(paths["out_json"], "w", encoding="utf-8") as f:
        json.dump(report, f, indent=2)
    logger.info("Wrote audit report to %s", len(supplement_rows) and paths["out_json"])

    return supplement_rows, report


def main(argv: Optional[Sequence[str]] = None) -> None:
    args = parse_args(argv)
    paths = resolve_default_paths(args)

    logger.info("Inventory path: %s", paths["inventory"])
    logger.info("GTDB r232 metadata: %s", paths["gtdb"])
    logger.info("CheckM2 refined path: %s", paths["cm2"])
    logger.info("MGnify comparators: %s", paths["mgn_comp"])
    logger.info("MGnify GTDB-Tk: %s", paths["mgn_tk"])
    logger.info("Output TSV: %s", paths["out_tsv"])
    logger.info("Output JSON: %s", paths["out_json"])

    run_pipeline(
        paths,
        verify_checksums=not args.skip_canonical_checksum_verification,
        enforce_484_counts=True,
    )


if __name__ == "__main__":
    main()
