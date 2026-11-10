#!/usr/bin/env python3
"""Extract bounded raw (MAPQ>=0) vs filtered (MAPQ>=20) coverage QC across all 14 samples x 7 focal viral targets.

Fail-closed implementation contract (EXP3-RAW-COVERAGE-WORKER-CONTRACT.txt & EXP3-RAW-COVERAGE-ROOT-REVIEW.txt):
* verify exact 14-sample scope from provenance/stage4_competitive_pipeline_submission.json;
* require actual non-empty stage4_reference_validation_receipt.json and REFERENCE_VALIDATION_DONE sentinel,
  verify sentinel receipt_sha256 binding, completed_verified/schema_version==3, and per-sample reference_binding hashes;
* verify per-sample MAPPING_DONE cryptographic binding to <sample>_mapping_receipt.json, require all producer
  mapping_metrics keys (never defaulting missing metrics to 0), and verify SHA256 integrity of raw BAM,
  filtered BAM, filtered coverage TSV, and focal evaluation TSV;
* run actual `samtools coverage --ff 3844 -q 0 -r <reference_contig> <raw_bam>` for the 7 exact
  viral target reference headers in each of the 14 samples (98 rows), using identical flag exclusion
  mask 3844 (UNMAP,SECONDARY,QCFAIL,DUP,SUPPLEMENTARY) as the primary MAPQ>=20 filter;
* cross-verify `samtools coverage --ff 3844 -q 20 -r <reference_contig> <filtered_bam>` against
  the receipted `focal_evaluations` and `<sample>_filtered_coverage.tsv`;
* join existing Stage 2 vOTU metadata (`results/stage2_target_viral_identity_table.tsv`) as
  context only (never as a demonstrated cause of read-mapping ambiguity; `is_votu_cluster_rep=0`
  is normal for a non-representative target whose `cluster_votu_rep_id` points elsewhere, and
  exact-sequence aliases can differ in representative flags);
* preserve primary 10x / 50% / MAPQ>=20 gates and analysis_role exclusions (`raw_reaches_10x_50pct_qc_only`
  is diagnostic QC only, never an eligibility override);
* distinguish high raw support lost at MAPQ>=20 (exact-reference attribution unresolved; compatible with
  competitive ambiguity or weak alignment discrimination) from low total read support and from
  measured zero assigned reads (`measured_zero_assigned_reads`, never proved biological absence;
  missing is never converted to zero).
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import json
import math
import os
import platform
import shutil
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

import pandas as pd

from stage4_competitive_common import (
     SHA256_HEX,
    atomic_json,
    atomic_text,
    load_analysis_pairs,
    require_columns,
    stream_sha256,
)
from stage4_map_sample_competitive import COVERAGE_COLUMNS, compute_coverage_metrics


def probe_tool_version(cmd: list[str]) -> str:
    result = subprocess.run(cmd, capture_output=True, text=True, errors="replace")
    if result.returncode != 0:
        raise RuntimeError(f"tool version probe failed: {cmd}: {result.stderr}")
    lines = (result.stdout or result.stderr).strip().splitlines()
    if not lines:
        raise RuntimeError(f"empty tool version output: {cmd}")
    return lines[0].strip()

BRANCH_ROOT = Path(os.environ.get("EXP3_ROOT", "/scratch/project/k10018/eq-mechanism-20261109/exp3_viruses"))
RES_DIR = BRANCH_ROOT / "results"
MAP_DIR = RES_DIR / "stage4_competitive_mapping"
REF_DIR = RES_DIR / "stage4_competitive_reference"
OUT_DIR = RES_DIR / "stage4_raw_vs_mapq20"

FILTER_FLAGS_INT = 3844
FILTER_FLAGS_NAMES = ["UNMAP", "SECONDARY", "QCFAIL", "DUP", "SUPPLEMENTARY"]
EXPECTED_MAPPING_FILTER_STRING = "MAPQ>=20; exclude UNMAP,SECONDARY,QCFAIL,DUP,SUPPLEMENTARY"
REQUIRED_MAPPING_METRICS_KEYS = {
    "raw_total_records",
    "raw_primary_mapped_records",
    "filtered_high_quality_mapq20_reads",
    "filtered_fraction_of_primary",
    "coverage_references_verified",
    "filter",
}

AUTHORITATIVE_9SR2_BENCHMARK = {
    "sample": "9Sr2",
    "entity_id": "9Sr2::9Sr2_284237",
    "reference_contig": "VIRUS_TARGET__9Sr2::9Sr2_284237",
    "length_bp": 56873,
    "raw_numreads": 12065,
    "raw_covbases": 56860,
    "raw_coverage_samtools_pct": 99.9771,
    "raw_meandepth": 30.919,
    "raw_meanbaseq": 35.6,
    "raw_meanmapq": 1.51,
    "filtered_numreads": 181,
    "filtered_covbases": 2124,
    "filtered_horizontal_breadth_pct": 3.7346368223937545,
    "filtered_meandepth": 0.469854,
}

STAGE2_REQUIRED_COLUMNS = [
    "target_id",
    "target_seq_sha256",
    "annotation_provenance",
    "cluster_votu_rep_id",
    "is_votu_cluster_rep",
    "target_namespaced_contig_id",
    "viral_identity_status",
    "checkv_quality",
    "checkv_viral_genes",
    "checkv_host_genes",
    "virus_hallmarks",
]


def require_valid_sha256_hex(value: Any, label: str) -> str:
    if not isinstance(value, str) or len(value) != 64 or not (set(value) <= SHA256_HEX):
        raise ValueError(f"invalid SHA256 for {label}: {value!r}")
    return value


def require_strict_int(value: Any, label: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise ValueError(f"expected integer for {label}, got {value!r}")
    return value


def load_expected_samples(submission_path: Path) -> tuple[list[str], dict[str, Any]]:
    if not submission_path.exists() or submission_path.stat().st_size == 0:
        raise FileNotFoundError(f"pipeline submission provenance missing/empty: {submission_path}")
    payload = json.loads(submission_path.read_text())
    samples = payload.get("samples")
    if not isinstance(samples, list) or not samples:
        raise ValueError(f"submission receipt missing non-empty samples list: {submission_path}")
    cleaned = [str(s).strip() for s in samples]
    if any(not s for s in cleaned) or len(cleaned) != len(set(cleaned)):
        raise ValueError(f"invalid or duplicate sample names in submission receipt: {cleaned}")
    if len(cleaned) != 14:
        raise ValueError(f"expected exact 14-sample scope in submission receipt, got {len(cleaned)}: {cleaned}")
    return cleaned, payload


def verify_reference_validation_contract(
    ref_dir: Path,
    expected_pair_ledger_sha256: str,
    expected_boundary_config_sha256: str,
) -> dict[str, Any]:
    """Require actual non-empty reference validation receipt + sentinel and validate bindings."""
    val_receipt_path = ref_dir / "stage4_reference_validation_receipt.json"
    val_done_path = ref_dir / "REFERENCE_VALIDATION_DONE"
    if not val_receipt_path.is_file() or val_receipt_path.stat().st_size == 0:
        raise FileNotFoundError(f"reference validation receipt missing/empty: {val_receipt_path}")
    if not val_done_path.is_file() or val_done_path.stat().st_size == 0:
        raise FileNotFoundError(f"REFERENCE_VALIDATION_DONE sentinel missing/empty: {val_done_path}")

    observed_val_sha = stream_sha256(val_receipt_path)
    done_payload = json.loads(val_done_path.read_text())
    if done_payload.get("receipt_sha256") != observed_val_sha:
        raise ValueError(
            f"REFERENCE_VALIDATION_DONE receipt_sha256 mismatch: "
            f"{done_payload.get('receipt_sha256')} != {observed_val_sha}"
        )

    val_receipt = json.loads(val_receipt_path.read_text())
    if val_receipt.get("status") != "completed_verified" or val_receipt.get("schema_version") != 3:
        raise ValueError(f"invalid reference validation receipt status/schema: {val_receipt_path}")
    if val_receipt.get("stage") != "stage4_validate_competitive_reference":
        raise ValueError(f"unexpected stage in reference validation receipt: {val_receipt.get('stage')!r}")

    if val_receipt.get("pair_ledger", {}).get("sha256") != expected_pair_ledger_sha256:
        raise ValueError("reference validation receipt pair_ledger.sha256 mismatch")
    if val_receipt.get("boundary_config", {}).get("sha256") != expected_boundary_config_sha256:
        raise ValueError("reference validation receipt boundary_config.sha256 mismatch")

    ref_files = val_receipt.get("reference_files", {})
    ref_fna_sha = require_valid_sha256_hex(
        ref_files.get("reference_fna", {}).get("sha256"),
        "reference_files.reference_fna.sha256",
    )
    ref_fai_sha = require_valid_sha256_hex(
        ref_files.get("reference_fai", {}).get("sha256"),
        "reference_files.reference_fai.sha256",
    )
    crosswalk_sha = require_valid_sha256_hex(
        ref_files.get("crosswalk_tsv", {}).get("sha256"),
        "reference_files.crosswalk_tsv.sha256",
    )
    build_receipt_sha = require_valid_sha256_hex(
        val_receipt.get("build_receipt", {}).get("sha256"),
        "build_receipt.sha256",
    )
    index_list = val_receipt.get("bowtie2_index", {}).get("files")
    if not isinstance(index_list, list) or len(index_list) != 6:
        raise ValueError(f"expected 6 bowtie2_index.files in reference validation receipt, got {index_list!r}")
    index_map: dict[str, str] = {}
    for entry in index_list:
        name = entry.get("name")
        if not isinstance(name, str) or not name:
            raise ValueError(f"invalid bowtie2 index entry name: {entry!r}")
        index_map[name] = require_valid_sha256_hex(entry.get("sha256"), f"bowtie2_index[{name}].sha256")
    if len(index_map) != 6:
        raise ValueError("duplicate bowtie2 index names in reference validation receipt")

    checks = val_receipt.get("checks", {})
    if not isinstance(checks, dict) or not checks or not all(v is True for v in checks.values()):
        raise ValueError(f"reference validation receipt checks not all True: {checks!r}")

    return {
        "validation_receipt_path": str(val_receipt_path),
        "validation_receipt_sha256": observed_val_sha,
        "validation_done_sha256": stream_sha256(val_done_path),
        "reference_fna_sha256": ref_fna_sha,
        "reference_fai_sha256": ref_fai_sha,
        "crosswalk_sha256": crosswalk_sha,
        "build_receipt_sha256": build_receipt_sha,
        "bowtie2_index_sha256": index_map,
    }


def validate_sample_mapping_metrics(map_metrics: Any, sample: str) -> dict[str, Any]:
    """Require producer's actual mapping_metrics keys and finite values; never default missing to zero."""
    if not isinstance(map_metrics, dict):
        raise ValueError(f"mapping_metrics missing or not a dict in {sample} receipt")
    missing_keys = sorted(REQUIRED_MAPPING_METRICS_KEYS - set(map_metrics))
    if missing_keys:
        raise ValueError(f"mapping_metrics in {sample} receipt missing required keys: {missing_keys}; missing is not zero")

    raw_total = require_strict_int(map_metrics["raw_total_records"], f"{sample}.raw_total_records")
    raw_primary = require_strict_int(map_metrics["raw_primary_mapped_records"], f"{sample}.raw_primary_mapped_records")
    filtered_hq = require_strict_int(
        map_metrics["filtered_high_quality_mapq20_reads"], f"{sample}.filtered_high_quality_mapq20_reads"
    )
    cov_refs = require_strict_int(
        map_metrics["coverage_references_verified"], f"{sample}.coverage_references_verified"
    )
    if not (raw_total >= raw_primary >= filtered_hq >= 0):
        raise ValueError(
            f"invalid record hierarchy in {sample} mapping_metrics: "
            f"raw_total={raw_total}, raw_primary={raw_primary}, filtered_hq={filtered_hq}"
        )
    if cov_refs <= 0:
        raise ValueError(f"invalid coverage_references_verified in {sample}: {cov_refs}")
    if map_metrics["filter"] != EXPECTED_MAPPING_FILTER_STRING:
        raise ValueError(
            f"unexpected mapping_metrics.filter in {sample}: "
            f"{map_metrics['filter']!r} != {EXPECTED_MAPPING_FILTER_STRING!r}"
        )
    frac = map_metrics["filtered_fraction_of_primary"]
    if raw_primary == 0:
        if frac is not None:
            raise ValueError(f"expected None filtered_fraction_of_primary when raw_primary==0 in {sample}, got {frac!r}")
        frac_value = math.nan
    else:
        if isinstance(frac, bool) or not isinstance(frac, (int, float)) or not math.isfinite(float(frac)):
            raise ValueError(f"non-finite filtered_fraction_of_primary in {sample}: {frac!r}")
        frac_value = float(frac)
        expected_frac = filtered_hq / raw_primary
        if not math.isclose(frac_value, expected_frac, rel_tol=1e-9, abs_tol=1e-9):
            raise ValueError(
                f"filtered_fraction_of_primary arithmetic mismatch in {sample}: {frac_value} != {expected_frac}"
            )
    return {
        "raw_total_records": raw_total,
        "raw_primary_mapped_records": raw_primary,
        "filtered_high_quality_mapq20_reads": filtered_hq,
        "filtered_fraction_of_primary": frac_value,
        "coverage_references_verified": cov_refs,
        "filter": map_metrics["filter"],
    }


def load_stage2_votu_context(stage2_path: Path, pairs: pd.DataFrame) -> dict[str, dict[str, Any]]:
    """Load existing vOTU metadata from stage2_target_viral_identity_table.tsv as context only."""
    if not stage2_path.exists() or stage2_path.stat().st_size == 0:
        raise FileNotFoundError(f"Stage 2 viral identity table missing/empty: {stage2_path}")
    frame = pd.read_csv(stage2_path, sep="\t", dtype=str, keep_default_na=False)
    require_columns(frame, STAGE2_REQUIRED_COLUMNS, "stage2_target_viral_identity_table")
    frame = frame.apply(lambda col: col.str.strip())
    by_target_id = {row["target_id"]: row for row in frame.to_dict(orient="records")}

    context_by_sha: dict[str, dict[str, Any]] = {}
    for target_sha, group in pairs.groupby("target_seq_sha256", sort=True):
        aliases: set[str] = set()
        for alias_list in group["target_aliases"]:
            if isinstance(alias_list, str):
                aliases.update(a.strip() for a in alias_list.split(";") if a.strip())
            else:
                aliases.update(str(a).strip() for a in alias_list if str(a).strip())
        sorted_aliases = sorted(aliases)
        canonical_ids = sorted(set(group["canonical_annotated_target_id"].astype(str)))
        if len(canonical_ids) != 1:
            raise ValueError(f"conflicting canonical_annotated_target_id for {target_sha}: {canonical_ids}")
        canonical_id = canonical_ids[0]
        if canonical_id not in by_target_id:
            raise ValueError(f"canonical target {canonical_id} absent from Stage 2 identity table")
        canon_row = by_target_id[canonical_id]
        if canon_row["target_seq_sha256"] != target_sha:
            raise ValueError(
                f"Stage 2 SHA mismatch for {canonical_id}: {canon_row['target_seq_sha256']} != {target_sha}"
            )
        alias_entries: list[str] = []
        alias_records: list[dict[str, Any]] = []
        for alias in sorted_aliases:
            if alias not in by_target_id:
                raise ValueError(f"target alias {alias} absent from Stage 2 identity table")
            arow = by_target_id[alias]
            if arow["target_seq_sha256"] != target_sha:
                raise ValueError(f"Stage 2 SHA mismatch for alias {alias}: {arow['target_seq_sha256']} != {target_sha}")
            is_rep = int(arow["is_votu_cluster_rep"])
            if is_rep not in (0, 1):
                raise ValueError(f"invalid is_votu_cluster_rep for {alias}: {arow['is_votu_cluster_rep']!r}")
            rep_id = arow["cluster_votu_rep_id"]
            if not rep_id:
                raise ValueError(f"empty cluster_votu_rep_id for {alias}")
            alias_entries.append(f"{alias}:rep={rep_id}:is_rep={is_rep}")
            alias_records.append({
                "target_id": alias,
                "cluster_votu_rep_id": rep_id,
                "is_votu_cluster_rep": is_rep,
                "viral_identity_status": arow["viral_identity_status"],
                "checkv_quality": arow["checkv_quality"],
                "checkv_viral_genes": int(arow["checkv_viral_genes"]),
                "checkv_host_genes": int(arow["checkv_host_genes"]),
                "virus_hallmarks": int(arow["virus_hallmarks"]),
            })
        canon_is_rep = int(canon_row["is_votu_cluster_rep"])
        canon_rep_id = canon_row["cluster_votu_rep_id"]
        context_by_sha[target_sha] = {
            "canonical_annotated_target_id": canonical_id,
            "stage2_cluster_votu_rep_id": canon_rep_id,
            "stage2_is_votu_cluster_rep": canon_is_rep,
            "stage2_alias_votu_metadata": ";".join(alias_entries),
            "stage2_alias_records": alias_records,
            "stage2_viral_identity_status": canon_row["viral_identity_status"],
            "stage2_checkv_quality": canon_row["checkv_quality"],
            "stage2_checkv_viral_genes": int(canon_row["checkv_viral_genes"]),
            "stage2_checkv_host_genes": int(canon_row["checkv_host_genes"]),
            "stage2_virus_hallmarks": int(canon_row["virus_hallmarks"]),
            "votu_metadata_role": (
                "context_only_not_demonstrated_cause_of_ambiguity; "
                "is_votu_cluster_rep=0_is_normal_nonrepresentative_member"
            ),
        }
    return context_by_sha


def build_pair_metadata_by_sha(pairs: pd.DataFrame) -> dict[str, dict[str, Any]]:
    meta: dict[str, dict[str, Any]] = {}
    for target_sha, group in pairs.groupby("target_seq_sha256", sort=True):
        aliases: set[str] = set()
        for values in group["target_aliases"]:
            if isinstance(values, str):
                aliases.update(v.strip() for v in values.split(";") if v.strip())
            else:
                aliases.update(str(v).strip() for v in values if str(v).strip())
        virus_samples: set[str] = set()
        for values in group["virus_samples"]:
            if isinstance(values, str):
                virus_samples.update(v.strip() for v in values.split(";") if v.strip())
            else:
                virus_samples.update(str(v).strip() for v in values if str(v).strip())
        roles = set(group["analysis_role"].astype(str))
        statuses = set(group["individual_virality_status"].astype(str))
        if len(roles) != 1 or len(statuses) != 1:
            raise ValueError(f"conflicting role/status for {target_sha}: roles={roles}, statuses={statuses}")
        meta[target_sha] = {
            "target_seq_sha256": target_sha,
            "entity_id": str(group.iloc[0]["target_id"]),
            "canonical_annotated_target_id": str(group.iloc[0]["canonical_annotated_target_id"]),
            "target_aliases": sorted(aliases),
            "host_mags": sorted(set(group["host_mag"].astype(str))),
            "host_phyla": sorted(set(group["host_phylum"].astype(str))),
            "host_samples": sorted(set(group["host_sample"].astype(str))),
            "virus_samples": sorted(virus_samples),
            "analysis_role": next(iter(roles)),
            "individual_virality_status": next(iter(statuses)),
        }
    if len(meta) != 7:
        raise ValueError(f"expected exactly 7 unique viral target SHAs, got {len(meta)}")
    return meta


def classify_recruitment_qc(
    raw_numreads: int,
    raw_meandepth: float,
    raw_breadth_pct: float,
    filtered_numreads: int,
    filtered_meandepth: float,
    filtered_breadth_pct: float,
    qualifies_primary: int,
    analysis_role: str,
    source_category_uniquely_viral: bool,
) -> tuple[str, str, str]:
    """Deterministically classify raw-vs-MAPQ20 recruitment without causal competitor overclaims."""
    if raw_numreads < 0 or filtered_numreads < 0 or filtered_numreads > raw_numreads:
        raise ValueError(f"invalid read counts: raw={raw_numreads}, filtered={filtered_numreads}")
    if raw_numreads == 0:
        return (
            "measured_zero_assigned_reads",
            "measured_zero_assigned_reads_to_reference_not_proved_biological_absence",
            (
                "Measured zero primary alignments assigned to this reference contig in this sample "
                "(not missing/untested, and not proved biological absence)."
            ),
        )
    if filtered_meandepth >= 10.0 and filtered_breadth_pct >= 50.0:
        if qualifies_primary == 1:
            return (
                "primary_mapq20_10x_50pct_eligible",
                "estimable_primary_mapq20_10x_50pct",
                (
                    "Meets primary MAPQ>=20 >=10x depth and >=50% breadth gates; supports competitive "
                    "read recruitment only, not active infection or replication."
                ),
            )
        return (
            "mapq20_10x_50pct_reached_excluded_by_role_or_ownership_gate",
            f"excluded_from_primary_instrain_by_role_{analysis_role}_or_ownership_gate",
            (
                "Reaches >=10x depth and >=50% breadth at MAPQ>=20, but excluded from primary inStrain "
                f"eligibility by analysis_role={analysis_role} or ownership gate."
            ),
        )
    if raw_meandepth >= 10.0 and raw_breadth_pct >= 50.0:
        return (
            "high_raw_support_10x_50pct_lost_at_mapq20",
            "high_raw_support_lost_at_mapq20_exact_reference_attribution_unresolved",
            (
                "High raw support (>=10x depth and >=50% breadth at MAPQ>=0) lost at MAPQ>=20; "
                "exact-reference attribution unresolved; compatible with competitive ambiguity or "
                "weak alignment discrimination (not sequence absence, unique exact-strain presence, "
                "or active infection)."
            ),
        )
    if filtered_numreads < raw_numreads:
        return (
            "low_total_read_support_with_mapq20_attrition",
            "low_total_sequence_support_with_mapq20_attrition_attribution_unresolved",
            (
                "Sub-threshold total raw read support (<10x depth or <50% breadth at MAPQ>=0) "
                "with additional read loss at MAPQ>=20; exact-reference attribution of lost reads "
                "unresolved; compatible with competitive ambiguity or weak alignment discrimination "
                "(not proved biological absence)."
            ),
        )
    return (
        "low_total_read_support_mapq20_retained",
        "low_total_sequence_support_below_10x_50pct_gate",
        (
            "Sub-threshold total read support (<10x depth or <50% breadth) with all assigned reads "
            "retained at MAPQ>=20; reflects low total sequence support, not proved biological absence."
        ),
    )


def run_samtools_region_coverage(
    bam_path: Path, reference_contig: str, min_mapq: int, expected_length: int
) -> tuple[dict[str, Any], str, list[str]]:
    cmd = [
        "samtools",
        "coverage",
        "--ff",
        str(FILTER_FLAGS_INT),
        "-q",
        str(min_mapq),
        "-r",
        reference_contig,
        str(bam_path),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"samtools coverage failed ({result.returncode}) for {bam_path} {reference_contig}: {result.stderr}"
        )
    frame = compute_coverage_metrics(result.stdout, {reference_contig: expected_length})
    if len(frame) != 1:
        raise ValueError(f"expected 1 coverage row for {reference_contig}, got {len(frame)}")
    row = frame.iloc[0].to_dict()
    data_lines = [line for line in result.stdout.splitlines() if line.strip() and not line.startswith("#")]
    if len(data_lines) != 1:
        raise ValueError(f"expected 1 non-header tabular line from samtools coverage, got {len(data_lines)}")
    return row, data_lines[0], cmd


def verify_and_extract_sample(
    sample: str,
    map_dir: Path,
    pair_meta_by_sha: dict[str, dict[str, Any]],
    stage2_context_by_sha: dict[str, dict[str, Any]],
    expected_pair_ledger_sha256: str,
    expected_boundary_config_sha256: str,
    expected_ref_validation: dict[str, Any],
) -> dict[str, Any]:
    sample_dir = map_dir / sample
    receipt_path = sample_dir / f"{sample}_mapping_receipt.json"
    done_path = sample_dir / "MAPPING_DONE"
    if not sample_dir.is_dir():
        raise FileNotFoundError(f"sample mapping directory missing: {sample_dir}; missing is not zero")
    if not receipt_path.is_file() or receipt_path.stat().st_size == 0:
        raise FileNotFoundError(f"mapping receipt missing/empty for {sample}: {receipt_path}")
    if not done_path.is_file() or done_path.stat().st_size == 0:
        raise FileNotFoundError(f"MAPPING_DONE sentinel missing/empty for {sample}: {done_path}")

    receipt_sha256 = stream_sha256(receipt_path)
    done_sha256 = stream_sha256(done_path)
    done_payload = json.loads(done_path.read_text())
    if done_payload.get("receipt_sha256") != receipt_sha256:
        raise ValueError(
            f"MAPPING_DONE receipt_sha256 mismatch for {sample}: "
            f"{done_payload.get('receipt_sha256')} != {receipt_sha256}"
        )

    receipt = json.loads(receipt_path.read_text())
    if receipt.get("status") != "completed_verified" or receipt.get("schema_version") != 3:
        raise ValueError(f"invalid mapping receipt status/schema for {sample}: {receipt_path}")
    if receipt.get("stage") != "stage4_map_sample_competitive":
        raise ValueError(f"unexpected receipt stage for {sample}: {receipt.get('stage')}")
    if receipt.get("sample") != sample:
        raise ValueError(f"sample mismatch in receipt {receipt_path}: {receipt.get('sample')} != {sample}")

    inputs = receipt.get("inputs", {})
    if inputs.get("pair_ledger_sha256") != expected_pair_ledger_sha256:
        raise ValueError(f"pair_ledger_sha256 mismatch in {sample} receipt")
    if inputs.get("boundary_config_sha256") != expected_boundary_config_sha256:
        raise ValueError(f"boundary_config_sha256 mismatch in {sample} receipt")

    ref_binding = receipt.get("reference_binding", {})
    for bind_key in ["validation_receipt_sha256", "reference_fna_sha256", "crosswalk_sha256", "build_receipt_sha256"]:
        if ref_binding.get(bind_key) != expected_ref_validation[bind_key]:
            raise ValueError(
                f"reference_binding.{bind_key} mismatch for {sample}: "
                f"{ref_binding.get(bind_key)} != {expected_ref_validation[bind_key]}"
            )
    if ref_binding.get("bowtie2_index_sha256") != expected_ref_validation["bowtie2_index_sha256"]:
        raise ValueError(f"reference_binding.bowtie2_index_sha256 mismatch for {sample}")

    validated_metrics = validate_sample_mapping_metrics(receipt.get("mapping_metrics"), sample)
    exec_env = receipt.get("execution_environment", {})
    slurm_job_id = str(exec_env.get("slurm_job_id", "")).strip()
    slurm_task_id = str(exec_env.get("slurm_array_task_id", "")).strip()
    if not slurm_job_id or not slurm_task_id:
        raise ValueError(f"missing slurm_job_id/slurm_array_task_id in {sample} execution_environment")

    outputs = receipt.get("outputs", {})
    verified_output_hashes: dict[str, str] = {}
    for path_key, sha_key in [
        ("raw_bam", "raw_bam_sha256"),
        ("filtered_bam", "filtered_bam_sha256"),
        ("coverage_tsv", "coverage_tsv_sha256"),
        ("focal_evaluation_tsv", "focal_evaluation_tsv_sha256"),
    ]:
        if path_key not in outputs or sha_key not in outputs:
            raise ValueError(f"missing output key {path_key}/{sha_key} in {sample} receipt")
        file_path = Path(outputs[path_key])
        if not file_path.is_file() or file_path.stat().st_size == 0:
            raise FileNotFoundError(f"receipted output missing/empty for {sample}: {file_path}")
        observed_sha = stream_sha256(file_path)
        if observed_sha != outputs[sha_key]:
            raise ValueError(
                f"receipted output hash mismatch for {sample} {path_key} ({file_path}): "
                f"{observed_sha} != {outputs[sha_key]}"
            )
        verified_output_hashes[sha_key] = observed_sha

    raw_bam = Path(outputs["raw_bam"])
    filtered_bam = Path(outputs["filtered_bam"])
    if "raw_bam_bytes" not in outputs:
        raise ValueError(f"missing raw_bam_bytes in {sample} receipt outputs")
    if raw_bam.stat().st_size != require_strict_int(outputs["raw_bam_bytes"], f"{sample}.raw_bam_bytes"):
        raise ValueError(f"raw_bam_bytes mismatch for {sample}: {raw_bam.stat().st_size} != {outputs['raw_bam_bytes']}")

    for bam in [raw_bam, filtered_bam]:
        bai = Path(str(bam) + ".bai")
        if not bai.is_file() or bai.stat().st_size == 0:
            raise FileNotFoundError(f"BAM index missing/empty for {sample}: {bai}")
        qc = subprocess.run(["samtools", "quickcheck", str(bam)], capture_output=True, text=True)
        if qc.returncode != 0:
            raise RuntimeError(f"samtools quickcheck failed for {bam}: {qc.stderr}")

    focal_evals = receipt.get("focal_evaluations", [])
    virus_evals = [row for row in focal_evals if row.get("entity_type") == "virus_target"]
    if len(virus_evals) != 7:
        raise ValueError(f"expected 7 virus_target rows in {sample} receipt, found {len(virus_evals)}")
    by_sha = {str(row["entity_sha256"]): row for row in virus_evals}
    if set(by_sha) != set(pair_meta_by_sha):
        raise ValueError(
            f"virus_target SHA set mismatch in {sample} receipt: "
            f"missing={sorted(set(pair_meta_by_sha) - set(by_sha))}, "
            f"extra={sorted(set(by_sha) - set(pair_meta_by_sha))}"
        )

    eval_tsv = pd.read_csv(outputs["focal_evaluation_tsv"], sep="\t")
    eval_virus_tsv = eval_tsv[eval_tsv["entity_type"] == "virus_target"]
    if len(eval_virus_tsv) != 7 or set(eval_virus_tsv["entity_sha256"].astype(str)) != set(pair_meta_by_sha):
        raise ValueError(f"focal_evaluation_tsv virus_target mismatch for {sample}")

    rows: list[dict[str, Any]] = []
    raw_tabular_lines: list[str] = []
    filtered_tabular_lines: list[str] = []
    commands_run: list[dict[str, Any]] = []

    for target_sha in sorted(pair_meta_by_sha):
        pmeta = pair_meta_by_sha[target_sha]
        vctx = stage2_context_by_sha[target_sha]
        feval = by_sha[target_sha]

        ref_contig = str(feval["reference_contig"])
        length_bp = int(feval["length_bp"])
        if length_bp <= 0:
            raise ValueError(f"invalid length_bp for {sample} {ref_contig}: {length_bp}")
        if str(feval["analysis_role"]) != pmeta["analysis_role"]:
            raise ValueError(
                f"analysis_role mismatch for {sample} {target_sha}: "
                f"{feval['analysis_role']} != {pmeta['analysis_role']}"
            )
        if str(feval["individual_status"]) != pmeta["individual_virality_status"]:
            raise ValueError(
                f"individual_status mismatch for {sample} {target_sha}: "
                f"{feval['individual_status']} != {pmeta['individual_virality_status']}"
            )

        raw_cov, raw_line, raw_cmd = run_samtools_region_coverage(raw_bam, ref_contig, 0, length_bp)
        filt_cov, filt_line, filt_cmd = run_samtools_region_coverage(filtered_bam, ref_contig, 20, length_bp)
        commands_run.append({"sample": sample, "target_sha256": target_sha, "mode": "raw_mapq0_ff3844", "command": raw_cmd})
        commands_run.append({"sample": sample, "target_sha256": target_sha, "mode": "filtered_mapq20_ff3844", "command": filt_cmd})
        raw_tabular_lines.append(f"{sample}\t{raw_line}")
        filtered_tabular_lines.append(f"{sample}\t{filt_line}")

        receipt_filt_reads = int(feval["filtered_reads_mapped"])
        receipt_filt_breadth = float(feval["horizontal_breadth_pct"])
        receipt_filt_depth = float(feval["mean_depth"])

        if int(filt_cov["numreads"]) != receipt_filt_reads:
            raise ValueError(
                f"filtered numreads mismatch for {sample} {ref_contig}: "
                f"samtools={filt_cov['numreads']} != receipt={receipt_filt_reads}"
            )
        recomputed_filt_breadth = 100.0 * int(filt_cov["covbases"]) / length_bp
        if not math.isclose(recomputed_filt_breadth, receipt_filt_breadth, rel_tol=1e-9, abs_tol=1e-9):
            raise ValueError(
                f"filtered breadth mismatch for {sample} {ref_contig}: "
                f"samtools={recomputed_filt_breadth} != receipt={receipt_filt_breadth}"
            )
        if not math.isclose(float(filt_cov["meandepth"]), receipt_filt_depth, rel_tol=1e-7, abs_tol=1e-7):
            raise ValueError(
                f"filtered meandepth mismatch for {sample} {ref_contig}: "
                f"samtools={filt_cov['meandepth']} != receipt={receipt_filt_depth}"
            )

        raw_reads = int(raw_cov["numreads"])
        raw_covbases = int(raw_cov["covbases"])
        raw_cov_pct = float(raw_cov["coverage"])
        raw_breadth_pct = 100.0 * raw_covbases / length_bp
        raw_depth = float(raw_cov["meandepth"])
        raw_baseq = float(raw_cov["meanbaseq"])
        raw_mapq = float(raw_cov["meanmapq"])

        filt_reads = int(filt_cov["numreads"])
        filt_covbases = int(filt_cov["covbases"])
        filt_cov_pct = float(filt_cov["coverage"])
        filt_breadth_pct = receipt_filt_breadth
        filt_depth = receipt_filt_depth
        filt_baseq = float(filt_cov["meanbaseq"])
        filt_mapq = float(filt_cov["meanmapq"])

        if filt_reads > raw_reads:
            raise ValueError(f"filtered reads ({filt_reads}) > raw reads ({raw_reads}) for {sample} {ref_contig}")
        if filt_covbases > raw_covbases:
            raise ValueError(
                f"filtered covbases ({filt_covbases}) > raw covbases ({raw_covbases}) for {sample} {ref_contig}"
            )
        if filt_depth - raw_depth > 1e-6:
            raise ValueError(f"filtered meandepth ({filt_depth}) > raw meandepth ({raw_depth}) for {sample} {ref_contig}")

        lost_reads = raw_reads - filt_reads
        retained_read_fraction = (filt_reads / raw_reads) if raw_reads > 0 else math.nan
        lost_read_fraction = (lost_reads / raw_reads) if raw_reads > 0 else math.nan
        lost_covbases = raw_covbases - filt_covbases
        breadth_loss_pct_points = raw_breadth_pct - filt_breadth_pct
        retained_breadth_fraction = (filt_covbases / raw_covbases) if raw_covbases > 0 else math.nan
        meandepth_loss = raw_depth - filt_depth
        retained_meandepth_fraction = (filt_depth / raw_depth) if raw_depth > 0 else math.nan

        raw_reaches_10x_50pct = int(raw_depth >= 10.0 and raw_breadth_pct >= 50.0)
        qualifies_primary = int(feval["qualifies_10x_50pct_instrain"])
        source_cat_unique = bool(feval["uniquely_attributable"])

        sample_roles: list[str] = []
        is_host_source = sample in pmeta["host_samples"]
        is_viral_source = sample in pmeta["virus_samples"]
        if is_host_source:
            sample_roles.append("host_source")
        if is_viral_source:
            sample_roles.append("viral_or_alias_source")
        if sample == "10Dr2":
            sample_roles.append("predetermined_unrelated_control")
        if not sample_roles:
            sample_roles.append("predetermined_coassembly_alias_or_copositive_sample")

        qc_class, non_est_reason, sci_qual = classify_recruitment_qc(
            raw_numreads=raw_reads,
            raw_meandepth=raw_depth,
            raw_breadth_pct=raw_breadth_pct,
            filtered_numreads=filt_reads,
            filtered_meandepth=filt_depth,
            filtered_breadth_pct=filt_breadth_pct,
            qualifies_primary=qualifies_primary,
            analysis_role=pmeta["analysis_role"],
            source_category_uniquely_viral=source_cat_unique,
        )

        rows.append({
            "sample": sample,
            "sample_roles": ";".join(sample_roles),
            "is_host_source_sample": is_host_source,
            "is_viral_or_alias_source_sample": is_viral_source,
            "entity_id": str(feval["entity_id"]),
            "canonical_annotated_target_id": pmeta["canonical_annotated_target_id"],
            "target_aliases": str(feval["target_aliases"]),
            "host_mags": str(feval["host_mags"]),
            "host_samples": ";".join(pmeta["host_samples"]),
            "virus_samples": ";".join(pmeta["virus_samples"]),
            "analysis_role": pmeta["analysis_role"],
            "individual_virality_status": pmeta["individual_virality_status"],
            "reference_contig": ref_contig,
            "reference_orientation_to_ledger": str(feval["reference_orientation_to_ledger"]),
            "raw_oriented_target_sha256": target_sha,
            "canonical_sha256": str(feval["canonical_sha256"]),
            "length_bp": length_bp,
            "startpos": int(raw_cov["startpos"]),
            "endpos": int(raw_cov["endpos"]),
            "source_category_uniquely_viral": source_cat_unique,
            "source_category_ambiguity_status": str(feval["ambiguity_status"]),
            "source_ownership_explanation": (
                "source_category_uniquely_viral reflects exact-sequence crosswalk category ownership only, "
                "never per-read mapping uniqueness across competing references"
            ),
            "stage2_cluster_votu_rep_id": vctx["stage2_cluster_votu_rep_id"],
            "stage2_is_votu_cluster_rep": vctx["stage2_is_votu_cluster_rep"],
            "stage2_alias_votu_metadata": vctx["stage2_alias_votu_metadata"],
            "stage2_checkv_quality": vctx["stage2_checkv_quality"],
            "stage2_checkv_viral_genes": vctx["stage2_checkv_viral_genes"],
            "stage2_checkv_host_genes": vctx["stage2_checkv_host_genes"],
            "stage2_virus_hallmarks": vctx["stage2_virus_hallmarks"],
            "votu_metadata_role": vctx["votu_metadata_role"],
            "raw_filter_flags_excluded": FILTER_FLAGS_INT,
            "raw_min_mapq": 0,
            "raw_numreads": raw_reads,
            "raw_covbases": raw_covbases,
            "raw_coverage_samtools_pct": raw_cov_pct,
            "raw_horizontal_breadth_pct": raw_breadth_pct,
            "raw_meandepth": raw_depth,
            "raw_meanbaseq": raw_baseq,
            "raw_meanmapq": raw_mapq,
            "raw_reaches_10x_50pct_qc_only": raw_reaches_10x_50pct,
            "filtered_filter_flags_excluded": FILTER_FLAGS_INT,
            "filtered_min_mapq": 20,
            "filtered_numreads": filt_reads,
            "unique_filtered_reads_mapped": int(feval["unique_filtered_reads_mapped"]),
            "ambiguous_filtered_reads_mapped": int(feval["ambiguous_filtered_reads_mapped"]),
            "filtered_covbases": filt_covbases,
            "filtered_coverage_samtools_pct": filt_cov_pct,
            "filtered_horizontal_breadth_pct": filt_breadth_pct,
            "filtered_meandepth": filt_depth,
            "filtered_meanbaseq": filt_baseq,
            "filtered_meanmapq": filt_mapq,
            "qualifies_primary_10x_50pct_instrain": qualifies_primary,
            "primary_instrain_action": str(feval["instrain_action"]),
            "lost_to_mapq20_reads": lost_reads,
            "retained_read_fraction": retained_read_fraction,
            "lost_to_mapq20_read_fraction": lost_read_fraction,
            "lost_to_mapq20_covbases": lost_covbases,
            "breadth_loss_pct_points": breadth_loss_pct_points,
            "retained_breadth_fraction": retained_breadth_fraction,
            "meandepth_loss": meandepth_loss,
            "retained_meandepth_fraction": retained_meandepth_fraction,
            "recruitment_qc_class": qc_class,
            "non_estimability_reason": non_est_reason,
            "scientific_qualification": sci_qual,
            "mapping_receipt_sha256": receipt_sha256,
            "raw_bam_sha256": verified_output_hashes["raw_bam_sha256"],
            "filtered_bam_sha256": verified_output_hashes["filtered_bam_sha256"],
        })

    sample_summary = {
        "sample": sample,
        "slurm_job_id": slurm_job_id,
        "slurm_array_task_id": slurm_task_id,
        "raw_total_records": validated_metrics["raw_total_records"],
        "raw_primary_mapped_records": validated_metrics["raw_primary_mapped_records"],
        "filtered_high_quality_mapq20_reads": validated_metrics["filtered_high_quality_mapq20_reads"],
        "filtered_fraction_of_primary": validated_metrics["filtered_fraction_of_primary"],
        "coverage_references_verified": validated_metrics["coverage_references_verified"],
        "primary_filter_spec": validated_metrics["filter"],
        "focal7_raw_reads_sum": int(sum(r["raw_numreads"] for r in rows)),
        "focal7_filtered_reads_sum": int(sum(r["filtered_numreads"] for r in rows)),
        "focal7_lost_to_mapq20_reads_sum": int(sum(r["lost_to_mapq20_reads"] for r in rows)),
        "focal7_n_targets_raw_reads_gt0": int(sum(r["raw_numreads"] > 0 for r in rows)),
        "focal7_n_targets_raw_reads_eq0": int(sum(r["raw_numreads"] == 0 for r in rows)),
        "focal7_n_targets_filtered_reads_gt0": int(sum(r["filtered_numreads"] > 0 for r in rows)),
        "focal7_n_targets_raw_10x_50pct_qc_only": int(sum(r["raw_reaches_10x_50pct_qc_only"] == 1 for r in rows)),
        "focal7_n_targets_primary_10x_50pct_eligible": int(sum(r["qualifies_primary_10x_50pct_instrain"] == 1 for r in rows)),
        "mapping_receipt_sha256": receipt_sha256,
        "done_sha256": done_sha256,
        "raw_bam_sha256": verified_output_hashes["raw_bam_sha256"],
        "filtered_bam_sha256": verified_output_hashes["filtered_bam_sha256"],
        "coverage_tsv_sha256": verified_output_hashes["coverage_tsv_sha256"],
        "focal_evaluation_tsv_sha256": verified_output_hashes["focal_evaluation_tsv_sha256"],
    }

    return {
        "sample": sample,
        "rows": rows,
        "raw_tabular_lines": raw_tabular_lines,
        "filtered_tabular_lines": filtered_tabular_lines,
        "sample_summary": sample_summary,
        "commands": commands_run,
    }


def summarize_by_target(frame: pd.DataFrame) -> pd.DataFrame:
    target_rows: list[dict[str, Any]] = []
    for target_sha, group in frame.groupby("raw_oriented_target_sha256", sort=True):
        first = group.iloc[0]
        source_rows = group[group["is_host_source_sample"] | group["is_viral_or_alias_source_sample"]]
        source_summary = ";".join(
            f"{r['sample']}:raw_reads={r['raw_numreads']},raw_depth={r['raw_meandepth']:.4f}x,"
            f"raw_breadth={r['raw_horizontal_breadth_pct']:.4f}%,raw_mapq={r['raw_meanmapq']:.2f},"
            f"filt_reads={r['filtered_numreads']},filt_depth={r['filtered_meandepth']:.4f}x,"
            f"filt_breadth={r['filtered_horizontal_breadth_pct']:.4f}%,class={r['recruitment_qc_class']}"
            for r in source_rows.to_dict(orient="records")
        )
        max_raw_idx = group["raw_meandepth"].astype(float).idxmax()
        max_raw_row = group.loc[max_raw_idx]
        max_filt_idx = group["filtered_meandepth"].astype(float).idxmax()
        max_filt_row = group.loc[max_filt_idx]
        target_rows.append({
            "entity_id": first["entity_id"],
            "canonical_annotated_target_id": first["canonical_annotated_target_id"],
            "target_aliases": first["target_aliases"],
            "raw_oriented_target_sha256": target_sha,
            "canonical_sha256": first["canonical_sha256"],
            "reference_contig": first["reference_contig"],
            "length_bp": int(first["length_bp"]),
            "host_mags": first["host_mags"],
            "host_samples": first["host_samples"],
            "virus_samples": first["virus_samples"],
            "analysis_role": first["analysis_role"],
            "individual_virality_status": first["individual_virality_status"],
            "source_category_uniquely_viral": bool(first["source_category_uniquely_viral"]),
            "stage2_cluster_votu_rep_id": first["stage2_cluster_votu_rep_id"],
            "stage2_is_votu_cluster_rep": int(first["stage2_is_votu_cluster_rep"]),
            "stage2_alias_votu_metadata": first["stage2_alias_votu_metadata"],
            "n_samples_evaluated": len(group),
            "n_samples_measured_zero_raw_reads": int((group["raw_numreads"] == 0).sum()),
            "n_samples_raw_reads_gt0": int((group["raw_numreads"] > 0).sum()),
            "n_samples_filtered_reads_gt0": int((group["filtered_numreads"] > 0).sum()),
            "n_samples_raw_reaches_10x_50pct_qc_only": int((group["raw_reaches_10x_50pct_qc_only"] == 1).sum()),
            "n_samples_high_raw_lost_at_mapq20": int(
                (group["recruitment_qc_class"] == "high_raw_support_10x_50pct_lost_at_mapq20").sum()
            ),
            "n_samples_primary_mapq20_10x_50pct_eligible": int(
                (group["qualifies_primary_10x_50pct_instrain"] == 1).sum()
            ),
            "total_raw_reads_across_14_samples": int(group["raw_numreads"].sum()),
            "total_filtered_reads_across_14_samples": int(group["filtered_numreads"].sum()),
            "total_lost_to_mapq20_reads_across_14_samples": int(group["lost_to_mapq20_reads"].sum()),
            "max_raw_meandepth_sample": str(max_raw_row["sample"]),
            "max_raw_meandepth": float(max_raw_row["raw_meandepth"]),
            "max_raw_horizontal_breadth_pct": float(max_raw_row["raw_horizontal_breadth_pct"]),
            "max_raw_meanmapq_at_peak_depth": float(max_raw_row["raw_meanmapq"]),
            "max_filtered_meandepth_sample": str(max_filt_row["sample"]),
            "max_filtered_meandepth": float(max_filt_row["filtered_meandepth"]),
            "max_filtered_horizontal_breadth_pct": float(max_filt_row["filtered_horizontal_breadth_pct"]),
            "source_sample_metrics": source_summary if source_summary else "no_direct_fastq_for_coassembly_source_sample",
        })
    return pd.DataFrame(target_rows)


def verify_authoritative_9sr2_row(frame: pd.DataFrame, verify_benchmark: bool = True) -> dict[str, Any]:
    sub = frame[
        (frame["sample"] == AUTHORITATIVE_9SR2_BENCHMARK["sample"])
        & (frame["reference_contig"] == AUTHORITATIVE_9SR2_BENCHMARK["reference_contig"])
    ]
    if sub.empty:
        if verify_benchmark:
            raise ValueError("authoritative 9Sr2 target row absent from 98-row coverage table")
        return {"verified": False, "reason": "9Sr2 not in synthetic fixture scope"}
    if len(sub) != 1:
        raise ValueError(f"expected 1 row for 9Sr2 target, got {len(sub)}")
    row = sub.iloc[0]
    if verify_benchmark:
        for int_key in ["length_bp", "raw_numreads", "raw_covbases", "filtered_numreads", "filtered_covbases"]:
            if int(row[int_key]) != int(AUTHORITATIVE_9SR2_BENCHMARK[int_key]):
                raise ValueError(
                    f"9Sr2 benchmark mismatch for {int_key}: {row[int_key]} != {AUTHORITATIVE_9SR2_BENCHMARK[int_key]}"
                )
        for float_key in [
            "raw_coverage_samtools_pct",
            "raw_meandepth",
            "raw_meanbaseq",
            "raw_meanmapq",
            "filtered_horizontal_breadth_pct",
            "filtered_meandepth",
        ]:
            if not math.isclose(float(row[float_key]), float(AUTHORITATIVE_9SR2_BENCHMARK[float_key]), rel_tol=1e-6, abs_tol=1e-6):
                raise ValueError(
                    f"9Sr2 benchmark mismatch for {float_key}: {row[float_key]} != {AUTHORITATIVE_9SR2_BENCHMARK[float_key]}"
                )
        if row["recruitment_qc_class"] != "high_raw_support_10x_50pct_lost_at_mapq20":
            raise ValueError(f"unexpected 9Sr2 recruitment_qc_class: {row['recruitment_qc_class']}")
        if int(row["qualifies_primary_10x_50pct_instrain"]) != 0:
            raise ValueError("9Sr2 raw coverage must never override primary MAPQ>=20 inStrain eligibility")
        if int(row["stage2_is_votu_cluster_rep"]) != 0 or str(row["stage2_cluster_votu_rep_id"]) != "55Sr3::55Sr3_129891":
            raise ValueError("9Sr2 Stage 2 vOTU context metadata mismatch")
    return {
        "verified": True,
        "sample": "9Sr2",
        "entity_id": str(row["entity_id"]),
        "reference_contig": str(row["reference_contig"]),
        "length_bp": int(row["length_bp"]),
        "raw_numreads": int(row["raw_numreads"]),
        "raw_covbases": int(row["raw_covbases"]),
        "raw_coverage_samtools_pct": float(row["raw_coverage_samtools_pct"]),
        "raw_horizontal_breadth_pct": float(row["raw_horizontal_breadth_pct"]),
        "raw_meandepth": float(row["raw_meandepth"]),
        "raw_meanbaseq": float(row["raw_meanbaseq"]),
        "raw_meanmapq": float(row["raw_meanmapq"]),
        "filtered_numreads": int(row["filtered_numreads"]),
        "filtered_covbases": int(row["filtered_covbases"]),
        "filtered_horizontal_breadth_pct": float(row["filtered_horizontal_breadth_pct"]),
        "filtered_meandepth": float(row["filtered_meandepth"]),
        "lost_to_mapq20_reads": int(row["lost_to_mapq20_reads"]),
        "retained_read_fraction": float(row["retained_read_fraction"]),
        "recruitment_qc_class": str(row["recruitment_qc_class"]),
        "stage2_cluster_votu_rep_id": str(row["stage2_cluster_votu_rep_id"]),
        "stage2_is_votu_cluster_rep": int(row["stage2_is_votu_cluster_rep"]),
    }


def run_raw_vs_mapq20_qc(
    map_dir: Path = MAP_DIR,
    ref_dir: Path = REF_DIR,
    pair_ledger: Path = RES_DIR / "stage3_tier1_pairs_individual_evidence_audit.tsv",
    boundary_config: Path = BRANCH_ROOT / "inputs/stage4_boundary_sensitive_intervals.json",
    stage2_identity_table: Path = RES_DIR / "stage2_target_viral_identity_table.tsv",
    submission_receipt: Path = BRANCH_ROOT / "provenance/stage4_competitive_pipeline_submission.json",
    out_dir: Path = OUT_DIR,
    threads: int = 14,
    verify_9sr2_benchmark: bool = True,
) -> dict[str, Any]:
    start = time.time()
    if threads < 1:
        raise ValueError("threads must be positive")
    if shutil.which("samtools") is None:
        raise RuntimeError("required tool unavailable: samtools")
    samtools_version = probe_tool_version(["samtools", "--version"])

    expected_samples, submission_payload = load_expected_samples(submission_receipt)
    existing_receipt_dirs = sorted(
        path.parent.name for path in map_dir.glob("*/*_mapping_receipt.json")
    )
    missing_samples = [s for s in expected_samples if not (map_dir / s / "MAPPING_DONE").is_file()]
    if missing_samples:
        raise FileNotFoundError(
            f"cannot complete 14-sample raw-vs-MAPQ20 QC: missing MAPPING_DONE for {missing_samples}; missing is not zero"
        )
    extra_samples = sorted(set(existing_receipt_dirs) - set(expected_samples))
    if extra_samples:
        raise ValueError(f"unexpected sample mapping receipts outside 14-sample scope: {extra_samples}")

    pairs, reviewed_pairs, boundary_pairs = load_analysis_pairs(pair_ledger, boundary_config)
    pair_meta_by_sha = build_pair_metadata_by_sha(pairs)
    stage2_context_by_sha = load_stage2_votu_context(stage2_identity_table, pairs)

    pair_ledger_sha = stream_sha256(pair_ledger)
    boundary_config_sha = stream_sha256(boundary_config)
    stage2_table_sha = stream_sha256(stage2_identity_table)
    submission_sha = stream_sha256(submission_receipt)

    expected_ref_validation = verify_reference_validation_contract(
        ref_dir=ref_dir,
        expected_pair_ledger_sha256=pair_ledger_sha,
        expected_boundary_config_sha256=boundary_config_sha,
    )

    out_dir.mkdir(parents=True, exist_ok=True)
    done_sentinel = out_dir / "STAGE4_RAW_VS_MAPQ20_DONE"
    done_sentinel.unlink(missing_ok=True)

    worker_count = min(threads, len(expected_samples))
    with ThreadPoolExecutor(max_workers=worker_count) as pool:
        futures = [
            pool.submit(
                verify_and_extract_sample,
                sample,
                map_dir,
                pair_meta_by_sha,
                stage2_context_by_sha,
                pair_ledger_sha,
                boundary_config_sha,
                expected_ref_validation,
            )
            for sample in expected_samples
        ]
        sample_results = [future.result() for future in futures]

    all_rows: list[dict[str, Any]] = []
    raw_tab_lines = ["sample\t#" + "\t".join(COVERAGE_COLUMNS)]
    filt_tab_lines = ["sample\t#" + "\t".join(COVERAGE_COLUMNS)]
    sample_summaries: list[dict[str, Any]] = []
    all_commands: list[dict[str, Any]] = []

    for res in sample_results:
        all_rows.extend(res["rows"])
        raw_tab_lines.extend(res["raw_tabular_lines"])
        filt_tab_lines.extend(res["filtered_tabular_lines"])
        sample_summaries.append(res["sample_summary"])
        all_commands.extend(res["commands"])

    if len(all_rows) != 98:
        raise ValueError(f"expected exactly 98 sample x target rows (14 x 7), got {len(all_rows)}")

    frame = pd.DataFrame(all_rows)
    target_summary_df = summarize_by_target(frame)
    sample_summary_df = pd.DataFrame(sample_summaries)

    auth_9sr2 = verify_authoritative_9sr2_row(frame, verify_benchmark=verify_9sr2_benchmark)

    rows_tsv = out_dir / "stage4_raw_vs_mapq20_focal_coverage_98rows.tsv"
    raw_tab_tsv = out_dir / "stage4_raw_vs_mapq20_samtools_raw_tabular_98rows.tsv"
    filt_tab_tsv = out_dir / "stage4_raw_vs_mapq20_samtools_filtered_tabular_98rows.tsv"
    target_tsv = out_dir / "stage4_raw_vs_mapq20_target_summary.tsv"
    sample_tsv = out_dir / "stage4_raw_vs_mapq20_sample_summary.tsv"
    summary_json = out_dir / "stage4_raw_vs_mapq20_summary.json"

    frame.to_csv(rows_tsv, sep="\t", index=False)
    atomic_text(raw_tab_tsv, "\n".join(raw_tab_lines) + "\n")
    atomic_text(filt_tab_tsv, "\n".join(filt_tab_lines) + "\n")
    target_summary_df.to_csv(target_tsv, sep="\t", index=False)
    sample_summary_df.to_csv(sample_tsv, sep="\t", index=False)

    class_counts_all = frame["recruitment_qc_class"].value_counts().to_dict()
    source_mask = frame["is_host_source_sample"] | frame["is_viral_or_alias_source_sample"]
    class_counts_source = frame.loc[source_mask, "recruitment_qc_class"].value_counts().to_dict()

    source_evaluations = frame.loc[
        source_mask,
        [
            "sample",
            "sample_roles",
            "entity_id",
            "reference_contig",
            "analysis_role",
            "length_bp",
            "raw_numreads",
            "raw_covbases",
            "raw_horizontal_breadth_pct",
            "raw_meandepth",
            "raw_meanmapq",
            "raw_reaches_10x_50pct_qc_only",
            "filtered_numreads",
            "filtered_covbases",
            "filtered_horizontal_breadth_pct",
            "filtered_meandepth",
            "filtered_meanmapq",
            "qualifies_primary_10x_50pct_instrain",
            "lost_to_mapq20_reads",
            "retained_read_fraction",
            "recruitment_qc_class",
            "non_estimability_reason",
            "stage2_cluster_votu_rep_id",
            "stage2_is_votu_cluster_rep",
        ],
    ].to_dict(orient="records")

    summary_payload = {
        "status": "completed_verified",
        "schema_version": 1,
        "stage": "stage4_raw_vs_mapq20_focal_coverage",
        "timestamp_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "elapsed_seconds": time.time() - start,
        "execution_environment": {
            "hostname": platform.node(),
            "kernel": platform.release(),
            "slurm_job_id": os.environ.get("SLURM_JOB_ID", "local"),
            "python_version": sys.version.split()[0],
            "tmpdir": os.environ.get("TMPDIR", ""),
            "tool_versions": {
                "samtools": samtools_version,
            },
        },
        "scope_counts": {
            "n_samples_expected_and_verified": len(expected_samples),
            "samples": expected_samples,
            "n_unique_viral_targets": len(pair_meta_by_sha),
            "n_primary_viral_targets": int((target_summary_df["analysis_role"] == "primary").sum()),
            "n_sensitivity_unresolved_virality_targets": int(
                (target_summary_df["analysis_role"] == "sensitivity_unresolved_virality").sum()
            ),
            "n_boundary_sensitive_provisional_targets": int(
                (target_summary_df["analysis_role"] == "boundary_sensitive_provisional").sum()
            ),
            "n_reviewed_host_target_sha_pairs": int(len(reviewed_pairs)),
            "n_boundary_host_target_sha_pairs": int(len(boundary_pairs)),
            "n_total_logical_pairs": int(len(pairs)),
            "n_sample_target_rows": int(len(frame)),
            "n_rows_measured_zero_raw_reads": int((frame["raw_numreads"] == 0).sum()),
            "n_rows_raw_reads_gt0": int((frame["raw_numreads"] > 0).sum()),
            "n_rows_filtered_reads_gt0": int((frame["filtered_numreads"] > 0).sum()),
            "n_rows_raw_reaches_10x_50pct_qc_only": int((frame["raw_reaches_10x_50pct_qc_only"] == 1).sum()),
            "n_rows_high_raw_lost_at_mapq20": int(
                (frame["recruitment_qc_class"] == "high_raw_support_10x_50pct_lost_at_mapq20").sum()
            ),
            "n_rows_primary_mapq20_10x_50pct_eligible": int(
                (frame["qualifies_primary_10x_50pct_instrain"] == 1).sum()
            ),
        },
        "filter_specification": {
            "exclusion_flags_int": FILTER_FLAGS_INT,
            "exclusion_flags_hex": hex(FILTER_FLAGS_INT),
            "exclusion_flags_names": FILTER_FLAGS_NAMES,
            "raw_coverage_command_template": "samtools coverage --ff 3844 -q 0 -r <reference_contig> <raw_bam>",
            "filtered_coverage_command_template": "samtools coverage --ff 3844 -q 20 -r <reference_contig> <filtered_bam>",
            "isolated_comparison_variable": "minimum MAPQ threshold only (-q 0 vs -q 20) with identical flag mask 3844",
            "primary_eligibility_gate": (
                "mean_depth >= 10.0 and horizontal_breadth_pct >= 50.0 and "
                "source_category_uniquely_viral and analysis_role == 'primary' at MAPQ >= 20"
            ),
            "raw_10x_50pct_gate_role": "diagnostic_qc_only_never_overrides_primary_eligibility",
        },
        "authoritative_9Sr2_verification": auth_9sr2,
        "recruitment_qc_class_counts_all98": {k: int(v) for k, v in sorted(class_counts_all.items())},
        "recruitment_qc_class_counts_source_samples": {k: int(v) for k, v in sorted(class_counts_source.items())},
        "source_sample_focal_evaluations": source_evaluations,
        "votu_metadata_context_summary": {
            "role": "context_only_not_demonstrated_cause_of_ambiguity",
            "semantics_note": (
                "is_votu_cluster_rep is a property of the target row in stage2_target_viral_identity_table.tsv. "
                "is_votu_cluster_rep=0 with a populated cluster_votu_rep_id is normal for a non-representative "
                "member (e.g., 9Sr2::9Sr2_284237 points to 55Sr3::55Sr3_129891), and exact-sequence aliases can "
                "differ in whether a given identifier was chosen as representative (e.g., 48Sr2::48Sr2_190170 "
                "has is_votu_cluster_rep=0 pointing to 48PRr2::48PRr2_1185466 while its exact sequence alias "
                "coasm_48_S::k127_1210458 has is_votu_cluster_rep=1)."
            ),
            "targets": stage2_context_by_sha,
        },
        "scientific_interpretation_and_limits": {
            "high_raw_lost_at_mapq20_vs_low_total_support": (
                "High raw support lost at MAPQ>=20 establishes mapping uncertainty and leaves exact-reference "
                "attribution unresolved; it is compatible with competitive ambiguity or weak alignment discrimination, "
                "and does not by itself identify which competing reference caused the loss or prove that every "
                "low-MAPQ read has a related competing alignment (e.g., 9Sr2::9Sr2_284237 in 9Sr2: 12,065 raw reads, "
                "30.919x depth, 99.9771% breadth, mean MAPQ 1.51 -> 181 MAPQ>=20 reads, 0.469854x depth, 3.734637% breadth)."
            ),
            "raw_assigned_reads_limit": (
                "Raw (MAPQ>=0) assigned primary alignments are NOT evidence of unique exact-strain presence, "
                "active infection, replication, or which competing reference genome is the true biological source."
            ),
            "measured_zero_limit": (
                "Measured zero assigned reads to a reference contig (measured_zero_assigned_reads) is NOT "
                "proved biological absence of the virus in the sample or habitat."
            ),
            "source_category_vs_per_read_uniqueness": (
                "source_category_uniquely_viral indicates only that exact sequence duplicates in the reference "
                "crosswalk belong to viral categories; it does NOT imply per-read mapping uniqueness."
            ),
            "primary_thresholds_preserved": (
                "Primary >=10x depth, >=50% breadth, MAPQ>=20, and analysis_role=='primary' gates are preserved "
                "unchanged; raw_reaches_10x_50pct_qc_only is diagnostic QC only."
            ),
        },
        "inputs": {
            "submission_receipt": str(submission_receipt),
            "submission_receipt_sha256": submission_sha,
            "job_validate_ref": str(submission_payload.get("job_validate_ref", "")),
            "job_map_array": str(submission_payload.get("job_map_array", "")),
            "pair_ledger": str(pair_ledger),
            "pair_ledger_sha256": pair_ledger_sha,
            "boundary_config": str(boundary_config),
            "boundary_config_sha256": boundary_config_sha,
            "stage2_identity_table": str(stage2_identity_table),
            "stage2_identity_table_sha256": stage2_table_sha,
            "reference_validation": expected_ref_validation,
            "script_path": str(Path(__file__).resolve()),
            "script_sha256": stream_sha256(Path(__file__).resolve()),
            "common_script_sha256": stream_sha256(Path(__file__).with_name("stage4_competitive_common.py")),
            "mapping_script_sha256": stream_sha256(Path(__file__).with_name("stage4_map_sample_competitive.py")),
            "sample_mapping_receipts": {
                s["sample"]: {
                    "mapping_receipt_sha256": s["mapping_receipt_sha256"],
                    "done_sha256": s["done_sha256"],
                    "raw_bam_sha256": s["raw_bam_sha256"],
                    "filtered_bam_sha256": s["filtered_bam_sha256"],
                    "coverage_tsv_sha256": s["coverage_tsv_sha256"],
                    "focal_evaluation_tsv_sha256": s["focal_evaluation_tsv_sha256"],
                }
                for s in sample_summaries
            },
        },
        "commands_executed_count": len(all_commands),
        "outputs": {
            "focal_coverage_98rows_tsv": str(rows_tsv),
            "focal_coverage_98rows_tsv_bytes": rows_tsv.stat().st_size,
            "focal_coverage_98rows_tsv_sha256": stream_sha256(rows_tsv),
            "samtools_raw_tabular_98rows_tsv": str(raw_tab_tsv),
            "samtools_raw_tabular_98rows_tsv_bytes": raw_tab_tsv.stat().st_size,
            "samtools_raw_tabular_98rows_tsv_sha256": stream_sha256(raw_tab_tsv),
            "samtools_filtered_tabular_98rows_tsv": str(filt_tab_tsv),
            "samtools_filtered_tabular_98rows_tsv_bytes": filt_tab_tsv.stat().st_size,
            "samtools_filtered_tabular_98rows_tsv_sha256": stream_sha256(filt_tab_tsv),
            "target_summary_tsv": str(target_tsv),
            "target_summary_tsv_bytes": target_tsv.stat().st_size,
            "target_summary_tsv_sha256": stream_sha256(target_tsv),
            "sample_summary_tsv": str(sample_tsv),
            "sample_summary_tsv_bytes": sample_tsv.stat().st_size,
            "sample_summary_tsv_sha256": stream_sha256(sample_tsv),
        },
    }
    atomic_json(summary_json, summary_payload)
    summary_sha = stream_sha256(summary_json)
    atomic_text(done_sentinel, json.dumps({"summary_sha256": summary_sha}) + "\n")
    print(f"[Completed verified] 98-row raw-vs-MAPQ20 focal coverage QC: {summary_json}", flush=True)
    return summary_payload


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--map-dir", type=Path, default=MAP_DIR)
    parser.add_argument("--ref-dir", type=Path, default=REF_DIR)
    parser.add_argument("--pair-ledger", type=Path, default=RES_DIR / "stage3_tier1_pairs_individual_evidence_audit.tsv")
    parser.add_argument("--boundary-config", type=Path, default=BRANCH_ROOT / "inputs/stage4_boundary_sensitive_intervals.json")
    parser.add_argument("--stage2-identity-table", type=Path, default=RES_DIR / "stage2_target_viral_identity_table.tsv")
    parser.add_argument("--submission-receipt", type=Path, default=BRANCH_ROOT / "provenance/stage4_competitive_pipeline_submission.json")
    parser.add_argument("--out-dir", type=Path, default=OUT_DIR)
    parser.add_argument("--threads", type=int, default=14)
    args = parser.parse_args()
    run_raw_vs_mapq20_qc(
        map_dir=args.map_dir,
        ref_dir=args.ref_dir,
        pair_ledger=args.pair_ledger,
        boundary_config=args.boundary_config,
        stage2_identity_table=args.stage2_identity_table,
        submission_receipt=args.submission_receipt,
        out_dir=args.out_dir,
        threads=args.threads,
        verify_9sr2_benchmark=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
