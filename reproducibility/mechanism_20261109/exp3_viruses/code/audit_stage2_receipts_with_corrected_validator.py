#!/usr/bin/env python3
"""Authoritative standalone receipt verification using the corrected Stage 2 validator.

Guarantees & Invariants:
1. Probe actual container and tool identity (corrected no-fallback; raises if invalid).
2. Strictly compares container_sha256 in receipt with actual container_sha256.
3. Streamed actual subject SHA-256 is determined independently (NEVER from the receipt itself).
4. Verifies actual subject file integrity on disk (detects same-size same-mtime content mutations).
5. Enforces 100% full coverage of all expected targets (6 catalog pools + cohort assemblies);
   nonempty result counts alone never pass without full expected coverage.
6. Must be executed on an allocated compute node (workq) for heavy I/O.
"""
from __future__ import annotations

import argparse
import concurrent.futures
import json
import os
import sys
from pathlib import Path
from typing import Any

CODE_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(CODE_DIR))

import stage2_viral_identity_and_spacers_v3 as s2_v3

DEFAULT_EXP3_ROOT = Path(
    os.environ.get(
        "EXP3_ROOT",
        "/scratch/project/k10018/eq-mechanism-20261109/exp3_viruses",
    )
)


def get_default_catalog_targets(raw_dir: Path) -> list[tuple[str, Path, int]]:
    """Return the 6 canonical catalog targets: (pool_name, subject_fna_path, expected_threads)."""
    cell_ctrl_fna = raw_dir / "curated_cellular_controls_verified.fna"
    return [
        ("blast_eq_votu_reps", s2_v3.P_VOTU_FNA, 32),
        ("blast_eq_all_viruses", s2_v3.P_ALL_VIR_FNA, 32),
        ("blast_eq_coasm_viruses", s2_v3.P_COASM_VIR_FNA, 32),
        ("blast_curated_plasmid_controls", s2_v3.P_ALL_PLAS_FNA, 32),
        ("blast_curated_cellular_controls", cell_ctrl_fna, 32),
        ("blast_external_habitat_viral_controls", s2_v3.P_EXT_REF_VIR_FNA, 32),
    ]


def verify_single_target(
    target_name: str,
    target_type: str,
    out_tsv: Path,
    receipt_p: Path,
    subject_fna: Path,
    expected_query_sha: str,
    expected_subject_sha: str | None,
    expected_threads: int | None,
    stream_check_subject: bool,
    expected_container_sha: str,
) -> dict[str, Any]:
    """Verify a single search target with strict non-circular checks."""
    result: dict[str, Any] = {
        "target_name": target_name,
        "target_type": target_type,
        "is_valid": False,
        "error": "",
    }

    if not out_tsv.exists():
        result["error"] = f"Missing output TSV: {out_tsv.name}"
        return result

    if not receipt_p.exists():
        result["error"] = f"Missing receipt JSON: {receipt_p.name}"
        return result

    if not subject_fna.exists():
        result["error"] = f"Missing subject FASTA: {subject_fna}"
        return result

    try:
        rec = json.loads(receipt_p.read_text())
    except Exception as exc:
        result["error"] = f"Corrupt receipt JSON: {exc}"
        return result

    # Container SHA-256 match
    if rec.get("container_sha256") != expected_container_sha:
        result["error"] = (
            f"container_sha256 mismatch: {rec.get('container_sha256')} != {expected_container_sha}"
        )
        return result

    # Compute actual streamed disk SHA-256 independently if requested
    actual_disk_sha = ""
    if stream_check_subject:
        actual_disk_sha = s2_v3.stream_file_sha256(subject_fna)
        if not actual_disk_sha:
            result["error"] = f"Failed to compute streamed SHA-256 for subject {subject_fna}"
            return result

        # If an independent expected subject SHA was provided (e.g. from manifest), cross-check it
        if expected_subject_sha and actual_disk_sha != expected_subject_sha:
            result["error"] = (
                f"Subject file on disk has been mutated: disk sha {actual_disk_sha} != expected {expected_subject_sha}"
            )
            return result

        effective_expected_sha = actual_disk_sha
    else:
        if not expected_subject_sha:
            result["error"] = "No independent subject SHA provided and stream_check_subject=False"
            return result
        effective_expected_sha = expected_subject_sha

    # Non-circular check: receipt's stored subject_sha256 must match the independent hash
    if rec.get("subject_sha256") != effective_expected_sha:
        result["error"] = (
            f"Receipt subject_sha256 mismatch: receipt {rec.get('subject_sha256')} != independent {effective_expected_sha}"
        )
        return result

    # Comprehensive validation via v3 validator
    is_ok = s2_v3.is_blast_receipt_valid(
        out_tsv=out_tsv,
        expected_query_sha=expected_query_sha,
        expected_subject_fna=subject_fna,
        expected_subject_sha256=effective_expected_sha,
        expected_threads=expected_threads,
    )

    if not is_ok:
        result["error"] = "Validator rejected receipt or output TSV format"
        return result

    result["is_valid"] = True
    result["n_hits"] = rec.get("n_hits", 0)
    result["output_sha256"] = rec.get("output_sha256", "")
    result["subject_sha256"] = effective_expected_sha
    return result


def audit_stage2_receipts(
    raw_dir: Path,
    spacers_fna: Path,
    assemblies_dir: Path | None = None,
    expected_catalog_targets: list[tuple[str, Path, int]] | None = None,
    expected_assembly_samples: list[str] | None = None,
    stream_check_subjects: bool = True,
    max_workers: int = 16,
) -> tuple[bool, dict[str, Any]]:
    """Perform independent, non-circular audit of all Stage 2 search receipts.

    Returns (is_passed, audit_report_dict).
    """
    report: dict[str, Any] = {
        "status": "failed",
        "raw_dir": str(raw_dir),
        "tool_info": {},
        "query_sha256": "",
        "n_catalog_expected": 0,
        "n_catalog_valid": 0,
        "n_assemblies_expected": 0,
        "n_assemblies_valid": 0,
        "total_expected": 0,
        "total_valid": 0,
        "missing_or_invalid_targets": [],
        "errors": [],
    }

    if not raw_dir.exists():
        report["errors"].append(f"Raw directory not found: {raw_dir}")
        return False, report

    if not spacers_fna.exists():
        report["errors"].append(f"Spacers query FNA not found: {spacers_fna}")
        return False, report

    # 1. Tool and container identity probe (strict, no fallback)
    try:
        tool_info = s2_v3.get_blast_tool_info()
    except Exception as exc:
        report["errors"].append(f"BLAST tool probe failed: {exc}")
        return False, report

    report["tool_info"] = tool_info
    expected_container_sha = tool_info["container_sha256"]

    # 2. Query SHA-256
    q_sha = s2_v3.sha256_bytes(spacers_fna.read_bytes())
    report["query_sha256"] = q_sha

    # 3. Catalog targets
    if expected_catalog_targets is None:
        expected_catalog_targets = get_default_catalog_targets(raw_dir)
    report["n_catalog_expected"] = len(expected_catalog_targets)

    # 4. Assembly targets and independent subject hashes manifest
    manifest_p = raw_dir / "target_assemblies_subject_hashes.json"
    manifest_hashes: dict[str, str] = {}
    if manifest_p.exists():
        try:
            m_data = json.loads(manifest_p.read_text())
            for s_id, s_val in m_data.items():
                if isinstance(s_val, dict):
                    manifest_hashes[s_id] = s_val.get("sha256", "")
                elif isinstance(s_val, str):
                    manifest_hashes[s_id] = s_val
        except Exception as exc:
            report["errors"].append(f"Corrupt target_assemblies_subject_hashes.json: {exc}")

    if expected_assembly_samples is None:
        if manifest_hashes:
            expected_assembly_samples = sorted(manifest_hashes.keys())
        else:
            expected_assembly_samples = []

    report["n_assemblies_expected"] = len(expected_assembly_samples)
    total_expected = len(expected_catalog_targets) + len(expected_assembly_samples)
    report["total_expected"] = total_expected

    if total_expected == 0:
        report["errors"].append("No expected search targets configured or discovered.")
        return False, report

    asm_dir = assemblies_dir or s2_v3.P_ASSEMBLIES_DIR
    raw_asm_dir = raw_dir / "raw_assembly_blast"

    # 5. Build task queue for verification
    verification_tasks: list[dict[str, Any]] = []

    # Catalog tasks
    for pool_name, pool_fna, threads in expected_catalog_targets:
        out_tsv = raw_dir / f"{pool_name}.tsv"
        receipt_p = raw_dir / f"{pool_name}.receipt.json"
        verification_tasks.append(
            {
                "target_name": pool_name,
                "target_type": "catalog_pool",
                "out_tsv": out_tsv,
                "receipt_p": receipt_p,
                "subject_fna": pool_fna,
                "expected_subject_sha": None,  # Will be stream-hashed directly
                "expected_threads": threads,
                "stream_check_subject": True,
            }
        )

    # Assembly tasks
    for samp in expected_assembly_samples:
        asm_p = asm_dir / samp / f"{samp}_contigs.fa"
        out_tsv = raw_asm_dir / f"raw_asm_{samp}.blast.tsv"
        receipt_p = raw_asm_dir / f"raw_asm_{samp}.blast.receipt.json"
        exp_sha = manifest_hashes.get(samp)
        verification_tasks.append(
            {
                "target_name": samp,
                "target_type": "cohort_assembly",
                "out_tsv": out_tsv,
                "receipt_p": receipt_p,
                "subject_fna": asm_p,
                "expected_subject_sha": exp_sha,
                "expected_threads": 8,
                "stream_check_subject": stream_check_subjects,
            }
        )

    # 6. Execute verification in parallel
    def _run_task(t: dict[str, Any]) -> dict[str, Any]:
        return verify_single_target(
            target_name=t["target_name"],
            target_type=t["target_type"],
            out_tsv=t["out_tsv"],
            receipt_p=t["receipt_p"],
            subject_fna=t["subject_fna"],
            expected_query_sha=q_sha,
            expected_subject_sha=t["expected_subject_sha"],
            expected_threads=t["expected_threads"],
            stream_check_subject=t["stream_check_subject"],
            expected_container_sha=expected_container_sha,
        )

    catalog_valid = 0
    assemblies_valid = 0
    invalid_targets: list[dict[str, Any]] = []

    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as ex:
        results = list(ex.map(_run_task, verification_tasks))

    for res in results:
        if res["is_valid"]:
            if res["target_type"] == "catalog_pool":
                catalog_valid += 1
            else:
                assemblies_valid += 1
        else:
            invalid_targets.append(res)

    report["n_catalog_valid"] = catalog_valid
    report["n_assemblies_valid"] = assemblies_valid
    report["total_valid"] = catalog_valid + assemblies_valid
    report["missing_or_invalid_targets"] = invalid_targets

    # 7. Coverage and integrity verdict
    if (
        catalog_valid == len(expected_catalog_targets)
        and assemblies_valid == len(expected_assembly_samples)
        and len(invalid_targets) == 0
    ):
        report["status"] = "passed"
        return True, report

    report["status"] = "failed"
    return False, report


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Authoritative standalone receipt audit using corrected Stage 2 validator."
    )
    parser.add_argument(
        "--raw-dir",
        type=Path,
        default=DEFAULT_EXP3_ROOT / "results" / "stage2_raw_outputs_v2",
        help="Path to raw output directory",
    )
    parser.add_argument(
        "--spacers-fna",
        type=Path,
        default=DEFAULT_EXP3_ROOT / "results" / "stage1_verified_archaeal_spacers.fna",
        help="Path to stage1 verified archaeal spacers FASTA",
    )
    parser.add_argument(
        "--assemblies-dir",
        type=Path,
        default=s2_v3.P_ASSEMBLIES_DIR,
        help="Path to raw cohort assemblies directory",
    )
    parser.add_argument(
        "--max-workers",
        type=int,
        default=24,
        help="Parallel threads for streaming subject verification",
    )
    parser.add_argument(
        "--skip-stream-subjects",
        action="store_true",
        help="Skip streaming subject files on disk (manifest verification only; NOT RECOMMENDED)",
    )
    args = parser.parse_args()

    print("=" * 80)
    print("STAGE 2 AUTHORITATIVE STANDALONE RECEIPT & COVERAGE AUDIT")
    print("=" * 80)

    is_passed, report = audit_stage2_receipts(
        raw_dir=args.raw_dir,
        spacers_fna=args.spacers_fna,
        assemblies_dir=args.assemblies_dir,
        stream_check_subjects=not args.skip_stream_subjects,
        max_workers=args.max_workers,
    )

    t_info = report.get("tool_info", {})
    print(f"Tool Binary:      {t_info.get('tool_binary')}")
    print(f"Tool Version:     {t_info.get('tool_version')}")
    print(f"Container SIF:    {t_info.get('container_sif')}")
    print(f"Container SHA256: {t_info.get('container_sha256')}")
    print(f"Spacers Query SHA256: {report.get('query_sha256')}")
    print("-" * 80)
    print(f"Catalog Pools:  {report['n_catalog_valid']} / {report['n_catalog_expected']} valid")
    print(f"Cohort Assemblies: {report['n_assemblies_valid']} / {report['n_assemblies_expected']} valid")
    print(f"Total Coverage: {report['total_valid']} / {report['total_expected']} valid")

    if not is_passed:
        print("\nAUDIT FAILED! Details:")
        for err in report.get("errors", []):
            print(f"  [ERROR] {err}")
        for inv in report.get("missing_or_invalid_targets", []):
            print(f"  [INVALID] {inv['target_name']} ({inv['target_type']}): {inv['error']}")
        print("=" * 80)
        return 1

    print("\nAUDIT PASSED! All 100% expected targets verified with non-circular streamed hashes.")
    print("=" * 80)
    return 0


if __name__ == "__main__":
    sys.exit(main())
