#!/usr/bin/env python3
"""Aggregate only cryptographically receipted competitive-mapping samples.

Every logical reviewed host+raw-SHA pair is evaluated in every completed mapped
sample.  Missing rows/errors are never converted to biological zeros.  inStrain
status comes from actual profile results, not eligibility alone.
"""
from __future__ import annotations

import argparse
import json
import os
import time
from pathlib import Path
from typing import Any

import pandas as pd

from stage4_competitive_common import atomic_json, atomic_text, load_analysis_pairs, stream_sha256

BRANCH_ROOT = Path(os.environ.get("EXP3_ROOT", "/scratch/project/k10018/eq-mechanism-20261109/exp3_viruses"))
RES_DIR = BRANCH_ROOT / "results"
MAP_DIR = RES_DIR / "stage4_competitive_mapping"


def one_row(frame: pd.DataFrame, **conditions: Any) -> pd.Series:
    selected = frame
    for column, value in conditions.items():
        selected = selected[selected[column] == value]
    if len(selected) != 1:
        raise ValueError(f"expected exactly one evaluation row for {conditions}, found={len(selected)}; missing is not zero")
    return selected.iloc[0]


def verify_sample_receipt(receipt_path: Path) -> dict[str, Any]:
    receipt = json.loads(receipt_path.read_text())
    if receipt.get("status") != "completed_verified" or receipt.get("schema_version") != 3:
        raise ValueError(f"invalid mapping receipt status/schema: {receipt_path}")
    sample = receipt.get("sample")
    if not sample or receipt_path.parent.name != sample:
        raise ValueError(f"sample/path mismatch in mapping receipt: {receipt_path}")
    done = receipt_path.parent / "MAPPING_DONE"
    if not done.exists():
        raise ValueError(f"MAPPING_DONE missing for {sample}")
    done_payload = json.loads(done.read_text())
    if done_payload.get("receipt_sha256") != stream_sha256(receipt_path):
        raise ValueError(f"MAPPING_DONE does not match current receipt for {sample}")
    outputs = receipt["outputs"]
    for path_key, sha_key in [
        ("coverage_tsv", "coverage_tsv_sha256"),
        ("focal_evaluation_tsv", "focal_evaluation_tsv_sha256"),
        ("raw_bam", "raw_bam_sha256"),
        ("filtered_bam", "filtered_bam_sha256"),
    ]:
        path = Path(outputs[path_key])
        if not path.exists() or path.stat().st_size == 0:
            raise FileNotFoundError(f"receipted output missing/empty: {path}")
        if stream_sha256(path) != outputs[sha_key]:
            raise ValueError(f"receipted output hash mismatch: {path}")
    evaluation = pd.read_csv(outputs["focal_evaluation_tsv"], sep="\t")
    if evaluation.empty or set(evaluation["sample"].astype(str)) != {sample}:
        raise ValueError(f"invalid focal evaluation table for {sample}")
    receipt["_evaluation_frame"] = evaluation
    return receipt


def main() -> int:
    start = time.time()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--map-dir", type=Path, default=MAP_DIR)
    parser.add_argument("--pair-ledger", type=Path, default=RES_DIR / "stage3_tier1_pairs_individual_evidence_audit.tsv")
    parser.add_argument("--boundary-config", type=Path, default=BRANCH_ROOT / "inputs/stage4_boundary_sensitive_intervals.json")
    parser.add_argument("--out-prefix", type=Path, default=RES_DIR / "stage4_competitive_mapping_all")
    args = parser.parse_args()
    (RES_DIR / "STAGE4_COMPETITIVE_MAPPING_DONE").unlink(missing_ok=True)
    receipt_paths = sorted(args.map_dir.glob("*/*_mapping_receipt.json"))
    if not receipt_paths:
        raise FileNotFoundError(f"no mapping receipts in {args.map_dir}")
    receipts = [verify_sample_receipt(path) for path in receipt_paths]
    samples = [receipt["sample"] for receipt in receipts]
    if len(samples) != len(set(samples)):
        raise ValueError("duplicate sample mapping receipts")
    pairs, reviewed_pairs, boundary_pairs = load_analysis_pairs(args.pair_ledger, args.boundary_config)

    evaluations = pd.concat([receipt["_evaluation_frame"] for receipt in receipts], ignore_index=True)
    evaluation_path = args.out_prefix.with_name(args.out_prefix.name + "_contigs_evaluation.tsv")
    pair_path = args.out_prefix.with_name(args.out_prefix.name + "_pairs_summary.tsv")
    boundary_path = args.out_prefix.with_name(args.out_prefix.name + "_boundary_intervals.tsv")
    summary_path = args.out_prefix.with_name(args.out_prefix.name + "_pairs_summary.json")
    evaluations.to_csv(evaluation_path, sep="\t", index=False)
    boundary_evaluations = evaluations[evaluations["entity_type"] == "boundary_interval"].copy()
    expected_interval_rows = len(boundary_pairs) * 3 * len(receipts)
    if len(boundary_evaluations) != expected_interval_rows:
        raise ValueError(f"boundary interval coverage incomplete: {len(boundary_evaluations)} != {expected_interval_rows}")
    boundary_evaluations.to_csv(boundary_path, sep="\t", index=False)

    rows: list[dict[str, Any]] = []
    for pair in pairs.to_dict(orient="records"):
        for receipt in receipts:
            sample = receipt["sample"]
            frame = receipt["_evaluation_frame"]
            host = one_row(frame, sample=sample, entity_type="host_mag", entity_id=pair["host_mag"])
            virus = one_row(frame, sample=sample, entity_type="virus_target", entity_sha256=pair["target_seq_sha256"])
            if virus["analysis_role"] != pair["analysis_role"]:
                raise ValueError(f"analysis role mismatch for {pair['target_seq_sha256']} in {sample}")
            profile = next(
                (result for result in receipt["instrain"]["results"] if result["target_raw_sha256"] == pair["target_seq_sha256"]),
                None,
            )
            qualifies = bool(int(virus["qualifies_10x_50pct_instrain"]))
            if bool(profile) != (qualifies and pair["analysis_role"] == "primary"):
                raise ValueError(f"inStrain eligibility/result mismatch for {pair['target_seq_sha256']} in {sample}")
            roles = []
            if sample == pair["host_sample"]:
                roles.append("host_source")
            if sample in pair["virus_samples"]:
                roles.append("viral_or_alias_source")
            if sample == "10Dr2":
                roles.append("predetermined_unrelated_control")
            if not roles:
                roles.append("predetermined_coassembly_alias_or_copositive_sample")
            rows.append({
                "host_mag": pair["host_mag"],
                "host_phylum": pair["host_phylum"],
                "target_virus": pair["target_id"],
                "target_aliases": ";".join(pair["target_aliases"]),
                "target_seq_sha256": pair["target_seq_sha256"],
                "analysis_role": pair["analysis_role"],
                "individual_virality_status": pair["individual_virality_status"],
                "sample": sample,
                "sample_roles": ";".join(roles),
                "host_total_mean_depth": float(host["mean_depth"]),
                "host_total_breadth_pct": float(host["horizontal_breadth_pct"]),
                "host_unique_mean_depth": float(host["unique_mean_depth"]),
                "host_unique_breadth_pct": float(host["unique_horizontal_breadth_pct"]),
                "host_unique_reads": int(host["unique_filtered_reads_mapped"]),
                "host_ambiguous_reads": int(host["ambiguous_filtered_reads_mapped"]),
                "virus_mean_depth": float(virus["mean_depth"]),
                "virus_breadth_pct": float(virus["horizontal_breadth_pct"]),
                "virus_uniquely_attributable": bool(virus["uniquely_attributable"]),
                "virus_unique_reads": int(virus["unique_filtered_reads_mapped"]),
                "virus_ambiguous_reads": int(virus["ambiguous_filtered_reads_mapped"]),
                "virus_ambiguity_status": virus["ambiguity_status"],
                "qualifies_10x_50pct_instrain": int(qualifies),
                "instrain_status": profile["status"] if profile else "not_profiled_not_eligible",
                "instrain_callable_positions": profile["callable_positions"] if profile else pd.NA,
                "instrain_analysis_mode": profile["analysis_mode"] if profile else pd.NA,
            })
    compiled = pd.DataFrame(rows)
    compiled.to_csv(pair_path, sep="\t", index=False)
    profile_count = sum(len(receipt["instrain"]["results"]) for receipt in receipts)
    summary = {
        "status": "completed_verified",
        "schema_version": 3,
        "stage": "stage4_aggregate_competitive_mapping",
        "timestamp_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "elapsed_seconds": time.time() - start,
        "n_samples_mapped": len(receipts),
        "samples": samples,
        "n_primary_host_target_sha_pairs": int((reviewed_pairs["analysis_role"] == "primary").sum()),
        "n_sensitivity_host_target_sha_pairs": int((reviewed_pairs["analysis_role"] != "primary").sum()),
        "n_boundary_sensitive_provisional_pairs": len(boundary_pairs),
        "n_boundary_interval_sample_rows": len(boundary_evaluations),
        "n_pair_sample_combinations": len(compiled),
        "n_qualifying_unique_viral_sha_sample_combinations": int(
            evaluations.loc[evaluations["entity_type"] == "virus_target", "qualifies_10x_50pct_instrain"].sum()
        ),
        "n_instrain_profiles_completed": profile_count,
        "interpretation_limits": {
            "supported": "competitive read recruitment and sequence association conditional on ownership ambiguity",
            "not_established": "active infection, replication, mortality, causality, or recombination",
            "unresolved_41Sr3": "sensitivity only; outside primary credible-virus analysis",
            "boundary_43PRr3_44Dr2": "provisional only; core/disputed extension/host flank reported separately; nested-competitor low coverage is not absence",
            "variant_note": "inStrain variants are sequence polymorphism, not recombination evidence",
        },
        "inputs": {
            "pair_ledger": str(args.pair_ledger), "pair_ledger_sha256": stream_sha256(args.pair_ledger),
            "boundary_config": str(args.boundary_config), "boundary_config_sha256": stream_sha256(args.boundary_config),
            "mapping_receipts": {receipt["sample"]: stream_sha256(args.map_dir / receipt["sample"] / f"{receipt['sample']}_mapping_receipt.json") for receipt in receipts},
            "aggregation_script_sha256": stream_sha256(Path(__file__)),
        },
        "outputs": {
            "evaluation_tsv": str(evaluation_path), "evaluation_tsv_sha256": stream_sha256(evaluation_path),
            "pairs_tsv": str(pair_path), "pairs_tsv_sha256": stream_sha256(pair_path),
            "boundary_intervals_tsv": str(boundary_path), "boundary_intervals_tsv_sha256": stream_sha256(boundary_path),
        },
    }
    atomic_json(summary_path, summary)
    atomic_text(RES_DIR / "STAGE4_COMPETITIVE_MAPPING_DONE", json.dumps({"summary_sha256": stream_sha256(summary_path)}) + "\n")
    print(f"[Aggregation completed verified] {summary_path}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
