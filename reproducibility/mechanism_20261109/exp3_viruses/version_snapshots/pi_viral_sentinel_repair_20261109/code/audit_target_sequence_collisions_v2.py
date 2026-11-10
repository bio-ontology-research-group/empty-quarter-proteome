#!/usr/bin/env python3
"""Authoritative Target Identity and Sequence Collision Audit (v2 Repair).

Resolves defect where raw assembly headers were extracted using namespaced IDs.
Extracts raw sequences using `raw_subject_id` and `source_sample` from `stage2_target_assemblies_manifest.tsv`,
and extracts catalog sequences using catalog headers.

Audits whether any target ID appears in multiple pools with conflicting sequences,
or if shared target IDs are identical duplicate content across pools.

Requires:
- Exactly all expected pools must be recovered.
- If any expected pool source is missing, marked as `unknown_incomplete_source_coverage`.
- Within-source repeated headers with conflicting sequences are detected and flagged.
- Sequence hashes are checked across all recovered source pools, against the identity table,
  and mapped through SHA-deduplication in `all_matched_primary_targets.fna`.

Outputs:
- `results/stage2_target_identity_collision_audit_v2.tsv`
- `results/stage2_target_identity_collision_summary_v2.json`
- `results/stage2_target_identity_collision_receipt_v2.json`
"""
from __future__ import annotations

import argparse
import collections
import hashlib
import json
import os
import sys
import time
from pathlib import Path
from typing import Any

import pandas as pd

CODE_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(CODE_DIR))

import stage2_viral_identity_and_spacers_v3 as s2_v3

BRANCH_ROOT = Path(
    os.environ.get(
        "EXP3_ROOT",
        "/scratch/project/k10018/eq-mechanism-20261109/exp3_viruses",
    )
)
RES_DIR = BRANCH_ROOT / "results"
RAW_DIR = RES_DIR / "stage2_raw_outputs_v2"


def sha256_bytes(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def stream_file_sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while chunk := f.read(1024 * 1024):
            h.update(chunk)
    return h.hexdigest()


def parse_fasta_for_headers(
    fasta_path: Path, needed_headers: set[str]
) -> tuple[dict[str, str], list[dict[str, Any]]]:
    """Scan a FASTA file once and extract sequences for needed headers.
    Detects duplicate headers with differing sequences.
    """
    if not fasta_path.exists():
        return {}, []

    results: dict[str, str] = {}
    within_source_conflicts: list[dict[str, Any]] = []
    curr_id = None
    curr_seq: list[str] = []

    with open(fasta_path, "r", encoding="ascii", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                if curr_id and curr_id in needed_headers:
                    seq_str = "".join(curr_seq)
                    if curr_id in results and results[curr_id] != seq_str:
                        within_source_conflicts.append({
                            "header": curr_id,
                            "file": str(fasta_path),
                            "len1": len(results[curr_id]),
                            "len2": len(seq_str),
                            "sha1": sha256_bytes(results[curr_id].encode("ascii")),
                            "sha2": sha256_bytes(seq_str.encode("ascii")),
                        })
                    else:
                        results[curr_id] = seq_str
                curr_id = line[1:].strip().split()[0]
                curr_seq = []
            else:
                curr_seq.append(line.strip())

        if curr_id and curr_id in needed_headers:
            seq_str = "".join(curr_seq)
            if curr_id in results and results[curr_id] != seq_str:
                within_source_conflicts.append({
                    "header": curr_id,
                    "file": str(fasta_path),
                    "len1": len(results[curr_id]),
                    "len2": len(seq_str),
                    "sha1": sha256_bytes(results[curr_id].encode("ascii")),
                    "sha2": sha256_bytes(seq_str.encode("ascii")),
                })
            else:
                results[curr_id] = seq_str

    return results, within_source_conflicts


def determine_target_status(
    missing_pools: list[str],
    shas_by_pool: dict[str, str],
    identity_table_sha: str,
    final_sha: str,
    target_conflicts: list[dict[str, Any]],
    final_conflicts: list[dict[str, Any]],
) -> str:
    """Truthful acceptance gate: no single defect can be hidden by identical pools."""
    if missing_pools:
        return "unknown_incomplete_source_coverage"
    if target_conflicts:
        return "conflict_within_source_repeated_header"
    unique_shas = set(shas_by_pool.values())
    if len(unique_shas) != 1:
        return "sequence_collision_across_pools"
    recovered_sha = next(iter(unique_shas))
    if len(identity_table_sha) != 64 or recovered_sha != identity_table_sha:
        return "conflict_with_identity_table"
    if final_conflicts or len(final_sha) != 64 or final_sha != recovered_sha:
        return "conflict_with_final_merged_fasta"
    return "identical_duplicate_across_pools"


def collision_free_verdict(audit_rows: list[dict[str, Any]], within_source_conflicts: list[dict[str, Any]], final_conflicts: list[dict[str, Any]]) -> bool:
    return (
        bool(audit_rows)
        and not within_source_conflicts
        and not final_conflicts
        and all(
            row["status"] == "identical_duplicate_across_pools"
            and row["has_complete_source_coverage"]
            and row["matches_identity_table"]
            and row["final_fna_sha"] == row["identity_table_sha"]
            for row in audit_rows
        )
    )


def main() -> int:
    t0 = time.time()
    scored_path = RES_DIR / "stage2_all_spacer_alignments_scored.tsv"
    if not scored_path.exists():
        print(f"ERROR: {scored_path} not found")
        return 1

    manifest_path = RES_DIR / "stage2_target_assemblies_manifest.tsv"
    if not manifest_path.exists():
        print(f"ERROR: {manifest_path} not found")
        return 1

    id_table_path = RES_DIR / "stage2_target_viral_identity_table.tsv"
    if not id_table_path.exists():
        print(f"ERROR: {id_table_path} not found")
        return 1

    df_scored = pd.read_csv(scored_path, sep="\t")
    print(f"[Audit v2] Total alignments in scored table: {len(df_scored)}")

    # 1. Identify all multi-pool targets
    pools_per_target = df_scored.groupby("target_id")["pool"].apply(lambda s: sorted(set(s)))
    total_unique_targets = len(pools_per_target)
    multi_pool_targets = {t_id: set(pools) for t_id, pools in pools_per_target.items() if len(pools) > 1}
    print(f"[Audit v2] Total unique targets: {total_unique_targets}")
    print(f"[Audit v2] Targets in multiple pools: {len(multi_pool_targets)}")

    multi_target_ids = set(multi_pool_targets.keys())
    df_multi = df_scored[df_scored["target_id"].isin(multi_target_ids)]

    # Load assembly manifest
    df_manifest = pd.read_csv(manifest_path, sep="\t")
    required_manifest = {"sample_id", "assembly_path", "assembly_exists", "subject_sha256", "receipt_verified"}
    if not required_manifest <= set(df_manifest.columns):
        raise ValueError(f"assembly manifest missing columns: {sorted(required_manifest-set(df_manifest.columns))}")
    asm_path_map = dict(zip(df_manifest["sample_id"], df_manifest["assembly_path"]))
    manifest_by_sample = df_manifest.set_index("sample_id", drop=False)
    source_paths_to_hash: set[Path] = set()

    # Catalog pool file paths
    catalog_pool_paths = {
        "eq_votu_reps": s2_v3.P_VOTU_FNA,
        "eq_all_viruses": s2_v3.P_ALL_VIR_FNA,
        "eq_coasm_viruses": s2_v3.P_COASM_VIR_FNA,
        "curated_plasmid_controls": s2_v3.P_ALL_PLAS_FNA,
        "curated_cellular_controls": RAW_DIR / "curated_cellular_controls_verified.fna",
        "external_habitat_viral_controls": s2_v3.P_EXT_REF_VIR_FNA,
    }

    # Store extracted sequences: pool_target_seqs[pool_name][target_id] = seq
    pool_target_seqs: dict[str, dict[str, str]] = collections.defaultdict(dict)
    all_within_source_conflicts: list[dict[str, Any]] = []

    # 2. Extract from Catalog Pools
    for pool_name, p_path in catalog_pool_paths.items():
        sub_pool = df_multi[df_multi["pool"] == pool_name]
        if sub_pool.empty:
            continue
        needed_target_ids = set(sub_pool["target_id"])
        print(f"[Audit v2] Extracting {len(needed_target_ids)} targets from catalog pool '{pool_name}' ({p_path})...")
        source_paths_to_hash.add(p_path)
        ext, conflicts = parse_fasta_for_headers(p_path, needed_target_ids)
        for conflict in conflicts:
            conflict["pool"] = pool_name
            conflict["target_ids"] = [conflict["header"]]
        all_within_source_conflicts.extend(conflicts)
        for t_id, seq in ext.items():
            pool_target_seqs[pool_name][t_id] = seq

    # 3. Extract from Raw Assemblies (grouped by sample so each assembly file is scanned ONCE)
    sub_raw = df_multi[df_multi["pool"] == "eq_raw_assemblies_arch_samples"]
    if not sub_raw.empty:
        print(f"[Audit v2] Extracting {sub_raw['target_id'].nunique()} targets from raw assemblies...")
        grouped_by_sample = sub_raw.groupby("source_sample")
        for samp, grp in grouped_by_sample:
            samp_str = str(samp)
            if samp_str not in asm_path_map:
                print(f"[ERROR] Sample {samp_str} not in assembly manifest!")
                continue
            asm_file = Path(asm_path_map[samp_str])
            source_paths_to_hash.add(asm_file)
            manifest_row = manifest_by_sample.loc[samp_str]
            if int(manifest_row["assembly_exists"]) != 1 or int(manifest_row["receipt_verified"]) != 1:
                raise ValueError(f"assembly source is not receipt verified: {samp_str}")
            # Map raw_subject_id to target_id without silently collapsing conflicts.
            raw_to_target: dict[str, str] = {}
            for raw_id, target_id in zip(grp["raw_subject_id"], grp["target_id"]):
                if raw_id in raw_to_target and raw_to_target[raw_id] != target_id:
                    raise ValueError(f"raw header maps to multiple target IDs: {samp_str}/{raw_id}")
                raw_to_target[raw_id] = target_id
            needed_raw = set(raw_to_target.keys())

            ext_raw, conflicts = parse_fasta_for_headers(asm_file, needed_raw)
            for conflict in conflicts:
                conflict["pool"] = "eq_raw_assemblies_arch_samples"
                conflict["target_ids"] = [f"{samp_str}::{conflict['header']}"]
            all_within_source_conflicts.extend(conflicts)

            for raw_id, seq in ext_raw.items():
                expected_target_id = raw_to_target[raw_id]
                # Apply stage2 namespacing transform: f"{samp_str}::{raw_id}"
                namespaced_id = f"{samp_str}::{raw_id}"
                if namespaced_id != expected_target_id:
                    raise ValueError(f"namespacing/source identity mismatch: {namespaced_id} vs {expected_target_id}")
                pool_target_seqs["eq_raw_assemblies_arch_samples"][expected_target_id] = seq

    # 4. Extract from all_matched_primary_targets.fna and index SHA-256 for alias resolution
    final_fna = RAW_DIR / "all_matched_primary_targets.fna"
    final_fna_seqs: dict[str, str] = {}
    final_fna_shas: dict[str, str] = {}
    sha_to_final_headers: dict[str, list[str]] = collections.defaultdict(list)
    final_fna_conflicts: list[dict[str, Any]] = []

    if final_fna.exists():
        source_paths_to_hash.add(final_fna)
        print(f"[Audit v2] Indexing {final_fna.name}...")
        curr_id = None
        curr_seq: list[str] = []
        with open(final_fna, "r", encoding="ascii", errors="replace") as fh:
            for line in fh:
                if line.startswith(">"):
                    if curr_id:
                        s = "".join(curr_seq)
                        s_sha = sha256_bytes(s.encode("ascii"))
                        if curr_id in final_fna_shas and final_fna_shas[curr_id] != s_sha:
                            final_fna_conflicts.append({"header": curr_id, "sha1": final_fna_shas[curr_id], "sha2": s_sha})
                        final_fna_seqs[curr_id] = s
                        final_fna_shas[curr_id] = s_sha
                        if curr_id not in sha_to_final_headers[s_sha]:
                            sha_to_final_headers[s_sha].append(curr_id)
                    curr_id = line[1:].strip().split()[0]
                    curr_seq = []
                else:
                    curr_seq.append(line.strip())
            if curr_id:
                s = "".join(curr_seq)
                s_sha = sha256_bytes(s.encode("ascii"))
                if curr_id in final_fna_shas and final_fna_shas[curr_id] != s_sha:
                    final_fna_conflicts.append({"header": curr_id, "sha1": final_fna_shas[curr_id], "sha2": s_sha})
                final_fna_seqs[curr_id] = s
                final_fna_shas[curr_id] = s_sha
                if curr_id not in sha_to_final_headers[s_sha]:
                    sha_to_final_headers[s_sha].append(curr_id)

    # 5. Load Identity Table
    df_id = pd.read_csv(id_table_path, sep="\t")
    id_table_shas = dict(zip(df_id["target_id"].astype(str), df_id["target_seq_sha256"].astype(str)))

    # 6. Evaluate all 144 multi-pool targets
    audit_rows: list[dict[str, Any]] = []
    n_identical_duplicates = 0
    n_sequence_collisions = 0
    n_incomplete_coverage = 0
    n_identity_table_mismatch = 0
    n_final_fasta_mismatch = 0
    n_within_source_target_conflicts = 0

    for t_id in sorted(multi_pool_targets.keys()):
        expected_pools = sorted(multi_pool_targets[t_id])
        recovered_pools = sorted([p for p in expected_pools if t_id in pool_target_seqs.get(p, {})])
        missing_pools = sorted(set(expected_pools) - set(recovered_pools))

        shas_by_pool: dict[str, str] = {}
        lengths_by_pool: dict[str, int] = {}
        for p in recovered_pools:
            seq = pool_target_seqs[p][t_id]
            shas_by_pool[p] = sha256_bytes(seq.encode("ascii"))
            lengths_by_pool[p] = len(seq)

        unique_shas = set(shas_by_pool.values())
        id_tbl_sha = id_table_shas.get(t_id, "missing")

        final_sha = final_fna_shas.get(t_id, "missing")
        final_alias_match = False
        if final_sha == "missing" and len(unique_shas) == 1:
            rec_sha = next(iter(unique_shas))
            if rec_sha in sha_to_final_headers:
                final_sha = rec_sha
                final_alias_match = True

        target_conflicts = [c for c in all_within_source_conflicts if t_id in c.get("target_ids", [])]
        target_final_conflicts = [c for c in final_fna_conflicts if c["header"] == t_id]
        status = determine_target_status(
            missing_pools, shas_by_pool, id_tbl_sha, final_sha, target_conflicts, target_final_conflicts,
        )
        if status == "unknown_incomplete_source_coverage":
            n_incomplete_coverage += 1
        elif status == "identical_duplicate_across_pools":
            n_identical_duplicates += 1
        elif status == "sequence_collision_across_pools":
            n_sequence_collisions += 1
        elif status == "conflict_with_identity_table":
            n_identity_table_mismatch += 1
        elif status == "conflict_with_final_merged_fasta":
            n_final_fasta_mismatch += 1
        elif status == "conflict_within_source_repeated_header":
            n_within_source_target_conflicts += 1

        audit_rows.append({
            "target_id": t_id,
            "expected_pools": ";".join(expected_pools),
            "recovered_pools": ";".join(recovered_pools),
            "missing_pools": ";".join(missing_pools),
            "n_expected_pools": len(expected_pools),
            "n_recovered_pools": len(recovered_pools),
            "has_complete_source_coverage": (len(missing_pools) == 0),
            "shas_by_pool": json.dumps(shas_by_pool),
            "lengths_by_pool": json.dumps(lengths_by_pool),
            "unique_sha_count": len(unique_shas),
            "final_fna_sha": final_sha,
            "final_fna_dedup_alias_resolved": final_alias_match,
            "identity_table_sha": id_tbl_sha,
            "matches_identity_table": (len(unique_shas) == 1 and id_tbl_sha in unique_shas),
            "within_source_conflicts": json.dumps(target_conflicts),
            "final_fasta_conflicts": json.dumps(target_final_conflicts),
            "status": status,
        })

    out_audit_df = pd.DataFrame(audit_rows)
    out_audit_tsv_v2 = RES_DIR / "stage2_target_identity_collision_audit_v2.tsv"
    out_audit_df.to_csv(out_audit_tsv_v2, sep="\t", index=False)

    # Also update baseline path (preserving previous as v1_flawed)
    baseline_tsv = RES_DIR / "stage2_target_identity_collision_audit.tsv"
    if baseline_tsv.exists():
        v1_tsv = RES_DIR / "stage2_target_identity_collision_audit_v1_flawed.tsv"
        if not v1_tsv.exists():
            v1_tsv.write_bytes(baseline_tsv.read_bytes())
    out_audit_df.to_csv(baseline_tsv, sep="\t", index=False)

    summary = {
        "audit_version": "v2_repaired_source_extraction",
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "total_targets_evaluated": total_unique_targets,
        "n_multi_pool_targets": len(multi_pool_targets),
        "n_complete_source_coverage": len(multi_pool_targets) - n_incomplete_coverage,
        "n_incomplete_source_coverage": n_incomplete_coverage,
        "n_identical_duplicates": n_identical_duplicates,
        "n_sequence_collisions": n_sequence_collisions,
        "n_identity_table_mismatches": n_identity_table_mismatch,
        "n_final_fasta_mismatches": n_final_fasta_mismatch,
        "n_targets_with_within_source_conflicts": n_within_source_target_conflicts,
        "n_within_source_repeated_headers": len(all_within_source_conflicts),
        "within_source_conflicts": all_within_source_conflicts,
        "final_fasta_conflicts": final_fna_conflicts,
        "collision_free": collision_free_verdict(audit_rows, all_within_source_conflicts, final_fna_conflicts),
    }

    out_summary_json_v2 = RES_DIR / "stage2_target_identity_collision_summary_v2.json"
    out_summary_json_v2.write_text(json.dumps(summary, indent=2) + "\n")

    baseline_json = RES_DIR / "stage2_target_identity_collision_summary.json"
    if baseline_json.exists():
        v1_json = RES_DIR / "stage2_target_identity_collision_summary_v1_flawed.json"
        if not v1_json.exists():
            v1_json.write_bytes(baseline_json.read_bytes())
    baseline_json.write_text(json.dumps(summary, indent=2) + "\n")

    # Bind every original source to its actual bytes.  Raw assemblies additionally
    # must agree with the SHA in their receipt-verified assembly manifest.
    source_hashes: dict[str, str] = {}
    for source_path in sorted(source_paths_to_hash, key=str):
        if not source_path.exists():
            raise FileNotFoundError(f"original source missing during receipt hashing: {source_path}")
        observed_sha = stream_file_sha256(source_path)
        source_hashes[str(source_path)] = observed_sha
    for sample, assembly_path in asm_path_map.items():
        path = Path(assembly_path)
        if path in source_paths_to_hash:
            expected_sha = str(manifest_by_sample.loc[sample, "subject_sha256"])
            if source_hashes[str(path)] != expected_sha:
                raise ValueError(f"assembly SHA differs from manifest for {sample}: {path}")

    # Cryptographic Receipt
    receipt = {
        "status": "completed" if summary["collision_free"] else "verified_with_findings",
        "stage": "stage2_target_identity_collision_audit_v2",
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "elapsed_seconds": round(time.time() - t0, 2),
        "inputs": {
            "scored_alignments": str(scored_path),
            "scored_alignments_sha256": stream_file_sha256(scored_path),
            "assemblies_manifest": str(manifest_path),
            "assemblies_manifest_sha256": stream_file_sha256(manifest_path),
            "identity_table": str(id_table_path),
            "identity_table_sha256": stream_file_sha256(id_table_path),
            "original_source_sha256": source_hashes,
            "original_source_hashes_verified_against_assembly_manifest": True,
        },
        "outputs": {
            "audit_tsv_v2": str(out_audit_tsv_v2),
            "audit_tsv_v2_sha256": stream_file_sha256(out_audit_tsv_v2),
            "summary_json_v2": str(out_summary_json_v2),
            "summary_json_v2_sha256": stream_file_sha256(out_summary_json_v2),
        },
        "metrics": summary,
    }
    receipt_p = RES_DIR / "stage2_target_identity_collision_receipt_v2.json"
    receipt_p.write_text(json.dumps(receipt, indent=2) + "\n")

    print("=" * 80)
    print("TARGET IDENTITY AND SEQUENCE COLLISION AUDIT REPORT (v2 REPAIRED)")
    print("=" * 80)
    print(f"Total Unique Targets:         {total_unique_targets}")
    print(f"Multi-Pool Shared Targets:    {len(multi_pool_targets)}")
    print(f"Complete Source Coverage:     {len(multi_pool_targets) - n_incomplete_coverage} / {len(multi_pool_targets)}")
    print(f"Incomplete Source Coverage:   {n_incomplete_coverage}")
    print(f"Identical Duplicates:         {n_identical_duplicates} (100% sequence equivalence)")
    print(f"Sequence Collisions:          {n_sequence_collisions}")
    print(f"Identity Table Mismatches:    {n_identity_table_mismatch}")
    print(f"Within-Source Conflicting IDs:{len(all_within_source_conflicts)}")
    print(f"Collision-Free Verdict:       {summary['collision_free']}")
    print(f"Receipt written to:           {receipt_p}")
    print("=" * 80)

    return 0 if summary["collision_free"] else 1


if __name__ == "__main__":
    sys.exit(main())
