#!/usr/bin/env python3
"""Authoritative Target Identity and Sequence Collision Audit.

Audits whether any target ID appears in multiple pools with conflicting sequences,
or if shared target IDs are identical duplicate content across pools.

Produces:
1. `results/stage2_target_identity_collision_audit.tsv`:
   Row-by-row comparison for all targets appearing in multiple pools.
2. `results/stage2_target_identity_collision_summary.json`:
   Summary metrics on shared targets, identical duplicates, and collisions.
"""
from __future__ import annotations

import collections
import hashlib
import json
import os
import sys
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


def parse_fasta_for_targets(fasta_path: Path, target_ids: set[str]) -> dict[str, str]:
    """Extract sequences for a specific set of target IDs from a FASTA file."""
    if not fasta_path.exists():
        return {}
    results: dict[str, str] = {}
    curr_id = None
    curr_seq = []
    with open(fasta_path, "r", encoding="ascii", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                if curr_id and curr_id in target_ids:
                    results[curr_id] = "".join(curr_seq)
                curr_id = line[1:].strip().split()[0]
                curr_seq = []
            else:
                curr_seq.append(line.strip())
        if curr_id and curr_id in target_ids:
            results[curr_id] = "".join(curr_seq)
    return results


def main() -> int:
    scored_path = RES_DIR / "stage2_all_spacer_alignments_scored.tsv"
    if not scored_path.exists():
        print(f"ERROR: {scored_path} not found")
        return 1

    df = pd.read_csv(scored_path, sep="\t")
    print(f"[Audit] Scored alignments: {len(df)}")
    target_pools: dict[str, set[str]] = collections.defaultdict(set)
    target_samples: dict[str, set[str]] = collections.defaultdict(set)

    for _, r in df.iterrows():
        t_id = str(r["target_id"])
        pool = str(r["pool"])
        samp = str(r.get("source_sample", ""))
        target_pools[t_id].add(pool)
        if samp and samp != "nan":
            target_samples[t_id].add(samp)

    total_unique_targets = len(target_pools)
    multi_pool_targets = {t_id: pools for t_id, pools in target_pools.items() if len(pools) > 1}
    print(f"[Audit] Unique targets: {total_unique_targets}")
    print(f"[Audit] Targets in multiple pools: {len(multi_pool_targets)}")

    # Catalog pool file paths
    catalog_pool_paths = {
        "eq_votu_reps": s2_v3.P_VOTU_FNA,
        "eq_all_viruses": s2_v3.P_ALL_VIR_FNA,
        "eq_coasm_viruses": s2_v3.P_COASM_VIR_FNA,
        "curated_plasmid_controls": s2_v3.P_ALL_PLAS_FNA,
        "curated_cellular_controls": RAW_DIR / "curated_cellular_controls_verified.fna",
        "external_habitat_viral_controls": s2_v3.P_EXT_REF_VIR_FNA,
    }

    # Extract sequences per pool for multi-pool targets
    pool_target_seqs: dict[str, dict[str, str]] = collections.defaultdict(dict)
    all_multi_ids = set(multi_pool_targets.keys())

    for pool_name, p_path in catalog_pool_paths.items():
        needed_for_pool = {t for t in all_multi_ids if pool_name in target_pools[t]}
        if needed_for_pool:
            print(f"[Audit] Extracting {len(needed_for_pool)} shared targets from {pool_name} ({p_path.name})...")
            extracted = parse_fasta_for_targets(p_path, needed_for_pool)
            pool_target_seqs[pool_name].update(extracted)

    # For raw assemblies, extract from the specific sample assembly
    asm_needed = {t for t in all_multi_ids if "eq_raw_assemblies_arch_samples" in target_pools[t]}
    if asm_needed:
        print(f"[Audit] Extracting {len(asm_needed)} shared targets from assembly files...")
        for t_id in asm_needed:
            samps = target_samples.get(t_id, set())
            for s in samps:
                asm_fa = s2_v3.P_ASSEMBLIES_DIR / s / f"{s}_contigs.fa"
                if asm_fa.exists():
                    ext = parse_fasta_for_targets(asm_fa, {t_id})
                    if t_id in ext:
                        pool_target_seqs["eq_raw_assemblies_arch_samples"][t_id] = ext[t_id]
                        break

    # Read the final assembled all_matched_primary_targets.fna
    final_fna = RAW_DIR / "all_matched_primary_targets.fna"
    final_seqs = {}
    if final_fna.exists():
        final_seqs = parse_fasta_for_targets(final_fna, all_multi_ids)

    # Read stage2_target_viral_identity_table.tsv
    vdf = pd.read_csv(RES_DIR / "stage2_target_viral_identity_table.tsv", sep="\t")
    id_table_shas = dict(zip(vdf["target_id"].astype(str), vdf["target_seq_sha256"].astype(str)))

    # Audit each shared target
    audit_rows: list[dict[str, Any]] = []
    n_identical_duplicates = 0
    n_collisions = 0
    n_single_source_recovered = 0

    for t_id, pools in sorted(multi_pool_targets.items()):
        shas_by_pool: dict[str, str] = {}
        lengths_by_pool: dict[str, int] = {}
        for p in pools:
            seq = pool_target_seqs.get(p, {}).get(t_id)
            if seq:
                shas_by_pool[p] = s2_v3.sha256_bytes(seq.encode("ascii"))
                lengths_by_pool[p] = len(seq)

        unique_shas = set(shas_by_pool.values())
        final_sha = s2_v3.sha256_bytes(final_seqs[t_id].encode("ascii")) if t_id in final_seqs else "missing"
        table_sha = id_table_shas.get(t_id, "missing")

        if len(unique_shas) == 1:
            status = "identical_duplicate_across_pools"
            n_identical_duplicates += 1
        elif len(unique_shas) > 1:
            status = "sequence_collision_across_pools"
            n_collisions += 1
        else:
            status = "single_pool_recovered"
            n_single_source_recovered += 1

        audit_rows.append(
            {
                "target_id": t_id,
                "pools": ";".join(sorted(pools)),
                "n_pools": len(pools),
                "shas_by_pool": json.dumps(shas_by_pool),
                "lengths_by_pool": json.dumps(lengths_by_pool),
                "unique_sha_count": len(unique_shas),
                "final_fna_sha": final_sha,
                "identity_table_sha": table_sha,
                "matches_final_fna": (final_sha in unique_shas) if unique_shas else False,
                "status": status,
            }
        )

    out_audit_tsv = RES_DIR / "stage2_target_identity_collision_audit.tsv"
    pd.DataFrame(audit_rows).to_csv(out_audit_tsv, sep="\t", index=False)

    summary = {
        "total_targets_evaluated": total_unique_targets,
        "n_multi_pool_targets": len(multi_pool_targets),
        "n_identical_duplicates": n_identical_duplicates,
        "n_sequence_collisions": n_collisions,
        "n_single_source_recovered": n_single_source_recovered,
        "collision_free": (n_collisions == 0),
    }

    out_summary_json = RES_DIR / "stage2_target_identity_collision_summary.json"
    out_summary_json.write_text(json.dumps(summary, indent=2) + "\n")

    print("=" * 80)
    print("TARGET IDENTITY AND SEQUENCE COLLISION AUDIT REPORT")
    print("=" * 80)
    print(f"Total Unique Targets:      {total_unique_targets}")
    print(f"Multi-Pool Shared Targets: {len(multi_pool_targets)}")
    print(f"Identical Duplicates:      {n_identical_duplicates} (100% sequence equivalence)")
    print(f"Sequence Collisions:       {n_collisions}")
    print(f"Single-Source Recovered:   {n_single_source_recovered}")
    print(f"Collision-Free Verdict:    {summary['collision_free']}")
    print("=" * 80)

    return 0 if n_collisions == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
