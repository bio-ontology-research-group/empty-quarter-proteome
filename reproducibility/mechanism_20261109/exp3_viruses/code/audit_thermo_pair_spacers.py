#!/usr/bin/env python3
"""Forensic audit of Thermoplasmatota host 41PRr2__single_SemiBin_236 vs V60PRr2::V60PRr2_715976.

Implements exact schema join against stage1 verified archaeal spacers:
Columns: unique_spacer_id, mag_id, is_arch591_rep, canonical_spacer_seq_rc,
         representative_spacer_seq, spacer_len, primary_namespaced_contig_id,
         primary_array_start, primary_array_end, all_array_loci, repeat_seq.

Deduplicates triplicated pool alignments (eq_votu_reps, eq_all_viruses,
eq_raw_assemblies_arch_samples) by (unique_spacer_id, protospacer_start,
protospacer_end, target_seq_sha256).

Fails explicitly on missing or empty sequence join.
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path
from typing import Any

import pandas as pd

BRANCH_ROOT = Path(
    os.environ.get(
        "EXP3_ROOT",
        "/scratch/project/k10018/eq-mechanism-20261109/exp3_viruses",
    )
)
RES_DIR = BRANCH_ROOT / "results"
RAW_DIR = RES_DIR / "stage2_raw_outputs_v2"


def audit_thermo_spacers(
    scored_p: Path,
    spc_p: Path,
    out_tsv: Path | None = None,
    out_json: Path | None = None,
) -> dict[str, Any]:
    """Audit Thermoplasmatota pair spacers with explicit schema and deduplication."""
    if not scored_p.exists():
        raise FileNotFoundError(f"Scored alignments file not found: {scored_p}")
    if not spc_p.exists():
        raise FileNotFoundError(f"Verified spacers file not found: {spc_p}")

    df = pd.read_csv(scored_p, sep="\t")
    spc_df = pd.read_csv(spc_p, sep="\t")

    # Validate stage1 schema
    required_spc_cols = {
        "unique_spacer_id",
        "mag_id",
        "canonical_spacer_seq_rc",
        "representative_spacer_seq",
        "spacer_len",
        "primary_namespaced_contig_id",
        "primary_array_start",
        "primary_array_end",
    }
    missing_cols = required_spc_cols - set(spc_df.columns)
    if missing_cols:
        raise ValueError(f"stage1 spacers table missing required schema columns: {sorted(missing_cols)}")

    spc_meta = spc_df.set_index("unique_spacer_id").to_dict(orient="index")

    # Filter for Thermoplasmatota pair
    m = df[
        (df["host_mag"] == "41PRr2__single_SemiBin_236")
        & (df["target_id"] == "V60PRr2::V60PRr2_715976")
    ].sort_values(["passes_strict_spacer_gate", "sstart"], ascending=[False, True])

    raw_match_count = len(m)
    if raw_match_count == 0:
        raise ValueError("No alignments found for 41PRr2__single_SemiBin_236 vs V60PRr2::V60PRr2_715976")

    # Deduplicate triplicated alignments across pools by (unique_spacer_id, sstart, send, target_seq_sha256)
    dedup_rows: list[dict[str, Any]] = []
    seen_dedup_keys: set[tuple[str, int, int, str]] = set()
    observed_pools: set[str] = set()

    for _, r in m.iterrows():
        sid = str(r["unique_spacer_id"])
        sstart = int(r["sstart"])
        send = int(r["send"])
        sha = str(r.get("target_seq_sha256", ""))
        pool = str(r.get("pool", ""))
        observed_pools.add(pool)

        k = (sid, sstart, send, sha)
        if k in seen_dedup_keys:
            continue
        seen_dedup_keys.add(k)

        # Explicit join on metadata
        if sid not in spc_meta:
            raise ValueError(f"CRITICAL: Spacer {sid} from alignments table not found in verified stage 1 spacers!")

        meta = spc_meta[sid]
        rep_seq = str(meta.get("representative_spacer_seq", "")).strip()
        rc_seq = str(meta.get("canonical_spacer_seq_rc", "")).strip()
        spc_len = int(meta.get("spacer_len", 0))
        contig_id = str(meta.get("primary_namespaced_contig_id", "")).strip()
        arr_start = int(meta.get("primary_array_start", 0))
        arr_end = int(meta.get("primary_array_end", 0))

        # Explicit validation: non-empty sequence and length > 0
        if not rep_seq or rep_seq == "nan" or spc_len <= 0:
            raise ValueError(f"CRITICAL: Spacer {sid} has invalid/empty sequence '{rep_seq}' or length {spc_len}!")

        lo = min(sstart, send)
        hi = max(sstart, send)

        row = {
            "unique_spacer_id": sid,
            "representative_spacer_seq": rep_seq,
            "canonical_spacer_seq_rc": rc_seq,
            "spacer_len": spc_len,
            "primary_namespaced_contig_id": contig_id,
            "primary_array_start": arr_start,
            "primary_array_end": arr_end,
            "protospacer_start": sstart,
            "protospacer_end": send,
            "protospacer_span_bp": hi - lo + 1,
            "effective_mismatches": int(r["effective_mismatches"]),
            "blast_pident": float(r["blast_pident"]),
            "evalue": float(r["evalue"]),
            "bitscore": float(r["bitscore"]),
            "match_tier": str(r["match_tier"]),
            "passes_strict_spacer_gate": bool(r["passes_strict_spacer_gate"]),
            "passes_exploratory_spacer_gate": bool(r.get("passes_exploratory_spacer_gate", False)),
            "btop": str(r["btop"]),
            "target_seq_sha256": sha,
            "pools_observed": sorted(observed_pools),
        }
        dedup_rows.append(row)

    out_df = pd.DataFrame(dedup_rows)
    if out_tsv is not None:
        out_df.to_csv(out_tsv, sep="\t", index=False)

    # Evaluate unique sequences and non-overlapping protospacer loci
    unique_rep_seqs = sorted(set(out_df["representative_spacer_seq"]))
    unique_rc_seqs = sorted(set(out_df["canonical_spacer_seq_rc"]))

    protospacer_intervals = [
        (min(r["protospacer_start"], r["protospacer_end"]), max(r["protospacer_start"], r["protospacer_end"]))
        for _, r in out_df.iterrows()
    ]
    sorted_intervals = sorted(protospacer_intervals)
    non_overlapping_loci = 0
    prev_hi = -1
    for lo, hi in sorted_intervals:
        if lo > prev_hi:
            non_overlapping_loci += 1
            prev_hi = hi

    strict_rows = out_df[out_df["passes_strict_spacer_gate"] == True]
    exact_0mm = out_df[out_df["effective_mismatches"] == 0]

    summary = {
        "host_mag": "41PRr2__single_SemiBin_236",
        "host_phylum": "p__Thermoplasmatota",
        "host_order": "o__JACPAO01",
        "target_virus": "V60PRr2::V60PRr2_715976",
        "target_seq_sha256": dedup_rows[0]["target_seq_sha256"] if dedup_rows else "",
        "raw_alignment_records_across_pools": raw_match_count,
        "n_pools_observed": len(observed_pools),
        "deduplicated_alignment_count": len(out_df),
        "n_strict_spacers": len(strict_rows),
        "n_exact_0mm_spacers": len(exact_0mm),
        "n_unique_representative_sequences": len(unique_rep_seqs),
        "n_unique_canonical_rc_sequences": len(unique_rc_seqs),
        "all_spacer_sequences_unique": (len(unique_rep_seqs) == len(out_df)),
        "n_distinct_protospacer_loci": non_overlapping_loci,
        "protospacer_intervals": sorted_intervals,
        "protospacer_min_coord": min(lo for lo, _ in sorted_intervals),
        "protospacer_max_coord": max(hi for _, hi in sorted_intervals),
        "protospacer_span_total_bp": max(hi for _, hi in sorted_intervals) - min(lo for lo, _ in sorted_intervals) + 1,
        "multiple_distinct_loci_confirmed": (non_overlapping_loci >= 3),
        "sample_cooccurrence": {
            "host_sample": "41PRr2",
            "virus_sample": "V60PRr2",
            "same_sample": False,
            "biological_interpretation": "Cross-sample CRISPR immunity: host 41PRr2 possesses spacers against virus V60PRr2",
        },
    }

    if out_json is not None:
        out_json.write_text(json.dumps(summary, indent=2) + "\n")

    return summary


def main() -> int:
    scored_p = RES_DIR / "stage2_all_spacer_alignments_scored.tsv"
    spc_p = RES_DIR / "stage1_verified_archaeal_spacers.tsv"
    out_tsv = RES_DIR / "stage3_thermo_pair_spacers_forensic_audit.tsv"
    out_json = RES_DIR / "stage3_thermo_pair_spacers_forensic_summary.json"

    print("=" * 80)
    print("THERMOPLASMATOTA PAIR FORENSIC SPACER & PROTOSPACER AUDIT (v2)")
    print("=" * 80)

    summary = audit_thermo_spacers(scored_p, spc_p, out_tsv, out_json)

    print(f"Raw Alignments Across Pools: {summary['raw_alignment_records_across_pools']}")
    print(f"Deduplicated Target Alignments: {summary['deduplicated_alignment_count']}")
    print(f"Strict Spacers: {summary['n_strict_spacers']}")
    print(f"Exact 0mm Spacers: {summary['n_exact_0mm_spacers']}")
    print(f"Unique Spacer Sequences: {summary['n_unique_representative_sequences']} / {summary['deduplicated_alignment_count']}")
    print(f"Distinct Non-overlapping Protospacer Loci: {summary['n_distinct_protospacer_loci']}")
    print(f"Protospacer Coordinate Span: {summary['protospacer_min_coord']}..{summary['protospacer_max_coord']} ({summary['protospacer_span_total_bp']} bp)")
    print("=" * 80)
    return 0


if __name__ == "__main__":
    sys.exit(main())
