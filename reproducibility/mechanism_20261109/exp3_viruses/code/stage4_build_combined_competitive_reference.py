#!/usr/bin/env python3
"""Stage 4: Build Full Unified Combined Competitive Reference Database.

Constructs ONE comprehensive competitive reference containing:
1. Canonical 5,438 MAG sequences streamed from:
   /project/k10018/empty-quarter/data/mobilome/tars/refined_bins_all.tar
   strictly validated against branch inputs/genome_manifest.tsv (expected SHA-256).
2. Full competitive viral and mobile background:
   - /project/k10018/empty-quarter/data/mobilome/catalogues/votu_vclust_1.3.1_ani95_af85/votu_reps.fna
   - /project/k10018/empty-quarter/data/mobilome/mobilome_inputs/all_viruses.fna
   - /project/k10018/empty-quarter/data/mobilome/mobilome_inputs/all_plasmids.fna
   - /scratch/project/k10018/work/mobilome_hyp_20261104/e4/derep40/coasm_viruses.fna
   - stage2_raw_outputs_v2/all_matched_primary_targets.fna (ensures exact focal targets present)
   - stage2_raw_outputs_v2/curated_cellular_controls_verified.fna (verified decoys)

Deduplicates exact sequences (and reverse complements) while preserving all ownership aliases.
Explicitly identifies shared host-virus sequences as ambiguous so coverage is not falsely assigned.
Builds Bowtie2 index with cryptographic receipts.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
import tarfile
import time
from pathlib import Path
from typing import Any

import pandas as pd

BRANCH_ROOT = Path(
    os.environ.get(
        "EXP3_ROOT",
        "/scratch/project/k10018/eq-mechanism-20261109/exp3_viruses",
    )
)
IN_DIR = BRANCH_ROOT / "inputs"
RES_DIR = BRANCH_ROOT / "results"
RAW_DIR = RES_DIR / "stage2_raw_outputs_v2"

# Authoritative input paths
P_REFINED_TAR = Path("/project/k10018/empty-quarter/data/mobilome/tars/refined_bins_all.tar")
P_GENOME_MANIFEST = IN_DIR / "genome_manifest.tsv"
P_VOTU_REPS = Path("/project/k10018/empty-quarter/data/mobilome/catalogues/votu_vclust_1.3.1_ani95_af85/votu_reps.fna")
P_ALL_VIR = Path("/project/k10018/empty-quarter/data/mobilome/mobilome_inputs/all_viruses.fna")
P_ALL_PLAS = Path("/project/k10018/empty-quarter/data/mobilome/mobilome_inputs/all_plasmids.fna")
P_COASM_VIR = Path("/scratch/project/k10018/work/mobilome_hyp_20261104/e4/derep40/coasm_viruses.fna")
P_MATCHED_TARGETS = RAW_DIR / "all_matched_primary_targets.fna"
P_CELL_CONTROLS = RAW_DIR / "curated_cellular_controls_verified.fna"

REVCOMP_TAB = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq: str) -> str:
    return seq.translate(REVCOMP_TAB)[::-1]


def stream_sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while chunk := f.read(4 * 1024 * 1024):
            h.update(chunk)
    return h.hexdigest()


def run_cmd(cmd: list[str], desc: str) -> None:
    print(f"[CMD] {desc}: {' '.join(cmd)}", flush=True)
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        print(f"[ERROR] {desc} failed (exit {res.returncode}):\n{res.stderr}", flush=True)
        raise RuntimeError(f"{desc} failed: {res.stderr}")


def parse_fasta_bytes_simple(raw_bytes: bytes) -> list[tuple[str, str]]:
    """Parse FASTA records from in-memory bytes."""
    text = raw_bytes.decode("utf-8", errors="replace")
    records: list[tuple[str, str]] = []
    cur_hdr = ""
    cur_seq: list[str] = []
    for line in text.splitlines():
        if line.startswith(">"):
            if cur_hdr:
                records.append((cur_hdr, "".join(cur_seq).upper()))
            cur_hdr = line[1:].strip().split()[0]
            cur_seq = []
        else:
            cur_seq.append(line.strip())
    if cur_hdr:
        records.append((cur_hdr, "".join(cur_seq).upper()))
    return records


def build_full_competitive_reference(
    refined_tar_path: Path,
    manifest_tsv_path: Path,
    viral_sources: list[tuple[str, Path, str]],
    out_fna: Path,
    crosswalk_tsv: Path,
) -> dict[str, Any]:
    """Stream canonical MAGs from tar, ingest all viral/plasmid sources, deduplicate, and write crosswalk."""
    t0 = time.time()
    out_fna.parent.mkdir(parents=True, exist_ok=True)

    # 1. Load canonical 5,438 manifest
    print(f"[Reference Build] Loading genome manifest from {manifest_tsv_path}...", flush=True)
    m_df = pd.read_csv(manifest_tsv_path, sep="\t")
    required_cols = {"mag", "sha256"}
    if not required_cols.issubset(set(m_df.columns)):
        raise ValueError(f"Manifest missing required columns: {required_cols - set(m_df.columns)}")
    manifest_lookup: dict[str, str] = dict(zip(m_df["mag"].astype(str), m_df["sha256"].astype(str)))
    n_expected_mags = len(manifest_lookup)
    print(f"[Reference Build] Manifest specifies {n_expected_mags} canonical MAGs.", flush=True)

    # Tracking structures
    seen_canonical_sha: dict[str, str] = {}  # can_sha -> representative_header
    seq_category_tracker: dict[str, set[str]] = {}  # can_sha -> set of source_categories (e.g. {"host_mag", "virus_target"})
    crosswalk_rows: list[dict[str, Any]] = []

    n_mags_found = 0
    n_host_contigs = 0
    n_viral_contigs = 0
    n_plasmid_contigs = 0
    n_decoy_contigs = 0
    n_written_contigs = 0
    n_duplicate_seqs = 0

    with open(out_fna, "wt") as out_fh:
        # Step A: Stream canonical 5,438 MAGs from refined_bins_all.tar
        print(f"[Reference Build] Streaming canonical MAGs from {refined_tar_path}...", flush=True)
        if not refined_tar_path.exists():
            raise FileNotFoundError(f"Missing refined_bins_all.tar at {refined_tar_path}")

        found_mags: set[str] = set()
        with tarfile.open(refined_tar_path, "r:") as tf:
            for member in tf:
                if not member.isfile() or not member.name.endswith(".fa"):
                    continue
                bname = Path(member.name).name
                mag_id = bname[:-3]
                if mag_id not in manifest_lookup:
                    continue

                fobj = tf.extractfile(member)
                if fobj is None:
                    continue
                raw_bytes = fobj.read()
                actual_sha = hashlib.sha256(raw_bytes).hexdigest()
                expected_sha = manifest_lookup[mag_id]

                if actual_sha != expected_sha:
                    raise ValueError(
                        f"CRITICAL: SHA-256 mismatch for canonical MAG {mag_id} in {refined_tar_path}! "
                        f"Expected {expected_sha}, got {actual_sha}"
                    )

                found_mags.add(mag_id)
                n_mags_found += 1
                records = parse_fasta_bytes_simple(raw_bytes)

                for c_hdr, s in records:
                    n_host_contigs += 1
                    rc = revcomp(s)
                    can_s = s if s <= rc else rc
                    can_sha = hashlib.sha256(can_s.encode("ascii")).hexdigest()
                    full_hdr = f"HOST__{mag_id}__{c_hdr}"

                    if can_sha not in seq_category_tracker:
                        seq_category_tracker[can_sha] = set()
                    seq_category_tracker[can_sha].add("host_mag")

                    if can_sha in seen_canonical_sha:
                        n_duplicate_seqs += 1
                        rep_h = seen_canonical_sha[can_sha]
                        crosswalk_rows.append({
                            "original_header": full_hdr,
                            "source_category": "host_mag",
                            "entity_id": mag_id,
                            "contig_id": c_hdr,
                            "seq_length": len(s),
                            "canonical_sha256": can_sha,
                            "is_reference_representative": 0,
                            "representative_header": rep_h,
                            "dedup_reason": "duplicate_sequence_of_earlier_contig",
                        })
                    else:
                        seen_canonical_sha[can_sha] = full_hdr
                        out_fh.write(f">{full_hdr}\n{s}\n")
                        n_written_contigs += 1
                        crosswalk_rows.append({
                            "original_header": full_hdr,
                            "source_category": "host_mag",
                            "entity_id": mag_id,
                            "contig_id": c_hdr,
                            "seq_length": len(s),
                            "canonical_sha256": can_sha,
                            "is_reference_representative": 1,
                            "representative_header": full_hdr,
                            "dedup_reason": "primary_representative",
                        })

                if n_mags_found % 500 == 0:
                    print(f"  ...processed {n_mags_found}/{n_expected_mags} canonical MAGs ({n_written_contigs} contigs written)", flush=True)

        if len(found_mags) != n_expected_mags:
            missing = set(manifest_lookup.keys()) - found_mags
            raise ValueError(
                f"CRITICAL: Incomplete canonical MAG cohort in tar! Found {len(found_mags)}/{n_expected_mags}. "
                f"Sample missing: {sorted(list(missing))[:10]}"
            )
        print(f"[Reference Build] Successfully verified and ingested all {n_mags_found} canonical MAGs.", flush=True)

        # Step B: Ingest full competitive viral, plasmid, and decoy sources
        for cat_name, file_path, hdr_prefix in viral_sources:
            if not file_path.exists():
                print(f"[WARN] Source file missing: {file_path}", flush=True)
                continue
            print(f"[Reference Build] Ingesting {cat_name} from {file_path}...", flush=True)
            cur_hdr = ""
            cur_seq: list[str] = []
            source_contig_count = 0

            with open(file_path, "rt", errors="replace") as fh:
                for line in fh:
                    if line.startswith(">"):
                        if cur_hdr:
                            source_contig_count += 1
                            s = "".join(cur_seq).upper()
                            rc = revcomp(s)
                            can_s = s if s <= rc else rc
                            can_sha = hashlib.sha256(can_s.encode("ascii")).hexdigest()
                            full_hdr = f"{hdr_prefix}__{cur_hdr}"

                            if can_sha not in seq_category_tracker:
                                seq_category_tracker[can_sha] = set()
                            seq_category_tracker[can_sha].add(cat_name)

                            if cat_name in {"virus_target", "all_viruses", "votu_reps", "coasm_viruses"}:
                                n_viral_contigs += 1
                            elif cat_name == "all_plasmids":
                                n_plasmid_contigs += 1
                            elif cat_name == "cellular_controls":
                                n_decoy_contigs += 1

                            if can_sha in seen_canonical_sha:
                                n_duplicate_seqs += 1
                                rep_h = seen_canonical_sha[can_sha]
                                crosswalk_rows.append({
                                    "original_header": full_hdr,
                                    "source_category": cat_name,
                                    "entity_id": cur_hdr,
                                    "contig_id": cur_hdr,
                                    "seq_length": len(s),
                                    "canonical_sha256": can_sha,
                                    "is_reference_representative": 0,
                                    "representative_header": rep_h,
                                    "dedup_reason": "duplicate_sequence_of_earlier_contig",
                                })
                            else:
                                seen_canonical_sha[can_sha] = full_hdr
                                out_fh.write(f">{full_hdr}\n{s}\n")
                                n_written_contigs += 1
                                crosswalk_rows.append({
                                    "original_header": full_hdr,
                                    "source_category": cat_name,
                                    "entity_id": cur_hdr,
                                    "contig_id": cur_hdr,
                                    "seq_length": len(s),
                                    "canonical_sha256": can_sha,
                                    "is_reference_representative": 1,
                                    "representative_header": full_hdr,
                                    "dedup_reason": "primary_representative",
                                })
                        cur_hdr = line[1:].strip().split()[0]
                        cur_seq = []
                    else:
                        cur_seq.append(line.strip())
                if cur_hdr:
                    source_contig_count += 1
                    s = "".join(cur_seq).upper()
                    rc = revcomp(s)
                    can_s = s if s <= rc else rc
                    can_sha = hashlib.sha256(can_s.encode("ascii")).hexdigest()
                    full_hdr = f"{hdr_prefix}__{cur_hdr}"

                    if can_sha not in seq_category_tracker:
                        seq_category_tracker[can_sha] = set()
                    seq_category_tracker[can_sha].add(cat_name)

                    if cat_name in {"virus_target", "all_viruses", "votu_reps", "coasm_viruses"}:
                        n_viral_contigs += 1
                    elif cat_name == "all_plasmids":
                        n_plasmid_contigs += 1
                    elif cat_name == "cellular_controls":
                        n_decoy_contigs += 1

                    if can_sha in seen_canonical_sha:
                        n_duplicate_seqs += 1
                        rep_h = seen_canonical_sha[can_sha]
                        crosswalk_rows.append({
                            "original_header": full_hdr,
                            "source_category": cat_name,
                            "entity_id": cur_hdr,
                            "contig_id": cur_hdr,
                            "seq_length": len(s),
                            "canonical_sha256": can_sha,
                            "is_reference_representative": 0,
                            "representative_header": rep_h,
                            "dedup_reason": "duplicate_sequence_of_earlier_contig",
                        })
                    else:
                        seen_canonical_sha[can_sha] = full_hdr
                        out_fh.write(f">{full_hdr}\n{s}\n")
                        n_written_contigs += 1
                        crosswalk_rows.append({
                            "original_header": full_hdr,
                            "source_category": cat_name,
                            "entity_id": cur_hdr,
                            "contig_id": cur_hdr,
                            "seq_length": len(s),
                            "canonical_sha256": can_sha,
                            "is_reference_representative": 1,
                            "representative_header": full_hdr,
                            "dedup_reason": "primary_representative",
                        })
            print(f"  ...ingested {source_contig_count} contigs from {cat_name}", flush=True)

    # Step C: Assign Crosswalk Ambiguity Status
    print("[Reference Build] Assigning crosswalk ambiguity status...", flush=True)
    for r in crosswalk_rows:
        sha = r["canonical_sha256"]
        cats = seq_category_tracker.get(sha, set())
        has_host = "host_mag" in cats
        has_viral = bool(cats & {"virus_target", "all_viruses", "votu_reps", "coasm_viruses"})
        has_plasmid = "all_plasmids" in cats

        if has_host and (has_viral or has_plasmid):
            ambig = "shared_host_viral_ambiguous"
        elif has_host:
            ambig = "unambiguous_host"
        elif has_viral and has_plasmid:
            ambig = "shared_viral_plasmid"
        elif has_viral:
            ambig = "unambiguous_virus"
        elif has_plasmid:
            ambig = "unambiguous_plasmid"
        else:
            ambig = "unambiguous_decoy"
        r["ambiguity_status"] = ambig

    cw_df = pd.DataFrame(crosswalk_rows)
    cw_df.to_csv(crosswalk_tsv, sep="\t", index=False)
    print(f"[Reference Build] Crosswalk saved to {crosswalk_tsv} ({len(cw_df)} records)", flush=True)

    stats = {
        "elapsed_seconds": round(time.time() - t0, 2),
        "n_canonical_mags_verified": n_mags_found,
        "n_host_contigs_raw": n_host_contigs,
        "n_viral_contigs_raw": n_viral_contigs,
        "n_plasmid_contigs_raw": n_plasmid_contigs,
        "n_decoy_contigs_raw": n_decoy_contigs,
        "n_total_input_contigs": len(crosswalk_rows),
        "n_duplicate_sequences_deduplicated": n_duplicate_seqs,
        "n_unique_representative_contigs_written": n_written_contigs,
        "n_ambiguous_host_viral_sequences": int((cw_df["ambiguity_status"] == "shared_host_viral_ambiguous").sum()),
        "reference_fna_bytes": out_fna.stat().st_size,
    }
    return stats


def main() -> int:
    t0 = time.time()
    parser = argparse.ArgumentParser(description="Build Full Unified Competitive Reference")
    parser.add_argument("--threads", type=int, default=64, help="Threads for bowtie2-build")
    args = parser.parse_args()

    ref_dir = RES_DIR / "stage4_competitive_reference"
    ref_dir.mkdir(parents=True, exist_ok=True)
    out_fna = ref_dir / "unified_competitive_reference.fna"
    cw_tsv = ref_dir / "stage4_reference_sequence_crosswalk.tsv"
    bt2_prefix = ref_dir / "unified_competitive_reference"

    viral_sources = [
        ("virus_target", P_MATCHED_TARGETS, "VIRUS_TARGET"),
        ("votu_reps", P_VOTU_REPS, "VIRUS_VOTU"),
        ("all_viruses", P_ALL_VIR, "VIRUS_ALL"),
        ("coasm_viruses", P_COASM_VIR, "VIRUS_COASM"),
        ("all_plasmids", P_ALL_PLAS, "PLASMID"),
        ("cellular_controls", P_CELL_CONTROLS, "DECOY"),
    ]

    print("=" * 80)
    print("STAGE 4: BUILD FULL UNIFIED COMPETITIVE REFERENCE")
    print("Scope: 5,438 Canonical MAGs + Full Viral/Plasmid Catalog + Decoys")
    print("=" * 80)

    stats = build_full_competitive_reference(
        refined_tar_path=P_REFINED_TAR,
        manifest_tsv_path=P_GENOME_MANIFEST,
        viral_sources=viral_sources,
        out_fna=out_fna,
        crosswalk_tsv=cw_tsv,
    )
    print("[Reference Build] Compilation summary:\n", json.dumps(stats, indent=2))

    # Bowtie2 build
    print(f"[Reference Build] Building Bowtie2 index using {args.threads} threads...", flush=True)
    bt2_cmd = [
        "bowtie2-build",
        "--threads", str(args.threads),
        str(out_fna),
        str(bt2_prefix),
    ]
    run_cmd(bt2_cmd, "Bowtie2 build unified competitive reference index")

    # Cryptographic index receipt
    idx_files = sorted(ref_dir.glob("unified_competitive_reference.*.bt2*"))
    idx_receipt = {f.name: {"size": f.stat().st_size, "sha256": stream_sha256(f)} for f in idx_files}

    receipt = {
        "status": "completed",
        "stage": "stage4_build_full_unified_competitive_reference",
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "elapsed_seconds": round(time.time() - t0, 2),
        "build_stats": stats,
        "index_prefix": str(bt2_prefix),
        "index_files": idx_receipt,
    }
    receipt_p = ref_dir / "unified_competitive_reference.receipt.json"
    receipt_p.write_text(json.dumps(receipt, indent=2) + "\n")
    (ref_dir / "REFERENCE_BUILD_DONE").touch()

    print("[Reference Build] Finished successfully. Receipt written to", receipt_p)
    return 0


if __name__ == "__main__":
    sys.exit(main())
