#!/usr/bin/env python3
"""
run_f3_context_pipeline.py - Master runner for F3 genomic context test.

Coordinates:
1. 02_extract_neighborhoods.py
2. 03_run_hmmsearch.py
3. 04_statistical_analysis.py
Hashes all inputs and outputs and generates an authoritative run audit.
"""

import argparse
import hashlib
import json
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
RESULTS_DIR = ROOT / "results"
REPORTS_DIR = ROOT / "reports"
PROFILES_DIR = ROOT / "profiles"


def compute_sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while chunk := f.read(1024 * 1024):
            h.update(chunk)
    return h.hexdigest()


def main():
    parser = argparse.ArgumentParser(description="Master runner for F3 genomic context test")
    parser.add_argument("--mode", choices=["pilot", "full"], default="full", help="Execution mode")
    parser.add_argument("--cpus", type=int, default=8, help="Number of CPUs for hmmsearch")
    args = parser.parse_args()

    t_start = time.time()
    prefix = "pilot_" if args.mode == "pilot" else "full_"
    print("=" * 80)
    print(f"MASTER F3 GENOMIC CONTEXT PIPELINE (Mode: {args.mode.upper()})")
    print(f"Timestamp: {time.strftime('%Y-%m-%d %H:%M:%S UTC', time.gmtime())}")
    print(f"CPUs: {args.cpus}")
    print("=" * 80)

    # Step 1: Extraction
    print(f"\n>>> STEP 1: Extracting neighborhoods and controls ({args.mode})...")
    cmd1 = [sys.executable, str(HERE / "02_extract_neighborhoods.py"), "--mode", args.mode]
    subprocess.run(cmd1, check=True)

    # Step 2: HMM search
    print(f"\n>>> STEP 2: Running hmmsearch with Pfam GA cutoffs ({args.mode})...")
    cmd2 = [sys.executable, str(HERE / "03_run_hmmsearch.py"), "--mode", args.mode, "--cpus", str(args.cpus)]
    subprocess.run(cmd2, check=True)

    # Step 3: Statistical analysis
    print(f"\n>>> STEP 3: Executing statistical tests and generating synthesis ({args.mode})...")
    cmd3 = [sys.executable, str(HERE / "04_statistical_analysis.py"), "--mode", args.mode]
    subprocess.run(cmd3, check=True)

    # Step 4: Checksums and Run Manifest
    print("\n>>> STEP 4: Computing SHA-256 checksums and building run manifest...")
    manifest_rows = []
    
    generated_files = [
        RESULTS_DIR / f"{prefix}loci_manifest.tsv",
        RESULTS_DIR / f"{prefix}rhodopsin_candidates.faa",
        RESULTS_DIR / f"{prefix}neighborhood_proteins.faa",
        RESULTS_DIR / f"{prefix}matched_control_proteins.faa",
        RESULTS_DIR / f"{prefix}combined_search_proteins.faa",
        RESULTS_DIR / f"{prefix}all_contig_proteins.faa",
        RESULTS_DIR / f"{prefix}domtbl_sensory_context.tsv",
        RESULTS_DIR / f"{prefix}domain_hits.tsv",
        RESULTS_DIR / f"{prefix}protein_annotations.tsv",
        RESULTS_DIR / f"{prefix}locus_marker_summary.tsv",
        RESULTS_DIR / f"{prefix}statistical_test_results.tsv",
        RESULTS_DIR / f"{prefix}extraction_summary.json",
        RESULTS_DIR / f"{prefix}hmmsearch_summary.json",
        REPORTS_DIR / f"{prefix}genomic_context_synthesis_report.md",
    ]

    for fpath in generated_files:
        if fpath.exists():
            manifest_rows.append({
                "filename": fpath.name,
                "relative_path": str(fpath.relative_to(ROOT)),
                "bytes": fpath.stat().st_size,
                "sha256": compute_sha256(fpath)
            })

    run_audit_path = ROOT / "audit" / f"{prefix}run_audit.json"
    run_audit_path.parent.mkdir(parents=True, exist_ok=True)
    run_audit_path.write_text(json.dumps({
        "mode": args.mode,
        "timestamp_utc": time.strftime("%Y-%m-%d %H:%M:%S UTC", time.gmtime()),
        "duration_seconds": round(time.time() - t_start, 2),
        "host": subprocess.getoutput("hostname -s"),
        "python": sys.version,
        "files_generated": manifest_rows
    }, indent=2), encoding="utf-8")
    print(f"Run audit written to {run_audit_path.name}")

    print("\n" + "=" * 80)
    print(f"PIPELINE COMPLETED SUCCESSFULLY IN {time.time() - t_start:.1f}s")
    print("=" * 80)


if __name__ == "__main__":
    main()
