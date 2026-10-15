#!/usr/bin/env python3
"""
03_run_hmmsearch.py - Run hmmsearch with Pfam GA cutoffs and annotate neighbor/control proteins.

Features:
1. Runs hmmsearch with --cut_ga using the source-verified f3_sensory_context.hmm database.
2. Parses domtblout strictly for Gathering-threshold compliant domain hits.
3. Classifies hits into sensory and signaling categories (TCS Kinase, TCS Regulator, Transducer/Linker, Sensor PAS/GAF, Second Messenger).
4. Emits detailed domain_hits.tsv and protein_annotations.tsv.
"""

import argparse
import csv
import hashlib
import json
import os
import subprocess
import sys
import time
from pathlib import Path
from typing import Dict, List, Set

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
RESULTS_DIR = ROOT / "results"
PROFILES_DIR = ROOT / "profiles"

# Category mapping for models
CATEGORY_MAP = {
    "HisKA": "TCS_Kinase",
    "HATPase_c": "TCS_Kinase",
    "HisKA_2": "TCS_Kinase",
    "HisKA_3": "TCS_Kinase",
    "Response_reg": "TCS_Regulator",
    "Trans_reg_C": "TCS_Regulator",
    "MCPsignal": "Transducer_Linker",
    "HAMP": "Transducer_Linker",
    "CheW": "Transducer_Linker",
    "PAS": "Sensor_Input",
    "PAS_3": "Sensor_Input",
    "PAS_4": "Sensor_Input",
    "PAS_9": "Sensor_Input",
    "GAF": "Sensor_Input",
    "GAF_2": "Sensor_Input",
    "GGDEF": "Second_Messenger",
    "EAL": "Second_Messenger",
    "Bac_rhodopsin": "Rhodopsin_Control",
    "Heliorhodopsin": "Rhodopsin_Control",
}

SENSORY_CATEGORIES = {
    "TCS_Kinase",
    "TCS_Regulator",
    "Transducer_Linker",
    "Sensor_Input",
    "Second_Messenger",
}


def compute_sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while chunk := f.read(1024 * 1024):
            h.update(chunk)
    return h.hexdigest()


def parse_domtblout(domtbl_path: Path) -> List[dict]:
    """Parse HMMER domtblout format."""
    hits = []
    with open(domtbl_path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            target_name = parts[0]
            target_len = int(parts[2])
            query_name = parts[3]
            query_acc = parts[4]
            query_len = int(parts[5])
            full_evalue = float(parts[6])
            full_score = float(parts[7])
            domain_idx = int(parts[9])
            domain_total = int(parts[10])
            c_evalue = float(parts[11])
            i_evalue = float(parts[12])
            domain_score = float(parts[13])
            hmm_from = int(parts[15])
            hmm_to = int(parts[16])
            ali_from = int(parts[17])
            ali_to = int(parts[18])
            
            category = CATEGORY_MAP.get(query_name, "Other")
            
            hits.append({
                "protein_id": target_name,
                "protein_len": target_len,
                "model_name": query_name,
                "model_acc": query_acc,
                "category": category,
                "domain_score": domain_score,
                "domain_ievalue": i_evalue,
                "full_score": full_score,
                "full_evalue": full_evalue,
                "domain_idx": domain_idx,
                "domain_total": domain_total,
                "hmm_from": hmm_from,
                "hmm_to": hmm_to,
                "ali_from": ali_from,
                "ali_to": ali_to,
            })
    return hits


def main():
    parser = argparse.ArgumentParser(description="Run hmmsearch for sensory/signaling context markers")
    parser.add_argument("--mode", choices=["pilot", "full"], default="full", help="Execution mode: pilot or full")
    parser.add_argument("--hmm-profile", type=Path, default=PROFILES_DIR / "f3_sensory_context.hmm", help="HMM database file")
    parser.add_argument("--cpus", type=int, default=8, help="Number of CPU threads for hmmsearch")
    args = parser.parse_args()

    t0 = time.time()
    prefix = "pilot_" if args.mode == "pilot" else "full_"
    target_faa = RESULTS_DIR / f"{prefix}combined_search_proteins.faa"
    domtbl_out = RESULTS_DIR / f"{prefix}domtbl_sensory_context.tsv"
    hmm_log = RESULTS_DIR / f"{prefix}hmmsearch.log"

    print("=" * 80)
    print(f"RUNNING HMMSEARCH FOR F3 SENSORY CONTEXT (Mode: {args.mode.upper()})")
    print(f"HMM Profile: {args.hmm_profile}")
    print(f"Target FASTA: {target_faa}")
    print(f"CPUs: {args.cpus}")
    print("=" * 80)

    assert args.hmm_profile.exists(), f"Missing HMM profile: {args.hmm_profile}"
    assert target_faa.exists(), f"Missing target FASTA: {target_faa}"

    # Run hmmsearch with --cut_ga
    cmd = [
        "hmmsearch",
        "--cut_ga",
        "--cpu", str(args.cpus),
        "--domtblout", str(domtbl_out),
        "-o", str(hmm_log),
        str(args.hmm_profile),
        str(target_faa)
    ]
    print(f"Executing: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)
    print(f"hmmsearch completed. Output written to {domtbl_out.name}")

    # Parse domtblout
    hits = parse_domtblout(domtbl_out)
    print(f"\nParsed {len(hits)} domain hits passing Pfam GA cutoffs.")

    # Write domain hits TSV
    domain_hits_tsv = RESULTS_DIR / f"{prefix}domain_hits.tsv"
    if hits:
        fieldnames = list(hits[0].keys())
        with open(domain_hits_tsv, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            writer.writerows(hits)
    else:
        with open(domain_hits_tsv, "w", newline="", encoding="utf-8") as f:
            f.write("protein_id\tmodel_name\tcategory\tdomain_score\tdomain_ievalue\n")
    print(f"Wrote domain hits TSV: {domain_hits_tsv.name}")

    # Aggregate hits by protein
    protein_hits: Dict[str, List[dict]] = {}
    for h in hits:
        protein_hits.setdefault(h["protein_id"], []).append(h)

    # Read target fasta to get full list of scanned proteins
    all_scanned_ids = []
    with open(target_faa, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith(">"):
                pid = line[1:].split()[0]
                all_scanned_ids.append(pid)

    protein_annotations = []
    for pid in all_scanned_ids:
        pHits = protein_hits.get(pid, [])
        models_hit = sorted(list({h["model_name"] for h in pHits}))
        cats_hit = sorted(list({h["category"] for h in pHits}))
        
        has_tcs_kinase = int("TCS_Kinase" in cats_hit)
        has_tcs_regulator = int("TCS_Regulator" in cats_hit)
        has_transducer_linker = int("Transducer_Linker" in cats_hit)
        has_sensor_input = int("Sensor_Input" in cats_hit)
        has_second_messenger = int("Second_Messenger" in cats_hit)
        has_rhodopsin_control = int("Rhodopsin_Control" in cats_hit)
        
        is_sensory_signaling = int(any(c in SENSORY_CATEGORIES for c in cats_hit))
        
        protein_annotations.append({
            "protein_id": pid,
            "has_hit": int(len(pHits) > 0),
            "n_domains": len(pHits),
            "models_hit": ";".join(models_hit),
            "categories_hit": ";".join(cats_hit),
            "is_sensory_signaling": is_sensory_signaling,
            "has_tcs_kinase": has_tcs_kinase,
            "has_tcs_regulator": has_tcs_regulator,
            "has_transducer_linker": has_transducer_linker,
            "has_sensor_input": has_sensor_input,
            "has_second_messenger": has_second_messenger,
            "has_rhodopsin_control": has_rhodopsin_control,
            "best_domain_score": max((h["domain_score"] for h in pHits), default=0.0),
            "best_model_name": min(pHits, key=lambda x: x["domain_ievalue"])["model_name"] if pHits else ""
        })

    anno_tsv = RESULTS_DIR / f"{prefix}protein_annotations.tsv"
    with open(anno_tsv, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(protein_annotations[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(protein_annotations)
    print(f"Wrote protein annotations: {anno_tsv.name} ({len(protein_annotations)} proteins)")

    # Print summary counts
    n_signaling = sum(1 for p in protein_annotations if p["is_sensory_signaling"])
    n_kinase = sum(1 for p in protein_annotations if p["has_tcs_kinase"])
    n_reg = sum(1 for p in protein_annotations if p["has_tcs_regulator"])
    n_transducer = sum(1 for p in protein_annotations if p["has_transducer_linker"])
    n_sensor = sum(1 for p in protein_annotations if p["has_sensor_input"])
    n_sec_msg = sum(1 for p in protein_annotations if p["has_second_messenger"])
    n_rho = sum(1 for p in protein_annotations if p["has_rhodopsin_control"])

    print("\nAnnotation Breakdown across all scanned proteins:")
    print(f"  Total scanned: {len(protein_annotations)}")
    print(f"  Any sensory/signaling marker: {n_signaling}")
    print(f"  TCS Kinase: {n_kinase}")
    print(f"  TCS Regulator: {n_reg}")
    print(f"  Transducer/Linker (MCP/HAMP/CheW): {n_transducer}")
    print(f"  Sensor Input (PAS/GAF): {n_sensor}")
    print(f"  Second Messenger (GGDEF/EAL): {n_sec_msg}")
    print(f"  Rhodopsin control hits: {n_rho}")

    # Audit JSON
    audit_path = RESULTS_DIR / f"{prefix}hmmsearch_summary.json"
    audit_path.write_text(json.dumps({
        "mode": args.mode,
        "timestamp_utc": time.strftime("%Y-%m-%d %H:%M:%S UTC", time.gmtime()),
        "hmm_database": str(args.hmm_profile),
        "target_fasta": str(target_faa),
        "target_fasta_sha256": compute_sha256(target_faa),
        "total_domain_hits": len(hits),
        "total_scanned_proteins": len(protein_annotations),
        "annotated_proteins": {
            "sensory_signaling": n_signaling,
            "tcs_kinase": n_kinase,
            "tcs_regulator": n_reg,
            "transducer_linker": n_transducer,
            "sensor_input": n_sensor,
            "second_messenger": n_sec_msg,
            "rhodopsin_control": n_rho
        }
    }, indent=2), encoding="utf-8")
    print(f"\nDone in {time.time() - t0:.1f}s.")


if __name__ == "__main__":
    main()
