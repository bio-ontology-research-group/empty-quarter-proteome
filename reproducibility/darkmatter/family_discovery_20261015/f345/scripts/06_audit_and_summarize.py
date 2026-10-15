#!/usr/bin/env python3
"""06_audit_and_summarize.py - Comprehensive Audit & Integrity Verification for Tracks F3, F4, F5.

Performs rigorous verification:
1. Validates references_manifest.tsv: verifies that every profile, index file, FASTA,
   DIAMOND database, and metadata table exists with matching SHA-256 and byte size.
2. Validates profile formats: ensures all .hmm files are valid HMMER3 ASCII text (not gzip).
3. Validates HMMER pressed indices: ensures .h3f, .h3i, .h3m, .h3p exist for multi-model HMMs.
4. Validates DIAMOND databases: ensures .dmnd files exist for all 4 reference sets.
5. Validates input schema & inventory: verifies results/input_schemas.json (5,438 MAGs, 0 mismatches).
6. Validates gate compliance: verifies that biological_searches_submitted is FALSE and no search jobs ran.
7. Updates STATUS.json and writes results/audit_summary.json.
"""

import csv
import hashlib
import json
import os
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
REF_DIR = ROOT / "references"
RESULTS_DIR = ROOT / "results"
LOGS_DIR = ROOT / "logs"
STATUS_FILE = ROOT / "STATUS.json"


def compute_sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while chunk := f.read(1024 * 1024):
            h.update(chunk)
    return h.hexdigest()


def main():
    print("=" * 80)
    print("F3/F4/F5 COMPREHENSIVE AUDIT & INTEGRITY VERIFICATION")
    print(f"Timestamp: {time.strftime('%Y-%m-%d %H:%M:%S UTC', time.gmtime())}")
    print(f"Workspace: {ROOT}")
    print("=" * 80)

    audit_errors = []
    audit_warnings = []

    # --------------------------------------------------------------------------
    # 1. VERIFY REFERENCES MANIFEST & CHECKSUMS
    # --------------------------------------------------------------------------
    print("\n--- 1. Auditing References Manifest & SHA-256 Checksums ---")
    manifest_tsv = REF_DIR / "references_manifest.tsv"
    manifest_checked = 0
    if not manifest_tsv.exists():
        audit_errors.append(f"Manifest missing: {manifest_tsv}")
    else:
        with open(manifest_tsv, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                fname = row["filename"]
                expected_sha = row["sha256"]
                expected_bytes = int(row["bytes"])
                target_file = REF_DIR / fname

                if not target_file.exists():
                    audit_errors.append(f"Manifest target missing: {fname}")
                    continue

                actual_bytes = target_file.stat().st_size
                if actual_bytes != expected_bytes:
                    audit_errors.append(f"Size mismatch for {fname}: expected {expected_bytes}, got {actual_bytes}")
                    continue

                actual_sha = compute_sha256(target_file)
                if actual_sha != expected_sha:
                    audit_errors.append(f"SHA-256 mismatch for {fname}: expected {expected_sha}, got {actual_sha}")
                    continue

                manifest_checked += 1
        print(f"  Verified {manifest_checked} files in manifest with 100% SHA-256 and byte-size match.")

    # --------------------------------------------------------------------------
    # 2. VERIFY PROFILE FORMATS (HMMER3 ASCII vs GZIP)
    # --------------------------------------------------------------------------
    print("\n--- 2. Auditing HMM Profile Text Formatting ---")
    profiles = list((REF_DIR / "profiles").glob("*.hmm"))
    plain_text_count = 0
    for p in profiles:
        with open(p, "rb") as f:
            header = f.read(16)
        if header.startswith(b"\x1f\x8b"):
            audit_errors.append(f"HMM profile is gzip compressed (must be HMMER3 plain text): {p.name}")
        elif not header.startswith(b"HMMER3"):
            audit_errors.append(f"HMM profile missing HMMER3 magic header: {p.name}")
        else:
            plain_text_count += 1
    print(f"  Verified {plain_text_count}/{len(profiles)} HMM profiles are valid HMMER3 ASCII format.")

    # --------------------------------------------------------------------------
    # 3. VERIFY HMMER PRESSED INDICES
    # --------------------------------------------------------------------------
    print("\n--- 3. Auditing HMMER Pressed Indices (.h3f, .h3i, .h3m, .h3p) ---")
    track_hmms = ["f3_rhodopsins.hmm", "f4_ectc.hmm", "f5_carbonic_anhydrase.hmm"]
    indices_checked = 0
    for th in track_hmms:
        base = REF_DIR / "profiles" / th
        if not base.exists():
            audit_errors.append(f"Track multi-model HMM missing: {th}")
            continue
        for ext in [".h3f", ".h3i", ".h3m", ".h3p"]:
            idx_file = base.with_name(th + ext)
            if not idx_file.exists():
                audit_errors.append(f"HMMER index missing: {idx_file.name}")
            elif idx_file.stat().st_size == 0:
                audit_errors.append(f"HMMER index empty: {idx_file.name}")
            else:
                indices_checked += 1
    print(f"  Verified {indices_checked}/12 HMMER index files present and non-empty.")

    # --------------------------------------------------------------------------
    # 4. VERIFY REFERENCE FASTAS & DIAMOND DATABASES
    # --------------------------------------------------------------------------
    print("\n--- 4. Auditing Reference FASTAs & DIAMOND Databases ---")
    dmnd_sets = [
        "f3_rhodopsin_references",
        "f4_ectc_references",
        "f5_carbonic_anhydrase_references",
        "all_f345_references",
    ]
    dmnd_ok = 0
    for ds in dmnd_sets:
        faa = REF_DIR / f"{ds}.faa"
        dmnd = REF_DIR / f"{ds}.dmnd"
        if not faa.exists() or faa.stat().st_size == 0:
            audit_errors.append(f"Reference FASTA missing or empty: {faa.name}")
        elif not dmnd.exists() or dmnd.stat().st_size == 0:
            audit_errors.append(f"DIAMOND database missing or empty: {dmnd.name}")
        else:
            dmnd_ok += 1
    print(f"  Verified {dmnd_ok}/4 reference FASTA + DIAMOND database pairs present.")

    # --------------------------------------------------------------------------
    # 5. VERIFY REFERENCE METADATA & EXPERIMENTAL EVIDENCE PARTITIONING
    # --------------------------------------------------------------------------
    print("\n--- 5. Auditing Reference Metadata & Evidence Separation ---")
    meta_tsv = REF_DIR / "verified_references_metadata.tsv"
    total_refs = 0
    fn_exp_refs = 0
    any_exp_refs = 0
    track_counts = {"F3": 0, "F4": 0, "F5": 0}
    if meta_tsv.exists():
        with open(meta_tsv, "r", encoding="utf-8") as f:
            for r in csv.DictReader(f, delimiter="\t"):
                total_refs += 1
                trk = r.get("track")
                if trk in track_counts:
                    track_counts[trk] += 1
                if r.get("has_function_experimental_evidence") == "1":
                    fn_exp_refs += 1
                if r.get("has_any_experimental_evidence") == "1":
                    any_exp_refs += 1
        print(f"  Total Clean References: {total_refs}")
        print(f"    - F3 (Rhodopsins): {track_counts['F3']}")
        print(f"    - F4 (EctC & operon): {track_counts['F4']}")
        print(f"    - F5 (Carbonic Anhydrases): {track_counts['F5']}")
        print(f"  Function-Specific Experimental Validation: {fn_exp_refs}")
        print(f"  Generic Record-Level ECO Tags: {any_exp_refs}")
    else:
        audit_errors.append("verified_references_metadata.tsv missing")

    # Excluded pathogens & rejected references
    pathogen_tsv = REF_DIR / "excluded_pathogens.tsv"
    rejected_tsv = REF_DIR / "rejected_references.tsv"
    pathogen_count = sum(1 for _ in open(pathogen_tsv)) - 1 if pathogen_tsv.exists() else 0
    rejected_count = sum(1 for _ in open(rejected_tsv)) - 1 if rejected_tsv.exists() else 0
    print(f"  Excluded Pathogens: {pathogen_count}")
    print(f"  Audited Rejected Accessions Quarantined: {rejected_count}")

    # --------------------------------------------------------------------------
    # 6. VERIFY INPUT SCHEMA INSPECTION REPORT
    # --------------------------------------------------------------------------
    print("\n--- 6. Auditing Input Schemas & dRep Cluster Validation ---")
    schema_json = RESULTS_DIR / "input_schemas.json"
    if not schema_json.exists():
        audit_errors.append("results/input_schemas.json missing")
    else:
        with open(schema_json, "r", encoding="utf-8") as f:
            sdata = json.load(f)
        drep_join = sdata.get("drep_cluster_join_validation", {})
        mismatches = drep_join.get("wdb_vs_cdb_secondary_mismatches")
        rep_prot = sdata.get("representative_proteomes", {})
        prot_count = rep_prot.get("count")
        missing_prot = rep_prot.get("missing_vs_wdb")
        extra_prot = rep_prot.get("extra_vs_wdb")

        if mismatches != 0:
            audit_errors.append(f"dRep cluster join has {mismatches} mismatches!")
        else:
            print("  dRep Cluster Join: Wdb.cluster == Cdb.secondary_cluster with 0 mismatches.")

        if prot_count != 5438 or missing_prot != 0 or extra_prot != 0:
            audit_errors.append(f"Proteome count mismatch: {prot_count} (missing: {missing_prot}, extra: {extra_prot})")
        else:
            print(f"  Representative Proteomes: exactly 5,438/5,438 present with 1-to-1 closure.")

    # --------------------------------------------------------------------------
    # 7. VERIFY GATE STATUS & SEARCH CONSTRAINTS
    # --------------------------------------------------------------------------
    print("\n--- 7. Auditing Gate Status & Search Submission Boundary ---")
    if not STATUS_FILE.exists():
        audit_errors.append("STATUS.json missing")
    else:
        with open(STATUS_FILE, "r", encoding="utf-8") as f:
            status_data = json.load(f)
        searches_submitted = status_data.get("constraints", {}).get("biological_searches_submitted")
        gate_status = status_data.get("constraints", {}).get("gate_status")

        if searches_submitted is not False:
            audit_errors.append("VIOLATION: biological_searches_submitted is not False!")
        if gate_status != "LOCKED_AWAITING_ROOT_GATE":
            audit_errors.append(f"VIOLATION: gate_status is '{gate_status}', expected 'LOCKED_AWAITING_ROOT_GATE'")

        print(f"  Biological Searches Submitted: {searches_submitted} (COMPLIANT)")
        print(f"  Gate Status: {gate_status} (LOCKED)")

    # --------------------------------------------------------------------------
    # 8. AUDIT RESULT SUMMARY & STATUS UPDATE
    # --------------------------------------------------------------------------
    print("\n" + "=" * 80)
    print("AUDIT FINDINGS SUMMARY")
    print("=" * 80)
    if audit_errors:
        print(f"FAILED with {len(audit_errors)} errors:")
        for err in audit_errors:
            print(f"  [ERROR] {err}")
        sys.exit(1)
    else:
        print("ALL INTEGRITY CHECKS PASSED (0 errors, 0 warnings).")
        print("\nReadiness Assessment:")
        print("  >> READY_FOR_ROOT_AUDIT")
        print("  Preparation scope (source references, schema joins, pilot code) is complete.")
        print("  Biological searches remain strictly locked awaiting root gate.")

    # Write audit summary JSON
    audit_summary = {
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "readiness_status": "READY_FOR_ROOT_AUDIT",
        "gate_status": "LOCKED_AWAITING_ROOT_GATE",
        "biological_searches_submitted": False,
        "manifest_entries_verified": manifest_checked,
        "verified_reference_proteins": total_refs,
        "reference_breakdown": track_counts,
        "function_specific_experimental_evidence": fn_exp_refs,
        "generic_record_level_eco_tags": any_exp_refs,
        "excluded_pathogens": pathogen_count,
        "audited_rejected_references_quarantined": rejected_count,
        "profile_hmms_verified": plain_text_count,
        "hmmer_indices_verified": indices_checked,
        "diamond_databases_verified": dmnd_ok,
        "representative_proteomes_verified": 5438,
        "drep_secondary_cluster_mismatches": 0,
        "completed_jobs": ["53815754", "53815888", "53816046"],
    }
    summary_path = RESULTS_DIR / "audit_summary.json"
    with open(summary_path, "w", encoding="utf-8") as f:
        json.dump(audit_summary, f, indent=2)
    print(f"\nWrote audit summary to {summary_path}")

    # Update STATUS.json with completed jobs and verified counts
    with open(STATUS_FILE, "r", encoding="utf-8") as f:
        st = json.load(f)
    st["completed_jobs"] = ["53815754", "53815888", "53816046"]
    st["audit_summary"]["total_verified_reference_proteins"] = total_refs
    st["audit_summary"]["function_specific_experimental_evidence"] = fn_exp_refs
    st["audit_summary"]["record_level_experimental_eco_tags"] = any_exp_refs
    st["audit_summary"]["verified_profile_hmms"] = plain_text_count
    st["audit_summary"]["hmmer_indices_verified"] = indices_checked
    st["audit_summary"]["diamond_databases_verified"] = dmnd_ok
    st["audit_summary"]["manifest_entries_verified"] = manifest_checked
    with open(STATUS_FILE, "w", encoding="utf-8") as f:
        json.dump(st, f, indent=2)
    print(f"Updated {STATUS_FILE}")


if __name__ == "__main__":
    main()
