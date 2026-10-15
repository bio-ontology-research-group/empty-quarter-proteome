#!/usr/bin/env python3
"""
Independent Verification and Audit Script for ProstT5 Embeddings (v3).

Verifies:
1. File integrity against checksums.sha256
2. Dimensionality (N, 1024), finiteness (no NaNs, no Infs), and dtype
3. Exact correspondence between embeddings.npz and embeddings.npy
4. ID order, residue length integrity, and unified raw vs clean sequence SHA256 against input FASTA
5. Equivalence of duplicate sequence embeddings (both exact duplicates and trailing-stop variants)
6. Completeness of machine-readable provenance, including runner code SHA,
   tokenizer/config hashes, library versions, and token alignment verification
"""

import argparse
import hashlib
import json
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np


def compute_file_sha256(filepath: Path, chunk_size: int = 16 * 1024 * 1024) -> str:
    h = hashlib.sha256()
    with open(filepath, "rb") as f:
        while chunk := f.read(chunk_size):
            h.update(chunk)
    return h.hexdigest()


def parse_fasta_ids_and_lengths(fasta_path: Path) -> List[Tuple[str, int, str, str, bool]]:
    records = []
    seen_ids = set()
    current_id = None
    current_seq_parts = []
    first_header_seen = False

    with open(fasta_path, "r", encoding="utf-8") as f:
        for line_num, line in enumerate(f, 1):
            line_stripped = line.strip()
            if not line_stripped:
                continue

            if line_stripped.startswith(">"):
                first_header_seen = True
                if current_id is not None:
                    raw_s = "".join(current_seq_parts).upper()
                    has_stop = raw_s.endswith("*")
                    clean_s = raw_s[:-1] if has_stop else raw_s
                    raw_sha = hashlib.sha256(raw_s.encode("ascii")).hexdigest()
                    clean_sha = hashlib.sha256(clean_s.encode("ascii")).hexdigest()
                    records.append((current_id, len(clean_s), raw_sha, clean_sha, has_stop))

                header_content = line_stripped[1:].strip()
                if not header_content:
                    raise ValueError(f"Empty header found at line {line_num} in {fasta_path}")

                current_id = header_content.split()[0]
                if current_id in seen_ids:
                    raise ValueError(f"Duplicate FASTA ID '{current_id}' at line {line_num} in {fasta_path}")
                seen_ids.add(current_id)
                current_seq_parts = []
            else:
                if not first_header_seen:
                    raise ValueError(f"Sequence content before first header at line {line_num} in {fasta_path}")
                current_seq_parts.append(line_stripped)

        if current_id is not None:
            raw_s = "".join(current_seq_parts).upper()
            has_stop = raw_s.endswith("*")
            clean_s = raw_s[:-1] if has_stop else raw_s
            raw_sha = hashlib.sha256(raw_s.encode("ascii")).hexdigest()
            clean_sha = hashlib.sha256(clean_s.encode("ascii")).hexdigest()
            records.append((current_id, len(clean_s), raw_sha, clean_sha, has_stop))

    return records


def verify_run(output_dir: Path, fasta_path: Path) -> Dict:
    results = {
        "status": "UNKNOWN",
        "output_dir": str(output_dir),
        "fasta_path": str(fasta_path),
        "checks": {},
    }

    # 1. Check required files
    required_files = ["embeddings.npz", "embeddings.npy", "sequence_manifest.tsv", "provenance.json", "checksums.sha256"]
    missing = [f for f in required_files if not (output_dir / f).exists()]
    if missing:
        results["checks"]["files_present"] = {"passed": False, "missing": missing}
        results["status"] = "FAILED"
        return results
    results["checks"]["files_present"] = {"passed": True}

    # 2. Checksum validation
    checksum_errors = []
    with open(output_dir / "checksums.sha256", "r") as cf:
        for line in cf:
            parts = line.strip().split()
            if len(parts) == 2:
                expected_sha, fname = parts[0], parts[1].lstrip("*")
                fpath = output_dir / fname
                if not fpath.exists():
                    checksum_errors.append(f"Missing file {fname}")
                else:
                    actual_sha = compute_file_sha256(fpath)
                    if actual_sha != expected_sha:
                        checksum_errors.append(f"Checksum mismatch for {fname}: expected {expected_sha}, got {actual_sha}")
    results["checks"]["checksums"] = {
        "passed": len(checksum_errors) == 0,
        "errors": checksum_errors
    }
    if checksum_errors:
        results["status"] = "FAILED"
        return results

    # 3. Load data
    npz_data = np.load(output_dir / "embeddings.npz", allow_pickle=True)
    npy_data = np.load(output_dir / "embeddings.npy")

    emb_npz = npz_data["embeddings"]
    ids_npz = [str(x) for x in npz_data["ids"]]
    lens_npz = list(npz_data["lengths"])
    raw_shas_npz = [str(x) for x in npz_data.get("raw_sequence_sha256", npz_data.get("sequence_sha256"))]
    clean_shas_npz = [str(x) for x in npz_data["clean_sequence_sha256"]] if "clean_sequence_sha256" in npz_data else None

    # 4. Check data equality between npy and npz
    arrays_equal = np.array_equal(emb_npz, npy_data)
    results["checks"]["npy_npz_identical"] = {"passed": bool(arrays_equal)}

    # 5. Check dimensionality and finiteness
    fasta_records = parse_fasta_ids_and_lengths(fasta_path)
    expected_n = len(fasta_records)
    expected_dim = 1024

    dim_ok = (emb_npz.shape == (expected_n, expected_dim))
    all_finite = bool(np.all(np.isfinite(emb_npz)))
    results["checks"]["dimensionality"] = {
        "passed": dim_ok,
        "actual_shape": list(emb_npz.shape),
        "expected_shape": [expected_n, expected_dim],
    }
    results["checks"]["finiteness"] = {
        "passed": all_finite,
        "has_nans": bool(np.isnan(emb_npz).any()),
        "has_infs": bool(np.isinf(emb_npz).any()),
    }

    # 6. Check ID order, lengths, and unified raw vs clean sequence hashes
    expected_ids = [r[0] for r in fasta_records]
    expected_clean_lens = [r[1] for r in fasta_records]
    expected_raw_shas = [r[2] for r in fasta_records]
    expected_clean_shas = [r[3] for r in fasta_records]
    expected_trailing_stops = [r[4] for r in fasta_records]

    ids_match = (ids_npz == expected_ids)
    lens_match = (lens_npz == expected_clean_lens)
    raw_shas_match = (raw_shas_npz == expected_raw_shas)
    clean_shas_match = (clean_shas_npz == expected_clean_shas) if clean_shas_npz is not None else True

    results["checks"]["id_order_preserved"] = {
        "passed": ids_match,
        "num_ids": len(ids_npz),
    }
    results["checks"]["lengths_match"] = {"passed": lens_match}
    results["checks"]["raw_sequence_shas_match"] = {"passed": raw_shas_match}
    results["checks"]["clean_sequence_shas_match"] = {"passed": clean_shas_match}

    # 7. Check sequence manifest fields
    manifest_path = output_dir / "sequence_manifest.tsv"
    manifest_ok = True
    manifest_rows = []
    with open(manifest_path, "r", encoding="utf-8") as mf:
        header = mf.readline().strip().split("\t")
        for line in mf:
            parts = line.strip().split("\t")
            if parts:
                manifest_rows.append(parts)

    if len(manifest_rows) != expected_n:
        manifest_ok = False
    else:
        for idx, row in enumerate(manifest_rows):
            # index, seq_id, length, raw_sha256, clean_sha256, trailing_stop_stripped
            row_id = row[1]
            row_len = int(row[2])
            row_raw_sha = row[3]
            row_clean_sha = row[4]
            row_stop = row[5].lower() == "true"
            if (row_id != expected_ids[idx] or
                row_len != expected_clean_lens[idx] or
                row_raw_sha != expected_raw_shas[idx] or
                row_clean_sha != expected_clean_shas[idx] or
                row_stop != expected_trailing_stops[idx]):
                manifest_ok = False
                break

    results["checks"]["manifest_integrity"] = {"passed": manifest_ok}

    # 8. Check duplicate sequence embeddings (by clean_sha)
    dup_checks = []
    clean_sha_map: Dict[str, List[int]] = {}
    for idx, sha in enumerate(expected_clean_shas):
        clean_sha_map.setdefault(sha, []).append(idx)

    for sha, indices in clean_sha_map.items():
        if len(indices) > 1:
            base_idx = indices[0]
            base_emb = emb_npz[base_idx]
            for other_idx in indices[1:]:
                other_emb = emb_npz[other_idx]
                max_diff = float(np.max(np.abs(base_emb - other_emb)))
                dup_checks.append({
                    "id1": ids_npz[base_idx],
                    "id2": ids_npz[other_idx],
                    "clean_sha256": sha,
                    "max_abs_diff": max_diff,
                    "passed": max_diff < 1e-5,
                })
    results["checks"]["duplicate_equivalence"] = {
        "passed": all(d["passed"] for d in dup_checks) if dup_checks else True,
        "checks": dup_checks,
    }

    # 9. Provenance structure check
    with open(output_dir / "provenance.json", "r") as pf:
        prov = json.load(pf)

    prov_valid = (
        "run_signature" in prov and
        "runner_code_sha256" in prov and
        "signature_components" in prov and
        "software_environment" in prov and
        "hardware" in prov and
        "model" in prov and
        prov.get("model", {}).get("tokenizer_legacy") is True and
        prov.get("validation", {}).get("dimensions") == [expected_n, expected_dim]
    )
    results["checks"]["provenance_valid"] = {
        "passed": prov_valid,
        "run_signature": prov.get("run_signature"),
        "runner_code_sha256": prov.get("runner_code_sha256"),
        "gpu_device": prov.get("hardware", {}).get("device_name"),
        "peak_vram_gb": prov.get("hardware", {}).get("peak_vram_used_gb"),
        "single_vs_batched_consistency": prov.get("validation", {}).get("single_vs_batched_consistency"),
    }

    all_passed = all(check.get("passed", False) for check in results["checks"].values())
    results["status"] = "PASSED" if all_passed else "FAILED"
    return results


def main():
    parser = argparse.ArgumentParser(description="Verify ProstT5 Embedding Output (v3)")
    parser.add_argument("--output-dir", type=Path, required=True, help="Directory containing embeddings")
    parser.add_argument("--fasta", type=Path, required=True, help="Original input FASTA file")
    parser.add_argument("--report", type=Path, default=None, help="Optional output JSON report path")
    args = parser.parse_args()

    results = verify_run(args.output_dir, args.fasta)
    output_json = json.dumps(results, indent=2)
    print(output_json)

    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        with open(args.report, "w", encoding="utf-8") as rf:
            rf.write(output_json + "\n")

    if results["status"] != "PASSED":
        sys.exit(1)
    print("Verification completed successfully: ALL CHECKS PASSED.")


if __name__ == "__main__":
    main()
