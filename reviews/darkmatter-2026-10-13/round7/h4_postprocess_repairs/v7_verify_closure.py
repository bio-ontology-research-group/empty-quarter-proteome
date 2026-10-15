#!/usr/bin/env python3
"""v7_verify_closure.py — Independent code-data closure and audit verification
for Hypothesis 4 V7 full-cohort outputs.

Root-Audited Contract:
- Standalone execution defaults to strict verification mode (mode="strict"): fails closed
  immediately if v7_results.sha256.tsv is missing or if any required artifact is missing/unmatched.
- In strict mode, verifies all listed files including v7_closure_verification.json and
  v7_raw_search_hashes.sha256.tsv without rewriting the receipt.
- Provisional mode (mode="provisional") verifies data tables against provisional manifest,
  runs scientific closure assertions, and writes the verification receipt.
- Verifies exact finite-population manifest (5,851 models, 7,773 regions, 161 0-region models).
- Verifies exact sub-stratum population totals across all 6 strata.
- Verifies category sum closure across all strata and Tiers 0-4.
- Verifies Tier 0-4 monotonicity (T0 >= T1 >= T2 >= T3 >= T4) for remaining hits and qTM scores.
- Verifies whole-chain / crop exclusion symmetry (Tier 0 == Tier 1 for dark candidates).
- Verifies that undefined domains return null/None (no [0.0, 0.0] intervals).
- Audits CIGAR M operations against (qspan + tspan - alnlen) for all audited example hits.
"""

import os
import re
import csv
import json
import hashlib
import argparse
from datetime import datetime, timezone
from collections import Counter, defaultdict

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES_DIR = os.path.join(BASE_DIR, "results")
CIGAR_RE = re.compile(r"(\d+)([MID])")

REQUIRED_DATA_ARTIFACTS = [
    "v7_full5851_cropped_regions_manifest.tsv",
    "v7_region_evaluations.tsv",
    "v7_model_summaries.tsv",
    "v7_stratum_census.tsv",
    "v7_figure_data.tsv",
    "v7_recovery_summary.json",
    "v7_audited_example_hits.tsv",
    "v7_raw_search_hashes.sha256.tsv",
]


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1048576), b""):
            h.update(chunk)
    return h.hexdigest()


def verify_closure(res_dir=RES_DIR, mode="strict"):
    print(f"=== Running Independent V7 Code-Data Closure Verification (mode={mode}) in {res_dir} ===")

    # 1. Checksum manifest verification (Strictly required in both strict and provisional modes)
    sha_path = os.path.join(res_dir, "v7_results.sha256.tsv")
    assert os.path.exists(sha_path), f"FAIL-CLOSED: Missing required checksum manifest {sha_path}"

    # 2. Fail-closed check: all required data artifacts must exist on disk
    for fn in REQUIRED_DATA_ARTIFACTS:
        fpath = os.path.join(res_dir, fn)
        assert os.path.exists(fpath), f"FAIL-CLOSED: Required artifact {fn} missing from {res_dir}"
        assert os.path.getsize(fpath) > 0, f"FAIL-CLOSED: Required artifact {fn} is empty in {res_dir}"

    if mode == "strict":
        receipt_path = os.path.join(res_dir, "v7_closure_verification.json")
        assert os.path.exists(receipt_path), f"FAIL-CLOSED: Missing receipt {receipt_path} in strict mode"
        assert os.path.getsize(receipt_path) > 0, f"FAIL-CLOSED: Empty receipt {receipt_path}"

    expected_hashes = {}
    with open(sha_path) as fh:
        for line in fh:
            parts = line.strip().split()
            if len(parts) == 2:
                expected_hashes[parts[1]] = parts[0]

    assert "v7_results.sha256.tsv" not in expected_hashes, "FAIL: v7_results.sha256.tsv contains itself!"

    target_artifacts = REQUIRED_DATA_ARTIFACTS + (["v7_closure_verification.json"] if mode == "strict" else [])
    for fn in target_artifacts:
        assert fn in expected_hashes, f"FAIL-CLOSED: Required artifact {fn} not listed in {sha_path}"
        actual = sha256_file(os.path.join(res_dir, fn))
        assert actual == expected_hashes[fn], f"FAIL: SHA-256 mismatch for {fn}: {actual} != {expected_hashes[fn]}"
        print(f"  [PASS] SHA-256 verified: {fn} ({actual[:16]}...)")

    # In strict mode, verify every single file listed in the manifest matches disk
    if mode == "strict":
        for fn, expected_digest in expected_hashes.items():
            fpath = os.path.join(res_dir, fn)
            assert os.path.exists(fpath), f"FAIL-CLOSED: Manifest lists file {fn} which does not exist in {res_dir}"
            actual = sha256_file(fpath)
            assert actual == expected_digest, f"FAIL: SHA-256 mismatch for {fn}: {actual} != {expected_digest}"
        print(f"  [PASS] All {len(expected_hashes)} files in final manifest strictly verified against disk.")

    # 3. Load tables
    with open(os.path.join(res_dir, "v7_full5851_cropped_regions_manifest.tsv")) as fh:
        full_manifest = list(csv.DictReader(fh, delimiter="\t"))
    with open(os.path.join(res_dir, "v7_region_evaluations.tsv")) as fh:
        regions = list(csv.DictReader(fh, delimiter="\t"))
    with open(os.path.join(res_dir, "v7_model_summaries.tsv")) as fh:
        models = list(csv.DictReader(fh, delimiter="\t"))
    with open(os.path.join(res_dir, "v7_recovery_summary.json")) as fh:
        summary = json.load(fh)
    with open(os.path.join(res_dir, "v7_audited_example_hits.tsv")) as fh:
        examples = list(csv.DictReader(fh, delimiter="\t"))
    with open(os.path.join(res_dir, "v7_stratum_census.tsv")) as fh:
        census_rows = list(csv.DictReader(fh, delimiter="\t"))
    with open(os.path.join(res_dir, "v7_figure_data.tsv")) as fh:
        figure_rows = list(csv.DictReader(fh, delimiter="\t"))

    # 4. Verify full N=5,851 manifest closure
    full_by_model = defaultdict(list)
    for r in full_manifest:
        full_by_model[r["model_id"]].append(r)
        assert len(r["source_pdb_sha256"]) == 64, f"Invalid source_pdb_sha256 in manifest: {r}"
        if r["cropped_pdb"] != "NONE":
            assert len(r["crop_pdb_sha256"]) == 64, f"Invalid crop_pdb_sha256 in manifest: {r}"

    assert len(full_by_model) == 5851, f"Expected 5851 unique models in full manifest, got {len(full_by_model)}"
    full_sub_counts = Counter(rows[0]["sub_stratum"] for rows in full_by_model.values())
    expected_sub_counts = {
        "K_local_pdb": 849,
        "K_local_afdb_only": 87,
        "K_dual_none": 95,
        "Q_not_K_pdb": 1107,
        "Q_not_K_afdb_only": 700,
        "Controls": 3013,
    }
    assert full_sub_counts == expected_sub_counts, f"Full manifest sub-stratum count mismatch: {full_sub_counts}"
    full_no_reg = sum(1 for rows in full_by_model.values() if int(rows[0]["n_confident_regions"]) == 0)
    full_crops = sum(1 for r in full_manifest if r["cropped_pdb"] != "NONE")
    assert full_no_reg == 161 and full_crops == 7773 and len(full_manifest) == 7934, (
        f"Manifest counts mismatch: 0-reg={full_no_reg}, crops={full_crops}, total_rows={len(full_manifest)}"
    )
    print("  [PASS] Full finite-population (N=5,851) region-eligibility manifest closure (7,773 regions + 161 0-region models).")

    # 5. Verify evaluated cohort size and strata
    assert len(models) == 5851, f"Expected 5851 models, got {len(models)}"
    assert len(regions) == 7773, f"Expected 7773 regions, got {len(regions)}"
    sub_counts = Counter(m["sub_stratum"] for m in models)
    assert sub_counts == expected_sub_counts, f"Model summaries sub-stratum mismatch: {sub_counts}"
    print("  [PASS] Evaluated full census cohort closure (N=5,851 models, R=7,773 regions across 6 sub-strata).")

    # 6. Verify no-region exclusion, tier monotonicity, and whole/crop symmetry
    for m in models:
        n_reg = int(m["n_confident_regions"])
        is_ctrl = int(m["is_control"]) == 1
        for k in range(5):
            st = m[f"crop_t{k}_status"]
            tr = m[f"trans_t{k}"]
            if n_reg == 0:
                assert st == "NO_CONFIDENT_REGION", f"0-reg model {m['model_id']} has status {st}"
                assert tr in ("WC_UNMATCHED_NO_CONFIDENT_REGION", "WC_MATCHED_LOST_NO_CONFIDENT_REGION")
            else:
                assert st != "NO_CONFIDENT_REGION", f"Model {m['model_id']} with {n_reg} regions has {st}"

        for k in range(4):
            assert int(m[f"wc_t{k}_remaining_hits"]) >= int(m[f"wc_t{k+1}_remaining_hits"])
            assert float(m[f"wc_t{k}_qtm"]) + 1e-9 >= float(m[f"wc_t{k+1}_qtm"])
            assert float(m[f"crop_t{k}_max_qtm"]) + 1e-9 >= float(m[f"crop_t{k+1}_max_qtm"])

        if not is_ctrl:
            assert int(m["wc_t0_remaining_hits"]) == int(m["wc_t1_remaining_hits"])
            assert float(m["wc_t0_qtm"]) == float(m["wc_t1_qtm"])
            assert float(m["crop_t0_max_qtm"]) == float(m["crop_t1_max_qtm"])
            assert m["trans_t0"] == m["trans_t1"]

    reg_by_crop_id = {r["crop_id"]: r for r in regions}
    for r in regions:
        for k in range(4):
            assert int(r[f"t{k}_remaining_hits"]) >= int(r[f"t{k+1}_remaining_hits"])
            assert float(r[f"t{k}_qtm"]) + 1e-9 >= float(r[f"t{k+1}_qtm"])
        for k in (0, 1):
            assert int(r[f"t{k}_aligned_pairs"]) <= int(r[f"t{k}_qspan"]) <= int(r["reg_len"])
            assert float(r[f"t{k}_aligned_pair_qcov"]) <= float(r[f"t{k}_qspan_cov"]) + 1e-6

    print("  [PASS] No-region exclusion, Tier 0-4 monotonicity, and whole/crop exclusion symmetry verified.")

    # 7. Verify stratum summary totals in JSON
    for st_name, st_data in summary["stratum_summaries"].items():
        n_m = st_data["model_count"]
        n_r = st_data["total_confident_regions"]
        for k in range(5):
            tk = st_data["tiers"][f"tier_{k}"]
            assert sum(tk["wc_status_counts"].values()) == n_m
            assert sum(tk["crop_model_status_counts"].values()) == n_m
            assert sum(tk["transition_counts"].values()) == n_m
            assert sum(tk["region_status_counts"].values()) == n_r

            sp = tk.get("model_based_superpopulation_inferential_statistics", {})
            for ci_key in [
                "unconditional_crop_match_wilson_ci95",
                "region_eligible_crop_match_wilson_ci95",
                "wc_unmatched_recovery_unconditional_wilson_ci95",
                "wc_unmatched_recovery_eligible_wilson_ci95",
                "wc_matched_retained_wilson_ci95",
            ]:
                val = sp.get(ci_key)
                if val is not None:
                    assert val != [0.0, 0.0], f"Found misleading [0.0, 0.0] in {st_name} {k} {ci_key}"
                    assert val != [100.0, 100.0], f"Found misleading [100.0, 100.0] in {st_name} {k} {ci_key}"

    print("  [PASS] Stratum summary category totals and absence of misleading CIs verified across all tiers.")

    # 8. Verify primary census table and figure data contain NO sampling CIs
    expected_strata = set(expected_sub_counts) | {"K_local_all", "Q_not_K_all"}
    assert len(census_rows) == len(expected_strata) * 5
    assert {(r["stratum"], r["tier"]) for r in census_rows} == {(s, f"tier_{k}") for s in expected_strata for k in range(5)}
    census_fieldnames = list(census_rows[0].keys())
    for col in census_fieldnames:
        assert "ci95" not in col.lower(), f"Primary stratum census table must NOT contain CI columns: {col}"

    figure_fieldnames = list(figure_rows[0].keys())
    for col in figure_fieldnames:
        assert "ci95" not in col.lower(), f"Primary figure data table must NOT contain CI columns: {col}"

    print("  [PASS] Primary stratum census and figure data verified free of sampling-error CIs.")

    # 9. Independently verify representative hit rows in v7_audited_example_hits.tsv
    mod_by_id = {m["model_id"]: m for m in models}
    matched_top_reg = 0
    matched_top_wc = 0
    for ex in examples:
        q = ex["query"]
        qlen = int(ex["qlen"])
        alnlen = int(ex["alnlen"])
        qstart, qend = int(ex["qstart"]), int(ex["qend"])
        tstart, tend = int(ex["tstart"]), int(ex["tend"])
        qspan = qend - qstart + 1
        tspan = tend - tstart + 1
        cigar = ex["cigar"]
        cigar_m = sum(int(n) for n, op in CIGAR_RE.findall(cigar) if op == "M")
        cigar_tot = sum(int(n) for n, op in CIGAR_RE.findall(cigar))
        assert cigar_tot == alnlen, f"CIGAR length mismatch for {q}"
        assert cigar_m == qspan + tspan - alnlen == int(ex["aligned_pairs"])
        assert abs(round(cigar_m / qlen, 4) - float(ex["aligned_pair_qcov"])) < 1e-6
        assert abs(round(qspan / qlen, 4) - float(ex["qspan_cov"])) < 1e-6

        if q in reg_by_crop_id:
            rr = reg_by_crop_id[q]
            for k in range(5):
                if rr[f"t{k}_top_target"] == ex["target"] and abs(float(rr[f"t{k}_qtm"]) - round(float(ex["qtm"]), 4)) < 1e-4:
                    assert abs(float(rr[f"t{k}_aligned_pair_qcov"]) - float(ex["aligned_pair_qcov"])) < 1e-6
                    assert abs(float(rr[f"t{k}_qspan_cov"]) - float(ex["qspan_cov"])) < 1e-6
                    if k in (0, 1):
                        assert int(rr[f"t{k}_aligned_pairs"]) == cigar_m
                        assert int(rr[f"t{k}_qspan"]) == qspan
                    matched_top_reg += 1
                    break
        elif q in mod_by_id:
            mm = mod_by_id[q]
            for k in range(5):
                if abs(float(mm[f"wc_t{k}_qtm"]) - round(float(ex["qtm"]), 4)) < 1e-4 and abs(float(mm[f"wc_t{k}_fident"]) - round(float(ex["fident"]), 4)) < 1e-4:
                    assert abs(float(mm[f"wc_t{k}_aligned_pair_qcov"]) - float(ex["aligned_pair_qcov"])) < 1e-6
                    assert abs(float(mm[f"wc_t{k}_qspan_cov"]) - float(ex["qspan_cov"])) < 1e-6
                    matched_top_wc += 1
                    break

    assert matched_top_reg >= 20 and matched_top_wc >= 15, (
        f"Expected >=20 top reg and >=15 top wc rows, got {matched_top_reg}, {matched_top_wc}"
    )
    print(f"  [PASS] Audited {len(examples)} representative raw hit rows (CIGAR aligned pairs vs span; {matched_top_reg} top region, {matched_top_wc} top wc).")

    # 10. Write verification receipt only in provisional mode
    receipt_path = os.path.join(res_dir, "v7_closure_verification.json")
    if mode == "provisional":
        receipt = {
            "status": "PASS",
            "timestamp_utc": datetime.now(timezone.utc).isoformat(),
            "verified_models": len(models),
            "verified_cropped_regions": len(regions),
            "models_no_confident_region": full_no_reg,
            "sub_strata_verified": expected_sub_counts,
            "manifest_checksums_verified": True,
            "all_closure_assertions_passed": True,
            "verified_artifacts": REQUIRED_DATA_ARTIFACTS,
        }
        with open(receipt_path, "w") as fh:
            json.dump(receipt, fh, indent=2)
        print(f"  [PASS] Verification receipt written to {receipt_path}")
    else:
        with open(receipt_path) as fh:
            receipt = json.load(fh)
        assert receipt.get("status") == "PASS", f"Receipt status is not PASS in {receipt_path}"
        print(f"  [PASS] Verified existing closure verification receipt ({receipt_path}).")

    print(f"=== ALL INDEPENDENT V7 CODE-DATA CLOSURE CHECKS PASSED SUCCESSFULLY (mode={mode}) ===")
    return receipt


def main():
    parser = argparse.ArgumentParser(description="H4 V7 Independent Closure Verification")
    parser.add_argument("--mode", choices=["strict", "provisional"], default="strict", help="Verification mode: strict (default) or provisional")
    parser.add_argument("--dir", default=RES_DIR, help="Path to results directory")
    args = parser.parse_args()

    verify_closure(res_dir=args.dir, mode=args.mode)


if __name__ == "__main__":
    main()
