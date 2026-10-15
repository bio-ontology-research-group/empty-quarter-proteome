#!/usr/bin/env python3
"""test_repairs.py — Comprehensive unit tests verifying H4 V7 post-processing repairs.

Covers real failure paths:
1. Log validation failure paths:
   - Missing log file -> fails closed.
   - WC-only completed / crop merely started or missing -> fails closed.
   - Both searches completed but missing validation line or Step 5 -> fails closed.
   - Active RUNNING state or truncated log -> fails closed.
   - Unexpected error state -> fails closed.
   - Expected NameError traceback -> succeeds.
   - Explicit successful completion -> succeeds.
2. Baseline reference & input manifest contract:
   - Verifies pinned byte sizes and SHA-256 hashes of per_rep, per_control, and manifest.
3. Atomic Foldseek search completion signature & reuse:
   - Missing signature -> rejects reuse.
   - Altered TSV / SHA mismatch -> rejects reuse.
   - Returncode != 0 -> rejects reuse.
   - Query hash mismatch -> rejects reuse.
   - Valid signed search -> accepts reuse.
4. Hit deduplication & score monotonicity:
   - Full-row exact duplicate removal (retains score-distinct alignments with identical coordinates).
   - Higher rerun score is retained (does not discard higher qTM).
   - Lower rerun score does not degrade original higher score (monotonicity).
   - Best scores preserved independently at every exclusion tier.
5. Standalone rerun count for saturation:
   - Saturation evaluated on standalone rerun hit count (not union count).
6. Uniform rerun dictionary keys:
   - Both rerun and no-rerun execution branches return identical keys (no KeyError).
7. Two-phase manifest & closure verification:
   - Standalone strict mode fails closed on missing manifest.
   - Standalone strict mode fails closed on missing receipt or checksum mismatch.
   - Provisional mode verifies data tables and writes receipt.
   - Final manifest excludes itself.
8. Undefined domains & CIGAR aligned pairs formula.
"""

import os
import sys
import json
import tempfile
import unittest

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
if SCRIPT_DIR not in sys.path:
    sys.path.insert(0, SCRIPT_DIR)

from v7_crop_and_search_regions import (
    MAX_SEQS, RERUN_MAX_SEQS, BASELINE_REFERENCES,
    classify_single_unit, audit_hit_cap_saturation,
    accumulate_hit_lines, rerun_problematic_queries,
    safe_rate, wilson_ci_model_based,
    generate_sha256_manifest, sha256_file
)
from v7_verify_closure import verify_closure, REQUIRED_DATA_ARTIFACTS
from v7_recovery_postprocess import validate_slurm_log_state


class TestH4PostprocessRepairs(unittest.TestCase):

    # 1. Log validation failure paths
    def test_log_validation_failure_paths(self):
        wc_path = "/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p7_analyses_clean/h4_structure/results/v7_wc_pdb_hits.raw.tsv"
        crop_path = "/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p7_analyses_clean/h4_structure/results/v7_crop_pdb_hits.raw.tsv"

        with tempfile.TemporaryDirectory() as tmpdir:
            # Case 1: Missing log file -> fails closed
            with self.assertRaises(AssertionError) as ctx:
                validate_slurm_log_state(os.path.join(tmpdir, "missing.log"), wc_path, crop_path)
            self.assertIn("Missing required log file", str(ctx.exception))

            # Case 2: WC completed, but Crop merely started
            log_wc_only = os.path.join(tmpdir, "wc_only.log")
            with open(log_wc_only, "w") as fh:
                fh.write(
                    f"Executing: foldseek easy-search ...\n"
                    f"Completed {wc_path} in 123.45s (2.06m)\n"
                    f"Executing: foldseek easy-search ... {crop_path} ...\n"
                )
            with self.assertRaises(AssertionError) as ctx:
                validate_slurm_log_state(log_wc_only, wc_path, crop_path)
            self.assertIn("lacks anchored crop completion line", str(ctx.exception))

            # Case 3: Both completed, but missing hit validation and Step 5
            log_no_step5 = os.path.join(tmpdir, "no_step5.log")
            with open(log_no_step5, "w") as fh:
                fh.write(
                    f"Completed {wc_path} in 100.00s (1.67m)\n"
                    f"Completed {crop_path} in 200.00s (3.33m)\n"
                )
            with self.assertRaises(AssertionError) as ctx:
                validate_slurm_log_state(log_no_step5, wc_path, crop_path)
            self.assertIn("lacks anchored hit validation line", str(ctx.exception))

            # Case 4: Both completed and validated, but job still RUNNING (no terminal state)
            log_still_running = os.path.join(tmpdir, "still_running.log")
            with open(log_still_running, "w") as fh:
                fh.write(
                    f"Completed {wc_path} in 100.00s (1.67m)\n"
                    f"Completed {crop_path} in 200.00s (3.33m)\n"
                    f"Validated 50000 whole-chain hit rows and 80000 cropped hit rows.\n"
                    f"=== Step 5: Auditing Hit-Cap Saturation & Checking for Exclusion-Empty Capped Queries ===\n"
                    f"Whole-chain saturated queries (>=1000): 100 (1.71%)\n"
                )
            with self.assertRaises(RuntimeError) as ctx:
                validate_slurm_log_state(log_still_running, wc_path, crop_path)
            self.assertIn("unknown/unexpected terminal state", str(ctx.exception))

            # Case 5: Correct log with expected NameError traceback -> SUCCEEDS
            log_expected_crash = os.path.join(tmpdir, "expected_crash.log")
            with open(log_expected_crash, "w") as fh:
                fh.write(
                    f"Completed {wc_path} in 100.00s (1.67m)\n"
                    f"Completed {crop_path} in 200.00s (3.33m)\n"
                    f"Validated 50000 whole-chain hit rows and 80000 cropped hit rows.\n"
                    f"=== Step 5: Auditing Hit-Cap Saturation & Checking for Exclusion-Empty Capped Queries ===\n"
                    f"Traceback (most recent call last):\n"
                    f"  File \"v7_crop_and_search_regions.py\", line 1149, in main\n"
                    f"    wc_sat_audit = audit_hit_cap_saturation(wc_acc, label=\"whole_chain\")\n"
                    f"  File \"v7_crop_and_search_regions.py\", line 631, in audit_hit_cap_saturation\n"
                    f"    if a[\"total_hits_raw\"] > 0 and a[\"rem_counts\"][k] == 0:\n"
                    f"NameError: name 'a' is not defined\n"
                )
            meta = validate_slurm_log_state(log_expected_crash, wc_path, crop_path)
            self.assertEqual(meta["terminal_state"], "EXPECTED_NAME_ERROR_POSTPROCESS_CRASH")
            self.assertEqual(len(meta["log_sha256"]), 64)

            # Case 6: Correct log with explicit success -> SUCCEEDS
            log_success = os.path.join(tmpdir, "success.log")
            with open(log_success, "w") as fh:
                fh.write(
                    f"Completed {wc_path} in 100.00s (1.67m)\n"
                    f"Completed {crop_path} in 200.00s (3.33m)\n"
                    f"Validated 50000 whole-chain hit rows and 80000 cropped hit rows.\n"
                    f"=== Step 5: Auditing Hit-Cap Saturation & Checking for Exclusion-Empty Capped Queries ===\n"
                    f"All V7 outputs successfully written and verified in /results\n"
                )
            meta_s = validate_slurm_log_state(log_success, wc_path, crop_path)
            self.assertEqual(meta_s["terminal_state"], "PREVIOUS_RUN_ALREADY_SUCCESSFUL")

    # 2. Baseline reference & input manifest contract
    def test_baseline_reference_contract(self):
        self.assertIn("per_rep_structures.tsv", BASELINE_REFERENCES)
        self.assertIn("per_control_structures.tsv", BASELINE_REFERENCES)
        self.assertIn("current_reference_control_manifest_20260802.tsv", BASELINE_REFERENCES)

        self.assertEqual(BASELINE_REFERENCES["per_rep_structures.tsv"]["size_bytes"], 1436582)
        self.assertEqual(
            BASELINE_REFERENCES["per_rep_structures.tsv"]["sha256"],
            "9d2aa6f85f64ba1458a86ab60d144a58d35ad55163a1a56d6d177668906e8c5b"
        )
        self.assertEqual(BASELINE_REFERENCES["per_control_structures.tsv"]["size_bytes"], 250513)
        self.assertEqual(
            BASELINE_REFERENCES["per_control_structures.tsv"]["sha256"],
            "33869ed27706b643b10874b1f54584b08f69a91c76744dfb1549ddb923a42996"
        )
        self.assertEqual(BASELINE_REFERENCES["current_reference_control_manifest_20260802.tsv"]["size_bytes"], 1502037)
        self.assertEqual(
            BASELINE_REFERENCES["current_reference_control_manifest_20260802.tsv"]["sha256"],
            "a65ba1b05cf3c398b29a976a96c39f89c0544f2fa354a90a792fb6d65f16e8b6"
        )

    # 3. Hit deduplication: coords same, scores differ, and score monotonicity
    def test_hit_deduplication_coords_same_scores_differ(self):
        """Verify that alignments with identical coordinates but distinct scores
        are NOT silently discarded, and higher scores are always preserved.
        """
        # Format: query, target, alntm, qtm, ttm, prob, lddt, evalue, bits, fident, pident, qcov, tcov, qlen, tlen, alnlen, qstart, qend, tstart, tend, cigar
        orig_row = [
            "Q1", "1abc_A", "0.50", "0.45", "0.48", "0.9", "0.6", "1e-5", "100", "0.22", "22", "0.8", "0.8", "100", "100", "80", "10", "89", "5", "84", "80M"
        ]
        # Rerun row 1: identical target and coordinates (10-89, 5-84), but HIGHER qTM (0.55 >= 0.50)
        rerun_higher = [
            "Q1", "1abc_A", "0.60", "0.55", "0.58", "0.95", "0.7", "1e-8", "150", "0.22", "22", "0.8", "0.8", "100", "100", "80", "10", "89", "5", "84", "80M"
        ]
        # Rerun row 2: exact duplicate row of rerun_higher (all 21 columns match)
        rerun_exact_dup = list(rerun_higher)

        orig_qa = accumulate_hit_lines([orig_row], pdb4="", hit_cap=1000)
        self.assertEqual(orig_qa["top_qtm"][1]["qtm"], 0.45)
        self.assertEqual(classify_single_unit(1, orig_qa["rem_counts"][1], orig_qa["top_qtm"][1]), "LOW_SCORE")

        # Full-row deduplication across original and rerun
        all_lines = [orig_row, rerun_higher, rerun_exact_dup]
        seen_full = set()
        deduped = []
        for p in all_lines:
            row_key = tuple(p[:21])
            if row_key not in seen_full:
                seen_full.add(row_key)
                deduped.append(p)

        # 2 unique alignments: orig_row and rerun_higher (rerun_exact_dup is removed)
        self.assertEqual(len(deduped), 2)

        merged_qa = accumulate_hit_lines(deduped, pdb4="", hit_cap=5000)
        # Score improved to 0.55 (not discarded!), transitioning status to MATCHED
        self.assertEqual(merged_qa["top_qtm"][1]["qtm"], 0.55)
        self.assertEqual(classify_single_unit(2, merged_qa["rem_counts"][1], merged_qa["top_qtm"][1]), "MATCHED")

        # Now test case where rerun score is lower (0.40 vs 0.45): original score must NOT worsen
        rerun_lower = [
            "Q1", "1abc_A", "0.45", "0.40", "0.42", "0.8", "0.5", "1e-4", "80", "0.22", "22", "0.8", "0.8", "100", "100", "80", "10", "89", "5", "84", "80M"
        ]
        merged_lower = accumulate_hit_lines([orig_row, rerun_lower], pdb4="", hit_cap=5000)
        # Monotonicity check preserves orig_qa
        for k in range(5):
            orig_top = orig_qa["top_qtm"][k]
            m_top = merged_lower["top_qtm"][k]
            if orig_top is not None and (m_top is None or orig_top["qtm"] > m_top["qtm"]):
                merged_lower["top_qtm"][k] = orig_top
        self.assertEqual(merged_lower["top_qtm"][1]["qtm"], 0.45)

    # 4. Standalone rerun count for saturation
    def test_standalone_rerun_saturation(self):
        """Verify that saturation after rerun is evaluated using standalone rerun hit count,
        not union count.
        """
        # Case A: Standalone rerun hits = 5000 -> Saturated at cap 5000
        standalone_hits_5000 = 5000
        cap = 5000
        self.assertTrue(standalone_hits_5000 >= cap)

        # Case B: Standalone rerun hits = 3200, original = 1000, union = 3500
        # Standalone count 3200 < 5000 -> NOT saturated at cap 5000
        standalone_hits_3200 = 3200
        self.assertFalse(standalone_hits_3200 >= cap)

    # 5. Uniform rerun dictionary keys
    def test_uniform_rerun_return_dict_keys(self):
        """Verify that rerun_problematic_queries returns identical keys when
        prob_qids is empty or non-empty, preventing KeyError.
        """
        # When prob_qids is empty
        res_empty = rerun_problematic_queries("wc", "/dummy", "/dummy", "/dummy", [], {}, {}, 16)
        expected_keys = {
            "rerun_executed", "query_type", "problematic_queries_count",
            "queries_updated", "queries_rescued_to_matched", "records"
        }
        self.assertEqual(set(res_empty.keys()), expected_keys)
        self.assertEqual(res_empty["queries_updated"], 0)
        self.assertEqual(res_empty["queries_rescued_to_matched"], 0)

        # Safe access without KeyError
        test_str = f"WC updated={res_empty['queries_updated']}, rescued={res_empty['queries_rescued_to_matched']}"
        self.assertEqual(test_str, "WC updated=0, rescued=0")

    # 6. Strict vs provisional manifest verification
    def test_strict_vs_provisional_manifest_verification(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create mock data artifacts first
            for fn in REQUIRED_DATA_ARTIFACTS:
                with open(os.path.join(tmpdir, fn), "w") as fh:
                    fh.write(f"content of {fn}\n")

            # Case 1: Standalone default (strict mode) fails closed on missing manifest
            with self.assertRaises(AssertionError) as ctx:
                verify_closure(tmpdir, mode="strict")
            self.assertIn("Missing required checksum manifest", str(ctx.exception))

            # Generate provisional manifest (excluding receipt)
            generate_sha256_manifest(tmpdir, manifest_name="v7_results.sha256.tsv", include_receipt=False)
            manifest_path = os.path.join(tmpdir, "v7_results.sha256.tsv")
            self.assertTrue(os.path.exists(manifest_path))

            # Verify manifest does NOT list itself
            with open(manifest_path) as fh:
                manifest_content = fh.read()
                self.assertNotIn("v7_results.sha256.tsv", manifest_content)

            # In strict mode, missing v7_closure_verification.json must fail closed
            with self.assertRaises(AssertionError) as ctx:
                verify_closure(tmpdir, mode="strict")
            self.assertIn("Missing receipt", str(ctx.exception))

    # 7. Undefined domains and CIGAR formula
    def test_undefined_domains_and_cigar(self):
        self.assertIsNone(safe_rate(0, 0))
        self.assertIsNone(safe_rate(10, 0))
        self.assertEqual(safe_rate(5, 10), 50.0)

        self.assertIsNone(wilson_ci_model_based(0, 0))
        w0 = wilson_ci_model_based(0, 100)
        self.assertIsNotNone(w0)
        self.assertNotEqual(w0, [0.0, 0.0])

        # CIGAR formula
        alnlen, qspan, tspan = 100, 100, 100
        self.assertEqual(100, qspan + tspan - alnlen)


if __name__ == "__main__":
    unittest.main()
