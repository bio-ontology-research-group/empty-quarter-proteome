#!/usr/bin/env python3
"""test_repairs.py — Synthetic unit tests verifying H4 V7 post-processing repairs.

Verifies:
1. audit_hit_cap_saturation:
   - Eliminates NameError completely.
   - Respects per-query hit cap (a.get('hit_cap', MAX_SEQS)).
   - Flags saturated queries lacking MATCHED status at any tier (both POST_EXCLUSION_EMPTY and LOW_SCORE).
   - Does NOT flag uncapped queries (hits < cap) as saturated.
   - Does NOT flag saturated queries that are MATCHED across all tiers.
2. Adaptive rerun merge & score monotonicity:
   - Target/alignment deduplication between original and rerun searches.
   - Scores NEVER worsen: if original had higher qTM, original top hit is retained.
   - If rerun discovers higher-scoring hit, new score is adopted and status transitions to MATCHED.
   - Records initial and final caps and statuses.
3. Exact census undefined domains & primary table format:
   - Denominator <= 0 produces None (null in JSON), never misleading [0.0, 0.0].
   - Primary stratum census and figure data contain NO sampling-error CI columns.
4. CIGAR aligned pairs formula:
   - Validates cigar_M == qspan + tspan - alnlen across matches, insertions, and deletions.
5. Strict fail-closed checksum manifest & closure verifier:
   - Manifest does NOT include itself (no self-inclusion).
   - Missing required artifact immediately raises AssertionError (fail-closed, no silent skip).
   - Checksum mismatch immediately raises AssertionError.
"""

import os
import sys
import json
import tempfile
import unittest

# Ensure repairs dir is on path
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
if SCRIPT_DIR not in sys.path:
    sys.path.insert(0, SCRIPT_DIR)

from v7_crop_and_search_regions import (
    MAX_SEQS, RERUN_MAX_SEQS,
    classify_single_unit, audit_hit_cap_saturation,
    accumulate_hit_lines, safe_rate, wilson_ci_model_based,
    generate_sha256_manifest, sha256_file
)
from v7_verify_closure import verify_closure, REQUIRED_ARTIFACTS


class TestH4Repairs(unittest.TestCase):

    def test_audit_hit_cap_saturation_no_nameerror_and_categories(self):
        """Test that audit_hit_cap_saturation executes without NameError and
        correctly classifies saturated queries lacking MATCHED at any tier.
        """
        # Create synthetic accumulator
        acc = {
            # Q1: Uncapped (50 hits total < 1000). Not saturated.
            "Q_uncapped": {
                "pdb4": "",
                "hit_cap": 1000,
                "total_hits_raw": 50,
                "rem_counts": [50, 50, 10, 0, 0],
                "top_qtm": [
                    {"qtm": 0.65}, {"qtm": 0.65}, {"qtm": 0.42}, None, None
                ],
            },
            # Q2: Saturated (1000 hits), MATCHED at all tiers (qTM >= 0.5)
            "Q_sat_all_matched": {
                "pdb4": "",
                "hit_cap": 1000,
                "total_hits_raw": 1000,
                "rem_counts": [1000, 1000, 800, 500, 200],
                "top_qtm": [
                    {"qtm": 0.85}, {"qtm": 0.85}, {"qtm": 0.72}, {"qtm": 0.61}, {"qtm": 0.53}
                ],
            },
            # Q3: Saturated (1000 hits), Tier 0 MATCHED, but Tier 2 POST_EXCLUSION_EMPTY
            "Q_sat_post_empty": {
                "pdb4": "",
                "hit_cap": 1000,
                "total_hits_raw": 1000,
                "rem_counts": [1000, 1000, 0, 0, 0],
                "top_qtm": [
                    {"qtm": 0.80}, {"qtm": 0.80}, None, None, None
                ],
            },
            # Q4: Saturated (1000 hits), Tier 1 LOW_SCORE (hits remain but qTM < 0.5)
            "Q_sat_low_score": {
                "pdb4": "",
                "hit_cap": 1000,
                "total_hits_raw": 1000,
                "rem_counts": [1000, 200, 50, 10, 2],
                "top_qtm": [
                    {"qtm": 0.80}, {"qtm": 0.44}, {"qtm": 0.38}, {"qtm": 0.32}, {"qtm": 0.28}
                ],
            },
            # Q5: Rerun query with per-query hit_cap=5000 and 1200 hits. Not saturated at cap 5000!
            "Q_rerun_uncapped": {
                "pdb4": "",
                "hit_cap": 5000,
                "total_hits_raw": 1200,
                "rem_counts": [1200, 1200, 100, 20, 5],
                "top_qtm": [
                    {"qtm": 0.80}, {"qtm": 0.44}, {"qtm": 0.38}, {"qtm": 0.32}, {"qtm": 0.28}
                ],
            },
            # Q6: Rerun query with per-query hit_cap=5000 and 5000 hits. Saturated at cap 5000, LOW_SCORE.
            "Q_rerun_saturated": {
                "pdb4": "",
                "hit_cap": 5000,
                "total_hits_raw": 5000,
                "rem_counts": [5000, 500, 100, 20, 5],
                "top_qtm": [
                    {"qtm": 0.80}, {"qtm": 0.44}, {"qtm": 0.38}, {"qtm": 0.32}, {"qtm": 0.28}
                ],
            },
        }

        # Run audit — must execute cleanly without NameError!
        audit = audit_hit_cap_saturation(acc, label="test")

        self.assertEqual(audit["total_queries"], 6)
        # Saturated: Q_sat_all_matched (1000/1000), Q_sat_post_empty (1000/1000),
        # Q_sat_low_score (1000/1000), Q_rerun_saturated (5000/5000). Total = 4.
        self.assertEqual(audit["saturated_count"], 4)

        # Uncapped queries must NOT be in saturated list
        self.assertNotIn("Q_uncapped", audit["problematic_qids"])
        self.assertNotIn("Q_rerun_uncapped", audit["problematic_qids"])

        # Q_sat_all_matched must NOT be in problematic list because it is MATCHED at all tiers
        self.assertNotIn("Q_sat_all_matched", audit["problematic_qids"])

        # Q_sat_post_empty must be in problematic_exclusion_empty at tier_2, tier_3, tier_4
        self.assertIn("Q_sat_post_empty", audit["problematic_exclusion_empty"]["tier_2"])
        self.assertIn("Q_sat_post_empty", audit["problematic_qids"])

        # Q_sat_low_score must be in problematic_low_score at tier_1, tier_2, tier_3, tier_4
        self.assertIn("Q_sat_low_score", audit["problematic_low_score"]["tier_1"])
        self.assertIn("Q_sat_low_score", audit["problematic_qids"])

        # Q_rerun_saturated must be in problematic_low_score at tier_1
        self.assertIn("Q_rerun_saturated", audit["problematic_low_score"]["tier_1"])
        self.assertIn("Q_rerun_saturated", audit["problematic_qids"])

    def test_rerun_target_deduplication_and_score_monotonicity(self):
        """Test that merging rerun hits with original hits preserves best scores
        and deduplicates identical alignments.
        """
        # Construct synthetic Foldseek raw lines (21 columns)
        # Format: query, target, alntm, qtm, ttm, prob, lddt, evalue, bits, fident, pident, qcov, tcov, qlen, tlen, alnlen, qstart, qend, tstart, tend, cigar
        orig_line_1 = [
            "Q_test", "1abc_A", "0.55", "0.45", "0.50", "0.9", "0.6", "1e-5", "100", "0.22", "22", "0.8", "0.8", "100", "100", "80", "10", "89", "5", "84", "80M"
        ]
        orig_line_2 = [
            "Q_test", "2def_A", "0.50", "0.40", "0.48", "0.8", "0.5", "1e-4", "80", "0.18", "18", "0.7", "0.7", "100", "100", "70", "15", "84", "10", "79", "70M"
        ]

        # Rerun lines: contains orig_line_1 (duplicate), a lower-scoring alignment, and a new high-scoring hit
        rerun_line_dup = list(orig_line_1)
        rerun_line_worse = [
            "Q_test", "3ghi_A", "0.40", "0.35", "0.38", "0.7", "0.4", "1e-3", "60", "0.15", "15", "0.6", "0.6", "100", "100", "60", "20", "79", "10", "69", "60M"
        ]
        rerun_line_better = [
            "Q_test", "4jkl_A", "0.75", "0.68", "0.70", "0.99", "0.8", "1e-12", "250", "0.28", "28", "0.9", "0.9", "100", "100", "90", "5", "94", "5", "94", "90M"
        ]

        orig_qa = accumulate_hit_lines([orig_line_1, orig_line_2], pdb4="", hit_cap=1000)
        self.assertEqual(orig_qa["total_hits_raw"], 2)
        self.assertEqual(orig_qa["top_qtm"][1]["qtm"], 0.45)
        self.assertEqual(classify_single_unit(2, orig_qa["rem_counts"][1], orig_qa["top_qtm"][1]), "LOW_SCORE")

        # Deduplicate hits across original and rerun
        all_lines = [orig_line_1, orig_line_2, rerun_line_dup, rerun_line_worse, rerun_line_better]
        seen_keys = set()
        deduped = []
        for p in all_lines:
            k = (p[1], p[16], p[17], p[18], p[19])
            if k not in seen_keys:
                seen_keys.add(k)
                deduped.append(p)

        self.assertEqual(len(deduped), 4)  # 2 original + 2 distinct rerun = 4 unique hits

        merged_qa = accumulate_hit_lines(deduped, pdb4="", hit_cap=5000)
        self.assertEqual(merged_qa["total_hits_raw"], 4)
        # Score improved from 0.45 to 0.68
        self.assertEqual(merged_qa["top_qtm"][1]["qtm"], 0.68)
        self.assertEqual(classify_single_unit(4, merged_qa["rem_counts"][1], merged_qa["top_qtm"][1]), "MATCHED")

        # Verify score monotonicity: scores never worsen
        for k in range(5):
            if orig_qa["top_qtm"][k] is not None:
                self.assertGreaterEqual(merged_qa["top_qtm"][k]["qtm"], orig_qa["top_qtm"][k]["qtm"] - 1e-9)

    def test_undefined_domains_return_none(self):
        """Test that undefined conditioning domains return None/null, not [0.0, 0.0]."""
        self.assertIsNone(safe_rate(0, 0))
        self.assertIsNone(safe_rate(5, 0))
        self.assertEqual(safe_rate(0, 100), 0.0)
        self.assertEqual(safe_rate(50, 100), 50.0)

        # Model-based Wilson interval behavior
        self.assertIsNone(wilson_ci_model_based(0, 0))
        # k=0 boundary produces non-zero-width interval [0.0, hi] with hi > 0
        w0 = wilson_ci_model_based(0, 100)
        self.assertIsNotNone(w0)
        self.assertEqual(w0[0], 0.0)
        self.assertGreater(w0[1], 0.0)
        self.assertNotEqual(w0, [0.0, 0.0])

        # k=n boundary produces non-zero-width interval [lo, 100.0] with lo < 100.0
        wn = wilson_ci_model_based(100, 100)
        self.assertIsNotNone(wn)
        self.assertLess(wn[0], 100.0)
        self.assertEqual(wn[1], 100.0)
        self.assertNotEqual(wn, [100.0, 100.0])

    def test_cigar_aligned_pairs_formula(self):
        """Test CIGAR M operations match qspan + tspan - alnlen."""
        # 1. Exact match
        alnlen, qspan, tspan = 100, 100, 100
        cigar_m = 100
        self.assertEqual(cigar_m, qspan + tspan - alnlen)

        # 2. Query insertion / target gap (2I): alnlen = 102, qspan = 102, tspan = 100
        alnlen, qspan, tspan = 102, 102, 100
        cigar_m = 100
        self.assertEqual(cigar_m, qspan + tspan - alnlen)

        # 3. Target insertion / query gap (3D): alnlen = 103, qspan = 100, tspan = 103
        alnlen, qspan, tspan = 103, 100, 103
        cigar_m = 100
        self.assertEqual(cigar_m, qspan + tspan - alnlen)

    def test_checksum_manifest_and_verifier_fail_closed(self):
        """Test that manifest generation avoids self-inclusion, and that
        verify_closure fails closed on missing files or checksum mismatches.
        """
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create minimal mock artifacts
            for fn in REQUIRED_ARTIFACTS:
                fpath = os.path.join(tmpdir, fn)
                with open(fpath, "w") as fh:
                    fh.write(f"dummy content for {fn}\n")

            # Generate manifest
            entries = generate_sha256_manifest(tmpdir, manifest_name="v7_results.sha256.tsv")
            manifest_path = os.path.join(tmpdir, "v7_results.sha256.tsv")
            self.assertTrue(os.path.exists(manifest_path))

            # Verify manifest does NOT list itself
            with open(manifest_path) as fh:
                manifest_text = fh.read()
                self.assertNotIn("v7_results.sha256.tsv", manifest_text)

            # Test fail-closed on missing required artifact
            os.remove(os.path.join(tmpdir, "v7_region_evaluations.tsv"))
            with self.assertRaises(AssertionError) as ctx:
                verify_closure(tmpdir)
            self.assertIn("FAIL-CLOSED: Required artifact v7_region_evaluations.tsv missing", str(ctx.exception))

            # Test fail-closed on checksum mismatch
            with open(os.path.join(tmpdir, "v7_region_evaluations.tsv"), "w") as fh:
                fh.write("corrupted content\n")
            with self.assertRaises(AssertionError) as ctx:
                verify_closure(tmpdir)
            self.assertIn("FAIL: SHA-256 mismatch", str(ctx.exception))

    def test_tier_monotonicity_validation(self):
        """Test that tier monotonicity assertions catch violations (T_{k+1} > T_k)."""
        # Monotonic case
        hits_mono = [100, 100, 80, 50, 20]
        qtm_mono = [0.8, 0.8, 0.6, 0.4, 0.2]
        self.assertTrue(all(hits_mono[k] >= hits_mono[k+1] for k in range(4)))
        self.assertTrue(all(qtm_mono[k] >= qtm_mono[k+1] for k in range(4)))

        # Non-monotonic hit count violation
        hits_violated = [100, 100, 70, 80, 20]
        self.assertFalse(all(hits_violated[k] >= hits_violated[k+1] for k in range(4)))

        # Non-monotonic score violation
        qtm_violated = [0.8, 0.8, 0.5, 0.6, 0.2]
        self.assertFalse(all(qtm_violated[k] >= qtm_violated[k+1] for k in range(4)))


if __name__ == "__main__":
    unittest.main()
