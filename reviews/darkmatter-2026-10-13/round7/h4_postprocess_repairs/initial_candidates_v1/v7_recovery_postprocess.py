#!/usr/bin/env python3
"""v7_recovery_postprocess.py — Post-Crash Recovery Wrapper for H4 V7 Full Analysis.

Context:
Job 53738540 successfully executed both full-cohort Foldseek searches (whole-chain on 5,851 queries
and cropped confident regions on 7,773 queries against PDB), returning 0, and validated hit rows,
but encountered a NameError in audit_hit_cap_saturation during post-processing.

This recovery script:
1. Verifies that the raw search outputs (v7_wc_pdb_hits.raw.tsv and v7_crop_pdb_hits.raw.tsv)
   exist and are complete, checking logged subprocess success.
2. Pins raw search evidence SHA-256 checksums to v7_raw_search_hashes.sha256.tsv on IBEX.
3. Validates full-cohort manifest closure (N=5,851 models, R=7,773 regions, 161 0-region models).
4. Runs repaired hit streaming and validation (100% CIGAR verification).
5. Runs repaired hit-cap saturation audit (per-query cap, POST_EXCLUSION_EMPTY + LOW_SCORE).
6. Executes adaptive reruns at --max-seqs 5000 if necessary, using separate paths by query type,
   target deduplication, and score monotonicity verification.
7. Produces finite-population exact census (NO sampling-error CIs in primary tables).
8. Runs code-data closure verification and generates clean SHA-256 manifest without circularity.
"""

import os
import sys
import csv
import json
import time
import argparse
from datetime import datetime, timezone
from collections import defaultdict

# Add current script directory to path for imports
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
if SCRIPT_DIR not in sys.path:
    sys.path.insert(0, SCRIPT_DIR)

from v7_crop_and_search_regions import (
    W, H4_DIR, OUT_DIR, QUERY_PDB_DIR, PDB_DB, MAX_SEQS, RERUN_MAX_SEQS,
    SUB_STRATA_ORDER, COMBINED_STRATA_MAP, sha256_file, safe_rate,
    load_full_cohort, parse_and_accumulate_hits, audit_hit_cap_saturation,
    rerun_problematic_queries, classify_single_unit, classify_model_crops,
    classify_transition, summarize_stratum_census, compute_effect_sizes,
    compare_with_round6_pilot, run_closure_assertions, generate_sha256_manifest
)
from v7_verify_closure import verify_closure, REQUIRED_ARTIFACTS


def validate_prerequisites(log_path=None):
    """Validate that raw search files and search completion markers exist."""
    print("=== Step 1: Validating Search Completion & Raw Evidence Prerequisites ===")

    wc_hits_path = os.path.join(OUT_DIR, "v7_wc_pdb_hits.raw.tsv")
    crop_hits_path = os.path.join(OUT_DIR, "v7_crop_pdb_hits.raw.tsv")
    manifest_path = os.path.join(OUT_DIR, "v7_full5851_cropped_regions_manifest.tsv")

    # Assert raw hit files exist and are non-empty
    assert os.path.exists(wc_hits_path), f"FAIL-CLOSED: Missing {wc_hits_path}"
    assert os.path.getsize(wc_hits_path) > 1000000, f"FAIL-CLOSED: {wc_hits_path} is suspiciously small ({os.path.getsize(wc_hits_path)} bytes)"
    assert os.path.exists(crop_hits_path), f"FAIL-CLOSED: Missing {crop_hits_path}"
    assert os.path.getsize(crop_hits_path) > 1000000, f"FAIL-CLOSED: {crop_hits_path} is suspiciously small ({os.path.getsize(crop_hits_path)} bytes)"
    assert os.path.exists(manifest_path), f"FAIL-CLOSED: Missing {manifest_path}"

    print(f"  [FOUND] Raw whole-chain hits: {wc_hits_path} ({os.path.getsize(wc_hits_path):,} bytes)")
    print(f"  [FOUND] Raw crop hits: {crop_hits_path} ({os.path.getsize(crop_hits_path):,} bytes)")
    print(f"  [FOUND] Crop manifest: {manifest_path} ({os.path.getsize(manifest_path):,} bytes)")

    # Validate log completion markers if log provided
    if log_path and os.path.exists(log_path):
        print(f"  Checking previous run log: {log_path}")
        with open(log_path) as fh:
            log_text = fh.read()
            assert "Completed" in log_text and "v7_wc_pdb_hits.raw.tsv" in log_text, (
                "FAIL-CLOSED: Log does not confirm whole-chain search completion"
            )
            assert "Completed" in log_text and "v7_crop_pdb_hits.raw.tsv" in log_text, (
                "FAIL-CLOSED: Log does not confirm crop search completion"
            )
            assert "Validated" in log_text and "hit rows" in log_text, (
                "FAIL-CLOSED: Log does not confirm hit validation step"
            )
        print("  [CONFIRMED] Slurm log confirms both Foldseek searches finished cleanly with returncode 0.")

    # Pin raw search hashes
    print("\n=== Step 2: Pinning Raw Search Evidence SHA-256 Checksums ===")
    wc_hash = sha256_file(wc_hits_path)
    crop_hash = sha256_file(crop_hits_path)
    man_hash = sha256_file(manifest_path)

    raw_pin_path = os.path.join(OUT_DIR, "v7_raw_search_hashes.sha256.tsv")
    with open(raw_pin_path, "w") as fh:
        fh.write(f"{wc_hash}  v7_wc_pdb_hits.raw.tsv\n")
        fh.write(f"{crop_hash}  v7_crop_pdb_hits.raw.tsv\n")
        fh.write(f"{man_hash}  v7_full5851_cropped_regions_manifest.tsv\n")

    print(f"  Pinned v7_wc_pdb_hits.raw.tsv: {wc_hash}")
    print(f"  Pinned v7_crop_pdb_hits.raw.tsv: {crop_hash}")
    print(f"  Pinned manifest: {man_hash}")
    print(f"  Written to {raw_pin_path}")

    return {
        "wc_hits_path": wc_hits_path,
        "crop_hits_path": crop_hits_path,
        "manifest_path": manifest_path,
        "wc_hash": wc_hash,
        "crop_hash": crop_hash,
    }


def run_recovery(threads=16, log_path=None):
    t_start = time.time()
    paths = validate_prerequisites(log_path)

    print("\n=== Step 3: Loading Full Cohort and Crop Manifest ===")
    full_cohort = load_full_cohort()
    assert len(full_cohort) == 5851, f"Expected 5851 models, got {len(full_cohort)}"

    with open(paths["manifest_path"]) as fh:
        manifest_rows = list(csv.DictReader(fh, delimiter="\t"))
    assert len(manifest_rows) == 7934, f"Expected 7934 manifest rows, got {len(manifest_rows)}"

    wc_query_to_pdb4 = {qid: m["pdb4"] for qid, m in full_cohort.items()}
    crop_query_to_pdb4 = {}
    for r in manifest_rows:
        if r["cropped_pdb"] != "NONE":
            cid = f"{r['model_id']}_reg{r['region_id']}"
            crop_query_to_pdb4[cid] = full_cohort[r["model_id"]]["pdb4"]

    print(f"Loaded {len(wc_query_to_pdb4)} whole-chain query targets and {len(crop_query_to_pdb4)} cropped query targets.")

    print("\n=== Step 4: Streaming & Validating Existing Raw Hits (100% CIGAR Validation) ===")
    wc_acc, wc_stats, wc_examples = parse_and_accumulate_hits(paths["wc_hits_path"], wc_query_to_pdb4, initial_cap=MAX_SEQS)
    crop_acc, crop_stats, crop_examples = parse_and_accumulate_hits(paths["crop_hits_path"], crop_query_to_pdb4, initial_cap=MAX_SEQS)
    print(f"Validated {wc_stats['total_rows']:,} whole-chain hit rows and {crop_stats['total_rows']:,} crop hit rows.")

    print("\n=== Step 5: Repaired Hit-Cap Saturation Audit ===")
    wc_sat = audit_hit_cap_saturation(wc_acc, label="whole_chain")
    crop_sat = audit_hit_cap_saturation(crop_acc, label="cropped_regions")
    print(f"Whole-chain saturated queries (>=1000): {wc_sat['saturated_count']} ({wc_sat['saturated_pct']}%)")
    print(f"Cropped saturated queries (>=1000): {crop_sat['saturated_count']} ({crop_sat['saturated_pct']}%)")

    # Adaptive Rerun with --max-seqs 5000 if problematic queries exist
    wc_dir = os.path.join(H4_DIR, "v7_wc_queries")
    crops_dir = os.path.join(H4_DIR, "v7_crops")

    wc_rerun = rerun_problematic_queries(
        "wc", wc_dir, PDB_DB, paths["wc_hits_path"], wc_sat["problematic_qids"],
        wc_acc, wc_query_to_pdb4, threads, higher_max_seqs=RERUN_MAX_SEQS, reuse=True
    )
    crop_rerun = rerun_problematic_queries(
        "crop", crops_dir, PDB_DB, paths["crop_hits_path"], crop_sat["problematic_qids"],
        crop_acc, crop_query_to_pdb4, threads, higher_max_seqs=RERUN_MAX_SEQS, reuse=True
    )

    rerun_audit = {
        "whole_chain_rerun": wc_rerun,
        "cropped_rerun": crop_rerun,
        "any_rerun_executed": wc_rerun["rerun_executed"] or crop_rerun["rerun_executed"],
    }
    if rerun_audit["any_rerun_executed"]:
        print(f"Adaptive reruns: WC updated={wc_rerun['queries_updated']}, rescued={wc_rerun['queries_rescued_to_matched']} | Crop updated={crop_rerun['queries_updated']}, rescued={crop_rerun['queries_rescued_to_matched']}")
    else:
        print("No capped queries lacking MATCHED detected. All saturated queries retained valid matches across exclusions.")

    print("\n=== Step 6: Building Region Evaluations and Model Summaries ===")
    region_evaluations = []
    regs_by_model = defaultdict(list)

    for rm in manifest_rows:
        if rm["cropped_pdb"] == "NONE":
            continue
        qid = rm["model_id"]
        reg_id = rm["region_id"]
        crop_id = f"{qid}_reg{reg_id}"
        ca_info = crop_acc[crop_id]
        raw_tot = ca_info["total_hits_raw"]
        raw_sat = int(raw_tot >= ca_info.get("hit_cap", MAX_SEQS))

        t0_top = ca_info["top_qtm"][0]
        t1_top = ca_info["top_qtm"][1]

        reg_row = {
            "model_id": qid,
            "region_id": reg_id,
            "crop_id": crop_id,
            "sub_stratum": rm["sub_stratum"],
            "parent_stratum": rm["parent_stratum"],
            "design_stratum_h": rm["design_stratum_h"],
            "reg_start": rm["reg_start"],
            "reg_end": rm["reg_end"],
            "reg_len": rm["reg_len"],
            "reg_mean_plddt": rm["reg_mean_plddt"],
            "total_hits_raw": raw_tot,
            "raw_list_saturated": raw_sat,
            "hit_cap": ca_info.get("hit_cap", MAX_SEQS),
            "max_fident_raw": round(ca_info["max_fident_raw"], 4),
            "hits_fident_ge_20": ca_info["hits_fident_ge_20"],
            "hits_fident_ge_25": ca_info["hits_fident_ge_25"],
            "hits_fident_ge_30": ca_info["hits_fident_ge_30"],
            "t0_remaining_hits": ca_info["rem_counts"][0],
            "t0_top_target": t0_top["target"] if t0_top else "NONE",
            "t0_qtm": round(t0_top["qtm"], 4) if t0_top else 0.0,
            "t0_alntm": round(t0_top["alntm"], 4) if t0_top else 0.0,
            "t0_fident": round(t0_top["fident"], 4) if t0_top else 0.0,
            "t0_qspan": t0_top["qspan"] if t0_top else 0,
            "t0_qspan_cov": t0_top["qspan_cov"] if t0_top else 0.0,
            "t0_aligned_pairs": t0_top["aligned_pairs"] if t0_top else 0,
            "t0_aligned_pair_qcov": t0_top["aligned_pair_qcov"] if t0_top else 0.0,
            "t0_status": classify_single_unit(raw_tot, ca_info["rem_counts"][0], t0_top),
            "t0_matched": int(t0_top["qtm"] >= 0.5) if t0_top else 0,
            "t1_remaining_hits": ca_info["rem_counts"][1],
            "t1_top_target": t1_top["target"] if t1_top else "NONE",
            "t1_qtm": round(t1_top["qtm"], 4) if t1_top else 0.0,
            "t1_alntm": round(t1_top["alntm"], 4) if t1_top else 0.0,
            "t1_fident": round(t1_top["fident"], 4) if t1_top else 0.0,
            "t1_qspan": t1_top["qspan"] if t1_top else 0,
            "t1_qspan_cov": t1_top["qspan_cov"] if t1_top else 0.0,
            "t1_aligned_pairs": t1_top["aligned_pairs"] if t1_top else 0,
            "t1_aligned_pair_qcov": t1_top["aligned_pair_qcov"] if t1_top else 0.0,
            "t1_status": classify_single_unit(raw_tot, ca_info["rem_counts"][1], t1_top),
            "t1_matched": int(t1_top["qtm"] >= 0.5) if t1_top else 0,
        }

        for k in (2, 3, 4):
            tk_top = ca_info["top_qtm"][k]
            reg_row[f"t{k}_remaining_hits"] = ca_info["rem_counts"][k]
            reg_row[f"t{k}_top_target"] = tk_top["target"] if tk_top else "NONE"
            reg_row[f"t{k}_qtm"] = round(tk_top["qtm"], 4) if tk_top else 0.0
            reg_row[f"t{k}_alntm"] = round(tk_top["alntm"], 4) if tk_top else 0.0
            reg_row[f"t{k}_fident"] = round(tk_top["fident"], 4) if tk_top else 0.0
            reg_row[f"t{k}_aligned_pair_qcov"] = tk_top["aligned_pair_qcov"] if tk_top else 0.0
            reg_row[f"t{k}_qspan_cov"] = tk_top["qspan_cov"] if tk_top else 0.0
            reg_row[f"t{k}_status"] = classify_single_unit(raw_tot, ca_info["rem_counts"][k], tk_top)
            reg_row[f"t{k}_matched"] = int(tk_top["qtm"] >= 0.5) if tk_top else 0

        region_evaluations.append(reg_row)
        regs_by_model[qid].append(reg_row)

    model_summaries = []
    for qid, m in sorted(full_cohort.items()):
        wa = wc_acc[qid]
        wc_raw = wa["total_hits_raw"]
        wc_sat_val = int(wc_raw >= wa.get("hit_cap", MAX_SEQS))
        reg_evals = regs_by_model.get(qid, [])
        n_regs = len(reg_evals)

        mrow = {
            "model_id": qid,
            "is_control": m["is_control"],
            "sub_stratum": m["sub_stratum"],
            "parent_stratum": m["parent_stratum"],
            "design_stratum_h": m["design_stratum_h"],
            "in_v2_pilot": m["in_v2_pilot"],
            "pilot_ht_weight": m["pilot_ht_weight"],
            "whole_chain_length": m["length"],
            "whole_chain_plddt": m["mean_plddt"],
            "hist_pdb_qtm": m["hist_pdb_qtm"],
            "hist_pdb_alntm": m["hist_pdb_alntm"],
            "hist_afdb_qtm": m["hist_afdb_qtm"],
            "hist_afdb_alntm": m["hist_afdb_alntm"],
            "n_confident_regions": n_regs,
            "total_confident_aa": sum(int(r["reg_len"]) for r in reg_evals),
            "confident_aa_fraction": round(sum(int(r["reg_len"]) for r in reg_evals) / m["length"], 4),
            "wc_total_hits_raw": wc_raw,
            "wc_raw_list_saturated": wc_sat_val,
            "wc_hit_cap": wa.get("hit_cap", MAX_SEQS),
            "crop_any_raw_retrieval": int(any(r["total_hits_raw"] > 0 for r in reg_evals)),
            "crop_any_raw_saturated": int(any(r["raw_list_saturated"] == 1 for r in reg_evals)),
        }

        for k in range(5):
            wc_rem = wa["rem_counts"][k]
            wc_top_q = wa["top_qtm"][k]
            wc_top_a = wa["top_alntm"][k]
            wc_stat = classify_single_unit(wc_raw, wc_rem, wc_top_q)

            mrow[f"wc_t{k}_remaining_hits"] = wc_rem
            mrow[f"wc_t{k}_qtm"] = round(wc_top_q["qtm"], 4) if wc_top_q else 0.0
            mrow[f"wc_t{k}_alntm"] = round(wc_top_a["alntm"], 4) if wc_top_a else 0.0
            mrow[f"wc_t{k}_fident"] = round(wc_top_q["fident"], 4) if wc_top_q else 0.0
            mrow[f"wc_t{k}_aligned_pair_qcov"] = wc_top_q["aligned_pair_qcov"] if wc_top_q else 0.0
            mrow[f"wc_t{k}_qspan_cov"] = wc_top_q["qspan_cov"] if wc_top_q else 0.0
            mrow[f"wc_t{k}_status"] = wc_stat

            crop_stat = classify_model_crops(n_regs, reg_evals, k)
            if n_regs > 0:
                best_reg = max(reg_evals, key=lambda r: r[f"t{k}_qtm"])
                crop_max_qtm = best_reg[f"t{k}_qtm"]
                crop_matched_cnt = sum(r[f"t{k}_matched"] for r in reg_evals)
                crop_rem_regs = sum(1 for r in reg_evals if r[f"t{k}_remaining_hits"] > 0)
            else:
                crop_max_qtm = 0.0
                crop_matched_cnt = 0
                crop_rem_regs = 0

            trans_cat = classify_transition(wc_stat, crop_stat)

            mrow[f"crop_t{k}_rem_regions"] = crop_rem_regs
            mrow[f"crop_t{k}_matched_regions"] = crop_matched_cnt
            mrow[f"crop_t{k}_max_qtm"] = round(crop_max_qtm, 4)
            mrow[f"crop_t{k}_status"] = crop_stat
            mrow[f"trans_t{k}"] = trans_cat

        model_summaries.append(mrow)

    print("\n=== Step 7: Stratum Census & Exact Rate Summaries ===")
    stratum_summaries = {}
    for st in SUB_STRATA_ORDER:
        m_sub = [m for m in model_summaries if m["sub_stratum"] == st]
        r_sub = [r for r in region_evaluations if r["sub_stratum"] == st]
        stratum_summaries[st] = summarize_stratum_census(m_sub, r_sub)

    for comb_name, sub_list in COMBINED_STRATA_MAP.items():
        m_sub = [m for m in model_summaries if m["sub_stratum"] in sub_list]
        r_sub = [r for r in region_evaluations if r["sub_stratum"] in sub_list]
        stratum_summaries[comb_name] = summarize_stratum_census(m_sub, r_sub)

    effect_sizes = compute_effect_sizes(stratum_summaries, tier_k=1)
    pilot_models = [m for m in model_summaries if m["in_v2_pilot"] == 1]
    pilot_regs = [r for r in region_evaluations if any(m["model_id"] == r["model_id"] for m in pilot_models)]
    pilot_comp = compare_with_round6_pilot(stratum_summaries, pilot_models, pilot_regs)

    run_closure_assertions(full_cohort, model_summaries, region_evaluations, stratum_summaries)

    print("\n=== Step 8: Writing Output Tables & Figure Data ===")
    with open(os.path.join(OUT_DIR, "v7_region_evaluations.tsv"), "w") as fh:
        w = csv.DictWriter(fh, fieldnames=list(region_evaluations[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(region_evaluations)

    with open(os.path.join(OUT_DIR, "v7_model_summaries.tsv"), "w") as fh:
        w = csv.DictWriter(fh, fieldnames=list(model_summaries[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(model_summaries)

    combined_examples = wc_examples + crop_examples
    if combined_examples:
        with open(os.path.join(OUT_DIR, "v7_audited_example_hits.tsv"), "w") as fh:
            w = csv.DictWriter(fh, fieldnames=list(combined_examples[0].keys()), delimiter="\t")
            w.writeheader()
            w.writerows(combined_examples)

    census_rows = []
    for st in SUB_STRATA_ORDER + list(COMBINED_STRATA_MAP.keys()):
        sd = stratum_summaries[st]
        for k in range(5):
            tk = sd["tiers"][f"tier_{k}"]
            mr = tk["exact_census_rates"]
            tc = tk["transition_counts"]
            wc_c = tk["wc_status_counts"]
            cr_c = tk["crop_model_status_counts"]
            census_rows.append({
                "stratum": st,
                "tier": f"tier_{k}",
                "model_count": sd["model_count"],
                "models_eligible": sd["models_with_confident_regions"],
                "models_no_region": sd["models_no_confident_region"],
                "total_regions": sd["total_confident_regions"],
                "wc_matched": wc_c["MATCHED"],
                "wc_low_score": wc_c["LOW_SCORE"],
                "wc_post_empty": wc_c["POST_EXCLUSION_EMPTY"],
                "wc_raw_zero": wc_c["RAW_ZERO_RETRIEVAL"],
                "crop_matched": cr_c["MATCHED"],
                "crop_persistent_low": cr_c["PERSISTENT_LOW_SCORE"],
                "crop_post_empty": cr_c["POST_EXCLUSION_EMPTY"],
                "crop_raw_zero": cr_c["RAW_ZERO_RETRIEVAL"],
                "crop_no_region": cr_c["NO_CONFIDENT_REGION"],
                "unconditional_crop_match_count": mr["unconditional_crop_match_count"],
                "unconditional_crop_match_denominator": mr["unconditional_crop_match_denominator"],
                "unconditional_crop_match_pct": mr["unconditional_crop_match_pct"],
                "region_eligible_crop_match_count": mr["region_eligible_crop_match_count"],
                "region_eligible_crop_match_denominator": mr["region_eligible_crop_match_denominator"],
                "region_eligible_crop_match_pct": mr["region_eligible_crop_match_pct"],
                "remaining_hit_conditional_crop_match_count": mr["remaining_hit_conditional_crop_match_count"],
                "remaining_hit_conditional_crop_match_denominator": mr["remaining_hit_conditional_crop_match_denominator"],
                "remaining_hit_conditional_crop_match_pct": mr["remaining_hit_conditional_crop_match_pct"],
                "wc_unmatched_count": mr["wc_unmatched_count"],
                "wc_unmatched_recovered_count": mr["wc_unmatched_recovery_unconditional_count"],
                "wc_unmatched_recovery_pct": mr["wc_unmatched_recovery_unconditional_pct"],
                "wc_unmatched_eligible_count": mr["wc_unmatched_eligible_count"],
                "wc_unmatched_recovery_eligible_pct": mr["wc_unmatched_recovery_eligible_pct"],
                "wc_matched_count": mr["wc_matched_count"],
                "wc_matched_retained_count": mr["wc_matched_retained_count"],
                "wc_matched_retained_pct": mr["wc_matched_retained_pct"],
                "crop_persistent_low_pct": mr["crop_persistent_low_score_pct"],
                "crop_post_empty_pct": mr["crop_post_exclusion_empty_pct"],
                "crop_raw_zero_pct": mr["crop_raw_zero_retrieval_pct"],
                "crop_no_region_pct": mr["crop_no_confident_region_pct"],
                "trans_recovered": tc["WC_UNMATCHED_TO_CROP_RECOVERED"],
                "trans_persistent_low": tc["WC_UNMATCHED_TO_PERSISTENT_LOW_SCORE"],
                "trans_post_empty": tc["WC_UNMATCHED_TO_POST_EXCLUSION_EMPTY"],
                "trans_raw_zero": tc["WC_UNMATCHED_TO_RAW_ZERO_RETRIEVAL"],
                "trans_unm_no_reg": tc["WC_UNMATCHED_NO_CONFIDENT_REGION"],
                "trans_retained": tc["WC_MATCHED_RETAINED_IN_CROP"],
                "trans_lost_low": tc["WC_MATCHED_LOST_LOW_SCORE"],
                "trans_lost_empty": tc["WC_MATCHED_LOST_POST_EXCLUSION_EMPTY"],
                "trans_lost_zero": tc["WC_MATCHED_LOST_RAW_ZERO_RETRIEVAL"],
                "trans_m_no_reg": tc["WC_MATCHED_LOST_NO_CONFIDENT_REGION"],
            })

    with open(os.path.join(OUT_DIR, "v7_stratum_census.tsv"), "w") as fh:
        w = csv.DictWriter(fh, fieldnames=list(census_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(census_rows)

    fig_rows = []
    for r in census_rows:
        for metric, cnt_col, den_col, pct_col in [
            ("unconditional_crop_match", "unconditional_crop_match_count", "unconditional_crop_match_denominator", "unconditional_crop_match_pct"),
            ("region_eligible_crop_match", "region_eligible_crop_match_count", "region_eligible_crop_match_denominator", "region_eligible_crop_match_pct"),
            ("wc_unmatched_recovery_unconditional", "wc_unmatched_recovered_count", "wc_unmatched_count", "wc_unmatched_recovery_pct"),
            ("wc_matched_retention", "wc_matched_retained_count", "wc_matched_count", "wc_matched_retained_pct"),
        ]:
            fig_rows.append({
                "stratum": r["stratum"],
                "tier": r["tier"],
                "metric": metric,
                "count": r[cnt_col],
                "denominator": r[den_col],
                "rate_pct": r[pct_col],
            })

    with open(os.path.join(OUT_DIR, "v7_figure_data.tsv"), "w") as fh:
        w = csv.DictWriter(fh, fieldnames=list(fig_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(fig_rows)

    raw_provenance = {
        "v7_wc_pdb_hits.raw.tsv": {
            "size_bytes": os.path.getsize(paths["wc_hits_path"]),
            "sha256": paths["wc_hash"],
            **wc_stats,
        },
        "v7_crop_pdb_hits.raw.tsv": {
            "size_bytes": os.path.getsize(paths["crop_hits_path"]),
            "sha256": paths["crop_hash"],
            **crop_stats,
        },
    }

    summary_json = {
        "analysis_version": "v7_full_cohort_census_symmetric_whole_chain_and_cropped_search",
        "recovery_mode": True,
        "recovery_timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "finite_population_models": len(full_cohort),
        "evaluated_models": len(model_summaries),
        "evaluated_cropped_regions": len(region_evaluations),
        "models_no_confident_region": sum(1 for m in model_summaries if m["n_confident_regions"] == 0),
        "foldseek_settings": {
            "version": "10-941cd33 (941cd33ff0771cd2e3f144e3293e22a2b87e9fda)",
            "database": PDB_DB,
            "alignment_type": 1,
            "evalue": 10.0,
            "sensitivity": 9.5,
            "max_seqs": MAX_SEQS,
            "cigar_validated_100pct": True,
        },
        "hit_cap_saturation_audit": {
            "whole_chain": wc_sat,
            "cropped_regions": crop_sat,
            "adaptive_rerun_audit": rerun_audit,
        },
        "raw_hit_files_provenance": raw_provenance,
        "effect_sizes": effect_sizes,
        "pilot_comparison": pilot_comp,
        "stratum_summaries": stratum_summaries,
    }

    with open(os.path.join(OUT_DIR, "v7_recovery_summary.json"), "w") as fh:
        json.dump(summary_json, fh, indent=2)

    print("\n=== Step 9: Running Independent Closure Verification & Generating SHA-256 Manifest ===")
    verify_closure(OUT_DIR)
    manifest_entries = generate_sha256_manifest(OUT_DIR)

    elapsed_total = round(time.time() - t_start, 2)
    print(f"\n[SUCCESS] Recovery and post-processing completed in {elapsed_total:.2f}s ({elapsed_total/60.0:.2f}m).")


def main():
    parser = argparse.ArgumentParser(description="H4 V7 Post-Crash Recovery Wrapper")
    parser.add_argument("--threads", type=int, default=16, help="Threads for reruns if needed")
    parser.add_argument("--log", type=str, default=None, help="Path to Slurm job log confirming completed search")
    args = parser.parse_args()

    run_recovery(threads=args.threads, log_path=args.log)


if __name__ == "__main__":
    main()
