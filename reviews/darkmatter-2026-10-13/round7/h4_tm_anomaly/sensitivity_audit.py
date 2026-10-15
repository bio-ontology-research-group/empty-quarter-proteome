#!/usr/bin/env python3
"""sensitivity_audit.py — H4 Foldseek TM-Score Range Anomaly Audit & Read-Only Exclusion Sensitivity.

Performs:
1. Complete full-row numerical and coordinate census across all 4 raw TSV result files.
2. Exact extraction and categorization of all anomalous rows (alntmscore > 1.00001, etc.).
3. Dual-arm recomputation of H4 structural hit accumulation, adaptive rerun merging,
   model classifications, cropped region evaluations, whole/crop transitions, and 8-strata x 5-tier census:
   - Arm A (Baseline): All rows included.
   - Arm B (Exclusion Sensitivity): Exclude any row where alntmscore > 1.00001 (or any invalid score).
4. Direct verification of Arm A against published v7 results (zero-discrepancy assertion).
5. Granular comparison of Arm B vs Arm A across all per-query top scores, unit statuses,
   crop statuses, transitions, and stratum census counts.
"""

import csv
import collections
import hashlib
import json
import math
import os
import re
import sys
import time
from collections import defaultdict
from pathlib import Path

# Insert scripts_recovery_v3 to reuse vetted H4 postprocessing functions
H4_DIR = Path("/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p7_analyses_clean/h4_structure")
SCRIPT_DIR = H4_DIR / "scripts_recovery_v3"
sys.path.insert(0, str(SCRIPT_DIR))

from v7_crop_and_search_regions import (
    load_full_cohort,
    classify_single_unit,
    classify_model_crops,
    classify_transition,
    summarize_stratum_census,
    SUB_STRATA_ORDER,
    COMBINED_STRATA_MAP,
    MAX_SEQS,
    RERUN_MAX_SEQS,
    CIGAR_RE,
    sha256_file,
)

RESULTS_DIR = H4_DIR / "results"
OUT_DIR = H4_DIR / "tm_anomaly"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def parse_and_accumulate_dual_stream(hit_tsv, query_to_pdb4, initial_cap=MAX_SEQS):
    """Single-pass streaming parser that accumulates into both Arm A (all rows)
    and Arm B (filtered rows excluding invalid scores), while computing complete
    column statistics across 100% of rows.
    """
    print(f"Streaming and auditing: {hit_tsv.name}...")
    acc_A = {}
    acc_B = {}
    for qid, pdb4 in query_to_pdb4.items():
        base_qa = {
            "pdb4": pdb4,
            "hit_cap": initial_cap,
            "total_hits_raw": 0,
            "max_fident_raw": 0.0,
            "hits_fident_ge_20": 0,
            "hits_fident_ge_25": 0,
            "hits_fident_ge_30": 0,
            "rem_counts": [0, 0, 0, 0, 0],
            "top_qtm": [None, None, None, None, None],
            "top_alntm": [None, None, None, None, None],
        }
        acc_A[qid] = json.loads(json.dumps(base_qa))
        acc_B[qid] = json.loads(json.dumps(base_qa))

    stats = {
        "file": hit_tsv.name,
        "rows": 0,
        "queries_with_hits": set(),
        "ranges": {
            "alntm": [float("inf"), float("-inf")],
            "qtm": [float("inf"), float("-inf")],
            "ttm": [float("inf"), float("-inf")],
            "prob": [float("inf"), float("-inf")],
            "lddt": [float("inf"), float("-inf")],
            "evalue": [float("inf"), float("-inf")],
            "bits": [float("inf"), float("-inf")],
            "fident": [float("inf"), float("-inf")],
            "pident": [float("inf"), float("-inf")],
            "qcov": [float("inf"), float("-inf")],
            "tcov": [float("inf"), float("-inf")],
        },
        "anomalies": {
            "alntm_gt_1": 0,
            "alntm_lt_0": 0,
            "qtm_gt_1": 0,
            "qtm_lt_0": 0,
            "ttm_gt_1": 0,
            "ttm_lt_0": 0,
            "fident_out_of_range": 0,
            "lddt_out_of_range": 0,
            "evalue_lt_0": 0,
            "cigar_length_mismatch": 0,
            "cigar_m_mismatch": 0,
            "coordinate_invalid": 0,
            "non_finite_numeric": 0,
            "total_invalid_score_rows": 0,
        },
    }
    anomalous_rows = []

    with open(hit_tsv, "r", encoding="utf-8") as fh:
        for line_idx, line in enumerate(fh, start=1):
            p = line.rstrip("\r\n").split("\t")
            if len(p) < 21:
                continue
            stats["rows"] += 1
            q = p[0][:-4] if p[0].endswith(".pdb") else p[0]
            target = p[1]
            stats["queries_with_hits"].add(q)

            try:
                alntm = float(p[2])
                qtm = float(p[3])
                ttm = float(p[4])
                prob = float(p[5])
                lddt = float(p[6])
                evalue = float(p[7])
                bits = float(p[8])
                fident = float(p[9])
                pident = float(p[10])
                fs_qcov = float(p[11])
                fs_tcov = float(p[12])
                qlen = int(p[13])
                tlen = int(p[14])
                alnlen = int(p[15])
                qstart = int(p[16])
                qend = int(p[17])
                tstart = int(p[18])
                tend = int(p[19])
                cigar = p[20]
            except ValueError:
                stats["anomalies"]["non_finite_numeric"] += 1
                continue

            # Check finiteness
            for val_name, v in [
                ("alntm", alntm), ("qtm", qtm), ("ttm", ttm), ("prob", prob),
                ("lddt", lddt), ("evalue", evalue), ("bits", bits),
                ("fident", fident), ("pident", pident), ("qcov", fs_qcov), ("tcov", fs_tcov)
            ]:
                if not math.isfinite(v):
                    stats["anomalies"]["non_finite_numeric"] += 1
                else:
                    if v < stats["ranges"][val_name][0]: stats["ranges"][val_name][0] = v
                    if v > stats["ranges"][val_name][1]: stats["ranges"][val_name][1] = v

            # Range checks
            row_has_invalid_score = False
            if alntm > 1.00001:
                stats["anomalies"]["alntm_gt_1"] += 1
                row_has_invalid_score = True
            elif alntm < 0.0:
                stats["anomalies"]["alntm_lt_0"] += 1
                row_has_invalid_score = True

            if qtm > 1.00001:
                stats["anomalies"]["qtm_gt_1"] += 1
                row_has_invalid_score = True
            elif qtm < 0.0:
                stats["anomalies"]["qtm_lt_0"] += 1
                row_has_invalid_score = True

            if ttm > 1.00001:
                stats["anomalies"]["ttm_gt_1"] += 1
                row_has_invalid_score = True
            elif ttm < 0.0:
                stats["anomalies"]["ttm_lt_0"] += 1
                row_has_invalid_score = True

            if not (0.0 <= fident <= 1.00001):
                stats["anomalies"]["fident_out_of_range"] += 1
            if not (0.0 <= lddt <= 1.00001):
                stats["anomalies"]["lddt_out_of_range"] += 1
            if evalue < 0.0:
                stats["anomalies"]["evalue_lt_0"] += 1

            # Coordinates and CIGAR
            qspan = qend - qstart + 1
            tspan = tend - tstart + 1
            if not (1 <= qstart <= qend <= qlen and 1 <= tstart <= tend <= tlen):
                stats["anomalies"]["coordinate_invalid"] += 1

            cigar_m = 0
            cigar_tot = 0
            for length_str, op in CIGAR_RE.findall(cigar):
                lval = int(length_str)
                cigar_tot += lval
                if op == "M":
                    cigar_m += lval

            expected_m = qspan + tspan - alnlen
            if cigar_tot != alnlen:
                stats["anomalies"]["cigar_length_mismatch"] += 1
            if cigar_m != expected_m:
                stats["anomalies"]["cigar_m_mismatch"] += 1

            if row_has_invalid_score:
                stats["anomalies"]["total_invalid_score_rows"] += 1
                anomalous_rows.append({
                    "file": hit_tsv.name,
                    "line": line_idx,
                    "query": q,
                    "target": target,
                    "alntm": alntm,
                    "qtm": qtm,
                    "ttm": ttm,
                    "fident": fident,
                    "qlen": qlen,
                    "tlen": tlen,
                    "alnlen": alnlen,
                    "qspan": qspan,
                    "tspan": tspan,
                    "cigar": cigar,
                })

            if q not in acc_A:
                continue

            row_dict = {
                "query": q,
                "target": target,
                "alntm": alntm,
                "qtm": qtm,
                "ttm": ttm,
                "evalue": evalue,
                "bits": bits,
                "fident": fident,
                "qlen": qlen,
                "tlen": tlen,
                "alnlen": alnlen,
                "qstart": qstart,
                "qend": qend,
                "tstart": tstart,
                "tend": tend,
                "qspan": qspan,
                "tspan": tspan,
                "aligned_pairs": cigar_m,
                "qspan_cov": round(qspan / qlen, 4),
                "aligned_pair_qcov": round(cigar_m / qlen, 4),
                "cigar": cigar,
            }

            # Feed Arm A (all rows)
            _accumulate_single_row(acc_A[q], row_dict)

            # Feed Arm B (exclude invalid score rows)
            if not row_has_invalid_score:
                _accumulate_single_row(acc_B[q], row_dict)

    stats["queries_with_hits"] = len(stats["queries_with_hits"])
    return acc_A, acc_B, stats, anomalous_rows


def _accumulate_single_row(qa, r):
    fident = r["fident"]
    qtm = r["qtm"]
    alntm = r["alntm"]
    target = r["target"]
    pdb4 = qa["pdb4"]

    qa["total_hits_raw"] += 1
    if fident > qa["max_fident_raw"]:
        qa["max_fident_raw"] = fident
    if fident >= 0.20:
        qa["hits_fident_ge_20"] += 1
    if fident >= 0.25:
        qa["hits_fident_ge_25"] += 1
    if fident >= 0.30:
        qa["hits_fident_ge_30"] += 1

    # Tier 0: Unfiltered
    qa["rem_counts"][0] += 1
    if qa["top_qtm"][0] is None or qtm > qa["top_qtm"][0]["qtm"]:
        qa["top_qtm"][0] = r
    if qa["top_alntm"][0] is None or alntm > qa["top_alntm"][0]["alntm"]:
        qa["top_alntm"][0] = r

    # Tier 1: Source-PDB excluded
    if pdb4 and target[:4].lower() == pdb4:
        return

    qa["rem_counts"][1] += 1
    if qa["top_qtm"][1] is None or qtm > qa["top_qtm"][1]["qtm"]:
        qa["top_qtm"][1] = r
    if qa["top_alntm"][1] is None or alntm > qa["top_alntm"][1]["alntm"]:
        qa["top_alntm"][1] = r

    # Tiers 2, 3, 4: fident < 0.30, 0.25, 0.20
    for t_idx, thresh in ((2, 0.30), (3, 0.25), (4, 0.20)):
        if fident < thresh:
            qa["rem_counts"][t_idx] += 1
            if qa["top_qtm"][t_idx] is None or qtm > qa["top_qtm"][t_idx]["qtm"]:
                qa["top_qtm"][t_idx] = r
            if qa["top_alntm"][t_idx] is None or alntm > qa["top_alntm"][t_idx]["alntm"]:
                qa["top_alntm"][t_idx] = r


def merge_reruns_arm(query_type, prob_qids, acc, query_to_pdb4, orig_tsv, rerun_tsv, exclude_invalid=False):
    """Replicates exact rerun merging from v7_crop_and_search_regions.py with
    optional score exclusion for Arm B.
    """
    prob_set = set(prob_qids)
    if not prob_set:
        return

    orig_lines = defaultdict(list)
    with open(orig_tsv) as fh:
        for line in fh:
            p = line.rstrip("\r\n").split("\t")
            if len(p) >= 21:
                q = p[0][:-4] if p[0].endswith(".pdb") else p[0]
                if q in prob_set:
                    if exclude_invalid:
                        alntm, qtm, ttm = float(p[2]), float(p[3]), float(p[4])
                        if not (0.0 <= alntm <= 1.00001 and 0.0 <= qtm <= 1.00001 and 0.0 <= ttm <= 1.00001):
                            continue
                    orig_lines[q].append(p)

    rerun_lines = defaultdict(list)
    with open(rerun_tsv) as fh:
        for line in fh:
            p = line.rstrip("\r\n").split("\t")
            if len(p) >= 21:
                q = p[0][:-4] if p[0].endswith(".pdb") else p[0]
                if q in prob_set:
                    if exclude_invalid:
                        alntm, qtm, ttm = float(p[2]), float(p[3]), float(p[4])
                        if not (0.0 <= alntm <= 1.00001 and 0.0 <= qtm <= 1.00001 and 0.0 <= ttm <= 1.00001):
                            continue
                    rerun_lines[q].append(p)

    for qid in prob_qids:
        orig_qa = acc[qid]
        merged_lines = []
        seen_full_rows = set()

        for p in orig_lines.get(qid, []):
            rk = tuple(p[:21])
            if rk not in seen_full_rows:
                seen_full_rows.add(rk)
                merged_lines.append(p)

        rerun_raw = rerun_lines.get(qid, [])
        standalone_rerun_count = len(rerun_raw)
        for p in rerun_raw:
            rk = tuple(p[:21])
            if rk not in seen_full_rows:
                seen_full_rows.add(rk)
                merged_lines.append(p)

        # Accumulate merged
        merged_qa = {
            "pdb4": query_to_pdb4.get(qid, ""),
            "hit_cap": RERUN_MAX_SEQS,
            "total_hits_raw": 0,
            "max_fident_raw": 0.0,
            "hits_fident_ge_20": 0,
            "hits_fident_ge_25": 0,
            "hits_fident_ge_30": 0,
            "rem_counts": [0, 0, 0, 0, 0],
            "top_qtm": [None, None, None, None, None],
            "top_alntm": [None, None, None, None, None],
        }
        for p in merged_lines:
            alntm = float(p[2])
            qtm = float(p[3])
            ttm = float(p[4])
            evalue = float(p[7])
            bits = float(p[8])
            fident = float(p[9])
            qlen = int(p[13])
            tlen = int(p[14])
            alnlen = int(p[15])
            qstart = int(p[16])
            qend = int(p[17])
            tstart = int(p[18])
            tend = int(p[19])
            cigar = p[20]
            qspan = qend - qstart + 1
            tspan = tend - tstart + 1
            cigar_m = sum(int(n) for n, op in CIGAR_RE.findall(cigar) if op == "M")

            r_dict = {
                "query": qid,
                "target": p[1],
                "alntm": alntm,
                "qtm": qtm,
                "ttm": ttm,
                "evalue": evalue,
                "bits": bits,
                "fident": fident,
                "qlen": qlen,
                "tlen": tlen,
                "alnlen": alnlen,
                "qstart": qstart,
                "qend": qend,
                "tstart": tstart,
                "tend": tend,
                "qspan": qspan,
                "tspan": tspan,
                "aligned_pairs": cigar_m,
                "qspan_cov": round(qspan / qlen, 4),
                "aligned_pair_qcov": round(cigar_m / qlen, 4),
                "cigar": cigar,
            }
            _accumulate_single_row(merged_qa, r_dict)

        merged_qa["standalone_search_hits_raw"] = standalone_rerun_count

        # Score monotonicity
        for k in range(5):
            orig_top = orig_qa["top_qtm"][k]
            merged_top = merged_qa["top_qtm"][k]
            if orig_top is not None:
                if merged_top is None or orig_top["qtm"] > merged_top["qtm"] + 1e-9:
                    merged_qa["top_qtm"][k] = orig_top

            orig_alntm = orig_qa["top_alntm"][k]
            merged_alntm = merged_qa["top_alntm"][k]
            if orig_alntm is not None:
                if merged_alntm is None or orig_alntm["alntm"] > merged_alntm["alntm"] + 1e-9:
                    merged_qa["top_alntm"][k] = orig_alntm

        acc[qid] = merged_qa


def build_pipeline_tables(full_cohort, manifest_rows, wc_acc, crop_acc):
    """Executes Step 7 and Step 8 of H4 pipeline to build region evaluations,
    model summaries, and stratum census table.
    """
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
        raw_sat = int(ca_info.get("standalone_search_hits_raw", raw_tot) >= ca_info.get("hit_cap", MAX_SEQS))

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
        wc_sat_val = int(wa.get("standalone_search_hits_raw", wc_raw) >= wa.get("hit_cap", MAX_SEQS))
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

    stratum_summaries = {}
    for st in SUB_STRATA_ORDER:
        m_sub = [m for m in model_summaries if m["sub_stratum"] == st]
        r_sub = [r for r in region_evaluations if r["sub_stratum"] == st]
        stratum_summaries[st] = summarize_stratum_census(m_sub, r_sub)

    for comb_name, sub_list in COMBINED_STRATA_MAP.items():
        m_sub = [m for m in model_summaries if m["sub_stratum"] in sub_list]
        r_sub = [r for r in region_evaluations if r["sub_stratum"] in sub_list]
        stratum_summaries[comb_name] = summarize_stratum_census(m_sub, r_sub)

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

    return model_summaries, region_evaluations, stratum_summaries, census_rows


def compare_with_baseline_disk(models_A, regs_A, census_A):
    """Asserts that Arm A reproduces the published v7 tables with 100% fidelity."""
    print("\n=== Validating Arm A against published v7 results ===")
    with open(RESULTS_DIR / "v7_model_summaries.tsv") as fh:
        base_models = list(csv.DictReader(fh, delimiter="\t"))
    with open(RESULTS_DIR / "v7_region_evaluations.tsv") as fh:
        base_regs = list(csv.DictReader(fh, delimiter="\t"))
    with open(RESULTS_DIR / "v7_stratum_census.tsv") as fh:
        base_census = list(csv.DictReader(fh, delimiter="\t"))

    assert len(models_A) == len(base_models), f"Model count mismatch: {len(models_A)} vs {len(base_models)}"
    assert len(regs_A) == len(base_regs), f"Region count mismatch: {len(regs_A)} vs {len(base_regs)}"
    assert len(census_A) == len(base_census), f"Census row count mismatch: {len(census_A)} vs {len(base_census)}"

    model_diffs = 0
    for a, b in zip(models_A, base_models):
        assert set(a.keys()) == set(b.keys()), f"Keys mismatch in model summary: {set(a.keys()) ^ set(b.keys())}"
        for k in b:
            val_a = "" if a[k] is None else str(a[k])
            val_b = "" if b[k] is None else str(b[k])
            if val_a != val_b:
                model_diffs += 1
                if model_diffs < 5:
                    print(f"  Model diff {a['model_id']}.{k}: ArmA='{a[k]}' vs Base='{b[k]}'")

    reg_diffs = 0
    for a, b in zip(regs_A, base_regs):
        assert set(a.keys()) == set(b.keys()), f"Keys mismatch in region evals: {set(a.keys()) ^ set(b.keys())}"
        for k in b:
            val_a = "" if a[k] is None else str(a[k])
            val_b = "" if b[k] is None else str(b[k])
            if val_a != val_b:
                reg_diffs += 1
                if reg_diffs < 5:
                    print(f"  Reg diff {a['crop_id']}.{k}: ArmA='{a[k]}' vs Base='{b[k]}'")

    census_diffs = 0
    for a, b in zip(census_A, base_census):
        assert set(a.keys()) == set(b.keys()), f"Keys mismatch in census: {set(a.keys()) ^ set(b.keys())}"
        for k in b:
            val_a = "" if a[k] is None else str(a[k])
            val_b = "" if b[k] is None else str(b[k])
            if val_a != val_b:
                census_diffs += 1
                if census_diffs < 5:
                    print(f"  Census diff {a['stratum']}.{a['tier']}.{k}: ArmA='{a[k]}' vs Base='{b[k]}'")

    print(f"Baseline fidelity check: model_diffs={model_diffs}, reg_diffs={reg_diffs}, census_diffs={census_diffs}")
    assert model_diffs == 0 and reg_diffs == 0 and census_diffs == 0, "Arm A does not exactly reproduce baseline!"
    print(" Arm A PERFECTLY reproduces 100% of baseline values across all tables.")


def compare_arms(models_A, regs_A, census_A, models_B, regs_B, census_B, wc_acc_A, wc_acc_B, crop_acc_A, crop_acc_B):
    """Compares Arm B (Filtered) against Arm A (Baseline) in complete detail."""
    print("\n=== Comparing Arm B (Exclusion Sensitivity) vs Arm A (Baseline) ===")
    diffs = {
        "wc_top_qtm_changes": [],
        "wc_top_alntm_changes": [],
        "wc_status_changes": [],
        "crop_top_qtm_changes": [],
        "crop_top_alntm_changes": [],
        "crop_status_changes": [],
        "model_crop_status_changes": [],
        "model_transition_changes": [],
        "stratum_census_changes": [],
    }

    # 1. Whole-Chain Query Level
    for qid in sorted(wc_acc_A.keys()):
        qa_A = wc_acc_A[qid]
        qa_B = wc_acc_B[qid]
        for k in range(5):
            top_q_A = qa_A["top_qtm"][k]
            top_q_B = qa_B["top_qtm"][k]
            qtm_A = top_q_A["qtm"] if top_q_A else 0.0
            qtm_B = top_q_B["qtm"] if top_q_B else 0.0
            if abs(qtm_A - qtm_B) > 1e-6:
                diffs["wc_top_qtm_changes"].append({
                    "query": qid, "tier": k, "baseline_qtm": qtm_A, "filtered_qtm": qtm_B, "delta": qtm_B - qtm_A
                })

            top_a_A = qa_A["top_alntm"][k]
            top_a_B = qa_B["top_alntm"][k]
            alntm_A = top_a_A["alntm"] if top_a_A else 0.0
            alntm_B = top_a_B["alntm"] if top_a_B else 0.0
            if abs(alntm_A - alntm_B) > 1e-6:
                diffs["wc_top_alntm_changes"].append({
                    "query": qid, "tier": k, "baseline_alntm": alntm_A, "filtered_alntm": alntm_B, "delta": alntm_B - alntm_A
                })

            stat_A = classify_single_unit(qa_A["total_hits_raw"], qa_A["rem_counts"][k], top_q_A)
            stat_B = classify_single_unit(qa_B["total_hits_raw"], qa_B["rem_counts"][k], top_q_B)
            if stat_A != stat_B:
                diffs["wc_status_changes"].append({
                    "query": qid, "tier": k, "baseline_status": stat_A, "filtered_status": stat_B
                })

    # 2. Cropped Regions Level
    for cid in sorted(crop_acc_A.keys()):
        ca_A = crop_acc_A[cid]
        ca_B = crop_acc_B[cid]
        for k in range(5):
            top_q_A = ca_A["top_qtm"][k]
            top_q_B = ca_B["top_qtm"][k]
            qtm_A = top_q_A["qtm"] if top_q_A else 0.0
            qtm_B = top_q_B["qtm"] if top_q_B else 0.0
            if abs(qtm_A - qtm_B) > 1e-6:
                diffs["crop_top_qtm_changes"].append({
                    "crop_id": cid, "tier": k, "baseline_qtm": qtm_A, "filtered_qtm": qtm_B, "delta": qtm_B - qtm_A
                })

            top_a_A = ca_A["top_alntm"][k]
            top_a_B = ca_B["top_alntm"][k]
            alntm_A = top_a_A["alntm"] if top_a_A else 0.0
            alntm_B = top_a_B["alntm"] if top_a_B else 0.0
            if abs(alntm_A - alntm_B) > 1e-6:
                diffs["crop_top_alntm_changes"].append({
                    "crop_id": cid, "tier": k, "baseline_alntm": alntm_A, "filtered_alntm": alntm_B, "delta": alntm_B - alntm_A
                })

            stat_A = classify_single_unit(ca_A["total_hits_raw"], ca_A["rem_counts"][k], top_q_A)
            stat_B = classify_single_unit(ca_B["total_hits_raw"], ca_B["rem_counts"][k], top_q_B)
            if stat_A != stat_B:
                diffs["crop_status_changes"].append({
                    "crop_id": cid, "tier": k, "baseline_status": stat_A, "filtered_status": stat_B
                })

    # 3. Model Summaries Level
    for mA, mB in zip(models_A, models_B):
        qid = mA["model_id"]
        for k in range(5):
            if mA[f"crop_t{k}_status"] != mB[f"crop_t{k}_status"]:
                diffs["model_crop_status_changes"].append({
                    "model_id": qid, "tier": k, "baseline": mA[f"crop_t{k}_status"], "filtered": mB[f"crop_t{k}_status"]
                })
            if mA[f"trans_t{k}"] != mB[f"trans_t{k}"]:
                diffs["model_transition_changes"].append({
                    "model_id": qid, "tier": k, "baseline": mA[f"trans_t{k}"], "filtered": mB[f"trans_t{k}"]
                })

    # 4. Stratum Census Level
    for cA, cB in zip(census_A, census_B):
        st = cA["stratum"]
        tier = cA["tier"]
        for k in cA:
            if str(cA[k]) != str(cB[k]):
                diffs["stratum_census_changes"].append({
                    "stratum": st, "tier": tier, "field": k, "baseline": cA[k], "filtered": cB[k]
                })

    print(f"Summary of Sensitivity Differences:")
    print(f"  Whole-Chain Top qTM changes:        {len(diffs['wc_top_qtm_changes'])}")
    print(f"  Whole-Chain Top alnTM changes:      {len(diffs['wc_top_alntm_changes'])}")
    print(f"  Whole-Chain Status changes:         {len(diffs['wc_status_changes'])}")
    print(f"  Crop Region Top qTM changes:        {len(diffs['crop_top_qtm_changes'])}")
    print(f"  Crop Region Top alnTM changes:      {len(diffs['crop_top_alntm_changes'])}")
    print(f"  Crop Region Status changes:         {len(diffs['crop_status_changes'])}")
    print(f"  Model-level Crop Status changes:    {len(diffs['model_crop_status_changes'])}")
    print(f"  Model Whole/Crop Transition changes:{len(diffs['model_transition_changes'])}")
    print(f"  Stratum Census Cell changes:        {len(diffs['stratum_census_changes'])}")

    return diffs


def main():
    t0 = time.time()
    print("=== Starting H4 Foldseek TM-Score Anomaly Sensitivity Audit ===")
    full_cohort = load_full_cohort()
    manifest_rows = list(csv.DictReader(open(RESULTS_DIR / "v7_full5851_cropped_regions_manifest.tsv"), delimiter="\t"))
    for row in manifest_rows:
        for key in ("region_id", "reg_start", "reg_end", "reg_len", "n_confident_regions"):
            row[key] = int(row[key])
        row["reg_mean_plddt"] = float(row["reg_mean_plddt"])

    wc_q2pdb4 = {qid: m["pdb4"] for qid, m in full_cohort.items()}
    crop_q2pdb4 = {}
    for r in manifest_rows:
        if r["cropped_pdb"] != "NONE":
            cid = f"{r['model_id']}_reg{r['region_id']}"
            crop_q2pdb4[cid] = full_cohort[r["model_id"]]["pdb4"]

    # Step 1: Scan and accumulate WC hits
    wc_acc_A, wc_acc_B, wc_stats, wc_anom = parse_and_accumulate_dual_stream(
        RESULTS_DIR / "v7_wc_pdb_hits.raw.tsv", wc_q2pdb4, initial_cap=MAX_SEQS
    )

    # Step 2: Scan and accumulate Crop hits
    crop_acc_A, crop_acc_B, crop_stats, crop_anom = parse_and_accumulate_dual_stream(
        RESULTS_DIR / "v7_crop_pdb_hits.raw.tsv", crop_q2pdb4, initial_cap=MAX_SEQS
    )

    # Step 3: Scan and accumulate Rerun hits
    wc_rerun_acc_A, wc_rerun_acc_B, wc_rerun_stats, wc_rerun_anom = parse_and_accumulate_dual_stream(
        RESULTS_DIR / "v7_wc_rerun_hits.raw.tsv", wc_q2pdb4, initial_cap=RERUN_MAX_SEQS
    )
    crop_rerun_acc_A, crop_rerun_acc_B, crop_rerun_stats, crop_rerun_anom = parse_and_accumulate_dual_stream(
        RESULTS_DIR / "v7_crop_rerun_hits.raw.tsv", crop_q2pdb4, initial_cap=RERUN_MAX_SEQS
    )

    # Step 4: Saturation audit for problematic queries
    wc_sat_A = audit_hit_cap_saturation_prob(wc_acc_A)
    crop_sat_A = audit_hit_cap_saturation_prob(crop_acc_A)
    wc_sat_B = audit_hit_cap_saturation_prob(wc_acc_B)
    crop_sat_B = audit_hit_cap_saturation_prob(crop_acc_B)

    # Step 5: Merge reruns
    print("\nMerging adaptive reruns for Arm A and Arm B...")
    merge_reruns_arm("wc", wc_sat_A, wc_acc_A, wc_q2pdb4,
                     RESULTS_DIR / "v7_wc_pdb_hits.raw.tsv", RESULTS_DIR / "v7_wc_rerun_hits.raw.tsv", exclude_invalid=False)
    merge_reruns_arm("crop", crop_sat_A, crop_acc_A, crop_q2pdb4,
                     RESULTS_DIR / "v7_crop_pdb_hits.raw.tsv", RESULTS_DIR / "v7_crop_rerun_hits.raw.tsv", exclude_invalid=False)

    merge_reruns_arm("wc", wc_sat_B, wc_acc_B, wc_q2pdb4,
                     RESULTS_DIR / "v7_wc_pdb_hits.raw.tsv", RESULTS_DIR / "v7_wc_rerun_hits.raw.tsv", exclude_invalid=True)
    merge_reruns_arm("crop", crop_sat_B, crop_acc_B, crop_q2pdb4,
                     RESULTS_DIR / "v7_crop_pdb_hits.raw.tsv", RESULTS_DIR / "v7_crop_rerun_hits.raw.tsv", exclude_invalid=True)

    # Step 6: Build tables
    print("\nBuilding tables for Arm A...")
    models_A, regs_A, strata_A, census_A = build_pipeline_tables(full_cohort, manifest_rows, wc_acc_A, crop_acc_A)

    print("\nBuilding tables for Arm B...")
    models_B, regs_B, strata_B, census_B = build_pipeline_tables(full_cohort, manifest_rows, wc_acc_B, crop_acc_B)

    # Step 7: Baseline fidelity check
    compare_with_baseline_disk(models_A, regs_A, census_A)

    # Step 8: Sensitivity comparison
    diffs = compare_arms(models_A, regs_A, census_A, models_B, regs_B, census_B,
                         wc_acc_A, wc_acc_B, crop_acc_A, crop_acc_B)

    # Step 9: Collate anomaly census
    all_anomalous_rows = wc_anom + crop_anom + wc_rerun_anom + crop_rerun_anom
    print(f"\nTotal invalid-score rows across all files: {len(all_anomalous_rows)}")

    # Add stratum annotations to anomalous rows
    for r in all_anomalous_rows:
        qid = r["query"].split("_reg")[0]
        m = full_cohort.get(qid, {})
        r["sub_stratum"] = m.get("sub_stratum", "UNKNOWN")
        r["is_control"] = m.get("is_control", -1)

    queries_affected = collections.Counter()
    strata_affected = collections.Counter()
    controls_vs_dark = collections.Counter()
    for r in all_anomalous_rows:
        queries_affected[r["query"]] += 1
        strata_affected[r["sub_stratum"]] += 1
        controls_vs_dark["control" if r["is_control"] == 1 else "dark_matter"] += 1

    census_summary = {
        "timestamp_utc_unix": time.time(),
        "elapsed_seconds": round(time.time() - t0, 2),
        "total_files_audited": 4,
        "total_rows_audited": wc_stats["rows"] + crop_stats["rows"] + wc_rerun_stats["rows"] + crop_rerun_stats["rows"],
        "total_invalid_score_rows": len(all_anomalous_rows),
        "invalid_score_fraction": len(all_anomalous_rows) / (wc_stats["rows"] + crop_stats["rows"] + wc_rerun_stats["rows"] + crop_rerun_stats["rows"]),
        "file_breakdown": {
            "v7_wc_pdb_hits.raw.tsv": wc_stats,
            "v7_crop_pdb_hits.raw.tsv": crop_stats,
            "v7_wc_rerun_hits.raw.tsv": wc_rerun_stats,
            "v7_crop_rerun_hits.raw.tsv": crop_rerun_stats,
        },
        "query_impact": {
            "distinct_queries_with_invalid_scores": len(queries_affected),
            "by_cohort": dict(controls_vs_dark),
            "by_stratum": dict(strata_affected),
        },
        "score_violation_types": {
            "alntm_gt_1": sum(s["anomalies"]["alntm_gt_1"] for s in [wc_stats, crop_stats, wc_rerun_stats, crop_rerun_stats]),
            "alntm_lt_0": sum(s["anomalies"]["alntm_lt_0"] for s in [wc_stats, crop_stats, wc_rerun_stats, crop_rerun_stats]),
            "qtm_gt_1": sum(s["anomalies"]["qtm_gt_1"] for s in [wc_stats, crop_stats, wc_rerun_stats, crop_rerun_stats]),
            "qtm_lt_0": sum(s["anomalies"]["qtm_lt_0"] for s in [wc_stats, crop_stats, wc_rerun_stats, crop_rerun_stats]),
            "ttm_gt_1": sum(s["anomalies"]["ttm_gt_1"] for s in [wc_stats, crop_stats, wc_rerun_stats, crop_rerun_stats]),
            "ttm_lt_0": sum(s["anomalies"]["ttm_lt_0"] for s in [wc_stats, crop_stats, wc_rerun_stats, crop_rerun_stats]),
            "coordinate_invalid": sum(s["anomalies"]["coordinate_invalid"] for s in [wc_stats, crop_stats, wc_rerun_stats, crop_rerun_stats]),
            "cigar_length_mismatch": sum(s["anomalies"]["cigar_length_mismatch"] for s in [wc_stats, crop_stats, wc_rerun_stats, crop_rerun_stats]),
            "cigar_m_mismatch": sum(s["anomalies"]["cigar_m_mismatch"] for s in [wc_stats, crop_stats, wc_rerun_stats, crop_rerun_stats]),
        },
    }

    # Write tm_anomaly_census.json
    with open(OUT_DIR / "tm_anomaly_census.json", "w") as fh:
        # Convert sets to ints for json serialization
        for f, s in census_summary["file_breakdown"].items():
            s["queries_with_hits"] = int(s["queries_with_hits"])
        json.dump(census_summary, fh, indent=2)

    # Write tm_anomaly_queries.tsv
    with open(OUT_DIR / "tm_anomaly_queries.tsv", "w") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["file", "line", "query", "target", "alntm", "qtm", "ttm", "fident", "qlen", "tlen", "alnlen", "qspan", "tspan", "is_control", "sub_stratum", "cigar"])
        for r in all_anomalous_rows:
            w.writerow([
                r["file"], r["line"], r["query"], r["target"], r["alntm"], r["qtm"], r["ttm"], r["fident"],
                r["qlen"], r["tlen"], r["alnlen"], r["qspan"], r["tspan"], r["is_control"], r["sub_stratum"], r["cigar"]
            ])

    # Write exclusion_sensitivity_diff.json
    with open(OUT_DIR / "exclusion_sensitivity_diff.json", "w") as fh:
        json.dump(diffs, fh, indent=2)

    # Write exclusion_sensitivity_census.tsv
    with open(OUT_DIR / "exclusion_sensitivity_census.tsv", "w") as fh:
        w = csv.DictWriter(fh, fieldnames=list(census_B[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(census_B)

    # Write SHA256 manifest
    manifest_lines = []
    for out_f in sorted(OUT_DIR.glob("*")):
        if out_f.name != "manifest.sha256" and out_f.is_file():
            h = sha256_file(str(out_f))
            manifest_lines.append(f"{h}  {out_f.name}\n")
    with open(OUT_DIR / "manifest.sha256", "w") as fh:
        fh.writelines(manifest_lines)

    print("\n=== Audit & Sensitivity Run Complete ===")
    print(f"Total time elapsed: {time.time() - t0:.2f} seconds.")


def audit_hit_cap_saturation_prob(acc):
    """Returns list of problematic qids from accumulator dict."""
    saturated = [qid for qid, a in acc.items() if a.get("standalone_search_hits_raw", a["total_hits_raw"]) >= a.get("hit_cap", MAX_SEQS)]
    unm_set = set()
    for qid in saturated:
        a = acc[qid]
        for k in range(5):
            st = classify_single_unit(a["total_hits_raw"], a["rem_counts"][k], a["top_qtm"][k])
            if st in ("POST_EXCLUSION_EMPTY", "LOW_SCORE"):
                unm_set.add(qid)
                break
    return sorted(list(unm_set))


if __name__ == "__main__":
    main()
