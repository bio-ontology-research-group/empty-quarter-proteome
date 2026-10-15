#!/usr/bin/env python3
"""01_audit_exact_hits_and_leakage.py — Audit exact Foldseek hit rows, disconnected maxima,
and PDB control leakage across all 2,838 dark candidates and 3,013 controls.

Objectives:
1. Exact-hit rows vs Disconnected Maxima:
   Quantify how often `best_qtm` and `best_alntm` refer to different target structures.
   Compute exact single-target query coverage (qcov = alnlen / qlen) from the top-scoring hit row.
2. Retrieval Decomposition:
   Decompose candidates and controls into:
   - No retrieval (0 returned rows from PDB/AFDB)
   - Retrieved with score failure (all hits TM < 0.5)
   - Retrieved with match (at least one hit TM >= 0.5)
3. PDB Control Overlap & Leakage Audit:
   For every control KF..., extract its parent PDB structure ID from the manifest.
   Identify self-hits (same 4-letter PDB accession).
   Measure control error floor before and after excluding self-hits.
4. Pfam Domain Architecture Intersection:
   Join with per_protein_pfam.tsv to evaluate multi-domain Pfam architecture among
   low-qTM vs fold-matched dark models.
"""

import os
import sys
import csv
import json
import math
from collections import defaultdict

W = "/ibex/scratch/projects/c2014/rob/metagenomics"
OUT_DIR = os.path.join(W, "06_darkmatter/p6_hypotheses_clean/h4_structure/results")

PER_REP = os.path.join(W, "04_structure/00_structure_clean/05_validate/per_rep_structures.tsv")
PER_CTRL = os.path.join(W, "04_structure/00_structure_clean/05_validate/per_control_structures.tsv")
MANIFEST = os.path.join(W, "rebuild_20260729/structural_calibration/current_reference_control_manifest_20260802.tsv")
PDB_HITS = os.path.join(W, "04_structure/00_structure_clean/05_validate/screen/pdb_hits.raw.tsv")
AFDB_HITS = os.path.join(W, "04_structure/00_structure_clean/05_validate/screen/afdb_hits.raw.tsv")
PFAM_PER_PROT = os.path.join(W, "04_structure/pfam_candidates_clean/results/per_protein_pfam.tsv")


def wilson_ci(k, n, z=1.959964):
    if n == 0:
        return (0.0, 0.0)
    p = k / n
    denom = 1 + z * z / n
    center = (p + z * z / (2 * n)) / denom
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / denom
    return (max(0.0, center - half), min(1.0, center + half))


def load_manifest(path):
    controls = {}
    with open(path) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            bid = row["blind_id"]
            pdb_id = row["pdb_structure_id"]
            pdb4 = pdb_id.split("-")[0].lower() if pdb_id else ""
            controls[bid] = {
                "blind_id": bid,
                "target_id": row["target_id"],
                "control_mode": row["control_mode"],
                "pdb_structure_id": pdb_id,
                "pdb4": pdb4,
                "control_length": int(row["control_length"]),
                "fold_mean_plddt": float(row["fold_mean_plddt"]),
            }
    return controls


def load_gated_cohorts():
    candidates = {}
    with open(PER_REP) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["confidence_qualified"] == "1":
                candidates[r["id"]] = {
                    "id": r["id"],
                    "length": int(r["length"]),
                    "mean_plddt": float(r["mean_plddt"]),
                    "best_pdb_qtm": float(r["best_pdb_qtm"]),
                    "best_pdb_alntm": float(r["best_pdb_alntm"]),
                    "best_afdb_qtm": float(r["best_afdb_qtm"]),
                    "best_afdb_alntm": float(r["best_afdb_alntm"]),
                    "conf_dual_novel": int(r["conf_qualified_dual_novel"]),
                }

    controls = {}
    with open(PER_CTRL) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["confidence_qualified"] == "1":
                controls[r["id"]] = {
                    "id": r["id"],
                    "length": int(r["length"]),
                    "mean_plddt": float(r["mean_plddt"]),
                    "best_pdb_qtm": float(r["best_pdb_qtm"]),
                    "best_pdb_alntm": float(r["best_pdb_alntm"]),
                    "best_afdb_qtm": float(r["best_afdb_qtm"]),
                    "best_afdb_alntm": float(r["best_afdb_alntm"]),
                    "conf_dual_novel": int(r["conf_qualified_dual_novel"]),
                }
    return candidates, controls


def stream_hits(hit_path, query_filter=None):
    """Yield parsed rows from raw hit tsv:
    query, target, alntm, qtm, ttm, prob, lddt, evalue, bits, qcov, tcov, qlen, tlen, alnlen
    """
    with open(hit_path) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 14:
                continue
            q = p[0][:-4] if p[0].endswith(".pdb") else p[0]
            if query_filter and q not in query_filter:
                continue
            yield {
                "query": q,
                "target": p[1],
                "alntm": float(p[2]),
                "qtm": float(p[3]),
                "ttm": float(p[4]),
                "prob": float(p[5]),
                "lddt": float(p[6]),
                "evalue": float(p[7]),
                "bits": float(p[8]),
                "qcov": float(p[9]),
                "tcov": float(p[10]),
                "qlen": int(p[11]),
                "tlen": int(p[12]),
                "alnlen": int(p[13]),
            }


def main():
    print("Step 1: Loading manifests and gated cohorts...")
    ctrl_meta = load_manifest(MANIFEST)
    cands, ctrls = load_gated_cohorts()
    print(f"Loaded {len(cands)} gate-passing candidates and {len(ctrls)} gate-passing controls.")

    target_queries = set(cands.keys()) | set(ctrls.keys())

    # Collect per-query hits
    # query -> list of hits from PDB and AFDB
    pdb_hits_by_q = defaultdict(list)
    afdb_hits_by_q = defaultdict(list)

    print("Step 2: Processing PDB raw hits...")
    for h in stream_hits(PDB_HITS, target_queries):
        pdb_hits_by_q[h["query"]].append(h)
    print(f"Queries with at least 1 PDB hit: {len(pdb_hits_by_q)}")

    print("Step 3: Processing AFDB raw hits...")
    for h in stream_hits(AFDB_HITS, target_queries):
        afdb_hits_by_q[h["query"]].append(h)
    print(f"Queries with at least 1 AFDB hit: {len(afdb_hits_by_q)}")

    # Audit Disconnected Maxima & Exact Single-Target Hits
    # For each candidate and control, analyze:
    # 1. Did best_qtm and best_alntm come from the same target in PDB? In AFDB? Across both?
    # 2. What is the single-target coverage of the best hit?
    records = []
    
    # Stratification sets:
    # K_local: qTM < 0.5, alnTM >= 0.5
    # K_dual_none: qTM < 0.5, alnTM < 0.5
    # Q_not_K: qTM >= 0.5
    # Controls: KF...
    
    diff_target_pdb_count = 0
    diff_target_afdb_count = 0
    total_with_multi_pdb_hits = 0
    total_with_multi_afdb_hits = 0

    print("Step 4: Evaluating single-target vs disconnected maxima...")
    for qid in sorted(target_queries):
        is_ctrl = qid.startswith("KF")
        info = ctrls[qid] if is_ctrl else cands[qid]
        phits = pdb_hits_by_q.get(qid, [])
        ahits = afdb_hits_by_q.get(qid, [])
        
        # PDB top hits
        best_pdb_q_row = max(phits, key=lambda x: x["qtm"]) if phits else None
        best_pdb_a_row = max(phits, key=lambda x: x["alntm"]) if phits else None
        same_target_pdb = (best_pdb_q_row["target"] == best_pdb_a_row["target"]) if (best_pdb_q_row and best_pdb_a_row) else None
        if len(phits) > 1:
            total_with_multi_pdb_hits += 1
            if same_target_pdb is False:
                diff_target_pdb_count += 1

        # AFDB top hits
        best_afdb_q_row = max(ahits, key=lambda x: x["qtm"]) if ahits else None
        best_afdb_a_row = max(ahits, key=lambda x: x["alntm"]) if ahits else None
        same_target_afdb = (best_afdb_q_row["target"] == best_afdb_a_row["target"]) if (best_afdb_q_row and best_afdb_a_row) else None
        if len(ahits) > 1:
            total_with_multi_afdb_hits += 1
            if same_target_afdb is False:
                diff_target_afdb_count += 1

        # Global best hit by bits or TM across all hits
        all_hits = phits + ahits
        top_hit_bits = max(all_hits, key=lambda x: x["bits"]) if all_hits else None
        top_hit_alntm = max(all_hits, key=lambda x: x["alntm"]) if all_hits else None
        top_hit_qtm = max(all_hits, key=lambda x: x["qtm"]) if all_hits else None

        records.append({
            "id": qid,
            "is_control": int(is_ctrl),
            "length": info["length"],
            "mean_plddt": info["mean_plddt"],
            "n_pdb_hits": len(phits),
            "n_afdb_hits": len(ahits),
            "has_retrieval": int(len(all_hits) > 0),
            "same_target_pdb": int(same_target_pdb) if same_target_pdb is not None else -1,
            "same_target_afdb": int(same_target_afdb) if same_target_afdb is not None else -1,
            "top_alntm_target": top_hit_alntm["target"] if top_hit_alntm else "NONE",
            "top_alntm_val": top_hit_alntm["alntm"] if top_hit_alntm else 0.0,
            "top_alntm_qcov": top_hit_alntm["qcov"] if top_hit_alntm else 0.0,
            "top_alntm_alnlen": top_hit_alntm["alnlen"] if top_hit_alntm else 0,
            "top_qtm_target": top_hit_qtm["target"] if top_hit_qtm else "NONE",
            "top_qtm_val": top_hit_qtm["qtm"] if top_hit_qtm else 0.0,
            "top_qtm_qcov": top_hit_qtm["qcov"] if top_hit_qtm else 0.0,
            "same_global_target": int(top_hit_alntm["target"] == top_hit_qtm["target"]) if (top_hit_alntm and top_hit_qtm) else -1
        })

    # Summary of disconnected maxima
    print(f"\n--- Disconnected Maxima Audit ---")
    print(f"PDB: {diff_target_pdb_count} / {total_with_multi_pdb_hits} queries with multiple hits "
          f"({100.0 * diff_target_pdb_count / max(1, total_with_multi_pdb_hits):.2f}%) have best_qtm and best_alntm pointing to DIFFERENT targets.")
    print(f"AFDB: {diff_target_afdb_count} / {total_with_multi_afdb_hits} queries with multiple hits "
          f"({100.0 * diff_target_afdb_count / max(1, total_with_multi_afdb_hits):.2f}%) have best_qtm and best_alntm pointing to DIFFERENT targets.")

    # Step 5: Control Leakage Audit (Self-PDB Exclusion)
    print("\nStep 5: Auditing Control Leakage and Self-Hit Exclusion...")
    ctrl_results = []
    
    # We evaluate controls under:
    # 1. Unfiltered (as shipped)
    # 2. Excluding exact PDB self-hits (matching 4-letter PDB code of parent structure)
    # 3. PDB-only vs AFDB-only vs Dual
    for cid, cinfo in ctrls.items():
        meta = ctrl_meta.get(cid, {})
        pdb4 = meta.get("pdb4", "")
        phits = pdb_hits_by_q.get(cid, [])
        ahits = afdb_hits_by_q.get(cid, [])

        # Unfiltered maxima
        u_pdb_qtm = max([h["qtm"] for h in phits], default=0.0)
        u_pdb_atm = max([h["alntm"] for h in phits], default=0.0)
        u_afdb_qtm = max([h["qtm"] for h in ahits], default=0.0)
        u_afdb_atm = max([h["alntm"] for h in ahits], default=0.0)
        u_dual_qtm_novel = int(u_pdb_qtm < 0.5 and u_afdb_qtm < 0.5)
        u_dual_atm_novel = int(u_pdb_atm < 0.5 and u_afdb_atm < 0.5)

        # Identify self-hits in PDB
        self_hits_pdb = [h for h in phits if h["target"][:4].lower() == pdb4] if pdb4 else []
        non_self_hits_pdb = [h for h in phits if h["target"][:4].lower() != pdb4] if pdb4 else phits

        # Excluded maxima (PDB self-hits removed)
        ex_pdb_qtm = max([h["qtm"] for h in non_self_hits_pdb], default=0.0)
        ex_pdb_atm = max([h["alntm"] for h in non_self_hits_pdb], default=0.0)
        # Note: AFDB entries are AlphaFold UniProt models (e.g. AF-P12345-F1).
        # We also evaluate PDB-alone failure rate to isolate exact crystallographic leakage.
        ex_pdb_only_qtm_novel = int(ex_pdb_qtm < 0.5)
        u_pdb_only_qtm_novel = int(u_pdb_qtm < 0.5)

        ctrl_results.append({
            "blind_id": cid,
            "pdb_parent": meta.get("pdb_structure_id", ""),
            "pdb4": pdb4,
            "control_mode": meta.get("control_mode", ""),
            "length": cinfo["length"],
            "plddt": cinfo["mean_plddt"],
            "n_pdb_hits": len(phits),
            "n_self_hits_pdb": len(self_hits_pdb),
            "top_self_qtm": max([h["qtm"] for h in self_hits_pdb], default=0.0),
            "top_non_self_qtm": ex_pdb_qtm,
            "top_non_self_alntm": ex_pdb_atm,
            "unfiltered_pdb_qtm": u_pdb_qtm,
            "unfiltered_dual_qtm_novel": u_dual_qtm_novel,
            "unfiltered_dual_atm_novel": u_dual_atm_novel,
            "unfiltered_pdb_only_qtm_novel": u_pdb_only_qtm_novel,
            "excluded_pdb_only_qtm_novel": ex_pdb_only_qtm_novel,
            "excluded_dual_qtm_novel": int(ex_pdb_qtm < 0.5 and u_afdb_qtm < 0.5),
        })

    # Summary of control self-hit leakage
    n_ctrl = len(ctrl_results)
    controls_with_self_hit = sum(1 for r in ctrl_results if r["n_self_hits_pdb"] > 0)
    print(f"Total gate-passing controls: {n_ctrl}")
    print(f"Controls with exact PDB 4-letter self-hit in PDB search: {controls_with_self_hit} / {n_ctrl} "
          f"({100.0 * controls_with_self_hit / n_ctrl:.2f}%)")

    u_dual_q_fail = sum(r["unfiltered_dual_qtm_novel"] for r in ctrl_results)
    u_pdb_q_fail = sum(r["unfiltered_pdb_only_qtm_novel"] for r in ctrl_results)
    ex_pdb_q_fail = sum(r["excluded_pdb_only_qtm_novel"] for r in ctrl_results)
    ex_dual_q_fail = sum(r["excluded_dual_qtm_novel"] for r in ctrl_results)

    print(f"Unfiltered dual qTM < 0.5 floor: {u_dual_q_fail} / {n_ctrl} ({100.0 * u_dual_q_fail / n_ctrl:.2f}%)")
    print(f"Unfiltered PDB-only qTM < 0.5 floor: {u_pdb_q_fail} / {n_ctrl} ({100.0 * u_pdb_q_fail / n_ctrl:.2f}%)")
    print(f"PDB-only qTM < 0.5 AFTER SELF-PDB EXCLUSION: {ex_pdb_q_fail} / {n_ctrl} ({100.0 * ex_pdb_q_fail / n_ctrl:.2f}%)")
    print(f"Dual qTM < 0.5 with PDB self-exclusion (AFDB untouched): {ex_dual_q_fail} / {n_ctrl} ({100.0 * ex_dual_q_fail / n_ctrl:.2f}%)")

    # Step 6: Candidate Retrieval vs Score Decomposition
    print("\nStep 6: Retrieval Decomposition...")
    # Candidates
    cand_no_ret = sum(1 for qid in cands if len(pdb_hits_by_q.get(qid, [])) == 0 and len(afdb_hits_by_q.get(qid, [])) == 0)
    cand_ret_dual_qtm_fail = sum(1 for qid, c in cands.items() if (len(pdb_hits_by_q.get(qid, [])) > 0 or len(afdb_hits_by_q.get(qid, [])) > 0) and c["conf_dual_novel"] == 1)
    cand_ret_match = len(cands) - cand_no_ret - cand_ret_dual_qtm_fail

    # Controls
    ctrl_no_ret = sum(1 for qid in ctrls if len(pdb_hits_by_q.get(qid, [])) == 0 and len(afdb_hits_by_q.get(qid, [])) == 0)
    ctrl_ret_dual_qtm_fail = sum(1 for qid, c in ctrls.items() if (len(pdb_hits_by_q.get(qid, [])) > 0 or len(afdb_hits_by_q.get(qid, [])) > 0) and c["conf_dual_novel"] == 1)
    ctrl_ret_match = len(ctrls) - ctrl_no_ret - ctrl_ret_dual_qtm_fail

    print(f"Candidates (N={len(cands)}):")
    print(f"  Zero Retrieval (no hits in either ref): {cand_no_ret} ({100.0 * cand_no_ret / len(cands):.2f}%)")
    print(f"  Retrieved Low-Score (all hits qTM < 0.5): {cand_ret_dual_qtm_fail} ({100.0 * cand_ret_dual_qtm_fail / len(cands):.2f}%)")
    print(f"  Retrieved Match (at least one hit qTM >= 0.5): {cand_ret_match} ({100.0 * cand_ret_match / len(cands):.2f}%)")

    print(f"Controls (N={len(ctrls)}):")
    print(f"  Zero Retrieval (no hits in either ref): {ctrl_no_ret} ({100.0 * ctrl_no_ret / len(ctrls):.2f}%)")
    print(f"  Retrieved Low-Score (all hits qTM < 0.5): {ctrl_ret_dual_qtm_fail} ({100.0 * ctrl_ret_dual_qtm_fail / len(ctrls):.2f}%)")
    print(f"  Retrieved Match (at least one hit qTM >= 0.5): {ctrl_ret_match} ({100.0 * ctrl_ret_match / len(ctrls):.2f}%)")

    # Step 7: Pfam Domain Architecture Intersection
    print("\nStep 7: Intersecting with Pfam Domain Architecture...")
    pfam_by_prot = {}
    if os.path.exists(PFAM_PER_PROT):
        with open(PFAM_PER_PROT) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                pfam_by_prot[r["protein"]] = r

    k_local_ids = [qid for qid, c in cands.items() if c["conf_dual_novel"] == 1 and (c["best_pdb_alntm"] >= 0.5 or c["best_afdb_alntm"] >= 0.5)]
    k_dual_none_ids = [qid for qid, c in cands.items() if c["conf_dual_novel"] == 1 and c["best_pdb_alntm"] < 0.5 and c["best_afdb_alntm"] < 0.5]
    q_not_k_ids = [qid for qid, c in cands.items() if c["conf_dual_novel"] == 0]

    def pfam_stats(qids, label):
        n_tot = len(qids)
        with_pfam = [qid for qid in qids if qid in pfam_by_prot and int(pfam_by_prot[qid]["n_dom_resolved"]) > 0]
        multi_dom = [qid for qid in with_pfam if int(pfam_by_prot[qid]["n_dom_resolved"]) > 1]
        with_gap80 = [qid for qid in qids if qid in pfam_by_prot and pfam_by_prot[qid]["has_gap80"] == "1"]
        print(f"Set: {label} (N={n_tot}):")
        print(f"  With Pfam domain: {len(with_pfam)} ({100.0 * len(with_pfam) / n_tot:.2f}%)")
        print(f"  Multi-domain Pfam: {len(multi_dom)} ({100.0 * len(multi_dom) / n_tot:.2f}%)")
        print(f"  Has >=80aa uncovered gap: {len(with_gap80)} ({100.0 * len(with_gap80) / n_tot:.2f}%)")

    pfam_stats(k_local_ids, "K_local (qTM < 0.5, alnTM >= 0.5)")
    pfam_stats(k_dual_none_ids, "K_dual_none (qTM < 0.5, alnTM < 0.5)")
    pfam_stats(q_not_k_ids, "Q_not_K (qTM >= 0.5)")

    # Save outputs to results directory
    os.makedirs(OUT_DIR, exist_ok=True)
    
    summary = {
        "candidate_counts": {
            "total_conf_qualified": len(cands),
            "k_local_alntm_ge_05": len(k_local_ids),
            "k_dual_none_alntm_lt_05": len(k_dual_none_ids),
            "q_not_k_fold_matched": len(q_not_k_ids),
            "zero_retrieval_count": cand_no_ret,
            "retrieved_low_score_count": cand_ret_dual_qtm_fail,
            "retrieved_match_count": cand_ret_match,
        },
        "control_leakage_audit": {
            "total_conf_controls": n_ctrl,
            "controls_with_pdb_self_hit": controls_with_self_hit,
            "pct_controls_with_pdb_self_hit": round(100.0 * controls_with_self_hit / n_ctrl, 2),
            "unfiltered_dual_qtm_floor": round(100.0 * u_dual_q_fail / n_ctrl, 2),
            "unfiltered_pdb_only_floor": round(100.0 * u_pdb_q_fail / n_ctrl, 2),
            "excluded_pdb_only_floor": round(100.0 * ex_pdb_q_fail / n_ctrl, 2),
            "excluded_dual_floor_pdb_self_excl": round(100.0 * ex_dual_q_fail / n_ctrl, 2),
            "ci95_excluded_pdb_only_floor": [round(100 * x, 2) for x in wilson_ci(ex_pdb_q_fail, n_ctrl)],
        },
        "disconnected_maxima_audit": {
            "multi_hit_pdb_queries": total_with_multi_pdb_hits,
            "diff_target_pdb_queries": diff_target_pdb_count,
            "pct_diff_target_pdb": round(100.0 * diff_target_pdb_count / max(1, total_with_multi_pdb_hits), 2),
            "multi_hit_afdb_queries": total_with_multi_afdb_hits,
            "diff_target_afdb_queries": diff_target_afdb_count,
            "pct_diff_target_afdb": round(100.0 * diff_target_afdb_count / max(1, total_with_multi_afdb_hits), 2),
        }
    }

    with open(os.path.join(OUT_DIR, "audit_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)

    with open(os.path.join(OUT_DIR, "control_leakage_audit.tsv"), "w") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(ctrl_results[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(ctrl_results)

    with open(os.path.join(OUT_DIR, "exact_hit_vs_maxima.tsv"), "w") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(records[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(records)

    print(f"\nAudit complete. Outputs written to {OUT_DIR}")


if __name__ == "__main__":
    main()
