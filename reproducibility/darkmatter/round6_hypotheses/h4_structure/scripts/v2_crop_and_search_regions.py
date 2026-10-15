#!/usr/bin/env python3
"""v2_crop_and_search_regions.py — Outcome-independent confident region cropping,
Foldseek search, and reference-conditional recovery calibration.

Methodological specifications (addressing ROOT_AUDIT_NOTES.md):
1. Unit of Analysis: Confident structural regions (contiguous segments >=50 aa with
   smoothed pLDDT >= 0.70). These are explicitly defined as confident regions,
   NOT asserted as validated evolutionary/functional domains.
2. Cropping: ATOM records corresponding to [reg_start, reg_end] are extracted into
   dedicated PDB files, preserving exact residue sequence numbers, coordinates, and chain IDs.
3. Cropped Search: Foldseek easy-search is executed directly on the cropped region models
   against PDB (pdb_20260730/foldseek_pdb), obtaining region-specific qTM, alnTM,
   alignment coordinates (qstart, qend), and sequence identity (fident).
4. Recovery Estimand: Within-model change in retrieval and score between whole-chain
   normalization and region-specific normalization:
   - For low-whole-query-TM models (qTM_wc < 0.5): does cropping restore match (qTM_crop >= 0.5)?
   - Decomposes outcomes into:
     * Zero-retrieval (0 hits returned)
     * Region recovery (qTM_crop >= 0.5 on at least one confident region)
     * Persistent low score (all hits qTM_crop < 0.5 on all confident regions)
5. Multi-Level Control Exclusion:
   - Unfiltered (including parent PDB structure)
   - Source-PDB exclusion (excluding same 4-letter PDB code)
   - Sequence-identity exclusion sensitivity (excluding hits with fident >= 0.30, 0.25, 0.20)
   - Explicitly documents residual homologs, search truncation, and limitations.
6. Internal Dark Comparator: Uses fold-matched dark models (Q_not_K: qTM_wc >= 0.5)
   as an internal comparator alongside PDB controls.
"""

import os
import sys
import csv
import json
import math
import random
import argparse
import subprocess
from collections import defaultdict

W = "/ibex/scratch/projects/c2014/rob/metagenomics"
H4_DIR = os.path.join(W, "06_darkmatter/p6_hypotheses_clean/h4_structure")
OUT_DIR = os.path.join(H4_DIR, "results")

PER_REP = os.path.join(W, "04_structure/00_structure_clean/05_validate/per_rep_structures.tsv")
PER_CTRL = os.path.join(W, "04_structure/00_structure_clean/05_validate/per_control_structures.tsv")
MANIFEST = os.path.join(W, "rebuild_20260729/structural_calibration/current_reference_control_manifest_20260802.tsv")
QUERY_PDB_DIR = os.path.join(W, "04_structure/00_structure_clean/05_validate/query_pdb")
PDB_DB = os.path.join(W, "rebuild_20260729/structural_calibration/references/pdb_20260730/foldseek_pdb")

PILOT_SEED = 20261013


def norm_plddt(val):
    """Normalize pLDDT to 0.0-1.0 scale if on 0-100 scale."""
    return val / 100.0 if val > 1.5 else val


def parse_pdb_ca(pdb_path):
    """Extract CA atoms: list of dicts with res_num, res_name, plddt, coord."""
    ca_atoms = []
    if not os.path.exists(pdb_path):
        return ca_atoms
    with open(pdb_path) as fh:
        for line in fh:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                try:
                    rnum = int(line[22:26].strip())
                    rname = line[17:20].strip()
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    b = float(line[60:66].strip())
                    ca_atoms.append({
                        "res_num": rnum,
                        "res_name": rname,
                        "plddt": norm_plddt(b),
                        "coord": (x, y, z)
                    })
                except (ValueError, IndexError):
                    continue
    return ca_atoms


def identify_confident_regions(ca_atoms, window=15, plddt_thresh=0.70, min_len=50, max_bridge=5):
    """Outcome-independent identification of confident regions:
    - Symmetric rolling-window smoothed pLDDT (window=15)
    - Threshold: smoothed pLDDT >= 0.70
    - Bridge small dips (<= 5 residues with mean pLDDT >= 0.55)
    - Minimum length threshold: >= 50 residues
    Returns: list of dicts {reg_id, start_res, end_res, length, mean_plddt}
    """
    L = len(ca_atoms)
    if L == 0:
        return []

    plddts = [a["plddt"] for a in ca_atoms]
    hw = window // 2
    smoothed = []
    for i in range(L):
        s = max(0, i - hw)
        e = min(L, i + hw + 1)
        smoothed.append(sum(plddts[s:e]) / (e - s))

    is_conf = [s >= plddt_thresh for s in smoothed]

    # Initial blocks
    blocks = []
    in_block = False
    start_idx = 0
    for i, c in enumerate(is_conf):
        if c and not in_block:
            in_block = True
            start_idx = i
        elif not c and in_block:
            in_block = False
            blocks.append((start_idx, i - 1))
    if in_block:
        blocks.append((start_idx, L - 1))

    # Bridge small dips
    merged = []
    for b in blocks:
        if not merged:
            merged.append(list(b))
        else:
            prev_e = merged[-1][1]
            gap = b[0] - prev_e - 1
            if gap <= max_bridge:
                gap_mean = sum(plddts[prev_e + 1: b[0]]) / max(1, gap)
                if gap_mean >= 0.55:
                    merged[-1][1] = b[1]
                else:
                    merged.append(list(b))
            else:
                merged.append(list(b))

    # Apply length threshold >= 50 aa
    regions = []
    for s_idx, e_idx in merged:
        reg_len = e_idx - s_idx + 1
        if reg_len >= min_len:
            reg_mean_pl = sum(plddts[s_idx: e_idx + 1]) / reg_len
            regions.append({
                "reg_id": len(regions) + 1,
                "start_idx": s_idx,
                "end_idx": e_idx,
                "start_res": ca_atoms[s_idx]["res_num"],
                "end_res": ca_atoms[e_idx]["res_num"],
                "length": reg_len,
                "mean_plddt": round(reg_mean_pl, 4)
            })
    return regions


def write_cropped_pdb(src_pdb, dst_pdb, start_res, end_res, source_id, reg_id):
    """Extract all ATOM records with residue numbers in [start_res, end_res]
    preserving exact residue numbers, chain IDs, and coordinates.
    """
    n_atoms = 0
    with open(src_pdb) as fh_in, open(dst_pdb, "w") as fh_out:
        fh_out.write(f"REMARK 999 CROPPED CONFIDENT REGION\n")
        fh_out.write(f"REMARK 999 SOURCE {source_id} REGION {reg_id} RESIDUES {start_res} TO {end_res}\n")
        for line in fh_in:
            if line.startswith("ATOM"):
                try:
                    rnum = int(line[22:26].strip())
                    if start_res <= rnum <= end_res:
                        fh_out.write(line)
                        n_atoms += 1
                except ValueError:
                    continue
        fh_out.write("TER\nEND\n")
    return n_atoms


def load_cohorts():
    """Load candidates and controls from frozen gate tables."""
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

    ctrl_meta = {}
    with open(MANIFEST) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            ctrl_meta[r["blind_id"]] = r

    controls = {}
    with open(PER_CTRL) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["confidence_qualified"] == "1":
                cid = r["id"]
                cm = ctrl_meta.get(cid, {})
                pdb_id = cm.get("pdb_structure_id", "")
                pdb4 = pdb_id.split("-")[0].lower() if pdb_id else ""
                controls[cid] = {
                    "id": cid,
                    "length": int(r["length"]),
                    "mean_plddt": float(r["mean_plddt"]),
                    "best_pdb_qtm": float(r["best_pdb_qtm"]),
                    "best_pdb_alntm": float(r["best_pdb_alntm"]),
                    "best_afdb_qtm": float(r["best_afdb_qtm"]),
                    "best_afdb_alntm": float(r["best_afdb_alntm"]),
                    "conf_dual_novel": int(r["conf_qualified_dual_novel"]),
                    "pdb_parent": pdb_id,
                    "pdb4": pdb4,
                    "control_mode": cm.get("control_mode", "")
                }
    return candidates, controls


def select_cohort(mode="pilot"):
    """Select cohort: either deterministic pilot (N=300) or full cohort (N=5,851)."""
    candidates, controls = load_cohorts()

    if mode == "full":
        cohort = {}
        for qid, c in candidates.items():
            if c["conf_dual_novel"] == 1 and (c["best_pdb_alntm"] >= 0.5 or c["best_afdb_alntm"] >= 0.5):
                st = "K_local"
            elif c["conf_dual_novel"] == 1:
                st = "K_dual_none"
            else:
                st = "Q_not_K"
            c["stratum"] = st
            cohort[qid] = c
        for cid, ctrl in controls.items():
            ctrl["stratum"] = "Controls"
            cohort[cid] = ctrl
        return cohort

    # Pilot selection (N=300)
    random.seed(PILOT_SEED)
    k_local = [c for c in candidates.values() if c["conf_dual_novel"] == 1 and (c["best_pdb_alntm"] >= 0.5 or c["best_afdb_alntm"] >= 0.5)]
    k_dual_none = [c for c in candidates.values() if c["conf_dual_novel"] == 1 and c["best_pdb_alntm"] < 0.5 and c["best_afdb_alntm"] < 0.5]
    q_not_k = [c for c in candidates.values() if c["conf_dual_novel"] == 0]
    ctrl_list = list(controls.values())

    k_local.sort(key=lambda x: x["id"])
    k_dual_none.sort(key=lambda x: x["id"])
    q_not_k.sort(key=lambda x: x["id"])
    ctrl_list.sort(key=lambda x: x["id"])

    # Stratify K_local by length quartiles (25 each)
    bins = [(100, 200), (200, 300), (300, 400), (400, 10000)]
    k_local_sample = []
    for lo, hi in bins:
        sub = [r for r in k_local if lo <= r["length"] < hi]
        k_local_sample.extend(random.sample(sub, min(25, len(sub))))

    k_dual_sample = random.sample(k_dual_none, min(50, len(k_dual_none)))
    q_not_k_sample = random.sample(q_not_k, min(75, len(q_not_k)))
    ctrl_sample = random.sample(ctrl_list, min(75, len(ctrl_list)))

    pilot = {}
    for r in k_local_sample:
        r["stratum"] = "K_local"
        pilot[r["id"]] = r
    for r in k_dual_sample:
        r["stratum"] = "K_dual_none"
        pilot[r["id"]] = r
    for r in q_not_k_sample:
        r["stratum"] = "Q_not_K"
        pilot[r["id"]] = r
    for r in ctrl_sample:
        r["stratum"] = "Controls"
        pilot[r["id"]] = r

    return pilot


def main():
    parser = argparse.ArgumentParser(description="Crop confident regions and evaluate recovery.")
    parser.add_argument("--mode", choices=["pilot", "full"], default="pilot", help="Cohort size: pilot (N=300) or full (N=5,851)")
    parser.add_argument("--threads", type=int, default=16, help="Foldseek search threads")
    args = parser.parse_args()

    print(f"=== Step 1: Loading Cohort (Mode: {args.mode}) ===")
    cohort = select_cohort(args.mode)
    print(f"Selected {len(cohort)} total models across strata:")
    strata_counts = defaultdict(int)
    for m in cohort.values():
        strata_counts[m["stratum"]] += 1
    for st, count in sorted(strata_counts.items()):
        print(f"  {st}: {count}")

    # Step 2: Crop Outcome-Independent Confident Regions
    print("\n=== Step 2: Cropping Outcome-Independent Confident Regions (>=50 aa, pLDDT >= 0.70) ===")
    crops_dir = os.path.join(H4_DIR, "crops")
    os.makedirs(crops_dir, exist_ok=True)

    model_region_manifest = []
    region_pdb_files = []
    models_with_no_conf_region = 0

    for qid, m in sorted(cohort.items()):
        src_pdb = os.path.join(QUERY_PDB_DIR, f"{qid}.pdb")
        ca = parse_pdb_ca(src_pdb)
        if not ca:
            print(f"Warning: could not parse CA atoms for {qid}")
            continue

        regions = identify_confident_regions(ca, window=15, plddt_thresh=0.70, min_len=50, max_bridge=5)
        if not regions:
            models_with_no_conf_region += 1
            model_region_manifest.append({
                "model_id": qid,
                "stratum": m["stratum"],
                "whole_chain_length": len(ca),
                "whole_chain_mean_plddt": m["mean_plddt"],
                "n_confident_regions": 0,
                "region_id": 0,
                "reg_start": 0,
                "reg_end": 0,
                "reg_len": 0,
                "reg_mean_plddt": 0.0,
                "cropped_pdb": "NONE"
            })
            continue

        for reg in regions:
            crop_filename = f"{qid}_reg{reg['reg_id']}.pdb"
            crop_path = os.path.join(crops_dir, crop_filename)
            write_cropped_pdb(src_pdb, crop_path, reg["start_res"], reg["end_res"], qid, reg["reg_id"])
            region_pdb_files.append(crop_path)
            model_region_manifest.append({
                "model_id": qid,
                "stratum": m["stratum"],
                "whole_chain_length": len(ca),
                "whole_chain_mean_plddt": m["mean_plddt"],
                "n_confident_regions": len(regions),
                "region_id": reg["reg_id"],
                "reg_start": reg["start_res"],
                "reg_end": reg["end_res"],
                "reg_len": reg["length"],
                "reg_mean_plddt": reg["mean_plddt"],
                "cropped_pdb": crop_filename
            })

    print(f"Total cropped confident regions created: {len(region_pdb_files)} (from {len(cohort) - models_with_no_conf_region} models)")
    print(f"Models with 0 confident regions (>=50 aa, pLDDT >= 0.70): {models_with_no_conf_region}")

    # Write region manifest
    os.makedirs(OUT_DIR, exist_ok=True)
    with open(os.path.join(OUT_DIR, "v2_cropped_regions_manifest.tsv"), "w") as fh:
        w = csv.DictWriter(fh, fieldnames=list(model_region_manifest[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(model_region_manifest)

    # Step 3: Foldseek Search of Cropped Regions Against PDB
    print("\n=== Step 3: Foldseek Search of Cropped Regions against PDB ===")
    crop_fs_out = os.path.join(OUT_DIR, "v2_crop_pdb_hits.tsv")
    tmp_out = crop_fs_out + ".tmp"
    tmp_scratch = os.path.join(H4_DIR, "fs_crop_tmp")
    os.makedirs(tmp_scratch, exist_ok=True)

    fmt = "query,target,alntmscore,qtmscore,ttmscore,prob,lddt,evalue,bits,fident,pident,qcov,tcov,qlen,tlen,alnlen,qstart,qend,tstart,tend"
    cmd = [
        "foldseek", "easy-search",
        crops_dir, PDB_DB, tmp_out, tmp_scratch,
        "--alignment-type", "1", "-e", "10", "-s", "9.5", "--max-seqs", "50",
        "--threads", str(args.threads),
        "--format-output", fmt
    ]
    print(f"Executing: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)
    os.rename(tmp_out, crop_fs_out)
    print(f"Cropped search complete: {crop_fs_out}")

    # Step 4: Parse Cropped Hits and Within-Model Recovery
    print("\n=== Step 4: Evaluating Within-Model Recovery and Sequence-Identity Exclusions ===")
    hits_by_crop = defaultdict(list)
    with open(crop_fs_out) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 20:
                continue
            cq = p[0][:-4] if p[0].endswith(".pdb") else p[0]
            hits_by_crop[cq].append({
                "query": cq,
                "target": p[1],
                "alntm": float(p[2]),
                "qtm": float(p[3]),
                "ttm": float(p[4]),
                "prob": float(p[5]),
                "lddt": float(p[6]),
                "evalue": float(p[7]),
                "bits": float(p[8]),
                "fident": float(p[9]),
                "pident": float(p[10]),
                "qcov": float(p[11]),
                "tcov": float(p[12]),
                "qlen": int(p[13]),
                "tlen": int(p[14]),
                "alnlen": int(p[15]),
                "qstart": int(p[16]),
                "qend": int(p[17]),
                "qalnlen": int(p[17]) - int(p[16]) + 1
            })

    # Evaluate each region under multiple exclusion tiers:
    # Tier 0: Unfiltered
    # Tier 1: Source-PDB excluded (for controls)
    # Tier 2: Sequence identity < 0.30 (sensitivity)
    # Tier 3: Sequence identity < 0.25
    # Tier 4: Sequence identity < 0.20
    region_evaluations = []
    model_evaluations = defaultdict(list)

    for rm in model_region_manifest:
        if rm["cropped_pdb"] == "NONE":
            continue
        crop_id = f"{rm['model_id']}_reg{rm['region_id']}"
        hits = hits_by_crop.get(crop_id, [])
        m_info = cohort[rm["model_id"]]
        is_ctrl = rm["stratum"] == "Controls"
        pdb4 = m_info.get("pdb4", "")

        # Tier 0: Unfiltered top hit
        t0_top = max(hits, key=lambda x: x["qtm"]) if hits else None

        # Tier 1: Source PDB excluded
        t1_hits = [h for h in hits if h["target"][:4].lower() != pdb4] if (is_ctrl and pdb4) else hits
        t1_top = max(t1_hits, key=lambda x: x["qtm"]) if t1_hits else None

        # Tier 2: Sequence identity < 0.30
        t2_hits = [h for h in t1_hits if h["fident"] < 0.30]
        t2_top = max(t2_hits, key=lambda x: x["qtm"]) if t2_hits else None

        # Tier 3: Sequence identity < 0.25
        t3_hits = [h for h in t1_hits if h["fident"] < 0.25]
        t3_top = max(t3_hits, key=lambda x: x["qtm"]) if t3_hits else None

        # Tier 4: Sequence identity < 0.20
        t4_hits = [h for h in t1_hits if h["fident"] < 0.20]
        t4_top = max(t4_hits, key=lambda x: x["qtm"]) if t4_hits else None

        reg_eval = {
            "model_id": rm["model_id"],
            "region_id": rm["region_id"],
            "crop_id": crop_id,
            "stratum": rm["stratum"],
            "reg_len": rm["reg_len"],
            "reg_mean_plddt": rm["reg_mean_plddt"],
            "total_hits_raw": len(hits),
            "has_retrieval": int(len(hits) > 0),
            # Tier 0 Unfiltered
            "t0_top_target": t0_top["target"] if t0_top else "NONE",
            "t0_qtm": t0_top["qtm"] if t0_top else 0.0,
            "t0_alntm": t0_top["alntm"] if t0_top else 0.0,
            "t0_fident": t0_top["fident"] if t0_top else 0.0,
            "t0_qalnlen": t0_top["qalnlen"] if t0_top else 0,
            "t0_qcov": round(t0_top["qalnlen"] / rm["reg_len"], 4) if t0_top else 0.0,
            "t0_matched": int(t0_top["qtm"] >= 0.5) if t0_top else 0,
            # Tier 1 Source-PDB Excluded
            "t1_top_target": t1_top["target"] if t1_top else "NONE",
            "t1_qtm": t1_top["qtm"] if t1_top else 0.0,
            "t1_alntm": t1_top["alntm"] if t1_top else 0.0,
            "t1_fident": t1_top["fident"] if t1_top else 0.0,
            "t1_qcov": round(t1_top["qalnlen"] / rm["reg_len"], 4) if t1_top else 0.0,
            "t1_matched": int(t1_top["qtm"] >= 0.5) if t1_top else 0,
            # Tier 2 SeqId < 0.30
            "t2_qtm": t2_top["qtm"] if t2_top else 0.0,
            "t2_matched": int(t2_top["qtm"] >= 0.5) if t2_top else 0,
            # Tier 3 SeqId < 0.25
            "t3_qtm": t3_top["qtm"] if t3_top else 0.0,
            "t3_matched": int(t3_top["qtm"] >= 0.5) if t3_top else 0,
            # Tier 4 SeqId < 0.20
            "t4_qtm": t4_top["qtm"] if t4_top else 0.0,
            "t4_matched": int(t4_top["qtm"] >= 0.5) if t4_top else 0,
        }
        region_evaluations.append(reg_eval)
        model_evaluations[rm["model_id"]].append(reg_eval)

    # Within-model summary
    model_summaries = []
    for qid, m in sorted(cohort.items()):
        reg_evals = model_evaluations.get(qid, [])
        n_regs = len(reg_evals)
        
        # Whole-chain reference metrics (from frozen tables)
        wc_qtm = m["best_pdb_qtm"]
        wc_alntm = m["best_pdb_alntm"]
        wc_matched = int(wc_qtm >= 0.5)

        # Region-level recovery (Tier 1: source-PDB excluded)
        if n_regs > 0:
            reg_t1_max_qtm = max(r["t1_qtm"] for r in reg_evals)
            reg_t1_matched_count = sum(r["t1_matched"] for r in reg_evals)
            reg_has_retrieval = any(r["has_retrieval"] == 1 for r in reg_evals)
            reg_t2_max_qtm = max(r["t2_qtm"] for r in reg_evals)
            reg_t2_matched_count = sum(r["t2_matched"] for r in reg_evals)
        else:
            reg_t1_max_qtm = 0.0
            reg_t1_matched_count = 0
            reg_has_retrieval = False
            reg_t2_max_qtm = 0.0
            reg_t2_matched_count = 0

        # Outcome-conditional transition category:
        # Category 1: "Recovered by Cropping" (wc_qtm < 0.5 -> reg_t1_max_qtm >= 0.5)
        # Category 2: "Persistent Low Score" (wc_qtm < 0.5 and reg_t1_max_qtm < 0.5 and reg_has_retrieval)
        # Category 3: "Zero Retrieval" (no hits returned on regions)
        # Category 4: "Whole-Chain Baseline Match" (wc_qtm >= 0.5)
        if wc_qtm < 0.5:
            if not reg_has_retrieval:
                trans_cat = "ZERO_RETRIEVAL"
            elif reg_t1_max_qtm >= 0.5:
                trans_cat = "RECOVERED_AFTER_CROPPING"
            else:
                trans_cat = "PERSISTENT_LOW_SCORE"
        else:
            trans_cat = "WHOLE_CHAIN_MATCHED"

        model_summaries.append({
            "model_id": qid,
            "stratum": m["stratum"],
            "whole_chain_length": m["length"],
            "whole_chain_plddt": m["mean_plddt"],
            "whole_chain_pdb_qtm": wc_qtm,
            "whole_chain_pdb_alntm": wc_alntm,
            "whole_chain_matched": wc_matched,
            "n_confident_regions": n_regs,
            "region_has_retrieval": int(reg_has_retrieval),
            "region_t1_max_qtm": round(reg_t1_max_qtm, 4),
            "region_t1_matched_count": reg_t1_matched_count,
            "region_t2_max_qtm_seqid30": round(reg_t2_max_qtm, 4),
            "region_t2_matched_count_seqid30": reg_t2_matched_count,
            "transition_category": trans_cat
        })

    # Summary Statistics across strata
    print("\n=== SUMMARY RESULTS: Within-Model Region Recovery ===")
    summary_by_stratum = {}
    for st in ["K_local", "K_dual_none", "Q_not_K", "Controls"]:
        sub = [m for m in model_summaries if m["stratum"] == st]
        if not sub:
            continue
        n_tot = len(sub)
        n_rec = sum(1 for m in sub if m["transition_category"] == "RECOVERED_AFTER_CROPPING")
        n_pers = sum(1 for m in sub if m["transition_category"] == "PERSISTENT_LOW_SCORE")
        n_zero = sum(1 for m in sub if m["transition_category"] == "ZERO_RETRIEVAL")
        n_wc = sum(1 for m in sub if m["transition_category"] == "WHOLE_CHAIN_MATCHED")

        # Region match rate (Tier 1: source PDB excluded)
        st_regs = [r for r in region_evaluations if r["stratum"] == st]
        n_regs_tot = len(st_regs)
        n_regs_t1_match = sum(r["t1_matched"] for r in st_regs)
        n_regs_t2_match = sum(r["t2_matched"] for r in st_regs)
        n_regs_t3_match = sum(r["t3_matched"] for r in st_regs)

        summary_by_stratum[st] = {
            "model_count": n_tot,
            "transition_recovered_after_cropping": n_rec,
            "pct_recovered": round(100.0 * n_rec / n_tot, 2),
            "transition_persistent_low_score": n_pers,
            "pct_persistent_low_score": round(100.0 * n_pers / n_tot, 2),
            "transition_zero_retrieval": n_zero,
            "pct_zero_retrieval": round(100.0 * n_zero / n_tot, 2),
            "whole_chain_matched": n_wc,
            "total_regions": n_regs_tot,
            "regions_t1_source_excl_matched": n_regs_t1_match,
            "pct_regions_t1_matched": round(100.0 * n_regs_t1_match / max(1, n_regs_tot), 2),
            "regions_t2_seqid30_matched": n_regs_t2_match,
            "pct_regions_t2_matched": round(100.0 * n_regs_t2_match / max(1, n_regs_tot), 2),
            "regions_t3_seqid25_matched": n_regs_t3_match,
            "pct_regions_t3_matched": round(100.0 * n_regs_t3_match / max(1, n_regs_tot), 2),
        }

        print(f"\nStratum: {st} (N={n_tot} models, N_regions={n_regs_tot}):")
        print(f"  Models Recovered after Cropping (qTM_wc < 0.5 -> qTM_crop >= 0.5): {n_rec} ({100.0*n_rec/n_tot:.1f}%)")
        print(f"  Models with Persistent Low Score (all regions qTM_crop < 0.5): {n_pers} ({100.0*n_pers/n_tot:.1f}%)")
        print(f"  Models with Zero Retrieval on Regions: {n_zero} ({100.0*n_zero/n_tot:.1f}%)")
        print(f"  Regions Matched (Tier 1 Source-PDB Excluded): {n_regs_t1_match} / {n_regs_tot} ({100.0*n_regs_t1_match/max(1, n_regs_tot):.1f}%)")
        print(f"  Regions Matched (Tier 2 SeqId < 30% Sensitivity): {n_regs_t2_match} / {n_regs_tot} ({100.0*n_regs_t2_match/max(1, n_regs_tot):.1f}%)")
        print(f"  Regions Matched (Tier 3 SeqId < 25% Sensitivity): {n_regs_t3_match} / {n_regs_tot} ({100.0*n_regs_t3_match/max(1, n_regs_tot):.1f}%)")

    # Save output tables
    with open(os.path.join(OUT_DIR, "v2_region_evaluations.tsv"), "w") as fh:
        w = csv.DictWriter(fh, fieldnames=list(region_evaluations[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(region_evaluations)

    with open(os.path.join(OUT_DIR, "v2_model_summaries.tsv"), "w") as fh:
        w = csv.DictWriter(fh, fieldnames=list(model_summaries[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(model_summaries)

    summary_json = {
        "analysis_version": "v2_cropped_regions",
        "cohort_mode": args.mode,
        "pilot_seed": PILOT_SEED if args.mode == "pilot" else None,
        "region_definition": {
            "min_length_aa": 50,
            "min_smoothed_plddt": 0.70,
            "smoothing_window_aa": 15,
            "max_bridge_gap_aa": 5
        },
        "stratum_summaries": summary_by_stratum
    }
    with open(os.path.join(OUT_DIR, "v2_recovery_summary.json"), "w") as fh:
        json.dump(summary_json, fh, indent=2)

    print(f"\nAll v2 results successfully written to {OUT_DIR}")


if __name__ == "__main__":
    main()
