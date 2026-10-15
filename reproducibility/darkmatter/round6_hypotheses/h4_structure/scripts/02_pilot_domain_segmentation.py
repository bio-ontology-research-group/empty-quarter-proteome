#!/usr/bin/env python3
"""02_pilot_domain_segmentation.py — Deterministic pilot of low-query-TM dark models,
fold-matched dark models, and PDB controls under independent domain segmentation.

Key components:
1. Deterministic cohort selection (N=250):
   - 100 K_local (qTM < 0.5, alnTM >= 0.5) stratified by length
   - 25 K_dual_none (qTM < 0.5, alnTM < 0.5)
   - 50 Q_not_K (qTM >= 0.5)
   - 75 KF positive controls (length and pLDDT matched)
2. Independent Domain Segmentation (outcome-blind):
   - Parse CA coordinates and B-factor pLDDT from PDB models.
   - Smooth pLDDT with 15-aa rolling window.
   - Partition into:
     * Confident domain segments: smoothed pLDDT >= 0.70, length >= 40 aa.
     * Low-confidence / disordered segments: pLDDT < 0.70 (tails, linkers).
   - Evaluate C-alpha contact compactness (distance < 10 A).
   - Map Pfam domain envelopes from pfam_domtbl.txt.
3. Residue-coordinate Foldseek search against PDB (pdb_20260730/foldseek_pdb):
   - Extract exact alignment intervals [qstart, qend].
4. Domain-level calibration and coverage mapping:
   - For each independently defined domain, measure exact overlap with top hits.
   - Classify domain status: matched (cov >= 0.70, alnTM >= 0.5) vs unmatched vs partial.
   - Quantify why low-qTM models failed whole-query TM:
     * Case A: Multi-domain with 1 matched domain and 1 unmatched confident domain.
     * Case B: Single folded domain + long disordered/low-confidence tail(s) (denominator inflation).
     * Case C: Flexible multi-domain linker reorientation (both domains match separately).
     * Case D: True whole-fold divergence (no confident domain matched).
5. Control Self-Exclusion calibration:
   - Recompute domain and chain match rates for controls with self-PDB chains excluded.
"""

import os
import sys
import csv
import json
import math
import random
import subprocess
from collections import defaultdict

W = "/ibex/scratch/projects/c2014/rob/metagenomics"
OUT_DIR = os.path.join(W, "06_darkmatter/p6_hypotheses_clean/h4_structure/results")
SCRIPTS_DIR = os.path.join(W, "06_darkmatter/p6_hypotheses_clean/h4_structure/scripts")

PER_REP = os.path.join(W, "04_structure/00_structure_clean/05_validate/per_rep_structures.tsv")
PER_CTRL = os.path.join(W, "04_structure/00_structure_clean/05_validate/per_control_structures.tsv")
MANIFEST = os.path.join(W, "rebuild_20260729/structural_calibration/current_reference_control_manifest_20260802.tsv")
QUERY_PDB_DIR = os.path.join(W, "04_structure/00_structure_clean/05_validate/query_pdb")
DOMTBL = os.path.join(W, "04_structure/pfam_candidates_clean/results/pfam_domtbl.txt")
PDB_DB = os.path.join(W, "rebuild_20260729/structural_calibration/references/pdb_20260730/foldseek_pdb")

PILOT_SEED = 20261013


def parse_pdb(pdb_path):
    """Extract CA atoms: residue index (1-based), resname, pLDDT, (x, y, z)."""
    ca_residues = []
    with open(pdb_path) as fh:
        for line in fh:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                try:
                    res_num = int(line[22:26].strip())
                    res_name = line[17:20].strip()
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    b = float(line[60:66].strip())
                    # normalize plddt if on 0-100 scale
                    plddt = b / 100.0 if b > 1.5 else b
                    ca_residues.append({
                        "res_num": res_num,
                        "res_name": res_name,
                        "plddt": plddt,
                        "coord": (x, y, z)
                    })
                except (ValueError, IndexError):
                    continue
    return ca_residues


def segment_domains(ca_residues, window=15, plddt_thresh=0.70, min_dom_len=40):
    """Independent domain segmentation:
    1. Compute smoothed rolling-window pLDDT.
    2. Identify continuous high-confidence blocks (pLDDT >= 0.70).
    3. Allow bridging of small dips (length <= 10 if mean pLDDT >= 0.60).
    4. Return segmented domain intervals [start_res, end_res] and disordered segments.
    """
    L = len(ca_residues)
    if L == 0:
        return [], []

    plddts = [r["plddt"] for r in ca_residues]
    smoothed = []
    hw = window // 2
    for i in range(L):
        s = max(0, i - hw)
        e = min(L, i + hw + 1)
        smoothed.append(sum(plddts[s:e]) / (e - s))

    # Binary mask of confident residues
    is_conf = [s >= plddt_thresh for s in smoothed]

    # Group into initial contiguous blocks
    raw_blocks = []
    in_block = False
    start = 0
    for i, conf in enumerate(is_conf):
        if conf and not in_block:
            in_block = True
            start = i
        elif not conf and in_block:
            in_block = False
            raw_blocks.append((start, i - 1))
    if in_block:
        raw_blocks.append((start, L - 1))

    # Merge blocks separated by short dips (<= 10 residues)
    merged_blocks = []
    for b in raw_blocks:
        if not merged_blocks:
            merged_blocks.append(list(b))
        else:
            prev_end = merged_blocks[-1][1]
            gap = b[0] - prev_end - 1
            if gap <= 10:
                # check mean pLDDT in gap
                gap_mean = sum(plddts[prev_end + 1: b[0]]) / max(1, gap)
                if gap_mean >= 0.55:
                    merged_blocks[-1][1] = b[1]
                else:
                    merged_blocks.append(list(b))
            else:
                merged_blocks.append(list(b))

    # Filter by minimum domain length
    domains = []
    disordered = []
    prev_idx = 0
    for s, e in merged_blocks:
        # any gap before this domain is disordered/linker
        if s > prev_idx:
            disordered.append((ca_residues[prev_idx]["res_num"], ca_residues[s - 1]["res_num"], s - prev_idx))
        length = e - s + 1
        if length >= min_dom_len:
            dom_plddt = sum(plddts[s:e + 1]) / length
            domains.append({
                "dom_id": len(domains) + 1,
                "start": ca_residues[s]["res_num"],
                "end": ca_residues[e]["res_num"],
                "length": length,
                "mean_plddt": dom_plddt,
            })
            prev_idx = e + 1
        else:
            # treat small block as disordered/structural peptide
            disordered.append((ca_residues[s]["res_num"], ca_residues[e]["res_num"], length))
            prev_idx = e + 1

    if prev_idx < L:
        disordered.append((ca_residues[prev_idx]["res_num"], ca_residues[L - 1]["res_num"], L - prev_idx))

    return domains, disordered


def load_pfam_domtbl(domtbl_path):
    """Parse pfam_domtbl: protein -> list of {pfam_acc, pfam_id, env_from, env_to, score}."""
    pfam_hits = defaultdict(list)
    if not os.path.exists(domtbl_path):
        return pfam_hits
    with open(domtbl_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split()
            if len(f) < 22:
                continue
            target, qacc, qname = f[0], f[4], f[3]
            dscore = float(f[13])
            env_from, env_to = int(f[19]), int(f[20])
            acc = qacc.split(".")[0] if qacc != "-" else qname
            pfam_hits[target].append({
                "pfam_acc": acc,
                "pfam_id": qname,
                "env_from": env_from,
                "env_to": env_to,
                "score": dscore
            })
    return pfam_hits


def select_pilot_cohort():
    """Deterministically select 250 models across 4 strata."""
    random.seed(PILOT_SEED)

    # Load candidates
    k_local = []
    k_dual_none = []
    q_not_k = []
    with open(PER_REP) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["confidence_qualified"] != "1":
                continue
            qid = r["id"]
            L = int(r["length"])
            pl = float(r["mean_plddt"])
            pq = float(r["best_pdb_qtm"])
            pa = float(r["best_pdb_alntm"])
            aq = float(r["best_afdb_qtm"])
            aa = float(r["best_afdb_alntm"])
            dual_q_novel = (pq < 0.5 and aq < 0.5)
            dual_a_novel = (pa < 0.5 and aa < 0.5)

            rec = {"id": qid, "length": L, "mean_plddt": pl, "best_pdb_qtm": pq,
                   "best_pdb_alntm": pa, "best_afdb_qtm": aq, "best_afdb_alntm": aa}
            if dual_q_novel and not dual_a_novel:
                k_local.append(rec)
            elif dual_q_novel and dual_a_novel:
                k_dual_none.append(rec)
            elif not dual_q_novel:
                q_not_k.append(rec)

    # Load controls
    ctrl_meta = {}
    with open(MANIFEST) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            ctrl_meta[r["blind_id"]] = r

    controls = []
    with open(PER_CTRL) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["confidence_qualified"] != "1":
                continue
            cid = r["id"]
            controls.append({
                "id": cid,
                "length": int(r["length"]),
                "mean_plddt": float(r["mean_plddt"]),
                "best_pdb_qtm": float(r["best_pdb_qtm"]),
                "best_pdb_alntm": float(r["best_pdb_alntm"]),
                "best_afdb_qtm": float(r["best_afdb_qtm"]),
                "best_afdb_alntm": float(r["best_afdb_alntm"]),
                "pdb_parent": ctrl_meta.get(cid, {}).get("pdb_structure_id", ""),
                "pdb4": ctrl_meta.get(cid, {}).get("pdb_structure_id", "").split("-")[0].lower(),
                "control_mode": ctrl_meta.get(cid, {}).get("control_mode", "")
            })

    # Sort deterministically
    k_local.sort(key=lambda x: x["id"])
    k_dual_none.sort(key=lambda x: x["id"])
    q_not_k.sort(key=lambda x: x["id"])
    controls.sort(key=lambda x: x["id"])

    # Sample K_local across length bins: [100, 200), [200, 300), [300, 400), [400, inf)
    bins = [(100, 200), (200, 300), (300, 400), (400, 10000)]
    k_local_strat = []
    for lo, hi in bins:
        sub = [r for r in k_local if lo <= r["length"] < hi]
        sampled = random.sample(sub, min(25, len(sub)))
        k_local_strat.extend(sampled)

    k_dual_none_sample = random.sample(k_dual_none, min(25, len(k_dual_none)))
    q_not_k_sample = random.sample(q_not_k, min(50, len(q_not_k)))
    controls_sample = random.sample(controls, min(75, len(controls)))

    pilot = {
        "K_local": k_local_strat,
        "K_dual_none": k_dual_none_sample,
        "Q_not_K": q_not_k_sample,
        "Controls": controls_sample,
    }
    return pilot


def run_foldseek_pilot(query_list, query_dir, out_tsv, threads=16):
    """Run Foldseek easy-search with full alignment coordinates."""
    # Write queries to temporary directory with symlinks
    tmp_qdir = os.path.join(OUT_DIR, "pilot_query_pdb")
    os.makedirs(tmp_qdir, exist_ok=True)
    for qid in query_list:
        src = os.path.join(QUERY_PDB_DIR, f"{qid}.pdb")
        dst = os.path.join(tmp_qdir, f"{qid}.pdb")
        if not os.path.exists(dst) and os.path.exists(src):
            os.symlink(src, dst)

    tmp_out = out_tsv + ".tmp"
    tmp_scratch = os.path.join(OUT_DIR, "fs_tmp")
    os.makedirs(tmp_scratch, exist_ok=True)

    fmt = "query,target,alntmscore,qtmscore,ttmscore,prob,lddt,evalue,bits,qcov,tcov,qlen,tlen,alnlen,qstart,qend,tstart,tend"
    cmd = [
        "foldseek", "easy-search",
        tmp_qdir, PDB_DB, tmp_out, tmp_scratch,
        "--alignment-type", "1", "-e", "10", "-s", "9.5", "--max-seqs", "50",
        "--threads", str(threads),
        "--format-output", fmt
    ]
    print(f"Running Foldseek command: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)
    os.rename(tmp_out, out_tsv)
    print(f"Foldseek search completed: {out_tsv}")


def main():
    print("=== Step 1: Deterministic Pilot Sampling ===")
    pilot = select_pilot_cohort()
    all_pilot_models = []
    pilot_meta = {}
    for stratum, models in pilot.items():
        print(f"Stratum {stratum}: {len(models)} models")
        for m in models:
            m["stratum"] = stratum
            all_pilot_models.append(m["id"])
            pilot_meta[m["id"]] = m

    # Step 2: Independent Domain Segmentation
    print("\n=== Step 2: Outcome-Blind Domain Segmentation ===")
    pfam_hits = load_pfam_domtbl(DOMTBL)
    segmentation_records = []
    domain_records = []

    for qid in all_pilot_models:
        pdb_path = os.path.join(QUERY_PDB_DIR, f"{qid}.pdb")
        if not os.path.exists(pdb_path):
            print(f"Warning: {pdb_path} not found")
            continue
        ca = parse_pdb(pdb_path)
        doms, disordered = segment_domains(ca)
        p_pfam = pfam_hits.get(qid, [])

        total_dom_len = sum(d["length"] for d in doms)
        total_dis_len = sum(d[2] for d in disordered)
        L = len(ca)

        segmentation_records.append({
            "id": qid,
            "stratum": pilot_meta[qid]["stratum"],
            "length": L,
            "mean_plddt": pilot_meta[qid]["mean_plddt"],
            "n_domains": len(doms),
            "domain_covered_len": total_dom_len,
            "domain_covered_frac": round(total_dom_len / max(1, L), 4),
            "disordered_len": total_dis_len,
            "disordered_frac": round(total_dis_len / max(1, L), 4),
            "n_pfam_domains": len(p_pfam),
        })

        for d in doms:
            # check Pfam overlap
            pfam_cov = 0
            for pf in p_pfam:
                ov = max(0, min(d["end"], pf["env_to"]) - max(d["start"], pf["env_from"]) + 1)
                pfam_cov = max(pfam_cov, ov)

            domain_records.append({
                "id": qid,
                "stratum": pilot_meta[qid]["stratum"],
                "dom_id": d["dom_id"],
                "start": d["start"],
                "end": d["end"],
                "length": d["length"],
                "mean_plddt": round(d["mean_plddt"], 4),
                "pfam_overlap_aa": pfam_cov,
                "pfam_overlap_frac": round(pfam_cov / d["length"], 4)
            })

    print(f"Segmented {len(all_pilot_models)} models into {len(domain_records)} structural domains.")

    # Step 3: Run Foldseek with Exact Coordinates
    print("\n=== Step 3: Foldseek Alignment with Exact Coordinates ===")
    fs_pilot_out = os.path.join(OUT_DIR, "pilot_pdb_hits_coords.tsv")
    run_foldseek_pilot(all_pilot_models, QUERY_PDB_DIR, fs_pilot_out)

    # Step 4: Map Alignment Intervals to Domain Boundaries
    print("\n=== Step 4: Exact Alignment Mapping to Domains ===")
    hits_by_query = defaultdict(list)
    with open(fs_pilot_out) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 18:
                continue
            q = p[0][:-4] if p[0].endswith(".pdb") else p[0]
            hits_by_query[q].append({
                "target": p[1],
                "alntm": float(p[2]),
                "qtm": float(p[3]),
                "ttm": float(p[4]),
                "prob": float(p[5]),
                "bits": float(p[8]),
                "qcov": float(p[9]),
                "qlen": int(p[11]),
                "alnlen": int(p[13]),
                "qstart": int(p[14]),
                "qend": int(p[15]),
            })

    # Evaluate domain coverage and match status
    # For each domain: compute max coverage by any single hit with alnTM >= 0.5
    calibrated_domains = []
    model_diagnostics = []

    # Map domains by query id
    doms_by_qid = defaultdict(list)
    for dr in domain_records:
        doms_by_qid[dr["id"]].append(dr)

    for qid in all_pilot_models:
        q_meta = pilot_meta[qid]
        q_doms = doms_by_qid.get(qid, [])
        q_hits = hits_by_query.get(qid, [])
        is_ctrl = q_meta["stratum"] == "Controls"
        pdb4 = q_meta.get("pdb4", "")

        # If control, separate self-hits from non-self hits
        if is_ctrl and pdb4:
            valid_hits = [h for h in q_hits if h["target"][:4].lower() != pdb4]
            self_hits = [h for h in q_hits if h["target"][:4].lower() == pdb4]
        else:
            valid_hits = q_hits
            self_hits = []

        # Global best hits
        best_hit_alntm = max([h["alntm"] for h in valid_hits], default=0.0)
        best_hit_qtm = max([h["qtm"] for h in valid_hits], default=0.0)
        self_hit_qtm = max([h["qtm"] for h in self_hits], default=0.0)

        # Domain level analysis
        n_matched_doms = 0
        n_unmatched_doms = 0
        dom_details = []

        for d in q_doms:
            d_len = d["length"]
            d_s = d["start"]
            d_e = d["end"]

            # Find best hit covering this domain
            best_cov = 0.0
            best_cov_alntm = 0.0
            best_cov_target = "NONE"

            for h in valid_hits:
                ov = max(0, min(d_e, h["qend"]) - max(d_s, h["qstart"]) + 1)
                cov = ov / d_len
                if cov > best_cov:
                    best_cov = cov
                    best_cov_alntm = h["alntm"]
                    best_cov_target = h["target"]

            # Domain is matched if covered >= 50% by an alignment with alnTM >= 0.5
            is_matched = int(best_cov >= 0.50 and best_cov_alntm >= 0.50)
            if is_matched:
                n_matched_doms += 1
            else:
                n_unmatched_doms += 1

            cal_d = dict(d)
            cal_d.update({
                "best_hit_target": best_cov_target,
                "best_hit_cov": round(best_cov, 4),
                "best_hit_alntm": round(best_cov_alntm, 4),
                "is_matched": is_matched,
            })
            calibrated_domains.append(cal_d)
            dom_details.append(cal_d)

        # Diagnostic classification for low-qTM candidates
        # Case A: Multi-domain with matched + unmatched confident domains
        # Case B: Single folded domain + disordered tail (inflated denominator)
        # Case C: Multi-domain with all domains matched (flexible linker reorientation)
        # Case D: True whole-protein novel/divergent fold (no domain matched)
        # Case E: Fold-matched (Q_not_K or Controls)
        diag_class = "UNKNOWN"
        if q_meta["stratum"] in ("K_local", "K_dual_none"):
            if len(q_doms) > 1 and n_matched_doms > 0 and n_unmatched_doms > 0:
                diag_class = "CASE_A_MULTI_DOMAIN_HYBRID"
            elif len(q_doms) == 1 and n_matched_doms == 1:
                diag_class = "CASE_B_DISORDERED_TAIL_INFLATION"
            elif len(q_doms) > 1 and n_matched_doms == len(q_doms):
                diag_class = "CASE_C_LINKER_REORIENTATION_ALL_MATCHED"
            elif n_matched_doms == 0:
                diag_class = "CASE_D_WHOLE_MODEL_UNMATCHED"
            else:
                diag_class = "CASE_OTHER"
        else:
            diag_class = f"REFERENCE_{q_meta['stratum']}"

        model_diagnostics.append({
            "id": qid,
            "stratum": q_meta["stratum"],
            "length": q_meta["length"],
            "mean_plddt": q_meta["mean_plddt"],
            "n_domains": len(q_doms),
            "n_matched_domains": n_matched_doms,
            "n_unmatched_domains": n_unmatched_doms,
            "diag_class": diag_class,
            "valid_best_pdb_alntm": round(best_hit_alntm, 4),
            "valid_best_pdb_qtm": round(best_hit_qtm, 4),
            "self_hit_pdb_qtm": round(self_hit_qtm, 4) if is_ctrl else 0.0,
            "has_retrieval": int(len(valid_hits) > 0),
        })

    # Summary Statistics by Stratum
    print("\n=== Summary of Model Diagnostics ===")
    diag_counts = defaultdict(lambda: defaultdict(int))
    for md in model_diagnostics:
        diag_counts[md["stratum"]][md["diag_class"]] += 1

    for stratum, counts in diag_counts.items():
        print(f"\nStratum: {stratum} (Total={sum(counts.values())})")
        for cls, n in sorted(counts.items()):
            print(f"  {cls}: {n} ({100.0 * n / sum(counts.values()):.1f}%)")

    # Save Results
    os.makedirs(OUT_DIR, exist_ok=True)
    with open(os.path.join(OUT_DIR, "pilot_model_diagnostics.tsv"), "w") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(model_diagnostics[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(model_diagnostics)

    with open(os.path.join(OUT_DIR, "pilot_domain_calibration.tsv"), "w") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(calibrated_domains[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(calibrated_domains)

    with open(os.path.join(OUT_DIR, "pilot_segmentation_summary.tsv"), "w") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(segmentation_records[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(segmentation_records)

    # Save summary JSON
    summary_out = {
        "pilot_size": len(all_pilot_models),
        "pilot_seed": PILOT_SEED,
        "stratum_counts": {s: len(m) for s, m in pilot.items()},
        "diagnostic_breakdown": {s: dict(c) for s, c in diag_counts.items()},
    }
    with open(os.path.join(OUT_DIR, "pilot_summary.json"), "w") as fh:
        json.dump(summary_out, fh, indent=2)

    print(f"\nPilot analysis complete! Outputs written to {OUT_DIR}")


if __name__ == "__main__":
    main()
