#!/usr/bin/env python3
"""v3_crop_and_search_regions.py — Outcome-Independent Confident Region Cropping,
Same-Run Whole-Chain & Cropped Foldseek Search (--max-seqs 1000), Symmetric
Multi-Tier Homology Exclusion, Horvitz-Thompson Design-Weighted Estimation,
and Exact CIGAR Aligned-Pair Calibration.

Repairs implemented for h4_review1 (C1-C5) and ROOT_AUDIT_NOTES.md:
- C1 (Control Hit-List Saturation vs Low Score): Executes Foldseek with --max-seqs 1000
  and records raw hit counts, saturation flags (total_hits_raw >= 1000), and per-tier
  remaining hit counts (t{k}_remaining_hits). Explicitly separates POST_EXCLUSION_EMPTY
  from PERSISTENT_LOW_SCORE / LOW_SCORE and RAW_ZERO_RETRIEVAL, and compares V3
  (--max-seqs 1000) directly against V2 (--max-seqs 50) on the identical 413 regions.
- C2 (No Eligible Region vs Zero Retrieval): Explicitly separates NO_CONFIDENT_REGION
  (n_confident_regions == 0, never searched as a crop) from RAW_ZERO_RETRIEVAL
  (n_confident_regions >= 1 and 0 raw hits returned), plus retains the full N=5,851
  finite-population region-eligibility census (v3_full5851_cropped_regions_manifest.tsv).
- C3 (Whole-Chain vs Cropped Symmetry & Reference-Consistent Strata): Executes both
  whole-chain and cropped-region searches in the same run against the identical
  PDB reference (pdb_20260730/foldseek_pdb) with identical parameters (--alignment-type 1
  -e 10 -s 9.5 --max-seqs 1000) and identical 21-column format including fident and cigar.
  Applies Tiers 0-4 symmetrically to both whole-chain and cropped searches, partitions
  reference-consistent strata (K_local_pdb, K_local_afdb_only, K_dual_none, Q_not_K_pdb,
  Q_not_K_afdb_only, Controls), and tracks both recovery (for whole-chain unmatched) and
  retention/loss (for whole-chain matched).
- C4 (Explicit Selection Weights & Design-Weighted Finite-Cohort Estimators): Records exact
  stratum/length-bin sampling probabilities pi_h = n_h / N_h and Horvitz-Thompson inverse-
  probability weights w_h = N_h / n_h (sum of weights = 5,851.0 across the N=300 stratified
  pilot), computes Horvitz-Thompson / Hajek design-weighted finite-population proportions
  and totals with finite-population-corrected (FPC) standard errors and 95% CIs across all
  strata and tiers, and validates HT estimates against the exact N=5,851 region census.
- C5 (Exact Aligned Residue Pairs from CIGAR & Statistical Reporting): Parses exact
  aligned residue pairs (M operations in cigar, asserted equal to qspan + tspan - alnlen
  on 100% of rows) distinct from query coordinate span (qend - qstart + 1), uses the
  parent protein model as the primary statistical unit across all tiers, repairs the
  3,013-control paired self-hit exclusion test (b=0, c=42, exact McNemar p=4.55e-13),
  and computes mean sequence identity strictly across retrieved hits.
"""

import os
import re
import sys
import csv
import json
import math
import time
import random
import shutil
import hashlib
import argparse
import subprocess
from collections import defaultdict, Counter

W = "/ibex/scratch/projects/c2014/rob/metagenomics"
H4_DIR = os.path.join(W, "06_darkmatter/p6_hypotheses_clean/h4_structure")
OUT_DIR = os.path.join(H4_DIR, "results")

PER_REP = os.path.join(W, "04_structure/00_structure_clean/05_validate/per_rep_structures.tsv")
PER_CTRL = os.path.join(W, "04_structure/00_structure_clean/05_validate/per_control_structures.tsv")
MANIFEST = os.path.join(W, "rebuild_20260729/structural_calibration/current_reference_control_manifest_20260802.tsv")
QUERY_PDB_DIR = os.path.join(W, "04_structure/00_structure_clean/05_validate/query_pdb")
PDB_DB = os.path.join(W, "rebuild_20260729/structural_calibration/references/pdb_20260730/foldseek_pdb")

PILOT_SEED = 20261013
MAX_SEQS = 1000
CIGAR_RE = re.compile(r"(\d+)([MID])")

SUB_STRATA_ORDER = [
    "K_local_pdb",
    "K_local_afdb_only",
    "K_dual_none",
    "Q_not_K_pdb",
    "Q_not_K_afdb_only",
    "Controls",
]

COMBINED_STRATA_MAP = {
    "K_local_all": ["K_local_pdb", "K_local_afdb_only"],
    "Q_not_K_all": ["Q_not_K_pdb", "Q_not_K_afdb_only"],
}


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1048576), b""):
            h.update(chunk)
    return h.hexdigest()


def wilson_ci(k, n, z=1.959964):
    if n == 0:
        return (0.0, 0.0)
    p = k / n
    denom = 1.0 + z * z / n
    center = (p + z * z / (2.0 * n)) / denom
    half = z * math.sqrt(p * (1.0 - p) / n + z * z / (4.0 * n * n)) / denom
    return (round(100.0 * max(0.0, center - half), 2), round(100.0 * min(1.0, center + half), 2))


def ht_domain_estimate(models_all, domain_pred, outcome_pred, z=1.959964):
    """Compute Horvitz-Thompson / Hajek domain ratio estimate, estimated population totals,
    and Taylor-linearized finite-population-corrected standard error across design strata.
    """
    by_h = defaultdict(list)
    for m in models_all:
        by_h[m["design_stratum_h"]].append(m)

    w_num = 0.0
    w_den = 0.0
    for m in models_all:
        w = m["pilot_ht_weight"]
        d = 1.0 if domain_pred(m) else 0.0
        y = 1.0 if (d > 0.0 and outcome_pred(m)) else 0.0
        w_den += w * d
        w_num += w * y

    if w_den <= 0.0:
        return {
            "est_pop_domain_N": 0.0,
            "est_pop_outcome_Y": 0.0,
            "ht_pct": 0.0,
            "fpc_se_pct": 0.0,
            "ci95": [0.0, 0.0],
        }

    p_hat = w_num / w_den
    var_total_linear = 0.0
    for h_key, m_list in by_h.items():
        n_h = len(m_list)
        if n_h <= 1:
            continue
        N_h = m_list[0]["pilot_pop_n"]
        f_h = n_h / N_h
        fpc = max(0.0, 1.0 - f_h)
        z_vals = []
        for m in m_list:
            d = 1.0 if domain_pred(m) else 0.0
            y = 1.0 if (d > 0.0 and outcome_pred(m)) else 0.0
            z_vals.append(d * (y - p_hat))
        z_bar = sum(z_vals) / n_h
        s2_z = sum((zv - z_bar) ** 2 for zv in z_vals) / (n_h - 1)
        var_total_linear += (N_h ** 2) * fpc * (s2_z / n_h)

    var_p = var_total_linear / (w_den ** 2)
    se_p = math.sqrt(max(0.0, var_p))
    lo = max(0.0, p_hat - z * se_p)
    hi = min(1.0, p_hat + z * se_p)
    return {
        "est_pop_domain_N": round(w_den, 2),
        "est_pop_outcome_Y": round(w_num, 2),
        "ht_pct": round(100.0 * p_hat, 2),
        "fpc_se_pct": round(100.0 * se_p, 2),
        "ci95": [round(100.0 * lo, 2), round(100.0 * hi, 2)],
    }


def norm_plddt(val):
    return val / 100.0 if val > 1.5 else val


def parse_pdb_ca(pdb_path):
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
                        "coord": (x, y, z),
                    })
                except (ValueError, IndexError):
                    continue
    return ca_atoms


def identify_confident_regions(ca_atoms, window=15, plddt_thresh=0.70, min_len=50, max_bridge=5):
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
                "mean_plddt": round(reg_mean_pl, 4),
            })
    return regions


def write_cropped_pdb(src_pdb, dst_pdb, start_res, end_res, source_id, reg_id):
    n_atoms = 0
    h = hashlib.sha256()
    with open(src_pdb) as fh_in, open(dst_pdb, "wb") as fh_out:
        h1 = f"REMARK 999 CROPPED CONFIDENT REGION\n".encode("ascii")
        h2 = f"REMARK 999 SOURCE {source_id} REGION {reg_id} RESIDUES {start_res} TO {end_res}\n".encode("ascii")
        fh_out.write(h1)
        fh_out.write(h2)
        h.update(h1)
        h.update(h2)
        for line in fh_in:
            if line.startswith("ATOM"):
                try:
                    rnum = int(line[22:26].strip())
                    if start_res <= rnum <= end_res:
                        bline = line.encode("ascii")
                        fh_out.write(bline)
                        h.update(bline)
                        n_atoms += 1
                except ValueError:
                    continue
        tail = b"TER\nEND\n"
        fh_out.write(tail)
        h.update(tail)
    return n_atoms, h.hexdigest()


def load_full_cohort():
    """Load all 2,838 confidence-qualified candidates and 3,013 controls,
    assigning reference-consistent sub-strata and reproducing V2 pilot flags/weights.
    """
    candidates = {}
    with open(PER_REP) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["confidence_qualified"] == "1":
                qid = r["id"]
                pq = float(r["best_pdb_qtm"])
                pa = float(r["best_pdb_alntm"])
                aq = float(r["best_afdb_qtm"])
                aa = float(r["best_afdb_alntm"])
                cdn = int(r["conf_qualified_dual_novel"])
                if cdn == 1:
                    if pa >= 0.5:
                        sub_st = "K_local_pdb"
                        par_st = "K_local_all"
                    elif aa >= 0.5:
                        sub_st = "K_local_afdb_only"
                        par_st = "K_local_all"
                    else:
                        sub_st = "K_dual_none"
                        par_st = "K_dual_none"
                else:
                    if pq >= 0.5:
                        sub_st = "Q_not_K_pdb"
                        par_st = "Q_not_K_all"
                    else:
                        sub_st = "Q_not_K_afdb_only"
                        par_st = "Q_not_K_all"
                candidates[qid] = {
                    "id": qid,
                    "is_control": 0,
                    "sub_stratum": sub_st,
                    "parent_stratum": par_st,
                    "length": int(r["length"]),
                    "mean_plddt": float(r["mean_plddt"]),
                    "hist_pdb_qtm": pq,
                    "hist_pdb_alntm": pa,
                    "hist_afdb_qtm": aq,
                    "hist_afdb_alntm": aa,
                    "conf_dual_novel": cdn,
                    "pdb_parent": "",
                    "pdb4": "",
                    "control_mode": "",
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
                    "is_control": 1,
                    "sub_stratum": "Controls",
                    "parent_stratum": "Controls",
                    "length": int(r["length"]),
                    "mean_plddt": float(r["mean_plddt"]),
                    "hist_pdb_qtm": float(r["best_pdb_qtm"]),
                    "hist_pdb_alntm": float(r["best_pdb_alntm"]),
                    "hist_afdb_qtm": float(r["best_afdb_qtm"]),
                    "hist_afdb_alntm": float(r["best_afdb_alntm"]),
                    "conf_dual_novel": int(r["conf_qualified_dual_novel"]),
                    "pdb_parent": pdb_id,
                    "pdb4": pdb4,
                    "control_mode": cm.get("control_mode", ""),
                }

    # Reproduce exact V2 pilot membership and Horvitz-Thompson weights for all 7 sampling cells
    random.seed(PILOT_SEED)
    k_local = sorted([c for c in candidates.values() if c["parent_stratum"] == "K_local_all"], key=lambda x: x["id"])
    k_dual_none = sorted([c for c in candidates.values() if c["sub_stratum"] == "K_dual_none"], key=lambda x: x["id"])
    q_not_k = sorted([c for c in candidates.values() if c["parent_stratum"] == "Q_not_K_all"], key=lambda x: x["id"])
    ctrl_list = sorted(list(controls.values()), key=lambda x: x["id"])

    bins = [(100, 200), (200, 300), (300, 400), (400, 10000)]
    pilot_weights = {}
    all_cell_meta = {}
    for lo, hi in bins:
        sub = [r for r in k_local if lo <= r["length"] < hi]
        sampled = random.sample(sub, min(25, len(sub)))
        w = len(sub) / len(sampled)
        h_name = f"K_local_len_{lo}_{hi}"
        for r in sub:
            all_cell_meta[r["id"]] = (h_name, f"[{lo},{hi})", len(sub), len(sampled), w)
        for s in sampled:
            pilot_weights[s["id"]] = (h_name, f"[{lo},{hi})", len(sub), len(sampled), w)

    for grp_name, grp, target_n in [
        ("K_dual_none_SRS", k_dual_none, 50),
        ("Q_not_K_SRS", q_not_k, 75),
        ("Controls_SRS", ctrl_list, 75),
    ]:
        sampled = random.sample(grp, min(target_n, len(grp)))
        w = len(grp) / len(sampled)
        for r in grp:
            all_cell_meta[r["id"]] = (grp_name, "SRS", len(grp), len(sampled), w)
        for s in sampled:
            pilot_weights[s["id"]] = (grp_name, "SRS", len(grp), len(sampled), w)

    cohort = {}
    for d in (candidates, controls):
        for qid, m in d.items():
            h_name, bname, npop, nsamp, w = all_cell_meta[qid]
            m["design_stratum_h"] = h_name
            m["pilot_bin"] = bname
            m["pilot_pop_n"] = npop
            m["pilot_samp_n"] = nsamp
            m["pilot_sampling_prob"] = round(nsamp / npop, 8)
            m["pilot_ht_weight"] = round(w, 6)
            m["in_v2_pilot"] = 1 if qid in pilot_weights else 0
            cohort[qid] = m

    return cohort


def run_foldseek_search(query_dir, db_path, out_tsv, tmp_dir, threads, max_seqs, reuse=False):
    if reuse and os.path.exists(out_tsv) and os.path.getsize(out_tsv) > 0:
        print(f"Reusing existing search output: {out_tsv}")
        return 0.0
    if os.path.exists(tmp_dir):
        shutil.rmtree(tmp_dir)
    os.makedirs(tmp_dir, exist_ok=True)
    tmp_out = out_tsv + ".tmp"
    if os.path.exists(tmp_out):
        os.remove(tmp_out)
    fmt = "query,target,alntmscore,qtmscore,ttmscore,prob,lddt,evalue,bits,fident,pident,qcov,tcov,qlen,tlen,alnlen,qstart,qend,tstart,tend,cigar"
    cmd = [
        "foldseek", "easy-search",
        query_dir, db_path, tmp_out, tmp_dir,
        "--alignment-type", "1", "-e", "10", "-s", "9.5", "--max-seqs", str(max_seqs),
        "--threads", str(threads),
        "--format-output", fmt,
    ]
    print(f"Executing: {' '.join(cmd)}")
    t0 = time.time()
    subprocess.run(cmd, check=True)
    os.rename(tmp_out, out_tsv)
    elapsed = round(time.time() - t0, 2)
    print(f"Completed {out_tsv} in {elapsed:.2f}s")
    if os.path.exists(tmp_dir):
        shutil.rmtree(tmp_dir)
    return elapsed


def parse_and_accumulate_hits(hit_tsv, query_to_pdb4):
    """Stream raw Foldseek hits (21 columns with cigar), validate exact CIGAR aligned pairs
    on 100% of rows, and accumulate Tier 0-4 summaries per query in O(N_queries) memory.
    """
    acc = {}
    for qid, pdb4 in query_to_pdb4.items():
        acc[qid] = {
            "pdb4": pdb4,
            "total_hits_raw": 0,
            "max_fident_raw": 0.0,
            "hits_fident_ge_20": 0,
            "hits_fident_ge_25": 0,
            "hits_fident_ge_30": 0,
            "rem_counts": [0, 0, 0, 0, 0],
            "top_qtm": [None, None, None, None, None],
            "top_alntm": [None, None, None, None, None],
        }

    total_rows = 0
    rows_with_q_ins = 0
    rows_with_t_ins = 0
    cand_rows_fident_ge20 = 0
    cand_rows_fident_ge25 = 0
    cand_rows_fident_ge30 = 0
    audited_examples = []

    with open(hit_tsv) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 21:
                continue
            total_rows += 1
            q = p[0][:-4] if p[0].endswith(".pdb") else p[0]
            if q not in acc:
                continue
            qa = acc[q]
            target = p[1]
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

            assert 0.0 <= fident <= 1.0, f"fident out of [0,1] scale: {fident}"

            qspan = qend - qstart + 1
            tspan = tend - tstart + 1
            cigar_m = 0
            cigar_tot = 0
            for length_str, op in CIGAR_RE.findall(cigar):
                lval = int(length_str)
                cigar_tot += lval
                if op == "M":
                    cigar_m += lval

            expected_m = qspan + tspan - alnlen
            assert cigar_tot == alnlen, f"CIGAR total {cigar_tot} != alnlen {alnlen} for {q} -> {target}"
            assert cigar_m == expected_m, f"CIGAR M {cigar_m} != qspan+tspan-alnlen {expected_m} for {q} -> {target}"

            if alnlen > tspan:
                rows_with_q_ins += 1
            if alnlen > qspan:
                rows_with_t_ins += 1

            qa["total_hits_raw"] += 1
            if fident > qa["max_fident_raw"]:
                qa["max_fident_raw"] = fident
            if fident >= 0.20:
                qa["hits_fident_ge_20"] += 1
                if not qa["pdb4"]:
                    cand_rows_fident_ge20 += 1
            if fident >= 0.25:
                qa["hits_fident_ge_25"] += 1
                if not qa["pdb4"]:
                    cand_rows_fident_ge25 += 1
            if fident >= 0.30:
                qa["hits_fident_ge_30"] += 1
                if not qa["pdb4"]:
                    cand_rows_fident_ge30 += 1

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

            # Tier 0: Unfiltered
            qa["rem_counts"][0] += 1
            if qa["top_qtm"][0] is None or qtm > qa["top_qtm"][0]["qtm"]:
                qa["top_qtm"][0] = row_dict
            if qa["top_alntm"][0] is None or alntm > qa["top_alntm"][0]["alntm"]:
                qa["top_alntm"][0] = row_dict

            # Tier 1: Source-PDB excluded
            pdb4 = qa["pdb4"]
            if pdb4 and target[:4].lower() == pdb4:
                continue

            qa["rem_counts"][1] += 1
            if qa["top_qtm"][1] is None or qtm > qa["top_qtm"][1]["qtm"]:
                qa["top_qtm"][1] = row_dict
            if qa["top_alntm"][1] is None or alntm > qa["top_alntm"][1]["alntm"]:
                qa["top_alntm"][1] = row_dict

            # Tiers 2, 3, 4: fident < 0.30, 0.25, 0.20
            for t_idx, thresh in ((2, 0.30), (3, 0.25), (4, 0.20)):
                if fident < thresh:
                    qa["rem_counts"][t_idx] += 1
                    if qa["top_qtm"][t_idx] is None or qtm > qa["top_qtm"][t_idx]["qtm"]:
                        qa["top_qtm"][t_idx] = row_dict
                    if qa["top_alntm"][t_idx] is None or alntm > qa["top_alntm"][t_idx]["alntm"]:
                        qa["top_alntm"][t_idx] = row_dict

            if len(audited_examples) < 60:
                if (not pdb4 and fident >= 0.20) or (qspan - cigar_m >= 18 and len(audited_examples) < 30):
                    audited_examples.append(row_dict)

    priority_queries = {
        "25Sr1::25Sr1_59902_2_reg1",
        "27Dr3::27Dr3_493525_1_reg1",
        "V34Dr2::V34Dr2_233315_1_reg1",
        "KF000066_reg1",
        "KF000126_reg1",
        "V34Dr2::V34Dr2_233315_1",
        "V12Sr1::V12Sr1_267603_7",
        "41PRr3::41PRr3_702372_18",
    }
    seen_keys = {(r["query"], r["target"], r["qstart"], r["tstart"]) for r in audited_examples}
    for qid in sorted(acc.keys()):
        if qid in priority_queries or len(audited_examples) < 110:
            for t_idx in (0, 1, 2, 3, 4):
                top_r = acc[qid]["top_qtm"][t_idx]
                if top_r is not None:
                    rk = (top_r["query"], top_r["target"], top_r["qstart"], top_r["tstart"])
                    if rk not in seen_keys:
                        audited_examples.append(top_r)
                        seen_keys.add(rk)

    stats = {
        "total_rows": total_rows,
        "rows_with_query_insertions_target_gaps": rows_with_q_ins,
        "pct_rows_with_query_insertions": round(100.0 * rows_with_q_ins / max(1, total_rows), 2),
        "rows_with_target_insertions_query_gaps": rows_with_t_ins,
        "pct_rows_with_target_insertions": round(100.0 * rows_with_t_ins / max(1, total_rows), 2),
        "candidate_rows_fident_ge_20": cand_rows_fident_ge20,
        "candidate_rows_fident_ge_25": cand_rows_fident_ge25,
        "candidate_rows_fident_ge_30": cand_rows_fident_ge30,
    }
    return acc, stats, audited_examples


def classify_single_unit(total_raw, rem_count, top_hit):
    if total_raw == 0:
        return "RAW_ZERO_RETRIEVAL"
    if rem_count == 0:
        return "POST_EXCLUSION_EMPTY"
    if top_hit is not None and top_hit["qtm"] >= 0.5:
        return "MATCHED"
    return "LOW_SCORE"


def classify_model_crops(n_regs, reg_evals, tier_k):
    if n_regs == 0:
        return "NO_CONFIDENT_REGION"
    if all(r["total_hits_raw"] == 0 for r in reg_evals):
        return "RAW_ZERO_RETRIEVAL"
    if all(r[f"t{tier_k}_remaining_hits"] == 0 for r in reg_evals):
        return "POST_EXCLUSION_EMPTY"
    if any(r[f"t{tier_k}_status"] == "MATCHED" for r in reg_evals):
        return "MATCHED"
    return "PERSISTENT_LOW_SCORE"


def classify_transition(wc_status, crop_status):
    if wc_status != "MATCHED":
        if crop_status == "MATCHED":
            return "WC_UNMATCHED_TO_CROP_RECOVERED"
        if crop_status == "PERSISTENT_LOW_SCORE":
            return "WC_UNMATCHED_TO_PERSISTENT_LOW_SCORE"
        if crop_status == "POST_EXCLUSION_EMPTY":
            return "WC_UNMATCHED_TO_POST_EXCLUSION_EMPTY"
        if crop_status == "RAW_ZERO_RETRIEVAL":
            return "WC_UNMATCHED_TO_RAW_ZERO_RETRIEVAL"
        if crop_status == "NO_CONFIDENT_REGION":
            return "WC_UNMATCHED_NO_CONFIDENT_REGION"
    else:
        if crop_status == "MATCHED":
            return "WC_MATCHED_RETAINED_IN_CROP"
        if crop_status == "PERSISTENT_LOW_SCORE":
            return "WC_MATCHED_LOST_LOW_SCORE"
        if crop_status == "POST_EXCLUSION_EMPTY":
            return "WC_MATCHED_LOST_POST_EXCLUSION_EMPTY"
        if crop_status == "RAW_ZERO_RETRIEVAL":
            return "WC_MATCHED_LOST_RAW_ZERO_RETRIEVAL"
        if crop_status == "NO_CONFIDENT_REGION":
            return "WC_MATCHED_LOST_NO_CONFIDENT_REGION"
    raise ValueError(f"Unhandled transition ({wc_status}, {crop_status})")


def summarize_stratum(models_all, models_sub, regs_sub, domain_pred):
    n_models = len(models_sub)
    n_regs = len(regs_sub)
    n_no_reg = sum(1 for m in models_sub if m["n_confident_regions"] == 0)
    n_elig = n_models - n_no_reg

    wc_raw_zero = sum(1 for m in models_sub if m["wc_total_hits_raw"] == 0)
    wc_raw_sat = sum(1 for m in models_sub if m["wc_raw_list_saturated"] == 1)
    reg_raw_zero = sum(1 for r in regs_sub if r["total_hits_raw"] == 0)
    reg_raw_sat = sum(1 for r in regs_sub if r["raw_list_saturated"] == 1)

    t1_ret_regs = [r for r in regs_sub if r["t1_remaining_hits"] > 0]
    mean_t1_fident_ret_only = (
        round(100.0 * sum(r["t1_fident"] for r in t1_ret_regs) / len(t1_ret_regs), 2)
        if t1_ret_regs else 0.0
    )
    max_t1_fident_ret_only = (
        round(100.0 * max(r["t1_fident"] for r in t1_ret_regs), 2)
        if t1_ret_regs else 0.0
    )
    regs_t1_fident_ge20 = sum(1 for r in t1_ret_regs if r["t1_fident"] >= 0.20)
    regs_t1_fident_ge25 = sum(1 for r in t1_ret_regs if r["t1_fident"] >= 0.25)
    regs_t1_fident_ge30 = sum(1 for r in t1_ret_regs if r["t1_fident"] >= 0.30)

    mean_t1_aligned_pair_qcov = (
        round(100.0 * sum(r["t1_aligned_pair_qcov"] for r in t1_ret_regs) / len(t1_ret_regs), 2)
        if t1_ret_regs else 0.0
    )
    mean_t1_qspan_cov = (
        round(100.0 * sum(r["t1_qspan_cov"] for r in t1_ret_regs) / len(t1_ret_regs), 2)
        if t1_ret_regs else 0.0
    )

    tiers_summary = {}
    for k in range(5):
        wc_stat_counts = Counter(m[f"wc_t{k}_status"] for m in models_sub)
        wc_unmatched_models = [m for m in models_sub if m[f"wc_t{k}_status"] != "MATCHED"]
        wc_matched_models = [m for m in models_sub if m[f"wc_t{k}_status"] == "MATCHED"]

        crop_stat_counts = Counter(m[f"crop_t{k}_status"] for m in models_sub)
        trans_counts = Counter(m[f"trans_t{k}"] for m in models_sub)

        reg_stat_counts = Counter(r[f"t{k}_status"] for r in regs_sub)
        n_reg_rem = sum(1 for r in regs_sub if r[f"t{k}_remaining_hits"] > 0)
        n_reg_match = reg_stat_counts["MATCHED"]

        n_mod_match = crop_stat_counts["MATCHED"]
        n_mod_rem = sum(1 for m in models_sub if m[f"crop_t{k}_status"] in ("MATCHED", "PERSISTENT_LOW_SCORE"))

        n_wc_unm = len(wc_unmatched_models)
        n_wc_unm_elig = sum(1 for m in wc_unmatched_models if m["n_confident_regions"] > 0)
        n_wc_unm_rec = trans_counts["WC_UNMATCHED_TO_CROP_RECOVERED"]

        n_wc_m = len(wc_matched_models)
        n_wc_m_ret = trans_counts["WC_MATCHED_RETAINED_IN_CROP"]

        # Horvitz-Thompson design-weighted estimates for this stratum/domain at Tier k
        ht_uncond_crop_match = ht_domain_estimate(
            models_all,
            domain_pred,
            lambda m, tk=k: m[f"crop_t{tk}_status"] == "MATCHED",
        )
        ht_elig_crop_match = ht_domain_estimate(
            models_all,
            lambda m: domain_pred(m) and m["n_confident_regions"] > 0,
            lambda m, tk=k: m[f"crop_t{tk}_status"] == "MATCHED",
        )
        ht_rem_crop_match = ht_domain_estimate(
            models_all,
            lambda m, tk=k: domain_pred(m) and m[f"crop_t{tk}_status"] in ("MATCHED", "PERSISTENT_LOW_SCORE"),
            lambda m, tk=k: m[f"crop_t{tk}_status"] == "MATCHED",
        )
        ht_wc_unm_rec_uncond = ht_domain_estimate(
            models_all,
            lambda m, tk=k: domain_pred(m) and m[f"wc_t{tk}_status"] != "MATCHED",
            lambda m, tk=k: m[f"trans_t{tk}"] == "WC_UNMATCHED_TO_CROP_RECOVERED",
        )
        ht_wc_unm_rec_elig = ht_domain_estimate(
            models_all,
            lambda m, tk=k: domain_pred(m) and m[f"wc_t{tk}_status"] != "MATCHED" and m["n_confident_regions"] > 0,
            lambda m, tk=k: m[f"trans_t{tk}"] == "WC_UNMATCHED_TO_CROP_RECOVERED",
        )
        ht_crop_pers_low = ht_domain_estimate(
            models_all,
            domain_pred,
            lambda m, tk=k: m[f"crop_t{tk}_status"] == "PERSISTENT_LOW_SCORE",
        )
        ht_crop_post_empty = ht_domain_estimate(
            models_all,
            domain_pred,
            lambda m, tk=k: m[f"crop_t{tk}_status"] == "POST_EXCLUSION_EMPTY",
        )
        ht_crop_raw_zero = ht_domain_estimate(
            models_all,
            domain_pred,
            lambda m, tk=k: m[f"crop_t{tk}_status"] == "RAW_ZERO_RETRIEVAL",
        )
        ht_crop_no_reg = ht_domain_estimate(
            models_all,
            domain_pred,
            lambda m, tk=k: m[f"crop_t{tk}_status"] == "NO_CONFIDENT_REGION",
        )

        tiers_summary[f"tier_{k}"] = {
            "wc_status_counts": {
                "MATCHED": wc_stat_counts["MATCHED"],
                "LOW_SCORE": wc_stat_counts["LOW_SCORE"],
                "POST_EXCLUSION_EMPTY": wc_stat_counts["POST_EXCLUSION_EMPTY"],
                "RAW_ZERO_RETRIEVAL": wc_stat_counts["RAW_ZERO_RETRIEVAL"],
            },
            "crop_model_status_counts": {
                "MATCHED": crop_stat_counts["MATCHED"],
                "PERSISTENT_LOW_SCORE": crop_stat_counts["PERSISTENT_LOW_SCORE"],
                "POST_EXCLUSION_EMPTY": crop_stat_counts["POST_EXCLUSION_EMPTY"],
                "RAW_ZERO_RETRIEVAL": crop_stat_counts["RAW_ZERO_RETRIEVAL"],
                "NO_CONFIDENT_REGION": crop_stat_counts["NO_CONFIDENT_REGION"],
            },
            "transition_counts": {
                "WC_UNMATCHED_TO_CROP_RECOVERED": trans_counts["WC_UNMATCHED_TO_CROP_RECOVERED"],
                "WC_UNMATCHED_TO_PERSISTENT_LOW_SCORE": trans_counts["WC_UNMATCHED_TO_PERSISTENT_LOW_SCORE"],
                "WC_UNMATCHED_TO_POST_EXCLUSION_EMPTY": trans_counts["WC_UNMATCHED_TO_POST_EXCLUSION_EMPTY"],
                "WC_UNMATCHED_TO_RAW_ZERO_RETRIEVAL": trans_counts["WC_UNMATCHED_TO_RAW_ZERO_RETRIEVAL"],
                "WC_UNMATCHED_NO_CONFIDENT_REGION": trans_counts["WC_UNMATCHED_NO_CONFIDENT_REGION"],
                "WC_MATCHED_RETAINED_IN_CROP": trans_counts["WC_MATCHED_RETAINED_IN_CROP"],
                "WC_MATCHED_LOST_LOW_SCORE": trans_counts["WC_MATCHED_LOST_LOW_SCORE"],
                "WC_MATCHED_LOST_POST_EXCLUSION_EMPTY": trans_counts["WC_MATCHED_LOST_POST_EXCLUSION_EMPTY"],
                "WC_MATCHED_LOST_RAW_ZERO_RETRIEVAL": trans_counts["WC_MATCHED_LOST_RAW_ZERO_RETRIEVAL"],
                "WC_MATCHED_LOST_NO_CONFIDENT_REGION": trans_counts["WC_MATCHED_LOST_NO_CONFIDENT_REGION"],
            },
            "model_rates": {
                "unconditional_crop_match_pct": round(100.0 * n_mod_match / max(1, n_models), 2),
                "unconditional_crop_match_ci95": list(wilson_ci(n_mod_match, n_models)),
                "region_eligible_crop_match_pct": round(100.0 * n_mod_match / max(1, n_elig), 2),
                "region_eligible_crop_match_ci95": list(wilson_ci(n_mod_match, n_elig)),
                "remaining_hit_conditional_crop_match_pct": round(100.0 * n_mod_match / max(1, n_mod_rem), 2),
                "remaining_hit_conditional_crop_match_ci95": list(wilson_ci(n_mod_match, n_mod_rem)),
                "wc_unmatched_count": n_wc_unm,
                "wc_unmatched_recovery_unconditional_pct": round(100.0 * n_wc_unm_rec / max(1, n_wc_unm), 2),
                "wc_unmatched_recovery_unconditional_ci95": list(wilson_ci(n_wc_unm_rec, n_wc_unm)),
                "wc_unmatched_eligible_count": n_wc_unm_elig,
                "wc_unmatched_recovery_eligible_pct": round(100.0 * n_wc_unm_rec / max(1, n_wc_unm_elig), 2),
                "wc_unmatched_recovery_eligible_ci95": list(wilson_ci(n_wc_unm_rec, n_wc_unm_elig)),
                "wc_matched_count": n_wc_m,
                "wc_matched_retained_pct": round(100.0 * n_wc_m_ret / max(1, n_wc_m), 2),
            },
            "horvitz_thompson_design_weighted": {
                "unconditional_crop_match": ht_uncond_crop_match,
                "region_eligible_crop_match": ht_elig_crop_match,
                "remaining_hit_conditional_crop_match": ht_rem_crop_match,
                "wc_unmatched_recovery_unconditional": ht_wc_unm_rec_uncond,
                "wc_unmatched_recovery_eligible": ht_wc_unm_rec_elig,
                "crop_persistent_low_score": ht_crop_pers_low,
                "crop_post_exclusion_empty": ht_crop_post_empty,
                "crop_raw_zero_retrieval": ht_crop_raw_zero,
                "crop_no_confident_region": ht_crop_no_reg,
            },
            "region_status_counts": {
                "MATCHED": reg_stat_counts["MATCHED"],
                "LOW_SCORE": reg_stat_counts["LOW_SCORE"],
                "POST_EXCLUSION_EMPTY": reg_stat_counts["POST_EXCLUSION_EMPTY"],
                "RAW_ZERO_RETRIEVAL": reg_stat_counts["RAW_ZERO_RETRIEVAL"],
            },
            "region_rates": {
                "unconditional_match_pct": round(100.0 * n_reg_match / max(1, n_regs), 2),
                "remaining_hit_conditional_match_pct": round(100.0 * n_reg_match / max(1, n_reg_rem), 2),
                "remaining_hit_regions": n_reg_rem,
            },
        }

    return {
        "model_count": n_models,
        "sum_horvitz_thompson_weights": round(sum(m["pilot_ht_weight"] for m in models_sub), 2),
        "models_with_confident_regions": n_elig,
        "models_no_confident_region": n_no_reg,
        "pct_no_confident_region": round(100.0 * n_no_reg / max(1, n_models), 2),
        "wc_raw_zero_retrieval_models": wc_raw_zero,
        "wc_raw_saturated_1000_models": wc_raw_sat,
        "total_confident_regions": n_regs,
        "mean_regions_per_model": round(n_regs / max(1, n_models), 3),
        "mean_regions_per_eligible_model": round(n_regs / max(1, n_elig), 3),
        "mean_region_length_aa": round(sum(r["reg_len"] for r in regs_sub) / max(1, n_regs), 2),
        "mean_region_plddt": round(sum(r["reg_mean_plddt"] for r in regs_sub) / max(1, n_regs), 4),
        "reg_raw_zero_retrieval_count": reg_raw_zero,
        "reg_raw_saturated_1000_count": reg_raw_sat,
        "t1_retrieved_regions_count": len(t1_ret_regs),
        "mean_t1_fident_retrieved_only_pct": mean_t1_fident_ret_only,
        "max_t1_fident_retrieved_only_pct": max_t1_fident_ret_only,
        "t1_regions_fident_ge_20": regs_t1_fident_ge20,
        "t1_regions_fident_ge_25": regs_t1_fident_ge25,
        "t1_regions_fident_ge_30": regs_t1_fident_ge30,
        "mean_t1_aligned_pair_qcov_retrieved_pct": mean_t1_aligned_pair_qcov,
        "mean_t1_qspan_cov_retrieved_pct": mean_t1_qspan_cov,
        "tiers": tiers_summary,
    }


def summarize_full_5851_region_manifest(full_manifest_path, full_cohort):
    """Summarize exact N=5,851 region-eligibility census from v3_full5851_cropped_regions_manifest.tsv."""
    if not os.path.exists(full_manifest_path):
        return None
    rows = list(csv.DictReader(open(full_manifest_path), delimiter="\t"))
    by_model = defaultdict(list)
    for r in rows:
        by_model[r["model_id"]].append(r)

    assert len(by_model) == 5851, f"Expected 5851 models in full manifest, got {len(by_model)}"
    census = {}
    for st in SUB_STRATA_ORDER + list(COMBINED_STRATA_MAP.keys()):
        if st in COMBINED_STRATA_MAP:
            sub_ids = [qid for qid, m in full_cohort.items() if m["sub_stratum"] in COMBINED_STRATA_MAP[st]]
        else:
            sub_ids = [qid for qid, m in full_cohort.items() if m["sub_stratum"] == st]
        n_m = len(sub_ids)
        n_no_reg = 0
        n_regs_tot = 0
        for qid in sub_ids:
            mrows = by_model[qid]
            n_r = int(mrows[0]["n_confident_regions"])
            if n_r == 0:
                n_no_reg += 1
            else:
                n_regs_tot += n_r
        census[st] = {
            "population_models": n_m,
            "models_with_confident_regions": n_m - n_no_reg,
            "models_no_confident_region": n_no_reg,
            "pct_no_confident_region": round(100.0 * n_no_reg / max(1, n_m), 2),
            "total_confident_regions": n_regs_tot,
            "mean_regions_per_model": round(n_regs_tot / max(1, n_m), 3),
        }
    return {
        "manifest_file": os.path.basename(full_manifest_path),
        "manifest_sha256": sha256_file(full_manifest_path),
        "total_manifest_rows": len(rows),
        "total_unique_models": len(by_model),
        "total_models_with_confident_regions": sum(1 for qid in by_model if int(by_model[qid][0]["n_confident_regions"]) > 0),
        "total_models_no_confident_region": sum(1 for qid in by_model if int(by_model[qid][0]["n_confident_regions"]) == 0),
        "total_cropped_regions": sum(1 for r in rows if int(r["n_confident_regions"]) > 0),
        "by_stratum": census,
    }


def compare_v2_vs_v3_regions(v2_reg_path, region_evaluations):
    """Compare V2 (--max-seqs 50) vs V3 (--max-seqs 1000) across the identical 413 regions."""
    if not os.path.exists(v2_reg_path):
        return None
    v2_rows = {f"{r['model_id']}_reg{r['region_id']}": r for r in csv.DictReader(open(v2_reg_path), delimiter="\t")}
    v3_rows = {r["crop_id"]: r for r in region_evaluations}
    if set(v2_rows.keys()) != set(v3_rows.keys()):
        return None

    ctrl_ids = [cid for cid, r in v3_rows.items() if r["sub_stratum"] == "Controls"]
    comp = {"n_control_regions": len(ctrl_ids), "tiers": {}}
    for k, label in ((1, "tier_1_source_excl"), (2, "tier_2_fident_lt_30"), (3, "tier_3_fident_lt_25"), (4, "tier_4_fident_lt_20")):
        v2_matched = sum(1 for cid in ctrl_ids if int(v2_rows[cid][f"t{k}_matched"]) == 1)
        v2_empty = sum(1 for cid in ctrl_ids if int(v2_rows[cid]["total_hits_raw"]) > 0 and float(v2_rows[cid][f"t{k}_qtm"]) == 0.0)
        v3_matched = sum(1 for cid in ctrl_ids if v3_rows[cid][f"t{k}_matched"] == 1)
        v3_empty = sum(1 for cid in ctrl_ids if v3_rows[cid][f"t{k}_status"] == "POST_EXCLUSION_EMPTY")
        v3_low = sum(1 for cid in ctrl_ids if v3_rows[cid][f"t{k}_status"] == "LOW_SCORE")
        v3_zero = sum(1 for cid in ctrl_ids if v3_rows[cid][f"t{k}_status"] == "RAW_ZERO_RETRIEVAL")
        rescued = sum(
            1 for cid in ctrl_ids
            if int(v2_rows[cid][f"t{k}_matched"]) == 0 and v3_rows[cid][f"t{k}_matched"] == 1
        )
        comp["tiers"][label] = {
            "v2_max_seqs_50_matched": v2_matched,
            "v2_max_seqs_50_match_pct": round(100.0 * v2_matched / len(ctrl_ids), 2),
            "v2_max_seqs_50_post_exclusion_empty": v2_empty,
            "v3_max_seqs_1000_matched": v3_matched,
            "v3_max_seqs_1000_match_pct": round(100.0 * v3_matched / len(ctrl_ids), 2),
            "v3_max_seqs_1000_post_exclusion_empty": v3_empty,
            "v3_max_seqs_1000_low_score": v3_low,
            "v3_max_seqs_1000_raw_zero": v3_zero,
            "v3_remaining_hit_conditional_match_pct": round(
                100.0 * v3_matched / max(1, v3_matched + v3_low), 2
            ),
            "regions_rescued_by_max_seqs_1000": rescued,
        }
    return comp


def compute_historical_3013_control_leakage_audit(leakage_path):
    """Exact paired McNemar / binomial audit on the 3,013 controls in control_leakage_audit.tsv."""
    if not os.path.exists(leakage_path):
        return None
    rows = list(csv.DictReader(open(leakage_path), delimiter="\t"))
    n_tot = len(rows)
    t0_unm = 0
    t1_unm = 0
    b_t0_unm_t1_m = 0
    c_t0_m_t1_unm = 0
    for r in rows:
        u_qtm = float(r["unfiltered_pdb_qtm"])
        ns_qtm = float(r["top_non_self_qtm"])
        t0_m = u_qtm >= 0.5
        t1_m = ns_qtm >= 0.5
        if not t0_m:
            t0_unm += 1
        if not t1_m:
            t1_unm += 1
        if (not t0_m) and t1_m:
            b_t0_unm_t1_m += 1
        elif t0_m and (not t1_m):
            c_t0_m_t1_unm += 1
    assert b_t0_unm_t1_m == 0
    exact_mcnemar_p = 2.0 * (0.5 ** c_t0_m_t1_unm)
    return {
        "n_controls": n_tot,
        "tier0_unfiltered_unmatched": t0_unm,
        "tier0_unfiltered_unmatched_pct": round(100.0 * t0_unm / n_tot, 2),
        "tier1_source_excl_unmatched": t1_unm,
        "tier1_source_excl_unmatched_pct": round(100.0 * t1_unm / n_tot, 2),
        "paired_discordant_b_fail_to_pass": b_t0_unm_t1_m,
        "paired_discordant_c_pass_to_fail": c_t0_m_t1_unm,
        "exact_paired_mcnemar_binomial_p": exact_mcnemar_p,
    }


def run_closure_assertions(mode, active_cohort, model_summaries, region_evaluations, stratum_summaries):
    """Run strict programmatic code-data closure assertions."""
    expected_n = 300 if mode == "pilot" else 5851
    assert len(active_cohort) == expected_n, f"Expected {expected_n} models, got {len(active_cohort)}"
    assert len(model_summaries) == expected_n, f"Expected {expected_n} model summaries, got {len(model_summaries)}"

    if mode == "pilot":
        tot_w = sum(m["pilot_ht_weight"] for m in model_summaries)
        assert abs(tot_w - 5851.0) < 0.01, f"Expected HT weight sum 5851.0, got {tot_w}"

    # 1. No-region exclusion and tier monotonicity per model
    for m in model_summaries:
        n_reg = m["n_confident_regions"]
        is_ctrl = m["is_control"] == 1
        for k in range(5):
            if n_reg == 0:
                assert m[f"crop_t{k}_status"] == "NO_CONFIDENT_REGION"
                assert m[f"trans_t{k}"] in ("WC_UNMATCHED_NO_CONFIDENT_REGION", "WC_MATCHED_LOST_NO_CONFIDENT_REGION")
            else:
                assert m[f"crop_t{k}_status"] != "NO_CONFIDENT_REGION"

        for k in range(4):
            assert m[f"wc_t{k}_remaining_hits"] >= m[f"wc_t{k+1}_remaining_hits"]
            assert m[f"wc_t{k}_qtm"] + 1e-9 >= m[f"wc_t{k+1}_qtm"]
            assert m[f"crop_t{k}_max_qtm"] + 1e-9 >= m[f"crop_t{k+1}_max_qtm"]

        if not is_ctrl:
            assert m["wc_t0_remaining_hits"] == m["wc_t1_remaining_hits"]
            assert m["wc_t0_qtm"] == m["wc_t1_qtm"]
            assert m["crop_t0_max_qtm"] == m["crop_t1_max_qtm"]
            assert m["trans_t0"] == m["trans_t1"]

    # 2. Region monotonicity and aligned-pair <= qspan <= reg_len
    for r in region_evaluations:
        for k in range(4):
            assert r[f"t{k}_remaining_hits"] >= r[f"t{k+1}_remaining_hits"]
            assert r[f"t{k}_qtm"] + 1e-9 >= r[f"t{k+1}_qtm"]
        for k in (0, 1):
            assert r[f"t{k}_aligned_pairs"] <= r[f"t{k}_qspan"] <= r["reg_len"]
            assert r[f"t{k}_aligned_pair_qcov"] <= r[f"t{k}_qspan_cov"] + 1e-6

    # 3. Stratum category total closure
    for st_name, st_data in stratum_summaries.items():
        n_m = st_data["model_count"]
        n_r = st_data["total_confident_regions"]
        for k in range(5):
            tk = st_data["tiers"][f"tier_{k}"]
            assert sum(tk["wc_status_counts"].values()) == n_m
            assert sum(tk["crop_model_status_counts"].values()) == n_m
            assert sum(tk["transition_counts"].values()) == n_m
            assert sum(tk["region_status_counts"].values()) == n_r
            assert tk["crop_model_status_counts"]["NO_CONFIDENT_REGION"] == st_data["models_no_confident_region"]

    print("ALL V3 CODE-DATA CLOSURE ASSERTIONS PASSED SUCCESSFULLY.")


def main():
    parser = argparse.ArgumentParser(description="H4 V3 Cropped & Whole-Chain Search Pipeline")
    parser.add_argument("--mode", choices=["pilot", "full"], default="pilot", help="Execution mode: pilot (N=300 stratified with HT weights + full 5,851 region census) or full (N=5,851)")
    parser.add_argument("--threads", type=int, default=16, help="Foldseek search threads")
    parser.add_argument("--max-seqs", type=int, default=MAX_SEQS, help="Foldseek --max-seqs prefilter cap")
    parser.add_argument("--reuse-search", action="store_true", help="Reuse existing v3 raw Foldseek outputs if present")
    args = parser.parse_args()

    print(f"=== Step 1: Loading Full Finite Cohort (N=5,851) and Selecting Mode={args.mode} ===")
    full_cohort = load_full_cohort()
    if args.mode == "pilot":
        active_cohort = {qid: m for qid, m in full_cohort.items() if m["in_v2_pilot"] == 1}
        wc_dir = os.path.join(H4_DIR, "v3_wc_pilot_queries")
        crops_dir = os.path.join(H4_DIR, "v3_pilot_crops")
    else:
        active_cohort = full_cohort
        wc_dir = os.path.join(H4_DIR, "v3_wc_queries")
        crops_dir = os.path.join(H4_DIR, "v3_crops")

    print(f"Active cohort size: {len(active_cohort)} models (out of {len(full_cohort)} finite population).")
    for st in SUB_STRATA_ORDER:
        cnt = sum(1 for m in active_cohort.values() if m["sub_stratum"] == st)
        pop_cnt = sum(1 for m in full_cohort.values() if m["sub_stratum"] == st)
        print(f"  {st:20s}: sample={cnt:4d} / pop={pop_cnt:4d}")

    print("\n=== Step 2: Preparing Whole-Chain Symlinks & Cropping Confident Regions (>=50 aa, pLDDT >= 0.70) ===")
    if os.path.exists(wc_dir):
        shutil.rmtree(wc_dir)
    if os.path.exists(crops_dir):
        shutil.rmtree(crops_dir)
    os.makedirs(wc_dir, exist_ok=True)
    os.makedirs(crops_dir, exist_ok=True)
    os.makedirs(OUT_DIR, exist_ok=True)

    model_region_manifest = []
    wc_query_to_pdb4 = {}
    crop_query_to_pdb4 = {}
    models_no_region = 0
    total_crops = 0

    for qid, m in sorted(active_cohort.items()):
        src_pdb = os.path.join(QUERY_PDB_DIR, f"{qid}.pdb")
        real_src = os.path.realpath(src_pdb)
        src_hash = sha256_file(real_src)
        m["source_pdb_sha256"] = src_hash

        wc_dst = os.path.join(wc_dir, f"{qid}.pdb")
        os.symlink(real_src, wc_dst)
        wc_query_to_pdb4[qid] = m["pdb4"]

        ca = parse_pdb_ca(real_src)
        assert len(ca) == m["length"], f"CA count {len(ca)} != manifest length {m['length']} for {qid}"

        regions = identify_confident_regions(ca, window=15, plddt_thresh=0.70, min_len=50, max_bridge=5)
        m["n_confident_regions"] = len(regions)
        m["total_confident_aa"] = sum(r["length"] for r in regions)
        m["confident_aa_fraction"] = round(m["total_confident_aa"] / m["length"], 4)

        if not regions:
            models_no_region += 1
            model_region_manifest.append({
                "model_id": qid,
                "sub_stratum": m["sub_stratum"],
                "parent_stratum": m["parent_stratum"],
                "design_stratum_h": m["design_stratum_h"],
                "pilot_ht_weight": m["pilot_ht_weight"],
                "whole_chain_length": m["length"],
                "whole_chain_mean_plddt": m["mean_plddt"],
                "n_confident_regions": 0,
                "region_id": 0,
                "start_idx_0based": -1,
                "end_idx_0based": -1,
                "reg_start": 0,
                "reg_end": 0,
                "reg_len": 0,
                "reg_mean_plddt": 0.0,
                "cropped_pdb": "NONE",
                "source_pdb_sha256": src_hash,
                "crop_pdb_sha256": "NONE",
            })
            continue

        for reg in regions:
            total_crops += 1
            crop_id = f"{qid}_reg{reg['reg_id']}"
            crop_filename = f"{crop_id}.pdb"
            crop_path = os.path.join(crops_dir, crop_filename)
            _, crop_hash = write_cropped_pdb(real_src, crop_path, reg["start_res"], reg["end_res"], qid, reg["reg_id"])
            crop_query_to_pdb4[crop_id] = m["pdb4"]
            model_region_manifest.append({
                "model_id": qid,
                "sub_stratum": m["sub_stratum"],
                "parent_stratum": m["parent_stratum"],
                "design_stratum_h": m["design_stratum_h"],
                "pilot_ht_weight": m["pilot_ht_weight"],
                "whole_chain_length": m["length"],
                "whole_chain_mean_plddt": m["mean_plddt"],
                "n_confident_regions": len(regions),
                "region_id": reg["reg_id"],
                "start_idx_0based": reg["start_idx"],
                "end_idx_0based": reg["end_idx"],
                "reg_start": reg["start_res"],
                "reg_end": reg["end_res"],
                "reg_len": reg["length"],
                "reg_mean_plddt": reg["mean_plddt"],
                "cropped_pdb": crop_filename,
                "source_pdb_sha256": src_hash,
                "crop_pdb_sha256": crop_hash,
            })

    print(f"Created {total_crops} cropped region PDBs from {len(active_cohort) - models_no_region} models.")
    print(f"Models with 0 confident regions (>=50 aa, smoothed pLDDT >= 0.70): {models_no_region}")

    manifest_out = os.path.join(OUT_DIR, "v3_cropped_regions_manifest.tsv")
    with open(manifest_out, "w") as fh:
        w = csv.DictWriter(fh, fieldnames=list(model_region_manifest[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(model_region_manifest)

    print(f"\n=== Step 3: Same-Run Foldseek Searches (Whole-Chain & Cropped Regions, --max-seqs {args.max_seqs}) ===")
    wc_hits_path = os.path.join(OUT_DIR, "v3_wc_pdb_hits.raw.tsv")
    crop_hits_path = os.path.join(OUT_DIR, "v3_crop_pdb_hits.raw.tsv")
    wc_tmp = os.path.join(H4_DIR, "v3_fs_wc_tmp")
    crop_tmp = os.path.join(H4_DIR, "v3_fs_crop_tmp")

    wc_elapsed = run_foldseek_search(wc_dir, PDB_DB, wc_hits_path, wc_tmp, args.threads, args.max_seqs, reuse=args.reuse_search)
    crop_elapsed = run_foldseek_search(crops_dir, PDB_DB, crop_hits_path, crop_tmp, args.threads, args.max_seqs, reuse=args.reuse_search)

    print("\n=== Step 4: Streaming & Validating Whole-Chain and Cropped Hits (CIGAR & Tiers 0-4) ===")
    wc_acc, wc_hit_stats, wc_examples = parse_and_accumulate_hits(wc_hits_path, wc_query_to_pdb4)
    crop_acc, crop_hit_stats, crop_examples = parse_and_accumulate_hits(crop_hits_path, crop_query_to_pdb4)
    print(f"Validated {wc_hit_stats['total_rows']} whole-chain hit rows and {crop_hit_stats['total_rows']} cropped hit rows.")

    region_evaluations = []
    regs_by_model = defaultdict(list)

    for rm in model_region_manifest:
        if rm["cropped_pdb"] == "NONE":
            continue
        qid = rm["model_id"]
        reg_id = rm["region_id"]
        crop_id = f"{qid}_reg{reg_id}"
        ca_info = crop_acc[crop_id]
        raw_tot = ca_info["total_hits_raw"]
        raw_sat = int(raw_tot >= args.max_seqs)

        t0_top = ca_info["top_qtm"][0]
        t1_top = ca_info["top_qtm"][1]

        reg_row = {
            "model_id": qid,
            "region_id": reg_id,
            "crop_id": crop_id,
            "sub_stratum": rm["sub_stratum"],
            "parent_stratum": rm["parent_stratum"],
            "design_stratum_h": rm["design_stratum_h"],
            "pilot_ht_weight": rm["pilot_ht_weight"],
            "reg_start": rm["reg_start"],
            "reg_end": rm["reg_end"],
            "reg_len": rm["reg_len"],
            "reg_mean_plddt": rm["reg_mean_plddt"],
            "total_hits_raw": raw_tot,
            "raw_list_saturated": raw_sat,
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
    for qid, m in sorted(active_cohort.items()):
        wa = wc_acc[qid]
        wc_raw = wa["total_hits_raw"]
        wc_sat = int(wc_raw >= args.max_seqs)
        reg_evals = regs_by_model.get(qid, [])
        n_regs = len(reg_evals)

        mrow = {
            "model_id": qid,
            "is_control": m["is_control"],
            "sub_stratum": m["sub_stratum"],
            "parent_stratum": m["parent_stratum"],
            "design_stratum_h": m["design_stratum_h"],
            "in_v2_pilot": m["in_v2_pilot"],
            "pilot_bin": m["pilot_bin"],
            "pilot_pop_n": m["pilot_pop_n"],
            "pilot_samp_n": m["pilot_samp_n"],
            "pilot_sampling_prob": m["pilot_sampling_prob"],
            "pilot_ht_weight": m["pilot_ht_weight"],
            "whole_chain_length": m["length"],
            "whole_chain_plddt": m["mean_plddt"],
            "hist_pdb_qtm": m["hist_pdb_qtm"],
            "hist_pdb_alntm": m["hist_pdb_alntm"],
            "hist_afdb_qtm": m["hist_afdb_qtm"],
            "hist_afdb_alntm": m["hist_afdb_alntm"],
            "n_confident_regions": n_regs,
            "total_confident_aa": m["total_confident_aa"],
            "confident_aa_fraction": m["confident_aa_fraction"],
            "wc_total_hits_raw": wc_raw,
            "wc_raw_list_saturated": wc_sat,
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
        stratum_summaries[st] = summarize_stratum(
            model_summaries, m_sub, r_sub, lambda m, s=st: m["sub_stratum"] == s
        )

    for comb_name, sub_list in COMBINED_STRATA_MAP.items():
        m_sub = [m for m in model_summaries if m["sub_stratum"] in sub_list]
        r_sub = [r for r in region_evaluations if r["sub_stratum"] in sub_list]
        stratum_summaries[comb_name] = summarize_stratum(
            model_summaries, m_sub, r_sub, lambda m, sl=sub_list: m["sub_stratum"] in sl
        )

    run_closure_assertions(args.mode, active_cohort, model_summaries, region_evaluations, stratum_summaries)

    # Full 5,851 region manifest census summary
    full_5851_manifest_path = os.path.join(OUT_DIR, "v3_full5851_cropped_regions_manifest.tsv")
    full_5851_census = summarize_full_5851_region_manifest(full_5851_manifest_path, full_cohort)

    # V2 (--max-seqs 50) vs V3 (--max-seqs 1000) paired comparison on the 413 regions
    v2_reg_path = os.path.join(OUT_DIR, "v2_region_evaluations.tsv")
    v2_vs_v3_comp = compare_v2_vs_v3_regions(v2_reg_path, region_evaluations)

    # Historical 3,013-control paired self-hit exclusion audit
    hist_leakage_audit = compute_historical_3013_control_leakage_audit(
        os.path.join(OUT_DIR, "control_leakage_audit.tsv")
    )

    # Whole-chain control self-hit paired exclusion audit at --max-seqs 1000
    ctrl_models = [m for m in model_summaries if m["sub_stratum"] == "Controls"]
    b_t0_fail_t1_pass = sum(1 for m in ctrl_models if m["wc_t0_status"] != "MATCHED" and m["wc_t1_status"] == "MATCHED")
    c_t0_pass_t1_fail = sum(1 for m in ctrl_models if m["wc_t0_status"] == "MATCHED" and m["wc_t1_status"] != "MATCHED")
    assert b_t0_fail_t1_pass == 0, "Removing self-hits cannot increase qTM"

    # Write output files
    with open(os.path.join(OUT_DIR, "v3_region_evaluations.tsv"), "w") as fh:
        w = csv.DictWriter(fh, fieldnames=list(region_evaluations[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(region_evaluations)

    with open(os.path.join(OUT_DIR, "v3_model_summaries.tsv"), "w") as fh:
        w = csv.DictWriter(fh, fieldnames=list(model_summaries[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(model_summaries)

    combined_examples = wc_examples + crop_examples
    if combined_examples:
        with open(os.path.join(OUT_DIR, "v3_audited_example_hits.tsv"), "w") as fh:
            w = csv.DictWriter(fh, fieldnames=list(combined_examples[0].keys()), delimiter="\t")
            w.writeheader()
            w.writerows(combined_examples)

    wc_raw_hash = sha256_file(wc_hits_path)
    crop_raw_hash = sha256_file(crop_hits_path)

    summary_json = {
        "analysis_version": "v3_symmetric_whole_chain_and_cropped_search",
        "cohort_mode": args.mode,
        "finite_population_models": len(full_cohort),
        "evaluated_models": len(model_summaries),
        "evaluated_cropped_regions": len(region_evaluations),
        "sum_horvitz_thompson_weights": round(sum(m["pilot_ht_weight"] for m in model_summaries), 2),
        "foldseek_settings": {
            "version": "10-941cd33 (941cd33ff0771cd2e3f144e3293e22a2b87e9fda)",
            "database": PDB_DB,
            "alignment_type": 1,
            "evalue": 10.0,
            "sensitivity": 9.5,
            "max_seqs": args.max_seqs,
            "cigar_validated_100pct": True,
            "wc_search_walltime_sec": wc_elapsed,
            "crop_search_walltime_sec": crop_elapsed,
        },
        "full_cohort_scaling_benchmark_job_53704591": {
            "slurm_job_id": 53704591,
            "node": "cn506-11-l",
            "cpus": 16,
            "full_cohort_cropping_completed": True,
            "full_cohort_models_cropped": 5851,
            "full_cohort_regions_created": 7773,
            "full_cohort_no_region_models": 161,
            "wc_prefilter_walltime_sec": 146.5,
            "wc_prefilter_passed_seqs_per_query": 990,
            "wc_structurealign_walltime_sec": 32.9,
            "wc_tmalign_observed_progress": "9/65 (13.85%) in 702s (11m 42s)",
            "projected_full_wc_tmalign_min": 84.5,
            "projected_full_wc_plus_crop_total_min": 183.0,
            "exceeds_1h_16cpu_allocation_limit": True,
        },
        "raw_hit_files_provenance": {
            "v3_wc_pdb_hits.raw.tsv": {
                "size_bytes": os.path.getsize(wc_hits_path),
                "sha256": wc_raw_hash,
                **wc_hit_stats,
            },
            "v3_crop_pdb_hits.raw.tsv": {
                "size_bytes": os.path.getsize(crop_hits_path),
                "sha256": crop_raw_hash,
                **crop_hit_stats,
            },
        },
        "historical_3013_control_paired_leakage_audit": hist_leakage_audit,
        "v3_control_wc_paired_self_exclusion": {
            "n_controls": len(ctrl_models),
            "t0_unfiltered_unmatched": sum(1 for m in ctrl_models if m["wc_t0_status"] != "MATCHED"),
            "t1_source_excl_unmatched": sum(1 for m in ctrl_models if m["wc_t1_status"] != "MATCHED"),
            "paired_discordant_b_fail_to_pass": b_t0_fail_t1_pass,
            "paired_discordant_c_pass_to_fail": c_t0_pass_t1_fail,
        },
        "v2_max_seqs_50_vs_v3_max_seqs_1000_control_regions": v2_vs_v3_comp,
        "full_5851_region_eligibility_census": full_5851_census,
        "stratum_summaries": stratum_summaries,
    }

    with open(os.path.join(OUT_DIR, "v3_recovery_summary.json"), "w") as fh:
        json.dump(summary_json, fh, indent=2)

    print("\n=== SUMMARY OF V3 RESULTS (Tier 1: Source-PDB Excluded) ===")
    for st in SUB_STRATA_ORDER + list(COMBINED_STRATA_MAP.keys()):
        sd = stratum_summaries[st]
        t1 = sd["tiers"]["tier_1"]
        mr = t1["model_rates"]
        ht = t1["horvitz_thompson_design_weighted"]
        print(
            f"{st:20s} | n={sd['model_count']:3d} (Elig={sd['models_with_confident_regions']:3d}, "
            f"0Reg={sd['models_no_confident_region']:2d}) | "
            f"Unw_CropMatch={mr['unconditional_crop_match_pct']:6.2f}% | "
            f"HT_CropMatch={ht['unconditional_crop_match']['ht_pct']:6.2f}% "
            f"(95%CI {ht['unconditional_crop_match']['ci95']}) | "
            f"WC_Unm={mr['wc_unmatched_count']:3d} -> Rec={t1['transition_counts']['WC_UNMATCHED_TO_CROP_RECOVERED']:3d} "
            f"(Unw={mr['wc_unmatched_recovery_unconditional_pct']:6.2f}%, HT={ht['wc_unmatched_recovery_unconditional']['ht_pct']:6.2f}%)"
        )

    print(f"\nAll V3 outputs successfully written and verified in {OUT_DIR}")


if __name__ == "__main__":
    main()
