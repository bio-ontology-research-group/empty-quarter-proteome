#!/usr/bin/env python3
"""v7_crop_and_search_regions.py — Full-Cohort Confident Region Cropping, Same-Run Whole-Chain
& Cropped Foldseek Searches (--max-seqs 1000), Symmetric Multi-Tier Homology Exclusion,
Hit-Cap Saturation Audit with Adaptive Rerun (--max-seqs 5000), and Exact Finite-Population
Census for all N=5,851 Models (7,773 Regions).

Repaired Version (Round 7 Root-Audited):
- Fixes NameError in audit_hit_cap_saturation (uses acc[qid]).
- Evaluates per-query hit cap and flags capped queries lacking MATCHED status at any tier
  (covering both POST_EXCLUSION_EMPTY and LOW_SCORE).
- Separates whole-chain and crop rerun queries, raw hit files, and temporary directories.
- Full-row exact duplicate removal for reruns, retaining score-distinct alignments deterministically;
  preserves best scores independently at every exclusion tier (score monotonicity).
- Evaluates saturation using standalone rerun count (not union count), with post-rerun final audit.
- Returns uniform keys across rerun and no-rerun execution paths.
- Enforces atomic completion signatures on Foldseek searches (query hashes, DB hashes, settings,
  raw SHA-256, subprocess returncode 0); refuses unsigned or partial results.
- Validates baseline reference identity and input hashes.
- Exact finite-population census without misleading sampling-error CIs in primary tables.
  Undefined domains return null/None (no [0.0, 0.0] intervals).
- Optional superpopulation model inferential stats explicitly labeled with independence assumptions.
- Generates clean SHA-256 manifest excluding itself.
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
from datetime import datetime, timezone
from collections import defaultdict, Counter

W = "/ibex/scratch/projects/c2014/rob/metagenomics"
H4_DIR = os.path.join(W, "06_darkmatter/p7_analyses_clean/h4_structure")
OUT_DIR = os.path.join(H4_DIR, "results")

PER_REP = os.path.join(W, "04_structure/00_structure_clean/05_validate/per_rep_structures.tsv")
PER_CTRL = os.path.join(W, "04_structure/00_structure_clean/05_validate/per_control_structures.tsv")
MANIFEST = os.path.join(W, "rebuild_20260729/structural_calibration/current_reference_control_manifest_20260802.tsv")
QUERY_PDB_DIR = os.path.join(W, "04_structure/00_structure_clean/05_validate/query_pdb")
PDB_DB = os.path.join(W, "rebuild_20260729/structural_calibration/references/pdb_20260730/foldseek_pdb")

PILOT_SEED = 20261013
MAX_SEQS = 1000
RERUN_MAX_SEQS = 5000
CIGAR_RE = re.compile(r"(\d+)([MID])")

BASELINE_REFERENCES = {
    "per_rep_structures.tsv": {
        "path": PER_REP,
        "size_bytes": 1436582,
        "sha256": "9d2aa6f85f64ba1458a86ab60d144a58d35ad55163a1a56d6d177668906e8c5b",
    },
    "per_control_structures.tsv": {
        "path": PER_CTRL,
        "size_bytes": 250513,
        "sha256": "33869ed27706b643b10874b1f54584b08f69a91c76744dfb1549ddb923a42996",
    },
    "current_reference_control_manifest_20260802.tsv": {
        "path": MANIFEST,
        "size_bytes": 1502037,
        "sha256": "a65ba1b05cf3c398b29a976a96c39f89c0544f2fa354a90a792fb6d65f16e8b6",
    },
}

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


def validate_baseline_references():
    """Validate referenced baseline tables against pinned size and SHA-256 hashes."""
    results = {}
    for name, ref in BASELINE_REFERENCES.items():
        p = ref["path"]
        if os.path.exists(p):
            sz = os.path.getsize(p)
            digest = sha256_file(p)
            assert sz == ref["size_bytes"], f"Size mismatch for {name}: {sz} != {ref['size_bytes']}"
            assert digest == ref["sha256"], f"SHA-256 mismatch for {name}: {digest} != {ref['sha256']}"
            results[name] = {"status": "VERIFIED", "size_bytes": sz, "sha256": digest}
        else:
            results[name] = {"status": "NOT_ACCESSIBLE_LOCALLY", "path": p}
    return results


def safe_rate(k, n):
    """Compute rate as percentage, returning None if n <= 0 (undefined domain)."""
    if n <= 0:
        return None
    return round(100.0 * k / n, 2)


def wilson_ci_model_based(k, n, z=1.959964):
    """Wilson score interval under superpopulation Bernoulli independence assumptions.
    Returns None if n <= 0. For k==0, returns [0.0, hi] with hi > 0.
    For k==n, returns [lo, 100.0] with lo < 100.0.
    """
    if n <= 0:
        return None
    p = k / n
    denom = 1.0 + z * z / n
    center = (p + z * z / (2.0 * n)) / denom
    half = z * math.sqrt(p * (1.0 - p) / n + z * z / (4.0 * n * n)) / denom
    lo = round(100.0 * max(0.0, center - half), 2)
    hi = round(100.0 * min(1.0, center + half), 2)
    return [lo, hi]


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
    assigning reference-consistent sub-strata and retaining round6 pilot flags for comparison.
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


def run_foldseek_search_atomic(query_dir, db_path, out_tsv, tmp_dir, threads, max_seqs, reuse=False):
    """Run Foldseek search with atomic signature generation and validation.
    When reuse=True: strictly validates query hashes, settings, and signed raw SHA-256.
    Refuses unsigned or partial files.
    """
    sig_path = out_tsv + ".signature.json"

    # Compute query hashes for requested queries
    query_files = sorted([f for f in os.listdir(query_dir) if f.endswith(".pdb")])
    assert len(query_files) > 0, f"No query PDB files in {query_dir}"
    query_hashes = {}
    for qf in query_files:
        qid = qf[:-4]
        query_hashes[qid] = sha256_file(os.path.join(query_dir, qf))

    # Check reuse condition with atomic signature
    if reuse:
        if os.path.exists(sig_path) and os.path.exists(out_tsv) and os.path.getsize(out_tsv) > 0:
            try:
                with open(sig_path) as fh:
                    sig = json.load(fh)
                assert sig.get("status") == "COMPLETED", "Signature status not COMPLETED"
                assert sig.get("returncode") == 0, "Signature returncode not 0"
                assert sig.get("settings", {}).get("max_seqs") == max_seqs, f"Signature max_seqs {sig.get('settings', {}).get('max_seqs')} != {max_seqs}"
                actual_sha = sha256_file(out_tsv)
                assert actual_sha == sig.get("raw_sha256"), f"Raw SHA mismatch: {actual_sha} != {sig.get('raw_sha256')}"
                assert sig.get("query_count") == len(query_files), f"Query count mismatch: {sig.get('query_count')} != {len(query_files)}"
                # Validate exact query IDs
                for qid, q_hash in query_hashes.items():
                    assert sig.get("query_id_hashes", {}).get(qid) == q_hash, f"Query hash mismatch for {qid}"
                print(f"Verified atomic completion signature for existing search output: {out_tsv}")
                return sig.get("elapsed_sec", 0.0)
            except Exception as e:
                print(f"Reuse validation failed ({e}); re-running search.")
        else:
            print(f"No valid atomic signature found for {out_tsv}; executing search fresh.")

    # Fresh execution
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
    raw_hash = sha256_file(out_tsv)

    # Get DB component metadata if available
    db_ver_path = db_path + ".version"
    db_ver_hash = sha256_file(db_ver_path) if os.path.exists(db_ver_path) else "n/a"

    signature = {
        "status": "COMPLETED",
        "out_tsv": out_tsv,
        "raw_sha256": raw_hash,
        "raw_size_bytes": os.path.getsize(out_tsv),
        "query_count": len(query_files),
        "query_id_hashes": query_hashes,
        "db_path": db_path,
        "db_version_sha256": db_ver_hash,
        "command": cmd,
        "settings": {
            "alignment_type": 1,
            "evalue": 10.0,
            "sensitivity": 9.5,
            "max_seqs": max_seqs,
            "threads": threads,
        },
        "returncode": 0,
        "elapsed_sec": elapsed,
        "completed_at_utc": datetime.now(timezone.utc).isoformat(),
    }
    tmp_sig = sig_path + ".tmp"
    with open(tmp_sig, "w") as fh:
        json.dump(signature, fh, indent=2)
    os.rename(tmp_sig, sig_path)

    print(f"Completed {out_tsv} in {elapsed:.2f}s ({elapsed/60.0:.2f}m)")
    if os.path.exists(tmp_dir):
        shutil.rmtree(tmp_dir)
    return elapsed


def parse_and_accumulate_hits(hit_tsv, query_to_pdb4, initial_cap=MAX_SEQS):
    """Stream raw Foldseek hits (21 columns with cigar), validate exact CIGAR aligned pairs
    on 100% of rows, and accumulate Tier 0-4 summaries per query in O(N_queries) memory.
    """
    acc = {}
    for qid, pdb4 in query_to_pdb4.items():
        acc[qid] = {
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

            if len(audited_examples) < 100:
                if (not pdb4 and fident >= 0.20) or (qspan - cigar_m >= 18 and len(audited_examples) < 50):
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
        if qid in priority_queries or len(audited_examples) < 250:
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


def audit_hit_cap_saturation(acc, label="query"):
    """Audit hit-cap saturation against each query's actual hit cap and identify
    saturated queries lacking MATCHED status at any tier (both POST_EXCLUSION_EMPTY
    and LOW_SCORE).
    """
    total_q = len(acc)
    # Check saturation against per-query hit cap
    saturated = [qid for qid, a in acc.items() if a["total_hits_raw"] >= a.get("hit_cap", MAX_SEQS)]

    empty_by_tier = {}
    low_score_by_tier = {}
    unmatched_by_tier = {}

    for k in range(5):
        empty_k = []
        low_k = []
        unm_k = []
        for qid in saturated:
            a = acc[qid]  # Explicitly index acc[qid] to avoid NameError
            status_k = classify_single_unit(a["total_hits_raw"], a["rem_counts"][k], a["top_qtm"][k])
            if status_k == "POST_EXCLUSION_EMPTY":
                empty_k.append(qid)
                unm_k.append(qid)
            elif status_k == "LOW_SCORE":
                low_k.append(qid)
                unm_k.append(qid)
        if empty_k:
            empty_by_tier[f"tier_{k}"] = sorted(empty_k)
        if low_k:
            low_score_by_tier[f"tier_{k}"] = sorted(low_k)
        if unm_k:
            unmatched_by_tier[f"tier_{k}"] = sorted(unm_k)

    all_prob_qids = sorted(list(set(
        qid for qids in unmatched_by_tier.values() for qid in qids
    )))

    return {
        "label": label,
        "total_queries": total_q,
        "saturated_count": len(saturated),
        "saturated_pct": round(100.0 * len(saturated) / max(1, total_q), 2),
        "problematic_exclusion_empty": empty_by_tier,
        "problematic_low_score": low_score_by_tier,
        "problematic_all_unmatched": unmatched_by_tier,
        "problematic_qids": all_prob_qids,
        "has_problematic_queries": len(all_prob_qids) > 0,
    }


def accumulate_hit_lines(lines, pdb4, hit_cap):
    """Accumulate raw Foldseek hit lines for a single query into an accumulator dict,
    validating CIGAR strings on every line.
    """
    qa = {
        "pdb4": pdb4,
        "hit_cap": hit_cap,
        "total_hits_raw": 0,
        "max_fident_raw": 0.0,
        "hits_fident_ge_20": 0,
        "hits_fident_ge_25": 0,
        "hits_fident_ge_30": 0,
        "rem_counts": [0, 0, 0, 0, 0],
        "top_qtm": [None, None, None, None, None],
        "top_alntm": [None, None, None, None, None],
    }
    for p in lines:
        if len(p) < 21:
            continue
        q = p[0][:-4] if p[0].endswith(".pdb") else p[0]
        target = p[1]
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

        qa["total_hits_raw"] += 1
        if fident > qa["max_fident_raw"]:
            qa["max_fident_raw"] = fident
        if fident >= 0.20:
            qa["hits_fident_ge_20"] += 1
        if fident >= 0.25:
            qa["hits_fident_ge_25"] += 1
        if fident >= 0.30:
            qa["hits_fident_ge_30"] += 1

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

        # Tier 0
        qa["rem_counts"][0] += 1
        if qa["top_qtm"][0] is None or qtm > qa["top_qtm"][0]["qtm"]:
            qa["top_qtm"][0] = row_dict
        if qa["top_alntm"][0] is None or alntm > qa["top_alntm"][0]["alntm"]:
            qa["top_alntm"][0] = row_dict

        # Tier 1
        if pdb4 and target[:4].lower() == pdb4:
            continue
        qa["rem_counts"][1] += 1
        if qa["top_qtm"][1] is None or qtm > qa["top_qtm"][1]["qtm"]:
            qa["top_qtm"][1] = row_dict
        if qa["top_alntm"][1] is None or alntm > qa["top_alntm"][1]["alntm"]:
            qa["top_alntm"][1] = row_dict

        # Tiers 2, 3, 4
        for t_idx, thresh in ((2, 0.30), (3, 0.25), (4, 0.20)):
            if fident < thresh:
                qa["rem_counts"][t_idx] += 1
                if qa["top_qtm"][t_idx] is None or qtm > qa["top_qtm"][t_idx]["qtm"]:
                    qa["top_qtm"][t_idx] = row_dict
                if qa["top_alntm"][t_idx] is None or alntm > qa["top_alntm"][t_idx]["alntm"]:
                    qa["top_alntm"][t_idx] = row_dict

    return qa


def rerun_problematic_queries(
    query_type,
    query_dir,
    db_path,
    orig_hits_tsv,
    prob_qids,
    acc,
    query_to_pdb4,
    threads,
    higher_max_seqs=RERUN_MAX_SEQS,
    reuse=False,
):
    """Rerun Foldseek with higher --max-seqs for problematic capped queries lacking MATCHED.

    Root-Audited Protocol:
    1. Separate query directory, raw hit file, and scratch directory by query type.
    2. Atomic completion signature validation on rerun searches.
    3. Full-row exact duplicate removal, retaining score-distinct alignments deterministically;
       preserves best scores independently at every exclusion tier.
    4. Saturation evaluates standalone rerun count against higher_max_seqs (not union count).
    5. Uniform return dictionary keys across rerun and no-rerun execution branches.
    """
    if not prob_qids:
        return {
            "rerun_executed": False,
            "query_type": query_type,
            "problematic_queries_count": 0,
            "queries_updated": 0,
            "queries_rescued_to_matched": 0,
            "records": [],
        }

    rerun_dir = os.path.join(H4_DIR, f"v7_{query_type}_rerun_queries")
    rerun_out = os.path.join(OUT_DIR, f"v7_{query_type}_rerun_hits.raw.tsv")
    rerun_tmp = os.path.join(H4_DIR, f"v7_{query_type}_rerun_fs_tmp")

    if not reuse or not os.path.exists(rerun_dir):
        if os.path.exists(rerun_dir):
            shutil.rmtree(rerun_dir)
        os.makedirs(rerun_dir, exist_ok=True)
        for qid in prob_qids:
            src = os.path.join(query_dir, f"{qid}.pdb")
            assert os.path.exists(src), f"Query source missing for rerun: {src}"
            os.symlink(os.path.realpath(src), os.path.join(rerun_dir, f"{qid}.pdb"))

    # Execute search with atomic completion signature
    run_foldseek_search_atomic(rerun_dir, db_path, rerun_out, rerun_tmp, threads, higher_max_seqs, reuse=reuse)

    # Read original hit lines for affected queries
    prob_set = set(prob_qids)
    orig_lines_by_q = defaultdict(list)
    if os.path.exists(orig_hits_tsv):
        with open(orig_hits_tsv) as fh:
            for line in fh:
                p = line.rstrip("\n").split("\t")
                if len(p) >= 21:
                    q = p[0][:-4] if p[0].endswith(".pdb") else p[0]
                    if q in prob_set:
                        orig_lines_by_q[q].append(p)

    # Read rerun hit lines for affected queries
    rerun_lines_by_q = defaultdict(list)
    if os.path.exists(rerun_out):
        with open(rerun_out) as fh:
            for line in fh:
                p = line.rstrip("\n").split("\t")
                if len(p) >= 21:
                    q = p[0][:-4] if p[0].endswith(".pdb") else p[0]
                    if q in prob_set:
                        rerun_lines_by_q[q].append(p)

    query_records = []
    rescued_count = 0
    updated_count = 0

    for qid in prob_qids:
        orig_qa = acc[qid]
        init_cap = orig_qa.get("hit_cap", MAX_SEQS)
        init_raw = orig_qa["total_hits_raw"]
        init_statuses = {
            f"tier_{k}": classify_single_unit(init_raw, orig_qa["rem_counts"][k], orig_qa["top_qtm"][k])
            for k in range(5)
        }
        init_top_qtm = {
            f"tier_{k}": round(orig_qa["top_qtm"][k]["qtm"], 4) if orig_qa["top_qtm"][k] else 0.0
            for k in range(5)
        }

        # Full-row exact duplicate removal: two rows are deduplicated ONLY if identical across all columns
        seen_full_rows = set()
        merged_lines = []
        for p in orig_lines_by_q.get(qid, []):
            row_key = tuple(p[:21])
            if row_key not in seen_full_rows:
                seen_full_rows.add(row_key)
                merged_lines.append(p)

        rerun_raw_lines = rerun_lines_by_q.get(qid, [])
        standalone_rerun_count = len(rerun_raw_lines)
        for p in rerun_raw_lines:
            row_key = tuple(p[:21])
            if row_key not in seen_full_rows:
                seen_full_rows.add(row_key)
                merged_lines.append(p)

        # Accumulate merged hits
        merged_qa = accumulate_hit_lines(merged_lines, query_to_pdb4.get(qid, ""), hit_cap=higher_max_seqs)

        # Score monotonicity: ensure scores never worsen compared to original at every tier
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

        final_raw = merged_qa["total_hits_raw"]
        final_statuses = {
            f"tier_{k}": classify_single_unit(final_raw, merged_qa["rem_counts"][k], merged_qa["top_qtm"][k])
            for k in range(5)
        }
        final_top_qtm = {
            f"tier_{k}": round(merged_qa["top_qtm"][k]["qtm"], 4) if merged_qa["top_qtm"][k] else 0.0
            for k in range(5)
        }

        improved = any(init_statuses[f"tier_{k}"] != "MATCHED" and final_statuses[f"tier_{k}"] == "MATCHED" for k in range(5))
        if improved:
            rescued_count += 1
        if final_raw > init_raw or any(final_top_qtm[f"tier_{k}"] > init_top_qtm[f"tier_{k}"] for k in range(5)):
            updated_count += 1

        acc[qid] = merged_qa

        # Standalone rerun saturation check (not union count)
        query_records.append({
            "query_id": qid,
            "initial_cap": init_cap,
            "final_cap": higher_max_seqs,
            "initial_total_hits_raw": init_raw,
            "standalone_rerun_hits_raw": standalone_rerun_count,
            "merged_unique_hits_raw": final_raw,
            "initial_statuses": init_statuses,
            "final_statuses": final_statuses,
            "initial_top_qtm": init_top_qtm,
            "final_top_qtm": final_top_qtm,
            "score_monotonicity_verified": all(final_top_qtm[f"tier_{k}"] >= init_top_qtm[f"tier_{k}"] - 1e-9 for k in range(5)),
            "status_improved_to_matched": improved,
            "remaining_saturated_at_final_cap": standalone_rerun_count >= higher_max_seqs,
        })

    if os.path.exists(rerun_tmp):
        shutil.rmtree(rerun_tmp)

    return {
        "rerun_executed": True,
        "query_type": query_type,
        "problematic_queries_count": len(prob_qids),
        "queries_updated": updated_count,
        "queries_rescued_to_matched": rescued_count,
        "records": query_records,
    }


def summarize_stratum_census(models_sub, regs_sub):
    """Compute exact finite-population stratum census metrics across all Tiers 0-4.
    Undefined domains return None (null in JSON). Primary census does NOT contain
    sampling-error Wilson intervals.
    """
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
        if t1_ret_regs else None
    )
    max_t1_fident_ret_only = (
        round(100.0 * max(r["t1_fident"] for r in t1_ret_regs), 2)
        if t1_ret_regs else None
    )
    regs_t1_fident_ge20 = sum(1 for r in t1_ret_regs if r["t1_fident"] >= 0.20)
    regs_t1_fident_ge25 = sum(1 for r in t1_ret_regs if r["t1_fident"] >= 0.25)
    regs_t1_fident_ge30 = sum(1 for r in t1_ret_regs if r["t1_fident"] >= 0.30)

    mean_t1_aligned_pair_qcov = (
        round(100.0 * sum(r["t1_aligned_pair_qcov"] for r in t1_ret_regs) / len(t1_ret_regs), 2)
        if t1_ret_regs else None
    )
    mean_t1_qspan_cov = (
        round(100.0 * sum(r["t1_qspan_cov"] for r in t1_ret_regs) / len(t1_ret_regs), 2)
        if t1_ret_regs else None
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
            "exact_census_rates": {
                "unconditional_crop_match_count": n_mod_match,
                "unconditional_crop_match_denominator": n_models,
                "unconditional_crop_match_pct": safe_rate(n_mod_match, n_models),
                "region_eligible_crop_match_count": n_mod_match,
                "region_eligible_crop_match_denominator": n_elig,
                "region_eligible_crop_match_pct": safe_rate(n_mod_match, n_elig),
                "remaining_hit_conditional_crop_match_count": n_mod_match,
                "remaining_hit_conditional_crop_match_denominator": n_mod_rem,
                "remaining_hit_conditional_crop_match_pct": safe_rate(n_mod_match, n_mod_rem),
                "wc_unmatched_count": n_wc_unm,
                "wc_unmatched_recovery_unconditional_count": n_wc_unm_rec,
                "wc_unmatched_recovery_unconditional_pct": safe_rate(n_wc_unm_rec, n_wc_unm),
                "wc_unmatched_eligible_count": n_wc_unm_elig,
                "wc_unmatched_recovery_eligible_count": n_wc_unm_rec,
                "wc_unmatched_recovery_eligible_pct": safe_rate(n_wc_unm_rec, n_wc_unm_elig),
                "wc_matched_count": n_wc_m,
                "wc_matched_retained_count": n_wc_m_ret,
                "wc_matched_retained_pct": safe_rate(n_wc_m_ret, n_wc_m),
                "crop_persistent_low_score_pct": safe_rate(crop_stat_counts["PERSISTENT_LOW_SCORE"], n_models),
                "crop_post_exclusion_empty_pct": safe_rate(crop_stat_counts["POST_EXCLUSION_EMPTY"], n_models),
                "crop_raw_zero_retrieval_pct": safe_rate(crop_stat_counts["RAW_ZERO_RETRIEVAL"], n_models),
                "crop_no_confident_region_pct": safe_rate(crop_stat_counts["NO_CONFIDENT_REGION"], n_models),
            },
            "region_status_counts": {
                "MATCHED": reg_stat_counts["MATCHED"],
                "LOW_SCORE": reg_stat_counts["LOW_SCORE"],
                "POST_EXCLUSION_EMPTY": reg_stat_counts["POST_EXCLUSION_EMPTY"],
                "RAW_ZERO_RETRIEVAL": reg_stat_counts["RAW_ZERO_RETRIEVAL"],
            },
            "region_rates": {
                "total_regions": n_regs,
                "matched_regions": n_reg_match,
                "unconditional_match_pct": safe_rate(n_reg_match, n_regs),
                "remaining_hit_regions": n_reg_rem,
                "remaining_hit_conditional_match_pct": safe_rate(n_reg_match, n_reg_rem),
            },
            "model_based_superpopulation_inferential_statistics": {
                "disclaimer": "Wilson intervals assume a hypothetical superpopulation Bernoulli model under independence assumptions, not finite-population uncertainty.",
                "unconditional_crop_match_wilson_ci95": wilson_ci_model_based(n_mod_match, n_models),
                "region_eligible_crop_match_wilson_ci95": wilson_ci_model_based(n_mod_match, n_elig),
                "remaining_hit_conditional_crop_match_wilson_ci95": wilson_ci_model_based(n_mod_match, n_mod_rem),
                "wc_unmatched_recovery_unconditional_wilson_ci95": wilson_ci_model_based(n_wc_unm_rec, n_wc_unm),
                "wc_unmatched_recovery_eligible_wilson_ci95": wilson_ci_model_based(n_wc_unm_rec, n_wc_unm_elig),
                "wc_matched_retained_wilson_ci95": wilson_ci_model_based(n_wc_m_ret, n_wc_m),
            },
        }

    return {
        "model_count": n_models,
        "models_with_confident_regions": n_elig,
        "models_no_confident_region": n_no_reg,
        "pct_no_confident_region": safe_rate(n_no_reg, n_models),
        "wc_raw_zero_retrieval_models": wc_raw_zero,
        "wc_raw_saturated_models": wc_raw_sat,
        "total_confident_regions": n_regs,
        "mean_regions_per_model": round(n_regs / max(1, n_models), 3),
        "mean_regions_per_eligible_model": round(n_regs / max(1, n_elig), 3),
        "mean_region_length_aa": round(sum(r["reg_len"] for r in regs_sub) / max(1, n_regs), 2) if regs_sub else None,
        "mean_region_plddt": round(sum(r["reg_mean_plddt"] for r in regs_sub) / max(1, n_regs), 4) if regs_sub else None,
        "reg_raw_zero_retrieval_count": reg_raw_zero,
        "reg_raw_saturated_count": reg_raw_sat,
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


def compute_effect_sizes(strata_summaries, tier_k=1):
    """Compute main-paper-ready effect sizes between candidate strata and controls."""
    tk = f"tier_{tier_k}"
    k_loc = strata_summaries["K_local_pdb"]["tiers"][tk]["exact_census_rates"]
    k_all = strata_summaries["K_local_all"]["tiers"][tk]["exact_census_rates"]
    k_none = strata_summaries["K_dual_none"]["tiers"][tk]["exact_census_rates"]
    ctrl = strata_summaries["Controls"]["tiers"][tk]["exact_census_rates"]

    def calc_2x2_stats(k1, n1, k2, n2):
        if n1 <= 0 or n2 <= 0:
            return None
        p1 = k1 / n1
        p2 = k2 / n2
        diff_pct = round(100.0 * (p1 - p2), 2)
        se_diff = math.sqrt((p1 * (1 - p1) / n1) + (p2 * (1 - p2) / n2))
        diff_ci95 = [round(100.0 * (p1 - p2 - 1.959964 * se_diff), 2), round(100.0 * (p1 - p2 + 1.959964 * se_diff), 2)]
        rr = round(p1 / max(1e-9, p2), 4) if p2 > 0 else None

        table = [[k1, n1 - k1], [k2, n2 - k2]]
        try:
            import scipy.stats
            odds_ratio, p_fisher = scipy.stats.fisher_exact(table)
            chi2, p_chi2, _, _ = scipy.stats.chi2_contingency(table, correction=True)
            p_fisher_str = float(f"{p_fisher:.6e}")
            p_chi2_str = float(f"{p_chi2:.6e}")
            or_val = round(odds_ratio, 4) if not math.isinf(odds_ratio) else "inf"
        except ImportError:
            or_val = "n/a"
            p_fisher_str = "n/a"
            p_chi2_str = "n/a"

        return {
            "group1_rate_pct": round(100.0 * p1, 2),
            "group2_rate_pct": round(100.0 * p2, 2),
            "rate_difference_percentage_points": diff_pct,
            "risk_ratio": rr,
            "odds_ratio": or_val,
            "model_based_superpopulation_inferential_stats": {
                "disclaimer": "Assumes independent Bernoulli observations from a superpopulation, not finite-census sampling error.",
                "rate_difference_ci95": diff_ci95,
                "fisher_exact_p": p_fisher_str,
                "chi2_p": p_chi2_str,
            }
        }

    loc_vs_none_uncond = calc_2x2_stats(
        k_loc["unconditional_crop_match_count"], k_loc["unconditional_crop_match_denominator"],
        k_none["unconditional_crop_match_count"], k_none["unconditional_crop_match_denominator"]
    )
    loc_vs_none_rec = calc_2x2_stats(
        k_loc["wc_unmatched_recovery_unconditional_count"], k_loc["wc_unmatched_count"],
        k_none["wc_unmatched_recovery_unconditional_count"], k_none["wc_unmatched_count"]
    )
    all_vs_none_uncond = calc_2x2_stats(
        k_all["unconditional_crop_match_count"], k_all["unconditional_crop_match_denominator"],
        k_none["unconditional_crop_match_count"], k_none["unconditional_crop_match_denominator"]
    )
    all_vs_none_rec = calc_2x2_stats(
        k_all["wc_unmatched_recovery_unconditional_count"], k_all["wc_unmatched_count"],
        k_none["wc_unmatched_recovery_unconditional_count"], k_none["wc_unmatched_count"]
    )

    return {
        "tier_evaluated": tk,
        "k_local_pdb_vs_k_dual_none_unconditional": loc_vs_none_uncond,
        "k_local_pdb_vs_k_dual_none_wc_unmatched_recovery": loc_vs_none_rec,
        "k_local_all_vs_k_dual_none_unconditional": all_vs_none_uncond,
        "k_local_all_vs_k_dual_none_wc_unmatched_recovery": all_vs_none_rec,
    }


def compare_with_round6_pilot(full_strata, pilot_models, pilot_regs):
    """Compare the full N=5,851 census against round6 pilot estimates."""
    pilot_strata = {}
    for st in SUB_STRATA_ORDER:
        m_sub = [m for m in pilot_models if m["sub_stratum"] == st]
        r_sub = [r for r in pilot_regs if r["sub_stratum"] == st]
        pilot_strata[st] = summarize_stratum_census(m_sub, r_sub)

    for comb_name, sub_list in COMBINED_STRATA_MAP.items():
        m_sub = [m for m in pilot_models if m["sub_stratum"] in sub_list]
        r_sub = [r for r in pilot_regs if r["sub_stratum"] in sub_list]
        pilot_strata[comb_name] = summarize_stratum_census(m_sub, r_sub)

    comparison = {}
    for st in SUB_STRATA_ORDER + list(COMBINED_STRATA_MAP.keys()):
        fc_t1 = full_strata[st]["tiers"]["tier_1"]["exact_census_rates"]
        pt_t1 = pilot_strata[st]["tiers"]["tier_1"]["exact_census_rates"]

        comparison[st] = {
            "population_models": full_strata[st]["model_count"],
            "pilot_sample_models": pilot_strata[st]["model_count"],
            "full_census_uncond_crop_match_pct": fc_t1["unconditional_crop_match_pct"],
            "pilot_unweighted_uncond_crop_match_pct": pt_t1["unconditional_crop_match_pct"],
            "difference_full_minus_pilot_pct": (
                round(fc_t1["unconditional_crop_match_pct"] - pt_t1["unconditional_crop_match_pct"], 2)
                if fc_t1["unconditional_crop_match_pct"] is not None and pt_t1["unconditional_crop_match_pct"] is not None
                else None
            ),
            "full_census_wc_unmatched_recovery_pct": fc_t1["wc_unmatched_recovery_unconditional_pct"],
            "pilot_wc_unmatched_recovery_pct": pt_t1["wc_unmatched_recovery_unconditional_pct"],
        }
    return comparison


def run_closure_assertions(full_cohort, model_summaries, region_evaluations, stratum_summaries):
    """Run strict programmatic code-data closure assertions for full cohort."""
    assert len(full_cohort) == 5851, f"Expected 5851 models, got {len(full_cohort)}"
    assert len(model_summaries) == 5851, f"Expected 5851 model summaries, got {len(model_summaries)}"
    assert len(region_evaluations) == 7773, f"Expected 7773 cropped regions, got {len(region_evaluations)}"

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

    for r in region_evaluations:
        for k in range(4):
            assert r[f"t{k}_remaining_hits"] >= r[f"t{k+1}_remaining_hits"]
            assert r[f"t{k}_qtm"] + 1e-9 >= r[f"t{k+1}_qtm"]
        for k in (0, 1):
            assert r[f"t{k}_aligned_pairs"] <= r[f"t{k}_qspan"] <= r["reg_len"]
            assert r[f"t{k}_aligned_pair_qcov"] <= r[f"t{k}_qspan_cov"] + 1e-6

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

            sp = tk["model_based_superpopulation_inferential_statistics"]
            for ci_key in [
                "unconditional_crop_match_wilson_ci95",
                "region_eligible_crop_match_wilson_ci95",
                "wc_unmatched_recovery_unconditional_wilson_ci95",
                "wc_unmatched_recovery_eligible_wilson_ci95",
                "wc_matched_retained_wilson_ci95",
            ]:
                val = sp[ci_key]
                if val is not None:
                    assert val != [0.0, 0.0], f"Found misleading [0.0, 0.0] interval in {st_name} {k} {ci_key}"
                    assert val != [100.0, 100.0], f"Found misleading [100.0, 100.0] interval in {st_name} {k} {ci_key}"

    print("ALL V7 CODE-DATA CLOSURE ASSERTIONS PASSED SUCCESSFULLY.")


def generate_sha256_manifest(res_dir, manifest_name="v7_results.sha256.tsv", include_receipt=True):
    """Generate SHA-256 manifest of results, atomically excluding the manifest itself."""
    entries = []
    for fn in sorted(os.listdir(res_dir)):
        if (fn.startswith("v7_") and (fn.endswith(".tsv") or fn.endswith(".json"))) and fn != manifest_name and not fn.endswith(".tmp"):
            if not include_receipt and fn == "v7_closure_verification.json":
                continue
            fpath = os.path.join(res_dir, fn)
            if os.path.isfile(fpath):
                entries.append((sha256_file(fpath), fn))

    tmp_manifest = os.path.join(res_dir, manifest_name + ".tmp")
    with open(tmp_manifest, "w") as fh:
        for digest, fn in entries:
            fh.write(f"{digest}  {fn}\n")
    os.rename(tmp_manifest, os.path.join(res_dir, manifest_name))
    print(f"Generated SHA-256 manifest with {len(entries)} entries at {manifest_name}")
    return entries


def main():
    parser = argparse.ArgumentParser(description="H4 V7 Full-Cohort Cropped & Whole-Chain Search Pipeline")
    parser.add_argument("--mode", choices=["full", "pilot"], default="full", help="Execution mode: full (N=5,851) or pilot (N=300)")
    parser.add_argument("--threads", type=int, default=32, help="Foldseek search threads")
    parser.add_argument("--max-seqs", type=int, default=MAX_SEQS, help="Foldseek --max-seqs prefilter cap")
    parser.add_argument("--reuse-search", action="store_true", help="Reuse existing v7 raw Foldseek outputs if verified")
    args = parser.parse_args()

    # Step 0: Validate baseline reference identity
    print("=== Step 0: Validating Baseline Reference Identity & Inputs ===")
    ref_status = validate_baseline_references()
    for ref_name, st in ref_status.items():
        print(f"  {ref_name}: {st['status']}")

    print(f"\n=== Step 1: Loading Full Finite Cohort (N=5,851) in Mode={args.mode} ===")
    full_cohort = load_full_cohort()
    if args.mode == "pilot":
        active_cohort = {qid: m for qid, m in full_cohort.items() if m["in_v2_pilot"] == 1}
        wc_dir = os.path.join(H4_DIR, "v7_wc_pilot_queries")
        crops_dir = os.path.join(H4_DIR, "v7_pilot_crops")
    else:
        active_cohort = full_cohort
        wc_dir = os.path.join(H4_DIR, "v7_wc_queries")
        crops_dir = os.path.join(H4_DIR, "v7_crops")

    print(f"Active cohort size: {len(active_cohort)} models (out of {len(full_cohort)} finite population).")
    for st in SUB_STRATA_ORDER:
        cnt = sum(1 for m in active_cohort.values() if m["sub_stratum"] == st)
        pop_cnt = sum(1 for m in full_cohort.values() if m["sub_stratum"] == st)
        print(f"  {st:20s}: sample={cnt:4d} / pop={pop_cnt:4d}")

    print("\n=== Step 2: Preparing Whole-Chain Symlinks & Cropping Confident Regions (>=50 aa, pLDDT >= 0.70) ===")
    if not args.reuse_search or not os.path.exists(wc_dir):
        if os.path.exists(wc_dir):
            shutil.rmtree(wc_dir)
        os.makedirs(wc_dir, exist_ok=True)
    if not args.reuse_search or not os.path.exists(crops_dir):
        if os.path.exists(crops_dir):
            shutil.rmtree(crops_dir)
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
        if not os.path.exists(wc_dst):
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
            crop_query_to_pdb4[crop_id] = m["pdb4"]
            if not os.path.exists(crop_path):
                _, crop_hash = write_cropped_pdb(real_src, crop_path, reg["start_res"], reg["end_res"], qid, reg["reg_id"])
            else:
                crop_hash = sha256_file(crop_path)
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

    print(f"Created/verified {total_crops} cropped region PDBs from {len(active_cohort) - models_no_region} models.")
    print(f"Models with 0 confident regions (>=50 aa, smoothed pLDDT >= 0.70): {models_no_region}")

    manifest_out = os.path.join(OUT_DIR, "v7_full5851_cropped_regions_manifest.tsv" if args.mode == "full" else "v7_pilot_cropped_regions_manifest.tsv")
    with open(manifest_out, "w") as fh:
        w = csv.DictWriter(fh, fieldnames=list(model_region_manifest[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(model_region_manifest)

    print(f"\n=== Step 3: Same-Run Foldseek Searches (Whole-Chain & Cropped Regions, --max-seqs {args.max_seqs}, --threads {args.threads}) ===")
    wc_hits_path = os.path.join(OUT_DIR, "v7_wc_pdb_hits.raw.tsv")
    crop_hits_path = os.path.join(OUT_DIR, "v7_crop_pdb_hits.raw.tsv")
    wc_tmp = os.path.join(H4_DIR, "v7_fs_wc_tmp")
    crop_tmp = os.path.join(H4_DIR, "v7_fs_crop_tmp")

    wc_elapsed = run_foldseek_search_atomic(wc_dir, PDB_DB, wc_hits_path, wc_tmp, args.threads, args.max_seqs, reuse=args.reuse_search)
    crop_elapsed = run_foldseek_search_atomic(crops_dir, PDB_DB, crop_hits_path, crop_tmp, args.threads, args.max_seqs, reuse=args.reuse_search)

    print("\n=== Step 4: Streaming & Validating Whole-Chain and Cropped Hits (CIGAR & Tiers 0-4) ===")
    wc_acc, wc_hit_stats, wc_examples = parse_and_accumulate_hits(wc_hits_path, wc_query_to_pdb4, initial_cap=args.max_seqs)
    crop_acc, crop_hit_stats, crop_examples = parse_and_accumulate_hits(crop_hits_path, crop_query_to_pdb4, initial_cap=args.max_seqs)
    print(f"Validated {wc_hit_stats['total_rows']} whole-chain hit rows and {crop_hit_stats['total_rows']} cropped hit rows.")

    print("\n=== Step 5: Auditing Hit-Cap Saturation & Checking for Exclusion-Empty Capped Queries ===")
    wc_sat_audit_initial = audit_hit_cap_saturation(wc_acc, label="whole_chain_initial")
    crop_sat_audit_initial = audit_hit_cap_saturation(crop_acc, label="cropped_regions_initial")
    print(f"Whole-chain saturated queries (>=1000): {wc_sat_audit_initial['saturated_count']} ({wc_sat_audit_initial['saturated_pct']}%)")
    print(f"Cropped-region saturated queries (>=1000): {crop_sat_audit_initial['saturated_count']} ({crop_sat_audit_initial['saturated_pct']}%)")

    # Adaptive Rerun with --max-seqs 5000 if problematic capped queries exist
    wc_rerun_res = rerun_problematic_queries(
        "wc", wc_dir, PDB_DB, wc_hits_path, wc_sat_audit_initial["problematic_qids"],
        wc_acc, wc_query_to_pdb4, args.threads, higher_max_seqs=RERUN_MAX_SEQS, reuse=args.reuse_search
    )
    crop_rerun_res = rerun_problematic_queries(
        "crop", crops_dir, PDB_DB, crop_hits_path, crop_sat_audit_initial["problematic_qids"],
        crop_acc, crop_query_to_pdb4, args.threads, higher_max_seqs=RERUN_MAX_SEQS, reuse=args.reuse_search
    )

    rerun_audit = {
        "whole_chain_rerun": wc_rerun_res,
        "cropped_rerun": crop_rerun_res,
        "any_rerun_executed": wc_rerun_res["rerun_executed"] or crop_rerun_res["rerun_executed"],
    }
    if rerun_audit["any_rerun_executed"]:
        print(f"Adaptive reruns complete: WC updated={wc_rerun_res['queries_updated']}, rescued={wc_rerun_res['queries_rescued_to_matched']} | Crop updated={crop_rerun_res['queries_updated']}, rescued={crop_rerun_res['queries_rescued_to_matched']}")
    else:
        print("No capped queries lacking MATCHED detected. All saturated queries retained valid matches across exclusions.")

    # Record post-rerun final saturation audit
    wc_sat_audit_final = audit_hit_cap_saturation(wc_acc, label="whole_chain_final")
    crop_sat_audit_final = audit_hit_cap_saturation(crop_acc, label="cropped_regions_final")

    print("\n=== Step 6: Building Region Evaluations and Model Summaries across Tiers 0-4 ===")
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
        raw_sat = int(raw_tot >= ca_info.get("hit_cap", args.max_seqs))

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
            "hit_cap": ca_info.get("hit_cap", args.max_seqs),
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
        wc_sat = int(wc_raw >= wa.get("hit_cap", args.max_seqs))
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
            "total_confident_aa": m["total_confident_aa"],
            "confident_aa_fraction": m["confident_aa_fraction"],
            "wc_total_hits_raw": wc_raw,
            "wc_raw_list_saturated": wc_sat,
            "wc_hit_cap": wa.get("hit_cap", args.max_seqs),
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

    wc_raw_hash = sha256_file(wc_hits_path)
    crop_raw_hash = sha256_file(crop_hits_path)

    raw_provenance = {
        "v7_wc_pdb_hits.raw.tsv": {
            "size_bytes": os.path.getsize(wc_hits_path),
            "sha256": wc_raw_hash,
            **wc_hit_stats,
        },
        "v7_crop_pdb_hits.raw.tsv": {
            "size_bytes": os.path.getsize(crop_hits_path),
            "sha256": crop_raw_hash,
            **crop_hit_stats,
        },
    }
    if os.path.exists(os.path.join(OUT_DIR, "v7_wc_rerun_hits.raw.tsv")):
        p = os.path.join(OUT_DIR, "v7_wc_rerun_hits.raw.tsv")
        raw_provenance["v7_wc_rerun_hits.raw.tsv"] = {
            "size_bytes": os.path.getsize(p),
            "sha256": sha256_file(p),
        }
    if os.path.exists(os.path.join(OUT_DIR, "v7_crop_rerun_hits.raw.tsv")):
        p = os.path.join(OUT_DIR, "v7_crop_rerun_hits.raw.tsv")
        raw_provenance["v7_crop_rerun_hits.raw.tsv"] = {
            "size_bytes": os.path.getsize(p),
            "sha256": sha256_file(p),
        }

    summary_json = {
        "analysis_version": "v7_full_cohort_census_symmetric_whole_chain_and_cropped_search",
        "cohort_mode": args.mode,
        "baseline_reference_verification": ref_status,
        "finite_population_models": len(full_cohort),
        "evaluated_models": len(model_summaries),
        "evaluated_cropped_regions": len(region_evaluations),
        "models_no_confident_region": models_no_region,
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
        "hit_cap_saturation_audit": {
            "initial": {
                "whole_chain": wc_sat_audit_initial,
                "cropped_regions": crop_sat_audit_initial,
            },
            "final": {
                "whole_chain": wc_sat_audit_final,
                "cropped_regions": crop_sat_audit_final,
            },
            "adaptive_rerun_audit": rerun_audit,
        },
        "raw_hit_files_provenance": raw_provenance,
        "effect_sizes": effect_sizes,
        "pilot_comparison": pilot_comp,
        "stratum_summaries": stratum_summaries,
    }

    with open(os.path.join(OUT_DIR, "v7_recovery_summary.json"), "w") as fh:
        json.dump(summary_json, fh, indent=2)

    print("\n=== SUMMARY OF V7 FULL CENSUS RESULTS (Tier 1: Source-PDB Excluded) ===")
    for st in SUB_STRATA_ORDER + list(COMBINED_STRATA_MAP.keys()):
        sd = stratum_summaries[st]
        t1 = sd["tiers"]["tier_1"]
        mr = t1["exact_census_rates"]
        print(
            f"{st:20s} | N={sd['model_count']:4d} (Elig={sd['models_with_confident_regions']:4d}, "
            f"0Reg={sd['models_no_confident_region']:3d}) | "
            f"CropMatch={mr['unconditional_crop_match_count']:4d}/{sd['model_count']} = {mr['unconditional_crop_match_pct']}% | "
            f"EligMatch={mr['region_eligible_crop_match_pct']}% | "
            f"WC_Unm={mr['wc_unmatched_count']:4d} -> Rec={mr['wc_unmatched_recovery_unconditional_count']:4d} "
            f"({mr['wc_unmatched_recovery_unconditional_pct']}%)"
        )

    # Manifest and closure verification sequence (Root-Audited):
    # 1. Generate provisional manifest (excluding closure verification receipt)
    print("\n=== Step 9: Provisional Manifest Generation & Closure Verification ===")
    generate_sha256_manifest(OUT_DIR, manifest_name="v7_results.sha256.tsv", include_receipt=False)

    # 2. Run closure verification in provisional mode (writes v7_closure_verification.json)
    from v7_verify_closure import verify_closure
    verify_closure(OUT_DIR, mode="provisional")

    # 3. Generate final manifest including receipt
    generate_sha256_manifest(OUT_DIR, manifest_name="v7_results.sha256.tsv", include_receipt=True)

    # 4. Strict check all files against final manifest
    verify_closure(OUT_DIR, mode="strict")

    print(f"\nAll V7 outputs successfully written, verified, and sealed in {OUT_DIR}")


if __name__ == "__main__":
    main()
