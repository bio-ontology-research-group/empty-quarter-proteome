#!/usr/bin/env python3
"""
scripts/08_f3_rhodopsin_full.py
Track F3: Microbial & Heliorhodopsin Full Discovery across:
1. All 5,438 Empty Quarter representative proteomes (cohort.faa, 20,451,986 proteins)
2. Harmonized Environmental SMAG cohort (smag_pool_canonical.faa)
3. Harmonized Environmental GEM cohort (gem_harmonized_pool_canonical_seed20261014.faa)

Strict Provenance & Audit Resolutions (Argon Section 4 & Root Requirements):
- Reuses verified pilot profile collection: references/profiles/f3_rhodopsins.hmm
  (PF01036, PF18761, NF038020, PF00001, PF00002).
- Never infer complete CDS from length>=180 or heuristic. Keep UNKNOWN distinct.
  For EQ MAGs: parse partial=00 from Prodigal header.
  For SMAG/GEM: look up partial from source coordinate mapping TSVs.
  Only partial == "00" is flagged is_complete_gene = 1.
- GTDB-Tk R232 taxonomies: loads both gtdbtk.bac120.summary.tsv and gtdbtk.ar53.summary.tsv.
- One-way DIAMOND blastp (--ultra-sensitive) against references/all_f345_references.dmnd.
  Candidate FASTA ID is bare protein_id. Query and target coverage calculated from coordinate spans.
- MMseqs2 clustering (prefix_cluster.tsv) via envs/mmseqs2/bin/mmseqs. Fail-closed.
- Sequence-level deduplication across cohorts preserving cross-source memberships.
- ProstT5 inputs strictly separated by actual family (Type-1 vs Heliorhodopsin), including curated references.
- Retains all rejected and ambiguous hits with explicit reasons.
- Emits results to results/full_v1/f3_rhodopsin/ and freezes executed code.
"""

import csv
import gzip
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent.parent
REFS_DIR = BASE_DIR / "references"
PROFILES_DIR = REFS_DIR / "profiles"
COHORTS_DIR = BASE_DIR / "cohorts"
RESULTS_DIR = BASE_DIR / "results" / "full_v1"
F3_DIR = RESULTS_DIR / "f3_rhodopsin"
PROSTT5_DIR = RESULTS_DIR / "prostt5_inputs"
FROZEN_DIR = RESULTS_DIR / "frozen_code"

for d in [F3_DIR, PROSTT5_DIR, FROZEN_DIR]:
    d.mkdir(parents=True, exist_ok=True)

# Profile constants (matching pilot exactly)
F3_HMM = PROFILES_DIR / "f3_rhodopsins.hmm"
TYPE1_ACCS = {"PF01036"}
HELIO_ACCS = {"PF18761", "NF038020"}
GPCR_ACCS = {"PF00001", "PF00002"}

# References & Tools
ALL_REFS_DMND = REFS_DIR / "all_f345_references.dmnd"
ALL_REFS_FAA = REFS_DIR / "all_f345_references.faa"
REFS_META_TSV = REFS_DIR / "verified_references_metadata.tsv"

MMSEQS_BIN = Path("/ibex/scratch/projects/c2014/rob/miniconda3/envs/mmseqs2/bin/mmseqs")
DIAMOND_BIN = Path("/ibex/scratch/projects/c2014/rob/miniconda3/envs/eggnog/bin/diamond")

W = Path("/ibex/scratch/projects/c2014/rob/metagenomics")
WDB_CSV = W / "05_refine/_drep/out/data_tables/Wdb.csv"
EQ_COHORT_FAA = W / "06_darkmatter/family_discovery_20261015_clean/root_audit/f1_full/cohort.faa"
EQ_COHORT_TSV = W / "06_darkmatter/family_discovery_20261015_clean/root_audit/f1_full/cohort.tsv"

H3_DIR = W / "06_darkmatter/p7_analyses_clean/h3_references/results"
SMAG_CANONICAL = H3_DIR / "smag_pool_canonical.faa"
SMAG_COORD_MAP = H3_DIR / "smag_protein_coordinate_mapping.tsv"
GEM_CANONICAL = H3_DIR / "gem_harmonized_pool_canonical_seed20261014.faa"
GEM_COORD_MAP = H3_DIR / "gem_harmonized_protein_coordinate_mapping_seed20261014.tsv"

GTDB_BAC = W / "05_refine/_gtdbtk_r232/gtdbtk.bac120.summary.tsv"
GTDB_AR = W / "05_refine/_gtdbtk_r232/gtdbtk.ar53.summary.tsv"


def load_metadata():
    print("Loading Wdb.csv and GTDB-Tk R232 taxonomies...")
    assert WDB_CSV.exists(), f"CRITICAL: Missing Wdb.csv: {WDB_CSV}"
    mag_to_cluster = {}
    with open(WDB_CSV, "r", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            mag = r["genome"].removesuffix(".fa").removesuffix(".fasta")
            mag_to_cluster[mag] = r.get("cluster", "Unknown")

    tax_map = {}
    for p in [GTDB_BAC, GTDB_AR]:
        assert p.exists(), f"CRITICAL: Missing GTDB-Tk summary: {p}"
        with open(p, "r", encoding="utf-8") as f:
            for r in csv.DictReader(f, delimiter="\t"):
                ug = r.get("user_genome", "").removesuffix(".fa").removesuffix(".fasta")
                cls = r.get("classification", "").strip()
                if ug and cls:
                    tax_map[ug] = cls

    print(f"Loaded {len(mag_to_cluster)} clusters and {len(tax_map)} GTDB-Tk taxonomies.")
    return mag_to_cluster, tax_map


def parse_domtblout(domtbl_path):
    hits = []
    with open(domtbl_path, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split()
            if len(parts) < 23:
                continue
            target_name = parts[0]
            target_len = int(parts[2])
            profile_name = parts[3]
            raw_acc = parts[4]
            profile_acc = raw_acc.split(".")[0]
            qlen = int(parts[5])
            full_evalue = float(parts[6])
            full_score = float(parts[7])
            dom_ievalue = float(parts[12])
            dom_score = float(parts[13])
            hmm_from = int(parts[15])
            hmm_to = int(parts[16])
            ali_from = int(parts[17])
            ali_to = int(parts[18])
            hmm_cov = (hmm_to - hmm_from + 1) / float(qlen) if qlen > 0 else 0.0

            hits.append({
                "target_name": target_name,
                "target_len": target_len,
                "profile_name": profile_name,
                "profile_acc": profile_acc,
                "profile_acc_raw": raw_acc,
                "qlen": qlen,
                "full_evalue": full_evalue,
                "full_score": full_score,
                "dom_ievalue": dom_ievalue,
                "dom_score": dom_score,
                "hmm_from": hmm_from,
                "hmm_to": hmm_to,
                "ali_from": ali_from,
                "ali_to": ali_to,
                "hmm_coverage": hmm_cov,
            })
    return hits


def evaluate_f3_gates(hits_by_protein):
    candidates = {}
    rejections = []

    for prot_id, prot_hits in hits_by_protein.items():
        type1_hits = [h for h in prot_hits if h["profile_acc"] in TYPE1_ACCS]
        helio_hits = [h for h in prot_hits if h["profile_acc"] in HELIO_ACCS]
        gpcr_hits = [h for h in prot_hits if h["profile_acc"] in GPCR_ACCS]

        best_type1 = max(type1_hits, key=lambda x: x["dom_score"]) if type1_hits else None
        best_helio = max(helio_hits, key=lambda x: x["dom_score"]) if helio_hits else None
        best_gpcr = max(gpcr_hits, key=lambda x: x["dom_score"]) if gpcr_hits else None

        score_type1 = best_type1["dom_score"] if best_type1 else 0.0
        score_helio = best_helio["dom_score"] if best_helio else 0.0
        score_gpcr = best_gpcr["dom_score"] if best_gpcr else 0.0

        target_score = max(score_type1, score_helio)
        best_target = best_type1 if score_type1 >= score_helio else best_helio

        if not best_target:
            rejections.append({
                "protein_id": prot_id,
                "reason": "Only_GPCR_control_hit",
                "target_score": 0.0,
                "gpcr_score": score_gpcr,
                "hmm_coverage": 0.0,
                "dom_ievalue": 1.0,
                "profile_acc": best_gpcr["profile_acc"] if best_gpcr else "None",
            })
            continue

        cov = best_target["hmm_coverage"]
        ie = best_target["dom_ievalue"]

        # Gate 1: dom_ievalue <= 1e-5
        if ie > 1e-5:
            rejections.append({
                "protein_id": prot_id,
                "reason": f"Failed_dom_ievalue_gate_{ie:.2e}>1e-5",
                "target_score": target_score,
                "gpcr_score": score_gpcr,
                "hmm_coverage": cov,
                "dom_ievalue": ie,
                "profile_acc": best_target["profile_acc"],
            })
            continue

        # Gate 2: domain coverage >= 0.70
        if cov < 0.70:
            rejections.append({
                "protein_id": prot_id,
                "reason": f"Failed_domain_coverage_gate_{cov:.2f}<0.70",
                "target_score": target_score,
                "gpcr_score": score_gpcr,
                "hmm_coverage": cov,
                "dom_ievalue": ie,
                "profile_acc": best_target["profile_acc"],
            })
            continue

        # Gate 3: Score discrimination vs GPCR controls (target >= gpcr + 20.0 bits)
        if (target_score - score_gpcr) < 20.0:
            rejections.append({
                "protein_id": prot_id,
                "reason": f"Excluded_GPCR_cross_reactivity_diff_{target_score - score_gpcr:.1f}<20b",
                "target_score": target_score,
                "gpcr_score": score_gpcr,
                "hmm_coverage": cov,
                "dom_ievalue": ie,
                "profile_acc": best_target["profile_acc"],
            })
            continue

        assigned_class = "Type-1_Microbial_Rhodopsin" if score_type1 >= score_helio else "Heliorhodopsin_Type-3"
        candidates[prot_id] = {
            "assigned_class": assigned_class,
            "best_hit": best_target,
            "best_type1_score": score_type1,
            "best_helio_score": score_helio,
            "best_gpcr_score": score_gpcr,
        }

    return candidates, rejections


def stream_sequences_from_faa(faa_path, target_ids):
    seqs = {}
    with open(faa_path, "r", encoding="utf-8") as f:
        curr_id = None
        curr_hdr = ""
        curr_seq = []
        for line in f:
            if line.startswith(">"):
                if curr_id and curr_id in target_ids:
                    seqs[curr_id] = (curr_hdr, "".join(curr_seq))
                hdr = line.strip()[1:]
                curr_id = hdr.split()[0]
                curr_hdr = hdr
                curr_seq = []
            else:
                curr_seq.append(line.strip())
        if curr_id and curr_id in target_ids:
            seqs[curr_id] = (curr_hdr, "".join(curr_seq))
    return seqs


def stream_h3_coordinates(mapping_path, target_ids):
    coords = {}
    with open(mapping_path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for r in reader:
            cid = r.get("canonical_seq_id", "")
            if cid in target_ids:
                coords[cid] = {
                    "mag": r.get("mag_id") or r.get("genome_id", "Unknown"),
                    "contig": r.get("contig_id", "Unknown"),
                    "start": r.get("start_coord", "0"),
                    "stop": r.get("end_coord", "0"),
                    "strand": r.get("strand", "+"),
                    "partial": r.get("partial", "unknown"),
                }
    return coords


def run_diamond_search(candidates_faa, out_tsv):
    print("Running one-way DIAMOND blastp (--ultra-sensitive) against all_f345_references.dmnd...")
    cmd = [
        str(DIAMOND_BIN), "blastp",
        "--query", str(candidates_faa),
        "--db", str(ALL_REFS_DMND),
        "--out", str(out_tsv),
        "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "evalue", "bitscore", "qstart", "qend", "sstart", "send", "qlen", "slen",
        "--max-target-seqs", "5",
        "--evalue", "1e-5",
        "--ultra-sensitive",
        "--threads", "4",
        "--quiet",
    ]
    subprocess.run(cmd, check=True)

    hits_by_q = {}
    if out_tsv.exists():
        with open(out_tsv, "r", encoding="utf-8") as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 12:
                    qid = parts[0]
                    qstart, qend = int(parts[6]), int(parts[7])
                    sstart, send = int(parts[8]), int(parts[9])
                    qlen, slen = int(parts[10]), int(parts[11])
                    qcov = (abs(qend - qstart) + 1) / float(qlen) if qlen > 0 else 0.0
                    scov = (abs(send - sstart) + 1) / float(slen) if slen > 0 else 0.0

                    rec = {
                        "sseqid": parts[1],
                        "pident": float(parts[2]),
                        "length": int(parts[3]),
                        "evalue": float(parts[4]),
                        "bitscore": float(parts[5]),
                        "qcov": round(qcov, 4),
                        "scov": round(scov, 4),
                    }
                    hits_by_q.setdefault(qid, []).append(rec)
    return hits_by_q


def run_mmseqs_clustering(candidates_faa, prefix_path, tmp_dir):
    print("Running MMseqs2 easy-linclust at 95% sequence identity...")
    cmd = [
        str(MMSEQS_BIN), "easy-linclust",
        str(candidates_faa),
        str(prefix_path),
        str(tmp_dir),
        "--min-seq-id", "0.95",
        "-c", "0.80",
        "--cov-mode", "0",
        "--threads", "4",
    ]
    subprocess.run(cmd, check=True)

    clu_tsv = Path(f"{prefix_path}_cluster.tsv")
    assert clu_tsv.exists(), f"CRITICAL: Expected MMseqs2 cluster file {clu_tsv} not found!"
    clu_map = {}
    with open(clu_tsv, "r", encoding="utf-8") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                rep, member = parts[0], parts[1]
                clu_map[member] = rep
    assert len(clu_map) > 0, "MMseqs2 cluster mapping is empty!"
    print(f"MMseqs2 clustering complete: {len(set(clu_map.values()))} clusters across {len(clu_map)} members.")
    return clu_map


def main():
    start_time = time.time()
    print("=" * 80)
    print("TRACK F3: MICROBIAL & HELIORHODOPSIN FULL DISCOVERY CENSUS")
    print(f"Timestamp: {time.strftime('%Y-%m-%d %H:%M:%S UTC', time.gmtime())}")
    print("Cohorts: 5,438 EQ MAGs, SMAG environmental, GEM environmental")
    print("=" * 80)

    # Freeze executed code BEFORE execution
    FROZEN_DIR.mkdir(parents=True, exist_ok=True)
    shutil.copy2(__file__, FROZEN_DIR / "executed_08_f3_full.py")

    # 1. Assert required inputs
    assert F3_HMM.exists(), f"CRITICAL: Missing {F3_HMM}"
    assert EQ_COHORT_FAA.exists(), f"CRITICAL: Missing EQ cohort FASTA: {EQ_COHORT_FAA}"
    assert SMAG_CANONICAL.exists(), f"CRITICAL: Missing SMAG FASTA: {SMAG_CANONICAL}"
    assert GEM_CANONICAL.exists(), f"CRITICAL: Missing GEM FASTA: {GEM_CANONICAL}"

    mag_to_cluster, tax_map = load_metadata()

    cohort_configs = [
        {
            "name": "EQ_MAG",
            "faa": EQ_COHORT_FAA,
            "type": "eq",
        },
        {
            "name": "Environmental_SMAG",
            "faa": SMAG_CANONICAL,
            "type": "smag",
            "coord_map": SMAG_COORD_MAP,
        },
        {
            "name": "Environmental_GEM",
            "faa": GEM_CANONICAL,
            "type": "gem",
            "coord_map": GEM_COORD_MAP,
        },
    ]

    all_candidates = []
    all_rejections = []
    cpus = os.environ.get("SLURM_CPUS_PER_TASK", "16")

    for cfg in cohort_configs:
        cname = cfg["name"]
        faa = cfg["faa"]
        print(f"\n--- Processing Cohort: {cname} ({faa.name}) ---")

        domtbl_out = F3_DIR / f"f3_domtbl_{cname}.tsv"
        hmm_log = F3_DIR / f"f3_hmmsearch_{cname}.log"

        if domtbl_out.exists() and domtbl_out.stat().st_size > 0:
            print(f"  Reusing existing hmmsearch domtblout: {domtbl_out}")
        else:
            hmm_cmd = [
                "hmmsearch",
                "--cpu", cpus,
                "-E", "1e-3",
                "--domtblout", str(domtbl_out),
                "-o", str(hmm_log),
                str(F3_HMM),
                str(faa),
            ]
            t0 = time.time()
            subprocess.run(hmm_cmd, check=True)
            print(f"  hmmsearch complete in {time.time() - t0:.1f}s.")

        raw_hits = parse_domtblout(domtbl_out)
        print(f"  Raw domain hits: {len(raw_hits)}")

        hits_by_prot = {}
        for h in raw_hits:
            hits_by_prot.setdefault(h["target_name"], []).append(h)

        cands, rejs = evaluate_f3_gates(hits_by_prot)
        for r in rejs:
            r["cohort"] = cname
        all_rejections.extend(rejs)
        print(f"  Proteins passing gates: {len(cands)}")

        if not cands:
            continue

        target_ids = set(cands.keys())
        print(f"  Extracting sequences for {len(target_ids)} candidates from {faa.name}...")
        seq_data = stream_sequences_from_faa(faa, target_ids)

        coord_data = {}
        if cfg["type"] in ("smag", "gem"):
            print(f"  Streaming coordinates from {cfg['coord_map'].name}...")
            coord_data = stream_h3_coordinates(cfg["coord_map"], target_ids)

        for pid, cinfo in cands.items():
            if pid not in seq_data:
                all_rejections.append({
                    "protein_id": pid,
                    "cohort": cname,
                    "reason": "Missing_in_FASTA_sequence_stream",
                    "target_score": cinfo["best_hit"]["dom_score"],
                    "gpcr_score": cinfo["best_gpcr_score"],
                    "hmm_coverage": cinfo["best_hit"]["hmm_coverage"],
                    "dom_ievalue": cinfo["best_hit"]["dom_ievalue"],
                    "profile_acc": cinfo["best_hit"]["profile_acc"],
                })
                continue

            hdr, seq = seq_data[pid]
            clean_seq = seq.rstrip("*")
            seq_sha = hashlib.sha256(clean_seq.encode("utf-8")).hexdigest()
            bh = cinfo["best_hit"]

            # Completeness & coordinate determination
            if cfg["type"] == "eq":
                # Parse Prodigal partial flag from header
                pmatch = re.search(r"partial=(\d{2})", hdr)
                partial_flag = pmatch.group(1) if pmatch else "unknown"
                is_comp = 1 if partial_flag == "00" else 0

                fields = hdr.split(" # ")
                start_p = fields[1] if len(fields) > 1 else "0"
                stop_p = fields[2] if len(fields) > 2 else "0"
                strand_p = fields[3] if len(fields) > 3 else "+"

                mag = pid.split("::")[0] if "::" in pid else pid.split("_")[0]
                contig = pid.split("::")[1].rsplit("_", 1)[0] if "::" in pid else pid.rsplit("_", 1)[0]
                sec_clust = mag_to_cluster.get(mag, "Unknown")
                tax = tax_map.get(mag, "Unclassified")
            else:
                cd = coord_data.get(pid, {})
                partial_flag = cd.get("partial", "unknown")
                is_comp = 1 if partial_flag == "00" else 0
                start_p = cd.get("start", "0")
                stop_p = cd.get("stop", "0")
                strand_p = cd.get("strand", "+")
                mag = cd.get("mag", pid.split("::")[0] if "::" in pid else "Unknown")
                contig = cd.get("contig", "Unknown")
                sec_clust = "Environmental_Cohort"
                tax = f"{cname}_Environmental"

            all_candidates.append({
                "protein_id": pid,
                "cohort": cname,
                "mag_or_genome": mag,
                "contig": contig,
                "start": start_p,
                "stop": stop_p,
                "strand": strand_p,
                "secondary_cluster_ani95": sec_clust,
                "gtdbtk_taxonomy": tax,
                "assigned_class": cinfo["assigned_class"],
                "best_profile_name": bh["profile_name"],
                "profile_accession": bh["profile_acc"],
                "domain_score": bh["dom_score"],
                "domain_ievalue": bh["dom_ievalue"],
                "hmm_coverage": round(bh["hmm_coverage"], 4),
                "best_type1_score": cinfo["best_type1_score"],
                "best_helio_score": cinfo["best_helio_score"],
                "gpcr_control_score": cinfo["best_gpcr_score"],
                "prodigal_partial": partial_flag,
                "is_complete_gene": is_comp,
                "length_aa": len(clean_seq),
                "seq_sha256": seq_sha,
                "sequence": clean_seq,
            })

    print(f"\nTotal F3 candidates across all cohorts: {len(all_candidates)}")
    total_complete = sum(c["is_complete_gene"] == 1 for c in all_candidates)
    print(f"Total complete coding genes (partial=00): {total_complete}")

    # Emit combined candidate FASTA with bare protein_id headers
    candidates_faa = F3_DIR / "f3_all_candidates.faa"
    with open(candidates_faa, "w", encoding="utf-8") as f_out:
        for c in all_candidates:
            f_out.write(f">{c['protein_id']}\n{c['sequence']}\n")

    # One-way DIAMOND blastp
    diamond_tsv = F3_DIR / "f3_diamond_hits.tsv"
    if diamond_tsv.exists() and diamond_tsv.stat().st_size > 0:
        print(f"Reusing existing DIAMOND hits from {diamond_tsv}...")
        diamond_hits = {}
        with open(diamond_tsv, "r", encoding="utf-8") as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 12:
                    qid = parts[0]
                    qstart, qend = int(parts[6]), int(parts[7])
                    sstart, send = int(parts[8]), int(parts[9])
                    qlen, slen = int(parts[10]), int(parts[11])
                    qcov = (abs(qend - qstart) + 1) / float(qlen) if qlen > 0 else 0.0
                    scov = (abs(send - sstart) + 1) / float(slen) if slen > 0 else 0.0
                    rec = {
                        "sseqid": parts[1],
                        "pident": float(parts[2]),
                        "length": int(parts[3]),
                        "evalue": float(parts[4]),
                        "bitscore": float(parts[5]),
                        "qcov": round(qcov, 4),
                        "scov": round(scov, 4),
                    }
                    diamond_hits.setdefault(qid, []).append(rec)
    else:
        diamond_hits = run_diamond_search(candidates_faa, diamond_tsv) if all_candidates else {}

    for c in all_candidates:
        dh = diamond_hits.get(c["protein_id"], [])
        if dh:
            top = dh[0]
            c["top_ref_id"] = top["sseqid"]
            c["top_ref_pident"] = top["pident"]
            c["top_ref_evalue"] = top["evalue"]
            c["top_ref_bitscore"] = top["bitscore"]
            c["top_ref_qcov"] = top["qcov"]
            c["top_ref_scov"] = top["scov"]
            c["top_competing_refs"] = ";".join([f"{h['sseqid']}(id={h['pident']}%,bs={h['bitscore']})" for h in dh[1:4]])
        else:
            c["top_ref_id"] = "NO_HIT"
            c["top_ref_pident"] = 0.0
            c["top_ref_evalue"] = 1.0
            c["top_ref_bitscore"] = 0.0
            c["top_ref_qcov"] = 0.0
            c["top_ref_scov"] = 0.0
            c["top_competing_refs"] = "None"

    # MMseqs2 95% clustering
    if all_candidates:
        clu_prefix = F3_DIR / "f3_clu95"
        clu_tsv = Path(f"{clu_prefix}_cluster.tsv")
        if clu_tsv.exists() and clu_tsv.stat().st_size > 0:
            print(f"Reusing existing MMseqs2 cluster TSV: {clu_tsv}...")
            clu_map = {}
            with open(clu_tsv, "r", encoding="utf-8") as f:
                for line in f:
                    parts = line.strip().split("\t")
                    if len(parts) >= 2:
                        clu_map[parts[1]] = parts[0]
        else:
            tmp_mmseqs = F3_DIR / "tmp_mmseqs"
            tmp_mmseqs.mkdir(exist_ok=True)
            clu_map = run_mmseqs_clustering(candidates_faa, clu_prefix, tmp_mmseqs)
        for c in all_candidates:
            c["mmseqs_cluster_rep_95"] = clu_map.get(c["protein_id"], c["protein_id"])
    else:
        for c in all_candidates:
            c["mmseqs_cluster_rep_95"] = c["protein_id"]

    # Sequence deduplication & cross-source membership tracking
    seq_to_members = {}
    for c in all_candidates:
        sha = c["seq_sha256"]
        if sha not in seq_to_members:
            seq_to_members[sha] = {
                "sequence": c["sequence"],
                "length_aa": c["length_aa"],
                "assigned_class": c["assigned_class"],
                "cohorts": set(),
                "protein_ids": [],
                "mags": set(),
                "taxonomies": set(),
                "has_complete": False,
            }
        seq_to_members[sha]["cohorts"].add(c["cohort"])
        seq_to_members[sha]["protein_ids"].append(c["protein_id"])
        seq_to_members[sha]["mags"].add(c["mag_or_genome"])
        seq_to_members[sha]["taxonomies"].add(c["gtdbtk_taxonomy"])
        if c["is_complete_gene"] == 1:
            seq_to_members[sha]["has_complete"] = True

    for c in all_candidates:
        sm = seq_to_members[c["seq_sha256"]]
        c["shared_across_cohorts"] = ";".join(sorted(sm["cohorts"]))

    # Write Candidate Summary TSV
    summary_tsv = F3_DIR / "f3_candidates_summary.tsv"
    fieldnames = [
        "protein_id", "cohort", "mag_or_genome", "contig", "start", "stop", "strand",
        "secondary_cluster_ani95", "gtdbtk_taxonomy", "assigned_class",
        "best_profile_name", "profile_accession", "domain_score", "domain_ievalue", "hmm_coverage",
        "best_type1_score", "best_helio_score", "gpcr_control_score",
        "top_ref_id", "top_ref_pident", "top_ref_evalue", "top_ref_bitscore", "top_ref_qcov", "top_ref_scov",
        "top_competing_refs", "mmseqs_cluster_rep_95", "shared_across_cohorts",
        "prodigal_partial", "is_complete_gene", "length_aa", "seq_sha256"
    ]
    with open(summary_tsv, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(all_candidates)
    print(f"Saved candidates summary to {summary_tsv}")

    # Write Rejections TSV
    rejections_tsv = F3_DIR / "f3_rejections_and_ambiguous.tsv"
    rej_fields = ["protein_id", "cohort", "reason", "target_score", "gpcr_score", "hmm_coverage", "dom_ievalue", "profile_acc"]
    with open(rejections_tsv, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=rej_fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(all_rejections)
    print(f"Saved rejections and ambiguous cases to {rejections_tsv}")

    # Prepare ProstT5 Inputs: separate by actual family and include curated references
    print("\nPreparing ProstT5 input packages with curated references...")
    ref_records = []
    if ALL_REFS_FAA.exists() and REFS_META_TSV.exists():
        ref_seqs = {}
        curr_acc = None
        curr_seq = []
        with open(ALL_REFS_FAA, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    if curr_acc:
                        ref_seqs[curr_acc] = "".join(curr_seq)
                    hdr = line[1:].split()[0]
                    curr_acc = hdr.split("|")[0]
                    curr_seq = []
                else:
                    curr_seq.append(line)
            if curr_acc:
                ref_seqs[curr_acc] = "".join(curr_seq)

        with open(REFS_META_TSV, "r", encoding="utf-8") as f:
            for r in csv.DictReader(f, delimiter="\t"):
                if r.get("track") == "F3":
                    acc = r["accession"]
                    if acc in ref_seqs:
                        r["sequence"] = ref_seqs[acc]
                        ref_records.append(r)

    families = [
        ("Type-1_Microbial_Rhodopsin", "prostt5_type-1_microbial_rhodopsin", "Type1_Microbial"),
        ("Heliorhodopsin_Type-3", "prostt5_heliorhodopsin_type-3", "Heliorhodopsin"),
    ]
    for fam, file_prefix, cat_match in families:
        fam_faa = PROSTT5_DIR / f"{file_prefix}.faa"
        fam_manifest_tsv = PROSTT5_DIR / f"{file_prefix}_manifest.tsv"

        fam_seqs = {sha: data for sha, data in seq_to_members.items() if data["assigned_class"] == fam}
        fam_refs = [r for r in ref_records if r.get("category") == cat_match]
        assert len(fam_refs) > 0, f"CRITICAL: 0 curated references matched for {fam}"

        complete_fam_seqs = {sha: data for sha, data in fam_seqs.items() if data["has_complete"]}
        partial_fam_seqs = {sha: data for sha, data in fam_seqs.items() if not data["has_complete"]}

        # 1. PRIMARY ProstT5 Package: Curated references + Complete (partial=00) candidates
        manifest_rows = []
        with open(fam_faa, "w", encoding="utf-8") as out_f:
            # 1a. Curated References
            for ref in fam_refs:
                rid = ref["accession"]
                seq = ref["sequence"]
                r_sha = hashlib.sha256(seq.encode("utf-8")).hexdigest()
                out_f.write(f">{rid}\n{seq}\n")
                manifest_rows.append({
                    "seq_id": rid,
                    "seq_sha256": r_sha,
                    "length_aa": len(seq),
                    "source_type": "Curated_Reference",
                    "member_count": 1,
                    "member_ids": rid,
                    "carrier_mag_or_genomes": ref.get("organism", "Reference"),
                    "taxonomies": ref.get("lineage", "Reference"),
                    "assigned_family": fam,
                    "is_complete_gene": 1,
                })

            # 1b. Complete Candidates (partial=00)
            for sha, sdata in complete_fam_seqs.items():
                rep_id = sdata["protein_ids"][0]
                out_f.write(f">{rep_id}\n{sdata['sequence']}\n")
                src_type = ";".join(sorted(sdata["cohorts"]))
                manifest_rows.append({
                    "seq_id": rep_id,
                    "seq_sha256": sha,
                    "length_aa": sdata["length_aa"],
                    "source_type": src_type,
                    "member_count": len(sdata["protein_ids"]),
                    "member_ids": ";".join(sdata["protein_ids"]),
                    "carrier_mag_or_genomes": ";".join(sorted(sdata["mags"])),
                    "taxonomies": ";".join(sorted(sdata["taxonomies"])),
                    "assigned_family": fam,
                    "is_complete_gene": 1,
                })

        m_fields = [
            "seq_id", "seq_sha256", "length_aa", "source_type", "member_count",
            "member_ids", "carrier_mag_or_genomes", "taxonomies", "assigned_family", "is_complete_gene"
        ]
        with open(fam_manifest_tsv, "w", newline="", encoding="utf-8") as mf:
            w = csv.DictWriter(mf, fieldnames=m_fields, delimiter="\t")
            w.writeheader()
            w.writerows(manifest_rows)

        # Validate PRIMARY FASTA ID closure
        fasta_ids = []
        with open(fam_faa, "r", encoding="utf-8") as f:
            for line in f:
                if line.startswith(">"):
                    fasta_ids.append(line[1:].strip().split()[0])
        manifest_ids = [m["seq_id"] for m in manifest_rows]
        assert fasta_ids == manifest_ids, f"ID closure failure for primary {fam}: {len(fasta_ids)} in FASTA vs {len(manifest_ids)} in manifest"
        assert len(set(manifest_ids)) == len(manifest_ids), f"Duplicate seq_ids in primary manifest for {fam}"

        # 2. SEPARATE Partial Candidates Package (retained, NOT silently dropped)
        part_faa = PROSTT5_DIR / f"{file_prefix}_partial_excluded.faa"
        part_manifest_tsv = PROSTT5_DIR / f"{file_prefix}_partial_excluded_manifest.tsv"
        part_rows = []
        with open(part_faa, "w", encoding="utf-8") as out_pf:
            for sha, sdata in partial_fam_seqs.items():
                rep_id = sdata["protein_ids"][0]
                out_pf.write(f">{rep_id}\n{sdata['sequence']}\n")
                src_type = ";".join(sorted(sdata["cohorts"]))
                part_rows.append({
                    "seq_id": rep_id,
                    "seq_sha256": sha,
                    "length_aa": sdata["length_aa"],
                    "source_type": src_type,
                    "member_count": len(sdata["protein_ids"]),
                    "member_ids": ";".join(sdata["protein_ids"]),
                    "carrier_mag_or_genomes": ";".join(sorted(sdata["mags"])),
                    "taxonomies": ";".join(sorted(sdata["taxonomies"])),
                    "assigned_family": fam,
                    "is_complete_gene": 0,
                    "exclusion_reason": "Incomplete_gene_fragment_partial_not_00",
                })
        p_fields = m_fields + ["exclusion_reason"]
        with open(part_manifest_tsv, "w", newline="", encoding="utf-8") as pf:
            w = csv.DictWriter(pf, fieldnames=p_fields, delimiter="\t")
            w.writeheader()
            w.writerows(part_rows)

        part_fasta_ids = []
        with open(part_faa, "r", encoding="utf-8") as f:
            for line in f:
                if line.startswith(">"):
                    part_fasta_ids.append(line[1:].strip().split()[0])
        part_manifest_ids = [m["seq_id"] for m in part_rows]
        assert part_fasta_ids == part_manifest_ids, f"ID closure failure for partial {fam}"

        print(f"Validated PRIMARY ProstT5 package for {fam}: {len(manifest_rows)} sequences ({len(fam_refs)} refs + {len(complete_fam_seqs)} complete candidates) to {fam_faa.name}")
        print(f"Retained {len(part_rows)} partial candidates in separate table {part_manifest_tsv.name}")

    # Step 7: Emit Metrics JSON
    metrics = {
        "track": "F3",
        "family": "Microbial_and_Heliorhodopsin",
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "total_candidates": len(all_candidates),
        "total_complete_candidates": total_complete,
        "by_cohort": {
            "EQ_MAG": sum(1 for c in all_candidates if c["cohort"] == "EQ_MAG"),
            "Environmental_SMAG": sum(1 for c in all_candidates if c["cohort"] == "Environmental_SMAG"),
            "Environmental_GEM": sum(1 for c in all_candidates if c["cohort"] == "Environmental_GEM"),
        },
        "by_class": {
            "Type-1_Microbial_Rhodopsin": sum(1 for c in all_candidates if c["assigned_class"] == "Type-1_Microbial_Rhodopsin"),
            "Heliorhodopsin_Type-3": sum(1 for c in all_candidates if c["assigned_class"] == "Heliorhodopsin_Type-3"),
        },
        "complete_by_class": {
            "Type-1_Microbial_Rhodopsin": sum(1 for c in all_candidates if c["assigned_class"] == "Type-1_Microbial_Rhodopsin" and c["is_complete_gene"] == 1),
            "Heliorhodopsin_Type-3": sum(1 for c in all_candidates if c["assigned_class"] == "Heliorhodopsin_Type-3" and c["is_complete_gene"] == 1),
        },
        "unique_sequences": len(seq_to_members),
        "total_rejections_and_ambiguous": len(all_rejections),
        "execution_time_seconds": round(time.time() - start_time, 2),
    }
    metrics_json = F3_DIR / "f3_full_metrics.json"
    with open(metrics_json, "w", encoding="utf-8") as f:
        json.dump(metrics, f, indent=2)
    print(f"\nSaved metrics to {metrics_json}")

    # Freeze code
    shutil.copy2(Path(__file__), FROZEN_DIR / "executed_08_f3_full.py")
    print(f"Code frozen to {FROZEN_DIR / 'executed_08_f3_full.py'}")
    print("=" * 80)
    print("TRACK F3 FULL DISCOVERY COMPLETE.")
    print("=" * 80)


if __name__ == "__main__":
    main()
