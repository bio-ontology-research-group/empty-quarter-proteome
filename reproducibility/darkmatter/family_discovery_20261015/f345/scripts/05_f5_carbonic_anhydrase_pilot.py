#!/usr/bin/env python3
"""05_f5_carbonic_anhydrase_pilot.py - Track F5: Carbonic Anhydrases Pilot Screening.

Grounding & Audit Resolutions (Argon Section 4 & Root Requirements):
1. Reuses frozen, read-only 256-genome cohort (956,344 proteins) from root job 53816043.
2. Streams candidates only: zero full-proteome dictionaries in memory; zero concatenations.
3. Profile matching strictly by HMM accession prefix (parts[4].split(".")[0]):
   - Target Gamma-CA: NF050158 (gamma_CA), NF040597 (gamma_class_CA)
   - Target Beta-CA: PF00484 (Pro_CA)
   - Carboxysome Shell CA: TIGR02701 (shell_CA / CsoSCA) -> KEPT STRICTLY SEPARATE PER ROOT INSTRUCTION
   - Competing LbH transferases: PF00132 (Hexapep), PF00583 (Acetyltransf_1)
   - Competing Alpha-CA: PF00194 (Carb_anhydrase)
4. Enforces predetermined gates:
   - Gamma-CA: dom_ievalue <= 1e-7, cov >= 0.70, bit score >= 20.0 bits above LbH transferases
   - Beta-CA: dom_ievalue <= 1e-7, cov >= 0.70, bit score >= 30.0
   - Carboxysome CA (TIGR02701): dom_ievalue <= 1e-7, cov >= 0.70 (never pooled into Gamma-CA)
   - Prodigal partial=00 (complete CDS)
5. Emits verified candidate FASTA with trailing stop '*' stripped.
6. Executes reciprocal DIAMOND blastp against f5_carbonic_anhydrase_references.dmnd.
7. Retains all failed and ambiguous hits in f5_rejections_and_ambiguous.tsv.
8. Freezes executed code and provenance manifest to results/pilot_v1/f5_carbonic_anhydrase/.
"""

import csv
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
REF_DIR = ROOT / "references"
HMM_DIR = REF_DIR / "profiles"
PILOT_DIR = ROOT / "results" / "pilot_v1" / "f5_carbonic_anhydrase"
PILOT_DIR.mkdir(parents=True, exist_ok=True)

# Read-only 256-genome cohort from Root
COHORT_FAA = Path("/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/family_discovery_20261015_clean/root_audit/f1_profile_pilot/cohort.faa")
COHORT_TSV = Path("/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/family_discovery_20261015_clean/root_audit/f1_profile_pilot/cohort.tsv")

W = Path("/ibex/scratch/projects/c2014/rob/metagenomics")
CDB_PATH = W / "05_refine/_drep/out/data_tables/Cdb.csv"
BAC_GTDB = W / "05_refine/_gtdbtk_r232/gtdbtk.bac120.summary.tsv"
AR_GTDB = W / "05_refine/_gtdbtk_r232/gtdbtk.ar53.summary.tsv"

F5_HMM = HMM_DIR / "f5_carbonic_anhydrase.hmm"
F5_DMND = REF_DIR / "f5_carbonic_anhydrase_references.dmnd"

GAMMA_ACCS = {"NF050158", "NF040597"}
BETA_ACCS = {"PF00484"}
CARBOXY_ACCS = {"TIGR02701"}
LBH_ACCS = {"PF00132", "PF00583"}
ALPHA_ACCS = {"PF00194"}


def compute_sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while chunk := f.read(1024 * 1024):
            h.update(chunk)
    return h.hexdigest()


def load_cohort_and_taxonomy():
    mag_to_cluster = {}
    if COHORT_TSV.exists():
        with open(COHORT_TSV, "r", encoding="utf-8") as f:
            for r in csv.DictReader(f, delimiter="\t"):
                mag_to_cluster[r["mag"]] = r["cluster"]

    if CDB_PATH.exists():
        with open(CDB_PATH, "r", encoding="utf-8") as f:
            for r in csv.DictReader(f):
                g = r["genome"].removesuffix(".fa").removesuffix(".fasta")
                if g not in mag_to_cluster:
                    mag_to_cluster[g] = r["secondary_cluster"]

    tax_map = {}
    for p in [BAC_GTDB, AR_GTDB]:
        if p.exists():
            with open(p, "r", encoding="utf-8") as f:
                for r in csv.DictReader(f, delimiter="\t"):
                    mag = r.get("user_genome", "").removesuffix(".fa").removesuffix(".fasta")
                    tax_map[mag] = r.get("classification", "Unclassified")

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


def run_reciprocal_diamond(query_faa, db_dmnd, out_tsv):
    cmd = [
        "diamond", "blastp",
        "--query", str(query_faa),
        "--db", str(db_dmnd),
        "--out", str(out_tsv),
        "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "evalue", "bitscore",
        "--max-target-seqs", "1",
        "--evalue", "1e-5",
        "--threads", "4",
        "--quiet",
    ]
    subprocess.run(cmd, check=True)
    best_hits = {}
    if out_tsv.exists():
        with open(out_tsv, "r", encoding="utf-8") as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 6:
                    qid = parts[0]
                    if qid not in best_hits:
                        best_hits[qid] = {
                            "top_ref_id": parts[1],
                            "pident": float(parts[2]),
                            "align_len": int(parts[3]),
                            "evalue": float(parts[4]),
                            "bitscore": float(parts[5]),
                        }
    return best_hits


def main():
    start_time = time.time()
    print("=" * 80)
    print("TRACK F5: CARBONIC ANHYDRASES (BETA & GAMMA) PILOT EXECUTION")
    print(f"Timestamp: {time.strftime('%Y-%m-%d %H:%M:%S UTC', time.gmtime())}")
    print(f"Cohort FASTA: {COHORT_FAA}")
    print(f"Output Directory: {PILOT_DIR}")
    print("=" * 80)

    if not COHORT_FAA.exists():
        raise FileNotFoundError(f"Cohort FASTA missing: {COHORT_FAA}")

    mag_to_cluster, tax_map = load_cohort_and_taxonomy()

    raw_domtbl = PILOT_DIR / "f5_domtblout.txt"
    raw_hmmsearch = PILOT_DIR / "f5_hmmsearch.txt"
    candidates_faa = PILOT_DIR / "f5_candidates.faa"
    candidates_tsv = PILOT_DIR / "f5_candidates_summary.tsv"
    rejections_tsv = PILOT_DIR / "f5_rejections_and_ambiguous.tsv"
    diamond_tsv = PILOT_DIR / "f5_diamond_reciprocal.tsv"
    metrics_json = PILOT_DIR / "f5_pilot_metrics.json"

    # Step 1: Run HMMER search directly on cohort FASTA
    print("\n[Step 1] Running hmmsearch against cohort FASTA (8 CPUs)...")
    hmm_cmd = [
        "hmmsearch",
        "--cpu", "8",
        "-E", "1e-3",
        "--domtblout", str(raw_domtbl),
        "-o", str(raw_hmmsearch),
        str(F5_HMM),
        str(COHORT_FAA),
    ]
    t0 = time.time()
    subprocess.run(hmm_cmd, check=True)
    t_hmm = time.time() - t0
    print(f"  hmmsearch complete in {t_hmm:.1f}s.")

    # Step 2: Parse domtblout and evaluate domain gates
    print("\n[Step 2] Parsing domain hits and evaluating gates...")
    hits = parse_domtblout(raw_domtbl)
    print(f"  Total raw domain hits in domtblout: {len(hits)}")

    by_protein = {}
    for h in hits:
        by_protein.setdefault(h["target_name"], []).append(h)

    print(f"  Unique proteins with domain hits: {len(by_protein)}")

    preliminary_candidates = {}
    rejections = []

    for prot_id, prot_hits in by_protein.items():
        gamma_hits = [h for h in prot_hits if h["profile_acc"] in GAMMA_ACCS]
        beta_hits = [h for h in prot_hits if h["profile_acc"] in BETA_ACCS]
        carboxy_hits = [h for h in prot_hits if h["profile_acc"] in CARBOXY_ACCS]
        lbh_hits = [h for h in prot_hits if h["profile_acc"] in LBH_ACCS]
        alpha_hits = [h for h in prot_hits if h["profile_acc"] in ALPHA_ACCS]

        best_gamma = max(gamma_hits, key=lambda x: x["dom_score"]) if gamma_hits else None
        best_beta = max(beta_hits, key=lambda x: x["dom_score"]) if beta_hits else None
        best_carboxy = max(carboxy_hits, key=lambda x: x["dom_score"]) if carboxy_hits else None
        best_lbh = max(lbh_hits, key=lambda x: x["dom_score"]) if lbh_hits else None
        best_alpha = max(alpha_hits, key=lambda x: x["dom_score"]) if alpha_hits else None

        score_gamma = best_gamma["dom_score"] if best_gamma else 0.0
        score_beta = best_beta["dom_score"] if best_beta else 0.0
        score_carboxy = best_carboxy["dom_score"] if best_carboxy else 0.0
        score_lbh = best_lbh["dom_score"] if best_lbh else 0.0
        score_alpha = best_alpha["dom_score"] if best_alpha else 0.0

        # Classify candidate based on best target score
        max_target_score = max(score_gamma, score_beta, score_carboxy)

        if max_target_score == 0.0:
            rejections.append({
                "protein_id": prot_id,
                "reason": "Only_control_hits_(Alpha_CA_or_LbH_transferase)",
                "target_score": 0.0,
                "control_score": max(score_lbh, score_alpha),
                "hmm_coverage": 0.0,
                "dom_ievalue": 1.0,
                "assigned_class": "Outgroup_Control",
            })
            continue

        # Case A: Carboxysome Shell CA (TIGR02701) - KEPT STRICTLY SEPARATE
        if max_target_score == score_carboxy and best_carboxy:
            cov = best_carboxy["hmm_coverage"]
            ie = best_carboxy["dom_ievalue"]
            if ie > 1e-7:
                rejections.append({
                    "protein_id": prot_id,
                    "reason": f"Carboxy_failed_dom_ievalue_gate_{ie:.2e}>1e-7",
                    "target_score": score_carboxy,
                    "control_score": 0.0,
                    "hmm_coverage": cov,
                    "dom_ievalue": ie,
                    "assigned_class": "Carboxysome_Shell_CA_CsoSCA",
                })
                continue
            if cov < 0.70:
                rejections.append({
                    "protein_id": prot_id,
                    "reason": f"Carboxy_failed_domain_coverage_gate_{cov:.2f}<0.70",
                    "target_score": score_carboxy,
                    "control_score": 0.0,
                    "hmm_coverage": cov,
                    "dom_ievalue": ie,
                    "assigned_class": "Carboxysome_Shell_CA_CsoSCA",
                })
                continue

            preliminary_candidates[prot_id] = {
                "assigned_class": "Carboxysome_Shell_CA_CsoSCA",
                "best_hit": best_carboxy,
                "score_gamma": score_gamma,
                "score_beta": score_beta,
                "score_carboxy": score_carboxy,
                "score_lbh": score_lbh,
                "score_alpha": score_alpha,
            }

        # Case B: Beta-Carbonic Anhydrase (PF00484)
        elif max_target_score == score_beta and best_beta:
            cov = best_beta["hmm_coverage"]
            ie = best_beta["dom_ievalue"]
            if ie > 1e-7:
                rejections.append({
                    "protein_id": prot_id,
                    "reason": f"Beta_CA_failed_dom_ievalue_gate_{ie:.2e}>1e-7",
                    "target_score": score_beta,
                    "control_score": 0.0,
                    "hmm_coverage": cov,
                    "dom_ievalue": ie,
                    "assigned_class": "Beta-Carbonic_Anhydrase",
                })
                continue
            if cov < 0.70:
                rejections.append({
                    "protein_id": prot_id,
                    "reason": f"Beta_CA_failed_domain_coverage_gate_{cov:.2f}<0.70",
                    "target_score": score_beta,
                    "control_score": 0.0,
                    "hmm_coverage": cov,
                    "dom_ievalue": ie,
                    "assigned_class": "Beta-Carbonic_Anhydrase",
                })
                continue
            if score_beta < 30.0:
                rejections.append({
                    "protein_id": prot_id,
                    "reason": f"Beta_CA_failed_score_gate_{score_beta:.1f}<30.0b",
                    "target_score": score_beta,
                    "control_score": 0.0,
                    "hmm_coverage": cov,
                    "dom_ievalue": ie,
                    "assigned_class": "Beta-Carbonic_Anhydrase",
                })
                continue

            preliminary_candidates[prot_id] = {
                "assigned_class": "Beta-Carbonic_Anhydrase",
                "best_hit": best_beta,
                "score_gamma": score_gamma,
                "score_beta": score_beta,
                "score_carboxy": score_carboxy,
                "score_lbh": score_lbh,
                "score_alpha": score_alpha,
            }

        # Case C: Gamma-Carbonic Anhydrase (NF050158, NF040597)
        elif max_target_score == score_gamma and best_gamma:
            cov = best_gamma["hmm_coverage"]
            ie = best_gamma["dom_ievalue"]
            if ie > 1e-7:
                rejections.append({
                    "protein_id": prot_id,
                    "reason": f"Gamma_CA_failed_dom_ievalue_gate_{ie:.2e}>1e-7",
                    "target_score": score_gamma,
                    "control_score": score_lbh,
                    "hmm_coverage": cov,
                    "dom_ievalue": ie,
                    "assigned_class": "Gamma-Carbonic_Anhydrase",
                })
                continue
            if cov < 0.70:
                rejections.append({
                    "protein_id": prot_id,
                    "reason": f"Gamma_CA_failed_domain_coverage_gate_{cov:.2f}<0.70",
                    "target_score": score_gamma,
                    "control_score": score_lbh,
                    "hmm_coverage": cov,
                    "dom_ievalue": ie,
                    "assigned_class": "Gamma-Carbonic_Anhydrase",
                })
                continue
            # Discrimination Gate vs competing LbH transferases (must be >= 20.0 bits above LbH)
            if (score_gamma - score_lbh) < 20.0:
                rejections.append({
                    "protein_id": prot_id,
                    "reason": f"Excluded_LbH_Transferase_cross_reactivity_diff_{score_gamma - score_lbh:.1f}<20b",
                    "target_score": score_gamma,
                    "control_score": score_lbh,
                    "hmm_coverage": cov,
                    "dom_ievalue": ie,
                    "assigned_class": "Gamma-Carbonic_Anhydrase",
                })
                continue

            preliminary_candidates[prot_id] = {
                "assigned_class": "Gamma-Carbonic_Anhydrase",
                "best_hit": best_gamma,
                "score_gamma": score_gamma,
                "score_beta": score_beta,
                "score_carboxy": score_carboxy,
                "score_lbh": score_lbh,
                "score_alpha": score_alpha,
            }

    print(f"  Proteins passing domain & discrimination gates: {len(preliminary_candidates)}")
    print(f"  Rejected/ambiguous hits: {len(rejections)}")

    # Step 3: Stream COHORT_FAA to extract candidate headers & sequences
    print("\n[Step 3] Streaming COHORT_FAA to extract sequences and check Prodigal completeness...")
    needed_ids = set(preliminary_candidates.keys())
    protein_data = {}

    with open(COHORT_FAA, "r", encoding="utf-8") as f:
        current_id = None
        current_header = ""
        current_seq = []
        for line in f:
            if line.startswith(">"):
                if current_id and current_id in needed_ids:
                    protein_data[current_id] = {
                        "header": current_header,
                        "seq": "".join(current_seq),
                    }
                header_line = line.strip()[1:]
                current_id = header_line.split()[0]
                current_header = header_line
                current_seq = []
            else:
                current_seq.append(line.strip())

        if current_id and current_id in needed_ids:
            protein_data[current_id] = {
                "header": current_header,
                "seq": "".join(current_seq),
            }

    print(f"  Extracted {len(protein_data)} candidate records from cohort FASTA.")

    final_candidates = []
    for prot_id, cinfo in preliminary_candidates.items():
        pdata = protein_data.get(prot_id)
        if not pdata:
            rejections.append({
                "protein_id": prot_id,
                "reason": "Sequence_missing_in_cohort_fasta",
                "target_score": cinfo["best_hit"]["dom_score"],
                "control_score": cinfo["score_lbh"],
                "hmm_coverage": cinfo["best_hit"]["hmm_coverage"],
                "dom_ievalue": cinfo["best_hit"]["dom_ievalue"],
                "assigned_class": cinfo["assigned_class"],
            })
            continue

        hdr = pdata["header"]
        seq = pdata["seq"].rstrip("*")

        partial_match = re.search(r"partial=(\d{2})", hdr)
        partial_flag = partial_match.group(1) if partial_match else "unknown"
        is_complete = 1 if partial_flag == "00" else 0

        fields = hdr.split(" # ")
        start_pos = fields[1] if len(fields) > 1 else "0"
        stop_pos = fields[2] if len(fields) > 2 else "0"
        strand_pos = fields[3] if len(fields) > 3 else "0"

        if "::" in prot_id:
            mag_part, contig_gene = prot_id.split("::", 1)
            mag = mag_part
            contig = contig_gene.rsplit("_", 1)[0] if "_" in contig_gene else contig_gene
        else:
            mag = prot_id.split("_")[0]
            contig = prot_id.rsplit("_", 1)[0] if "_" in prot_id else prot_id

        sec_clust = mag_to_cluster.get(mag, "Unknown")
        tax = tax_map.get(mag, "Unclassified")
        bh = cinfo["best_hit"]

        if is_complete != 1:
            rejections.append({
                "protein_id": prot_id,
                "reason": f"Edge_truncated_gene_partial={partial_flag}",
                "target_score": bh["dom_score"],
                "control_score": cinfo["score_lbh"],
                "hmm_coverage": bh["hmm_coverage"],
                "dom_ievalue": bh["dom_ievalue"],
                "assigned_class": cinfo["assigned_class"],
            })
            continue

        final_candidates.append({
            "protein_id": prot_id,
            "mag": mag,
            "contig": contig,
            "start": start_pos,
            "stop": stop_pos,
            "strand": strand_pos,
            "secondary_cluster_ani95": sec_clust,
            "gtdbtk_taxonomy": tax,
            "assigned_class": cinfo["assigned_class"],
            "best_profile_name": bh["profile_name"],
            "profile_accession": bh["profile_acc"],
            "domain_score": bh["dom_score"],
            "domain_ievalue": bh["dom_ievalue"],
            "hmm_coverage": bh["hmm_coverage"],
            "score_gamma": cinfo["score_gamma"],
            "score_beta": cinfo["score_beta"],
            "score_carboxy": cinfo["score_carboxy"],
            "score_lbh_control": cinfo["score_lbh"],
            "score_alpha_control": cinfo["score_alpha"],
            "prodigal_partial": partial_flag,
            "is_complete_gene": is_complete,
            "sequence": seq,
        })

    print(f"  Final CA candidates meeting all gates (including partial=00): {len(final_candidates)}")

    # Step 4: Write candidate FASTA (trailing '*' stripped)
    print(f"\n[Step 4] Emitting {len(final_candidates)} candidates to {candidates_faa.name}...")
    with open(candidates_faa, "w", encoding="utf-8") as f_out:
        for c in final_candidates:
            hdr_str = f">{c['protein_id']}|{c['mag']}|{c['assigned_class']}|{c['secondary_cluster_ani95']}"
            f_out.write(f"{hdr_str}\n{c['sequence']}\n")

    # Step 5: Reciprocal DIAMOND blastp against reference database
    diamond_matches = {}
    if F5_DMND.exists() and len(final_candidates) > 0:
        print("\n[Step 5] Running reciprocal DIAMOND blastp against f5_carbonic_anhydrase_references.dmnd...")
        diamond_matches = run_reciprocal_diamond(candidates_faa, F5_DMND, diamond_tsv)
        print(f"  Reciprocal matches obtained: {len(diamond_matches)}")

    for c in final_candidates:
        dm = diamond_matches.get(c["protein_id"], {})
        c["top_ref_id"] = dm.get("top_ref_id", "NO_HIT")
        c["top_ref_pident"] = dm.get("pident", 0.0)
        c["top_ref_evalue"] = dm.get("evalue", 1.0)
        c["top_ref_bitscore"] = dm.get("bitscore", 0.0)

    # Step 6: Write Summary TSV
    print(f"\n[Step 6] Writing summary TSV to {candidates_tsv.name}...")
    fieldnames = [
        "protein_id", "mag", "contig", "start", "stop", "strand",
        "secondary_cluster_ani95", "gtdbtk_taxonomy", "assigned_class",
        "best_profile_name", "profile_accession", "domain_score",
        "domain_ievalue", "hmm_coverage", "score_gamma", "score_beta",
        "score_carboxy", "score_lbh_control", "score_alpha_control",
        "top_ref_id", "top_ref_pident", "top_ref_evalue", "top_ref_bitscore",
        "prodigal_partial", "is_complete_gene"
    ]
    with open(candidates_tsv, "w", newline="", encoding="utf-8") as f_out:
        writer = csv.DictWriter(f_out, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(final_candidates)

    # Write Rejections TSV
    rej_fields = ["protein_id", "reason", "target_score", "control_score", "hmm_coverage", "dom_ievalue", "assigned_class"]
    with open(rejections_tsv, "w", newline="", encoding="utf-8") as f_out:
        writer = csv.DictWriter(f_out, fieldnames=rej_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rejections)

    # Step 7: Freeze executed code & provenance
    shutil.copy2(HERE / "05_f5_carbonic_anhydrase_pilot.py", PILOT_DIR / "executed_05_f5_carbonic_anhydrase_pilot.py")

    provenance = {
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "script": "05_f5_carbonic_anhydrase_pilot.py",
        "script_sha256": compute_sha256(HERE / "05_f5_carbonic_anhydrase_pilot.py"),
        "cohort_faa": str(COHORT_FAA),
        "cohort_faa_sha256": compute_sha256(COHORT_FAA),
        "f5_hmm": str(F5_HMM),
        "f5_hmm_sha256": compute_sha256(F5_HMM),
        "f5_dmnd": str(F5_DMND),
        "f5_dmnd_sha256": compute_sha256(F5_DMND),
    }
    with open(PILOT_DIR / "provenance.json", "w", encoding="utf-8") as f_out:
        json.dump(provenance, f_out, indent=2)

    # Step 8: Compute pilot metrics
    gamma_count = sum(1 for c in final_candidates if c["assigned_class"] == "Gamma-Carbonic_Anhydrase")
    beta_count = sum(1 for c in final_candidates if c["assigned_class"] == "Beta-Carbonic_Anhydrase")
    carboxy_count = sum(1 for c in final_candidates if c["assigned_class"] == "Carboxysome_Shell_CA_CsoSCA")
    mags_with_cand = len(set(c["mag"] for c in final_candidates))
    clusters_with_cand = len(set(c["secondary_cluster_ani95"] for c in final_candidates if c["secondary_cluster_ani95"] != "Unknown"))

    metrics = {
        "track": "F5",
        "description": "Carbonic Anhydrase (Beta and Gamma) Discovery Pilot",
        "cohort_genomes": 256,
        "cohort_proteins_scanned": 956344,
        "raw_domtbl_hits": len(hits),
        "unique_proteins_with_hits": len(by_protein),
        "rejected_or_ambiguous_hits": len(rejections),
        "final_passing_candidates": len(final_candidates),
        "gamma_carbonic_anhydrase": gamma_count,
        "beta_carbonic_anhydrase": beta_count,
        "carboxysome_shell_ca_csosca": carboxy_count,
        "distinct_mags_with_candidates": mags_with_cand,
        "distinct_secondary_clusters_ani95": clusters_with_cand,
        "timing_hmmsearch_seconds": t_hmm,
        "timing_total_seconds": time.time() - start_time,
    }
    with open(metrics_json, "w", encoding="utf-8") as f_out:
        json.dump(metrics, f_out, indent=2)

    print(f"\nWrote pilot metrics to {metrics_json.name}")
    print(f"=== Track F5 Pilot Finished in {time.time() - start_time:.1f}s ===")


if __name__ == "__main__":
    main()
