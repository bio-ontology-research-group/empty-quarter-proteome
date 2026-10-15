#!/usr/bin/env python3
"""03_f3_rhodopsin_pilot.py - Track F3: Rhodopsin & Heliorhodopsin Pilot Screening.

Grounding & Audit Resolutions (Argon Section 4 & Root Requirements):
1. Reuses frozen, read-only 256-genome cohort (956,344 proteins) from root job 53816043.
2. Streams candidates only: zero full-proteome dictionaries in memory; zero concatenations.
3. Profile matching strictly by HMM accession prefix (parts[4].split(".")[0]):
   - Target Type-1 Rhodopsin: PF01036 (Bac_rhodopsin)
   - Target Heliorhodopsin: PF18761 (Heliorhodopsin), NF038020 (Heliorhodopsin_HeR)
   - Competing GPCR controls: PF00001 (7tm_1), PF00002 (7tm_2)
4. Enforces predetermined gates:
   - dom_ievalue <= 1e-5
   - domain coverage >= 0.70
   - target bit score >= 20.0 bits above competing GPCR score
   - Prodigal partial=00 (complete CDS)
5. Emits verified candidate FASTA with trailing stop '*' stripped.
6. Executes reciprocal DIAMOND blastp against f3_rhodopsin_references.dmnd.
7. Retains all failed and ambiguous hits in f3_rejections_and_ambiguous.tsv.
8. Freezes executed code and provenance manifest to results/pilot_v1/f3_rhodopsin/.
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
PILOT_DIR = ROOT / "results" / "pilot_v1" / "f3_rhodopsin"
PILOT_DIR.mkdir(parents=True, exist_ok=True)

# Read-only 256-genome cohort from Root
COHORT_FAA = Path("/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/family_discovery_20261015_clean/root_audit/f1_profile_pilot/cohort.faa")
COHORT_TSV = Path("/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/family_discovery_20261015_clean/root_audit/f1_profile_pilot/cohort.tsv")

W = Path("/ibex/scratch/projects/c2014/rob/metagenomics")
CDB_PATH = W / "05_refine/_drep/out/data_tables/Cdb.csv"
BAC_GTDB = W / "05_refine/_gtdbtk_r232/gtdbtk.bac120.summary.tsv"
AR_GTDB = W / "05_refine/_gtdbtk_r232/gtdbtk.ar53.summary.tsv"

F3_HMM = HMM_DIR / "f3_rhodopsins.hmm"
F3_DMND = REF_DIR / "f3_rhodopsin_references.dmnd"

TYPE1_ACCS = {"PF01036"}
HELIO_ACCS = {"PF18761", "NF038020"}
GPCR_ACCS = {"PF00001", "PF00002"}


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
    print("TRACK F3: RHODOPSIN & HELIORHODOPSIN PILOT EXECUTION")
    print(f"Timestamp: {time.strftime('%Y-%m-%d %H:%M:%S UTC', time.gmtime())}")
    print(f"Cohort FASTA: {COHORT_FAA}")
    print(f"Output Directory: {PILOT_DIR}")
    print("=" * 80)

    if not COHORT_FAA.exists():
        raise FileNotFoundError(f"Cohort FASTA missing: {COHORT_FAA}")

    mag_to_cluster, tax_map = load_cohort_and_taxonomy()

    raw_domtbl = PILOT_DIR / "f3_domtblout.txt"
    raw_hmmsearch = PILOT_DIR / "f3_hmmsearch.txt"
    candidates_faa = PILOT_DIR / "f3_candidates.faa"
    candidates_tsv = PILOT_DIR / "f3_candidates_summary.tsv"
    rejections_tsv = PILOT_DIR / "f3_rejections_and_ambiguous.tsv"
    diamond_tsv = PILOT_DIR / "f3_diamond_reciprocal.tsv"
    metrics_json = PILOT_DIR / "f3_pilot_metrics.json"

    # Step 1: Run HMMER search directly on cohort FASTA
    print("\n[Step 1] Running hmmsearch against cohort FASTA (8 CPUs)...")
    hmm_cmd = [
        "hmmsearch",
        "--cpu", "8",
        "-E", "1e-3",
        "--domtblout", str(raw_domtbl),
        "-o", str(raw_hmmsearch),
        str(F3_HMM),
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

    # Group hits by protein ID
    by_protein = {}
    for h in hits:
        by_protein.setdefault(h["target_name"], []).append(h)

    print(f"  Unique proteins with domain hits: {len(by_protein)}")

    # Evaluate candidates and rejections
    preliminary_candidates = {}
    rejections = []

    for prot_id, prot_hits in by_protein.items():
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
        best_target_hit = best_type1 if score_type1 >= score_helio else best_helio

        if not best_target_hit:
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

        cov = best_target_hit["hmm_coverage"]
        ie = best_target_hit["dom_ievalue"]

        # Gate 1: Domain i-Evalue <= 1e-5
        if ie > 1e-5:
            rejections.append({
                "protein_id": prot_id,
                "reason": f"Failed_dom_ievalue_gate_{ie:.2e}>1e-5",
                "target_score": target_score,
                "gpcr_score": score_gpcr,
                "hmm_coverage": cov,
                "dom_ievalue": ie,
                "profile_acc": best_target_hit["profile_acc"],
            })
            continue

        # Gate 2: Domain coverage >= 0.70
        if cov < 0.70:
            rejections.append({
                "protein_id": prot_id,
                "reason": f"Failed_domain_coverage_gate_{cov:.2f}<0.70",
                "target_score": target_score,
                "gpcr_score": score_gpcr,
                "hmm_coverage": cov,
                "dom_ievalue": ie,
                "profile_acc": best_target_hit["profile_acc"],
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
                "profile_acc": best_target_hit["profile_acc"],
            })
            continue

        assigned_class = "Type-1_Microbial_Rhodopsin" if score_type1 >= score_helio else "Heliorhodopsin_Type-3"
        preliminary_candidates[prot_id] = {
            "assigned_class": assigned_class,
            "best_hit": best_target_hit,
            "best_type1_score": score_type1,
            "best_helio_score": score_helio,
            "best_gpcr_score": score_gpcr,
        }

    print(f"  Proteins passing HMM domain & discrimination gates: {len(preliminary_candidates)}")
    print(f"  Rejected/ambiguous hits: {len(rejections)}")

    # Step 3: Stream COHORT_FAA to extract candidate headers & sequences (check partial=00)
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
                "gpcr_score": cinfo["best_gpcr_score"],
                "hmm_coverage": cinfo["best_hit"]["hmm_coverage"],
                "dom_ievalue": cinfo["best_hit"]["dom_ievalue"],
                "profile_acc": cinfo["best_hit"]["profile_acc"],
            })
            continue

        hdr = pdata["header"]
        seq = pdata["seq"].rstrip("*")

        # Parse Prodigal completeness: partial=00
        partial_match = re.search(r"partial=(\d{2})", hdr)
        partial_flag = partial_match.group(1) if partial_match else "unknown"
        is_complete = 1 if partial_flag == "00" else 0

        # Parse coordinates
        fields = hdr.split(" # ")
        start_pos = fields[1] if len(fields) > 1 else "0"
        stop_pos = fields[2] if len(fields) > 2 else "0"
        strand_pos = fields[3] if len(fields) > 3 else "0"

        # Parse MAG and contig
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
                "gpcr_score": cinfo["best_gpcr_score"],
                "hmm_coverage": bh["hmm_coverage"],
                "dom_ievalue": bh["dom_ievalue"],
                "profile_acc": bh["profile_acc"],
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
            "best_type1_score": cinfo["best_type1_score"],
            "best_helio_score": cinfo["best_helio_score"],
            "gpcr_control_score": cinfo["best_gpcr_score"],
            "prodigal_partial": partial_flag,
            "is_complete_gene": is_complete,
            "sequence": seq,
        })

    print(f"  Final candidate proteins meeting all gates (including partial=00): {len(final_candidates)}")

    # Step 4: Write candidate FASTA (trailing '*' stripped)
    print(f"\n[Step 4] Emitting {len(final_candidates)} candidates to {candidates_faa.name}...")
    with open(candidates_faa, "w", encoding="utf-8") as f_out:
        for c in final_candidates:
            hdr_str = f">{c['protein_id']}|{c['mag']}|{c['assigned_class']}|{c['secondary_cluster_ani95']}"
            f_out.write(f"{hdr_str}\n{c['sequence']}\n")

    # Step 5: Reciprocal DIAMOND blastp against reference database
    diamond_matches = {}
    if F3_DMND.exists() and len(final_candidates) > 0:
        print("\n[Step 5] Running reciprocal DIAMOND blastp against f3_rhodopsin_references.dmnd...")
        diamond_matches = run_reciprocal_diamond(candidates_faa, F3_DMND, diamond_tsv)
        print(f"  Reciprocal matches obtained: {len(diamond_matches)}")

    for c in final_candidates:
        dm = diamond_matches.get(c["protein_id"], {})
        c["top_ref_id"] = dm.get("top_ref_id", "NO_HIT")
        c["top_ref_pident"] = dm.get("pident", 0.0)
        c["top_ref_evalue"] = dm.get("evalue", 1.0)
        c["top_ref_bitscore"] = dm.get("bitscore", 0.0)

    # Step 6: Write Candidates Summary TSV
    print(f"\n[Step 6] Writing summary TSV to {candidates_tsv.name}...")
    fieldnames = [
        "protein_id", "mag", "contig", "start", "stop", "strand",
        "secondary_cluster_ani95", "gtdbtk_taxonomy", "assigned_class",
        "best_profile_name", "profile_accession", "domain_score",
        "domain_ievalue", "hmm_coverage", "best_type1_score", "best_helio_score",
        "gpcr_control_score", "top_ref_id", "top_ref_pident", "top_ref_evalue", "top_ref_bitscore",
        "prodigal_partial", "is_complete_gene"
    ]
    with open(candidates_tsv, "w", newline="", encoding="utf-8") as f_out:
        writer = csv.DictWriter(f_out, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(final_candidates)

    # Write Rejections TSV
    rej_fields = ["protein_id", "reason", "target_score", "gpcr_score", "hmm_coverage", "dom_ievalue", "profile_acc"]
    with open(rejections_tsv, "w", newline="", encoding="utf-8") as f_out:
        writer = csv.DictWriter(f_out, fieldnames=rej_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rejections)

    # Step 7: Freeze executed code & provenance
    shutil.copy2(HERE / "03_f3_rhodopsin_pilot.py", PILOT_DIR / "executed_03_f3_pilot.py")
    
    provenance = {
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "script": "03_f3_rhodopsin_pilot.py",
        "script_sha256": compute_sha256(HERE / "03_f3_rhodopsin_pilot.py"),
        "cohort_faa": str(COHORT_FAA),
        "cohort_faa_sha256": compute_sha256(COHORT_FAA),
        "f3_hmm": str(F3_HMM),
        "f3_hmm_sha256": compute_sha256(F3_HMM),
        "f3_dmnd": str(F3_DMND),
        "f3_dmnd_sha256": compute_sha256(F3_DMND),
    }
    with open(PILOT_DIR / "provenance.json", "w", encoding="utf-8") as f_out:
        json.dump(provenance, f_out, indent=2)

    # Step 8: Compute pilot metrics
    type1_count = sum(1 for c in final_candidates if c["assigned_class"] == "Type-1_Microbial_Rhodopsin")
    helio_count = sum(1 for c in final_candidates if c["assigned_class"] == "Heliorhodopsin_Type-3")
    mags_with_cand = len(set(c["mag"] for c in final_candidates))
    clusters_with_cand = len(set(c["secondary_cluster_ani95"] for c in final_candidates if c["secondary_cluster_ani95"] != "Unknown"))

    metrics = {
        "track": "F3",
        "description": "Rhodopsin and Heliorhodopsin Discovery Pilot",
        "cohort_genomes": 256,
        "cohort_proteins_scanned": 956344,
        "raw_domtbl_hits": len(hits),
        "unique_proteins_with_hits": len(by_protein),
        "rejected_or_ambiguous_hits": len(rejections),
        "final_passing_candidates": len(final_candidates),
        "type1_microbial_rhodopsin": type1_count,
        "heliorhodopsin_type3": helio_count,
        "distinct_mags_with_candidates": mags_with_cand,
        "distinct_secondary_clusters_ani95": clusters_with_cand,
        "timing_hmmsearch_seconds": t_hmm,
        "timing_total_seconds": time.time() - start_time,
    }
    with open(metrics_json, "w", encoding="utf-8") as f_out:
        json.dump(metrics, f_out, indent=2)

    print(f"\nWrote pilot metrics to {metrics_json.name}")
    print(f"=== Track F3 Pilot Finished in {time.time() - start_time:.1f}s ===")


if __name__ == "__main__":
    main()
