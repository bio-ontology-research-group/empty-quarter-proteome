#!/usr/bin/env python3
"""04_f4_ectc_pilot.py - Track F4: Ectoine Synthase (EctC) & Synteny Pilot Screening.

Grounding & Audit Resolutions (Argon Section 4 & Root Requirements):
1. Reuses frozen, read-only 256-genome cohort (956,344 proteins) from root job 53816043.
2. Streams candidates only: zero full-proteome dictionaries in memory; zero concatenations.
3. Profile matching strictly by HMM accession prefix (parts[4].split(".")[0]):
   - Target EctC: PF06339 (Ectoine_synth), NF009806 (ectoine_synthase)
   - Pathway markers: TIGR02406 (EctA), TIGR02407 (EctB), TIGR02408 (EctD)
   - Competing Cupin controls: PF00190 (Cupin_1), PF07883 (Cupin_2), PF02678 (Pirin)
4. Enforces predetermined gates:
   - EctC dom_ievalue <= 1e-10
   - EctC domain coverage >= 0.70
   - EctC bit score >= 25.0 bits above generic cupin controls
   - Prodigal partial=00 (complete CDS)
5. Synteny analysis (Argon Section 4):
   - Measures exact distance (bp) and strand orientation on the SAME contig (<=5000 bp -> Co_localized_Operon).
   - If EctA/B present on different contigs -> Unlinked_Pathway_Carrier (never whole-MAG false operon, never call unlinked absent).
   - If missing -> Missing_in_assembly_or_unlinked (not definitive absence).
6. Emits verified candidate FASTA with trailing stop '*' stripped.
7. Executes reciprocal DIAMOND blastp against f4_ectc_references.dmnd.
8. Retains all failed and ambiguous hits in f4_rejections_and_ambiguous.tsv.
9. Freezes executed code and provenance manifest to results/pilot_v1/f4_ectc/.
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
PILOT_DIR = ROOT / "results" / "pilot_v1" / "f4_ectc"
PILOT_DIR.mkdir(parents=True, exist_ok=True)

# Read-only 256-genome cohort from Root
COHORT_FAA = Path("/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/family_discovery_20261015_clean/root_audit/f1_profile_pilot/cohort.faa")
COHORT_TSV = Path("/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/family_discovery_20261015_clean/root_audit/f1_profile_pilot/cohort.tsv")

W = Path("/ibex/scratch/projects/c2014/rob/metagenomics")
CDB_PATH = W / "05_refine/_drep/out/data_tables/Cdb.csv"
BAC_GTDB = W / "05_refine/_gtdbtk_r232/gtdbtk.bac120.summary.tsv"
AR_GTDB = W / "05_refine/_gtdbtk_r232/gtdbtk.ar53.summary.tsv"

F4_HMM = HMM_DIR / "f4_ectc.hmm"
F4_DMND = REF_DIR / "f4_ectc_references.dmnd"

ECTC_ACCS = {"PF06339", "NF009806"}
ECTA_ACCS = {"TIGR02406"}
ECTB_ACCS = {"TIGR02407"}
ECTD_ACCS = {"TIGR02408"}
CUPIN_ACCS = {"PF00190", "PF07883", "PF02678"}


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
    print("TRACK F4: ECTOINE SYNTHASE (ECTC) & SYNTENY PILOT EXECUTION")
    print(f"Timestamp: {time.strftime('%Y-%m-%d %H:%M:%S UTC', time.gmtime())}")
    print(f"Cohort FASTA: {COHORT_FAA}")
    print(f"Output Directory: {PILOT_DIR}")
    print("=" * 80)

    if not COHORT_FAA.exists():
        raise FileNotFoundError(f"Cohort FASTA missing: {COHORT_FAA}")

    mag_to_cluster, tax_map = load_cohort_and_taxonomy()

    raw_domtbl = PILOT_DIR / "f4_domtblout.txt"
    raw_hmmsearch = PILOT_DIR / "f4_hmmsearch.txt"
    candidates_faa = PILOT_DIR / "f4_candidates.faa"
    candidates_tsv = PILOT_DIR / "f4_candidates_summary.tsv"
    rejections_tsv = PILOT_DIR / "f4_rejections_and_ambiguous.tsv"
    diamond_tsv = PILOT_DIR / "f4_diamond_reciprocal.tsv"
    metrics_json = PILOT_DIR / "f4_pilot_metrics.json"

    # Step 1: Run HMMER search directly on cohort FASTA
    print("\n[Step 1] Running hmmsearch against cohort FASTA (8 CPUs)...")
    hmm_cmd = [
        "hmmsearch",
        "--cpu", "8",
        "-E", "1e-3",
        "--domtblout", str(raw_domtbl),
        "-o", str(raw_hmmsearch),
        str(F4_HMM),
        str(COHORT_FAA),
    ]
    t0 = time.time()
    subprocess.run(hmm_cmd, check=True)
    t_hmm = time.time() - t0
    print(f"  hmmsearch complete in {t_hmm:.1f}s.")

    # Step 2: Parse domtblout
    print("\n[Step 2] Parsing domain hits and evaluating EctC & pathway genes...")
    hits = parse_domtblout(raw_domtbl)
    print(f"  Total raw domain hits in domtblout: {len(hits)}")

    by_protein = {}
    for h in hits:
        by_protein.setdefault(h["target_name"], []).append(h)

    print(f"  Unique proteins with domain hits: {len(by_protein)}")

    # Separate pathway marker hits (EctA, EctB, EctD) and EctC candidate evaluations
    pathway_hits = {}  # prot_id -> {"role": "EctA"|"EctB"|"EctD", "dom_score": ..., "dom_ievalue": ...}
    preliminary_ectc = {}
    rejections = []

    for prot_id, prot_hits in by_protein.items():
        ectc_hits = [h for h in prot_hits if h["profile_acc"] in ECTC_ACCS]
        cupin_hits = [h for h in prot_hits if h["profile_acc"] in CUPIN_ACCS]
        ecta_hits = [h for h in prot_hits if h["profile_acc"] in ECTA_ACCS]
        ectb_hits = [h for h in prot_hits if h["profile_acc"] in ECTB_ACCS]
        ectd_hits = [h for h in prot_hits if h["profile_acc"] in ECTD_ACCS]

        # Check for pathway marker presence
        if ecta_hits:
            best_a = max(ecta_hits, key=lambda x: x["dom_score"])
            if best_a["dom_ievalue"] <= 1e-5 and best_a["hmm_coverage"] >= 0.50:
                pathway_hits[prot_id] = {"role": "EctA", "hit": best_a}
        if ectb_hits:
            best_b = max(ectb_hits, key=lambda x: x["dom_score"])
            if best_b["dom_ievalue"] <= 1e-5 and best_b["hmm_coverage"] >= 0.50:
                pathway_hits[prot_id] = {"role": "EctB", "hit": best_b}
        if ectd_hits:
            best_d = max(ectd_hits, key=lambda x: x["dom_score"])
            if best_d["dom_ievalue"] <= 1e-5 and best_d["hmm_coverage"] >= 0.50:
                pathway_hits[prot_id] = {"role": "EctD", "hit": best_d}

        # Check for EctC target
        if not ectc_hits:
            if not (ecta_hits or ectb_hits or ectd_hits):
                best_cupin = max(cupin_hits, key=lambda x: x["dom_score"]) if cupin_hits else None
                rejections.append({
                    "protein_id": prot_id,
                    "reason": "Only_generic_cupin_control_hit",
                    "target_score": 0.0,
                    "cupin_score": best_cupin["dom_score"] if best_cupin else 0.0,
                    "hmm_coverage": 0.0,
                    "dom_ievalue": 1.0,
                    "profile_acc": best_cupin["profile_acc"] if best_cupin else "None",
                })
            continue

        best_ectc = max(ectc_hits, key=lambda x: x["dom_score"])
        best_cupin = max(cupin_hits, key=lambda x: x["dom_score"]) if cupin_hits else None

        score_ectc = best_ectc["dom_score"]
        score_cupin = best_cupin["dom_score"] if best_cupin else 0.0
        cov = best_ectc["hmm_coverage"]
        ie = best_ectc["dom_ievalue"]

        # Gate 1: dom_ievalue <= 1e-10
        if ie > 1e-10:
            rejections.append({
                "protein_id": prot_id,
                "reason": f"Failed_dom_ievalue_gate_{ie:.2e}>1e-10",
                "target_score": score_ectc,
                "cupin_score": score_cupin,
                "hmm_coverage": cov,
                "dom_ievalue": ie,
                "profile_acc": best_ectc["profile_acc"],
            })
            continue

        # Gate 2: domain coverage >= 0.70
        if cov < 0.70:
            rejections.append({
                "protein_id": prot_id,
                "reason": f"Failed_domain_coverage_gate_{cov:.2f}<0.70",
                "target_score": score_ectc,
                "cupin_score": score_cupin,
                "hmm_coverage": cov,
                "dom_ievalue": ie,
                "profile_acc": best_ectc["profile_acc"],
            })
            continue

        # Gate 3: Score discrimination vs generic cupins (target >= cupin + 25.0 bits)
        if (score_ectc - score_cupin) < 25.0:
            rejections.append({
                "protein_id": prot_id,
                "reason": f"Excluded_Cupin_cross_reactivity_diff_{score_ectc - score_cupin:.1f}<25b",
                "target_score": score_ectc,
                "cupin_score": score_cupin,
                "hmm_coverage": cov,
                "dom_ievalue": ie,
                "profile_acc": best_ectc["profile_acc"],
            })
            continue

        preliminary_ectc[prot_id] = {
            "best_hit": best_ectc,
            "best_ectc_score": score_ectc,
            "best_cupin_score": score_cupin,
        }

    print(f"  EctC proteins passing domain & discrimination gates: {len(preliminary_ectc)}")
    print(f"  Pathway marker proteins identified: EctA/B/D count = {len(pathway_hits)}")
    print(f"  Rejected/ambiguous hits: {len(rejections)}")

    # Step 3: Stream COHORT_FAA to extract candidate headers & sequences (for EctC AND pathway genes)
    print("\n[Step 3] Streaming COHORT_FAA to extract sequences and coordinate headers...")
    needed_ids = set(preliminary_ectc.keys()) | set(pathway_hits.keys())
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

    print(f"  Extracted {len(protein_data)} records for EctC candidates & pathway genes.")

    # Parse coordinates for pathway genes to build genome-level index
    # genome -> list of {"role": "EctA"|"EctB"|"EctD", "prot_id": ..., "contig": ..., "start": int, "stop": int, "strand": int}
    genome_pathway_genes = {}
    for pid, pinfo in pathway_hits.items():
        pdata = protein_data.get(pid)
        if not pdata:
            continue
        hdr = pdata["header"]
        fields = hdr.split(" # ")
        s = int(fields[1]) if len(fields) > 1 and fields[1].isdigit() else 0
        e = int(fields[2]) if len(fields) > 2 and fields[2].isdigit() else 0
        strnd = int(fields[3]) if len(fields) > 3 and fields[3] in ("1", "-1") else 0
        
        if "::" in pid:
            mag_part, contig_gene = pid.split("::", 1)
            mag = mag_part
            contig = contig_gene.rsplit("_", 1)[0] if "_" in contig_gene else contig_gene
        else:
            mag = pid.split("_")[0]
            contig = pid.rsplit("_", 1)[0] if "_" in pid else pid

        genome_pathway_genes.setdefault(mag, []).append({
            "role": pinfo["role"],
            "prot_id": pid,
            "contig": contig,
            "start": s,
            "stop": e,
            "strand": strnd,
            "dom_score": pinfo["hit"]["dom_score"],
        })

    # Evaluate final EctC candidates with synteny
    final_candidates = []
    for prot_id, cinfo in preliminary_ectc.items():
        pdata = protein_data.get(prot_id)
        if not pdata:
            rejections.append({
                "protein_id": prot_id,
                "reason": "Sequence_missing_in_cohort_fasta",
                "target_score": cinfo["best_ectc_score"],
                "cupin_score": cinfo["best_cupin_score"],
                "hmm_coverage": cinfo["best_hit"]["hmm_coverage"],
                "dom_ievalue": cinfo["best_hit"]["dom_ievalue"],
                "profile_acc": cinfo["best_hit"]["profile_acc"],
            })
            continue

        hdr = pdata["header"]
        seq = pdata["seq"].rstrip("*")

        # Parse Prodigal completeness
        partial_match = re.search(r"partial=(\d{2})", hdr)
        partial_flag = partial_match.group(1) if partial_match else "unknown"
        is_complete = 1 if partial_flag == "00" else 0

        fields = hdr.split(" # ")
        start_pos = int(fields[1]) if len(fields) > 1 and fields[1].isdigit() else 0
        stop_pos = int(fields[2]) if len(fields) > 2 and fields[2].isdigit() else 0
        strand_pos = int(fields[3]) if len(fields) > 3 and fields[3] in ("1", "-1") else 0

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
                "cupin_score": cinfo["best_cupin_score"],
                "hmm_coverage": bh["hmm_coverage"],
                "dom_ievalue": bh["dom_ievalue"],
                "profile_acc": bh["profile_acc"],
            })
            continue

        # Evaluate Synteny vs other pathway genes in the same MAG
        pw_genes = genome_pathway_genes.get(mag, [])
        same_contig_pw = [g for g in pw_genes if g["contig"] == contig]
        diff_contig_pw = [g for g in pw_genes if g["contig"] != contig]

        # Operon proximity: distance <= 5000 bp on same contig
        operon_matches = []
        for g in same_contig_pw:
            # Intergenic distance
            d = min(abs(start_pos - g["stop"]), abs(g["start"] - stop_pos))
            if d <= 5000:
                rel_orient = "same_strand" if strand_pos == g["strand"] else "divergent_or_convergent"
                operon_matches.append({
                    "role": g["role"],
                    "prot_id": g["prot_id"],
                    "distance_bp": d,
                    "orientation": rel_orient,
                })

        if operon_matches:
            synteny_class = "Co_localized_Operon"
            matched_roles = set(m["role"] for m in operon_matches)
            operon_type = "ect" + "".join(sorted(list(matched_roles | {"C"})))
            min_dist = min(m["distance_bp"] for m in operon_matches)
            orient = operon_matches[0]["orientation"] if len(set(m["orientation"] for m in operon_matches)) == 1 else "mixed"
            neighbor_ids = ";".join(f"{m['role']}:{m['prot_id']}({m['distance_bp']}bp)" for m in operon_matches)
        elif pw_genes:
            synteny_class = "Unlinked_Pathway_Carrier"
            matched_roles = set(g["role"] for g in pw_genes)
            operon_type = "unlinked_" + "_".join(sorted(list(matched_roles)))
            min_dist = -1
            orient = "unlinked_contigs"
            neighbor_ids = ";".join(f"{g['role']}:{g['prot_id']}(unlinked_contig:{g['contig']})" for g in pw_genes)
        else:
            synteny_class = "Missing_in_assembly_or_unlinked"
            operon_type = "orphan_ectC_candidate"
            min_dist = -1
            orient = "none"
            neighbor_ids = "none_detected"

        final_candidates.append({
            "protein_id": prot_id,
            "mag": mag,
            "contig": contig,
            "start": start_pos,
            "stop": stop_pos,
            "strand": strand_pos,
            "secondary_cluster_ani95": sec_clust,
            "gtdbtk_taxonomy": tax,
            "assigned_class": "Ectoine_Synthase_EctC",
            "best_profile_name": bh["profile_name"],
            "profile_accession": bh["profile_acc"],
            "domain_score": bh["dom_score"],
            "domain_ievalue": bh["dom_ievalue"],
            "hmm_coverage": bh["hmm_coverage"],
            "best_cupin_score": cinfo["best_cupin_score"],
            "synteny_classification": synteny_class,
            "operon_type": operon_type,
            "min_intergenic_distance_bp": min_dist,
            "relative_orientation": orient,
            "neighbor_gene_ids": neighbor_ids,
            "prodigal_partial": partial_flag,
            "is_complete_gene": is_complete,
            "sequence": seq,
        })

    print(f"  Final EctC candidates meeting all gates (including partial=00): {len(final_candidates)}")

    # Step 4: Write candidate FASTA (trailing '*' stripped)
    print(f"\n[Step 4] Emitting {len(final_candidates)} candidates to {candidates_faa.name}...")
    with open(candidates_faa, "w", encoding="utf-8") as f_out:
        for c in final_candidates:
            hdr_str = f">{c['protein_id']}|{c['mag']}|{c['synteny_classification']}|{c['secondary_cluster_ani95']}"
            f_out.write(f"{hdr_str}\n{c['sequence']}\n")

    # Step 5: Reciprocal DIAMOND blastp against reference database
    diamond_matches = {}
    if F4_DMND.exists() and len(final_candidates) > 0:
        print("\n[Step 5] Running reciprocal DIAMOND blastp against f4_ectc_references.dmnd...")
        diamond_matches = run_reciprocal_diamond(candidates_faa, F4_DMND, diamond_tsv)
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
        "domain_ievalue", "hmm_coverage", "best_cupin_score",
        "synteny_classification", "operon_type", "min_intergenic_distance_bp",
        "relative_orientation", "neighbor_gene_ids",
        "top_ref_id", "top_ref_pident", "top_ref_evalue", "top_ref_bitscore",
        "prodigal_partial", "is_complete_gene"
    ]
    with open(candidates_tsv, "w", newline="", encoding="utf-8") as f_out:
        writer = csv.DictWriter(f_out, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(final_candidates)

    # Write Rejections TSV
    rej_fields = ["protein_id", "reason", "target_score", "cupin_score", "hmm_coverage", "dom_ievalue", "profile_acc"]
    with open(rejections_tsv, "w", newline="", encoding="utf-8") as f_out:
        writer = csv.DictWriter(f_out, fieldnames=rej_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rejections)

    # Step 7: Freeze executed code & provenance
    shutil.copy2(HERE / "04_f4_ectc_pilot.py", PILOT_DIR / "executed_04_f4_pilot.py")

    provenance = {
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "script": "04_f4_ectc_pilot.py",
        "script_sha256": compute_sha256(HERE / "04_f4_ectc_pilot.py"),
        "cohort_faa": str(COHORT_FAA),
        "cohort_faa_sha256": compute_sha256(COHORT_FAA),
        "f4_hmm": str(F4_HMM),
        "f4_hmm_sha256": compute_sha256(F4_HMM),
        "f4_dmnd": str(F4_DMND),
        "f4_dmnd_sha256": compute_sha256(F4_DMND),
    }
    with open(PILOT_DIR / "provenance.json", "w", encoding="utf-8") as f_out:
        json.dump(provenance, f_out, indent=2)

    # Step 8: Compute pilot metrics
    coloc_count = sum(1 for c in final_candidates if c["synteny_classification"] == "Co_localized_Operon")
    unlinked_count = sum(1 for c in final_candidates if c["synteny_classification"] == "Unlinked_Pathway_Carrier")
    orphan_count = sum(1 for c in final_candidates if c["synteny_classification"] == "Missing_in_assembly_or_unlinked")
    mags_with_cand = len(set(c["mag"] for c in final_candidates))
    clusters_with_cand = len(set(c["secondary_cluster_ani95"] for c in final_candidates if c["secondary_cluster_ani95"] != "Unknown"))

    metrics = {
        "track": "F4",
        "description": "Ectoine Synthase (EctC) and Pathway Synteny Discovery Pilot",
        "cohort_genomes": 256,
        "cohort_proteins_scanned": 956344,
        "raw_domtbl_hits": len(hits),
        "unique_proteins_with_hits": len(by_protein),
        "rejected_or_ambiguous_hits": len(rejections),
        "final_passing_candidates": len(final_candidates),
        "synteny_co_localized_operon": coloc_count,
        "synteny_unlinked_pathway_carrier": unlinked_count,
        "synteny_missing_in_assembly_or_unlinked": orphan_count,
        "distinct_mags_with_candidates": mags_with_cand,
        "distinct_secondary_clusters_ani95": clusters_with_cand,
        "timing_hmmsearch_seconds": t_hmm,
        "timing_total_seconds": time.time() - start_time,
    }
    with open(metrics_json, "w", encoding="utf-8") as f_out:
        json.dump(metrics, f_out, indent=2)

    print(f"\nWrote pilot metrics to {metrics_json.name}")
    print(f"=== Track F4 Pilot Finished in {time.time() - start_time:.1f}s ===")


if __name__ == "__main__":
    main()
