#!/usr/bin/env python3
"""02_f1_screen.py - CPF Family Identity & Diversity Census across 256 Deterministic Representatives.

Audit-Compliant Pilot Screen (Root & Argon Audit Gate, 2026-10-15):
1. Deterministic SHA256-ranked cohort of 256 representative proteomes from the 5,438 Wdb ANI95 winners:
   - Sorts 5,438 representatives alphabetically, computes SHA256(mag_name), sorts by hash, selects top 256.
   - Retains exact manifest in results_pilot_v3/f1_pilot_cohort_256.tsv before search.
2. Extracts and concatenates only the 256 pilot proteomes into scratch.
3. DIAMOND homology search against curated CPF references (E <= 1e-3, ultra-sensitive).
4. Dedicated HMMER hmmsearch with separate --domtblout and -o report files:
   - Exact domtblout schema: parts[3]=query_name, parts[4]=query_acc (PF00875, PF03441, PF04244, PF04055).
   - Domain score (parts[13]), domain evalue (parts[12]), and domain coverage.
5. Strict candidate gating:
   - is_cpf = 1 ONLY if domain E-value <= 1e-4 to a CPF profile (PF00875, PF03441, PF04244) and dom_score >= 25.0.
   - Candidates without CPF domain are rejected from confirmed CPF set.
   - Outgroup Radical SAM SplB (PF04055) is separated.
   - Broad family categories: CPD, 6-4, cryptochrome, SPL, SplB, CPF_unclassified (no heuristic fallbacks).
6. Evaluates Prodigal coding completeness (partial=00 is complete; else partial=10, 01, 11, or unknown).
7. Joins with GTDB-Tk R232 taxonomy and sample metadata.
8. Writes detailed f1_candidates.tsv, f1_subfamily_summary.tsv, and f1_complete_candidates.faa in results_pilot_v3/.
"""

import csv
import gzip
import hashlib
import os
import subprocess
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
REF_DIR = ROOT / "references"
RES_DIR = ROOT / "results_pilot_v3"
RAW_RES_DIR = RES_DIR / "raw_hits"
RES_DIR.mkdir(parents=True, exist_ok=True)
RAW_RES_DIR.mkdir(parents=True, exist_ok=True)

W = Path("/ibex/scratch/projects/c2014/rob/metagenomics")
REP_PROTEOMES_DIR = W / "06_darkmatter/drep_rep_proteomes"
DREP_WDB = W / "05_refine/_drep/out/data_tables/Wdb.csv"
GTDB_BAC = W / "05_refine/_gtdbtk_r232/gtdbtk.bac120.summary.tsv"
GTDB_AR = W / "05_refine/_gtdbtk_r232/gtdbtk.ar53.summary.tsv"
SAMPLE_META = W / "metadata/sample_metadata.tsv"

CPF_HMM = REF_DIR / "cpf_all.hmm"
CPF_DMND = REF_DIR / "cpf_references.dmnd"
REF_META_FILE = REF_DIR / "verified_references_metadata.tsv"

DIAMOND_BIN = Path("/ibex/scratch/projects/c2014/rob/miniconda3/envs/eggnog/bin/diamond")
if not DIAMOND_BIN.exists():
    DIAMOND_BIN = Path("diamond")

THREADS = int(os.environ.get("SLURM_CPUS_PER_TASK", "16"))
PILOT_COHORT_SIZE = 256


def load_taxonomy():
    tax_map = {}
    for p in [GTDB_BAC, GTDB_AR]:
        if not p.exists():
            continue
        with open(p) as f:
            for r in csv.DictReader(f, delimiter="\t"):
                mag = r["user_genome"]
                tax = r["classification"]
                parts = tax.split(";")
                d = p_ = c = o = f_ = g = s = ""
                for part in parts:
                    if part.startswith("d__"): d = part[3:]
                    elif part.startswith("p__"): p_ = part[3:]
                    elif part.startswith("c__"): c = part[3:]
                    elif part.startswith("o__"): o = part[3:]
                    elif part.startswith("f__"): f_ = part[3:]
                    elif part.startswith("g__"): g = part[3:]
                    elif part.startswith("s__"): s = part[3:]
                tax_map[mag] = {
                    "classification": tax,
                    "domain": d,
                    "phylum": p_,
                    "class": c,
                    "order": o,
                    "family": f_,
                    "genus": g,
                    "species": s,
                }
    return tax_map


def load_sample_metadata():
    s_map = {}
    if SAMPLE_META.exists():
        with open(SAMPLE_META) as f:
            for r in csv.DictReader(f, delimiter="\t"):
                s_map[r["sample"]] = {
                    "compartment": r.get("compartment", "Unknown"),
                    "trip": r.get("trip", ""),
                    "site": r.get("site", ""),
                }
    return s_map


def load_reference_subfamilies():
    subfam_map = {}
    if REF_META_FILE.exists():
        with open(REF_META_FILE) as f:
            for r in csv.DictReader(f, delimiter="\t"):
                if r["track"] == "F1_CPF":
                    subfam_map[r["accession"]] = {
                        "class_assignment": r["class_assignment"],
                        "subcategory": r["subcategory"],
                        "category": r["category"],
                        "protein_name": r["protein_name"],
                        "organism": r["organism"],
                        "ec_number": r["ec_number"],
                        "has_direct_assay": r.get("has_direct_assay_evidence", "0"),
                        "has_experimental": r.get("has_experimental_evidence", "0"),
                    }
    return subfam_map


def select_deterministic_pilot_cohort():
    """Deterministically select 256 representatives from 5,438 Wdb winners using SHA256 ranking."""
    all_reps = []
    with open(DREP_WDB) as f:
        for r in csv.DictReader(f):
            mag = r["genome"].replace(".fa", "").replace(".fasta", "")
            cluster = r["cluster"]
            all_reps.append((mag, cluster))

    assert len(all_reps) == 5438, f"Expected 5,438 Wdb representatives, found {len(all_reps)}!"
    all_reps.sort(key=lambda x: x[0])

    # Hash each genome name with SHA256 and sort
    hashed_reps = []
    for mag, cluster in all_reps:
        h = hashlib.sha256(mag.encode("utf-8")).hexdigest()
        hashed_reps.append((h, mag, cluster))

    hashed_reps.sort(key=lambda x: x[0])
    selected = hashed_reps[:PILOT_COHORT_SIZE]
    return selected


def main():
    t0 = time.time()
    print("=== [F1 Pilot Screen - Audit Compliant v3] Photolyase / Cryptochrome Diversity Pilot ===")
    print(f"Threads: {THREADS}")
    print(f"Cohort Size: {PILOT_COHORT_SIZE} deterministic representatives")
    print(f"Representative Proteomes Source: {REP_PROTEOMES_DIR}")
    print(f"Reference Metadata: {REF_META_FILE}")

    tax_map = load_taxonomy()
    meta_map = load_sample_metadata()
    subfam_map = load_reference_subfamilies()
    print(f"Loaded taxonomy for {len(tax_map)} MAGs, metadata for {len(meta_map)} samples, {len(subfam_map)} reference subfamilies.")
    assert len(subfam_map) > 0, "No reference subfamilies loaded! Check verified_references_metadata.tsv."

    # 1. Deterministic Pilot Cohort Selection
    print("\n--- Step 1: Deterministic Pilot Cohort Selection (SHA256 Ranking) ---")
    pilot_cohort = select_deterministic_pilot_cohort()
    cohort_mags = set(x[1] for x in pilot_cohort)
    print(f"Selected {len(pilot_cohort)} representative MAGs.")

    # 2. Extract and Concatenate Pilot Proteomes
    tmp_dir = ROOT / "tmp_f1"
    tmp_dir.mkdir(parents=True, exist_ok=True)
    pilot_faa = tmp_dir / "f1_pilot_256_reps.faa"

    cohort_manifest_tsv = RES_DIR / "f1_pilot_cohort_256.tsv"
    cohort_rows = []
    total_proteins_scanned = 0

    with open(pilot_faa, "w") as out_f:
        for rank, (h, mag, cluster) in enumerate(pilot_cohort, 1):
            faa_path = REP_PROTEOMES_DIR / f"{mag}.faa.gz"
            assert faa_path.exists(), f"Proteome archive not found for {mag} at {faa_path}!"
            mag_proteins = 0
            with gzip.open(faa_path, "rt") as in_f:
                for line in in_f:
                    if line.startswith(">"):
                        mag_proteins += 1
                        total_proteins_scanned += 1
                        pid = line[1:].strip().split()[0]
                        if not pid.startswith(mag_name + "::"):
                            line = f">{mag}::{line[1:]}"
                    out_f.write(line)

            sample = mag.split("__")[0] if "__" in mag else ""
            s_info = meta_map.get(sample, {})
            t_info = tax_map.get(mag, {})
            cohort_rows.append({
                "rank": rank,
                "mag": mag,
                "sha256_rank_hash": h,
                "secondary_cluster": cluster,
                "domain": t_info.get("domain", "Unknown"),
                "phylum": t_info.get("phylum", "Unknown"),
                "class": t_info.get("class", "Unknown"),
                "order": t_info.get("order", "Unknown"),
                "family": t_info.get("family", "Unknown"),
                "genus": t_info.get("genus", "Unknown"),
                "species": t_info.get("species", "Unknown"),
                "sample": sample,
                "compartment": s_info.get("compartment", "Unknown"),
                "total_proteins": mag_proteins,
            })

    with open(cohort_manifest_tsv, "w") as f:
        writer = csv.DictWriter(f, fieldnames=list(cohort_rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(cohort_rows)

    print(f"Wrote cohort manifest ({len(cohort_rows)} MAGs, {total_proteins_scanned} total proteins) to {cohort_manifest_tsv.name}.")
    print(f"Concatenated pilot proteomes to {pilot_faa.name}.")

    # 3. DIAMOND Search against Curated References
    print("\n--- Step 2: DIAMOND Homology Search Against Curated CPF Reference Set ---")
    dmnd_out = RAW_RES_DIR / "f1_pilot_diamond_hits.tsv"
    cmd = [
        str(DIAMOND_BIN), "blastp",
        "-d", str(CPF_DMND),
        "-q", str(pilot_faa),
        "-o", str(dmnd_out),
        "--ultra-sensitive",
        "-e", "1e-3",
        "--max-target-seqs", "5",
        "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
        "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen", "qcovhsp", "scovhsp",
        "--threads", str(THREADS),
    ]
    dmnd_start = time.time()
    subprocess.run(cmd, check=True)
    print(f"DIAMOND search finished in {time.time() - dmnd_start:.1f}s.")

    diamond_hits_by_query = defaultdict(list)
    with open(dmnd_out) as f:
        for line in f:
            parts = line.strip().split("\t")
            if not parts or len(parts) < 16:
                continue
            qid = parts[0]
            diamond_hits_by_query[qid].append({
                "ref_id": parts[1],
                "pident": float(parts[2]),
                "aln_len": int(parts[3]),
                "evalue": float(parts[10]),
                "bitscore": float(parts[11]),
                "qlen": int(parts[12]),
                "slen": int(parts[13]),
                "qcovhsp": float(parts[14]),
                "scovhsp": float(parts[15]),
            })

    print(f"Identified {len(diamond_hits_by_query)} unique query candidate hits.")

    # 4. Extract Candidate Sequences
    print("\n--- Step 3: Extracting Full Candidate Sequences & Parsing Prodigal Header ---")
    cand_seqs = {}
    cand_metadata = {}

    with open(pilot_faa) as f:
        cur_qid = None
        cur_hdr = None
        seq_parts = []
        for line in f:
            if line.startswith(">"):
                if cur_qid and cur_qid in diamond_hits_by_query:
                    cand_seqs[cur_qid] = "".join(seq_parts)
                hdr = line[1:].strip()
                cur_qid = hdr.split()[0]
                cur_hdr = hdr
                seq_parts = []
                if cur_qid in diamond_hits_by_query:
                    partial = "unknown"
                    start_bp = 0
                    stop_bp = 0
                    strand = 1
                    tokens = hdr.split(" # ")
                    if len(tokens) >= 4:
                        try:
                            start_bp = int(tokens[1])
                            stop_bp = int(tokens[2])
                            strand = int(tokens[3])
                        except ValueError:
                            pass
                    if len(tokens) >= 5:
                        for item in tokens[4].split(";"):
                            if item.startswith("partial="):
                                partial = item.split("=")[1]

                    mag_name = cur_qid.split("::")[0] if "::" in cur_qid else ""
                    cand_metadata[cur_qid] = {
                        "raw_header": hdr,
                        "mag": mag_name,
                        "start_bp": start_bp,
                        "stop_bp": stop_bp,
                        "strand": strand,
                        "partial": partial,
                        "is_complete": (partial == "00"),
                    }
            else:
                seq_parts.append(line.strip())
        if cur_qid and cur_qid in diamond_hits_by_query:
            cand_seqs[cur_qid] = "".join(seq_parts)

    cand_faa = RAW_RES_DIR / "f1_pilot_candidates.faa"
    with open(cand_faa, "w") as out_f:
        for qid, seq in cand_seqs.items():
            out_f.write(f">{qid}\n{seq}\n")
    print(f"Extracted {len(cand_seqs)} candidate sequences to {cand_faa.name}.")

    # 5. Dedicated HMMER hmmsearch
    print("\n--- Step 4: Running hmmsearch with Dedicated Outputs ---")
    hmm_domtbl = RAW_RES_DIR / "f1_cpf_hmm_domtbl.txt"
    hmm_report = RAW_RES_DIR / "f1_cpf_hmm_report.txt"
    cmd = [
        "hmmsearch",
        "--domtblout", str(hmm_domtbl),
        "-o", str(hmm_report),
        "--domE", "1e-4",
        "--cpu", str(THREADS),
        str(CPF_HMM),
        str(cand_faa),
    ]
    hmm_start = time.time()
    subprocess.run(cmd, check=True)
    print(f"hmmsearch completed in {time.time() - hmm_start:.1f}s.")

    # Parse domtblout: query_name is parts[3], query_acc is parts[4]
    hmm_hits = defaultdict(list)
    with open(hmm_domtbl) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 23:
                continue
            cand_id = parts[0]
            hmm_name = parts[3]
            hmm_acc = parts[4].split(".")[0]
            qlen = int(parts[5])
            full_eval = float(parts[6])
            full_score = float(parts[7])
            dom_eval = float(parts[12])
            dom_score = float(parts[13])
            hmm_from = int(parts[15])
            hmm_to = int(parts[16])
            coverage = (hmm_to - hmm_from + 1) / max(1, qlen)
            hmm_hits[cand_id].append({
                "hmm_acc": hmm_acc,
                "hmm_name": hmm_name,
                "full_evalue": full_eval,
                "full_score": full_score,
                "dom_evalue": dom_eval,
                "dom_score": dom_score,
                "hmm_coverage": coverage,
            })
    print(f"HMM matches recorded for {len(hmm_hits)} candidate sequences.")

    # 6. Biological Candidate Gating & Subfamily Assignment
    print("\n--- Step 5: Strict Gating & Broad Family Classification ---")
    final_candidates = []
    cpf_profile_accs = {"PF00875", "PF03441", "PF04244"}

    for qid, hit_list in diamond_hits_by_query.items():
        meta = cand_metadata.get(qid, {})
        seq = cand_seqs.get(qid, "")
        h_list = hmm_hits.get(qid, [])

        top_hit = hit_list[0]
        top_ref = top_hit["ref_id"]
        ref_entry = subfam_map.get(top_ref, {})
        assigned_broad = ref_entry.get("class_assignment", "CPF_unclassified")

        matched_hmms = sorted(set(h["hmm_acc"] for h in h_list))
        cpf_hmms = [h for h in h_list if h["hmm_acc"] in cpf_profile_accs and h["dom_evalue"] <= 1e-4]
        splb_hmms = [h for h in h_list if h["hmm_acc"] == "PF04055" and h["dom_evalue"] <= 1e-4]

        best_cpf_score = max((h["dom_score"] for h in cpf_hmms), default=0.0)
        best_cpf_eval = min((h["dom_evalue"] for h in cpf_hmms), default=1.0)
        best_cpf_cov = max((h["hmm_coverage"] for h in cpf_hmms), default=0.0)

        best_splb_score = max((h["dom_score"] for h in splb_hmms), default=0.0)

        # Strict Gating
        is_cpf = 0
        is_outgroup_splb = 0
        confidence_tier = "Rejected_No_Domain"

        if cpf_hmms:
            is_cpf = 1
            if best_cpf_score >= 50.0 and best_cpf_cov >= 0.50:
                confidence_tier = "High"
            elif best_cpf_score >= 25.0:
                confidence_tier = "Moderate"
            else:
                confidence_tier = "Low"
        elif splb_hmms and not cpf_hmms:
            is_outgroup_splb = 1
            assigned_broad = "SplB"
            confidence_tier = "SplB_Radical_SAM_Outgroup"
        else:
            confidence_tier = "Rejected_No_CPF_Domain"

        seq_len = len(seq)
        cys_count = seq.count("C")
        trp_count = seq.count("W")
        is_short_form = 1 if (seq_len < 420 and is_cpf == 1) else 0

        mag = meta.get("mag", "")
        sample = mag.split("__")[0] if "__" in mag else ""
        sample_info = meta_map.get(sample, {})
        tax_info = tax_map.get(mag, {})

        final_candidates.append({
            "candidate_id": qid,
            "mag": mag,
            "sample": sample,
            "compartment": sample_info.get("compartment", "Unknown"),
            "domain": tax_info.get("domain", "Unknown"),
            "phylum": tax_info.get("phylum", "Unknown"),
            "class": tax_info.get("class", "Unknown"),
            "order": tax_info.get("order", "Unknown"),
            "family": tax_info.get("family", "Unknown"),
            "genus": tax_info.get("genus", "Unknown"),
            "species": tax_info.get("species", "Unknown"),
            "length_aa": seq_len,
            "partial": meta.get("partial", "unknown"),
            "is_complete": 1 if meta.get("is_complete", False) else 0,
            "start_bp": meta.get("start_bp", 0),
            "stop_bp": meta.get("stop_bp", 0),
            "strand": meta.get("strand", 1),
            "top_reference_id": top_ref,
            "assigned_broad_family": assigned_broad,
            "confidence_tier": confidence_tier,
            "is_confirmed_cpf": is_cpf,
            "is_outgroup_splb": is_outgroup_splb,
            "is_short_form": is_short_form,
            "ref_pident": top_hit["pident"],
            "ref_evalue": top_hit["evalue"],
            "ref_bitscore": top_hit["bitscore"],
            "ref_qcovhsp": top_hit["qcovhsp"],
            "matched_hmms": ",".join(matched_hmms),
            "best_cpf_hmm_score": best_cpf_score,
            "best_cpf_hmm_eval": best_cpf_eval,
            "best_cpf_hmm_coverage": round(best_cpf_cov, 3),
            "cys_count": cys_count,
            "trp_count": trp_count,
            "sequence": seq,
        })

    # Write f1_candidates.tsv
    cand_tsv = RES_DIR / "f1_candidates.tsv"
    fieldnames = [
        "candidate_id", "mag", "sample", "compartment", "domain", "phylum", "class", "order",
        "family", "genus", "species", "length_aa", "partial", "is_complete", "start_bp", "stop_bp",
        "strand", "top_reference_id", "assigned_broad_family", "confidence_tier",
        "is_confirmed_cpf", "is_outgroup_splb", "is_short_form",
        "ref_pident", "ref_evalue", "ref_bitscore", "ref_qcovhsp",
        "matched_hmms", "best_cpf_hmm_score", "best_cpf_hmm_eval", "best_cpf_hmm_coverage",
        "cys_count", "trp_count"
    ]
    with open(cand_tsv, "w") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for c in final_candidates:
            writer.writerow(c)
    print(f"\nWrote {len(final_candidates)} candidates to {cand_tsv.name}.")

    # Write complete sequences FASTA
    complete_faa = RES_DIR / "f1_complete_candidates.faa"
    complete_count = 0
    with open(complete_faa, "w") as f:
        for c in final_candidates:
            if c["is_confirmed_cpf"] == 1 and c["is_complete"] == 1:
                complete_count += 1
                f.write(f">{c['candidate_id']} family={c['assigned_broad_family']} tier={c['confidence_tier']} len={c['length_aa']} phylum={c['phylum']}\n{c['sequence']}\n")
    print(f"Wrote {complete_count} complete confirmed CPF candidate sequences to {complete_faa.name}.")

    # Summary Table
    print("\n--- Step 6: Generating Summary Table ---")
    fam_counts = Counter()
    fam_complete = Counter()
    fam_phyla = defaultdict(set)
    fam_comps = defaultdict(Counter)

    for c in final_candidates:
        if c["is_confirmed_cpf"] == 1:
            fam = c["assigned_broad_family"]
            fam_counts[fam] += 1
            if c["is_complete"] == 1:
                fam_complete[fam] += 1
            fam_phyla[fam].add(c["phylum"])
            fam_comps[fam][c["compartment"]] += 1

    summary_tsv = RES_DIR / "f1_subfamily_summary.tsv"
    with open(summary_tsv, "w") as f:
        f.write("broad_family\ttotal_candidates\tcomplete_loci\tpct_complete\tn_phyla\ttop_phyla\tsurface_count\tdeep_count\tplantroot_count\n")
        for fam, tot in sorted(fam_counts.items(), key=lambda x: -x[1]):
            comp_cnt = fam_complete[fam]
            pct = 100.0 * comp_cnt / tot if tot else 0.0
            phyla_str = ",".join(sorted(list(fam_phyla[fam]))[:4])
            c_dict = fam_comps[fam]
            f.write(f"{fam}\t{tot}\t{comp_cnt}\t{pct:.1f}%\t{len(fam_phyla[fam])}\t{phyla_str}\t{c_dict['Surface']}\t{c_dict['Deep']}\t{c_dict['PlantRoot']}\n")
            print(f"  {fam:<25} total={tot:<4} complete={comp_cnt:<4} ({pct:.1f}%) phyla={len(fam_phyla[fam])} Surface={c_dict['Surface']} Deep={c_dict['Deep']}")

    n_confirmed = sum(1 for c in final_candidates if c["is_confirmed_cpf"] == 1)
    n_high = sum(1 for c in final_candidates if c["confidence_tier"] == "High")
    n_short = sum(1 for c in final_candidates if c["is_short_form"] == 1)
    print(f"\nPilot Summary across 256 representatives:")
    print(f"  Confirmed CPF loci (passing HMM domain gate): {n_confirmed} ({complete_count} complete, {n_high} high-confidence)")
    print(f"  Short photolyase candidates (<420 aa): {n_short}")
    print(f"=== F1 Screen Pipeline Complete in {time.time() - t0:.1f}s ===")


if __name__ == "__main__":
    main()
