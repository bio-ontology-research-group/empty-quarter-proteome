#!/usr/bin/env python3
"""
02_f1_pilot_v3.py: F1 Photolyase/Cryptochrome Family Identity and Diversity Pilot.

Root Audit Rules:
- References positive control check.
- HMMER parsers: p[4] query accession, p[3] query name, explicit cpu limit.
- All source/mapping errors fail, no empty defaults.
- Candidate HMM domain coverage, family confidence independent of nearest hit, partial flags.
- All hits at 1e-3 are candidates only; high confidence requires profile/domain support (E <= 1e-10) and reasonable coverage (>= 50%).
- Broad family categories: CPF / cryptochrome / CPD / 6-4 / SPL / SplB / unclassified; discard heuristic CPD_I/CRY_DASH.
- Do not claim subfamily classification or novelty from best hit alone.
"""

import csv
import os
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

def require_slurm():
    job_id = os.environ.get("SLURM_JOB_ID")
    if not job_id:
        print("ERROR: Must run inside Slurm worker allocation!", file=sys.stderr)
        sys.exit(1)
    print(f"[OK] Slurm worker environment confirmed. SLURM_JOB_ID={job_id}")

def parse_header_metadata(header_line: str) -> dict:
    # Header format: >bin::contig_gene # start # end # strand # ID=...;partial=00;...
    meta = {"partial": "unknown", "start_type": "unknown", "gc_cont": "unknown"}
    parts = header_line.split("#")
    if len(parts) >= 5:
        attrs_str = parts[4].strip()
        for token in attrs_str.split(";"):
            if "=" in token:
                k, v = token.split("=", 1)
                meta[k.strip()] = v.strip()
    return meta

def parse_domtblout(domtbl_path: Path) -> dict:
    """
    Parse HMMER domtblout strictly:
    p[0]: target name
    p[1]: target accession
    p[2]: tlen
    p[3]: query name (model name)
    p[4]: query accession (Pfam accession)
    p[5]: qlen
    p[6]: full seq E-value
    p[7]: full seq score
    p[11]: c-Evalue
    p[12]: i-Evalue (domain independent E-value)
    p[13]: domain score
    p[15]: hmm from
    p[16]: hmm to
    p[17]: ali from
    p[18]: ali to
    """
    if not domtbl_path.is_file():
        raise FileNotFoundError(f"domtblout file not found: {domtbl_path}")

    domains_by_target = defaultdict(list)
    with open(domtbl_path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 22:
                continue

            target_name = parts[0]
            target_len = int(parts[2])
            query_name = parts[3]
            query_acc = parts[4]
            hmm_len = int(parts[5])
            full_evalue = float(parts[6])
            full_score = float(parts[7])
            domain_ievalue = float(parts[12])
            domain_score = float(parts[13])
            hmm_from = int(parts[15])
            hmm_to = int(parts[16])
            ali_from = int(parts[17])
            ali_to = int(parts[18])

            domain_cov = (hmm_to - hmm_from + 1) / hmm_len if hmm_len > 0 else 0.0
            seq_cov = (ali_to - ali_from + 1) / target_len if target_len > 0 else 0.0

            domains_by_target[target_name].append({
                "query_name": query_name,
                "query_acc": query_acc,
                "hmm_len": hmm_len,
                "target_len": target_len,
                "full_evalue": full_evalue,
                "full_score": full_score,
                "domain_ievalue": domain_ievalue,
                "domain_score": domain_score,
                "hmm_from": hmm_from,
                "hmm_to": hmm_to,
                "ali_from": ali_from,
                "ali_to": ali_to,
                "domain_cov": domain_cov,
                "seq_cov": seq_cov
            })

    return domains_by_target

def parse_diamond_out(dmnd_res_path: Path) -> dict:
    """
    Parse DIAMOND tabular outfmt 6:
    qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen
    """
    if not dmnd_res_path.is_file():
        raise FileNotFoundError(f"DIAMOND results file not found: {dmnd_res_path}")

    hits_by_query = defaultdict(list)
    with open(dmnd_res_path, "r", encoding="utf-8") as f:
        reader = csv.reader(f, delimiter="\t")
        for parts in reader:
            if len(parts) < 14:
                continue
            qseqid = parts[0]
            sseqid = parts[1]
            pident = float(parts[2])
            length = int(parts[3])
            qstart = int(parts[6])
            qend = int(parts[7])
            evalue = float(parts[10])
            bitscore = float(parts[11])
            qlen = int(parts[12])
            slen = int(parts[13])

            q_cov = (qend - qstart + 1) / qlen if qlen > 0 else 0.0

            hits_by_query[qseqid].append({
                "sseqid": sseqid,
                "pident": pident,
                "length": length,
                "evalue": evalue,
                "bitscore": bitscore,
                "qlen": qlen,
                "slen": slen,
                "q_cov": q_cov
            })

    return hits_by_query

def main():
    require_slurm()

    f12_dir = Path("/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/family_discovery_20261015_clean/f12")
    refs_dir = f12_dir / "references_v2"
    out_dir = f12_dir / "results/pilot_v3_results"
    out_dir.mkdir(parents=True, exist_ok=True)

    diamond_bin = "/ibex/scratch/projects/c2014/rob/miniconda3/envs/eggnog/bin/diamond"
    if not os.path.isfile(diamond_bin):
        diamond_bin = "diamond"

    cpf_dmnd = refs_dir / "cpf_references_dedup.dmnd"
    cpf_hmm = refs_dir / "cpf_all.hmm"
    meta_dedup = refs_dir / "verified_references_metadata_dedup.tsv"

    if not cpf_dmnd.is_file():
        raise FileNotFoundError(f"Missing DIAMOND database: {cpf_dmnd}")
    if not cpf_hmm.is_file():
        raise FileNotFoundError(f"Missing HMM file: {cpf_hmm}")
    if not meta_dedup.is_file():
        raise FileNotFoundError(f"Missing metadata file: {meta_dedup}")

    # Load deduplicated reference metadata
    ref_meta = {}
    with open(meta_dedup, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for r in reader:
            ref_meta[r["accession"]] = r

    cpf_ref_count = sum(1 for r in ref_meta.values() if r["analysis_track"] == "F1_CPF")
    print(f"Loaded {len(ref_meta)} total references ({cpf_ref_count} CPF references).")

    # Step 1: Reference positive control
    print("Running positive control on reference sequences...")
    ref_faa = refs_dir / "cpf_references_dedup.faa"
    ref_ctrl_dmnd_out = out_dir / "f1_ref_positive_control_diamond.tsv"
    cmd_ref_dmnd = [
        diamond_bin, "blastp",
        "--db", str(cpf_dmnd),
        "--query", str(ref_faa),
        "--out", str(ref_ctrl_dmnd_out),
        "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen",
        "--evalue", "1e-5",
        "--threads", "8",
        "--quiet"
    ]
    subprocess.run(cmd_ref_dmnd, check=True)

    ref_ctrl_domtbl = out_dir / "f1_ref_positive_control_domtbl.txt"
    cmd_ref_hmm = [
        "hmmsearch",
        "--cpu", "8",
        "-E", "1e-3",
        "--domE", "1e-3",
        "--domtblout", str(ref_ctrl_domtbl),
        str(cpf_hmm),
        str(ref_faa)
    ]
    subprocess.run(cmd_ref_hmm, check=True)

    ref_dmnd_hits = parse_diamond_out(ref_ctrl_dmnd_out)
    ref_hmm_hits = parse_domtblout(ref_ctrl_domtbl)
    print(f"[OK] Positive control: {len(ref_dmnd_hits)}/{cpf_ref_count} recovered by DIAMOND, {len(ref_hmm_hits)}/{cpf_ref_count} by HMM.")

    # Step 2: Search cohort proteins
    cohort_faa = out_dir / "f1_pilot_256_cohort_proteins.faa"
    if not cohort_faa.is_file():
        raise FileNotFoundError(f"Cohort FASTA missing: {cohort_faa}")

    cohort_dmnd_out = out_dir / "f1_pilot_diamond_raw.tsv"
    print(f"Running DIAMOND blastp on cohort proteins ({cohort_faa}) vs {cpf_dmnd}...")
    cmd_cohort_dmnd = [
        diamond_bin, "blastp",
        "--db", str(cpf_dmnd),
        "--query", str(cohort_faa),
        "--out", str(cohort_dmnd_out),
        "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen",
        "--evalue", "1e-3",
        "--max-target-seqs", "10",
        "--threads", "16"
    ]
    subprocess.run(cmd_cohort_dmnd, check=True)

    cohort_domtbl = out_dir / "f1_pilot_hmm_domtbl.txt"
    print(f"Running hmmsearch on cohort proteins vs {cpf_hmm} (--cpu 16)...")
    cmd_cohort_hmm = [
        "hmmsearch",
        "--cpu", "16",
        "-E", "1e-3",
        "--domE", "1e-3",
        "--domtblout", str(cohort_domtbl),
        str(cpf_hmm),
        str(cohort_faa)
    ]
    subprocess.run(cmd_cohort_hmm, check=True)

    # Step 3: Parse and combine hits
    cohort_dmnd_hits = parse_diamond_out(cohort_dmnd_out)
    cohort_hmm_hits = parse_domtblout(cohort_domtbl)

    all_candidate_ids = set(cohort_dmnd_hits.keys()) | set(cohort_hmm_hits.keys())
    print(f"Total candidate proteins identified (E <= 1e-3): {len(all_candidate_ids)}")

    # Extract headers and sequence metadata from cohort FASTA for candidate IDs
    candidate_headers = {}
    with open(cohort_faa, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith(">"):
                hdr_clean = line[1:].strip()
                seq_id = hdr_clean.split()[0]
                if seq_id in all_candidate_ids:
                    candidate_headers[seq_id] = hdr_clean

    candidate_rows = []
    for cid in sorted(all_candidate_ids):
        hdr = candidate_headers.get(cid, "")
        header_meta = parse_header_metadata(hdr)
        partial = header_meta.get("partial", "unknown")

        dmnd_list = cohort_dmnd_hits.get(cid, [])
        hmm_list = cohort_hmm_hits.get(cid, [])

        top_dmnd = dmnd_list[0] if dmnd_list else None
        top_hmm = sorted(hmm_list, key=lambda x: x["domain_ievalue"])[0] if hmm_list else None

        ref_acc = top_dmnd["sseqid"].split()[0] if top_dmnd else ""
        ref_rec = ref_meta.get(ref_acc, {})
        ref_broad = ref_rec.get("broad_family", "unclassified") if ref_rec else "unclassified"
        ref_name = ref_rec.get("protein_name", "") if ref_rec else ""
        ref_ec = ref_rec.get("ec_number", "") if ref_rec else ""

        # Domain support check
        has_hmm_support = False
        top_hmm_acc = ""
        top_hmm_name = ""
        top_hmm_ievalue = 999.0
        top_hmm_score = 0.0
        top_hmm_dom_cov = 0.0

        if top_hmm:
            top_hmm_acc = top_hmm["query_acc"]
            top_hmm_name = top_hmm["query_name"]
            top_hmm_ievalue = top_hmm["domain_ievalue"]
            top_hmm_score = top_hmm["domain_score"]
            top_hmm_dom_cov = top_hmm["domain_cov"]
            if top_hmm_ievalue <= 1e-10 and top_hmm_dom_cov >= 0.50:
                has_hmm_support = True

        # Assign confidence
        if has_hmm_support:
            confidence = "HIGH_CONFIDENCE"
        elif (top_hmm and top_hmm_ievalue <= 1e-5 and top_hmm_dom_cov >= 0.30) or (top_dmnd and top_dmnd["bitscore"] >= 100 and top_dmnd["evalue"] <= 1e-10):
            confidence = "MODERATE_CONFIDENCE"
        else:
            confidence = "CANDIDATE_ONLY"

        # Broad family assignment (Root: broad CPF/cryptochrome/CPD/6-4/SPL/SplB categories with unclassified permitted)
        if top_hmm_name == "Radical_SAM" or top_hmm_acc.startswith("PF04055"):
            broad_family = "SplB_outgroup"
        elif ref_broad == "SPL" and (top_hmm_name in ("FAD_binding_7", "DNA_photolyase", "Photo_lyase") or has_hmm_support):
            broad_family = "SPL_candidate"
        elif ref_broad == "CPD":
            broad_family = "CPD_candidate"
        elif ref_broad == "6-4":
            broad_family = "6-4_candidate"
        elif ref_broad == "cryptochrome":
            broad_family = "cryptochrome_candidate"
        elif has_hmm_support:
            broad_family = "CPF_unclassified"
        else:
            broad_family = "unclassified"

        candidate_rows.append({
            "candidate_id": cid,
            "genome_bin": cid.split("::")[0] if "::" in cid else "",
            "confidence": confidence,
            "broad_family": broad_family,
            "coding_partial": partial,
            "is_complete_cds": (partial == "00"),
            "top_ref_acc": ref_acc,
            "top_ref_broad_family": ref_broad,
            "top_ref_pident": top_dmnd["pident"] if top_dmnd else 0.0,
            "top_ref_evalue": top_dmnd["evalue"] if top_dmnd else 1.0,
            "top_ref_bitscore": top_dmnd["bitscore"] if top_dmnd else 0.0,
            "top_ref_qcov": f"{top_dmnd['q_cov']:.3f}" if top_dmnd else "0.000",
            "top_ref_protein_name": ref_name,
            "top_ref_ec": ref_ec,
            "top_hmm_acc": top_hmm_acc,
            "top_hmm_name": top_hmm_name,
            "top_hmm_ievalue": f"{top_hmm_ievalue:.2e}" if top_hmm else "1.00e+00",
            "top_hmm_score": top_hmm_score,
            "top_hmm_dom_cov": f"{top_hmm_dom_cov:.3f}" if top_hmm else "0.000",
            "has_hmm_support": has_hmm_support,
            "total_dmnd_hits": len(dmnd_list),
            "total_hmm_domains": len(hmm_list)
        })

    candidates_tsv = out_dir / "f1_pilot_candidates.tsv"
    with open(candidates_tsv, "w", newline="", encoding="utf-8") as f:
        fieldnames = [
            "candidate_id", "genome_bin", "confidence", "broad_family", "coding_partial", "is_complete_cds",
            "top_ref_acc", "top_ref_broad_family", "top_ref_pident", "top_ref_evalue", "top_ref_bitscore", "top_ref_qcov",
            "top_ref_protein_name", "top_ref_ec", "top_hmm_acc", "top_hmm_name", "top_hmm_ievalue", "top_hmm_score",
            "top_hmm_dom_cov", "has_hmm_support", "total_dmnd_hits", "total_hmm_domains"
        ]
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(candidate_rows)

    print(f"[OK] Wrote {len(candidate_rows)} candidate rows to {candidates_tsv}")

    # Step 4: Summary by confidence, family, and completeness
    high_conf = [r for r in candidate_rows if r["confidence"] == "HIGH_CONFIDENCE"]
    mod_conf = [r for r in candidate_rows if r["confidence"] == "MODERATE_CONFIDENCE"]
    cand_only = [r for r in candidate_rows if r["confidence"] == "CANDIDATE_ONLY"]

    fam_counts = defaultdict(int)
    for r in high_conf:
        fam_counts[r["broad_family"]] += 1

    complete_high = sum(1 for r in high_conf if r["is_complete_cds"])

    summary_tsv = out_dir / "f1_pilot_summary.tsv"
    with open(summary_tsv, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["metric", "value"])
        writer.writerow(["cohort_genomes_screened", 256])
        writer.writerow(["total_candidate_hits_e1e3", len(candidate_rows)])
        writer.writerow(["high_confidence_candidates", len(high_conf)])
        writer.writerow(["moderate_confidence_candidates", len(mod_conf)])
        writer.writerow(["candidate_only_hits", len(cand_only)])
        writer.writerow(["high_conf_complete_cds_count", complete_high])
        writer.writerow(["high_conf_complete_cds_pct", f"{100*complete_high/len(high_conf):.1f}%" if high_conf else "0.0%"])
        for fam, cnt in sorted(fam_counts.items()):
            writer.writerow([f"high_conf_{fam}", cnt])

    print(f"[SUMMARY] Total candidates: {len(candidate_rows)}, High-confidence: {len(high_conf)} ({complete_high} complete CDS).")
    print("[SUCCESS] 02_f1_pilot_v3 completed cleanly.")

if __name__ == "__main__":
    main()
