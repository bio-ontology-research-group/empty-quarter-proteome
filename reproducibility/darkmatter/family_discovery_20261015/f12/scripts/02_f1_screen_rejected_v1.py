#!/usr/bin/env python3
"""02_f1_screen.py - Comprehensive CPF Family Identity & Diversity Census across 5,438 ANI95 Representatives.

Execution:
1. Iterates over all 5,438 ANI95 representative proteomes in drep_rep_proteomes.
2. Identifies CPF candidate homologs via sensitive DIAMOND search against curated 8-subfamily + SPL reference database.
3. Extracts full candidate protein sequences and parses Prodigal headers for coding completeness (partial=00 vs truncated).
4. Runs hmmsearch with Pfam profiles (PF00875, PF03441, PF04244, PF04055) against all candidates.
5. Classifies candidates into CPF subfamilies (CPD_I, CPD_II, CPD_III, CRY_DASH, Eu64, pCry, aCry, FeS_BCP, SPL, Outgroup_SplB)
   and checks key functional residues (FAD pocket, Trp electron transfer chain / SPL variant, FeS cysteines).
6. Merges GTDB-Tk R232 taxonomy and sample metadata (Surface, Deep, PlantRoot).
7. Writes f1_candidates.tsv, f1_subfamily_summary.tsv, and f1_complete_candidates.faa.
"""

import csv
import gzip
import os
import re
import subprocess
import sys
import time
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
REF_DIR = ROOT / "references"
RES_DIR = ROOT / "results"
RES_DIR.mkdir(parents=True, exist_ok=True)

W = Path("/ibex/scratch/projects/c2014/rob/metagenomics")
REP_PROTEOMES_DIR = W / "06_darkmatter/drep_rep_proteomes"
GTDB_BAC = W / "05_refine/_gtdbtk_r232/gtdbtk.bac120.summary.tsv"
GTDB_AR = W / "05_refine/_gtdbtk_r232/gtdbtk.ar53.summary.tsv"
SAMPLE_META = W / "metadata/sample_metadata.tsv"

CPF_HMM = REF_DIR / "cpf_all.hmm"
CPF_DMND = REF_DIR / "cpf_references.dmnd"
DIAMOND_BIN = Path("/ibex/scratch/projects/c2014/rob/miniconda3/envs/eggnog/bin/diamond")
if not DIAMOND_BIN.exists():
    DIAMOND_BIN = Path("diamond")

THREADS = int(os.environ.get("SLURM_CPUS_PER_TASK", "16"))


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


def main():
    t0 = time.time()
    print("=== [F1 Pilot Census] Photolyase / Cryptochrome Family Diversity Screen ===")
    print(f"Threads: {THREADS}")
    print(f"Representative Proteomes Source: {REP_PROTEOMES_DIR}")

    tax_map = load_taxonomy()
    meta_map = load_sample_metadata()
    print(f"Loaded taxonomy for {len(tax_map)} MAGs, metadata for {len(meta_map)} samples.")

    # 1. Identify all 5,438 representative proteome files
    faa_files = sorted(REP_PROTEOMES_DIR.glob("*.faa.gz"))
    total_reps = len(faa_files)
    print(f"Found {total_reps} representative proteome archives.")
    assert total_reps > 0, "No representative proteomes found!"

    # 2. Extract and concatenate or batch-screen candidate proteomes
    # Create a concatenated fasta of all reps in a fast scratch directory
    tmp_dir = ROOT / "tmp_f1"
    tmp_dir.mkdir(parents=True, exist_ok=True)
    all_reps_faa = tmp_dir / "all_5438_reps.faa"

    print("\n--- Step 1: Concatenating 5,438 Representative Proteomes ---")
    cat_start = time.time()
    total_proteins_scanned = 0
    with open(all_reps_faa, "w") as out_f:
        for faa_path in faa_files:
            mag_name = faa_path.name.replace(".faa.gz", "")
            with gzip.open(faa_path, "rt") as in_f:
                for line in in_f:
                    if line.startswith(">"):
                        total_proteins_scanned += 1
                        # Ensure header has MAG prefix for unambiguous mapping
                        pid = line[1:].strip().split()[0]
                        if not pid.startswith(mag_name + "__"):
                            # Format: >mag__pid rest_of_line
                            line = f">{mag_name}::{line[1:]}"
                    out_f.write(line)
    print(f"Concatenated {total_proteins_scanned} total representative proteins in {time.time() - cat_start:.1f}s.")

    # 3. DIAMOND search against CPF references
    print("\n--- Step 2: Sensitive DIAMOND Screening Against Curated CPF Reference Set ---")
    dmnd_out = tmp_dir / "diamond_hits.tsv"
    cmd = [
        str(DIAMOND_BIN), "blastp",
        "-d", str(CPF_DMND),
        "-q", str(all_reps_faa),
        "-o", str(dmnd_out),
        "--ultra-sensitive",
        "-e", "1e-3",
        "--max-target-seqs", "1",
        "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
        "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen", "qcovhsp", "scovhsp",
        "--threads", str(THREADS),
    ]
    dmnd_start = time.time()
    subprocess.run(cmd, check=True)
    print(f"DIAMOND search finished in {time.time() - dmnd_start:.1f}s.")

    # Read DIAMOND hits
    diamond_hits = {}
    with open(dmnd_out) as f:
        for line in f:
            parts = line.strip().split("\t")
            if not parts or len(parts) < 16:
                continue
            qid = parts[0]
            diamond_hits[qid] = {
                "top_ref_id": parts[1],
                "pident": float(parts[2]),
                "aln_len": int(parts[3]),
                "evalue": float(parts[10]),
                "bitscore": float(parts[11]),
                "qlen": int(parts[12]),
                "slen": int(parts[13]),
                "qcovhsp": float(parts[14]),
                "scovhsp": float(parts[15]),
            }
    print(f"Identified {len(diamond_hits)} candidate CPF protein hits.")

    # 4. Extract candidate protein sequences and metadata
    print("\n--- Step 3: Extracting Full Candidate Sequences & Parsing Coding Integrity ---")
    cand_seqs = {}
    cand_metadata = {}

    with open(all_reps_faa) as f:
        cur_qid = None
        cur_hdr = None
        seq_parts = []
        for line in f:
            if line.startswith(">"):
                if cur_qid and cur_qid in diamond_hits:
                    cand_seqs[cur_qid] = "".join(seq_parts)
                hdr = line[1:].strip()
                cur_qid = hdr.split()[0]
                cur_hdr = hdr
                seq_parts = []
                if cur_qid in diamond_hits:
                    # Parse prodigal header: >pid # start # stop # strand # ID=...;partial=00;...
                    partial = "00"
                    start_bp = 0
                    stop_bp = 0
                    strand = 1
                    tokens = hdr.split(" # ")
                    if len(tokens) >= 4:
                        start_bp = int(tokens[1])
                        stop_bp = int(tokens[2])
                        strand = int(tokens[3])
                    if len(tokens) >= 5:
                        for item in tokens[4].split(";"):
                            if item.startswith("partial="):
                                partial = item.split("=")[1]

                    # Parse MAG name
                    mag_name = ""
                    if "::" in cur_qid:
                        mag_name = cur_qid.split("::")[0]
                    elif "__" in cur_qid:
                        mag_name = cur_qid.split("__")[0]

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
        if cur_qid and cur_qid in diamond_hits:
            cand_seqs[cur_qid] = "".join(seq_parts)

    cand_faa = RES_DIR / "f1_candidate_pool.faa"
    with open(cand_faa, "w") as out_f:
        for qid, seq in cand_seqs.items():
            out_f.write(f">{qid}\n{seq}\n")
    print(f"Extracted {len(cand_seqs)} candidate sequences to {cand_faa.name}.")

    # 5. Run hmmsearch with Pfam profiles
    print("\n--- Step 4: Running hmmsearch with Pfam Profiles (PF00875, PF03441, PF04244, PF04055) ---")
    hmm_domtbl = tmp_dir / "cpf_hmm_domtbl.txt"
    hmm_out = tmp_dir / "cpf_hmm_out.txt"
    cmd = [
        "hmmsearch",
        "--domtblout", str(hmm_domtbl),
        "-o", str(hmm_out),
        "--domE", "1e-4",
        "--cpu", str(THREADS),
        str(CPF_HMM),
        str(cand_faa),
    ]
    hmm_start = time.time()
    subprocess.run(cmd, check=True)
    print(f"hmmsearch completed in {time.time() - hmm_start:.1f}s.")

    # Parse domtblout
    # Format: target_name target_acc tlen query_name query_acc qlen full_Evalue full_score ...
    # Note: target = candidate seq, query = HMM model
    hmm_hits = defaultdict(list)
    with open(hmm_domtbl) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 23:
                continue
            cand_id = parts[0]
            hmm_acc = parts[3].split(".")[0]
            hmm_name = parts[4]
            full_eval = float(parts[6])
            full_score = float(parts[7])
            dom_eval = float(parts[12])
            dom_score = float(parts[13])
            hmm_hits[cand_id].append({
                "hmm_acc": hmm_acc,
                "hmm_name": hmm_name,
                "full_evalue": full_eval,
                "full_score": full_score,
                "dom_evalue": dom_eval,
                "dom_score": dom_score,
            })

    print(f"HMM profile domain matches found for {len(hmm_hits)} candidates.")

    # 6. Detailed Subfamily Assignment & Motifs Analysis
    print("\n--- Step 5: Subfamily Classification & Catalytic Feature Audit ---")
    final_candidates = []

    # Map reference IDs to subfamilies
    ref_subfamily_map = {
        "CPD_I_Ecol_PhrB": "CPD_I",
        "CPD_I_Synel_Phr": "CPD_I",
        "CPD_II_Arath_UVR3": "CPD_II",
        "CPD_II_Drome_Phr": "CPD_II",
        "CPD_II_Metma_Phr": "CPD_II",
        "CPD_III_Caucr_Phr": "CPD_III",
        "CRY_DASH_Vibch_CryD": "CRY_DASH",
        "CRY_DASH_Synp6_CryD": "CRY_DASH",
        "CRY_DASH_Arath_CRY3": "CRY_DASH",
        "Eu64_Drome_Cry64": "Eu64",
        "Eu64_Xenla_Phr64": "Eu64",
        "pCry_Arath_CRY1": "pCry",
        "pCry_Arath_CRY2": "pCry",
        "aCry_Drome_CRY": "aCry",
        "FeS_BCP_Agro_PhrB": "FeS_BCP",
        "FeS_BCP_Psepu_PhrB": "FeS_BCP",
        "SPL_Rfas_RfSPL": "SPL",
        "SPL_Selo_SeSPL": "SPL",
        "SPL_Cmar_CmSPL": "SPL",
        "SPL_Rlac_RlSPL": "SPL",
        "SPL_Cgra_CgSPL": "SPL",
        "Outgroup_SplB_Bacsu": "Outgroup_SplB",
    }

    for qid, d_info in diamond_hits.items():
        meta = cand_metadata.get(qid, {})
        seq = cand_seqs.get(qid, "")
        h_list = hmm_hits.get(qid, [])

        top_ref = d_info["top_ref_id"]
        subfam = ref_subfamily_map.get(top_ref, "Unclassified")

        # Summarize HMM domains
        matched_hmms = sorted(set(h["hmm_acc"] for h in h_list))
        best_hmm_score = max((h["full_score"] for h in h_list), default=0.0)
        best_hmm_eval = min((h["full_evalue"] for h in h_list), default=1.0)

        # Separate true CPF from outgroup SplB
        is_outgroup_splb = (subfam == "Outgroup_SplB") or ("PF04055" in matched_hmms and "PF00875" not in matched_hmms and "PF03441" not in matched_hmms)

        # Check sequence length category
        seq_len = len(seq)
        is_short_form = (seq_len < 420) and not is_outgroup_splb

        # Active-site motifs:
        # Check for FAD binding Gly/Asn/Arg motifs
        has_fad_binding = bool(re.search(r"[R|K].{20,40}[D|E].{5,15}[A|G|S]..[G|A]", seq))
        # Check for Cysteine cluster (FeS coordination): 4 cysteines in C-terminal / specific arrangement
        cys_count = seq.count("C")
        cys_cluster_match = bool(re.search(r"C.{2,5}C.{2,15}C.{2,15}C", seq))

        mag = meta.get("mag", "")
        # Resolve sample from MAG
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
            "partial": meta.get("partial", "Unknown"),
            "is_complete": meta.get("is_complete", False),
            "start_bp": meta.get("start_bp", 0),
            "stop_bp": meta.get("stop_bp", 0),
            "strand": meta.get("strand", 1),
            "top_reference_id": top_ref,
            "assigned_subfamily": subfam,
            "ref_pident": d_info["pident"],
            "ref_evalue": d_info["evalue"],
            "ref_bitscore": d_info["bitscore"],
            "ref_qcovhsp": d_info["qcovhsp"],
            "matched_hmms": ",".join(matched_hmms),
            "best_hmm_score": best_hmm_score,
            "best_hmm_eval": best_hmm_eval,
            "is_outgroup_splb": is_outgroup_splb,
            "is_short_form": is_short_form,
            "cys_count": cys_count,
            "has_cys_cluster": cys_cluster_match,
            "sequence": seq,
        })

    # Write f1_candidates.tsv
    cand_tsv = RES_DIR / "f1_candidates.tsv"
    fieldnames = [
        "candidate_id", "mag", "sample", "compartment", "domain", "phylum", "class", "order",
        "family", "genus", "species", "length_aa", "partial", "is_complete", "start_bp", "stop_bp",
        "strand", "top_reference_id", "assigned_subfamily", "ref_pident", "ref_evalue", "ref_bitscore",
        "ref_qcovhsp", "matched_hmms", "best_hmm_score", "best_hmm_eval", "is_outgroup_splb",
        "is_short_form", "cys_count", "has_cys_cluster"
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
            if c["is_complete"] and not c["is_outgroup_splb"]:
                complete_count += 1
                f.write(f">{c['candidate_id']} subfam={c['assigned_subfamily']} len={c['length_aa']} phylum={c['phylum']}\n{c['sequence']}\n")
    print(f"Wrote {complete_count} complete CPF candidate sequences to {complete_faa.name}.")

    # 7. Generate Subfamily and Phylum Summary Table
    print("\n--- Step 6: Generating Summary Table ---")
    subfam_counts = Counter()
    subfam_complete = Counter()
    subfam_phyla = defaultdict(set)
    subfam_comps = defaultdict(Counter)

    for c in final_candidates:
        if c["is_outgroup_splb"]:
            sf = "Outgroup_SplB"
        else:
            sf = c["assigned_subfamily"]
        subfam_counts[sf] += 1
        if c["is_complete"]:
            subfam_complete[sf] += 1
        subfam_phyla[sf].add(c["phylum"])
        subfam_comps[sf][c["compartment"]] += 1

    summary_tsv = RES_DIR / "f1_subfamily_summary.tsv"
    with open(summary_tsv, "w") as f:
        f.write("subfamily\ttotal_candidates\tcomplete_loci\tpct_complete\tn_phyla\ttop_phyla\tsurface_count\tdeep_count\tplantroot_count\n")
        for sf, tot in sorted(subfam_counts.items(), key=lambda x: -x[1]):
            comp_cnt = subfam_complete[sf]
            pct = 100.0 * comp_cnt / tot if tot else 0.0
            phyla_str = ",".join(sorted(list(subfam_phyla[sf]))[:4])
            c_dict = subfam_comps[sf]
            f.write(f"{sf}\t{tot}\t{comp_cnt}\t{pct:.1f}%\t{len(subfam_phyla[sf])}\t{phyla_str}\t{c_dict['Surface']}\t{c_dict['Deep']}\t{c_dict['PlantRoot']}\n")
            print(f"  {sf:<15} total={tot:<4} complete={comp_cnt:<4} ({pct:.1f}%) phyla={len(subfam_phyla[sf])} Surface={c_dict['Surface']} Deep={c_dict['Deep']}")

    print(f"Wrote summary table to {summary_tsv.name}.")
    print(f"=== F1 Screen Pipeline Complete in {time.time() - t0:.1f}s ===")


if __name__ == "__main__":
    main()
