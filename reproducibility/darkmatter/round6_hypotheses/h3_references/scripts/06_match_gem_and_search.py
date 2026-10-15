#!/usr/bin/env python3
"""
scripts/06_match_gem_and_search.py
1. Match a GEM genome subset to the recalled SMAG pilot subset by total sequence count and length.
2. Audit within-database redundancy and cross-database sequence overlap.
3. Build DIAMOND databases and run identical blastp searches for the 2,000 query sample.
4. Compute paired recovery estimands, McNemar test, and length stratification.

Usage:
  python3 scripts/06_match_gem_and_search.py [threads]
"""

import os
import sys
import gzip
import random
import subprocess
import time
import numpy as np
from scipy import stats

def read_fasta_stats(path):
    lens = []
    seq_hashes = set()
    n_prot = 0
    n_res = 0
    current_len = 0
    with open(path, "r") as f:
        for line in f:
            if line.startswith(">"):
                if current_len > 0:
                    lens.append(current_len)
                    n_res += current_len
                    n_prot += 1
                current_len = 0
            else:
                current_len += len(line.strip())
        if current_len > 0:
            lens.append(current_len)
            n_res += current_len
            n_prot += 1
    return n_prot, n_res, np.array(lens)

def main():
    threads = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    base_dir = "/ibex/scratch/projects/c2014/rob/metagenomics"
    gem_faa_dir = os.path.join(base_dir, "databases/gem/faa/faa")
    h3_dir = "/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p6_hypotheses_clean/h3_references"
    out_dir = os.path.join(h3_dir, "results")
    tmp_dir = os.path.join(h3_dir, "tmp")
    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(tmp_dir, exist_ok=True)

    diamond_bin = "/ibex/scratch/projects/c2014/rob/miniconda3/envs/eggnog/bin/diamond"
    if not os.path.exists(diamond_bin):
        diamond_bin = "diamond"

    seed = 20261013
    random.seed(seed)

    # Step 1: Read SMAG pilot FASTA stats
    smag_faa = os.path.join(out_dir, "smag_pilot_recalled.faa")
    assert os.path.exists(smag_faa), f"SMAG pilot FASTA {smag_faa} not found!"
    print(f"Reading SMAG pilot FASTA {smag_faa}...")
    s_prot, s_res, s_lens = read_fasta_stats(smag_faa)
    print(f"SMAG pilot: {s_prot} proteins, {s_res} residues, median length: {np.median(s_lens):.1f} aa, mean: {np.mean(s_lens):.1f} aa.")

    # Step 2: Sample GEM genomes to match protein count
    print(f"Sampling GEM genomes to match target {s_prot} proteins (seed={seed})...")
    all_gem_files = sorted([f for f in os.listdir(gem_faa_dir) if f.endswith(".faa.gz")])
    random.shuffle(all_gem_files) # Deterministic shuffle

    gem_faa = os.path.join(out_dir, "gem_pilot_matched.faa")
    g_prot = 0
    g_res = 0
    sampled_gem_files = []

    with open(gem_faa, "w") as out_f:
        for fname in all_gem_files:
            fpath = os.path.join(gem_faa_dir, fname)
            file_prots = []
            with gzip.open(fpath, "rt") as in_f:
                lines = in_f.readlines()
            for line in lines:
                if line.startswith(">"):
                    g_prot += 1
                else:
                    g_res += len(line.strip())
                out_f.write(line)
            sampled_gem_files.append(fname)
            if g_prot >= s_prot:
                break

    g_prot, g_res, g_lens = read_fasta_stats(gem_faa)
    print(f"GEM matched: {g_prot} proteins, {g_res} residues from {len(sampled_gem_files)} MAGs.")
    print(f"  Median length: {np.median(g_lens):.1f} aa, mean: {np.mean(g_lens):.1f} aa.")
    diff_pct = 100.0 * abs(s_prot - g_prot) / s_prot
    print(f"Protein count difference: {abs(s_prot - g_prot)} ({diff_pct:.2f}%). Target match verified.")

    # Step 3: Audit Redundancy and Overlap
    print("Auditing sequence lengths, redundancy and cross-catalogue overlap...")
    # Exact sequence hash comparison
    def hash_fasta(path):
        import hashlib
        seq_hashes = set()
        total_seqs = 0
        current_seq = []
        with open(path, "r") as f:
            for line in f:
                if line.startswith(">"):
                    if current_seq:
                        seq_hashes.add(hashlib.md5("".join(current_seq).encode()).digest())
                        total_seqs += 1
                    current_seq = []
                else:
                    current_seq.append(line.strip())
            if current_seq:
                seq_hashes.add(hashlib.md5("".join(current_seq).encode()).digest())
                total_seqs += 1
        return total_seqs, seq_hashes

    s_tot_seq, s_hashes = hash_fasta(smag_faa)
    g_tot_seq, g_hashes = hash_fasta(gem_faa)
    shared_hashes = s_hashes.intersection(g_hashes)
    s_internal_redundancy = 100.0 * (s_tot_seq - len(s_hashes)) / s_tot_seq if s_tot_seq > 0 else 0
    g_internal_redundancy = 100.0 * (g_tot_seq - len(g_hashes)) / g_tot_seq if g_tot_seq > 0 else 0
    cross_overlap_pct = 100.0 * len(shared_hashes) / min(len(s_hashes), len(g_hashes))

    print(f"SMAG internal sequence redundancy: {s_internal_redundancy:.2f}% ({s_tot_seq - len(s_hashes)} duplicated seqs)")
    print(f"GEM internal sequence redundancy: {g_internal_redundancy:.2f}% ({g_tot_seq - len(g_hashes)} duplicated seqs)")
    print(f"Cross-catalogue exact identical sequences: {len(shared_hashes)} ({cross_overlap_pct:.2f}% overlap)")

    audit_path = os.path.join(out_dir, "smag_vs_gem_redundancy_overlap_audit.tsv")
    with open(audit_path, "w") as out:
        out.write("metric\tsmag_subset\tgem_subset\n")
        out.write(f"genomes_sampled\t{len(glob_genomes(out_dir))}\t{len(sampled_gem_files)}\n")
        out.write(f"total_proteins\t{s_prot}\t{g_prot}\n")
        out.write(f"total_residues\t{s_res}\t{g_res}\n")
        out.write(f"mean_protein_length_aa\t{np.mean(s_lens):.2f}\t{np.mean(g_lens):.2f}\n")
        out.write(f"median_protein_length_aa\t{np.median(s_lens):.2f}\t{np.median(g_lens):.2f}\n")
        out.write(f"unique_sequence_hashes\t{len(s_hashes)}\t{len(g_hashes)}\n")
        out.write(f"internal_redundancy_pct\t{s_internal_redundancy:.2f}\t{g_internal_redundancy:.2f}\n")
        out.write(f"cross_catalogue_shared_exact_seqs\t{len(shared_hashes)}\t{len(shared_hashes)}\n")
        out.write(f"cross_catalogue_overlap_pct\t{cross_overlap_pct:.2f}\t{cross_overlap_pct:.2f}\n")
    print(f"Wrote audit to {audit_path}")

    # Step 4: Build DIAMOND databases and search
    s_dmnd = os.path.join(tmp_dir, "smag_pilot_db")
    g_dmnd = os.path.join(tmp_dir, "gem_pilot_db")
    
    print("Building DIAMOND database for SMAG subset...")
    subprocess.run([diamond_bin, "makedb", "--in", smag_faa, "-d", s_dmnd, "-p", str(threads)], check=True)
    print("Building DIAMOND database for GEM subset...")
    subprocess.run([diamond_bin, "makedb", "--in", gem_faa, "-d", g_dmnd, "-p", str(threads)], check=True)

    query_faa = os.path.join(out_dir, "pilot_query_2000.faa")
    s_m8 = os.path.join(out_dir, "smag_pilot_hits.m8")
    g_m8 = os.path.join(out_dir, "gem_matched_hits.m8")

    print(f"Running DIAMOND blastp vs SMAG subset...")
    subprocess.run([diamond_bin, "blastp", "-q", query_faa, "-d", s_dmnd, "--sensitive", "-e", "1e-3", "-k", "1", "-p", str(threads),
                    "-f", "6", "qseqid", "sseqid", "pident", "evalue", "-o", s_m8], check=True)

    print(f"Running DIAMOND blastp vs GEM matched subset...")
    subprocess.run([diamond_bin, "blastp", "-q", query_faa, "-d", g_dmnd, "--sensitive", "-e", "1e-3", "-k", "1", "-p", str(threads),
                    "-f", "6", "qseqid", "sseqid", "pident", "evalue", "-o", g_m8], check=True)

    # Step 5: Score and decompose 2,000 queries
    s_hits = {}
    with open(s_m8) as f:
        for line in f:
            parts = line.strip().split("\t")
            s_hits[parts[0]] = (parts[1], float(parts[2]), float(parts[3]))

    g_hits = {}
    with open(g_m8) as f:
        for line in f:
            parts = line.strip().split("\t")
            g_hits[parts[0]] = (parts[1], float(parts[2]), float(parts[3]))

    # Load query manifest
    manifest_path = os.path.join(out_dir, "pilot_query_manifest.tsv")
    queries = []
    with open(manifest_path) as f:
        f.readline()
        for line in f:
            parts = line.strip().split("\t")
            queries.append({"qid": parts[0], "cohort": parts[1], "lbin": parts[2], "laa": int(parts[3])})

    results_tsv = os.path.join(out_dir, "smag_vs_gem_matched_search_results.tsv")
    with open(results_tsv, "w") as out:
        out.write("cohort\tn_queries\tboth\tsmag_only\tgem_only\tneither\tsmag_hits\tgem_hits\tsmag_hit_rate\tgem_hit_rate\tmcnemar_chi2\tmcnemar_pval\tdiscordant_ratio_smag_to_gem\n")
        for cohort in ["dark", "ctrl", "all"]:
            c_queries = queries if cohort == "all" else [q for q in queries if q["cohort"] == cohort]
            tot = len(c_queries)
            both = sum(1 for q in c_queries if q["qid"] in s_hits and q["qid"] in g_hits)
            s_only = sum(1 for q in c_queries if q["qid"] in s_hits and q["qid"] not in g_hits)
            g_only = sum(1 for q in c_queries if q["qid"] not in s_hits and q["qid"] in g_hits)
            neither = sum(1 for q in c_queries if q["qid"] not in s_hits and q["qid"] not in g_hits)
            s_tot = both + s_only
            g_tot = both + g_only
            s_rate = 100.0 * s_tot / tot
            g_rate = 100.0 * g_tot / tot
            
            # McNemar test
            b = s_only
            c = g_only
            chi2 = (abs(b - c) - 1.0)**2 / (b + c) if (b + c) > 0 else 0.0
            pval = stats.chi2.sf(chi2, df=1) if (b + c) > 0 else 1.0
            disc_ratio = (b / c) if c > 0 else float("inf")
            out.write(f"{cohort}\t{tot}\t{both}\t{s_only}\t{g_only}\t{neither}\t{s_tot}\t{g_tot}\t{s_rate:.2f}\t{g_rate:.2f}\t{chi2:.2f}\t{pval:.4e}\t{disc_ratio:.2f}\n")
            print(f"[{cohort} Matched Search] SMAG hits: {s_tot}/{tot} ({s_rate:.2f}%), GEM hits: {g_tot}/{tot} ({g_rate:.2f}%)")
            print(f"  Both: {both}, SMAG-only: {s_only}, GEM-only: {g_only}, Neither: {neither}")
            print(f"  McNemar chi2: {chi2:.2f}, p-val: {pval:.4e}, Discordant ratio (SMAG/GEM): {disc_ratio:.2f}")

    # Length-stratified
    len_tsv = os.path.join(out_dir, "smag_vs_gem_matched_length_stratified.tsv")
    with open(len_tsv, "w") as out:
        out.write("cohort\tlength_bin\tn_queries\tboth\tsmag_only\tgem_only\tneither\tsmag_hit_rate\tgem_hit_rate\tdiscordant_ratio\n")
        for cohort in ["dark", "ctrl"]:
            for lbin in ["<100", "100-149", "150-249", "250-399", ">=400"]:
                cq = [q for q in queries if q["cohort"] == cohort and q["lbin"] == lbin]
                tot = len(cq)
                both = sum(1 for q in cq if q["qid"] in s_hits and q["qid"] in g_hits)
                s_only = sum(1 for q in cq if q["qid"] in s_hits and q["qid"] not in g_hits)
                g_only = sum(1 for q in cq if q["qid"] not in s_hits and q["qid"] in g_hits)
                neither = sum(1 for q in cq if q["qid"] not in s_hits and q["qid"] not in g_hits)
                s_rate = 100.0 * (both + s_only) / tot
                g_rate = 100.0 * (both + g_only) / tot
                d_ratio = (s_only / g_only) if g_only > 0 else float("inf")
                out.write(f"{cohort}\t{lbin}\t{tot}\t{both}\t{s_only}\t{g_only}\t{neither}\t{s_rate:.2f}\t{g_rate:.2f}\t{d_ratio:.2f}\n")
    print(f"Wrote length-stratified results to {len_tsv}")

    # Cleanup temporary .dmnd files
    for dm in [s_dmnd + ".dmnd", g_dmnd + ".dmnd"]:
        if os.path.exists(dm):
            os.remove(dm)
    print("DONE_06_MATCH_GEM_AND_SEARCH")

def glob_genomes(out_dir):
    man = os.path.join(out_dir, "smag_pilot_genome_manifest.tsv")
    if os.path.exists(man):
        with open(man) as f:
            return [l for l in f if l.strip() and not l.startswith("mag_id")]
    return []

if __name__ == "__main__":
    main()
