#!/usr/bin/env python3
"""
scripts/03_pilot_subsample_search.py
Execute the empirical search pilot to test and prove the Top-1 Subsampling Fallacy:
1. Deterministically sample 2,500 GEM genomes (~5% of GEM; ~6M proteins) with seed 20261013.
2. Build a DIAMOND database from these 2,500 genomes.
3. Search the 2,000 length-stratified query subset (1,000 dark + 1,000 control) against the subsampled database.
4. Compare empirical hits vs post-hoc filtering of the full GEM top-1 hit table (gem_hits.m8).
5. Quantify the exact masking bias (undercount rate of post-hoc filtering).

Usage:
  python3 scripts/03_pilot_subsample_search.py [threads]
"""

import os
import sys
import gzip
import random
import subprocess
import time

def main():
    threads = int(sys.argv[1]) if len(sys.argv) > 1 else 16
    base_dir = "/ibex/scratch/projects/c2014/rob/metagenomics"
    gem_faa_dir = os.path.join(base_dir, "databases/gem/faa/faa")
    spec_dir = os.path.join(base_dir, "06_darkmatter/specificity_clean/results")
    h3_dir = "/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p6_hypotheses_clean/h3_references"
    out_dir = os.path.join(h3_dir, "results")
    os.makedirs(out_dir, exist_ok=True)
    tmp_dir = os.path.join(h3_dir, "tmp")
    os.makedirs(tmp_dir, exist_ok=True)

    diamond_bin = "/ibex/scratch/projects/c2014/rob/miniconda3/envs/eggnog/bin/diamond"
    if not os.path.exists(diamond_bin):
        diamond_bin = "diamond"

    seed = 20261013
    random.seed(seed)

    print(f"--- Pilot Subsample Search & Top-1 Bias Proof (Threads: {threads}) ---")
    t0 = time.time()

    # Step 1: Enumerate GEM genomes and select 2,500
    print("Step 1: Listing GEM genomes...")
    all_gem_files = [f for f in os.listdir(gem_faa_dir) if f.endswith(".faa.gz")]
    all_gem_files.sort()
    print(f"Total GEM genome files available: {len(all_gem_files)}")
    assert len(all_gem_files) == 52515, f"Expected 52515 genomes, found {len(all_gem_files)}"

    n_sample_genomes = 2500
    sampled_files = random.sample(all_gem_files, n_sample_genomes)
    sampled_files.sort()
    print(f"Sampled {len(sampled_files)} GEM genomes (seed={seed}).")

    # Record sampled genome list
    sample_list_path = os.path.join(out_dir, "pilot_gem_2500_genomes.txt")
    with open(sample_list_path, "w") as out:
        for f in sampled_files:
            out.write(f + "\n")

    # Step 2: Concatenate sampled FASTA files and collect protein ID set
    print("Step 2: Concatenating sampled FASTA and indexing subject IDs...")
    sub_faa = os.path.join(tmp_dir, "gem_sub2500.faa")
    retained_subjects = set()
    n_proteins = 0
    n_residues = 0

    with open(sub_faa, "w") as out_f:
        for fname in sampled_files:
            fpath = os.path.join(gem_faa_dir, fname)
            with gzip.open(fpath, "rt") as in_f:
                for line in in_f:
                    if line.startswith(">"):
                        n_proteins += 1
                        # extract subject ID (first token after >)
                        sid = line[1:].strip().split()[0]
                        retained_subjects.add(sid)
                    else:
                        n_residues += len(line.strip())
                    out_f.write(line)

    print(f"Compiled subsample FASTA: {n_proteins} proteins, {n_residues} residues from {n_sample_genomes} genomes.")
    print(f"Retained subject IDs in memory: {len(retained_subjects)}")

    # Step 3: Build DIAMOND database
    sub_dmnd = os.path.join(tmp_dir, "gem_sub2500")
    print("Step 3: Building DIAMOND database with makedb...")
    cmd_makedb = [
        diamond_bin, "makedb",
        "--in", sub_faa,
        "-d", sub_dmnd,
        "-p", str(threads)
    ]
    t_make = time.time()
    subprocess.run(cmd_makedb, check=True)
    print(f"DIAMOND makedb finished in {time.time() - t_make:.2f}s.")

    # Step 4: Run empirical DIAMOND search of 2,000 query sample
    query_faa = os.path.join(out_dir, "pilot_query_2000.faa")
    empirical_m8 = os.path.join(out_dir, "pilot_sub2500_hits.m8")
    print(f"Step 4: Running DIAMOND blastp of {query_faa} vs {sub_dmnd}...")
    cmd_blastp = [
        diamond_bin, "blastp",
        "-q", query_faa,
        "-d", sub_dmnd,
        "--sensitive",
        "--evalue", "1e-3",
        "-k", "1",
        "-p", str(threads),
        "-f", "6", "qseqid", "sseqid", "pident", "evalue",
        "-o", empirical_m8
    ]
    t_search = time.time()
    subprocess.run(cmd_blastp, check=True)
    print(f"DIAMOND blastp finished in {time.time() - t_search:.2f}s.")

    # Step 5: Load empirical hits
    empirical_hits = {}
    with open(empirical_m8, "r") as f:
        for line in f:
            parts = line.strip().split("\t")
            empirical_hits[parts[0]] = (parts[1], float(parts[2]), float(parts[3]))
    print(f"Empirical search returned hits for {len(empirical_hits)} / 2000 queries.")

    # Step 6: Load full GEM hits and evaluate post-hoc filtering
    print("Step 6: Evaluating post-hoc filtering from full GEM hits (gem_hits.m8)...")
    full_gem_m8 = os.path.join(spec_dir, "gem_hits.m8")
    full_hits = {}
    with open(full_gem_m8, "r") as f:
        for line in f:
            parts = line.strip().split("\t")
            full_hits[parts[0]] = (parts[1], float(parts[2]), float(parts[3]))

    # Load query manifest
    manifest_path = os.path.join(out_dir, "pilot_query_manifest.tsv")
    queries = []
    with open(manifest_path, "r") as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split("\t")
            queries.append({
                "query_id": parts[0],
                "cohort": parts[1],
                "length_bin": parts[2],
                "length_aa": int(parts[3])
            })

    # Step 7: Compare Empirical Search vs Post-hoc Filter
    # Categories:
    # 1. both_hit: Hit in empirical search AND post-hoc filter (top hit was in subsample)
    # 2. masked_hit: Hit in empirical search BUT missed by post-hoc filter (alternative hit exists!)
    # 3. both_absent: No hit in either
    # 4. impossible_posthoc_only: Post-hoc hit but no empirical hit (should be 0)

    results_table = []
    for q in queries:
        qid = q["query_id"]
        has_emp = qid in empirical_hits
        
        # Post-hoc rule: does full GEM have hit, and is hit subject in retained_subjects?
        has_posthoc = False
        posthoc_sid = None
        if qid in full_hits:
            sid = full_hits[qid][0]
            if sid in retained_subjects:
                has_posthoc = True
                posthoc_sid = sid

        if has_emp and has_posthoc:
            cat = "retained_top_hit"
        elif has_emp and not has_posthoc:
            cat = "masked_alternative_hit"
        elif not has_emp and not has_posthoc:
            cat = "concordant_absence"
        else:
            cat = "inconsistent_posthoc_only"

        emp_s = empirical_hits[qid][0] if has_emp else "NONE"
        emp_pid = empirical_hits[qid][1] if has_emp else 0.0
        emp_eval = empirical_hits[qid][2] if has_emp else 1.0

        full_s = full_hits[qid][0] if qid in full_hits else "NONE"
        full_pid = full_hits[qid][1] if qid in full_hits else 0.0
        full_eval = full_hits[qid][2] if qid in full_hits else 1.0

        results_table.append({
            "query_id": qid,
            "cohort": q["cohort"],
            "length_bin": q["length_bin"],
            "length_aa": q["length_aa"],
            "has_empirical": has_emp,
            "has_posthoc": has_posthoc,
            "category": cat,
            "empirical_subject": emp_s,
            "empirical_pident": emp_pid,
            "empirical_evalue": emp_eval,
            "full_gem_subject": full_s,
            "full_gem_pident": full_pid,
            "full_gem_evalue": full_eval
        })

    # Write detailed comparison table
    comp_tsv = os.path.join(out_dir, "pilot_query_comparison_details.tsv")
    with open(comp_tsv, "w") as out:
        cols = ["query_id", "cohort", "length_bin", "length_aa", "has_empirical", "has_posthoc", "category",
                "empirical_subject", "empirical_pident", "empirical_evalue", "full_gem_subject", "full_gem_pident", "full_gem_evalue"]
        out.write("\t".join(cols) + "\n")
        for row in results_table:
            out.write("\t".join(str(row[c]) for c in cols) + "\n")
    print(f"Wrote query comparison details to {comp_tsv}")

    # Summary statistics by cohort and length
    summary_tsv = os.path.join(out_dir, "pilot_top1_bias_summary.tsv")
    with open(summary_tsv, "w") as out:
        out.write("cohort\tlength_bin\tn_queries\tempirical_hits\tposthoc_hits\tmasked_hits\tposthoc_undercount_pct\tempirical_hit_rate\tposthoc_hit_rate\n")
        
        for cohort in ["dark", "ctrl", "all"]:
            c_rows = results_table if cohort == "all" else [r for r in results_table if r["cohort"] == cohort]
            
            # Overall for cohort
            n_tot = len(c_rows)
            n_emp = sum(1 for r in c_rows if r["has_empirical"])
            n_post = sum(1 for r in c_rows if r["has_posthoc"])
            n_masked = sum(1 for r in c_rows if r["category"] == "masked_alternative_hit")
            undercount_pct = (100.0 * n_masked / n_emp) if n_emp > 0 else 0.0
            emp_rate = 100.0 * n_emp / n_tot
            post_rate = 100.0 * n_post / n_tot
            out.write(f"{cohort}\tall\t{n_tot}\t{n_emp}\t{n_post}\t{n_masked}\t{undercount_pct:.2f}\t{emp_rate:.2f}\t{post_rate:.2f}\n")
            print(f"[{cohort} Overall] Empirical Hits: {n_emp}/{n_tot} ({emp_rate:.2f}%), Post-hoc Hits: {n_post}/{n_tot} ({post_rate:.2f}%)")
            print(f"  -> Masked Alternative Hits (Top-1 Bias): {n_masked} ({undercount_pct:.2f}% of true hits missed by post-hoc filter!)")

            # Per length bin
            for lbin in ["<100", "100-149", "150-249", "250-399", ">=400"]:
                b_rows = [r for r in c_rows if r["length_bin"] == lbin]
                bn_tot = len(b_rows)
                bn_emp = sum(1 for r in b_rows if r["has_empirical"])
                bn_post = sum(1 for r in b_rows if r["has_posthoc"])
                bn_masked = sum(1 for r in b_rows if r["category"] == "masked_alternative_hit")
                b_undercount = (100.0 * bn_masked / bn_emp) if bn_emp > 0 else 0.0
                b_emp_rate = 100.0 * bn_emp / bn_tot
                b_post_rate = 100.0 * bn_post / bn_tot
                out.write(f"{cohort}\t{lbin}\t{bn_tot}\t{bn_emp}\t{bn_post}\t{bn_masked}\t{b_undercount:.2f}\t{b_emp_rate:.2f}\t{b_post_rate:.2f}\n")

    print(f"Wrote summary to {summary_tsv}")

    # Clean up large temporary FASTA and diamond database
    print("Step 8: Cleaning up temporary intermediate files...")
    if os.path.exists(sub_faa):
        os.remove(sub_faa)
    if os.path.exists(sub_dmnd + ".dmnd"):
        os.remove(sub_dmnd + ".dmnd")
    print("Cleaned up temporary FASTA and .dmnd files.")

    print(f"DONE_03_PILOT_SEARCH in {time.time() - t0:.2f}s total.")

if __name__ == "__main__":
    main()
