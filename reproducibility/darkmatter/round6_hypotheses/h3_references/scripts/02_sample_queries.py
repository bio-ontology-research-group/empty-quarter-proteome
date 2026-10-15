#!/usr/bin/env python3
"""
scripts/02_sample_queries.py
Generate a deterministic length-stratified query subset of 2,000 family representatives
(1,000 dark and 1,000 matched control; 200 per length stratum) for the empirical search pilot.

Usage:
  python3 scripts/02_sample_queries.py
"""

import os
import sys
import random

def main():
    base_dir = "/ibex/scratch/projects/c2014/rob/metagenomics"
    spec_dir = os.path.join(base_dir, "06_darkmatter/specificity_clean/results")
    basis_dir = os.path.join(base_dir, "06_darkmatter/darkfam_basis_clean/results")
    out_dir = "/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p6_hypotheses_clean/h3_references/results"
    os.makedirs(out_dir, exist_ok=True)

    seed = 20261013
    random.seed(seed)

    print("Step 1: Reading query sequences and lengths...")
    query_faa = os.path.join(spec_dir, "query_all.faa")
    query_seqs = {}
    current_id = None
    current_lines = []
    
    with open(query_faa, "r") as f:
        for line in f:
            if line.startswith(">"):
                if current_id:
                    seq_str = "".join(current_lines)
                    query_seqs[current_id] = (seq_str, len(seq_str))
                current_id = line[1:].strip().split()[0]
                current_lines = []
            else:
                current_lines.append(line.strip())
        if current_id:
            seq_str = "".join(current_lines)
            query_seqs[current_id] = (seq_str, len(seq_str))
    print(f"Loaded {len(query_seqs)} query sequences.")

    print("Step 2: Loading dark and control lists...")
    with open(os.path.join(basis_dir, "FAMILIES_FINAL.txt")) as f:
        dark_fams = [line.strip() for line in f if line.strip()]
    with open(os.path.join(basis_dir, "CONTROL_FAMILIES.txt")) as f:
        ctrl_fams = [line.strip() for line in f if line.strip()]

    length_bins = [
        ("<100", lambda l: l < 100),
        ("100-149", lambda l: 100 <= l <= 149),
        ("150-249", lambda l: 150 <= l <= 249),
        ("250-399", lambda l: 250 <= l <= 399),
        (">=400", lambda l: l >= 400)
    ]
    
    sample_dark = []
    sample_ctrl = []
    manifest_rows = []

    for bin_name, bin_func in length_bins:
        # dark
        dark_pool = [q for q in dark_fams if q in query_seqs and bin_func(query_seqs[q][1])]
        dark_pool.sort() # Ensure deterministic order before shuffle
        random.seed(seed + len(bin_name))
        chosen_dark = random.sample(dark_pool, 200)
        sample_dark.extend(chosen_dark)
        for q in chosen_dark:
            manifest_rows.append((q, "dark", bin_name, query_seqs[q][1]))

        # ctrl
        ctrl_pool = [q for q in ctrl_fams if q in query_seqs and bin_func(query_seqs[q][1])]
        ctrl_pool.sort()
        random.seed(seed + 100 + len(bin_name))
        chosen_ctrl = random.sample(ctrl_pool, 200)
        sample_ctrl.extend(chosen_ctrl)
        for q in chosen_ctrl:
            manifest_rows.append((q, "ctrl", bin_name, query_seqs[q][1]))

        print(f"Stratum {bin_name}: Sampled 200 dark from {len(dark_pool)}, 200 ctrl from {len(ctrl_pool)}")

    print(f"Total sampled: {len(sample_dark)} dark, {len(sample_ctrl)} ctrl.")

    # Write FASTA
    out_faa = os.path.join(out_dir, "pilot_query_2000.faa")
    with open(out_faa, "w") as out:
        for q in sample_dark + sample_ctrl:
            seq = query_seqs[q][0]
            out.write(f">{q}\n{seq}\n")
    print(f"Wrote {out_faa}")

    # Write Manifest
    out_man = os.path.join(out_dir, "pilot_query_manifest.tsv")
    with open(out_man, "w") as out:
        out.write("query_id\tcohort\tlength_bin\tlength_aa\n")
        for q, cohort, bin_name, length_aa in manifest_rows:
            out.write(f"{q}\t{cohort}\t{bin_name}\t{length_aa}\n")
    print(f"Wrote {out_man}")
    print("DONE_02_SAMPLE_QUERIES")

if __name__ == "__main__":
    main()
