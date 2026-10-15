#!/usr/bin/env python3
"""01_denominator_caller_audit.py - Reconcile gene-caller basis, protein denominators,
and partial/complete gene structure across the 237 desert Thermoplasmatota representatives.

Audits:
1. Catalogue protein count (v2_reps.tsv) vs single-genome Pyrodigal count (v2_thermo_relatives_aai.tsv.gz).
2. RBH fraction denominator: min(N_a, N_b) vs N_desert vs N_nearest.
3. Protein completeness (partial=00 vs partial!=00) and length (>=100aa vs <100aa) from CheckM2 .faa calls.
4. Lineage-stratified comparisons between c__SW-10-69-26, c__E2, and orders.
"""

import os
import sys
import glob
import gzip
import pandas as pd
import numpy as np

def main():
    # Detect environment (IBEX workbench vs local)
    W_IBEX = "/ibex/scratch/projects/c2014/rob/metagenomics"
    if os.path.exists(W_IBEX):
        w_dir = W_IBEX
        faa_dir = os.path.join(w_dir, "05_refine", "_checkm2_refined", "protein_files")
        aai_path = os.path.join(w_dir, "06_darkmatter", "lineage_darkness_clean", "results", "v2_thermo_relatives_aai.tsv.gz")
        pairs_path = os.path.join(w_dir, "06_darkmatter", "lineage_darkness_clean", "results", "v2_thermo_relatives_pairs.tsv")
        reps_path = os.path.join(w_dir, "06_darkmatter", "lineage_darkness_clean", "results", "v2_reps.tsv")
        out_dir = os.path.join(w_dir, "06_darkmatter", "p6_hypotheses_clean", "h1_lineage", "results")
    else:
        # Local paths relative to repo root
        here = os.path.dirname(os.path.abspath(__file__))
        repo_root = os.path.abspath(os.path.join(here, "..", "..", "..", ".."))
        pairs_path = os.path.join(repo_root, "reviews", "darkmatter-2026-10-13", "round5", "reviewer_inputs",
                                  "06_darkmatter", "lineage_darkness_clean", "results", "v2_thermo_relatives_pairs.tsv")
        reps_path = os.path.join(repo_root, "reproducibility", "darkmatter", "inputs",
                                 "06_darkmatter", "lineage_darkness_clean", "results", "v2_reps.tsv.gz")
        aai_path = None
        faa_dir = None
        out_dir = os.path.abspath(os.path.join(here, "..", "results"))

    os.makedirs(out_dir, exist_ok=True)
    print(f"Running denominator and caller audit. Output dir: {out_dir}")

    # 1. Load pairs and reps
    pairs = pd.read_csv(pairs_path, sep="\t")
    print(f"Loaded pairs: {len(pairs)} rows")

    reps = pd.read_csv(reps_path, sep="\t")
    reps_thermo = reps[reps["phylum"].str.contains("Thermoplasmatota", na=False)].copy()
    print(f"Loaded reps: {len(reps_thermo)} Thermoplasmatota representatives")

    # 2. Load AAI table for raw single-genome counts
    prot_map = {}
    if aai_path and os.path.exists(aai_path):
        print("Reading AAI table for single-genome protein counts...")
        aai = pd.read_csv(aai_path, sep="\t")
        for _, row in aai.iterrows():
            ga, gb = row["genome_a"], row["genome_b"]
            if ga not in prot_map:
                prot_map[ga] = row["n_prot_a"]
            if gb not in prot_map:
                prot_map[gb] = row["n_prot_b"]
        print(f"Mapped {len(prot_map)} genomes from AAI table")

    # 3. Process CheckM2 .faa files for completeness & length
    faa_stats = {}
    if faa_dir and os.path.exists(faa_dir):
        print(f"Scanning CheckM2 .faa files in {faa_dir} for 237 representatives...")
        for g in pairs["genome"]:
            faa_file = os.path.join(faa_dir, f"{g}.faa")
            if not os.path.exists(faa_file):
                print(f"Warning: missing faa for {g}")
                continue
            n_tot = 0
            n_comp = 0
            n_ge100 = 0
            n_comp_ge100 = 0
            lengths = []
            with open(faa_file, "r") as f:
                cur_hdr = None
                cur_seq = []
                for line in f:
                    if line.startswith(">"):
                        if cur_hdr is not None:
                            slen = len("".join(cur_seq))
                            lengths.append(slen)
                            is_comp = "partial=00" in cur_hdr
                            if is_comp:
                                n_comp += 1
                            if slen >= 100:
                                n_ge100 += 1
                            if is_comp and slen >= 100:
                                n_comp_ge100 += 1
                        cur_hdr = line.strip()
                        cur_seq = []
                        n_tot += 1
                    else:
                        cur_seq.append(line.strip())
                if cur_hdr is not None:
                    slen = len("".join(cur_seq))
                    lengths.append(slen)
                    is_comp = "partial=00" in cur_hdr
                    if is_comp:
                        n_comp += 1
                    if slen >= 100:
                        n_ge100 += 1
                    if is_comp and slen >= 100:
                        n_comp_ge100 += 1
            faa_stats[g] = {
                "n_proteins_checkm2": n_tot,
                "n_complete": n_comp,
                "n_partial": n_tot - n_comp,
                "frac_complete": n_comp / n_tot if n_tot > 0 else np.nan,
                "n_ge100aa": n_ge100,
                "n_complete_ge100aa": n_comp_ge100,
                "frac_complete_ge100aa": n_comp_ge100 / n_tot if n_tot > 0 else np.nan,
                "median_protein_len": float(np.median(lengths)) if lengths else np.nan,
                "mean_protein_len": float(np.mean(lengths)) if lengths else np.nan
            }
        print(f"Parsed FASTA headers for {len(faa_stats)} representatives")

    # 4. Assemble merged audit table
    merged = pd.merge(pairs, reps_thermo[["genome", "completeness", "contamination", "length", "N50",
                                          "n_proteins", "n_dark", "dark_fraction", "resid_A"]], on="genome")
    
    merged["n_proteins_cat"] = merged["n_proteins"]
    merged["n_prot_aai"] = merged["genome"].map(prot_map)
    merged["n_prot_nearest"] = merged["aai_nearest"].map(prot_map)

    # Denominators
    merged["desert_is_smaller"] = merged["n_prot_aai"] < merged["n_prot_nearest"]
    merged["frac_rbh_of_desert_aai"] = merged["aai_n_rbh"] / merged["n_prot_aai"]
    merged["frac_rbh_of_desert_cat"] = merged["aai_n_rbh"] / merged["n_proteins_cat"]
    merged["frac_rbh_of_nearest"] = merged["aai_n_rbh"] / merged["n_prot_nearest"]

    # Incorporate CheckM2 stats if available
    if faa_stats:
        fs_df = pd.DataFrame.from_dict(faa_stats, orient="index").reset_index().rename(columns={"index": "genome"})
        merged = pd.merge(merged, fs_df, on="genome", how="left")

    # Save detailed per-genome audit table
    out_tsv = os.path.join(out_dir, "denominator_caller_audit.tsv")
    merged.to_csv(out_tsv, sep="\t", index=False)
    print(f"Saved audit table to {out_tsv}")

    # Generate summary report text
    report_lines = []
    report_lines.append("=" * 80)
    report_lines.append("DENOMINATOR AND CALLER AUDIT SUMMARY: 237 DESERT THERMOPLASMATOTA")
    report_lines.append("=" * 80)
    report_lines.append("")
    report_lines.append("1. TAXONOMIC BREAKDOWN")
    report_lines.append(f"Total representatives: {len(merged)}")
    for cls_name, cnt in merged['cls'].value_counts().items():
        report_lines.append(f"  Class {cls_name}: n = {cnt} ({cnt/len(merged):.1%})")
    report_lines.append("")

    report_lines.append("2. GENE-CALLER COMPARISON (Catalogue vs Single-Genome Pyrodigal)")
    if merged['n_prot_aai'].notna().any():
        diff = merged['n_prot_aai'] - merged['n_proteins_cat']
        ratio = merged['n_prot_aai'] / merged['n_proteins_cat']
        report_lines.append(f"  Catalogue n_proteins: mean = {merged['n_proteins_cat'].mean():.1f}, median = {merged['n_proteins_cat'].median():.1f}, range = [{merged['n_proteins_cat'].min()}, {merged['n_proteins_cat'].max()}]")
        report_lines.append(f"  Single-genome n_prot_aai: mean = {merged['n_prot_aai'].mean():.1f}, median = {merged['n_prot_aai'].median():.1f}, range = [{merged['n_prot_aai'].min():.0f}, {merged['n_prot_aai'].max():.0f}]")
        report_lines.append(f"  Difference (AAI - Cat): mean = {diff.mean():+.1f}, median = {diff.median():+.1f}, range = [{diff.min():+.0f}, {diff.max():+.0f}]")
        report_lines.append(f"  Ratio (AAI / Cat): mean = {ratio.mean():.4f}, median = {ratio.median():.4f}")
        report_lines.append(f"  Pearson correlation between callers: r = {merged['n_proteins_cat'].corr(merged['n_prot_aai']):.6f}")
    report_lines.append("")

    report_lines.append("3. RBH FRACTION DENOMINATOR DISTORTION")
    n_desert_larger = (merged['desert_is_smaller'] == False).sum()
    pct_desert_larger = n_desert_larger / len(merged) * 100
    report_lines.append(f"  Genomes where desert MAG is LARGER than nearest comparator: {n_desert_larger} / {len(merged)} ({pct_desert_larger:.1f}%)")
    report_lines.append(f"  Published aai_frac_rbh (divided by min(N_a, N_b)): median = {merged['aai_frac_rbh'].median():.4f}, IQR = [{merged['aai_frac_rbh'].quantile(0.25):.4f}, {merged['aai_frac_rbh'].quantile(0.75):.4f}]")
    if 'frac_rbh_of_desert_aai' in merged:
        report_lines.append(f"  Corrected frac_rbh_of_desert (divided by desert MAG N): median = {merged['frac_rbh_of_desert_aai'].median():.4f}, IQR = [{merged['frac_rbh_of_desert_aai'].quantile(0.25):.4f}, {merged['frac_rbh_of_desert_aai'].quantile(0.75):.4f}]")
        report_lines.append(f"  frac_rbh_of_nearest (divided by comparator N): median = {merged['frac_rbh_of_nearest'].median():.4f}, IQR = [{merged['frac_rbh_of_nearest'].quantile(0.25):.4f}, {merged['frac_rbh_of_nearest'].quantile(0.75):.4f}]")
    report_lines.append("")

    if 'frac_complete' in merged:
        report_lines.append("4. CALLABLE GENE STRUCTURE AND FRAGMENTATION SENSITIVITY")
        report_lines.append(f"  Complete gene fraction (partial=00): mean = {merged['frac_complete'].mean():.3f}, median = {merged['frac_complete'].median():.3f}, range = [{merged['frac_complete'].min():.3f}, {merged['frac_complete'].max():.3f}]")
        report_lines.append(f"  Complete & >=100aa fraction: mean = {merged['frac_complete_ge100aa'].mean():.3f}, median = {merged['frac_complete_ge100aa'].median():.3f}, range = [{merged['frac_complete_ge100aa'].min():.3f}, {merged['frac_complete_ge100aa'].max():.3f}]")
        report_lines.append(f"  Median protein length: mean across genomes = {merged['median_protein_len'].mean():.1f} aa, median = {merged['median_protein_len'].median():.1f} aa")
        
        # Compare classes
        sw = merged[merged['cls'] == 'SW-10-69-26']
        e2 = merged[merged['cls'] == 'E2']
        report_lines.append("  By Class:")
        report_lines.append(f"    c__SW-10-69-26 (n={len(sw)}): frac_complete median = {sw['frac_complete'].median():.3f}, frac_complete_ge100aa median = {sw['frac_complete_ge100aa'].median():.3f}")
        report_lines.append(f"    c__E2 (n={len(e2)}): frac_complete median = {e2['frac_complete'].median():.3f}, frac_complete_ge100aa median = {e2['frac_complete_ge100aa'].median():.3f}")
    report_lines.append("")
    report_lines.append("=" * 80)

    report_text = "\n".join(report_lines)
    print(report_text)
    out_txt = os.path.join(out_dir, "denominator_caller_audit_summary.txt")
    with open(out_txt, "w") as f:
        f.write(report_text + "\n")
    print(f"Saved summary report to {out_txt}")

if __name__ == "__main__":
    main()
