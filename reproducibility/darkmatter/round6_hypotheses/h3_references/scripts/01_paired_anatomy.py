#!/usr/bin/env python3
"""
scripts/01_paired_anatomy.py
Perform paired quadrant decomposition, McNemar symmetry test, length stratification,
and alignment quality comparison between SMAG and GEM hits across the 106,042 dark
and 106,042 control families.

Usage:
  python3 scripts/01_paired_anatomy.py
"""

import sys
import os
import math
import numpy as np
from scipy import stats

def main():
    base_dir = "/ibex/scratch/projects/c2014/rob/metagenomics"
    spec_dir = os.path.join(base_dir, "06_darkmatter/specificity_clean/results")
    basis_dir = os.path.join(base_dir, "06_darkmatter/darkfam_basis_clean/results")
    out_dir = "/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p6_hypotheses_clean/h3_references/results"
    os.makedirs(out_dir, exist_ok=True)

    print("Step 1: Reading query sequences and lengths...")
    query_faa = os.path.join(spec_dir, "query_all.faa")
    query_lens = {}
    current_id = None
    current_len = 0
    with open(query_faa, "r") as f:
        for line in f:
            if line.startswith(">"):
                if current_id:
                    query_lens[current_id] = current_len
                current_id = line[1:].strip().split()[0]
                current_len = 0
            else:
                current_len += len(line.strip())
        if current_id:
            query_lens[current_id] = current_len
    print(f"Loaded {len(query_lens)} query lengths.")

    print("Step 2: Loading family sets...")
    with open(os.path.join(basis_dir, "FAMILIES_FINAL.txt")) as f:
        dark_fams = [line.strip() for line in f if line.strip()]
    with open(os.path.join(basis_dir, "CONTROL_FAMILIES.txt")) as f:
        ctrl_fams = [line.strip() for line in f if line.strip()]
    print(f"Dark families: {len(dark_fams)}, Control families: {len(ctrl_fams)}")

    print("Step 3: Loading SMAG and GEM hits...")
    smag_hits = {}
    with open(os.path.join(spec_dir, "smag_hits.m8")) as f:
        for line in f:
            parts = line.strip().split("\t")
            # qseqid sseqid pident evalue
            smag_hits[parts[0]] = (parts[1], float(parts[2]), float(parts[3]))
            
    gem_hits = {}
    with open(os.path.join(spec_dir, "gem_hits.m8")) as f:
        for line in f:
            parts = line.strip().split("\t")
            gem_hits[parts[0]] = (parts[1], float(parts[2]), float(parts[3]))
    print(f"SMAG hits: {len(smag_hits)}, GEM hits: {len(gem_hits)}")

    # 4-quadrant decomposition
    quad_summary_path = os.path.join(out_dir, "paired_quadrant_summary.tsv")
    with open(quad_summary_path, "w") as out:
        out.write("cohort\ttotal\tboth\tsmag_only\tgem_only\tneither\tpct_both\tpct_smag_only\tpct_gem_only\tpct_neither\tmcnemar_chi2\tmcnemar_pvalue\todds_ratio_smag_vs_gem\n")
        
        for cohort_name, fam_list in [("dark", dark_fams), ("ctrl", ctrl_fams)]:
            tot = len(fam_list)
            both = 0
            smag_only = 0
            gem_only = 0
            neither = 0
            for q in fam_list:
                in_s = q in smag_hits
                in_g = q in gem_hits
                if in_s and in_g:
                    both += 1
                elif in_s and not in_g:
                    smag_only += 1
                elif not in_s and in_g:
                    gem_only += 1
                else:
                    neither += 1
            
            # McNemar test with continuity correction: (|b - c| - 1)^2 / (b + c)
            b = smag_only
            c = gem_only
            chi2 = (abs(b - c) - 1.0)**2 / (b + c) if (b + c) > 0 else 0.0
            pval = stats.chi2.sf(chi2, df=1) if (b + c) > 0 else 1.0
            odds_ratio = (b / c) if c > 0 else float("inf")
            
            p_both = 100.0 * both / tot
            p_smag = 100.0 * smag_only / tot
            p_gem = 100.0 * gem_only / tot
            p_neither = 100.0 * neither / tot
            
            out.write(f"{cohort_name}\t{tot}\t{both}\t{smag_only}\t{gem_only}\t{neither}\t{p_both:.2f}\t{p_smag:.2f}\t{p_gem:.2f}\t{p_neither:.2f}\t{chi2:.2f}\t{pval:.4e}\t{odds_ratio:.2f}\n")
            print(f"[{cohort_name}] Both: {both} ({p_both:.2f}%), SMAG-only: {smag_only} ({p_smag:.2f}%), GEM-only: {gem_only} ({p_gem:.2f}%), Neither: {neither} ({p_neither:.2f}%)")
            print(f"[{cohort_name}] McNemar chi2: {chi2:.2f}, p-val: {pval:.4e}, Odds Ratio (SMAG/GEM discordant): {odds_ratio:.2f}")

    print(f"Wrote quadrant summary to {quad_summary_path}")

    # Length stratification
    length_bins = [
        ("<100", lambda l: l < 100),
        ("100-149", lambda l: 100 <= l <= 149),
        ("150-249", lambda l: 150 <= l <= 249),
        ("250-399", lambda l: 250 <= l <= 399),
        (">=400", lambda l: l >= 400)
    ]
    
    len_summary_path = os.path.join(out_dir, "length_stratified_quadrants.tsv")
    with open(len_summary_path, "w") as out:
        out.write("cohort\tlength_bin\tn_fams\tboth\tsmag_only\tgem_only\tneither\tpct_both\tpct_smag_only\tpct_gem_only\tpct_neither\tsmag_hit_rate\tgem_hit_rate\tdiscordant_ratio_smag_to_gem\n")
        
        for cohort_name, fam_list in [("dark", dark_fams), ("ctrl", ctrl_fams)]:
            for bin_name, bin_func in length_bins:
                subset = [q for q in fam_list if bin_func(query_lens.get(q, 0))]
                tot = len(subset)
                if tot == 0:
                    continue
                both = sum(1 for q in subset if q in smag_hits and q in gem_hits)
                smag_only = sum(1 for q in subset if q in smag_hits and q not in gem_hits)
                gem_only = sum(1 for q in subset if q not in smag_hits and q in gem_hits)
                neither = sum(1 for q in subset if q not in smag_hits and q not in gem_hits)
                
                p_both = 100.0 * both / tot
                p_smag = 100.0 * smag_only / tot
                p_gem = 100.0 * gem_only / tot
                p_neither = 100.0 * neither / tot
                smag_hit_rate = 100.0 * (both + smag_only) / tot
                gem_hit_rate = 100.0 * (both + gem_only) / tot
                disc_ratio = (smag_only / gem_only) if gem_only > 0 else float("inf")
                
                out.write(f"{cohort_name}\t{bin_name}\t{tot}\t{both}\t{smag_only}\t{gem_only}\t{neither}\t{p_both:.2f}\t{p_smag:.2f}\t{p_gem:.2f}\t{p_neither:.2f}\t{smag_hit_rate:.2f}\t{gem_hit_rate:.2f}\t{disc_ratio:.2f}\n")

    print(f"Wrote length-stratified summary to {len_summary_path}")

    # Percent identity and E-value comparison for shared hits
    ident_summary_path = os.path.join(out_dir, "shared_hits_identity_comparison.tsv")
    with open(ident_summary_path, "w") as out:
        out.write("cohort\tn_shared\tsmag_pident_mean\tsmag_pident_median\tgem_pident_mean\tgem_pident_median\tpident_diff_mean\tpident_diff_median\twilcoxon_pident_stat\twilcoxon_pident_pval\tsmag_better_evalue_pct\tgem_better_evalue_pct\ttied_evalue_pct\n")
        
        for cohort_name, fam_list in [("dark", dark_fams), ("ctrl", ctrl_fams)]:
            shared = [q for q in fam_list if q in smag_hits and q in gem_hits]
            if not shared:
                continue
            
            s_pident = np.array([smag_hits[q][1] for q in shared])
            g_pident = np.array([gem_hits[q][1] for q in shared])
            s_eval = np.array([smag_hits[q][2] for q in shared])
            g_eval = np.array([gem_hits[q][2] for q in shared])
            
            diff = s_pident - g_pident
            
            # Wilcoxon signed-rank test
            w_res = stats.wilcoxon(s_pident, g_pident)
            
            smag_better_eval = np.sum(s_eval < g_eval)
            gem_better_eval = np.sum(g_eval < s_eval)
            tied_eval = np.sum(s_eval == g_eval)
            n_tot = len(shared)
            
            out.write(f"{cohort_name}\t{n_tot}\t{np.mean(s_pident):.2f}\t{np.median(s_pident):.2f}\t{np.mean(g_pident):.2f}\t{np.median(g_pident):.2f}\t{np.mean(diff):.2f}\t{np.median(diff):.2f}\t{w_res.statistic:.2f}\t{w_res.pvalue:.4e}\t{100.0*smag_better_eval/n_tot:.2f}\t{100.0*gem_better_eval/n_tot:.2f}\t{100.0*tied_eval/n_tot:.2f}\n")
            print(f"[{cohort_name} Shared N={n_tot}] SMAG pident median: {np.median(s_pident):.2f}%, GEM pident median: {np.median(g_pident):.2f}%, Diff median: {np.median(diff):.2f}%")
            print(f"[{cohort_name}] SMAG better e-val: {100.0*smag_better_eval/n_tot:.2f}%, GEM better e-val: {100.0*gem_better_eval/n_tot:.2f}%")

    print(f"Wrote identity comparison to {ident_summary_path}")
    print("DONE_01_PAIRED_ANATOMY")

if __name__ == "__main__":
    main()
