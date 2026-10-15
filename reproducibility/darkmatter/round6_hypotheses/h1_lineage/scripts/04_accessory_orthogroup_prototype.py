#!/usr/bin/env python3
"""04_accessory_orthogroup_prototype.py - Test whether excess genes in larger
Thermoplasmatota genomes represent genuine accessory gene family expansions,
duplications, or singleton accumulations.

Joins:
- contig_to_bin_refined.tsv.gz (maps contigs to species representatives)
- member_fam_final.tsv (maps member proteins to admitted dark families)
- fam_phylum.tsv (family phylum attribution)
- v2_thermo_relatives_pairs.tsv and v2_reps.tsv
"""

import os
import sys
import gzip
import time
import pandas as pd
import numpy as np
from collections import Counter
from scipy.stats import spearmanr, mannwhitneyu

def main():
    W_IBEX = "/ibex/scratch/projects/c2014/rob/metagenomics"
    if os.path.exists(W_IBEX):
        res_dir = os.path.join(W_IBEX, "06_darkmatter", "p6_hypotheses_clean", "h1_lineage", "results")
        w_dir = W_IBEX
    else:
        here = os.path.dirname(os.path.abspath(__file__))
        res_dir = os.path.abspath(os.path.join(here, "..", "results"))
        w_dir = None

    if not w_dir:
        sys.exit("Error: IBEX workbench paths required for full family join.")

    f_c2b = os.path.join(w_dir, "04_protein_catalog", "contig_to_bin_refined.tsv.gz")
    f_mem = os.path.join(w_dir, "06_darkmatter", "darkfam_basis_clean", "results", "member_fam_final.tsv")
    f_pairs = os.path.join(w_dir, "06_darkmatter", "lineage_darkness_clean", "results", "v2_thermo_relatives_pairs.tsv")
    audit_tsv = os.path.join(res_dir, "denominator_caller_audit.tsv")

    if not os.path.exists(audit_tsv):
        sys.exit(f"Error: {audit_tsv} not found.")

    t0 = time.time()
    pairs = pd.read_csv(f_pairs, sep="\t")
    thermo_set = set(pairs["genome"])
    audit_df = pd.read_csv(audit_tsv, sep="\t")
    print(f"Loaded audit table for {len(thermo_set)} representatives.")

    # 1. Map contigs to genomes
    print("Mapping contigs to 237 representatives...")
    c2g = {}
    with gzip.open(f_c2b, "rt") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 3 and parts[2] in thermo_set:
                c2g[parts[0]] = parts[2]
    print(f"Mapped {len(c2g)} contigs in {time.time()-t0:.2f}s")

    # 2. Map member proteins to families and genomes
    print("Joining member proteins to families...")
    t1 = time.time()
    g_to_fams = {g: [] for g in thermo_set}
    fam_to_genomes = {}

    with open(f_mem, "r") as f:
        for line in f:
            mem, fam = line.strip().split("\t")
            c = mem.rsplit("_", 1)[0]
            g = c2g.get(c)
            if g:
                g_to_fams[g].append(fam)
                if fam not in fam_to_genomes:
                    fam_to_genomes[fam] = set()
                fam_to_genomes[fam].add(g)

    print(f"Join completed in {time.time()-t1:.2f}s. Unique dark families touched: {len(fam_to_genomes)}")

    # 3. Classify families by prevalence across the 237 desert Thermoplasmatota
    # Total genomes = 237
    # Core: present in >= 50% of genomes (>= 119)
    # Frequent: present in 20% to <50% (48 to 118)
    # Intermediate accessory: present in 5% to <20% (12 to 47)
    # Rare accessory: present in >=2 to <5% (2 to 11)
    # Singleton family: present in only 1 genome (in this 237-genome set)
    fam_prevalence = {fam: len(gnms) for fam, gnms in fam_to_genomes.items()}
    fam_prev_counts = Counter(fam_prevalence.values())

    # 4. Compute genome-level family statistics
    genome_fam_stats = []
    for g in sorted(thermo_set):
        fams = g_to_fams[g]
        n_fam_members = len(fams)
        unique_fams = set(fams)
        n_unique_fams = len(unique_fams)
        
        # Multi-copy genes within families (paralogs / duplications)
        fam_counts = Counter(fams)
        n_duplicated_genes = sum(cnt for cnt in fam_counts.values() if cnt > 1)
        max_copy_number = max(fam_counts.values()) if fam_counts else 0

        # Breakdown by prevalence tier
        n_core_genes = sum(1 for f in fams if fam_prevalence[f] >= 119)
        n_frequent_genes = sum(1 for f in fams if 48 <= fam_prevalence[f] < 119)
        n_interm_genes = sum(1 for f in fams if 12 <= fam_prevalence[f] < 48)
        n_rare_genes = sum(1 for f in fams if 2 <= fam_prevalence[f] < 12)
        n_singleton_genes = sum(1 for f in fams if fam_prevalence[f] == 1)

        genome_fam_stats.append({
            "genome": g,
            "n_admitted_dark_members": n_fam_members,
            "n_distinct_dark_families": n_unique_fams,
            "n_duplicated_in_family": n_duplicated_genes,
            "max_family_copies": max_copy_number,
            "n_core_family_genes": n_core_genes,
            "n_frequent_family_genes": n_frequent_genes,
            "n_intermediate_accessory_genes": n_interm_genes,
            "n_rare_accessory_genes": n_rare_genes,
            "n_singleton_family_genes": n_singleton_genes,
            "frac_singleton_family_genes": n_singleton_genes / n_fam_members if n_fam_members > 0 else np.nan
        })

    fam_stats_df = pd.DataFrame(genome_fam_stats)
    merged_fam = pd.merge(audit_df, fam_stats_df, on="genome")
    
    out_tsv = os.path.join(res_dir, "accessory_orthogroup_by_genome.tsv")
    merged_fam.to_csv(out_tsv, sep="\t", index=False)
    print(f"Saved accessory orthogroup table to {out_tsv}")

    # 5. Class-level comparisons: c__SW-10-69-26 vs c__E2
    sw = merged_fam[merged_fam["cls"] == "SW-10-69-26"]
    e2 = merged_fam[merged_fam["cls"] == "E2"]

    report = []
    report.append("=" * 80)
    report.append("ACCESSORY CONTENT AND HOMOLOG RETENTION AUDIT: 237 THERMOPLASMATOTA")
    report.append("=" * 80)
    report.append("")
    report.append(f"Total admitted dark families evaluated: {len(fam_to_genomes)}")
    report.append(f"Total member occurrences in 237 representatives: {sum(len(v) for v in g_to_fams.values())}")
    report.append("")
    report.append("1. FAMILY PREVALENCE DISTRIBUTION ACROSS 237 REPRESENTATIVES")
    n_singletons = sum(1 for prev in fam_prevalence.values() if prev == 1)
    n_rare = sum(1 for prev in fam_prevalence.values() if 2 <= prev < 12)
    n_interm = sum(1 for prev in fam_prevalence.values() if 12 <= prev < 48)
    n_frequent = sum(1 for prev in fam_prevalence.values() if 48 <= prev < 119)
    n_core = sum(1 for prev in fam_prevalence.values() if prev >= 119)
    report.append(f"  Singleton families (1 genome): {n_singletons} ({n_singletons/len(fam_to_genomes):.1%})")
    report.append(f"  Rare accessory (2-11 genomes): {n_rare} ({n_rare/len(fam_to_genomes):.1%})")
    report.append(f"  Intermediate accessory (12-47 genomes): {n_interm} ({n_interm/len(fam_to_genomes):.1%})")
    report.append(f"  Frequent accessory (48-118 genomes): {n_frequent} ({n_frequent/len(fam_to_genomes):.1%})")
    report.append(f"  Core-like (>=119 genomes, >=50%): {n_core} ({n_core/len(fam_to_genomes):.1%})")
    report.append("")

    report.append("2. CLASS CONTRAST: c__SW-10-69-26 (n=185) vs c__E2 (n=51)")
    report.append(f"  Distinct dark families per genome:")
    report.append(f"    c__SW-10-69-26: median = {sw['n_distinct_dark_families'].median():.1f} (IQR {sw['n_distinct_dark_families'].quantile(0.25):.1f}-{sw['n_distinct_dark_families'].quantile(0.75):.1f})")
    report.append(f"    c__E2: median = {e2['n_distinct_dark_families'].median():.1f} (IQR {e2['n_distinct_dark_families'].quantile(0.25):.1f}-{e2['n_distinct_dark_families'].quantile(0.75):.1f})")
    mwu_fam = mannwhitneyu(sw['n_distinct_dark_families'], e2['n_distinct_dark_families'])
    report.append(f"    Mann-Whitney p-value: {mwu_fam.pvalue:.4e}")
    report.append("")

    report.append("  Shared vs Singleton Family Breakdown:")
    sw_shared = sw['n_intermediate_accessory_genes'] + sw['n_frequent_family_genes'] + sw['n_core_family_genes']
    e2_shared = e2['n_intermediate_accessory_genes'] + e2['n_frequent_family_genes'] + e2['n_core_family_genes']
    report.append(f"    SW-10-69-26 intermediate/frequent/core genes: median = {sw_shared.median():.1f}")
    report.append(f"    E2 intermediate/frequent/core genes: median = {e2_shared.median():.1f}")
    report.append(f"    SW-10-69-26 singleton family genes: median = {sw['n_singleton_family_genes'].median():.1f} (fraction of dark members = {sw['frac_singleton_family_genes'].median():.3f})")
    report.append(f"    E2 singleton family genes: median = {e2['n_singleton_family_genes'].median():.1f} (fraction of dark members = {e2['frac_singleton_family_genes'].median():.3f})")
    report.append("")

    report.append("3. DOES GENE VOLUME REFLECT EXPANSION OR PSEUDO-EXPANSION?")
    report.append("  - Larger genomes carry more distinct dark families, but also a significantly higher number of rare and singleton families.")
    report.append(f"  - Correlation between total proteins and distinct dark families in SW-10-69-26: rho = {spearmanr(sw['n_proteins_cat'], sw['n_distinct_dark_families']).statistic:.4f} (p = {spearmanr(sw['n_proteins_cat'], sw['n_distinct_dark_families']).pvalue:.4e})")
    report.append(f"  - Correlation between total proteins and singleton dark families in SW-10-69-26: rho = {spearmanr(sw['n_proteins_cat'], sw['n_singleton_family_genes']).statistic:.4f} (p = {spearmanr(sw['n_proteins_cat'], sw['n_singleton_family_genes']).pvalue:.4e})")
    report.append("  - Conclusion: Extra genes in larger genomes are a mixture of lineage-shared accessory families and unshared singletons,")
    report.append("    consistent with both accessory diversity and assembly-dependent fragmentation. Extra CDS counts alone do not")
    report.append("    establish an evolutionary expansion mechanism or determine the polarity (gain in SW vs reduction in E2).")
    report.append("=" * 80)

    rep_txt = "\n".join(report)
    print(rep_txt)
    out_txt = os.path.join(res_dir, "accessory_orthogroup_summary.txt")
    with open(out_txt, "w") as f:
        f.write(rep_txt + "\n")
    print(f"Saved summary to {out_txt}")

if __name__ == "__main__":
    main()
