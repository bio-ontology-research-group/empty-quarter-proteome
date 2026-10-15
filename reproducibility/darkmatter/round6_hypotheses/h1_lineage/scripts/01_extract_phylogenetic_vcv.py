#!/usr/bin/env python3
"""01_extract_phylogenetic_vcv.py - Parse ar53 FastTree Newick, validate tip mapping,
and construct the phylogenetic variance-covariance and correlation matrix for the
237 desert Thermoplasmatota representatives.

Executed using DendroPy in the verified drep environment:
/ibex/scratch/projects/c2014/rob/miniconda3/envs/drep/bin/python
"""

import os
import sys
import gzip
import time
import pandas as pd
import numpy as np

def main():
    W_IBEX = "/ibex/scratch/projects/c2014/rob/metagenomics"
    if os.path.exists(W_IBEX):
        w_dir = W_IBEX
        tree_path = os.path.join(w_dir, "06_darkmatter", "lineage_darkness_clean", "results", "thermo_tree.nwk")
        pairs_path = os.path.join(w_dir, "06_darkmatter", "lineage_darkness_clean", "results", "v2_thermo_relatives_pairs.tsv")
        out_dir = os.path.join(w_dir, "06_darkmatter", "p6_hypotheses_clean", "h1_lineage", "results")
    else:
        here = os.path.dirname(os.path.abspath(__file__))
        repo_root = os.path.abspath(os.path.join(here, "..", "..", "..", ".."))
        pairs_path = os.path.join(repo_root, "reviews", "darkmatter-2026-10-13", "round5", "reviewer_inputs",
                                  "06_darkmatter", "lineage_darkness_clean", "results", "v2_thermo_relatives_pairs.tsv")
        tree_path = None
        out_dir = os.path.abspath(os.path.join(here, "..", "results"))

    os.makedirs(out_dir, exist_ok=True)
    print(f"Extracting phylogenetic covariance matrix. Output: {out_dir}")

    # Load targets
    pairs = pd.read_csv(pairs_path, sep="\t")
    target_tips = sorted(list(pairs["genome"]))
    n_targets = len(target_tips)
    print(f"Target desert Thermoplasmatota representatives: {n_targets}")

    if not tree_path or not os.path.exists(tree_path):
        sys.exit(f"Error: tree file {tree_path} not accessible.")

    import dendropy
    t0 = time.time()
    print(f"Reading Newick tree from {tree_path}...")
    tree = dendropy.Tree.get(path=tree_path, schema="newick", preserve_underscores=True)
    all_tips = [leaf.taxon.label for leaf in tree.leaf_nodes() if leaf.taxon]
    print(f"Total tips in tree: {len(all_tips)} (read in {time.time()-t0:.2f}s)")

    # Validate tip mapping
    leaf_set = set(all_tips)
    matched = set(target_tips).intersection(leaf_set)
    print(f"Matched target tips: {len(matched)} / {n_targets}")
    if len(matched) < n_targets:
        missing = set(target_tips) - leaf_set
        sys.exit(f"Fatal error: missing tips in tree: {missing}")

    # Verify unique tip labels
    tip_counts = pd.Series(all_tips).value_counts()
    dups = tip_counts[tip_counts > 1]
    if len(dups) > 0:
        print(f"Warning: duplicated tip labels in tree: {len(dups)}")

    # Record tip mapping table
    tip_info = []
    for g in target_tips:
        tip_info.append({
            "genome": g,
            "in_tree": g in leaf_set,
            "tree_label": g
        })
    tip_df = pd.DataFrame(tip_info)
    tip_tsv = os.path.join(out_dir, "phylo_tree_tips.tsv")
    tip_df.to_csv(tip_tsv, sep="\t", index=False)
    print(f"Saved tip mapping to {tip_tsv}")

    # Root tree at midpoint for consistent root-to-tip shared distances
    t1 = time.time()
    tree.reroot_at_midpoint()
    print(f"Rerooted at midpoint in {time.time()-t1:.2f}s")

    # Prune tree to only the target 237 representatives
    t2 = time.time()
    tree.retain_taxa_with_labels(target_tips)
    tree.calc_node_root_distances()
    print(f"Pruned tree to {len(tree.leaf_nodes())} target leaves in {time.time()-t2:.2f}s")

    # Compute phylogenetic distance matrix
    t3 = time.time()
    pdm = tree.phylogenetic_distance_matrix()
    print(f"Computed cophenetic distance matrix in {time.time()-t3:.2f}s")

    # Construct variance-covariance matrix (C_ij = (d_Ri + d_Rj - d_ij) / 2)
    taxa = [tree.taxon_namespace.get_taxon(t) for t in target_tips]
    root_dist = np.array([tree.find_node_with_taxon_label(t).root_distance for t in target_tips])

    vcv = np.zeros((n_targets, n_targets))
    dist_mat = np.zeros((n_targets, n_targets))

    for i in range(n_targets):
        vcv[i, i] = root_dist[i]
        dist_mat[i, i] = 0.0
        for j in range(i + 1, n_targets):
            d_ij = pdm(taxa[i], taxa[j])
            dist_mat[i, j] = d_ij
            dist_mat[j, i] = d_ij
            c_ij = 0.5 * (root_dist[i] + root_dist[j] - d_ij)
            vcv[i, j] = c_ij
            vcv[j, i] = c_ij

    # Construct correlation matrix C_corr = VCV / sqrt(diag * diag')
    std_diag = np.sqrt(np.diag(vcv))
    corr_mat = vcv / np.outer(std_diag, std_diag)
    np.fill_diagonal(corr_mat, 1.0)
    # Ensure numerical symmetry and non-negative correlations
    corr_mat = 0.5 * (corr_mat + corr_mat.T)
    corr_mat = np.clip(corr_mat, 0.0, 1.0)

    # Save matrices as compressed TSV
    vcv_df = pd.DataFrame(vcv, index=target_tips, columns=target_tips)
    corr_df = pd.DataFrame(corr_mat, index=target_tips, columns=target_tips)
    dist_df = pd.DataFrame(dist_mat, index=target_tips, columns=target_tips)

    vcv_tsv = os.path.join(out_dir, "phylo_vcv.tsv.gz")
    corr_tsv = os.path.join(out_dir, "phylo_correlation.tsv.gz")
    dist_tsv = os.path.join(out_dir, "phylo_cophenetic_dist.tsv.gz")

    vcv_df.to_csv(vcv_tsv, sep="\t", compression="gzip", float_format="%.6f")
    corr_df.to_csv(corr_tsv, sep="\t", compression="gzip", float_format="%.6f")
    dist_df.to_csv(dist_tsv, sep="\t", compression="gzip", float_format="%.6f")

    print(f"Saved VCV matrix to {vcv_tsv}")
    print(f"Saved correlation matrix to {corr_tsv}")
    print(f"Saved cophenetic distance matrix to {dist_tsv}")

    # Summary report
    rep = []
    rep.append("=" * 80)
    rep.append("PHYLOGENETIC TREE AND COVARIANCE AUDIT: 237 THERMOPLASMATOTA")
    rep.append("=" * 80)
    rep.append(f"Tree source: {tree_path}")
    rep.append(f"Full tree tips: {len(all_tips)}")
    rep.append(f"Target representatives: {n_targets}")
    rep.append(f"Verified mapped tips: {len(matched)} / {n_targets} (100.0%)")
    rep.append("")
    rep.append("Root-to-Tip Distance Summary:")
    rep.append(f"  Min: {root_dist.min():.4f}, 25%: {np.percentile(root_dist, 25):.4f}, Median: {np.median(root_dist):.4f}, 75%: {np.percentile(root_dist, 75):.4f}, Max: {root_dist.max():.4f}")
    rep.append("")
    rep.append("Pairwise Cophenetic Distance Summary (Substitutions per site):")
    off_diag_dist = dist_mat[np.triu_indices_from(dist_mat, 1)]
    rep.append(f"  Min: {off_diag_dist.min():.4f}, 25%: {np.percentile(off_diag_dist, 25):.4f}, Median: {np.median(off_diag_dist):.4f}, 75%: {np.percentile(off_diag_dist, 75):.4f}, Max: {off_diag_dist.max():.4f}")
    rep.append("")
    rep.append("Phylogenetic Correlation Summary (Off-diagonal):")
    off_diag_corr = corr_mat[np.triu_indices_from(corr_mat, 1)]
    rep.append(f"  Min: {off_diag_corr.min():.4f}, 25%: {np.percentile(off_diag_corr, 25):.4f}, Median: {np.median(off_diag_corr):.4f}, 75%: {np.percentile(off_diag_corr, 75):.4f}, Max: {off_diag_corr.max():.4f}")
    rep.append("=" * 80)

    summary_text = "\n".join(rep)
    print(summary_text)
    with open(os.path.join(out_dir, "phylo_covariance_summary.txt"), "w") as f:
        f.write(summary_text + "\n")

if __name__ == "__main__":
    main()
