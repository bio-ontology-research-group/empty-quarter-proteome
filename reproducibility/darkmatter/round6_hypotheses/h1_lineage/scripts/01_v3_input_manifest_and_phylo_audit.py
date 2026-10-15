#!/usr/bin/env python3
"""01_v3_input_manifest_and_phylo_audit.py - Enumerate and checksum all raw inputs
read by the H1 v3 pipeline (eliminating directory/dynamic templates) and validate
phylogenetic covariance matrices under midpoint vs explicit outgroup rooting.

Executed using DendroPy in the verified drep environment:
/ibex/scratch/projects/c2014/rob/miniconda3/envs/drep/bin/python
"""

import os
import sys
import hashlib
import json
import time
import numpy as np
import pandas as pd


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    W_IBEX = "/ibex/scratch/projects/c2014/rob/metagenomics"
    if not os.path.exists(W_IBEX):
        sys.exit(f"Error: IBEX workbench {W_IBEX} required for full raw-input enumeration and tree audit.")

    here = os.path.dirname(os.path.abspath(__file__))
    h1_dir = os.path.abspath(os.path.join(here, ".."))
    out_dir = os.path.join(h1_dir, "results")
    os.makedirs(out_dir, exist_ok=True)

    t0 = time.time()
    print(f"[01_v3] Starting enumerated input manifest and phylogenetic VCV validation in {h1_dir}")

    # 1. Core input files on IBEX
    ld_res = os.path.join(W_IBEX, "06_darkmatter", "lineage_darkness_clean", "results")
    faa_dir = os.path.join(W_IBEX, "05_refine", "_checkm2_refined", "protein_files")

    core_inputs = [
        ("v2_thermo_relatives_pairs", os.path.join(ld_res, "v2_thermo_relatives_pairs.tsv"), "ibex",
         "237-row nearest relative table for desert Thermoplasmatota"),
        ("v2_reps", os.path.join(ld_res, "v2_reps.tsv"), "ibex",
         "5,438 species representatives with catalogue dark counts and CheckM2 quality metrics"),
        ("v2_thermo_relatives_aai", os.path.join(ld_res, "v2_thermo_relatives_aai.tsv.gz"), "ibex",
         "Pairwise RBH AAI and single-genome Pyrodigal protein counts"),
        ("v2_thermo_relatives_comparators", os.path.join(ld_res, "v2_thermo_relatives_comparators.tsv"), "ibex",
         "Manifest of 218 GTDB r232 species reps and 83 MGnify comparator genomes"),
        ("v2_thermo_relatives_gtdbtk_mgnify", os.path.join(ld_res, "v2_thermo_relatives_gtdbtk_mgnify.tsv"), "ibex",
         "GTDB-Tk r232 classifications of MGnify comparator genomes"),
        ("thermo_tree_nwk", os.path.join(ld_res, "thermo_tree.nwk"), "ibex",
         "FastTree LG+Gamma approximate-ML tree of 1,914 ar53 taxa (237 desert + 1,647 GTDB + 30 outgroup)"),
        ("thermo_tree_taxa", os.path.join(ld_res, "thermo_tree_taxa.tsv"), "ibex",
         "Taxon membership and GTDB lineage table for the 1,914 tips in thermo_tree.nwk"),
        ("contig_to_bin_refined", os.path.join(W_IBEX, "04_protein_catalog", "contig_to_bin_refined.tsv.gz"), "ibex",
         "Contig to refined bin mapping across the catalogue"),
        ("member_fam_final", os.path.join(W_IBEX, "06_darkmatter", "darkfam_basis_clean", "results", "member_fam_final.tsv"), "ibex",
         "Member protein to admitted dark family mapping"),
        ("fam_phylum", os.path.join(W_IBEX, "06_darkmatter", "darkfam_basis_clean", "results", "fam_phylum.tsv"), "ibex",
         "Phylum attribution for 78,727 dark families"),
    ]

    manifest_rows = []
    for input_id, fpath, host, desc in core_inputs:
        if not os.path.exists(fpath):
            sys.exit(f"Fatal error: required core input missing: {fpath}")
        sz = os.path.getsize(fpath)
        chksum = sha256_file(fpath)
        manifest_rows.append({
            "input_id": input_id,
            "file_path": fpath,
            "host": host,
            "access_status": "verified_readable",
            "size_bytes": sz,
            "sha256_checksum": chksum,
            "description": desc
        })

    # 2. Enumerate all 237 target genomes for CheckM2 .faa files
    pairs_path = os.path.join(ld_res, "v2_thermo_relatives_pairs.tsv")
    pairs = pd.read_csv(pairs_path, sep="\t")
    target_tips = sorted(list(pairs["genome"]))
    assert len(target_tips) == 237, f"Expected 237 target genomes, found {len(target_tips)}"
    assert len(set(target_tips)) == 237, "Duplicate genome IDs in v2_thermo_relatives_pairs.tsv"

    n_faa_present = 0
    missing_faa_genomes = []
    for g in target_tips:
        faa_path = os.path.join(faa_dir, f"{g}.faa")
        if os.path.exists(faa_path):
            sz = os.path.getsize(faa_path)
            chksum = sha256_file(faa_path)
            manifest_rows.append({
                "input_id": f"checkm2_faa__{g}",
                "file_path": faa_path,
                "host": "ibex",
                "access_status": "verified_readable",
                "size_bytes": sz,
                "sha256_checksum": chksum,
                "description": f"CheckM2 single-genome Pyrodigal protein FASTA for {g}"
            })
            n_faa_present += 1
        else:
            missing_faa_genomes.append(g)
            manifest_rows.append({
                "input_id": f"checkm2_faa__{g}",
                "file_path": faa_path,
                "host": "ibex",
                "access_status": "missing_in_checkm2_protein_files",
                "size_bytes": 0,
                "sha256_checksum": "NONE_MISSING_FILE",
                "description": f"Absent in _checkm2_refined/protein_files/ for {g} (excluded in S7 frac_complete sensitivity)"
            })

    print(f"[01_v3] Enumerated {len(core_inputs)} core inputs + {n_faa_present} readable CheckM2 .faa files + {len(missing_faa_genomes)} missing .faa files ({missing_faa_genomes})")

    manifest_df = pd.DataFrame(manifest_rows)
    manifest_tsv = os.path.join(out_dir, "v3_input_file_manifest.tsv")
    manifest_df.to_csv(manifest_tsv, sep="\t", index=False)
    # Also write directly to INPUTS.tsv so the top-level manifest has zero DYNAMIC/directory entries
    top_inputs_tsv = os.path.join(h1_dir, "INPUTS.tsv")
    manifest_df.to_csv(top_inputs_tsv, sep="\t", index=False)
    print(f"[01_v3] Saved 247-row enumerated manifest to {manifest_tsv} and {top_inputs_tsv}")

    # 3. Validate phylogenetic tree rooting (Midpoint vs Explicit Outgroup Rooting) and VCV properties
    import dendropy
    tree_path = os.path.join(ld_res, "thermo_tree.nwk")
    taxa_path = os.path.join(ld_res, "thermo_tree_taxa.tsv")
    taxa_df = pd.read_csv(taxa_path, sep="\t")
    outgroups = sorted(list(taxa_df[taxa_df["set"] == "outgroup"]["taxon"]))
    desert_from_taxa = sorted(list(taxa_df[taxa_df["set"] == "desert"]["taxon"]))
    assert desert_from_taxa == target_tips, "Mismatch between thermo_tree_taxa.tsv desert set and pairs table"
    assert len(outgroups) == 30, f"Expected 30 outgroup taxa, got {len(outgroups)}"

    # A. Midpoint-rooted extraction (as in 01_extract_phylogenetic_vcv.py)
    tree_mid = dendropy.Tree.get(path=tree_path, schema="newick", preserve_underscores=True)
    all_tips_mid = [leaf.taxon.label for leaf in tree_mid.leaf_nodes() if leaf.taxon]
    assert len(all_tips_mid) == 1914, f"Expected 1914 tips, got {len(all_tips_mid)}"
    assert len(set(all_tips_mid)) == 1914, "Duplicate tip labels in thermo_tree.nwk"
    assert set(target_tips).issubset(set(all_tips_mid)), "Missing target tips in thermo_tree.nwk"

    tree_mid.reroot_at_midpoint()
    tree_mid.retain_taxa_with_labels(target_tips)
    tree_mid.calc_node_root_distances()
    pdm_mid = tree_mid.phylogenetic_distance_matrix()
    taxa_mid = [tree_mid.taxon_namespace.get_taxon(t) for t in target_tips]
    rd_mid = np.array([tree_mid.find_node_with_taxon_label(t).root_distance for t in target_tips])

    n_t = len(target_tips)
    vcv_mid = np.zeros((n_t, n_t))
    dist_mid = np.zeros((n_t, n_t))
    for i in range(n_t):
        vcv_mid[i, i] = rd_mid[i]
        for j in range(i + 1, n_t):
            d_ij = pdm_mid(taxa_mid[i], taxa_mid[j])
            dist_mid[i, j] = dist_mid[j, i] = d_ij
            c_ij = 0.5 * (rd_mid[i] + rd_mid[j] - d_ij)
            vcv_mid[i, j] = vcv_mid[j, i] = c_ij

    std_mid = np.sqrt(np.diag(vcv_mid))
    corr_mid = vcv_mid / np.outer(std_mid, std_mid)
    np.fill_diagonal(corr_mid, 1.0)
    corr_mid = 0.5 * (corr_mid + corr_mid.T)
    corr_mid = np.clip(corr_mid, 0.0, 1.0)

    # B. Explicit Outgroup-rooted extraction (rooting on MRCA of the 30 Halobacteriota/Thermoproteota outgroups)
    tree_og = dendropy.Tree.get(path=tree_path, schema="newick", preserve_underscores=True)
    tree_og.is_rooted = True
    mrca_og = tree_og.mrca(taxon_labels=outgroups)
    tree_og.reroot_at_edge(mrca_og.edge, update_bipartitions=False)
    tree_og.retain_taxa_with_labels(target_tips)
    # Set stem above the ingroup crown node to 0.0 so root distances are measured from the desert Thermoplasmatota crown MRCA
    if tree_og.seed_node.edge.length is not None:
        tree_og.seed_node.edge.length = 0.0
    tree_og.calc_node_root_distances()
    pdm_og = tree_og.phylogenetic_distance_matrix()
    taxa_og = [tree_og.taxon_namespace.get_taxon(t) for t in target_tips]
    rd_og = np.array([tree_og.find_node_with_taxon_label(t).root_distance for t in target_tips])

    vcv_og = np.zeros((n_t, n_t))
    for i in range(n_t):
        vcv_og[i, i] = rd_og[i]
        for j in range(i + 1, n_t):
            d_ij = pdm_og(taxa_og[i], taxa_og[j])
            c_ij = 0.5 * (rd_og[i] + rd_og[j] - d_ij)
            vcv_og[i, j] = vcv_og[j, i] = c_ij

    std_og = np.sqrt(np.diag(vcv_og))
    corr_og = vcv_og / np.outer(std_og, std_og)
    np.fill_diagonal(corr_og, 1.0)
    corr_og = 0.5 * (corr_og + corr_og.T)
    corr_og = np.clip(corr_og, 0.0, 1.0)

    max_vcv_diff_mid_vs_og = float(np.max(np.abs(vcv_mid - vcv_og)))
    max_corr_diff_mid_vs_og = float(np.max(np.abs(corr_mid - corr_og)))
    assert max_vcv_diff_mid_vs_og < 1e-12, f"Midpoint vs outgroup VCV mismatch: {max_vcv_diff_mid_vs_og}"
    assert max_corr_diff_mid_vs_og < 1e-12, f"Midpoint vs outgroup Corr mismatch: {max_corr_diff_mid_vs_og}"

    # Verify saved phylo_correlation.tsv.gz and phylo_vcv.tsv.gz from Step 1
    saved_corr_df = pd.read_csv(os.path.join(out_dir, "phylo_correlation.tsv.gz"), sep="\t", index_col=0)
    saved_vcv_df = pd.read_csv(os.path.join(out_dir, "phylo_vcv.tsv.gz"), sep="\t", index_col=0)
    saved_corr = saved_corr_df.loc[target_tips, target_tips].values
    saved_vcv = saved_vcv_df.loc[target_tips, target_tips].values

    # Note saved files were formatted to %.6f
    assert np.max(np.abs(saved_corr - corr_mid)) < 1e-6, "Saved phylo_correlation.tsv.gz deviates from freshly computed matrix"
    assert np.max(np.abs(saved_vcv - vcv_mid)) < 1e-6, "Saved phylo_vcv.tsv.gz deviates from freshly computed matrix"

    eigs_corr = np.linalg.eigvalsh(saved_corr)
    eigs_vcv = np.linalg.eigvalsh(saved_vcv)
    assert eigs_corr.min() > 0, f"Correlation matrix not positive definite: min eig = {eigs_corr.min()}"
    assert eigs_vcv.min() > 0, f"Raw VCV matrix not positive definite: min eig = {eigs_vcv.min()}"

    off_diag_zero_pairs = int(np.sum(saved_vcv[np.triu_indices(n_t, 1)] == 0.0))
    total_pairs = int(n_t * (n_t - 1) // 2)

    val_info = {
        "dendropy_version": dendropy.__version__,
        "n_full_tree_tips": len(all_tips_mid),
        "n_outgroup_tips": len(outgroups),
        "n_target_desert_tips": n_t,
        "all_tips_unique": True,
        "mapped_fraction": 1.0,
        "max_abs_diff_vcv_midpoint_vs_outgroup_rooting": max_vcv_diff_mid_vs_og,
        "max_abs_diff_corr_midpoint_vs_outgroup_rooting": max_corr_diff_mid_vs_og,
        "crown_split_zero_covariance_pairs": off_diag_zero_pairs,
        "total_off_diagonal_pairs": total_pairs,
        "crown_split_zero_fraction": off_diag_zero_pairs / total_pairs,
        "root_to_tip_distance_min": float(rd_mid.min()),
        "root_to_tip_distance_q25": float(np.percentile(rd_mid, 25)),
        "root_to_tip_distance_median": float(np.median(rd_mid)),
        "root_to_tip_distance_q75": float(np.percentile(rd_mid, 75)),
        "root_to_tip_distance_max": float(rd_mid.max()),
        "corr_matrix_min_eigenvalue": float(eigs_corr.min()),
        "corr_matrix_max_eigenvalue": float(eigs_corr.max()),
        "corr_matrix_condition_number": float(eigs_corr.max() / eigs_corr.min()),
        "raw_vcv_min_eigenvalue": float(eigs_vcv.min()),
        "raw_vcv_max_eigenvalue": float(eigs_vcv.max()),
        "checkm2_faa_files_verified": n_faa_present,
        "checkm2_faa_files_missing": len(missing_faa_genomes),
        "checkm2_missing_genomes": missing_faa_genomes,
        "manifest_total_rows": len(manifest_df)
    }

    val_json = os.path.join(out_dir, "v3_phylo_vcv_validation.json")
    with open(val_json, "w") as f:
        json.dump(val_info, f, indent=2)

    lines = [
        "=" * 80,
        "V3 PHYLOGENETIC COVARIANCE AND INPUT MANIFEST VALIDATION REPORT",
        "=" * 80,
        f"DendroPy version: {dendropy.__version__}",
        f"Full ar53 tree tips: {len(all_tips_mid)} (all unique: True)",
        f"Target desert Thermoplasmatota tips mapped: {n_t} / {n_t} (100.0%)",
        f"Outgroup taxa (Halobacteriota + Thermoproteota): {len(outgroups)}",
        f"Midpoint vs Outgroup-rooted VCV max |diff|: {max_vcv_diff_mid_vs_og:.3e}",
        f"Midpoint vs Outgroup-rooted Correlation max |diff|: {max_corr_diff_mid_vs_og:.3e}",
        f"Crown split zero-covariance pairs (c__E2 vs c__SW-10-69-26 + c__Thermoplasmata): {off_diag_zero_pairs} / {total_pairs} ({off_diag_zero_pairs/total_pairs:.6f})",
        f"Correlation matrix C eigenvalues: min = {eigs_corr.min():.6f}, max = {eigs_corr.max():.6f}, cond = {eigs_corr.max()/eigs_corr.min():.2f}",
        f"Unnormalized VCV eigenvalues: min = {eigs_vcv.min():.6f}, max = {eigs_vcv.max():.6f}",
        f"Enumerated input manifest rows: {len(manifest_df)} (10 core files + {n_faa_present} readable CheckM2 .faa + {len(missing_faa_genomes)} missing CheckM2 .faa)",
        f"Missing CheckM2 .faa genomes (n=5): {', '.join(missing_faa_genomes)}",
        "=" * 80
    ]
    val_txt = os.path.join(out_dir, "v3_phylo_vcv_validation.txt")
    with open(val_txt, "w") as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    print(f"[01_v3] Completed in {time.time()-t0:.2f}s")


if __name__ == "__main__":
    main()
