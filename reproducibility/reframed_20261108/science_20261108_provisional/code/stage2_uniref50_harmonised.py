#!/usr/bin/env python3
"""Stage 2: Harmonised Direct UniRef50 (2026_01) Search Analysis across 538 Thermoplasmatota Genomes.

Implements Supervisor Review Item 8:
- Primary direct-search darkness matches active `main.tex` Methods:
  MMseqs2 against UniRef50 2026_01, sensitivity `-s 7.5`, `E <= 1e-3`, `--max-accept 10`,
  NO extra identity or coverage gate (`primary_direct_uniref50_dark = no returned hit`).
- Keeps inherited catalogue darkness (`archaeal_genes.tsv` `dark` / `dark_strict`) strictly
  distinct from `primary_direct_uniref50_dark` and outputs a full cross-tabulation.
- Retains stricter thresholds (`E <= 1e-5` and `pident < 30%`) only as clearly named
  sensitivity columns (`sensitivity_direct_uniref50_e1e5_dark`, `sensitivity_direct_uniref50_pident30_dark`).
- Deduplicates the 54 external genome aliases (`247` distinct external assemblies) for all
  statistical inference across EQ vs external comparators.
"""
from __future__ import annotations

import collections
import gzip
import json
import math
import os
import re
import sys
import tarfile
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common_science import (
    audit_orf_integrity,
    block_bootstrap_stat,
    build_stage2_per_genome_summary,
    classify_genome_quality,
    identify_mgnify_gtdb_duplicates,
    normalize_genus_groups,
    require_file,
    sha256_file,
    validate_primary_key,
    validated_merge,
)

P_ANNOT = Path("/project/k10018/work/eq-annot-20261030")
P_ARCH = Path("/scratch/project/k10018/eq-archaeal-function-20261106")
P_UNIREF = Path("/project/k10018/sw/db/eq-c2/uniref50")
P_OUT_ROOT = Path(os.environ.get("EQ_SCIENCE_OUT", "/scratch/project/k10018/eq-science-20261108"))
P_IN = P_OUT_ROOT / "inputs"
P_RES = P_OUT_ROOT / "results"


def parse_gtdb_lineage(tax_str: str) -> dict[str, str]:
    ranks = {"d": "", "p": "", "c": "", "o": "", "f": "", "g": "", "s": ""}
    if not isinstance(tax_str, str):
        return ranks
    for part in tax_str.split(";"):
        part = part.strip()
        if "__" in part:
            k, v = part.split("__", 1)
            if k in ranks:
                ranks[k] = f"{k}__{v}" if v else f"{k}__"
    return ranks


def verify_uniref50_release() -> dict[str, str]:
    ver_file = require_file(P_UNIREF / "uniref50.version", "uniref50.version")
    src_file = require_file(P_UNIREF / "uniref50.source", "uniref50.source")
    ver_txt = ver_file.read_text(encoding="utf-8", errors="replace")
    if "2026_01" not in ver_txt:
        raise ValueError(f"Expected UniRef50 Release 2026_01 in {ver_file}, got:\n{ver_txt}")
    return {
        "uniref50_version_path": str(ver_file),
        "uniref50_version_sha256": sha256_file(ver_file),
        "uniref50_source_sha256": sha256_file(src_file),
        "uniref50_version_text": " ".join(ver_txt.split()),
    }


def parse_538_proteomes_metadata(tar_path: Path) -> pd.DataFrame:
    """Extract per-gene length and Prodigal partial flags across all 538 proteomes in thermo_proteomes.tar.gz."""
    require_file(tar_path, "thermo_proteomes.tar.gz")
    rows = []
    partial_re = re.compile(r"partial=([01]{2})")
    with tarfile.open(tar_path, "r:gz") as tf:
        for member in tf:
            if not member.isfile() or not member.name.endswith(".faa"):
                continue
            genome_id = Path(member.name).stem
            f = tf.extractfile(member)
            if f is None:
                continue
            cur_gene: str | None = None
            cur_partial: str = "unknown"
            seq_len = 0
            for raw_line in f:
                line = raw_line.decode("utf-8", errors="replace").strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if cur_gene is not None:
                        rows.append((genome_id, cur_gene, f"{genome_id}|{cur_gene}", seq_len, cur_partial))
                    hdr = line[1:]
                    cur_gene = hdr.split()[0]
                    m = partial_re.search(hdr)
                    cur_partial = m.group(1) if m else "unknown"
                    seq_len = 0
                else:
                    seq_len += len(line.replace("*", ""))
            if cur_gene is not None:
                rows.append((genome_id, cur_gene, f"{genome_id}|{cur_gene}", seq_len, cur_partial))
    df = pd.DataFrame(rows, columns=["genome", "gene", "query_id", "aa_length", "prodigal_partial"])
    validate_primary_key(df, ["query_id"], "thermo538_proteomes_genes")
    return df


def parse_mmseqs_m8(m8_gz_path: Path) -> dict[str, tuple[float, float, str]]:
    """Parse MMseqs2 m8 output (query, target, pident, alnlen, evalue, bits) returning query_id -> (min_evalue, max_pident_at_e1e3, best_target).

    Note: m8 file was produced with `-s 7.5 -e 1e-3 --max-accept 10`, so every row in the file already satisfies `E <= 1e-3`.
    """
    require_file(m8_gz_path, "MMseqs2 UniRef50 m8 output")
    hits: dict[str, tuple[float, float, str]] = {}
    with gzip.open(m8_gz_path, "rt") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            qid, tgt, pid_s, _, ev_s, _ = parts[:6]
            ev = float(ev_s)
            pid = float(pid_s)
            if ev > 1e-3:
                continue
            prev = hits.get(qid)
            if prev is None:
                hits[qid] = (ev, pid, tgt)
            else:
                min_ev = min(prev[0], ev)
                max_pid = max(prev[1], pid)
                best_t = tgt if ev < prev[0] else prev[2]
                hits[qid] = (min_ev, max_pid, best_t)
    return hits


def main() -> None:
    m8_path = Path(sys.argv[1]) if len(sys.argv) > 1 else (P_RES / "thermo538_vs_uniref50_s75_e1e3.m8.gz")
    P_RES.mkdir(parents=True, exist_ok=True)

    print("[Stage 2] Verifying UniRef50 2026_01 database metadata and hashes...", flush=True)
    uniref_meta = verify_uniref50_release()
    tar_path = require_file(P_ANNOT / "results/c7_thermo_evol/thermo_proteomes.tar.gz", "thermo_proteomes.tar.gz")
    uniref_meta["thermo_proteomes_tar_sha256"] = sha256_file(tar_path)

    print("[Stage 2] Parsing 538 Thermoplasmatota proteomes from thermo_proteomes.tar.gz...", flush=True)
    prot_df = parse_538_proteomes_metadata(tar_path)
    print(f"[Stage 2] Loaded {len(prot_df)} total proteins across {prot_df['genome'].nunique()} genomes.", flush=True)

    print(f"[Stage 2] Parsing MMseqs2 UniRef50 search results from {m8_path}...", flush=True)
    m8_hits = parse_mmseqs_m8(m8_path)
    print(f"[Stage 2] Found {len(m8_hits)} queries with >=1 UniRef50 2026_01 hit at E <= 1e-3 (-s 7.5, max-accept 10).", flush=True)

    # Annotate primary direct-search darkness (NO hit at E <= 1e-3, NO extra identity/coverage gate)
    # and separate named sensitivities
    q_list = prot_df["query_id"].tolist()
    has_hit_e1e3 = np.zeros(len(q_list), dtype=bool)
    has_hit_e1e5 = np.zeros(len(q_list), dtype=bool)
    has_hit_pident30 = np.zeros(len(q_list), dtype=bool)
    min_ev_arr = np.full(len(q_list), np.nan, dtype="float64")
    max_pid_arr = np.full(len(q_list), np.nan, dtype="float64")

    for i, qid in enumerate(q_list):
        h = m8_hits.get(qid)
        if h is not None:
            has_hit_e1e3[i] = True
            min_ev_arr[i] = h[0]
            max_pid_arr[i] = h[1]
            if h[0] <= 1e-5:
                has_hit_e1e5[i] = True
            if h[1] >= 30.0:
                has_hit_pident30[i] = True

    prot_df["primary_direct_uniref50_dark"] = ~has_hit_e1e3
    prot_df["sensitivity_direct_uniref50_e1e5_dark"] = ~has_hit_e1e5
    prot_df["sensitivity_direct_uniref50_pident30_dark"] = ~has_hit_pident30
    prot_df["direct_uniref50_min_evalue"] = min_ev_arr
    prot_df["direct_uniref50_max_pident"] = max_pid_arr
    prot_df["len_bin"] = np.where(
        prot_df["aa_length"] < 100,
        "<100",
        np.where(prot_df["aa_length"] < 250, "100-249", ">=250"),
    )

    # =========================================================================
    # 1. Explicit Cross-Tabulation on the 237 EQ Genomes:
    #    Inherited Catalogue Darkness (`archaeal_genes.tsv`) vs Direct UniRef50 Search Darkness
    # =========================================================================
    print("[Stage 2] Computing cross-tabulation of inherited catalogue darkness vs direct UniRef50 search darkness on 237 EQ genomes...", flush=True)
    thermo_df = pd.read_csv(require_file(P_ARCH / "inputs/thermoplasmatota.tsv", "thermoplasmatota.tsv"), sep="\t")
    thermo_genomes_set = set(thermo_df["genome"])

    arch_genes_df = pd.read_csv(require_file(P_ARCH / "inputs/archaeal_genes.tsv", "archaeal_genes.tsv"), sep="\t", low_memory=False)
    arch_genes_df = arch_genes_df[arch_genes_df["genome"].isin(thermo_genomes_set)].copy().reset_index(drop=True)
    arch_genes_df["start_num"] = pd.to_numeric(arch_genes_df["start"], errors="coerce")
    arch_genes_df["end_num"] = pd.to_numeric(arch_genes_df["end"], errors="coerce")
    arch_genes_df = arch_genes_df.sort_values(["genome", "contig", "start_num", "end_num", "gene"]).reset_index(drop=True)
    arch_genes_df["contig_rank"] = arch_genes_df.groupby(["genome", "contig"]).cumcount() + 1
    arch_genes_df["contig_n_genes"] = arch_genes_df.groupby(["genome", "contig"])["gene"].transform("count")
    arch_genes_df["has_full_seq"] = True
    arch_genes_df["coords_valid"] = True
    arch_genes_df = audit_orf_integrity(
        arch_genes_df,
        partial_col="partial",
        internal_stop_col="internal_stop",
        length_col="length",
        contig_rank_col="contig_rank",
        contig_n_genes_col="contig_n_genes",
        contig_in_cat_col="contig_in_catalogue",
        has_full_seq_col="has_full_seq",
        coords_valid_col="coords_valid",
        edge_gene_buffer=2,
    )

    eq_prot_df = prot_df[prot_df["genome"].isin(thermo_genomes_set)].copy()
    eq_merged = validated_merge(
        arch_genes_df,
        eq_prot_df[
            [
                "genome",
                "gene",
                "primary_direct_uniref50_dark",
                "sensitivity_direct_uniref50_e1e5_dark",
                "sensitivity_direct_uniref50_pident30_dark",
                "direct_uniref50_min_evalue",
                "direct_uniref50_max_pident",
            ]
        ],
        on=["genome", "gene"],
        how="inner",
        validate="one_to_one",
        require_all_left=True,
        left_name="archaeal_genes_eq237",
        right_name="thermo_proteomes_eq237",
    )

    def _cat_dark_label(v: Any) -> str:
        if pd.isna(v):
            return "unlabelled_in_catalogue"
        fv = float(v)
        if fv == 1.0:
            return "inherited_catalogue_dark"
        if fv == 0.0:
            return "inherited_catalogue_known"
        return "unlabelled_in_catalogue"

    eq_merged["inherited_catalogue_status"] = [_cat_dark_label(x) for x in eq_merged["dark"]]
    eq_merged["inherited_catalogue_strict_status"] = [_cat_dark_label(x) for x in eq_merged["dark_strict"]]

    ct_rows = []
    for subset_name, sub_df in [
        ("all_777177_eq_genes", eq_merged),
        (
            "complete_interior_orfs_only",
            eq_merged[(eq_merged["orf_status"] == "complete") & (eq_merged["contig_edge_status"] == "interior")],
        ),
        (
            "identical_locus_in_catalogue",
            eq_merged[eq_merged["locus_class"] == "identical"],
        ),
        (
            "non_identical_or_recalled_locus",
            eq_merged[eq_merged["locus_class"] != "identical"],
        ),
    ]:
        for cat_status, grp in sub_df.groupby("inherited_catalogue_status", sort=False):
            n_tot = len(grp)
            n_dir_dark = int(grp["primary_direct_uniref50_dark"].sum())
            n_dir_known = n_tot - n_dir_dark
            n_e1e5_dark = int(grp["sensitivity_direct_uniref50_e1e5_dark"].sum())
            n_pid30_dark = int(grp["sensitivity_direct_uniref50_pident30_dark"].sum())
            ct_rows.append(
                {
                    "subset": subset_name,
                    "inherited_catalogue_status": cat_status,
                    "n_genes": n_tot,
                    "n_primary_direct_uniref50_dark_e1e3": n_dir_dark,
                    "n_primary_direct_uniref50_known_e1e3": n_dir_known,
                    "pct_primary_direct_uniref50_dark_e1e3": 100.0 * n_dir_dark / n_tot if n_tot else np.nan,
                    "n_sensitivity_direct_e1e5_dark": n_e1e5_dark,
                    "pct_sensitivity_direct_e1e5_dark": 100.0 * n_e1e5_dark / n_tot if n_tot else np.nan,
                    "n_sensitivity_direct_pident30_dark": n_pid30_dark,
                    "pct_sensitivity_direct_pident30_dark": 100.0 * n_pid30_dark / n_tot if n_tot else np.nan,
                }
            )
    pd.DataFrame(ct_rows).to_csv(P_RES / "task_b_uniref50_darkness_crosstab_eq237.tsv", sep="\t", index=False)

    # =========================================================================
    # 2. Harmonised Direct-Search Darkness Comparison: 237 EQ vs 247 Deduplicated External Comparators
    # =========================================================================
    print("[Stage 2] Running harmonised direct UniRef50 darkness comparison (237 EQ vs 247 deduplicated external comparators)...", flush=True)
    comp_df = pd.read_csv(require_file(P_IN / "v2_thermo_relatives_comparators.tsv", "v2_thermo_relatives_comparators.tsv"), sep="\t")
    skani_df = pd.read_csv(require_file(P_IN / "v2_thermo_relatives_skani.tsv.gz", "v2_thermo_relatives_skani.tsv.gz"), sep="\t")
    dup_ext_ids = identify_mgnify_gtdb_duplicates(comp_df, skani_df)
    comp_df["is_duplicate_of_gtdb_rep"] = comp_df["genome_id"].isin(dup_ext_ids)
    comp_dedup = comp_df[~comp_df["is_duplicate_of_gtdb_rep"]].copy().reset_index(drop=True)
    assert len(comp_dedup) == 247, f"Expected 247 deduplicated external comparators, got {len(comp_dedup)}"

    # Build unified 484-genome metadata table (237 EQ + 247 deduplicated external comparators)
    eq_meta = pd.DataFrame(
        {
            "genome": thermo_df["genome"],
            "cohort": "EQ_Thermoplasmatota_237",
            "lineage": thermo_df["gtdb_classification"],
            "completeness": pd.to_numeric(thermo_df["checkm_completeness"], errors="coerce"),
            "contamination": pd.to_numeric(thermo_df["contamination"], errors="coerce"),
        }
    )
    ext_meta = pd.DataFrame(
        {
            "genome": comp_dedup["genome_id"],
            "cohort": "External_Comparators_Deduplicated_247",
            "lineage": comp_dedup["lineage"],
            "completeness": pd.to_numeric(comp_dedup["completeness"], errors="coerce"),
            "contamination": pd.to_numeric(comp_dedup["contamination"], errors="coerce"),
        }
    )
    meta_484 = pd.concat([eq_meta, ext_meta], ignore_index=True)
    validate_primary_key(meta_484, ["genome"], "unified_484_genomes")

    tax_484 = [parse_gtdb_lineage(x) for x in meta_484["lineage"]]
    for rk, col_name in [("c", "tax_class"), ("o", "tax_order"), ("f", "tax_family"), ("g", "tax_genus")]:
        meta_484[col_name] = [d[rk] for d in tax_484]
    g_blocks_484, _ = normalize_genus_groups(meta_484["genome"].tolist(), meta_484["tax_genus"].tolist())
    meta_484["genus_singleton_block"] = g_blocks_484
    meta_484["quality_known"] = meta_484["completeness"].notna() & meta_484["contamination"].notna()
    meta_484["comp70_contam5_pass"] = (
        meta_484["quality_known"]
        & (meta_484["completeness"] >= 70.0)
        & (meta_484["contamination"] <= 5.0)
    )

    # Identify AAI >= 65% matched genome subset
    pairs_df = pd.read_csv(require_file(P_IN / "v2_thermo_relatives_pairs.tsv", "v2_thermo_relatives_pairs.tsv"), sep="\t")
    pairs_aai65 = pairs_df[pd.to_numeric(pairs_df["aai"], errors="coerce").fillna(0.0) >= 65.0]
    aai65_eq_genomes = set(pairs_aai65["genome"])
    aai65_ext_genomes = set(pairs_aai65["aai_nearest"].dropna()) | set(pairs_aai65["aai_nearest_gtdb"].dropna())

    # Join per-gene direct UniRef50 results onto the 484 deduplicated genomes
    prot_484 = prot_df[prot_df["genome"].isin(set(meta_484["genome"]))].copy()
    # Map EQ prodigal_partial from archaeal_genes_df (`partial == '00'`)
    eq_partial_map = dict(zip(arch_genes_df["gene"], arch_genes_df["orf_status"] == "complete"))
    is_complete_orf = []
    for g_id, gene_id, p_flag in zip(prot_484["genome"], prot_484["gene"], prot_484["prodigal_partial"]):
        if g_id in thermo_genomes_set:
            is_complete_orf.append(bool(eq_partial_map.get(gene_id, False)))
        else:
            is_complete_orf.append(p_flag == "00")
    prot_484["is_complete_orf"] = is_complete_orf

    prot_484 = validated_merge(
        prot_484,
        meta_484,
        on="genome",
        how="left",
        validate="many_to_one",
        require_all_left=True,
        left_name="prot_484",
        right_name="meta_484",
    )

    ext_man_df: pd.DataFrame | None = None
    for cand_man in (P_RES / "external_signalp_manifest.tsv", P_IN / "external_signalp_manifest.tsv"):
        if cand_man.is_file():
            ext_man_df = pd.read_csv(cand_man, sep="\t")
            break

    # Per-genome summary table across all 484 genomes (preserving unknown completeness/contamination as NaN)
    per_genome_df = build_stage2_per_genome_summary(
        prot_484,
        meta_484,
        expected_n_eq=237,
        expected_n_ext=247,
        ext_manifest_df=ext_man_df,
    )
    assert len(per_genome_df) == 484, f"Expected 484 genomes in per_genome_df, got {len(per_genome_df)}"
    per_genome_df.to_csv(P_RES / "task_b_uniref50_harmonised_per_genome_484.tsv", sep="\t", index=False)

    # Compare EQ vs External across matched strata with genus-singleton block bootstrap
    shared_orders = {
        o for o in per_genome_df[per_genome_df["cohort"] == "EQ_Thermoplasmatota_237"]["tax_order"] if o != "o__"
    } & {
        o
        for o in per_genome_df[per_genome_df["cohort"] == "External_Comparators_Deduplicated_247"]["tax_order"]
        if o != "o__"
    }
    strata_defs = [
        ("all_deduplicated_484_genomes", per_genome_df),
        (
            "comp70_contam5_genomes",
            per_genome_df[per_genome_df["quality_known"] & per_genome_df["comp70_contam5_pass"]],
        ),
        ("class_c__SW-10-69-26_only", per_genome_df[per_genome_df["tax_class"] == "c__SW-10-69-26"]),
        ("class_c__E2_only", per_genome_df[per_genome_df["tax_class"] == "c__E2"]),
        ("shared_taxonomic_orders_only", per_genome_df[per_genome_df["tax_order"].isin(shared_orders)]),
        (
            "aai_ge65_matched_relatives_only",
            per_genome_df[per_genome_df["genome"].isin(aai65_eq_genomes | aai65_ext_genomes)],
        ),
    ]

    harm_rows = []
    for s_name, s_df in strata_defs:
        eq_sub = s_df[s_df["cohort"] == "EQ_Thermoplasmatota_237"]
        ext_sub = s_df[s_df["cohort"] == "External_Comparators_Deduplicated_247"]
        if len(eq_sub) < 2 or len(ext_sub) < 2:
            continue

        for metric_col, gene_set_label in [
            ("pct_primary_direct_dark_all", "all_orfs"),
            ("pct_primary_direct_dark_complete", "complete_orfs_only"),
        ]:
            sub_valid = s_df[s_df[metric_col].notna()].copy()
            is_eq = (sub_valid["cohort"] == "EQ_Thermoplasmatota_237").to_numpy(dtype=bool)
            vals = sub_valid[metric_col].to_numpy(dtype="float64")

            def _diff_eq_minus_ext(idx: np.ndarray) -> float:
                m_eq = is_eq[idx]
                m_ext = ~m_eq
                if not m_eq.any() or not m_ext.any():
                    return np.nan
                return float(vals[idx][m_eq].mean() - vals[idx][m_ext].mean())

            boot = block_bootstrap_stat(
                sub_valid["genus_singleton_block"].tolist(),
                _diff_eq_minus_ext,
                n_boot=1000,
                seed=20261108,
            )
            harm_rows.append(
                {
                    "stratum": s_name,
                    "orf_subset": gene_set_label,
                    "darkness_definition": "primary_direct_uniref50_2026_01_s7.5_e1e-3_max10_no_id_cov_gate",
                    "harmonisation_scope": "database_and_search_rule_harmonised_UniRef50_2026_01_not_gene_caller_mode_harmonised",
                    "eq_gene_caller_provenance": "eq_representative_proteomes_metagenomic_mode_p_meta",
                    "ext_gene_caller_provenance": "c7b_orthofinder_sbatch_prodigal_p_single_fallback_p_meta",
                    "n_eq_genomes": int(is_eq.sum()),
                    "n_ext_dedup_genomes": int((~is_eq).sum()),
                    "n_eq_genus_blocks": int(sub_valid[is_eq]["genus_singleton_block"].nunique()),
                    "n_ext_genus_blocks": int(sub_valid[~is_eq]["genus_singleton_block"].nunique()),
                    "eq_mean_pct_dark": float(vals[is_eq].mean()),
                    "ext_mean_pct_dark": float(vals[~is_eq].mean()),
                    "diff_eq_minus_ext_pp": float(boot["estimate"]),
                    "diff_se_genus_blocked": float(boot["se"]),
                    "diff_ci95_low": float(boot["ci95_low"]),
                    "diff_ci95_high": float(boot["ci95_high"]),
                }
            )

    pd.DataFrame(harm_rows).to_csv(P_RES / "task_b_uniref50_harmonised_eq_vs_external.tsv", sep="\t", index=False)

    # Save per-gene direct UniRef50 darkness table across all 484 deduplicated genomes so IBEX external SignalP 6.0h
    # predictions can be joined directly by (genome, gene) / query_id without re-running MMseqs2
    prot_484[
        [
            "genome",
            "gene",
            "query_id",
            "cohort",
            "tax_class",
            "tax_order",
            "tax_family",
            "tax_genus",
            "genus_singleton_block",
            "completeness",
            "contamination",
            "quality_known",
            "comp70_contam5_pass",
            "aa_length",
            "len_bin",
            "prodigal_partial",
            "is_complete_orf",
            "primary_direct_uniref50_dark",
            "sensitivity_direct_uniref50_e1e5_dark",
            "sensitivity_direct_uniref50_pident30_dark",
            "direct_uniref50_min_evalue",
            "direct_uniref50_max_pident",
        ]
    ].to_csv(P_RES / "task_b_uniref50_direct_per_gene_484.tsv.gz", sep="\t", index=False, compression="gzip")

    eq_pg = per_genome_df[per_genome_df["cohort"] == "EQ_Thermoplasmatota_237"]
    ext_pg = per_genome_df[per_genome_df["cohort"] == "External_Comparators_Deduplicated_247"]
    ext_pg_mgyg = ext_pg[ext_pg["quality_known"]]
    ext_pg_gtdb = ext_pg[~ext_pg["quality_known"]]

    stage2_summary = {
        "stage": "stage2_uniref50_harmonised",
        "status": "completed",
        "harmonisation_scope": "database_and_search_rule_harmonised_UniRef50_2026_01_not_gene_caller_mode_harmonised",
        "gene_calling_provenance": {
            "eq_237_genomes": "extracted in scripts/c7b_orthofinder.sbatch from eq-annot-20261030/inputs/proteomes.faa (representative MAG proteomes called in metagenomic mode -p meta)",
            "external_247_dedup_comparators": "called in scripts/c7b_orthofinder.sbatch from eq-annot-20261030/inputs/thermo_comparators.tar using prodigal -p single -q (with fallback to prodigal -p meta -q), terminal * stripped",
            "external_genomic_fasta_archive_on_ibex": "/ibex/scratch/projects/c2014/rob/metagenomics/eq-annot-20261030/inputs/thermo_comparators.tar (or /ibex/scratch/projects/c2014/rob/metagenomics/v2_thermo_relatives/)",
            "confounding_warning": "Database and MMseqs2 search parameters are identical across all 484 genomes, but gene-calling mode differs (-p meta on EQ MAGs vs prodigal -p single on external comparators). Do not claim fully gene-call-harmonised darkness without recalling both cohorts under identical gene-caller settings.",
        },
        "uniref50_metadata": uniref_meta,
        "n_total_genomes_538": int(prot_df["genome"].nunique()),
        "n_deduplicated_genomes_484": int(len(per_genome_df)),
        "n_eq_genomes": int(len(eq_pg)),
        "n_ext_dedup_genomes": int(len(ext_pg)),
        "n_ext_dedup_genomes_quality_known_mgyg": int(len(ext_pg_mgyg)),
        "n_ext_dedup_genomes_quality_unknown_gtdb": int(len(ext_pg_gtdb)),
        "n_total_proteins_538": int(len(prot_df)),
        "n_proteins_with_e1e3_hit": int(len(m8_hits)),
        "n_total_proteins_484_deduplicated": int(len(prot_484)),
        "n_eq237_proteins_all": int(eq_pg["n_genes_all"].sum()),
        "n_eq237_proteins_complete_orfs": int(eq_pg["n_genes_complete"].sum()),
        "n_ext301_raw_proteins_all": int(len(prot_df) - len(eq_prot_df)),
        "n_ext247_dedup_proteins_all": int(ext_pg["n_genes_all"].sum()),
        "n_ext247_dedup_proteins_complete_orfs": int(ext_pg["n_genes_complete"].sum()),
        "n_ext29_mgyg_quality_known_proteins_all": int(ext_pg_mgyg["n_genes_all"].sum()),
        "n_ext218_gtdb_quality_unknown_proteins_all": int(ext_pg_gtdb["n_genes_all"].sum()),
        "eq237_pct_primary_direct_dark_all": float(100.0 * eq_merged["primary_direct_uniref50_dark"].mean()),
        "eq237_pct_inherited_catalogue_dark_labelled": float(
            100.0
            * (
                eq_merged[eq_merged["inherited_catalogue_status"] != "unlabelled_in_catalogue"][
                    "inherited_catalogue_status"
                ]
                == "inherited_catalogue_dark"
            ).mean()
        ),
    }
    (P_RES / "stage2_summary.json").write_text(json.dumps(stage2_summary, indent=2))
    print("[Stage 2] Completed successfully!", json.dumps(stage2_summary), flush=True)


if __name__ == "__main__":
    main()
