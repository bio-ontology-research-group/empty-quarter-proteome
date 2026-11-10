#!/usr/bin/env python3
"""Stage 2b: Join Canonical EQ Taxonomy, Stage 2 Direct UniRef50 Darkness, and IBEX External Comparator SignalP 6.0h Predictions.

Implements:
- `external-signalp-design.txt`
- `external-comparison-scope-review.txt`
- `external-signalp-final-checks.txt`
- `canonical-taxonomy-sensitivity-design.txt` (Stage 2b section):
  * Joins the 237 EQ Thermoplasmatota genomes to canonical `refined_representatives.tsv` (`mag == genome`)
    with SHA-256 (`dbe9db11362e77c4274543431ad40ee92ed6b181c2c0f7ddc055c4f5a32971d6`), 1:1 uniqueness,
    exact 237 membership, `p__Thermoplasmatota` phylum, and completeness/contamination checks (`42/237` strings differ).
  * Retains legacy analytical EQ taxonomy (`thermoplasmatota.tsv`) as an explicit sensitivity analysis
    (`taxonomy_scheme = 'legacy_analytical_taxonomy_sensitivity'`) alongside the primary canonical GTDB R232
    data-foundation taxonomy (`taxonomy_scheme = 'canonical_gtdb_r232_data_foundation'`).
  * Preserves external comparator taxonomy (`External_Comparators_Deduplicated_247`) unchanged from
    `thermo_comparators_301.tsv` (`external_taxonomy_provenance = 'unchanged_from_thermo_comparators_301_tsv'`).
  * Labels genus-or-unclassified-singleton blocks explicitly as resampling blocks, never as independent
    phylogenetic lineages.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common_science import (
    CANONICAL_REFINED_REPRESENTATIVES_SHA256,
    block_bootstrap_stat,
    compute_thermo_proteomes_tar_sha256_map,
    join_and_validate_canonical_eq_taxonomy,
    require_file,
    validate_full_stage2b_signalp_and_uniref50_provenance,
    validate_primary_key,
)
from stage1_audit_and_export import load_all_sp_tm
from stage1b_refined_fasta_and_corrections import len_bin_7, mantel_haenszel_with_rbg_ci

P_OUT_ROOT = Path(os.environ.get("EQ_SCIENCE_OUT", "/scratch/project/k10018/eq-science-20261108"))
P_IN = P_OUT_ROOT / "inputs"
P_RES = P_OUT_ROOT / "results"


def run_eq237_canonical_taxonomy_reconciliation(
    thermo_tsv_path: Path,
    canonical_tsv_path: Path,
    out_res_dir: Path,
) -> tuple[pd.DataFrame, dict[str, object]]:
    """Reconcile the 237 EQ Thermoplasmatota genomes against canonical `refined_representatives.tsv`."""
    eq_df = pd.read_csv(require_file(thermo_tsv_path, "thermoplasmatota.tsv"), sep="\t")
    joined, summary = join_and_validate_canonical_eq_taxonomy(
        eq_df,
        canonical_tsv_path,
        genome_col="genome",
        legacy_taxonomy_col="gtdb_classification",
        expected_sha256=CANONICAL_REFINED_REPRESENTATIVES_SHA256,
        expected_n_eq=237,
        require_canonical_diff_count=42,
    )
    out_res_dir.mkdir(parents=True, exist_ok=True)
    reconcil_p = out_res_dir / "task_b_stage2b_eq237_canonical_vs_legacy_taxonomy_reconciliation.tsv"
    summary_p = out_res_dir / "task_b_stage2b_eq237_canonical_vs_legacy_taxonomy_summary.json"
    joined.to_csv(reconcil_p, sep="\t", index=False)
    with open(summary_p, "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2, sort_keys=True)
    return joined, summary


def evaluate_stage2b_uniref50_and_signalp_models(
    df_484: pd.DataFrame,
    eq_reconciled: pd.DataFrame,
    tax_summary: dict[str, object],
    pairs_df: pd.DataFrame,
    *,
    has_signalp_joined: bool,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Evaluate EQ-vs-External direct UniRef50 darkness and (when present) SignalP 6.0h export across both canonical and legacy EQ taxonomy schemes."""
    # Map EQ 237 canonical and legacy order/genus/block onto df_484 without altering External 247 taxonomy
    can_order_map = dict(zip(eq_reconciled["genome"].astype(str), eq_reconciled["canonical_order_label"].astype(str)))
    can_genus_map = dict(zip(eq_reconciled["genome"].astype(str), eq_reconciled["canonical_genus_label"].astype(str)))
    can_block_map = dict(zip(eq_reconciled["genome"].astype(str), eq_reconciled["canonical_genus_or_singleton_block"].astype(str)))

    leg_order_map = dict(zip(eq_reconciled["genome"].astype(str), eq_reconciled["legacy_order_label"].astype(str)))
    leg_genus_map = dict(zip(eq_reconciled["genome"].astype(str), eq_reconciled["legacy_genus_label"].astype(str)))
    leg_block_map = dict(zip(eq_reconciled["genome"].astype(str), eq_reconciled["legacy_genus_or_singleton_block"].astype(str)))

    eq_genomes_in_df = set(df_484.loc[df_484["cohort"] == "EQ_Thermoplasmatota_237", "genome"].astype(str).unique())
    if eq_genomes_in_df != set(can_order_map.keys()):
        raise ValueError(
            f"EQ 237 genome membership mismatch in df_484 vs canonical reconciliation: "
            f"n_df={len(eq_genomes_in_df)}, n_reconcil={len(can_order_map)}"
        )

    pairs_aai65 = pairs_df[pd.to_numeric(pairs_df["aai"], errors="coerce").fillna(0.0) >= 65.0]
    aai65_genomes = (
        set(pairs_aai65["genome"].astype(str))
        | set(pairs_aai65["aai_nearest"].dropna().astype(str))
        | set(pairs_aai65["aai_nearest_gtdb"].dropna().astype(str))
    )

    darkness_rows: list[dict[str, object]] = []
    signalp_rows: list[dict[str, object]] = []

    for scheme_name, order_map, genus_map, block_map in [
        ("canonical_gtdb_r232_data_foundation", can_order_map, can_genus_map, can_block_map),
        ("legacy_analytical_taxonomy_sensitivity", leg_order_map, leg_genus_map, leg_block_map),
    ]:
        work = df_484.copy()
        is_eq_row = (work["cohort"] == "EQ_Thermoplasmatota_237").to_numpy(dtype=bool)
        g_str = work["genome"].astype(str)

        work["scheme_tax_order"] = np.where(is_eq_row, g_str.map(order_map), work["tax_order"].astype(str))
        work["scheme_tax_genus"] = np.where(is_eq_row, g_str.map(genus_map), work["tax_genus"].astype(str))
        work["scheme_resampling_block"] = np.where(is_eq_row, g_str.map(block_map), work["genus_singleton_block"].astype(str))

        eq_orders_named = {
            o for o in work.loc[is_eq_row, "scheme_tax_order"].unique() if str(o) != "o__"
        }
        ext_orders_named = {
            o for o in work.loc[~is_eq_row, "scheme_tax_order"].unique() if str(o) != "o__"
        }
        shared_named_orders = eq_orders_named & ext_orders_named

        if "completeness" in work.columns and "contamination" in work.columns:
            qual_pass = (
                work["comp70_contam5_pass"].astype(bool)
                & work["completeness"].notna()
                & work["contamination"].notna()
            )
        else:
            qual_pass = work["comp70_contam5_pass"].astype(bool)

        strata = [
            ("all_484_deduplicated_genomes", work),
            ("comp70_contam5_genomes", work[qual_pass]),
            ("shared_named_taxonomic_orders_only", work[work["scheme_tax_order"].isin(shared_named_orders)]),
            ("aai_ge65_matched_relatives_only", work[work["genome"].astype(str).isin(aai65_genomes)]),
        ]

        for s_name, s_df in strata:
            for orf_label, sub_df in [
                ("all_orfs", s_df),
                ("complete_orfs_only", s_df[s_df["is_complete_orf"]]),
            ]:
                eq_sub = sub_df[sub_df["cohort"] == "EQ_Thermoplasmatota_237"]
                ext_sub = sub_df[sub_df["cohort"] == "External_Comparators_Deduplicated_247"]
                if len(eq_sub) == 0 or len(ext_sub) == 0:
                    continue

                # 1. Direct UniRef50 darkness per-genome comparison under this taxonomy scheme
                g_dark = (
                    sub_df.groupby(["genome", "cohort", "scheme_resampling_block"], as_index=False, sort=False)
                    .agg(
                        n_genes=("primary_direct_uniref50_dark", "size"),
                        n_dark=("primary_direct_uniref50_dark", "sum"),
                    )
                )
                g_dark["pct_dark"] = 100.0 * g_dark["n_dark"] / g_dark["n_genes"]
                is_eq_gd = (g_dark["cohort"] == "EQ_Thermoplasmatota_237").to_numpy(dtype=bool)
                pct_vals = g_dark["pct_dark"].to_numpy(dtype=float)

                def _diff_pct_dark(idx: np.ndarray) -> float:
                    m_e = is_eq_gd[idx]
                    if not m_e.any() or (~m_e).sum() == 0:
                        return np.nan
                    return float(pct_vals[idx][m_e].mean() - pct_vals[idx][~m_e].mean())

                boot_dark = block_bootstrap_stat(
                    g_dark["scheme_resampling_block"].tolist(), _diff_pct_dark, n_boot=1000, seed=20261108
                )
                darkness_rows.append(
                    {
                        "taxonomy_scheme": scheme_name,
                        "canonical_refined_representatives_sha256": str(tax_summary["canonical_tsv_sha256"]),
                        "external_taxonomy_provenance": "unchanged_from_thermo_comparators_301_tsv",
                        "block_semantics_note": "genus_or_unclassified_singleton_resampling_blocks_not_independent_phylogenetic_lineages",
                        "stratum": s_name,
                        "orf_subset": orf_label,
                        "harmonisation_scope": "database_and_search_rule_harmonised_not_gene_caller_mode_harmonised",
                        "n_eq_genomes": int(is_eq_gd.sum()),
                        "n_ext_genomes": int((~is_eq_gd).sum()),
                        "n_eq_resampling_blocks": int(g_dark.loc[is_eq_gd, "scheme_resampling_block"].nunique()),
                        "n_ext_resampling_blocks": int(g_dark.loc[~is_eq_gd, "scheme_resampling_block"].nunique()),
                        "eq_mean_genome_pct_dark": float(pct_vals[is_eq_gd].mean()),
                        "ext_mean_genome_pct_dark": float(pct_vals[~is_eq_gd].mean()),
                        "diff_mean_genome_pct_dark_eq_minus_ext": float(boot_dark["estimate"]),
                        "diff_se_genus_or_singleton_block": float(boot_dark["se"]),
                        "diff_ci95_low": float(boot_dark["ci95_low"]),
                        "diff_ci95_high": float(boot_dark["ci95_high"]),
                    }
                )

                # 2. SignalP 6.0h export by direct UniRef50 darkness under this taxonomy scheme
                if has_signalp_joined:
                    sp_sub = sub_df[sub_df["has_signalp6"]]
                    eq_sp = sp_sub[sp_sub["cohort"] == "EQ_Thermoplasmatota_237"]
                    ext_sp = sp_sub[sp_sub["cohort"] == "External_Comparators_Deduplicated_247"]
                    if len(eq_sp) == 0 or len(ext_sp) == 0:
                        continue

                    mh_eq_gen = mantel_haenszel_with_rbg_ci(
                        eq_sp, ["genome", "len_bin7"], "primary_direct_uniref50_dark", "sp_any"
                    )
                    mh_ext_gen = mantel_haenszel_with_rbg_ci(
                        ext_sp, ["genome", "len_bin7"], "primary_direct_uniref50_dark", "sp_any"
                    )

                    g_agg = (
                        sp_sub.groupby(
                            ["genome", "cohort", "scheme_resampling_block", "primary_direct_uniref50_dark"],
                            sort=False,
                        )
                        .agg(n_tot=("sp_any", "size"), n_sp=("sp_any", "sum"))
                        .unstack("primary_direct_uniref50_dark", fill_value=0)
                    )
                    g_agg.columns = [f"{c}_{'dark' if d else 'known'}" for c, d in g_agg.columns]
                    g_agg = g_agg.reset_index()
                    for col_req in ("n_tot_dark", "n_sp_dark", "n_tot_known", "n_sp_known"):
                        if col_req not in g_agg.columns:
                            g_agg[col_req] = 0
                    a = g_agg["n_sp_dark"].to_numpy(dtype=float) + 0.5
                    b = (g_agg["n_tot_dark"] - g_agg["n_sp_dark"]).to_numpy(dtype=float) + 0.5
                    c = g_agg["n_sp_known"].to_numpy(dtype=float) + 0.5
                    d = (g_agg["n_tot_known"] - g_agg["n_sp_known"]).to_numpy(dtype=float) + 0.5
                    g_agg["log_or_sp_dark_vs_known"] = np.log((a * d) / (b * c))
                    is_eq_g = (g_agg["cohort"] == "EQ_Thermoplasmatota_237").to_numpy(dtype=bool)
                    lor_vals = g_agg["log_or_sp_dark_vs_known"].to_numpy(dtype=float)

                    def _diff_lor(idx: np.ndarray) -> float:
                        m_e = is_eq_g[idx]
                        if not m_e.any() or (~m_e).sum() == 0:
                            return np.nan
                        return float(lor_vals[idx][m_e].mean() - lor_vals[idx][~m_e].mean())

                    boot_sp = block_bootstrap_stat(
                        g_agg["scheme_resampling_block"].tolist(), _diff_lor, n_boot=1000, seed=20261108
                    )
                    signalp_rows.append(
                        {
                            "taxonomy_scheme": scheme_name,
                            "canonical_refined_representatives_sha256": str(tax_summary["canonical_tsv_sha256"]),
                            "external_taxonomy_provenance": "unchanged_from_thermo_comparators_301_tsv",
                            "block_semantics_note": "genus_or_unclassified_singleton_resampling_blocks_not_independent_phylogenetic_lineages",
                            "stratum": s_name,
                            "orf_subset": orf_label,
                            "harmonisation_scope": "database_search_and_SignalP6_0h_harmonised_not_gene_caller_mode_harmonised",
                            "eq_gene_caller_mode": "metagenomic_mode_p_meta",
                            "ext_gene_caller_mode": "prodigal_p_single_fallback_p_meta_in_c7b_orthofinder_sbatch",
                            "n_eq_genomes": int(eq_sp["genome"].nunique()),
                            "n_ext_genomes": int(ext_sp["genome"].nunique()),
                            "n_eq_resampling_blocks": int(g_agg.loc[is_eq_g, "scheme_resampling_block"].nunique()),
                            "n_ext_resampling_blocks": int(g_agg.loc[~is_eq_g, "scheme_resampling_block"].nunique()),
                            "eq_dark_sp_pct": float(100.0 * eq_sp[eq_sp["primary_direct_uniref50_dark"]]["sp_any"].mean()),
                            "eq_known_sp_pct": float(100.0 * eq_sp[~eq_sp["primary_direct_uniref50_dark"]]["sp_any"].mean()),
                            "eq_within_genome_len7_mh_or": float(mh_eq_gen["or_mh"]),
                            "eq_within_genome_len7_mh_ci95_low": float(mh_eq_gen["or_ci95_low"]),
                            "eq_within_genome_len7_mh_ci95_high": float(mh_eq_gen["or_ci95_high"]),
                            "ext_dark_sp_pct": float(100.0 * ext_sp[ext_sp["primary_direct_uniref50_dark"]]["sp_any"].mean()),
                            "ext_known_sp_pct": float(100.0 * ext_sp[~ext_sp["primary_direct_uniref50_dark"]]["sp_any"].mean()),
                            "ext_within_genome_len7_mh_or": float(mh_ext_gen["or_mh"]),
                            "ext_within_genome_len7_mh_ci95_low": float(mh_ext_gen["or_ci95_low"]),
                            "ext_within_genome_len7_mh_ci95_high": float(mh_ext_gen["or_ci95_high"]),
                            "diff_mean_genome_log_or_eq_minus_ext": float(boot_sp["estimate"]),
                            "diff_mean_genome_log_or_se_genus_blocked": float(boot_sp["se"]),
                            "diff_mean_genome_log_or_ci95_low": float(boot_sp["ci95_low"]),
                            "diff_mean_genome_log_or_ci95_high": float(boot_sp["ci95_high"]),
                            "ratio_of_ors_eq_over_ext": float(math.exp(float(boot_sp["estimate"]))),
                        }
                    )

    return pd.DataFrame(darkness_rows), pd.DataFrame(signalp_rows)


def main() -> None:
    ap = argparse.ArgumentParser(description="Stage 2b: Canonical EQ taxonomy reconciliation + Direct UniRef50 + External SignalP 6.0h join")
    ap.add_argument("--reconcile-taxonomy-only", action="store_true", help="Run only the 237-genome EQ canonical vs legacy taxonomy reconciliation and exit.")
    ap.add_argument("--allow-pending-signalp", action="store_true", help="Compute canonical-vs-legacy direct UniRef50 darkness models if external247_signalp6_predictions.tsv.gz is still pending on IBEX.")
    args = ap.parse_args()

    thermo_p = P_IN / "thermoplasmatota.tsv"
    canonical_p = P_IN / "refined_representatives.tsv"
    eq_reconciled, tax_summary = run_eq237_canonical_taxonomy_reconciliation(thermo_p, canonical_p, P_RES)
    print(
        f"[Stage 2b] Verified canonical refined_representatives.tsv (sha256={tax_summary['canonical_tsv_sha256']}): "
        f"{tax_summary['n_taxonomy_strings_differ']}/{tax_summary['n_eq_genomes']} EQ taxonomy strings differ "
        f"(canonical: {tax_summary['canonical_n_named_orders']} named orders + {tax_summary['canonical_n_unassigned_order_genomes']} unassigned, "
        f"{tax_summary['canonical_n_named_genera']} named genera + {tax_summary['canonical_n_unassigned_genus_genomes']} unassigned).",
        flush=True,
    )
    if args.reconcile_taxonomy_only:
        return

    direct_p = require_file(
        P_RES / "task_b_uniref50_direct_per_gene_484.tsv.gz",
        "task_b_uniref50_direct_per_gene_484.tsv.gz",
    )
    pairs_df = pd.read_csv(
        require_file(P_IN / "v2_thermo_relatives_pairs.tsv", "v2_thermo_relatives_pairs.tsv"),
        sep="\t",
    )

    print("[Stage 2b] Loading per-gene direct UniRef50 darkness across 484 genomes...", flush=True)
    df_484 = pd.read_csv(direct_p, sep="\t", low_memory=False)
    validate_primary_key(df_484, ["genome", "gene"], "df_484_direct_uniref50")
    n_eq_g_484 = int(df_484.loc[df_484["cohort"] == "EQ_Thermoplasmatota_237", "genome"].nunique())
    n_ext_g_484 = int(df_484.loc[df_484["cohort"] == "External_Comparators_Deduplicated_247", "genome"].nunique())
    if n_eq_g_484 != 237 or n_ext_g_484 != 247 or int(df_484["genome"].nunique()) != 484:
        raise ValueError(
            f"df_484 genome membership mismatch: expected 237 EQ + 247 External = 484 genomes, "
            f"got {n_eq_g_484} EQ + {n_ext_g_484} External = {df_484['genome'].nunique()}"
        )
    df_484["len_bin7"] = len_bin_7(df_484["aa_length"]).to_numpy()

    ext_sig_cand = P_RES / "external247_signalp6_predictions.tsv.gz"
    ext_man_cand = P_RES / "external_signalp_manifest.tsv"
    if not ext_man_cand.is_file() and (P_IN / "external_signalp_manifest.tsv").is_file():
        ext_man_cand = P_IN / "external_signalp_manifest.tsv"

    has_ext_sig = ext_sig_cand.is_file() and ext_sig_cand.stat().st_size > 0
    if not has_ext_sig and not args.allow_pending_signalp:
        require_file(ext_sig_cand, "external247_signalp6_predictions.tsv.gz")

    if has_ext_sig:
        require_file(ext_man_cand, "external_signalp_manifest.tsv")
        thermo_tar_p = require_file(P_IN / "thermo_proteomes.tar.gz", "thermo_proteomes.tar.gz")
        print("[Stage 2b] Loading EQ 237 SignalP 6.0h predictions...", flush=True)
        _, sp_eq = load_all_sp_tm(n_workers=16)

        print("[Stage 2b] Loading External 247 SignalP 6.0h predictions and manifest from IBEX...", flush=True)
        ext_sig_df = pd.read_csv(ext_sig_cand, sep="\t", low_memory=False)
        ext_man_df = pd.read_csv(ext_man_cand, sep="\t", low_memory=False)
        ext_genome_ids = set(
            df_484.loc[df_484["cohort"] == "External_Comparators_Deduplicated_247", "genome"].astype(str).unique()
        )
        print("[Stage 2b] Verifying per-genome .faa SHA-256 digests against direct-search source thermo_proteomes.tar.gz...", flush=True)
        expected_ext_sha_map = compute_thermo_proteomes_tar_sha256_map(thermo_tar_p, target_genome_ids=ext_genome_ids)

        df_484 = validate_full_stage2b_signalp_and_uniref50_provenance(
            df_484,
            ext_sig_df,
            ext_man_df,
            sp_eq,
            expected_ext_genome_faa_sha256=expected_ext_sha_map,
        )
        print(
            f"[Stage 2b] Verified exact (genome, gene) key sets, zero duplicates, and .faa SHA-256 provenance across all {len(df_484)} genes.",
            flush=True,
        )
    else:
        print(
            "[Stage 2b] external247_signalp6_predictions.tsv.gz is pending IBEX job 55372579; "
            "computing canonical-vs-legacy direct UniRef50 darkness models now.",
            flush=True,
        )

    darkness_df, signalp_df = evaluate_stage2b_uniref50_and_signalp_models(
        df_484,
        eq_reconciled,
        tax_summary,
        pairs_df,
        has_signalp_joined=has_ext_sig,
    )
    darkness_df.to_csv(
        P_RES / "task_b_stage2b_uniref50_eq_vs_external_canonical_and_legacy_taxonomy.tsv",
        sep="\t",
        index=False,
    )
    if has_ext_sig:
        signalp_df.to_csv(
            P_RES / "task_b_external_vs_eq_signalp6_by_direct_uniref50_darkness.tsv",
            sep="\t",
            index=False,
        )
    print("[Stage 2b] Completed successfully!", flush=True)


if __name__ == "__main__":
    main()
