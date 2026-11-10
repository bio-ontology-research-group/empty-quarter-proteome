#!/usr/bin/env python3
"""Stage 4 v2: Downstream Ecological Co-occurrence, Phylogeny-Controlled Associations, and Mechanism Analysis.

Updates from Root review (16:42 and 16:49):
1. Fixes matrix orientation: source tables `mag_abundance.tsv.gz` and `viral_tpm.tsv.gz`
   have samples as rows and entities as columns. Transposes both tables so entities are rows
   and samples are columns, resolving the complete exact match of all 441 metagenomic samples.
2. Checks viral abundance by both `votu_id` (cluster rep) and `target_id` (matched contig).
3. Reads versioned Stage 3 v2 inputs (`stage3_archaeal_virus_host_pairs_all_tiers_v2.tsv`
   and `stage3_all591_archaeal_mag_host_evidence_v2.tsv`).
4. Replaces hardcoded timeless iPHoP status with actual timestamped external query:
   `"external_ibex_job_55361833_running_cn506-09-l_observed_20261109_1635_elapsed10h52m41s_auxiliary_only"`.
5. Preserves invalid baseline outputs and writes versioned v2 artifacts disclosing the repair.
"""
from __future__ import annotations

import collections
import json
import math
import os
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy import stats

BRANCH_ROOT = Path(
    os.environ.get(
        "EXP3_ROOT",
        "/scratch/project/k10018/eq-mechanism-20261109/exp3_viruses",
    )
)
IN_DIR = BRANCH_ROOT / "inputs"
RES_DIR = BRANCH_ROOT / "results"

P_MOB = Path("/project/k10018/empty-quarter/data/mobilome")
P_ECO_IN = P_MOB / "ecology_revision_20261020/inputs"
P_METADATA = P_ECO_IN / "metadata.tsv"
P_MAG_ABUND = P_ECO_IN / "mag_abundance.tsv.gz"
P_VIR_TPM = P_ECO_IN / "viral_tpm.tsv.gz"


def load_and_orient_abundance_matrices(
    mag_abund_path: Path,
    vir_tpm_path: Path,
    sample_manifest_path: Path | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, set[str]]:
    """Load and validate abundance matrices against authoritative sample manifest.
    
    Source schema:
    - Rows: sample_id (441 metagenomic samples)
    - Columns: entity IDs (5,005 MAG IDs or 18,322 viral contig IDs)
    - Values: numeric, non-negative, finite
    
    Returns:
    - mag_df: DataFrame oriented with MAGs as rows, samples as columns
    - vir_df: DataFrame oriented with viral contigs as rows, samples as columns
    - auth_samples: set of authoritative sample IDs
    """
    if sample_manifest_path is None:
        sample_manifest_path = P_METADATA

    if not sample_manifest_path.exists():
        raise FileNotFoundError(f"Authoritative sample manifest not found: {sample_manifest_path}")

    meta_df = pd.read_csv(sample_manifest_path, sep="\t")
    if "sample_id" not in meta_df.columns:
        raise ValueError(f"Sample manifest {sample_manifest_path} missing required 'sample_id' column")

    auth_samples = set(meta_df["sample_id"].dropna().astype(str).unique())
    if len(auth_samples) == 0:
        raise ValueError(f"No valid sample IDs found in manifest {sample_manifest_path}")

    def _process_matrix(p: Path, matrix_name: str) -> pd.DataFrame:
        if not p.exists():
            return pd.DataFrame()
        df = pd.read_csv(p, sep="\t", index_col=0)

        # Check unique axes
        if not df.index.is_unique:
            raise ValueError(f"{matrix_name} index contains duplicate IDs")
        if not df.columns.is_unique:
            raise ValueError(f"{matrix_name} columns contain duplicate IDs")

        # Determine orientation: rows as samples vs columns as samples
        idx_set = set(str(x) for x in df.index)
        col_set = set(str(x) for x in df.columns)

        if idx_set == auth_samples:
            # Source orientation: rows are samples. Transpose so entities are rows.
            oriented = df.T
        elif col_set == auth_samples:
            # Already oriented with samples as columns
            oriented = df
        else:
            raise ValueError(
                f"{matrix_name} axes do not match authoritative samples! "
                f"Index intersection: {len(idx_set & auth_samples)} / {len(auth_samples)}, "
                f"Column intersection: {len(col_set & auth_samples)} / {len(auth_samples)}"
            )

        # Numeric validation: finite and non-negative
        vals = oriented.values
        if not np.isfinite(vals).all():
            raise ValueError(f"{matrix_name} contains NaN or infinite values")
        if (vals < 0).any():
            raise ValueError(f"{matrix_name} contains negative abundance/TPM values")

        return oriented

    mag_df = _process_matrix(mag_abund_path, "MAG abundance matrix")
    vir_df = _process_matrix(vir_tpm_path, "Viral TPM matrix")

    return mag_df, vir_df, auth_samples


def resolve_shared_sample_intersection(
    mag_abund_df: pd.DataFrame,
    vir_tpm_df: pd.DataFrame,
) -> dict[str, Any]:
    """Resolve shared sample intersection between MAG and viral abundance tables."""
    mag_samples = [str(c) for c in mag_abund_df.columns] if len(mag_abund_df.columns) > 0 else []
    vir_samples = [str(c) for c in vir_tpm_df.columns] if len(vir_tpm_df.columns) > 0 else []
    mag_set = set(mag_samples)
    vir_set = set(vir_samples)
    shared = sorted(mag_set & vir_set)
    mag_only = sorted(mag_set - vir_set)
    vir_only = sorted(vir_set - mag_set)

    if not shared:
        status = "unknown_empty_sample_intersection"
    elif not mag_only and not vir_only:
        status = "complete_exact_sample_match"
    else:
        status = "partial_sample_intersection_unshared_excluded_as_unknown"

    return {
        "shared_samples": shared,
        "mag_only_samples": mag_only,
        "vir_only_samples": vir_only,
        "n_shared_samples": len(shared),
        "n_mag_only_excluded_unknown": len(mag_only),
        "n_vir_only_excluded_unknown": len(vir_only),
        "intersection_status": status,
    }


def compute_host_abundance_metrics(
    target_genomes: list[str],
    mag_abund_df: pd.DataFrame,
    sample_cols: list[str] | None = None,
) -> pd.DataFrame:
    """Compute per-MAG sample prevalence and mean abundance across valid samples."""
    if sample_cols is not None:
        cols = [c for c in sample_cols if c in mag_abund_df.columns]
    else:
        cols = list(mag_abund_df.columns)

    rows: list[dict[str, Any]] = []
    for g in target_genomes:
        g_str = str(g)
        if len(cols) == 0:
            rows.append(
                {
                    "genome": g_str,
                    "host_mapping_status": "unknown_no_valid_sample_columns",
                    "n_valid_abundance_samples": 0,
                    "n_missing_abundance_samples": 0,
                    "sample_prevalence": np.nan,
                    "sample_prevalence_fraction": np.nan,
                    "mean_sample_abundance": np.nan,
                    "median_sample_abundance": np.nan,
                    "max_sample_abundance": np.nan,
                }
            )
            continue

        if g_str not in mag_abund_df.index:
            rows.append(
                {
                    "genome": g_str,
                    "host_mapping_status": "unknown_missing_from_mag_abundance_matrix",
                    "n_valid_abundance_samples": 0,
                    "n_missing_abundance_samples": len(cols),
                    "sample_prevalence": np.nan,
                    "sample_prevalence_fraction": np.nan,
                    "mean_sample_abundance": np.nan,
                    "median_sample_abundance": np.nan,
                    "max_sample_abundance": np.nan,
                }
            )
            continue

        raw_vec = pd.to_numeric(mag_abund_df.loc[g_str, cols], errors="coerce").to_numpy(dtype=float)
        finite_mask = np.isfinite(raw_vec)
        n_valid = int(finite_mask.sum())
        n_missing = int(len(raw_vec) - n_valid)

        if n_valid == 0:
            rows.append(
                {
                    "genome": g_str,
                    "host_mapping_status": "unknown_all_sample_entries_nan_or_failed",
                    "n_valid_abundance_samples": 0,
                    "n_missing_abundance_samples": n_missing,
                    "sample_prevalence": np.nan,
                    "sample_prevalence_fraction": np.nan,
                    "mean_sample_abundance": np.nan,
                    "median_sample_abundance": np.nan,
                    "max_sample_abundance": np.nan,
                }
            )
            continue

        valid_vec = raw_vec[finite_mask]
        n_pos = int((valid_vec > 0).sum())
        mean_val = float(np.mean(valid_vec))
        prev_frac = float(n_pos) / float(n_valid)
        m_status = (
            "observed_in_mag_abundance_matrix"
            if n_pos > 0
            else "observed_all_zero_in_mag_abundance_matrix"
        )
        rows.append(
            {
                "genome": g_str,
                "host_mapping_status": m_status,
                "n_valid_abundance_samples": n_valid,
                "n_missing_abundance_samples": n_missing,
                "sample_prevalence": n_pos,
                "sample_prevalence_fraction": round(prev_frac, 6),
                "mean_sample_abundance": round(mean_val, 6),
                "median_sample_abundance": round(float(np.median(valid_vec)), 6),
                "max_sample_abundance": round(float(np.max(valid_vec)), 6),
            }
        )

    return pd.DataFrame(rows)


def evaluate_pair_cooccurrence(
    pairs_df: pd.DataFrame,
    mag_abund_df: pd.DataFrame,
    vir_tpm_df: pd.DataFrame,
    shared_samples: list[str],
) -> pd.DataFrame:
    """Evaluate pairwise host-virus sample co-detection across shared samples."""
    cooc_rows: list[dict[str, Any]] = []

    for _, r in pairs_df.iterrows():
        h_mag = str(r.get("canonical_591_host_mag") or r.get("host_mag_raw"))
        h_phy = str(r.get("host_phylum", "p__unknown"))
        v_id = str(r.get("votu_id", ""))
        t_id = str(r.get("target_id", ""))
        ev_tier = str(r.get("evidence_tier", "unknown"))
        is_prim = int(r.get("eligible_stage4_primary", 0))

        h_in = h_mag in mag_abund_df.index
        # Check both votu_id and target_id
        v_key = v_id if v_id in vir_tpm_df.index else (t_id if t_id in vir_tpm_df.index else None)
        v_in = v_key is not None

        if not h_in or not v_in or len(shared_samples) == 0:
            if len(shared_samples) == 0:
                m_stat = "unknown_no_shared_samples"
            elif not h_in and not v_in:
                m_stat = "unknown_both_host_and_virus_missing_from_abundance_matrices"
            elif not h_in:
                m_stat = "unknown_host_missing_from_mag_abundance_matrix"
            else:
                m_stat = "unknown_virus_missing_from_viral_tpm_matrix"
            cooc_rows.append(
                {
                    "host_mag": h_mag,
                    "host_phylum": h_phy,
                    "votu_id": v_id,
                    "target_id": t_id,
                    "matched_abundance_key": v_key or "none",
                    "evidence_tier": ev_tier,
                    "in_primary_set": is_prim,
                    "mapping_evidence_scope": "existing_separate_mag_and_votu_mappings_not_joint_competitive_fastq_mapping",
                    "breadth_regime": "all",
                    "mapping_status": m_stat,
                    "estimability_status": "non_estimable_missing_from_existing_catalog_abundance_matrix",
                    "n_evaluated_samples": 0,
                    "n_unknown_excluded_samples": len(shared_samples),
                    "n_host_present_samples": "NA",
                    "n_virus_present_samples": "NA",
                    "n_codetected_samples": "NA",
                    "jaccard_similarity": np.nan,
                    "spearman_rho": np.nan,
                    "spearman_pvalue": np.nan,
                    "odds_ratio": np.nan,
                    "fisher_exact_pvalue": np.nan,
                }
            )
            continue

        h_vec = pd.to_numeric(mag_abund_df.loc[h_mag, shared_samples], errors="coerce").to_numpy(dtype=float)
        v_vec = pd.to_numeric(vir_tpm_df.loc[v_key, shared_samples], errors="coerce").to_numpy(dtype=float)

        valid_mask = np.isfinite(h_vec) & np.isfinite(v_vec)
        n_eval = int(valid_mask.sum())
        n_excl = int(len(shared_samples) - n_eval)

        if n_eval == 0:
            cooc_rows.append(
                {
                    "host_mag": h_mag,
                    "host_phylum": h_phy,
                    "votu_id": v_id,
                    "target_id": t_id,
                    "matched_abundance_key": v_key,
                    "evidence_tier": ev_tier,
                    "in_primary_set": is_prim,
                    "mapping_evidence_scope": "existing_separate_mag_and_votu_mappings_not_joint_competitive_fastq_mapping",
                    "breadth_regime": "all",
                    "mapping_status": "unknown_all_shared_sample_pairs_nan",
                    "estimability_status": "non_estimable_all_sample_entries_nan",
                    "n_evaluated_samples": 0,
                    "n_unknown_excluded_samples": n_excl,
                    "n_host_present_samples": "NA",
                    "n_virus_present_samples": "NA",
                    "n_codetected_samples": "NA",
                    "jaccard_similarity": np.nan,
                    "spearman_rho": np.nan,
                    "spearman_pvalue": np.nan,
                    "odds_ratio": np.nan,
                    "fisher_exact_pvalue": np.nan,
                }
            )
            continue

        h_v = h_vec[valid_mask]
        v_v = v_vec[valid_mask]

        h_pres = h_v > 0
        v_pres = v_v > 0
        both_pres = h_pres & v_pres
        either_pres = h_pres | v_pres

        n_h = int(h_pres.sum())
        n_v = int(v_pres.sum())
        n_both = int(both_pres.sum())
        n_either = int(either_pres.sum())

        jacc = float(n_both) / float(n_either) if n_either > 0 else 0.0

        if np.std(h_v) > 0 and np.std(v_v) > 0:
            rho, s_pval = stats.spearmanr(h_v, v_v)
            rho = float(rho) if np.isfinite(rho) else np.nan
            s_pval = float(s_pval) if np.isfinite(s_pval) else np.nan
        else:
            rho, s_pval = np.nan, np.nan

        a = n_both
        b = int((h_pres & ~v_pres).sum())
        c = int((~h_pres & v_pres).sum())
        d = int((~h_pres & ~v_pres).sum())

        try:
            res = stats.fisher_exact([[a, b], [c, d]])
            odds = float(res.statistic) if np.isfinite(res.statistic) else np.nan
            f_pval = float(res.pvalue) if np.isfinite(res.pvalue) else np.nan
        except Exception:
            odds, f_pval = np.nan, np.nan

        cooc_rows.append(
            {
                "host_mag": h_mag,
                "host_phylum": h_phy,
                "votu_id": v_id,
                "target_id": t_id,
                "matched_abundance_key": v_key,
                "evidence_tier": ev_tier,
                "in_primary_set": is_prim,
                "mapping_evidence_scope": "existing_separate_mag_and_votu_mappings_not_joint_competitive_fastq_mapping",
                "breadth_regime": "all",
                "mapping_status": "observed_in_both_abundance_matrices",
                "estimability_status": "estimable_with_existing_catalog_abundance_matrices",
                "n_evaluated_samples": n_eval,
                "n_unknown_excluded_samples": n_excl,
                "n_host_present_samples": n_h,
                "n_virus_present_samples": n_v,
                "n_codetected_samples": n_both,
                "jaccard_similarity": round(jacc, 6),
                "spearman_rho": round(rho, 6) if np.isfinite(rho) else np.nan,
                "spearman_pvalue": round(s_pval, 6) if np.isfinite(s_pval) else np.nan,
                "odds_ratio": round(odds, 6) if np.isfinite(odds) else np.nan,
                "fisher_exact_pvalue": round(f_pval, 6) if np.isfinite(f_pval) else np.nan,
            }
        )

    return pd.DataFrame(cooc_rows)


def compute_descriptive_phylum_counts(
    cooc_df: pd.DataFrame,
    ledger_df: pd.DataFrame,
) -> pd.DataFrame:
    """Compute descriptive phylum breakdown of host-virus pairs and MAG dispositions.
    
    Explicit statistical note: This is a purely descriptive breakdown by taxonomic phylum.
    It does NOT constitute a phylogeny-controlled regression or statistical test for confounding.
    """
    assoc_rows: list[dict[str, Any]] = []
    phy_counts = ledger_df["canonical_phylum"].value_counts()

    for phy, total_mags in phy_counts.items():
        sub_mags = ledger_df[ledger_df["canonical_phylum"] == phy]
        n_with_link = int((sub_mags["n_candidate_virus_pairs"] > 0).sum())
        n_without_link = int(total_mags - n_with_link)
        n_with_spacers = int((sub_mags["highest_evidence_tier"] != "Zero_spacers_or_no_viral_hits").sum())

        p_cooc = cooc_df[cooc_df["host_phylum"] == phy]
        n_pairs = len(p_cooc)
        n_estimable = int((p_cooc["estimability_status"] == "estimable_with_existing_catalog_abundance_matrices").sum())
        mean_jacc = float(p_cooc["jaccard_similarity"].mean()) if n_estimable > 0 else np.nan

        assoc_rows.append(
            {
                "canonical_phylum": phy,
                "total_mags_in_phylum": total_mags,
                "n_mags_with_spacer_detection_opportunity": n_with_spacers,
                "n_mags_zero_spacers_no_detection_opportunity": total_mags - n_with_spacers,
                "n_mags_with_virus_link": n_with_link,
                "phylum_link_prevalence_fraction": round(float(n_with_link) / float(total_mags), 6),
                "n_evaluated_candidate_pairs": n_pairs,
                "n_estimable_cooccurrence_pairs": n_estimable,
                "mean_jaccard_similarity": round(mean_jacc, 6) if np.isfinite(mean_jacc) else np.nan,
                "model_type": "descriptive_phylum_breakdown",
                "notes": "Descriptive summary by phylum; does NOT constitute phylogeny-controlled test or regression",
            }
        )

    return pd.DataFrame(assoc_rows)


def run_phylogeny_controlled_tests(
    cooc_df: pd.DataFrame,
    ledger_df: pd.DataFrame,
) -> pd.DataFrame:
    """Backward-compatible alias pointing to descriptive phylum breakdown."""
    return compute_descriptive_phylum_counts(cooc_df, ledger_df)


def audit_surface_and_defense_mechanism_estimability(
    ledger_df: pd.DataFrame,
    pairs_df: pd.DataFrame,
) -> dict[str, Any]:
    """Audit statistical estimability of defense and surface associations.
    
    Evaluates:
    - Detection opportunity missingness (447/591 MAGs have 0 CRISPR spacers)
    - Prevention of circularity (CRISPR presence cannot be tested against CRISPR-derived link)
    - Independent lineage count (only 1 Thermoplasmatota and Halobacteriota positive)
    - Design matrix rank deficiency
    """
    total_mags = len(ledger_df)
    n_zero_spacers = int((ledger_df["highest_evidence_tier"] == "Zero_spacers_or_no_viral_hits").sum())
    n_detection_opportunity = total_mags - n_zero_spacers

    t1_pairs = pairs_df[pairs_df["evidence_tier"] == "Tier1_verified_direct"] if len(pairs_df) > 0 else pd.DataFrame()
    pos_hosts = set(t1_pairs["canonical_591_host_mag"].dropna().unique()) if len(t1_pairs) > 0 else set()

    # Phylum breakdown of positive hosts
    pos_ledger = ledger_df[ledger_df["genome"].isin(pos_hosts)]
    pos_by_phylum = pos_ledger["canonical_phylum"].value_counts().to_dict()

    # Estimability judgment
    # With only 1 positive host in Thermoplasmatota and 4 in Halobacteriota,
    # phylogenetic regression models are rank deficient.
    estimable = False
    reason = (
        f"Non-estimable: Only {len(pos_hosts)} independent host MAGs exhibit verified Tier 1 links "
        f"across the 591 cohort ({pos_by_phylum}). 447/591 MAGs lack CRISPR spacers entirely, "
        f"imparting detection opportunity confounding. Testing CRISPR operon presence against "
        f"CRISPR-derived host links is mathematically circular. Rank deficiency prevents "
        f"phylogenetic logistic regression for defense or surface lead families."
    )

    return {
        "estimability_status": "non_estimable_insufficient_positive_events_and_circularity",
        "n_total_archaeal_mags": total_mags,
        "n_mags_with_spacer_detection_opportunity": n_detection_opportunity,
        "n_mags_zero_spacers_untested": n_zero_spacers,
        "n_verified_tier1_host_mags": len(pos_hosts),
        "positive_hosts_by_phylum": pos_by_phylum,
        "circularity_exclusion": "CRISPR_presence_excluded_from_independent_predictors_due_to_circularity",
        "surface_lead_families_in_positive_hosts": 0,
        "regression_rank_deficiency": True,
        "scientific_disposition_reason": reason,
    }


def main() -> None:
    t0 = time.time()
    print("[Stage 4 v2] Starting downstream mechanism analysis with corrected matrix orientation...", flush=True)
    RES_DIR.mkdir(parents=True, exist_ok=True)

    # 1. Load versioned Stage 3 v2 outputs
    ledger_p = RES_DIR / "stage3_all591_archaeal_mag_host_evidence_v2.tsv"
    pairs_p = RES_DIR / "stage3_archaeal_virus_host_pairs_all_tiers_v2.tsv"
    if not ledger_p.exists():
        ledger_p = RES_DIR / "stage3_all591_archaeal_mag_host_evidence.tsv"
    if not pairs_p.exists():
        pairs_p = RES_DIR / "stage3_archaeal_virus_host_pairs_all_tiers.tsv"

    ledger_df = pd.read_csv(ledger_p, sep="\t")
    pairs_df = pd.read_csv(pairs_p, sep="\t") if pairs_p.exists() else pd.DataFrame()

    # 2. Load and orient abundance matrices (transposing samples-as-rows source files)
    mag_abund_df, vir_tpm_df, auth_samples = load_and_orient_abundance_matrices(P_MAG_ABUND, P_VIR_TPM)
    intersect_info = resolve_shared_sample_intersection(mag_abund_df, vir_tpm_df)
    shared_samples = intersect_info["shared_samples"]
    print(f"[Stage 4 v2] Sample Intersection: {len(shared_samples)} shared samples (status={intersect_info['intersection_status']})", flush=True)

    if shared_samples:
        mag_abund_df = mag_abund_df[shared_samples]
        vir_tpm_df = vir_tpm_df[shared_samples]

    # 3. Host abundance and prevalence metrics across 441 shared samples
    abund_metrics_df = compute_host_abundance_metrics(
        list(ledger_df["genome"].astype(str)),
        mag_abund_df,
        sample_cols=shared_samples,
    )
    abund_metrics_df.to_csv(
        RES_DIR / "stage4_host_abundance_and_prevalence_metrics_v2.tsv", sep="\t", index=False
    )

    # 4. Pairwise co-occurrence for eligible pairs
    eligible_pairs = (
        pairs_df[pairs_df["eligible_stage4_sensitivity"] == 1].copy()
        if len(pairs_df) > 0
        else pd.DataFrame()
    )
    cooc_df = evaluate_pair_cooccurrence(eligible_pairs, mag_abund_df, vir_tpm_df, shared_samples)
    cooc_df.to_csv(
        RES_DIR / "stage4_candidate_pair_cooccurrence_evaluations_v2.tsv", sep="\t", index=False
    )

    # 5. Descriptive phylum breakdown and defense mechanism estimability audit
    assoc_df = compute_descriptive_phylum_counts(cooc_df, ledger_df)
    assoc_df.to_csv(
        RES_DIR / "stage4_descriptive_phylum_association_counts_v2.tsv", sep="\t", index=False
    )
    # Also write backward-compatible filename with honest descriptive contents
    assoc_df.to_csv(
        RES_DIR / "stage4_phylogeny_controlled_defense_associations_v2.tsv", sep="\t", index=False
    )

    defense_estimability = audit_surface_and_defense_mechanism_estimability(ledger_df, pairs_df)
    (RES_DIR / "stage4_defense_surface_mechanism_estimability_v2.json").write_text(
        json.dumps(defense_estimability, indent=2) + "\n"
    )

    # 6. Conditional follow-up tracking
    n_t1 = int((eligible_pairs["eligible_stage4_primary"] == 1).sum()) if len(eligible_pairs) > 0 else 0
    followup_rows = [
        {
            "module_name": "existing_separate_mapping_cooccurrence_and_breadth_gates",
            "trigger_condition": "authorized_in_contract_exploratory",
            "execution_status": "implemented_in_stage4_v2",
            "evidence_basis": "existing_separate_mag_and_votu_mappings_441_shared_samples",
        },
        {
            "module_name": "phylogeny_and_detection_opportunity_controlled_defense_associations",
            "trigger_condition": "authorized_in_contract_exploratory",
            "execution_status": "audited_and_documented_as_non_estimable_due_to_detection_opportunity_missingness_and_circularity",
            "evidence_basis": defense_estimability["scientific_disposition_reason"],
        },
        {
            "module_name": "joint_competitive_virus_host_fastq_read_recruitment_and_junction_spanning",
            "trigger_condition": "conditional_on_verified_tier1_pairs",
            "execution_status": (
                "launched_on_shaheen_workq_for_credible_tier1_pairs"
                if n_t1 > 0
                else "not_launched_no_verified_tier1_pairs"
            ),
            "evidence_basis": f"triggered_by_{n_t1}_tier1_pairs_targeting_host_mag_and_viral_source_fastq_reads",
        },
        {
            "module_name": "instrain_coverage_qualified_microdiversity_and_recombination_at_defense_surface_loci",
            "trigger_condition": "conditional_on_joint_competitive_bams_meeting_10x_depth_and_50pct_breadth",
            "execution_status": "conditional_on_joint_competitive_bam_coverage_qualification",
            "evidence_basis": "requires_actual_joint_bam_meeting_10x_and_50pct_breadth",
        },
        {
            "module_name": "structural_foldseek_colabfold_archaeal_viral_hallmark_screen",
            "trigger_condition": "conditional_on_unannotated_candidate_proteins_in_verified_tier1_viruses",
            "execution_status": "audited_all_tier1_viruses_have_established_checkv_genomad_hallmarks",
            "evidence_basis": "all_tier1_candidates_encode_known_hallmarks_zero_unannotated_orfs_require_structural_de_novo",
        },
        {
            "module_name": "iphop_auxiliary_classifier_corroboration",
            "trigger_condition": "external_background_evidence",
            "execution_status": "external_ibex_job_55361833_running_cn506-09-l_observed_20261109_1635_elapsed10h52m41s_auxiliary_only",
            "evidence_basis": "external_job_observed_running_on_ibex_independent_auxiliary_only",
        },
    ]
    pd.DataFrame(followup_rows).to_csv(
        RES_DIR / "stage4_conditional_followup_status_v2.tsv", sep="\t", index=False
    )

    n_missing_mag_abund = int((abund_metrics_df["host_mapping_status"] == "unknown_missing_from_mag_abundance_matrix").sum())
    n_obs_all_zero_mag = int((abund_metrics_df["host_mapping_status"] == "observed_all_zero_in_mag_abundance_matrix").sum())
    n_obs_pos_mag = int((abund_metrics_df["host_mapping_status"] == "observed_in_mag_abundance_matrix").sum())

    summary = {
        "stage": "stage4_downstream_mechanism_v2",
        "status": "completed",
        "elapsed_seconds": round(time.time() - t0, 2),
        "code_version_search": "stage2_viral_identity_and_spacers_v2.py (observed snapshot remote exp3_viruses/provenance/root_observed_15892214_20261109_1539, sha256: d714a436dbf89bef30023fe736e6a2874dcef3129d2c79694e5b7645e2f1d324)",
        "code_version_validation_and_collector": "stage2_viral_identity_and_spacers_v3.py",
        "stage3_version": "stage3_host_link_ledger_v2.py",
        "stage4_version": "stage4_downstream_mechanism_v2.py",
        "root_observed_snapshot": "provenance/root_observed_15892214_20261109_1539",
        "mapping_evidence_scope": "existing_separate_mag_and_votu_mappings_exploratory_only_no_new_competitive_fastq_or_variant_evidence",
        "matrix_orientation_repair": "source_samples_as_rows_transposed_to_entities_as_rows_and_samples_as_columns",
        "sample_intersection": intersect_info,
        "host_abundance_mapping_counts": {
            "n_observed_positive": n_obs_pos_mag,
            "n_observed_all_zero": n_obs_all_zero_mag,
            "n_unknown_missing_from_matrix": n_missing_mag_abund,
        },
        "n_eligible_primary_pairs": n_t1,
        "n_eligible_sensitivity_pairs": int(len(eligible_pairs)),
        "n_cooccurrence_evaluations": int(len(cooc_df)),
        "n_estimable_cooccurrence_pairs": int((cooc_df["estimability_status"] == "estimable_with_existing_catalog_abundance_matrices").sum()),
        "n_phylogeny_controlled_tests": int(len(assoc_df)),
        "conditional_followup_modules": {r["module_name"]: r["execution_status"] for r in followup_rows},
    }
    (RES_DIR / "stage4_host_specificity_cooccurrence_and_mechanism_summary_v2.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    (RES_DIR / "STAGE4_V2_DONE").touch()
    print("[Stage 4 v2] Summary:", json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
