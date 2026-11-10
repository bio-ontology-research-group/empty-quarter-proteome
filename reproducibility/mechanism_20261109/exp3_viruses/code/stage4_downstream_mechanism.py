#!/usr/bin/env python3
"""Stage 4: Downstream Host-Linked Mechanism Analyses Using Existing Pre-Computed
Mappings (with Breadth Gates & Matched Host Controls), Phylogeny/Detection-Opportunity
Controlled Defense/Surface Associations, and Explicit Conditional Follow-Up Ledger.

Implements all requirements from `EXP3-CONTRACT.txt`, `EXP3-RUNTIME-HANDOFF.txt`,
and `EXP3-FINAL-PROVENANCE-REVIEW.txt`:
1. Consumes ONLY eligible verified links from Stage 3 (`Tier1_verified_direct` primary;
   `Tier1` + `Tier2` sensitivity), excluding all `Tier3_*`, `Unknown_unverified_viral_identity`,
   `Unresolved_unexcised_provirus_candidate`, `External_habitat_viral_control_match_retained`,
   `Cellular_or_plasmid_excluded`, `Decoy_control_rejected`, and `Technical_self_crispr_array_masked`.
2. Missing MAG or viral abundances, unshared sample columns, and failed/missing alignment or
   caller rows are strictly UNKNOWN (`NaN` / `"NA"` + explicit `mapping_status`), NEVER coerced
   to `0` / true absence, and are excluded from prevalence, co-occurrence, control matching,
   and regression denominators. Observed all-zero rows (`observed_all_zero_in_mag_abundance_matrix`)
   are tracked distinctly from missing MAGs (`unknown_missing_from_mag_abundance_matrix`).
3. Explicitly labels co-occurrence / breadth analyses as exploratory analyses derived from
   pre-existing separate catalog mappings (`mapping_evidence_scope = "existing_separate_mag_and_votu_mappings_not_joint_competitive_fastq_mapping"`),
   never claiming competitive virus-to-host abundance ratios or variant evidence from separate mappings.
4. Explicitly records the pending / non-estimable status and exact reason for conditional
   follow-up branches (`joint_competitive_virus_host_fastq_mapping`,
   `instrain_recombination_and_variant_analysis`, `structural_foldseek_colabfold_hallmark_screen`,
   and `iphop_auxiliary_corroboration`) in `results/stage4_conditional_followup_status.tsv`.
"""
from __future__ import annotations

import collections
import gzip
import json
import os
import tarfile
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
P_MAG_ABUND = P_ECO_IN / "mag_abundance.tsv.gz"
P_VIR_TPM = P_ECO_IN / "viral_tpm.tsv.gz"
P_E1_TAR = P_MOB / "results_20261104/curated/e1_metrics.tar"
P_ARCH_IN = Path("/scratch/project/k10018/eq-archaeal-function-20261106/inputs")


def resolve_shared_sample_intersection(
    mag_abund_df: pd.DataFrame,
    vir_tpm_df: pd.DataFrame,
) -> dict[str, Any]:
    """Resolve shared sample intersection between MAG and viral abundance tables.
    Unshared samples are never imputed as zero; they are excluded as unknown.
    """
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
    genomes: list[str],
    mag_abund_df: pd.DataFrame,
    sample_cols: list[str] | None = None,
) -> pd.DataFrame:
    """Compute per-MAG sample prevalence and mean abundance while keeping missing MAGs
    and NaN sample entries strictly as `np.nan` + explicit `host_mapping_status`.

    Distinguishes:
    - `unknown_missing_from_mag_abundance_matrix`: genome absent from `mag_abund_df.index` -> `NaN`
    - `unknown_all_sample_entries_nan_or_failed`: genome row present but all evaluated sample entries are `NaN` -> `NaN`
    - `observed_all_zero_in_mag_abundance_matrix`: genome row present with finite `0.0` across all valid samples -> `0.0`
    - `observed_in_mag_abundance_matrix`: genome row present with >= 1 positive sample -> finite > 0
    """
    if sample_cols is not None:
        cols = [c for c in sample_cols if c in mag_abund_df.columns]
    else:
        cols = list(mag_abund_df.columns)

    rows: list[dict[str, Any]] = []
    for g in genomes:
        g_str = str(g)
        if mag_abund_df.empty or len(cols) == 0:
            rows.append(
                {
                    "genome": g_str,
                    "host_mapping_status": "unknown_no_shared_abundance_samples",
                    "n_valid_abundance_samples": 0,
                    "n_missing_abundance_samples": len(sample_cols) if sample_cols is not None else 0,
                    "sample_prevalence": np.nan,
                    "sample_prevalence_fraction": np.nan,
                    "mean_sample_abundance": np.nan,
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
                }
            )
            continue

        raw_series = pd.to_numeric(mag_abund_df.loc[g_str, cols], errors="coerce")
        raw_vec = raw_series.to_numpy(dtype=float)
        finite_mask = np.isfinite(raw_vec)
        n_valid = int(finite_mask.sum())
        n_missing = int((~finite_mask).sum())

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
                "sample_prevalence": float(n_pos),
                "sample_prevalence_fraction": prev_frac,
                "mean_sample_abundance": mean_val,
            }
        )

    return pd.DataFrame(rows)


def load_e1_breadth_masks(
    target_votus: set[str],
    e1_tar_path: Path = P_E1_TAR,
) -> tuple[dict[str, set[str]], dict[str, set[str]], set[str], set[str]]:
    """Load per-sample breadth >= 0.50 (e1 rule) and breadth >= 0.75 (strict) masks
    from `e1_metrics.tar` for the target vOTUs, along with `observed_samples` and
    `failed_samples` so missing/failed sample tables are treated as unknown rather than 0.
    """
    e1_pass: dict[str, set[str]] = collections.defaultdict(set)
    b75_pass: dict[str, set[str]] = collections.defaultdict(set)
    observed_samples: set[str] = set()
    failed_samples: set[str] = set()
    if not e1_tar_path.exists():
        return e1_pass, b75_pass, observed_samples, failed_samples

    with tarfile.open(e1_tar_path, "r:") as tf:
        for member in tf:
            if not member.isfile() or not member.name.endswith(".tsv.gz"):
                continue
            sample_id = Path(member.name).name.replace(".tsv.gz", "")
            fobj = tf.extractfile(member)
            if fobj is None:
                failed_samples.add(sample_id)
                continue
            try:
                with gzip.open(fobj, "rt") as gz:
                    hdr = gz.readline().rstrip("\n").split("\t")
                    if not hdr or hdr == [""]:
                        failed_samples.add(sample_id)
                        continue
                    col_idx = {c: i for i, c in enumerate(hdr)}
                    c_vid = col_idx.get("contig_id", col_idx.get("votu_id", 0))
                    c_br = col_idx.get("breadth", col_idx.get("covered_fraction", -1))
                    c_cov_bp = col_idx.get("covered_bases", -1)
                    if c_br < 0 and c_cov_bp < 0:
                        failed_samples.add(sample_id)
                        continue
                    for line in gz:
                        parts = line.rstrip("\n").split("\t")
                        if len(parts) <= c_vid:
                            continue
                        vid = parts[c_vid]
                        if target_votus and vid not in target_votus:
                            continue
                        br_str = parts[c_br] if (0 <= c_br < len(parts)) else ""
                        cov_str = parts[c_cov_bp] if (0 <= c_cov_bp < len(parts)) else ""
                        br = float(br_str) if br_str not in ("", "NA", "NaN", "nan") else np.nan
                        cov_bp = float(cov_str) if cov_str not in ("", "NA", "NaN", "nan") else np.nan
                        if (np.isfinite(br) and br >= 0.50) or (np.isfinite(cov_bp) and cov_bp >= 5000):
                            e1_pass[vid].add(sample_id)
                        if np.isfinite(br) and br >= 0.75:
                            b75_pass[vid].add(sample_id)
                    observed_samples.add(sample_id)
            except Exception:
                failed_samples.add(sample_id)
                observed_samples.discard(sample_id)

    return e1_pass, b75_pass, observed_samples, failed_samples


def evaluate_pair_cooccurrence(
    eligible_pairs: pd.DataFrame,
    ledger_df: pd.DataFrame,
    mag_abund_df: pd.DataFrame,
    vir_tpm_df: pd.DataFrame,
    shared_samples: list[str],
    e1_mask: dict[str, set[str]],
    b75_mask: dict[str, set[str]],
    e1_observed_samples: set[str] | None = None,
    e1_failed_samples: set[str] | None = None,
) -> pd.DataFrame:
    """Evaluate virus-host co-occurrence across shared samples while preserving unknown (`NA`)
    for missing MAGs, missing vOTUs, failed/missing alignment rows, and missing/failed sample metrics.
    """
    cooc_rows: list[dict[str, Any]] = []
    if len(eligible_pairs) == 0:
        return pd.DataFrame(cooc_rows)

    for _, pr in eligible_pairs.iterrows():
        h_mag = str(pr.get("canonical_591_host_mag", ""))
        v_id = str(pr.get("votu_id", ""))
        h_phy = str(pr.get("host_phylum", ""))
        ev_tier = str(pr.get("evidence_tier", ""))
        is_prim = int(pr.get("eligible_stage4_primary", 0) or 0)
        align_status = str(pr.get("alignment_status", pr.get("pair_disposition", "valid")))

        # Preserve unknown if the pair row itself represents a failed or missing alignment
        if "failed" in align_status.lower() or "missing" in align_status.lower() or "unknown" in ev_tier.lower():
            cooc_rows.append(
                {
                    "host_mag": h_mag,
                    "host_phylum": h_phy,
                    "votu_id": v_id,
                    "evidence_tier": ev_tier,
                    "in_primary_set": is_prim,
                    "mapping_evidence_scope": "existing_separate_mag_and_votu_mappings_not_joint_competitive_fastq_mapping",
                    "breadth_regime": "all",
                    "mapping_status": "unknown_failed_or_missing_alignment_row",
                    "estimability_status": "non_estimable_unknown_failed_or_missing_alignment_row",
                    "n_evaluated_samples": 0,
                    "n_unknown_excluded_samples": len(shared_samples),
                    "n_host_present_samples": "NA",
                    "n_virus_present_samples": "NA",
                    "n_codetected_samples": "NA",
                    "jaccard_observed": "NA",
                    "p_virus_given_host_observed": "NA",
                    "matched_control_mean_jaccard": "NA",
                    "delta_jaccard_vs_matched_unrelated_hosts": "NA",
                    "empirical_p_vs_matched_controls": "NA",
                    "competitive_vhr_status": "non_estimable_requires_joint_competitive_fastq_mapping",
                }
            )
            continue

        h_in = (not mag_abund_df.empty) and (h_mag in mag_abund_df.index)
        v_in = (not vir_tpm_df.empty) and (v_id in vir_tpm_df.index)
        if len(shared_samples) == 0 or not h_in or not v_in:
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
                    "jaccard_observed": "NA",
                    "p_virus_given_host_observed": "NA",
                    "matched_control_mean_jaccard": "NA",
                    "delta_jaccard_vs_matched_unrelated_hosts": "NA",
                    "empirical_p_vs_matched_controls": "NA",
                    "competitive_vhr_status": "non_estimable_requires_joint_competitive_fastq_mapping",
                }
            )
            continue

        h_raw = pd.to_numeric(mag_abund_df.loc[h_mag, shared_samples], errors="coerce").to_numpy(dtype=float)
        v_raw = pd.to_numeric(vir_tpm_df.loc[v_id, shared_samples], errors="coerce").to_numpy(dtype=float)
        base_valid_mask = np.isfinite(h_raw) & np.isfinite(v_raw)

        if int(base_valid_mask.sum()) == 0:
            cooc_rows.append(
                {
                    "host_mag": h_mag,
                    "host_phylum": h_phy,
                    "votu_id": v_id,
                    "evidence_tier": ev_tier,
                    "in_primary_set": is_prim,
                    "mapping_evidence_scope": "existing_separate_mag_and_votu_mappings_not_joint_competitive_fastq_mapping",
                    "breadth_regime": "all",
                    "mapping_status": "unknown_all_shared_sample_values_nan",
                    "estimability_status": "non_estimable_all_shared_sample_values_nan",
                    "n_evaluated_samples": 0,
                    "n_unknown_excluded_samples": len(shared_samples),
                    "n_host_present_samples": "NA",
                    "n_virus_present_samples": "NA",
                    "n_codetected_samples": "NA",
                    "jaccard_observed": "NA",
                    "p_virus_given_host_observed": "NA",
                    "matched_control_mean_jaccard": "NA",
                    "delta_jaccard_vs_matched_unrelated_hosts": "NA",
                    "empirical_p_vs_matched_controls": "NA",
                    "competitive_vhr_status": "non_estimable_requires_joint_competitive_fastq_mapping",
                }
            )
            continue

        # Candidate controls MUST be observed in `mag_abund_df` (`sample_prevalence.notna()`)
        # and have `Zero_spacers_detected` (never missing MAGs or failed/unscreened callers!)
        obs_ledger = ledger_df[
            ledger_df["sample_prevalence"].notna()
            & (ledger_df["genome"] != h_mag)
            & (ledger_df["genome"].isin(mag_abund_df.index))
        ].copy()
        cand_ctrls = obs_ledger[
            (obs_ledger["canonical_phylum"] == h_phy)
            & (obs_ledger["mag_host_disposition"] == "Zero_spacers_detected")
        ].copy()
        if len(cand_ctrls) < 5:
            cand_ctrls = obs_ledger[obs_ledger["mag_host_disposition"] == "Zero_spacers_detected"].copy()
        host_total_prev = int((h_raw[base_valid_mask] > 0).sum())
        cand_ctrls["prev_diff"] = (cand_ctrls["sample_prevalence"] - host_total_prev).abs()
        matched_ctrls = list(cand_ctrls.sort_values("prev_diff")["genome"].head(25))

        for regime in ("unrestricted_tpm_gt0", "e1_presence_rule", "strict_breadth75"):
            if regime == "unrestricted_tpm_gt0":
                regime_valid_mask = base_valid_mask.copy()
                v_pres_full = v_raw > 0
            else:
                # For breadth regimes, if e1_observed_samples is explicitly provided,
                # samples missing or failed in e1_metrics.tar are UNKNOWN (excluded, never 0!)
                if e1_observed_samples is not None:
                    if len(e1_observed_samples) == 0:
                        cooc_rows.append(
                            {
                                "host_mag": h_mag,
                                "host_phylum": h_phy,
                                "votu_id": v_id,
                                "evidence_tier": ev_tier,
                                "in_primary_set": is_prim,
                                "mapping_evidence_scope": "existing_separate_mag_and_votu_mappings_not_joint_competitive_fastq_mapping",
                                "breadth_regime": regime,
                                "mapping_status": "unknown_missing_e1_breadth_metrics",
                                "estimability_status": "non_estimable_missing_e1_breadth_metrics",
                                "n_evaluated_samples": 0,
                                "n_unknown_excluded_samples": len(shared_samples),
                                "n_host_present_samples": "NA",
                                "n_virus_present_samples": "NA",
                                "n_codetected_samples": "NA",
                                "jaccard_observed": "NA",
                                "p_virus_given_host_observed": "NA",
                                "matched_control_mean_jaccard": "NA",
                                "delta_jaccard_vs_matched_unrelated_hosts": "NA",
                                "empirical_p_vs_matched_controls": "NA",
                                "competitive_vhr_status": "non_estimable_requires_joint_competitive_fastq_mapping",
                            }
                        )
                        continue
                    e1_ok = np.array(
                        [
                            (s in e1_observed_samples)
                            and (e1_failed_samples is None or s not in e1_failed_samples)
                            for s in shared_samples
                        ],
                        dtype=bool,
                    )
                    regime_valid_mask = base_valid_mask & e1_ok
                else:
                    regime_valid_mask = base_valid_mask.copy()

                allowed = e1_mask.get(v_id, set()) if regime == "e1_presence_rule" else b75_mask.get(v_id, set())
                m_arr = np.array([s in allowed for s in shared_samples], dtype=bool)
                v_pres_full = (v_raw > 0) & m_arr

            n_eval = int(regime_valid_mask.sum())
            n_unk = int((~regime_valid_mask).sum())
            if n_eval == 0:
                cooc_rows.append(
                    {
                        "host_mag": h_mag,
                        "host_phylum": h_phy,
                        "votu_id": v_id,
                        "evidence_tier": ev_tier,
                        "in_primary_set": is_prim,
                        "mapping_evidence_scope": "existing_separate_mag_and_votu_mappings_not_joint_competitive_fastq_mapping",
                        "breadth_regime": regime,
                        "mapping_status": "unknown_zero_valid_samples_in_regime",
                        "estimability_status": "non_estimable_zero_valid_samples_in_regime",
                        "n_evaluated_samples": 0,
                        "n_unknown_excluded_samples": n_unk,
                        "n_host_present_samples": "NA",
                        "n_virus_present_samples": "NA",
                        "n_codetected_samples": "NA",
                        "jaccard_observed": "NA",
                        "p_virus_given_host_observed": "NA",
                        "matched_control_mean_jaccard": "NA",
                        "delta_jaccard_vs_matched_unrelated_hosts": "NA",
                        "empirical_p_vs_matched_controls": "NA",
                        "competitive_vhr_status": "non_estimable_requires_joint_competitive_fastq_mapping",
                    }
                )
                continue

            h_pres = (h_raw > 0)[regime_valid_mask]
            v_pres = v_pres_full[regime_valid_mask]

            n_h = int(h_pres.sum())
            n_v = int(v_pres.sum())
            n_co = int((h_pres & v_pres).sum())
            n_un = int((h_pres | v_pres).sum())

            if n_h == 0 or n_v == 0:
                cooc_rows.append(
                    {
                        "host_mag": h_mag,
                        "host_phylum": h_phy,
                        "votu_id": v_id,
                        "evidence_tier": ev_tier,
                        "in_primary_set": is_prim,
                        "mapping_evidence_scope": "existing_separate_mag_and_votu_mappings_not_joint_competitive_fastq_mapping",
                        "breadth_regime": regime,
                        "mapping_status": "observed_zero_detections_under_breadth_gate",
                        "estimability_status": "non_estimable_zero_detections_under_breadth_gate",
                        "n_evaluated_samples": n_eval,
                        "n_unknown_excluded_samples": n_unk,
                        "n_host_present_samples": n_h,
                        "n_virus_present_samples": n_v,
                        "n_codetected_samples": n_co,
                        "jaccard_observed": 0.0,
                        "p_virus_given_host_observed": 0.0,
                        "matched_control_mean_jaccard": "NA",
                        "delta_jaccard_vs_matched_unrelated_hosts": "NA",
                        "empirical_p_vs_matched_controls": "NA",
                        "competitive_vhr_status": "non_estimable_requires_joint_competitive_fastq_mapping",
                    }
                )
                continue

            j_obs = n_co / float(n_un) if n_un > 0 else 0.0
            p_vh = n_co / float(n_h) if n_h > 0 else 0.0

            ctrl_jaccards = []
            for cm in matched_ctrls:
                if cm in mag_abund_df.index:
                    c_raw = pd.to_numeric(mag_abund_df.loc[cm, shared_samples], errors="coerce").to_numpy(dtype=float)
                    c_valid = regime_valid_mask & np.isfinite(c_raw)
                    if int(c_valid.sum()) == 0:
                        continue
                    c_pres = (c_raw > 0)[c_valid]
                    v_pres_c = v_pres_full[c_valid]
                    c_un = int((c_pres | v_pres_c).sum())
                    c_co = int((c_pres & v_pres_c).sum())
                    ctrl_jaccards.append(c_co / float(c_un) if c_un > 0 else 0.0)

            mean_ctrl_j = float(np.mean(ctrl_jaccards)) if ctrl_jaccards else np.nan
            delta_j = (j_obs - mean_ctrl_j) if np.isfinite(mean_ctrl_j) else np.nan
            emp_p = (
                (1 + sum(1 for cj in ctrl_jaccards if cj >= j_obs)) / float(1 + len(ctrl_jaccards))
                if ctrl_jaccards
                else np.nan
            )

            cooc_rows.append(
                {
                    "host_mag": h_mag,
                    "host_phylum": h_phy,
                    "votu_id": v_id,
                    "evidence_tier": ev_tier,
                    "in_primary_set": is_prim,
                    "mapping_evidence_scope": "existing_separate_mag_and_votu_mappings_not_joint_competitive_fastq_mapping",
                    "breadth_regime": regime,
                    "mapping_status": "observed_in_both_host_and_viral_matrices",
                    "estimability_status": "estimable_from_existing_separate_mappings",
                    "n_evaluated_samples": n_eval,
                    "n_unknown_excluded_samples": n_unk,
                    "n_host_present_samples": n_h,
                    "n_virus_present_samples": n_v,
                    "n_codetected_samples": n_co,
                    "jaccard_observed": round(j_obs, 6),
                    "p_virus_given_host_observed": round(p_vh, 6),
                    "matched_control_mean_jaccard": round(mean_ctrl_j, 6) if np.isfinite(mean_ctrl_j) else "NA",
                    "delta_jaccard_vs_matched_unrelated_hosts": round(delta_j, 6) if np.isfinite(delta_j) else "NA",
                    "empirical_p_vs_matched_controls": round(emp_p, 6) if np.isfinite(emp_p) else "NA",
                    "competitive_vhr_status": "non_estimable_separate_mappings_not_comparable_for_competitive_vhr",
                }
            )

    return pd.DataFrame(cooc_rows)


def ols_residualize(y: np.ndarray, X: np.ndarray) -> np.ndarray:
    """Return OLS residuals of y on [1, X]."""
    n = len(y)
    if n == 0:
        return y
    X_des = np.column_stack([np.ones(n), X])
    beta, _, _, _ = np.linalg.lstsq(X_des, y, rcond=None)
    return y - X_des @ beta


def run_phylogeny_controlled_tests(merged: pd.DataFrame) -> pd.DataFrame:
    """Run phylogeny- and detection-opportunity-controlled associations while preserving
    `NaN` for missing MAG abundances, missing defense annotations, and failed/unknown
    CRISPR caller or alignment statuses (excluding them from regression denominators).
    """
    df = merged.copy()

    # Failed or unscreened CRISPR caller rows have UNKNOWN spacer/cas/link status (`np.nan`), NEVER 0!
    unknown_caller_mask = df["mag_host_disposition"].isin(["Failed_caller_execution", "Unknown_unscreened"])
    if "caller_execution_failed" in df.columns:
        unknown_caller_mask = unknown_caller_mask | (pd.to_numeric(df["caller_execution_failed"], errors="coerce") == 1)

    df["n_verified_spacers_clean"] = pd.to_numeric(df["n_verified_spacers"], errors="coerce")
    df.loc[unknown_caller_mask, "n_verified_spacers_clean"] = np.nan

    df["n_cas_operons_clean"] = pd.to_numeric(df["n_cas_operons"], errors="coerce")
    df.loc[unknown_caller_mask, "n_cas_operons_clean"] = np.nan

    df["has_verified_spacers"] = np.where(
        df["n_verified_spacers_clean"].isna(),
        np.nan,
        (df["n_verified_spacers_clean"] > 0).astype(float),
    )
    df["has_cas_operon"] = np.where(
        df["n_cas_operons_clean"].isna(),
        np.nan,
        (df["n_cas_operons_clean"] > 0).astype(float),
    )

    # For viral link sensitivity status, MAGs with failed/unscreened callers OR unresolved
    # unknown viral identity hits (`Unknown_unverified_viral_identity`) have UNKNOWN link status (`np.nan`), not 0!
    t1 = pd.to_numeric(df.get("n_tier1_verified_direct_links", 0), errors="coerce")
    t2 = pd.to_numeric(df.get("n_tier2_corroborated_links", 0), errors="coerce")
    unknown_link_mask = unknown_caller_mask | (df["mag_host_disposition"] == "Unknown_unverified_viral_identity")
    df["has_verified_viral_link_sens"] = np.where(
        unknown_link_mask | t1.isna() | t2.isna(),
        np.nan,
        ((t1 + t2) > 0).astype(float),
    )

    # Never `.fillna(0.0)` on `c8_defense_systems`! Missing annotation is `np.nan`.
    df["c8_def_num"] = pd.to_numeric(df.get("c8_defense_systems"), errors="coerce")
    df["mean_sample_abundance_num"] = pd.to_numeric(df.get("mean_sample_abundance"), errors="coerce")
    df["checkm2_completeness_num"] = pd.to_numeric(df.get("checkm2_completeness"), errors="coerce")
    df["refined_n_bases_num"] = pd.to_numeric(df.get("refined_n_bases"), errors="coerce")
    df["refined_n_contigs_num"] = pd.to_numeric(df.get("refined_n_contigs"), errors="coerce")

    assoc_rows: list[dict[str, Any]] = []
    for scope_name, scope_df in (
        ("all591_archaea", df),
        ("thermoplasmatota_237", df[df["canonical_phylum"] == "p__Thermoplasmatota"]),
    ):
        n_total = len(scope_df)
        for pred_col in ("has_verified_spacers", "has_cas_operon", "has_verified_viral_link_sens"):
            for resp_col, resp_clean in (
                ("c8_def_num", "c8_def_num"),
                ("n_cas_operons", "n_cas_operons_clean"),
                ("n_verified_spacers", "n_verified_spacers_clean"),
            ):
                if (pred_col == "has_verified_spacers" and resp_col == "n_verified_spacers") or (
                    pred_col == "has_cas_operon" and resp_col == "n_cas_operons"
                ):
                    continue

                miss_abund = scope_df["mean_sample_abundance_num"].isna()
                miss_pred_or_caller = scope_df[pred_col].isna() | scope_df[resp_clean].isna()
                miss_cov = (
                    scope_df["checkm2_completeness_num"].isna()
                    | scope_df["refined_n_bases_num"].isna()
                    | scope_df["refined_n_contigs_num"].isna()
                    | scope_df["canonical_order"].isna()
                )
                valid_mask = ~(miss_abund | miss_pred_or_caller | miss_cov)
                sub = scope_df[valid_mask].copy()

                n_eval = len(sub)
                n_miss_abund = int(miss_abund.sum())
                n_miss_caller_align = int((~miss_abund & miss_pred_or_caller).sum())
                n_miss_cov = int((~miss_abund & ~miss_pred_or_caller & miss_cov).sum())

                if n_eval < 10:
                    assoc_rows.append(
                        {
                            "scope": scope_name,
                            "n_scope_total_mags": n_total,
                            "n_mags": n_eval,
                            "n_missing_abundance_excluded": n_miss_abund,
                            "n_unknown_caller_or_alignment_excluded": n_miss_caller_align,
                            "n_missing_annotation_or_covariate_excluded": n_miss_cov,
                            "predictor": pred_col,
                            "n_predictor_positive": int((sub[pred_col] > 0).sum()) if n_eval > 0 else 0,
                            "response": resp_col,
                            "estimability_status": "non_estimable_insufficient_complete_observations",
                            "unadjusted_spearman_r": "NA",
                            "unadjusted_p": "NA",
                            "phylo_detection_adjusted_partial_r": "NA",
                            "phylo_detection_adjusted_p": "NA",
                        }
                    )
                    continue

                x_vec = sub[pred_col].to_numpy(dtype=float)
                y_vec = sub[resp_clean].to_numpy(dtype=float)
                n_pos = int((x_vec > 0).sum())

                if n_pos < 3 or np.std(x_vec) == 0 or np.std(y_vec) == 0:
                    assoc_rows.append(
                        {
                            "scope": scope_name,
                            "n_scope_total_mags": n_total,
                            "n_mags": n_eval,
                            "n_missing_abundance_excluded": n_miss_abund,
                            "n_unknown_caller_or_alignment_excluded": n_miss_caller_align,
                            "n_missing_annotation_or_covariate_excluded": n_miss_cov,
                            "predictor": pred_col,
                            "n_predictor_positive": n_pos,
                            "response": resp_col,
                            "estimability_status": "non_estimable_insufficient_positive_mags_or_variance",
                            "unadjusted_spearman_r": "NA",
                            "unadjusted_p": "NA",
                            "phylo_detection_adjusted_partial_r": "NA",
                            "phylo_detection_adjusted_p": "NA",
                        }
                    )
                    continue

                X_cov = np.column_stack(
                    [
                        sub["checkm2_completeness_num"].to_numpy(dtype=float),
                        np.log10(np.maximum(1.0, sub["refined_n_bases_num"].to_numpy(dtype=float))),
                        np.log10(np.maximum(1.0, sub["refined_n_contigs_num"].to_numpy(dtype=float))),
                        np.log10(sub["mean_sample_abundance_num"].to_numpy(dtype=float) + 1e-6),
                    ]
                )
                order_counts = sub["canonical_order"].value_counts()
                major_orders = [o for o, c in order_counts.items() if c >= 5]
                if len(major_orders) > 1:
                    order_dummies = np.column_stack(
                        [(sub["canonical_order"] == o).to_numpy(dtype=float) for o in major_orders[:-1]]
                    )
                    X_full = np.column_stack([X_cov, order_dummies])
                else:
                    X_full = X_cov

                r_raw, p_raw = stats.spearmanr(x_vec, y_vec)
                rx = ols_residualize(x_vec, X_full)
                ry = ols_residualize(y_vec, X_full)
                if np.std(rx) == 0 or np.std(ry) == 0:
                    r_adj, p_adj = np.nan, np.nan
                else:
                    r_adj, p_adj = stats.pearsonr(rx, ry)

                assoc_rows.append(
                    {
                        "scope": scope_name,
                        "n_scope_total_mags": n_total,
                        "n_mags": n_eval,
                        "n_missing_abundance_excluded": n_miss_abund,
                        "n_unknown_caller_or_alignment_excluded": n_miss_caller_align,
                        "n_missing_annotation_or_covariate_excluded": n_miss_cov,
                        "predictor": pred_col,
                        "n_predictor_positive": n_pos,
                        "response": resp_col,
                        "estimability_status": "estimable_exploratory_existing_mappings",
                        "unadjusted_spearman_r": round(float(r_raw), 6) if np.isfinite(r_raw) else "NA",
                        "unadjusted_p": round(float(p_raw), 6) if np.isfinite(p_raw) else "NA",
                        "phylo_detection_adjusted_partial_r": round(float(r_adj), 6) if np.isfinite(r_adj) else "NA",
                        "phylo_detection_adjusted_p": round(float(p_adj), 6) if np.isfinite(p_adj) else "NA",
                    }
                )

    return pd.DataFrame(assoc_rows)


def main() -> None:
    t0 = time.time()
    RES_DIR.mkdir(parents=True, exist_ok=True)

    ledger_df = pd.read_csv(RES_DIR / "stage3_all591_archaeal_mag_host_evidence.tsv", sep="\t")
    pairs_path = RES_DIR / "stage3_archaeal_virus_host_pairs_all_tiers.tsv"
    pairs_df = (
        pd.read_csv(pairs_path, sep="\t")
        if (pairs_path.exists() and pairs_path.stat().st_size > 10)
        else pd.DataFrame()
    )

    mag_abund_df = pd.read_csv(P_MAG_ABUND, sep="\t", index_col=0) if P_MAG_ABUND.exists() else pd.DataFrame()
    vir_tpm_df = pd.read_csv(P_VIR_TPM, sep="\t", index_col=0) if P_VIR_TPM.exists() else pd.DataFrame()

    intersect_info = resolve_shared_sample_intersection(mag_abund_df, vir_tpm_df)
    shared_samples = intersect_info["shared_samples"]
    if shared_samples:
        mag_abund_df = mag_abund_df[shared_samples]
        vir_tpm_df = vir_tpm_df[shared_samples]

    abund_metrics_df = compute_host_abundance_metrics(
        list(ledger_df["genome"].astype(str)),
        mag_abund_df,
        sample_cols=shared_samples,
    )
    ledger_df = ledger_df.merge(abund_metrics_df, on="genome", how="left")
    # Persist updated host mapping status & NaN-preserving abundance metrics
    ledger_df.to_csv(
        RES_DIR / "stage4_all591_host_abundance_and_mapping_status.tsv",
        sep="\t",
        index=False,
        na_rep="NA",
    )

    eligible_pairs = (
        pairs_df[pairs_df["eligible_stage4_sensitivity"] == 1].copy()
        if len(pairs_df) > 0
        else pd.DataFrame()
    )
    target_votus = set(eligible_pairs["votu_id"].astype(str)) if len(eligible_pairs) > 0 else set()
    e1_mask, b75_mask, e1_obs_samples, e1_fail_samples = load_e1_breadth_masks(target_votus)

    cooc_df = evaluate_pair_cooccurrence(
        eligible_pairs=eligible_pairs,
        ledger_df=ledger_df,
        mag_abund_df=mag_abund_df,
        vir_tpm_df=vir_tpm_df,
        shared_samples=shared_samples,
        e1_mask=e1_mask,
        b75_mask=b75_mask,
        e1_observed_samples=e1_obs_samples if P_E1_TAR.exists() else set(),
        e1_failed_samples=e1_fail_samples,
    )
    cooc_df.to_csv(
        RES_DIR / "stage4_virus_host_cooccurrence_and_controls.tsv",
        sep="\t",
        index=False,
        na_rep="NA",
    )

    ver_df = pd.read_csv(RES_DIR / "stage0_arch591_fasta_verification.tsv", sep="\t")
    merged = ledger_df.merge(
        ver_df[["genome", "refined_n_contigs", "refined_n_bases"]],
        on="genome",
        how="left",
    )

    arch_annot_path = P_ARCH_IN / "archaea.tsv"
    if arch_annot_path.exists():
        aadf = pd.read_csv(arch_annot_path, sep="\t")
        gcol = "genome" if "genome" in aadf.columns else aadf.columns[0]
        num_cols = [
            c
            for c in aadf.columns
            if c != gcol and c not in merged.columns and pd.api.types.is_numeric_dtype(aadf[c])
        ]
        if num_cols:
            merged = merged.merge(aadf[[gcol] + num_cols], left_on="genome", right_on=gcol, how="left")

    assoc_df = run_phylogeny_controlled_tests(merged)
    assoc_df.to_csv(
        RES_DIR / "stage4_surface_defense_phylogeny_controlled_tests.tsv",
        sep="\t",
        index=False,
        na_rep="NA",
    )

    # Explicit conditional follow-up status table (`EXP3-RUNTIME-HANDOFF.txt` & `EXP3-FINAL-PROVENANCE-REVIEW.txt`)
    n_t1 = int((eligible_pairs["eligible_stage4_primary"] == 1).sum()) if len(eligible_pairs) > 0 else 0
    followup_rows = [
        {
            "module": "existing_separate_mapping_cooccurrence_and_breadth_gates",
            "status": "implemented_in_stage4_exploratory_only",
            "launched_as_slurm_job": "yes_stage4",
            "evidence_scope": "existing_separate_mag_and_votu_mappings_not_joint_competitive_fastq_mapping",
            "reason_or_gate": "Exploratory evaluation of Jaccard and P(virus|host) across unrestricted, E1 breadth>=0.50, and strict breadth>=0.75 vs phylum/prevalence-matched unrelated archaeal host controls; missing MAGs/samples excluded as NA.",
        },
        {
            "module": "phylogeny_and_detection_opportunity_controlled_defense_associations",
            "status": "implemented_in_stage4_exploratory_only",
            "launched_as_slurm_job": "yes_stage4",
            "evidence_scope": "all591_archaea_and_237_thermoplasmatota_with_gtdb_order_and_detection_covariates",
            "reason_or_gate": "Controls for canonical GTDB order, CheckM2 completeness, log10(bases), log10(contigs), and log10(mean_abundance); excludes missing MAG abundances and failed/unknown callers as NA.",
        },
        {
            "module": "joint_competitive_virus_host_fastq_read_recruitment_and_junction_spanning",
            "status": "not_launched_conditional_on_verified_tier1_pairs",
            "launched_as_slurm_job": "no_conditional_not_launched",
            "evidence_scope": f"conditional_on_stage3_tier1_verified_pairs (current_n_tier1={n_t1})",
            "reason_or_gate": "Not launched: separate pre-existing MAG and vOTU mappings cannot be used to claim competitive host/virus abundance or read-junction spanning; targeted joint Bowtie2 competitive FASTQ mapping is conditional on verified Tier 1 pairs from Stage 3.",
        },
        {
            "module": "instrain_coverage_qualified_microdiversity_and_recombination_at_defense_surface_loci",
            "status": "not_launched_non_estimable_without_joint_competitive_bams_meeting_10x_50pct_breadth",
            "launched_as_slurm_job": "no_conditional_not_launched",
            "evidence_scope": "requires_competitive_bam_depth_ge_10x_and_breadth_ge_0.50",
            "reason_or_gate": "Not launched: allelic microdiversity (pN/pS, pi) and recombination at surface/defense loci cannot be estimated from summary abundance matrices and require competitive BAMs meeting >=10x depth and >=50% breadth.",
        },
        {
            "module": "structural_foldseek_colabfold_archaeal_viral_hallmark_screen",
            "status": "not_launched_conditional_on_stage2_unannotated_candidate_proteins",
            "launched_as_slurm_job": "no_conditional_not_launched",
            "evidence_scope": "sequence_level_fresh_checkv_v1.5_and_genomad_v1.12_implemented_in_stage2",
            "reason_or_gate": "Not launched: Stage 2 runs fresh CheckV + geNomad sequence/HMM annotation on all matched primary targets; structural Foldseek/ColabFold on unannotated dark viral ORFs is conditional on Stage 2 output.",
        },
        {
            "module": "iphop_auxiliary_classifier_corroboration",
            "status": "external_ibex_job_55361833_running_auxiliary_only",
            "launched_as_slurm_job": "external_ibex_55361833",
            "evidence_scope": "auxiliary_genus_level_prediction_never_species_level_proof",
            "reason_or_gate": "iPHoP RaFAH stage running on IBEX (job 55361833); per contract, classifier predictions are auxiliary corroboration only and never substitute for direct CRISPR/provirus evidence.",
        },
    ]
    pd.DataFrame(followup_rows).to_csv(
        RES_DIR / "stage4_conditional_followup_status.tsv", sep="\t", index=False
    )

    n_missing_mag_abund = int((abund_metrics_df["host_mapping_status"] == "unknown_missing_from_mag_abundance_matrix").sum())
    n_obs_all_zero_mag = int((abund_metrics_df["host_mapping_status"] == "observed_all_zero_in_mag_abundance_matrix").sum())
    n_obs_pos_mag = int((abund_metrics_df["host_mapping_status"] == "observed_in_mag_abundance_matrix").sum())

    summary = {
        "stage": "stage4_downstream_mechanism",
        "status": "completed",
        "elapsed_seconds": round(time.time() - t0, 2),
        "mapping_evidence_scope": "existing_separate_mag_and_votu_mappings_exploratory_only_no_new_competitive_fastq_or_variant_evidence",
        "sample_intersection": intersect_info,
        "host_abundance_mapping_counts": {
            "n_observed_positive": n_obs_pos_mag,
            "n_observed_all_zero": n_obs_all_zero_mag,
            "n_unknown_missing_from_matrix": n_missing_mag_abund,
        },
        "n_eligible_primary_pairs": n_t1,
        "n_eligible_sensitivity_pairs": int(len(eligible_pairs)),
        "n_cooccurrence_evaluations": int(len(cooc_df)),
        "n_phylogeny_controlled_tests": int(len(assoc_df)),
        "conditional_followup_modules": {r["module"]: r["status"] for r in followup_rows},
    }
    (RES_DIR / "stage4_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print("[Stage 4] Summary:", json.dumps(summary), flush=True)


if __name__ == "__main__":
    main()
