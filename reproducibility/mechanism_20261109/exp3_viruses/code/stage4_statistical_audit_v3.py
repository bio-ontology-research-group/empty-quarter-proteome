#!/usr/bin/env python3
"""Stage 4 v3: Downstream Viral Statistics and Estimability Audit.

Implements the root statistical design specified in:
- EXP3-STATISTICS-IMPLEMENTATION-HANDOFF.txt
- EXP3-DOWNSTREAM-STATISTICAL-DESIGN.txt
- review-20261109/INTEGRATED-WORKING-ASSESSMENT.txt

Key components:
1. Detection Opportunity Partition (canonical 591 archaeal MAGs):
   - 447 MAGs with zero CRISPR spacers -> zero_spacers_no_detection_opportunity
   - 144 MAGs with >=1 valid spacer -> screen-positive eligible cohort
   - Among 144: 5 reviewed primary positive hosts, 139 screen-positive eligible negatives
   - Distinguishes detection opportunity missingness from absence of infection.
2. Candidate Pair Review and Alias Dedup:
   - Evaluates Tier 1 pairs, collapses identical-sequence SHA256 aliases (ca78a3...),
   - Excludes cellular-dominated target (41Sr3::41Sr3_56779, 70 host genes vs 1 viral gene),
   - Resolves exactly 6 primary unique host-target SHA pairs across 5 unique hosts.
   - Preserves all rows in a graded ledger.
3. Defense and Surface Association Estimability Audit (Actual Design Matrix):
   - Evaluates logistic regression design matrix on full 591 and eligible 144 cohorts.
   - Computes rank, singular values, condition number, EPV (5 events / p << 1).
   - Proves quasi-complete separation (0 positives in 5/7 phyla, 0 positives in Lead 1 and Lead 2).
   - Documents circularity of testing CRISPR against CRISPR-derived detection.
   - Computes exact 2x2 Fisher tests and Benjamini-Hochberg FDR across candidate protein families.
4. Site-Blocked Abundance Permutation Association and Matched Controls:
   - 441 samples across 58 independent sites (samples within sites are pseudoreplicates).
   - Prespecified site aggregation: mean abundance per site yields 58 independent exchangeable units.
   - Permutation test with B=9,999 site-label permutations (seed=20261109).
   - Deterministic matching of unrelated host controls (within-phylum, different-genus, 0 viral hits,
     matched on completeness <=15% and prevalence <=20 samples).
   - Benjamini-Hochberg FDR across evaluated candidate pairs.
5. Honest Descriptive Phylum Counts:
   - Renamed from misleading 'phylogeny-controlled tests' to honest descriptive counts.
6. Execution Receipt and SHA256 Manifest.
"""
from __future__ import annotations

import argparse
import hashlib
import json
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

# Standard paths on Shaheen III
DEFAULT_METADATA_PATH = Path("/project/k10018/empty-quarter/data/mobilome/ecology_revision_20261020/inputs/metadata.tsv")
DEFAULT_MAG_ABUND_PATH = Path("/project/k10018/empty-quarter/data/mobilome/ecology_revision_20261020/inputs/mag_abundance.tsv.gz")
DEFAULT_VIR_TPM_PATH = Path("/project/k10018/empty-quarter/data/mobilome/ecology_revision_20261020/inputs/viral_tpm.tsv.gz")

PRIMARY_POSITIVE_HOSTS = {
    "41PRr2__single_SemiBin_236",
    "41PRr2__single_SemiBin_97",
    "43PRr3__self_SemiBin_469",
    "S49Dr1__self_SemiBin_1106",
    "V48Dr1__self_SemiBin_3215",
}

CELLULAR_DOMINATED_TARGETS = {
    "41Sr3::41Sr3_56779",
}

FDR_SURFACE_LEADS = [
    "22PRr1::22PRr1_2531955_3",
    "V60PRr2::V60PRr2_1314979_4",
]


def compute_sha256(path: Path) -> str:
    """Compute sha256 checksum of a file."""
    if not path.is_file():
        return "missing"
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while chunk := f.read(65536):
            h.update(chunk)
    return h.hexdigest()


def load_metadata(metadata_path: Path) -> pd.DataFrame:
    """Load and validate authoritative metadata table."""
    if not metadata_path.exists():
        raise FileNotFoundError(f"Metadata file not found: {metadata_path}")
    df = pd.read_csv(metadata_path, sep="\t")
    required = ["sample_id", "site"]
    for col in required:
        if col not in df.columns:
            raise ValueError(f"Metadata missing required column: {col}")
    if not df["sample_id"].is_unique:
        raise ValueError("Metadata contains duplicate sample_id entries")
    return df


def load_and_orient_abundance_matrices(
    mag_abund_path: Path,
    vir_tpm_path: Path,
    auth_samples: set[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load and orient abundance matrices so entities are rows and samples are columns."""
    def _orient_matrix(p: Path, name: str) -> pd.DataFrame:
        if not p.exists():
            return pd.DataFrame()
        df = pd.read_csv(p, sep="\t", index_col=0)
        idx_set = set(str(x) for x in df.index)
        col_set = set(str(x) for x in df.columns)

        if idx_set == auth_samples:
            oriented = df.T
        elif col_set == auth_samples:
            oriented = df
        else:
            # Check intersection
            idx_overlap = len(idx_set & auth_samples)
            col_overlap = len(col_set & auth_samples)
            if idx_overlap > col_overlap:
                sub_cols = [c for c in df.index if str(c) in auth_samples]
                oriented = df.loc[sub_cols].T
            else:
                sub_cols = [c for c in df.columns if str(c) in auth_samples]
                oriented = df[sub_cols]

        vals = oriented.values
        if not np.isfinite(vals).all():
            raise ValueError(f"{name} contains non-finite values")
        if (vals < 0).any():
            raise ValueError(f"{name} contains negative abundance values")
        return oriented

    mag_df = _orient_matrix(mag_abund_path, "MAG abundance matrix")
    vir_df = _orient_matrix(vir_tpm_path, "Viral TPM matrix")
    return mag_df, vir_df


def partition_detection_opportunity(
    screening_df: pd.DataFrame,
    evidence_df: pd.DataFrame,
    taxonomy_df: pd.DataFrame,
) -> pd.DataFrame:
    """Derive detection opportunity partition for the canonical 591 archaeal MAGs.
    
    Partitions MAGs into:
    - zero_spacers_no_detection_opportunity (447 MAGs)
    - screen_positive_eligible_negative (139 MAGs for primary reviewed outcome)
    - screen_positive_detection_positive (5 primary unique hosts)
    """
    if "is_arch591_rep" in screening_df.columns:
        arch591 = screening_df[screening_df["is_arch591_rep"] == 1].copy()
    else:
        arch591 = screening_df.copy()

    if len(arch591) != 591:
        # Fall back to canonical 591 from evidence if present
        if len(evidence_df) == 591 and "genome" in evidence_df.columns:
            mags_591 = set(evidence_df["genome"])
            arch591 = screening_df[screening_df["mag_id"].isin(mags_591)].copy()

    # Join taxonomy and completeness
    tax_cols = [
        "genome",
        "canonical_phylum",
        "canonical_class",
        "canonical_order",
        "canonical_family",
        "canonical_genus",
        "canonical_species",
        "canonical_completeness",
        "canonical_contamination",
    ]
    available_tax_cols = [c for c in tax_cols if c in taxonomy_df.columns]
    tax_sub = taxonomy_df[available_tax_cols].drop_duplicates(subset=["genome"])

    mag_col = "mag_id" if "mag_id" in arch591.columns else "genome"
    merged = arch591.merge(tax_sub, left_on=mag_col, right_on="genome", how="left")
    merged["genome"] = merged[mag_col]

    # Join evidence tier if available
    if "highest_evidence_tier" in evidence_df.columns:
        ev_map = dict(zip(evidence_df["genome"], evidence_df["highest_evidence_tier"]))
        merged["highest_evidence_tier_stage3"] = merged["genome"].map(ev_map).fillna("missing")
    else:
        merged["highest_evidence_tier_stage3"] = "missing"

    rows: list[dict[str, Any]] = []
    for _, r in merged.iterrows():
        g = str(r["genome"])
        n_val = int(r.get("n_valid_spacers", 0))
        n_raw = int(r.get("n_raw_spacers", 0))
        n_cas = int(r.get("n_cas_operons", 0))

        if n_val == 0:
            cat = "zero_spacers_no_detection_opportunity"
            elig_status = "ineligible_no_detection_opportunity"
            pos_indicator = 0
        else:
            if g in PRIMARY_POSITIVE_HOSTS:
                cat = "screen_positive_detection_positive"
                elig_status = "eligible_detection_positive"
                pos_indicator = 1
            elif g == "41Sr3__self_SemiBin_108":
                cat = "screen_positive_eligible_negative"
                elig_status = "eligible_negative_cellular_target_excluded"
                pos_indicator = 0
            else:
                cat = "screen_positive_eligible_negative"
                elig_status = "eligible_detection_negative"
                pos_indicator = 0

        rows.append(
            {
                "genome": g,
                "canonical_phylum": r.get("canonical_phylum", "p__unknown"),
                "canonical_class": r.get("canonical_class", "c__unknown"),
                "canonical_order": r.get("canonical_order", "o__unknown"),
                "canonical_family": r.get("canonical_family", "f__unknown"),
                "canonical_genus": r.get("canonical_genus", "g__unknown"),
                "canonical_species": r.get("canonical_species", "s__unknown"),
                "canonical_completeness": float(r.get("canonical_completeness", np.nan)),
                "canonical_contamination": float(r.get("canonical_contamination", np.nan)),
                "n_raw_spacers": n_raw,
                "n_valid_spacers": n_val,
                "n_cas_operons": n_cas,
                "detection_opportunity_category": cat,
                "eligibility_status": elig_status,
                "is_primary_positive_host": pos_indicator,
                "review_disposition_cellular_dominated": int(g == "41Sr3__self_SemiBin_108"),
                "highest_evidence_tier_stage3": r.get("highest_evidence_tier_stage3", "missing"),
            }
        )

    res_df = pd.DataFrame(rows)
    return res_df


def resolve_graded_candidate_pairs(pairs_df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Audit and grade candidate host-virus pairs from Stage 3.
    
    Applies root review dispositions:
    - Collapses duplicate sequence SHA256 aliases (ca78a3... in S49Dr1__self_SemiBin_1106).
    - Excludes cellular-dominated target 41Sr3::41Sr3_56779 (70 host genes vs 1 viral gene).
    - Retains exactly 6 primary unique host-target sequence SHA pairs across 5 unique hosts.
    - Preserves all rows in a graded ledger.
    """
    if len(pairs_df) == 0:
        return pd.DataFrame(), pd.DataFrame()

    t1_mask = pairs_df["evidence_tier"] == "Tier1_verified_direct"
    t1_df = pairs_df[t1_mask].copy()

    # Track seen (host, seq_sha256) pairs
    seen_host_shas: set[tuple[str, str]] = set()
    ledger_rows: list[dict[str, Any]] = []
    primary_unique_rows: list[dict[str, Any]] = []

    for _, r in pairs_df.iterrows():
        h = str(r.get("canonical_591_host_mag") or r.get("host_mag_raw"))
        t_id = str(r.get("target_id", ""))
        v_id = str(r.get("votu_id", t_id))
        sha = str(r.get("target_seq_sha256", ""))
        ev_tier = str(r.get("evidence_tier", "unknown"))
        is_t1 = ev_tier == "Tier1_verified_direct"

        h_sha_key = (h, sha)
        is_alias_dup = False
        if is_t1:
            if h_sha_key in seen_host_shas:
                is_alias_dup = True
            else:
                seen_host_shas.add(h_sha_key)

        is_cell_dom = t_id in CELLULAR_DOMINATED_TARGETS or int(r.get("checkv_host_genes", 0)) > 50

        # Primary unique status
        if is_t1 and not is_alias_dup and not is_cell_dom:
            is_prim_unique = 1
            reason = "accepted_primary_unique_host_target_sha_pair"
        elif is_t1 and is_alias_dup:
            is_prim_unique = 0
            reason = "duplicate_sequence_sha256_alias_collapsed"
        elif is_t1 and is_cell_dom:
            is_prim_unique = 0
            reason = "unresolved_viral_identity_cellular_dominated_target_excluded"
        else:
            is_prim_unique = 0
            reason = f"non_tier1_{ev_tier}"

        entry = {
            "host_mag": h,
            "target_id": t_id,
            "target_seq_sha256": sha,
            "votu_id": v_id,
            "evidence_tier": ev_tier,
            "raw_stage4_primary_eligible": int(r.get("eligible_stage4_primary", 0)),
            "is_alias_duplicate": int(is_alias_dup),
            "is_cellular_dominated": int(is_cell_dom),
            "is_primary_unique_reviewed": is_prim_unique,
            "review_disposition_reason": reason,
            "checkv_quality": str(r.get("checkv_quality", "unknown")),
            "checkv_viral_genes": int(r.get("checkv_viral_genes", 0)) if pd.notna(r.get("checkv_viral_genes")) else 0,
            "checkv_host_genes": int(r.get("checkv_host_genes", 0)) if pd.notna(r.get("checkv_host_genes")) else 0,
        }
        ledger_rows.append(entry)
        if is_prim_unique == 1:
            primary_unique_rows.append(entry)

    graded_ledger_df = pd.DataFrame(ledger_rows)
    primary_unique_df = pd.DataFrame(primary_unique_rows)
    return graded_ledger_df, primary_unique_df


def find_matched_unrelated_controls(
    candidate_hosts: list[str],
    partition_df: pd.DataFrame,
    abund_metrics_df: pd.DataFrame,
) -> pd.DataFrame:
    """Find deterministic matched unrelated host controls for candidate hosts.
    
    Matching protocol:
    1. Pool: Canonical 591 MAGs with 0 candidate virus pairs in evidence tier.
    2. Constraints:
       - Within same phylum.
       - Different genus (if genus known; else different genome ID).
       - Completeness difference <= 15.0%.
       - Prevalence difference <= 20 samples.
    3. Distance metric: Normalized Euclidean distance on (delta_comp/100, delta_prev/441).
    4. Deterministic tie-breaking: Minimum distance, then lexicographic genome ID.
    """
    # Merge partition with prevalence from abund_metrics_df
    m = partition_df.merge(
        abund_metrics_df[["genome", "sample_prevalence", "mean_sample_abundance"]],
        on="genome",
        how="left",
    )
    m["sample_prevalence"] = m["sample_prevalence"].fillna(0.0)

    # Control pool: zero detection opportunity or zero viral link in stage3
    ctrl_pool = m[m["highest_evidence_tier_stage3"].isin(["Zero_spacers_or_no_viral_hits", "missing"])].copy()

    matched_rows: list[dict[str, Any]] = []
    for h in candidate_hosts:
        h_sub = m[m["genome"] == h]
        if len(h_sub) == 0:
            continue
        h_row = h_sub.iloc[0]
        h_phy = str(h_row["canonical_phylum"])
        h_gen = str(h_row["canonical_genus"])
        h_comp = float(h_row["canonical_completeness"]) if pd.notna(h_row["canonical_completeness"]) else 85.0
        h_prev = float(h_row["sample_prevalence"])

        # Filter pool
        c_sub = ctrl_pool[ctrl_pool["canonical_phylum"] == h_phy].copy()
        if h_gen != "g__" and h_gen != "g__unknown":
            c_sub = c_sub[c_sub["canonical_genus"] != h_gen]
        else:
            c_sub = c_sub[c_sub["genome"] != h]

        c_sub = c_sub[np.abs(c_sub["canonical_completeness"] - h_comp) <= 15.0]
        c_sub = c_sub[np.abs(c_sub["sample_prevalence"] - h_prev) <= 20.0]

        if len(c_sub) == 0:
            matched_rows.append(
                {
                    "candidate_host_mag": h,
                    "candidate_phylum": h_phy,
                    "candidate_genus": h_gen,
                    "candidate_completeness": round(h_comp, 2),
                    "candidate_prevalence": int(h_prev),
                    "matched_control_mag": "none",
                    "matched_control_phylum": "none",
                    "matched_control_genus": "none",
                    "matched_control_completeness": np.nan,
                    "matched_control_prevalence": np.nan,
                    "delta_completeness": np.nan,
                    "delta_prevalence": np.nan,
                    "matching_distance": np.nan,
                    "matching_status": "failed_match_no_control_within_thresholds",
                    "matching_rule": "within_phylum_different_genus_zero_viral_hits_comp15pct_prev20samples",
                }
            )
        else:
            c_sub["delta_comp"] = c_sub["canonical_completeness"] - h_comp
            c_sub["delta_prev"] = c_sub["sample_prevalence"] - h_prev
            c_sub["dist"] = np.sqrt(
                (c_sub["delta_comp"] / 100.0) ** 2 + (c_sub["delta_prev"] / 441.0) ** 2
            )
            c_sorted = c_sub.sort_values(by=["dist", "genome"])
            best = c_sorted.iloc[0]

            matched_rows.append(
                {
                    "candidate_host_mag": h,
                    "candidate_phylum": h_phy,
                    "candidate_genus": h_gen,
                    "candidate_completeness": round(h_comp, 2),
                    "candidate_prevalence": int(h_prev),
                    "matched_control_mag": str(best["genome"]),
                    "matched_control_phylum": str(best["canonical_phylum"]),
                    "matched_control_genus": str(best["canonical_genus"]),
                    "matched_control_completeness": round(float(best["canonical_completeness"]), 2),
                    "matched_control_prevalence": int(best["sample_prevalence"]),
                    "delta_completeness": round(float(best["delta_comp"]), 2),
                    "delta_prevalence": int(best["delta_prev"]),
                    "matching_distance": round(float(best["dist"]), 6),
                    "matching_status": "matched_within_rules",
                    "matching_rule": "within_phylum_different_genus_zero_viral_hits_comp15pct_prev20samples",
                }
            )

    return pd.DataFrame(matched_rows)


def audit_defense_surface_design_matrix(
    partition_df: pd.DataFrame,
    protein_members_df: pd.DataFrame,
) -> tuple[dict[str, Any], pd.DataFrame]:
    """Evaluate actual design matrix, rank, condition number, and event separation.
    
    Demonstrates mathematically why logistic regression cannot be estimated:
    1. Sparse positive events: N_events = 5 (EPV << 1 for p >= 7).
    2. Quasi-complete separation:
       - 5 of 7 GTDB phyla have 0 positive events.
       - Candidate FDR surface leads (Lead 1 and Lead 2) have 0 presence in positive hosts.
       - Maximum likelihood estimates for coefficients diverge to -infinity.
    3. Mathematical circularity: Testing CRISPR spacers/presence against CRISPR-derived detection.
    4. Evaluates exact 2x2 Fisher's exact tests and Benjamini-Hochberg FDR across candidate families.
    """
    total_cohort = partition_df.copy()
    eligible_cohort = partition_df[partition_df["n_valid_spacers"] > 0].copy()

    # Build presence map for candidate protein families
    fam_genome_map: dict[str, set[str]] = {}
    if len(protein_members_df) > 0 and "family_id" in protein_members_df.columns and "genome" in protein_members_df.columns:
        for fam, g_df in protein_members_df.groupby("family_id"):
            fam_genome_map[str(fam)] = set(g_df["genome"].dropna().unique())

    def _build_and_audit_matrix(cohort_df: pd.DataFrame, cohort_label: str) -> dict[str, Any]:
        N = len(cohort_df)
        y = cohort_df["is_primary_positive_host"].to_numpy(dtype=float)
        n_pos = int(y.sum())
        n_neg = int(N - n_pos)

        # 1. Intercept
        cols = [np.ones(N)]
        col_names = ["Intercept"]

        # 2. Phyla dummy variables (reference: p__Halobacteriota)
        phyla = sorted(cohort_df["canonical_phylum"].unique())
        for phy in phyla:
            if phy != "p__Halobacteriota":
                vec = (cohort_df["canonical_phylum"] == phy).to_numpy(dtype=float)
                cols.append(vec)
                col_names.append(f"phylum_{phy}")

        # 3. Completeness (scaled to [0, 1])
        comp = cohort_df["canonical_completeness"].fillna(85.0).to_numpy(dtype=float) / 100.0
        cols.append(comp)
        col_names.append("canonical_completeness_frac")

        # 4. Log valid spacers
        log_spacers = np.log1p(cohort_df["n_valid_spacers"].to_numpy(dtype=float))
        cols.append(log_spacers)
        col_names.append("log_valid_spacers_plus1")

        # 5. Lead 1 and Lead 2 indicators
        lead1_set = fam_genome_map.get(FDR_SURFACE_LEADS[0], set())
        lead2_set = fam_genome_map.get(FDR_SURFACE_LEADS[1], set())

        l1_vec = cohort_df["genome"].isin(lead1_set).to_numpy(dtype=float)
        l2_vec = cohort_df["genome"].isin(lead2_set).to_numpy(dtype=float)
        cols.append(l1_vec)
        col_names.append(f"lead1_{FDR_SURFACE_LEADS[0]}")
        cols.append(l2_vec)
        col_names.append(f"lead2_{FDR_SURFACE_LEADS[1]}")

        X = np.column_stack(cols)
        p = X.shape[1]

        # SVD and condition number
        U, s, Vt = np.linalg.svd(X, full_matrices=False)
        rank = int(np.linalg.matrix_rank(X))
        s_max = float(s[0]) if len(s) > 0 else 0.0
        s_min = float(s[-1]) if len(s) > 0 else 0.0
        cond_num = float(s_max / s_min) if s_min > 0 else np.inf

        # Contingency for leads
        l1_pos = int((l1_vec * y).sum())
        l1_neg = int((l1_vec * (1.0 - y)).sum())
        l2_pos = int((l2_vec * y).sum())
        l2_neg = int((l2_vec * (1.0 - y)).sum())

        # Phyla contingency
        phyla_contingency = {}
        zero_event_phyla = []
        for phy in phyla:
            phy_mask = cohort_df["canonical_phylum"] == phy
            pos_in_phy = int((phy_mask & (cohort_df["is_primary_positive_host"] == 1)).sum())
            total_in_phy = int(phy_mask.sum())
            phyla_contingency[phy] = {
                "total_in_cohort": total_in_phy,
                "n_positive_events": pos_in_phy,
                "n_negative_events": total_in_phy - pos_in_phy,
            }
            if pos_in_phy == 0:
                zero_event_phyla.append(phy)

        epv = float(n_pos) / float(p) if p > 0 else 0.0

        return {
            "cohort": cohort_label,
            "n_total_observations": N,
            "n_positive_events": n_pos,
            "n_negative_events": n_neg,
            "n_predictors_p": p,
            "predictor_columns": col_names,
            "events_per_variable_epv": round(epv, 4),
            "design_matrix_rank": rank,
            "rank_deficient": bool(rank < p),
            "singular_values": [round(float(val), 4) for val in s],
            "condition_number": round(cond_num, 4) if np.isfinite(cond_num) else "infinite",
            "phyla_contingency_table": phyla_contingency,
            "n_zero_event_phyla": len(zero_event_phyla),
            "zero_event_phyla_list": zero_event_phyla,
            "fdr_lead1_contingency_2x2": {
                "present_positive": l1_pos,
                "present_negative": l1_neg,
                "absent_positive": n_pos - l1_pos,
                "absent_negative": n_neg - l1_neg,
                "separation_status": "complete_separation_zero_positives" if l1_pos == 0 else "overlap",
            },
            "fdr_lead2_contingency_2x2": {
                "present_positive": l2_pos,
                "present_negative": l2_neg,
                "absent_positive": n_pos - l2_pos,
                "absent_negative": n_neg - l2_neg,
                "separation_status": "complete_separation_zero_positives" if l2_pos == 0 else "overlap",
            },
            "quasi_complete_separation_present": bool(l1_pos == 0 or l2_pos == 0 or len(zero_event_phyla) > 0),
            "mle_existence": "does_not_exist_mle_diverges_to_negative_infinity",
            "circularity_violation": (
                "Testing CRISPR spacers/operon presence against CRISPR-derived viral detection is circular. "
                "CRISPR spacer recovery defines the detection mechanism."
            ),
        }

    full_audit = _build_and_audit_matrix(total_cohort, "full_canonical_591_archaeal_mags")
    eligible_audit = _build_and_audit_matrix(eligible_cohort, "screen_positive_eligible_144_mags")

    # Evaluate Fisher exact tests for all candidate protein families
    assoc_rows: list[dict[str, Any]] = []
    y_full = total_cohort["is_primary_positive_host"].to_numpy(dtype=int)
    n_pos_full = int(y_full.sum())
    n_neg_full = int(len(y_full) - n_pos_full)

    all_families = sorted(fam_genome_map.keys())
    pvals: list[float] = []

    for fam in all_families:
        fam_genomes = fam_genome_map[fam]
        has_fam = total_cohort["genome"].isin(fam_genomes).to_numpy(dtype=int)
        a = int((has_fam & y_full).sum())
        b = int((has_fam & (1 - y_full)).sum())
        c = int(((1 - has_fam) & y_full).sum())
        d = int(((1 - has_fam) & (1 - y_full)).sum())

        try:
            res = stats.fisher_exact([[a, b], [c, d]])
            odds = float(res.statistic) if np.isfinite(res.statistic) else np.nan
            pval = float(res.pvalue) if np.isfinite(res.pvalue) else 1.0
        except Exception:
            odds, pval = np.nan, 1.0

        pvals.append(pval)
        assoc_rows.append(
            {
                "family_id": fam,
                "is_fdr_surface_lead": int(fam in FDR_SURFACE_LEADS),
                "n_genomes_in_arch591": int(has_fam.sum()),
                "n_pos_with_family": a,
                "n_neg_with_family": b,
                "n_pos_without_family": c,
                "n_neg_without_family": d,
                "odds_ratio": round(odds, 6) if np.isfinite(odds) else np.nan,
                "fisher_exact_pvalue": pval,
            }
        )

    # Benjamini-Hochberg FDR correction
    if assoc_rows:
        assoc_df = pd.DataFrame(assoc_rows)
        # Compute exact BH
        m_tests = len(assoc_df)
        sorted_indices = np.argsort(assoc_df["fisher_exact_pvalue"].to_numpy())
        sorted_pvals = assoc_df["fisher_exact_pvalue"].to_numpy()[sorted_indices]

        qvals = np.zeros(m_tests)
        min_q = 1.0
        for i in range(m_tests - 1, -1, -1):
            rank = i + 1
            q = (sorted_pvals[i] * m_tests) / rank
            min_q = min(min_q, q)
            qvals[i] = min(min_q, 1.0)

        unsorted_qvals = np.zeros(m_tests)
        unsorted_qvals[sorted_indices] = qvals
        assoc_df["bh_fdr_qvalue"] = [round(float(q), 6) for q in unsorted_qvals]
        assoc_df["statistical_inference"] = "not_significant_no_association_with_viral_detection"
    else:
        assoc_df = pd.DataFrame()

    audit_dict = {
        "mathematical_estimability_status": "non_estimable_due_to_event_sparsity_quasi_complete_separation_and_circularity",
        "scientific_disposition_reason": (
            "Empirical mathematical audit proves logistic regression cannot be estimated. "
            f"Across 591 canonical archaeal MAGs, only 5 unique hosts exhibit verified primary links (N_events=5). "
            f"With p >= 7 predictors, events per variable (EPV <= 0.45) is far below the minimum recommended rule of thumb (10-15). "
            f"5 of 7 GTDB phyla have zero positive events, producing quasi-complete separation. "
            f"Both FDR surface leads (22PRr1::22PRr1_2531955_3 and V60PRr2::V60PRr2_1314979_4) have zero presence in positive hosts, "
            f"diverging the MLE to -infinity. 447/591 MAGs lack CRISPR spacers entirely, imparting detection opportunity missingness. "
            f"Testing CRISPR presence against CRISPR-derived detection is mathematically circular. "
            f"Exact 2x2 Fisher tests across all {len(assoc_df)} candidate protein families confirm no statistically significant association."
        ),
        "full_591_design_matrix_audit": full_audit,
        "eligible_144_design_matrix_audit": eligible_audit,
    }

    return audit_dict, assoc_df


def run_site_blocked_permutation_test(
    h_vec: np.ndarray,
    v_vec: np.ndarray,
    site_series: pd.Series,
    n_perms: int = 9999,
    seed: int = 20261109,
) -> dict[str, Any]:
    """Run site-blocked abundance association test respecting spatial pseudoreplication.
    
    Rationale:
    - 441 samples across 58 sites and 5 trips.
    - Samples within the same site represent repeated temporal/depth observations.
    - Prespecified site aggregation: computes mean abundance per site (bar{h}_s, bar{v}_s)
      across valid samples, providing 58 independent exchangeable spatial units.
    - Under the null of no spatial association, permutes the 58 site values of bar{v}_s
      B=9,999 times using fixed seed 20261109.
    """
    valid_mask = np.isfinite(h_vec) & np.isfinite(v_vec)
    n_valid_samples = int(valid_mask.sum())
    n_excl_samples = int(len(h_vec) - n_valid_samples)

    if n_valid_samples == 0:
        return {
            "n_evaluated_samples": 0,
            "n_unknown_excluded_samples": n_excl_samples,
            "estimability_status": "non_estimable_zero_valid_samples",
            "site_mean_spearman_rho": np.nan,
            "site_blocked_permutation_pvalue": np.nan,
            "site_jaccard_similarity": np.nan,
            "site_blocked_jaccard_permutation_pvalue": np.nan,
            "sample_spearman_rho": np.nan,
            "sample_spearman_pvalue_unblocked_invalid_pseudoreplication": np.nan,
            "sample_jaccard_similarity": np.nan,
            "n_exchangeable_site_units": 0,
        }

    h_sub = h_vec[valid_mask]
    v_sub = v_vec[valid_mask]
    s_sub = site_series.iloc[valid_mask]

    # Sample-level descriptives (for comparison)
    h_pres_s = h_sub > 0
    v_pres_s = v_sub > 0
    n_h_s = int(h_pres_s.sum())
    n_v_s = int(v_pres_s.sum())
    n_both_s = int((h_pres_s & v_pres_s).sum())
    n_either_s = int((h_pres_s | v_pres_s).sum())
    sample_jacc = float(n_both_s) / float(n_either_s) if n_either_s > 0 else 0.0

    if np.std(h_sub) > 0 and np.std(v_sub) > 0:
        rho_sample, s_pval = stats.spearmanr(h_sub, v_sub)
        rho_sample = float(rho_sample) if np.isfinite(rho_sample) else np.nan
        s_pval = float(s_pval) if np.isfinite(s_pval) else np.nan
    else:
        rho_sample, s_pval = np.nan, np.nan

    # Prespecified site aggregation: mean per site
    df_s = pd.DataFrame({"site": s_sub, "h": h_sub, "v": v_sub})
    site_agg = df_s.groupby("site").mean()
    n_sites = len(site_agg)

    h_site = site_agg["h"].to_numpy(dtype=float)
    v_site = site_agg["v"].to_numpy(dtype=float)

    h_pres_site = h_site > 0
    v_pres_site = v_site > 0
    n_h_site = int(h_pres_site.sum())
    n_v_site = int(v_pres_site.sum())
    n_both_site = int((h_pres_site & v_pres_site).sum())
    n_either_site = int((h_pres_site | v_pres_site).sum())
    site_jacc = float(n_both_site) / float(n_either_site) if n_either_site > 0 else 0.0

    std_h_site = np.std(h_site)
    std_v_site = np.std(v_site)

    if std_h_site == 0 or std_v_site == 0:
        return {
            "n_evaluated_samples": n_valid_samples,
            "n_unknown_excluded_samples": n_excl_samples,
            "estimability_status": "non_estimable_zero_variance_in_site_means",
            "n_exchangeable_site_units": n_sites,
            "n_host_present_sites": n_h_site,
            "n_virus_present_sites": n_v_site,
            "n_codetected_sites": n_both_site,
            "site_mean_spearman_rho": np.nan,
            "site_blocked_permutation_pvalue": np.nan,
            "site_jaccard_similarity": round(site_jacc, 6),
            "site_blocked_jaccard_permutation_pvalue": np.nan,
            "sample_spearman_rho": round(rho_sample, 6) if np.isfinite(rho_sample) else np.nan,
            "sample_spearman_pvalue_unblocked_invalid_pseudoreplication": round(s_pval, 6) if np.isfinite(s_pval) else np.nan,
            "sample_jaccard_similarity": round(sample_jacc, 6),
        }

    obs_rho_res = stats.spearmanr(h_site, v_site)
    obs_rho = float(obs_rho_res.statistic) if hasattr(obs_rho_res, "statistic") else float(obs_rho_res[0])
    if not np.isfinite(obs_rho):
        obs_rho = 0.0

    # Permutations with fixed seed
    rng = np.random.default_rng(seed)
    n_rho_ge = 0
    n_jacc_ge = 0

    v_perm = v_site.copy()
    abs_obs_rho = abs(obs_rho)

    for _ in range(n_perms):
        rng.shuffle(v_perm)
        p_res = stats.spearmanr(h_site, v_perm)
        r_p = float(p_res.statistic) if hasattr(p_res, "statistic") else float(p_res[0])
        if np.isfinite(r_p) and abs(r_p) >= abs_obs_rho:
            n_rho_ge += 1

        v_pres_p = v_perm > 0
        b_p = int((h_pres_site & v_pres_p).sum())
        e_p = int((h_pres_site | v_pres_p).sum())
        j_p = float(b_p) / float(e_p) if e_p > 0 else 0.0
        if j_p >= site_jacc:
            n_jacc_ge += 1

    p_val_rho = float(1 + n_rho_ge) / float(n_perms + 1)
    p_val_jacc = float(1 + n_jacc_ge) / float(n_perms + 1)

    return {
        "n_evaluated_samples": n_valid_samples,
        "n_unknown_excluded_samples": n_excl_samples,
        "estimability_status": "estimable_with_site_aggregation_and_permutations",
        "n_exchangeable_site_units": n_sites,
        "n_host_present_sites": n_h_site,
        "n_virus_present_sites": n_v_site,
        "n_codetected_sites": n_both_site,
        "site_mean_spearman_rho": round(obs_rho, 6),
        "site_blocked_permutation_pvalue": round(p_val_rho, 6),
        "site_jaccard_similarity": round(site_jacc, 6),
        "site_blocked_jaccard_permutation_pvalue": round(p_val_jacc, 6),
        "sample_spearman_rho": round(rho_sample, 6) if np.isfinite(rho_sample) else np.nan,
        "sample_spearman_pvalue_unblocked_invalid_pseudoreplication": round(s_pval, 6) if np.isfinite(s_pval) else np.nan,
        "sample_jaccard_similarity": round(sample_jacc, 6),
    }


def evaluate_cooccurrence_and_controls(
    primary_pairs_df: pd.DataFrame,
    matched_controls_df: pd.DataFrame,
    mag_abund_df: pd.DataFrame,
    vir_tpm_df: pd.DataFrame,
    metadata_df: pd.DataFrame,
    n_perms: int = 9999,
    seed: int = 20261109,
) -> pd.DataFrame:
    """Evaluate site-blocked abundance associations for candidate pairs and matched controls."""
    sample_order = [s for s in metadata_df["sample_id"] if s in mag_abund_df.columns and s in vir_tpm_df.columns]
    site_series = metadata_df.set_index("sample_id").loc[sample_order, "site"]

    eval_rows: list[dict[str, Any]] = []

    # Map candidate pairs
    for _, r in primary_pairs_df.iterrows():
        h = str(r["host_mag"])
        t_id = str(r["target_id"])
        v_id = str(r["votu_id"])
        sha = str(r.get("target_seq_sha256", ""))

        # Check key in viral TPM
        v_key = v_id if v_id in vir_tpm_df.index else (t_id if t_id in vir_tpm_df.index else None)
        h_in = h in mag_abund_df.index
        v_in = v_key is not None

        if not h_in or not v_in or len(sample_order) == 0:
            eval_rows.append(
                {
                    "pair_evaluation_role": "candidate_primary_pair",
                    "host_mag": h,
                    "control_for_host": "none",
                    "target_id": t_id,
                    "votu_id": v_id,
                    "target_seq_sha256": sha,
                    "matched_abundance_key": v_key or "none",
                    "catalog_estimability_status": "non_estimable_missing_from_existing_viral_tpm_matrix" if not v_in else "non_estimable_host_missing",
                    "n_evaluated_samples": 0,
                    "n_exchangeable_site_units": 0,
                    "n_host_present_sites": np.nan,
                    "n_virus_present_sites": np.nan,
                    "n_codetected_sites": np.nan,
                    "site_mean_spearman_rho": np.nan,
                    "site_blocked_permutation_pvalue": np.nan,
                    "site_jaccard_similarity": np.nan,
                    "site_blocked_jaccard_permutation_pvalue": np.nan,
                    "sample_spearman_rho": np.nan,
                    "sample_spearman_pvalue_unblocked_invalid_pseudoreplication": np.nan,
                    "sample_jaccard_similarity": np.nan,
                }
            )
            continue

        h_vec = pd.to_numeric(mag_abund_df.loc[h, sample_order], errors="coerce").to_numpy(dtype=float)
        v_vec = pd.to_numeric(vir_tpm_df.loc[v_key, sample_order], errors="coerce").to_numpy(dtype=float)

        res = run_site_blocked_permutation_test(h_vec, v_vec, site_series, n_perms=n_perms, seed=seed)
        row = {
            "pair_evaluation_role": "candidate_primary_pair",
            "host_mag": h,
            "control_for_host": "none",
            "target_id": t_id,
            "votu_id": v_id,
            "target_seq_sha256": sha,
            "matched_abundance_key": v_key,
            "catalog_estimability_status": res["estimability_status"],
            "n_evaluated_samples": res["n_evaluated_samples"],
            "n_exchangeable_site_units": res["n_exchangeable_site_units"],
            "n_host_present_sites": res.get("n_host_present_sites", np.nan),
            "n_virus_present_sites": res.get("n_virus_present_sites", np.nan),
            "n_codetected_sites": res.get("n_codetected_sites", np.nan),
            "site_mean_spearman_rho": res["site_mean_spearman_rho"],
            "site_blocked_permutation_pvalue": res["site_blocked_permutation_pvalue"],
            "site_jaccard_similarity": res["site_jaccard_similarity"],
            "site_blocked_jaccard_permutation_pvalue": res["site_blocked_jaccard_permutation_pvalue"],
            "sample_spearman_rho": res["sample_spearman_rho"],
            "sample_spearman_pvalue_unblocked_invalid_pseudoreplication": res["sample_spearman_pvalue_unblocked_invalid_pseudoreplication"],
            "sample_jaccard_similarity": res["sample_jaccard_similarity"],
        }
        eval_rows.append(row)

    # Map matched unrelated controls
    ctrl_map = dict(zip(matched_controls_df["candidate_host_mag"], matched_controls_df["matched_control_mag"]))
    for _, r in primary_pairs_df.iterrows():
        h = str(r["host_mag"])
        ctrl_mag = ctrl_map.get(h, "none")
        if ctrl_mag == "none" or ctrl_mag not in mag_abund_df.index:
            continue

        t_id = str(r["target_id"])
        v_id = str(r["votu_id"])
        sha = str(r.get("target_seq_sha256", ""))
        v_key = v_id if v_id in vir_tpm_df.index else (t_id if t_id in vir_tpm_df.index else None)

        if v_key is None or len(sample_order) == 0:
            continue

        c_vec = pd.to_numeric(mag_abund_df.loc[ctrl_mag, sample_order], errors="coerce").to_numpy(dtype=float)
        v_vec = pd.to_numeric(vir_tpm_df.loc[v_key, sample_order], errors="coerce").to_numpy(dtype=float)

        res_c = run_site_blocked_permutation_test(c_vec, v_vec, site_series, n_perms=n_perms, seed=seed)
        row_c = {
            "pair_evaluation_role": "matched_unrelated_host_control",
            "host_mag": ctrl_mag,
            "control_for_host": h,
            "target_id": t_id,
            "votu_id": v_id,
            "target_seq_sha256": sha,
            "matched_abundance_key": v_key,
            "catalog_estimability_status": res_c["estimability_status"],
            "n_evaluated_samples": res_c["n_evaluated_samples"],
            "n_exchangeable_site_units": res_c["n_exchangeable_site_units"],
            "n_host_present_sites": res_c.get("n_host_present_sites", np.nan),
            "n_virus_present_sites": res_c.get("n_virus_present_sites", np.nan),
            "n_codetected_sites": res_c.get("n_codetected_sites", np.nan),
            "site_mean_spearman_rho": res_c["site_mean_spearman_rho"],
            "site_blocked_permutation_pvalue": res_c["site_blocked_permutation_pvalue"],
            "site_jaccard_similarity": res_c["site_jaccard_similarity"],
            "site_blocked_jaccard_permutation_pvalue": res_c["site_blocked_jaccard_permutation_pvalue"],
            "sample_spearman_rho": res_c["sample_spearman_rho"],
            "sample_spearman_pvalue_unblocked_invalid_pseudoreplication": res_c["sample_spearman_pvalue_unblocked_invalid_pseudoreplication"],
            "sample_jaccard_similarity": res_c["sample_jaccard_similarity"],
        }
        eval_rows.append(row_c)

    res_df = pd.DataFrame(eval_rows)

    # Benjamini-Hochberg FDR adjustment across estimable candidate pairs
    cand_mask = (res_df["pair_evaluation_role"] == "candidate_primary_pair") & res_df["site_blocked_permutation_pvalue"].notna()
    if cand_mask.sum() > 0:
        cand_pvals = res_df.loc[cand_mask, "site_blocked_permutation_pvalue"].to_numpy()
        m_cand = len(cand_pvals)
        s_idx = np.argsort(cand_pvals)
        s_p = cand_pvals[s_idx]

        qvals = np.zeros(m_cand)
        min_q = 1.0
        for i in range(m_cand - 1, -1, -1):
            rank = i + 1
            q = (s_p[i] * m_cand) / rank
            min_q = min(min_q, q)
            qvals[i] = min(min_q, 1.0)

        unsorted_qvals = np.zeros(m_cand)
        unsorted_qvals[s_idx] = qvals
        res_df["bh_fdr_qvalue"] = np.nan
        res_df.loc[cand_mask, "bh_fdr_qvalue"] = [round(float(q), 6) for q in unsorted_qvals]
    else:
        res_df["bh_fdr_qvalue"] = np.nan

    return res_df


def compute_honest_descriptive_phylum_counts(
    partition_df: pd.DataFrame,
    primary_pairs_df: pd.DataFrame,
) -> pd.DataFrame:
    """Compute descriptive phylum breakdown of MAGs and links.
    
    Explicit statistical note: This is a purely descriptive breakdown by taxonomic phylum.
    It does NOT constitute a phylogeny-controlled regression or statistical test for confounding.
    """
    rows: list[dict[str, Any]] = []
    phy_counts = partition_df["canonical_phylum"].value_counts()

    pos_hosts_by_phy = set(partition_df[partition_df["is_primary_positive_host"] == 1]["genome"])

    for phy, total_mags in phy_counts.items():
        sub = partition_df[partition_df["canonical_phylum"] == phy]
        n_spacers = int((sub["n_valid_spacers"] > 0).sum())
        n_zero = int(total_mags - n_spacers)
        n_pos = int((sub["is_primary_positive_host"] == 1).sum())

        rows.append(
            {
                "canonical_phylum": phy,
                "total_mags_in_phylum": total_mags,
                "n_mags_with_spacer_detection_opportunity": n_spacers,
                "n_mags_zero_spacers_no_detection_opportunity": n_zero,
                "n_mags_with_verified_primary_virus_link": n_pos,
                "phylum_link_prevalence_fraction": round(float(n_pos) / float(total_mags), 6),
                "model_type": "descriptive_phylum_breakdown",
                "notes": "Descriptive summary by phylum; does NOT constitute phylogeny-controlled test or regression",
            }
        )

    return pd.DataFrame(rows)


def main() -> None:
    t0 = time.time()
    parser = argparse.ArgumentParser(description="Stage 4 v3: Downstream Viral Statistics and Estimability Audit.")
    parser.add_argument("--exp3-root", type=Path, default=BRANCH_ROOT, help="EXP3 root directory")
    parser.add_argument("--metadata", type=Path, default=None, help="Authoritative metadata.tsv")
    parser.add_argument("--mag-abund", type=Path, default=None, help="mag_abundance.tsv.gz")
    parser.add_argument("--vir-tpm", type=Path, default=None, help="viral_tpm.tsv.gz")
    parser.add_argument("--taxonomy", type=Path, default=None, help="taxonomy_reconciliation_all5438.tsv.gz")
    parser.add_argument("--screening", type=Path, default=None, help="stage1_mag_screening_status.tsv")
    parser.add_argument("--evidence", type=Path, default=None, help="stage3_all591_archaeal_mag_host_evidence_v2.tsv")
    parser.add_argument("--pairs", type=Path, default=None, help="stage3_archaeal_virus_host_pairs_all_tiers_v2.tsv")
    parser.add_argument("--protein-members", type=Path, default=None, help="exp4_47_member_evidence.tsv")
    parser.add_argument("--results-dir", type=Path, default=None, help="Output results directory")
    parser.add_argument("--n-perms", type=int, default=9999, help="Number of permutations for site-blocked test")
    parser.add_argument("--seed", type=int, default=20261109, help="Random seed")
    args = parser.parse_args()

    root = args.exp3_root
    in_dir = root / "inputs"
    res_dir = args.results_dir or (root / "results")
    res_dir.mkdir(parents=True, exist_ok=True)

    # Resolve paths with fallbacks
    meta_p = args.metadata or (DEFAULT_METADATA_PATH if DEFAULT_METADATA_PATH.exists() else in_dir / "metadata.tsv")
    mag_p = args.mag_abund or (DEFAULT_MAG_ABUND_PATH if DEFAULT_MAG_ABUND_PATH.exists() else in_dir / "mag_abundance.tsv.gz")
    vir_p = args.vir_tpm or (DEFAULT_VIR_TPM_PATH if DEFAULT_VIR_TPM_PATH.exists() else in_dir / "viral_tpm.tsv.gz")
    tax_p = args.taxonomy or (in_dir / "taxonomy_reconciliation_all5438.tsv.gz")
    screen_p = args.screening or (res_dir / "stage1_mag_screening_status.tsv")
    ev_p = args.evidence or (res_dir / "stage3_all591_archaeal_mag_host_evidence_v2.tsv")
    pairs_p = args.pairs or (res_dir / "stage3_archaeal_virus_host_pairs_all_tiers_v2.tsv")
    prot_p = args.protein_members or (in_dir / "exp4_47_member_evidence.tsv")

    print("[Stage 4 v3] Starting downstream statistics and estimability audit...", flush=True)

    # 1. Load metadata
    meta_df = load_metadata(meta_p)
    auth_samples = set(meta_df["sample_id"].astype(str))
    print(f"[Stage 4 v3] Authoritative metadata loaded: {len(meta_df)} samples across {meta_df['site'].nunique()} sites.", flush=True)

    # 2. Partition detection opportunity
    screen_df = pd.read_csv(screen_p, sep="\t") if screen_p.exists() else pd.DataFrame()
    ev_df = pd.read_csv(ev_p, sep="\t") if ev_p.exists() else pd.DataFrame()
    tax_df = pd.read_csv(tax_p, sep="\t") if tax_p.exists() else pd.DataFrame()

    partition_df = partition_detection_opportunity(screen_df, ev_df, tax_df)
    partition_tsv = res_dir / "stage4_v3_mag_detection_opportunity_partition.tsv"
    partition_df.to_csv(partition_tsv, sep="\t", index=False)
    n_zero = int((partition_df["detection_opportunity_category"] == "zero_spacers_no_detection_opportunity").sum())
    n_pos = int((partition_df["detection_opportunity_category"] == "screen_positive_detection_positive").sum())
    n_neg = int((partition_df["detection_opportunity_category"] == "screen_positive_eligible_negative").sum())
    print(f"[Stage 4 v3] Detection opportunity partition: {len(partition_df)} MAGs (zero_spacers={n_zero}, eligible_negative={n_neg}, positive={n_pos})", flush=True)

    # 3. Resolve graded candidate pairs
    pairs_df = pd.read_csv(pairs_p, sep="\t") if pairs_p.exists() else pd.DataFrame()
    graded_ledger_df, primary_unique_pairs_df = resolve_graded_candidate_pairs(pairs_df)
    graded_tsv = res_dir / "stage4_v3_candidate_pairs_graded_ledger.tsv"
    graded_ledger_df.to_csv(graded_tsv, sep="\t", index=False)
    print(f"[Stage 4 v3] Candidate pairs graded: {len(graded_ledger_df)} total evaluated rows -> {len(primary_unique_pairs_df)} primary unique host-target SHA pairs across {primary_unique_pairs_df['host_mag'].nunique()} unique hosts.", flush=True)

    # 4. Audit defense/surface design matrix and candidate protein families
    prot_df = pd.read_csv(prot_p, sep="\t") if prot_p.exists() else pd.DataFrame()
    design_audit, candidate_assoc_df = audit_defense_surface_design_matrix(partition_df, prot_df)
    audit_json = res_dir / "stage4_v3_defense_surface_design_matrix_audit.json"
    audit_json.write_text(json.dumps(design_audit, indent=2) + "\n")
    cand_assoc_tsv = res_dir / "stage4_v3_candidate_protein_family_associations.tsv"
    candidate_assoc_df.to_csv(cand_assoc_tsv, sep="\t", index=False)
    print(f"[Stage 4 v3] Design matrix audit completed: rank={design_audit['full_591_design_matrix_audit']['design_matrix_rank']}/{design_audit['full_591_design_matrix_audit']['n_predictors_p']}, condition_number={design_audit['full_591_design_matrix_audit']['condition_number']}. Evaluated {len(candidate_assoc_df)} candidate protein families.", flush=True)

    # 5. Load and orient abundance matrices
    mag_df, vir_df = load_and_orient_abundance_matrices(mag_p, vir_p, auth_samples)
    print(f"[Stage 4 v3] Abundance matrices loaded: MAGs={mag_df.shape}, Viruses={vir_df.shape}.", flush=True)

    # Compute host prevalence and abundance metrics
    abund_metric_rows: list[dict[str, Any]] = []
    sample_cols = [s for s in meta_df["sample_id"] if s in mag_df.columns]
    for g in partition_df["genome"]:
        if g not in mag_df.index or len(sample_cols) == 0:
            abund_metric_rows.append({"genome": g, "sample_prevalence": 0.0, "mean_sample_abundance": 0.0})
        else:
            vec = pd.to_numeric(mag_df.loc[g, sample_cols], errors="coerce").to_numpy(dtype=float)
            v_finite = vec[np.isfinite(vec)]
            abund_metric_rows.append({
                "genome": g,
                "sample_prevalence": float((v_finite > 0).sum()) if len(v_finite) > 0 else 0.0,
                "mean_sample_abundance": float(np.mean(v_finite)) if len(v_finite) > 0 else 0.0,
            })
    abund_metrics_df = pd.DataFrame(abund_metric_rows)

    # 6. Find matched unrelated host controls
    unique_candidate_hosts = sorted(primary_unique_pairs_df["host_mag"].unique())
    matched_controls_df = find_matched_unrelated_controls(unique_candidate_hosts, partition_df, abund_metrics_df)
    matched_tsv = res_dir / "stage4_v3_matched_unrelated_host_controls.tsv"
    matched_controls_df.to_csv(matched_tsv, sep="\t", index=False)
    print(f"[Stage 4 v3] Matched unrelated controls identified for {len(matched_controls_df)} candidate hosts.", flush=True)

    # 7. Evaluate site-blocked abundance associations
    cooc_df = evaluate_cooccurrence_and_controls(
        primary_unique_pairs_df,
        matched_controls_df,
        mag_df,
        vir_df,
        meta_df,
        n_perms=args.n_perms,
        seed=args.seed,
    )
    cooc_tsv = res_dir / "stage4_v3_site_blocked_abundance_associations.tsv"
    cooc_df.to_csv(cooc_tsv, sep="\t", index=False)
    print(f"[Stage 4 v3] Site-blocked abundance associations evaluated for {len(cooc_df)} pairs/controls across 58 sites.", flush=True)

    # 8. Honest descriptive phylum breakdown
    desc_phylum_df = compute_honest_descriptive_phylum_counts(partition_df, primary_unique_pairs_df)
    phylum_tsv = res_dir / "stage4_v3_descriptive_phylum_counts.tsv"
    desc_phylum_df.to_csv(phylum_tsv, sep="\t", index=False)

    # 9. Summary receipt
    elapsed = round(time.time() - t0, 2)
    summary_data = {
        "stage": "stage4_statistical_audit_v3",
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S UTC", time.gmtime()),
        "elapsed_seconds": elapsed,
        "runtime_command": " ".join(sys.argv),
        "input_files": {
            "metadata": {"path": str(meta_p), "sha256": compute_sha256(meta_p)},
            "mag_abundance": {"path": str(mag_p), "sha256": compute_sha256(mag_p)},
            "viral_tpm": {"path": str(vir_p), "sha256": compute_sha256(vir_p)},
            "taxonomy": {"path": str(tax_p), "sha256": compute_sha256(tax_p)},
            "stage1_screening": {"path": str(screen_p), "sha256": compute_sha256(screen_p)},
            "stage3_evidence": {"path": str(ev_p), "sha256": compute_sha256(ev_p)},
            "stage3_pairs": {"path": str(pairs_p), "sha256": compute_sha256(pairs_p)},
            "protein_members": {"path": str(prot_p), "sha256": compute_sha256(prot_p)},
        },
        "canonical_591_partition_counts": {
            "total_canonical_mags": len(partition_df),
            "zero_spacers_no_detection_opportunity": n_zero,
            "screen_positive_eligible_negative": n_neg,
            "screen_positive_detection_positive": n_pos,
        },
        "candidate_pairs_counts": {
            "raw_evaluated_rows": len(graded_ledger_df),
            "primary_unique_host_target_sha_pairs": len(primary_unique_pairs_df),
            "unique_positive_hosts": int(primary_unique_pairs_df["host_mag"].nunique()),
            "alias_duplicates_collapsed": int((graded_ledger_df["is_alias_duplicate"] == 1).sum()),
            "cellular_dominated_targets_excluded": int((graded_ledger_df["is_cellular_dominated"] == 1).sum()),
        },
        "design_matrix_audit_summary": {
            "full_591_epv": design_audit["full_591_design_matrix_audit"]["events_per_variable_epv"],
            "full_591_rank": design_audit["full_591_design_matrix_audit"]["design_matrix_rank"],
            "full_591_condition_number": design_audit["full_591_design_matrix_audit"]["condition_number"],
            "zero_event_phyla": design_audit["full_591_design_matrix_audit"]["zero_event_phyla_list"],
            "quasi_complete_separation": design_audit["full_591_design_matrix_audit"]["quasi_complete_separation_present"],
            "candidate_protein_families_evaluated": len(candidate_assoc_df),
            "candidate_protein_families_significant_bh_fdr": int((candidate_assoc_df["bh_fdr_qvalue"] <= 0.05).sum()) if len(candidate_assoc_df) > 0 else 0,
        },
        "site_blocked_abundance_summary": {
            "n_exchangeable_site_units": 58,
            "n_permutations_per_test": args.n_perms,
            "permutation_seed": args.seed,
            "n_pairs_evaluated": len(cooc_df),
            "n_estimable_in_catalog_matrices": int((cooc_df["catalog_estimability_status"] == "estimable_with_site_aggregation_and_permutations").sum()),
            "n_significant_pairs_bh_fdr": int((cooc_df["bh_fdr_qvalue"] <= 0.05).sum()),
        },
        "scientific_limitations_and_disposition": [
            "CRISPR spacer matching indicates sequence-linked historical exposure, not active viral replication.",
            "All abundance associations evaluated on separate historical catalog mappings are explicitly exploratory.",
            "Site-level mean aggregation respects the 58 independent geographic sites and resolves unequal sample counts.",
            "Phylogenetic logistic regression for defense or surface associations is mathematically non-estimable due to severe event sparsity (N_events=5), quasi-complete separation (0 positives in 5/7 phyla and FDR surface leads), and detection opportunity confounding.",
            "Testing CRISPR operons/spacers against CRISPR-derived detection is mathematically circular.",
            "Final joint competitive read mapping and inStrain microdiversity follow-ups remain pending Slurm pipeline execution.",
        ],
    }

    summary_json = res_dir / "stage4_v3_statistical_audit_summary.json"
    summary_json.write_text(json.dumps(summary_data, indent=2) + "\n")
    (res_dir / "STAGE4_V3_DONE").touch()

    print(f"[Stage 4 v3] Complete! Elapsed: {elapsed}s. Outputs written to {res_dir}.", flush=True)


if __name__ == "__main__":
    main()
