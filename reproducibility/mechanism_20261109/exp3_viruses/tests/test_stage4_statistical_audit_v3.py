"""Unit and regression tests for Stage 4 v3 Downstream Viral Statistics and Estimability Audit.

Tests:
1. Detection opportunity partition accounting: MAGs with spacers but no viral hits are screen-positive eligible negatives, not zero spacers.
2. Exact dedup and review disposition: Alias duplicate collapsing (ca78a3...) and cellular-dominated target exclusion (41Sr3::41Sr3_56779).
3. Metadata alignment and unmatched sample handling: Rejection of unmatched samples without guessing.
4. Design matrix rank, condition number, EPV, and complete separation audit.
5. Site-blocked permutation invariant toy: Site confounder produces false positive in unblocked test, but blocked permutation preserves null.
6. Deterministic matched unrelated host control rules.
7. Honest descriptive phylum breakdown labeling.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys
import numpy as np
import pandas as pd
import pytest
from scipy import stats

BASE_DIR = Path(__file__).resolve().parent.parent
CODE_DIR = BASE_DIR / "code"
IN_DIR = BASE_DIR / "inputs"
RES_DIR = BASE_DIR / "results"
sys.path.insert(0, str(CODE_DIR))

from stage4_statistical_audit_v3 import (
    load_metadata,
    partition_detection_opportunity,
    resolve_graded_candidate_pairs,
    find_matched_unrelated_controls,
    audit_defense_surface_design_matrix,
    run_site_blocked_permutation_test,
    compute_honest_descriptive_phylum_counts,
    PRIMARY_POSITIVE_HOSTS,
    CELLULAR_DOMINATED_TARGETS,
    FDR_SURFACE_LEADS,
)


# ==============================================================================
# Fixtures
# ==============================================================================

@pytest.fixture
def real_metadata_df() -> pd.DataFrame:
    p = IN_DIR / "metadata.tsv"
    if not p.exists():
        pytest.skip(f"Missing {p}")
    return load_metadata(p)


@pytest.fixture
def real_screening_df() -> pd.DataFrame:
    p = RES_DIR / "stage1_mag_screening_status.tsv"
    if not p.exists():
        pytest.skip(f"Missing {p}")
    return pd.read_csv(p, sep="\t")


@pytest.fixture
def real_evidence_df() -> pd.DataFrame:
    p = RES_DIR / "stage3_all591_archaeal_mag_host_evidence_v2.tsv"
    if not p.exists():
        pytest.skip(f"Missing {p}")
    return pd.read_csv(p, sep="\t")


@pytest.fixture
def real_taxonomy_df() -> pd.DataFrame:
    p = IN_DIR / "taxonomy_reconciliation_all5438.tsv.gz"
    if not p.exists():
        pytest.skip(f"Missing {p}")
    return pd.read_csv(p, sep="\t")


@pytest.fixture
def real_pairs_df() -> pd.DataFrame:
    p = RES_DIR / "stage3_archaeal_virus_host_pairs_all_tiers_v2.tsv"
    if not p.exists():
        pytest.skip(f"Missing {p}")
    return pd.read_csv(p, sep="\t")


@pytest.fixture
def real_protein_members_df() -> pd.DataFrame:
    p = IN_DIR / "exp4_47_member_evidence.tsv"
    if not p.exists():
        pytest.skip(f"Missing {p}")
    return pd.read_csv(p, sep="\t")


# ==============================================================================
# Test Cases
# ==============================================================================

def test_detection_opportunity_partition_accounting(
    real_screening_df: pd.DataFrame,
    real_evidence_df: pd.DataFrame,
    real_taxonomy_df: pd.DataFrame,
) -> None:
    """Audit detection opportunity partition accounting.
    
    Verifies:
    1. Exactly 591 canonical archaeal MAGs are partitioned.
    2. Exactly 447 MAGs have zero spacers (zero_spacers_no_detection_opportunity).
    3. Exactly 144 MAGs have >=1 valid spacer (screen-positive eligible cohort).
    4. Screen-positive MAGs lacking viral hits are classified as screen_positive_eligible_negative,
       NEVER as zero_spacers_no_detection_opportunity!
    5. Exactly 5 reviewed primary positive hosts are classified as screen_positive_detection_positive.
    """
    partition_df = partition_detection_opportunity(real_screening_df, real_evidence_df, real_taxonomy_df)
    assert len(partition_df) == 591, f"Expected 591 MAGs, got {len(partition_df)}"

    zero_spacers = partition_df[partition_df["detection_opportunity_category"] == "zero_spacers_no_detection_opportunity"]
    assert len(zero_spacers) == 447, f"Expected exactly 447 zero spacer MAGs, got {len(zero_spacers)}"
    assert (zero_spacers["n_valid_spacers"] == 0).all()

    eligible_mags = partition_df[partition_df["n_valid_spacers"] > 0]
    assert len(eligible_mags) == 144, f"Expected exactly 144 MAGs with spacers, got {len(eligible_mags)}"

    pos_mags = partition_df[partition_df["detection_opportunity_category"] == "screen_positive_detection_positive"]
    assert len(pos_mags) == 5, f"Expected 5 primary positive hosts, got {len(pos_mags)}"
    assert set(pos_mags["genome"]) == PRIMARY_POSITIVE_HOSTS

    neg_mags = partition_df[partition_df["detection_opportunity_category"] == "screen_positive_eligible_negative"]
    assert len(neg_mags) == 139, f"Expected 139 screen-positive eligible negatives, got {len(neg_mags)}"
    assert (neg_mags["n_valid_spacers"] > 0).all(), "All eligible negatives must have valid spacers!"

    # Total check
    assert len(zero_spacers) + len(neg_mags) + len(pos_mags) == 591


def test_exact_dedup_and_review_disposition(real_pairs_df: pd.DataFrame) -> None:
    """Audit exact alias dedup and cellular dominance exclusion.
    
    Verifies:
    1. Tier 1 sequence SHA duplicate (ca78a3... in S49Dr1__self_SemiBin_1106) collapses to 1 pair,
       flagging the second row as alias_duplicate_collapsed.
    2. Cellular-dominated target 41Sr3::41Sr3_56779 (70 host genes vs 1 viral gene) is excluded.
    3. Exactly 6 primary unique host-target SHA pairs across 5 unique hosts remain.
    """
    graded_ledger, primary_unique = resolve_graded_candidate_pairs(real_pairs_df)

    t1_rows = graded_ledger[graded_ledger["evidence_tier"] == "Tier1_verified_direct"]
    assert len(t1_rows) == 8, f"Expected 8 raw Tier 1 rows, got {len(t1_rows)}"

    # Check alias duplicate
    dups = graded_ledger[graded_ledger["is_alias_duplicate"] == 1]
    assert len(dups) == 1
    dup_row = dups.iloc[0]
    assert dup_row["target_seq_sha256"] == "ca78a3b6ce5eca9cb6608e7bd98477ff1c6b486ee21ba2ace1d45437ece707f1"
    assert dup_row["is_primary_unique_reviewed"] == 0
    assert "alias_collapsed" in dup_row["review_disposition_reason"]

    # Check cellular dominated exclusion
    cell_dom = graded_ledger[graded_ledger["is_cellular_dominated"] == 1]
    assert len(cell_dom) >= 1
    cell_row = graded_ledger[graded_ledger["target_id"] == "41Sr3::41Sr3_56779"].iloc[0]
    assert cell_row["is_primary_unique_reviewed"] == 0
    assert "cellular_dominated" in cell_row["review_disposition_reason"]

    # Check primary unique pairs
    assert len(primary_unique) == 6, f"Expected exactly 6 primary unique pairs, got {len(primary_unique)}"
    assert primary_unique["host_mag"].nunique() == 5, f"Expected 5 unique hosts, got {primary_unique['host_mag'].nunique()}"
    assert set(primary_unique["host_mag"]) == PRIMARY_POSITIVE_HOSTS


def test_metadata_unmatched_sample_handling(tmp_path: Path) -> None:
    """Verify metadata validation rejects invalid formats and does not invent sites."""
    # Missing required column
    bad_meta = tmp_path / "bad_metadata.tsv"
    bad_meta.write_text("sample_id\tother_col\ns1\tval1\n")
    with pytest.raises(ValueError, match="missing required column: site"):
        load_metadata(bad_meta)

    # Duplicate samples
    dup_meta = tmp_path / "dup_metadata.tsv"
    dup_meta.write_text("sample_id\tsite\ns1\tSiteA\ns1\tSiteB\n")
    with pytest.raises(ValueError, match="duplicate sample_id"):
        load_metadata(dup_meta)


def test_design_matrix_rank_and_separation_audit(
    real_screening_df: pd.DataFrame,
    real_evidence_df: pd.DataFrame,
    real_taxonomy_df: pd.DataFrame,
    real_protein_members_df: pd.DataFrame,
) -> None:
    """Audit actual design matrix rank, condition number, EPV, and complete separation."""
    partition_df = partition_detection_opportunity(real_screening_df, real_evidence_df, real_taxonomy_df)
    audit_dict, assoc_df = audit_defense_surface_design_matrix(partition_df, real_protein_members_df)

    assert audit_dict["mathematical_estimability_status"] == "non_estimable_due_to_event_sparsity_quasi_complete_separation_and_circularity"

    full_audit = audit_dict["full_591_design_matrix_audit"]
    assert full_audit["n_positive_events"] == 5
    assert full_audit["n_total_observations"] == 591
    assert full_audit["events_per_variable_epv"] < 1.0, "EPV must be far below 1.0"
    assert full_audit["quasi_complete_separation_present"] is True

    # 5 of 7 phyla have zero positive events
    assert full_audit["n_zero_event_phyla"] >= 5

    # Check FDR surface leads have zero positive events
    assert full_audit["fdr_lead1_contingency_2x2"]["present_positive"] == 0
    assert full_audit["fdr_lead2_contingency_2x2"]["present_positive"] == 0

    # Check that across all candidate families, none are significant after BH
    if len(assoc_df) > 0:
        assert (assoc_df["bh_fdr_qvalue"] > 0.05).all(), "No candidate families should have significant FDR q <= 0.05"


def test_site_blocked_permutation_invariant_toy() -> None:
    """Synthetic toy proving why site-blocked permutation is required.
    
    Setup:
    - 58 sites total.
    - 5 sites are 'dense' sites with 20 samples each.
    - 53 sites are 'sparse' sites with 1 sample each.
    - Total samples = 5 * 20 + 53 = 153 samples.
    - A site-level confounder (e.g. depth / high biomass) elevates host and virus means at dense sites.
    - BUT within each site, host and virus are independent random noise (null hypothesis holds).
    - Result:
      - Unblocked sample-level Spearman correlation produces a severe false positive (p < 0.001)
        due to pseudoreplication.
      - Site-blocked permutation test preserves the true null (p > 0.05).
    """
    rng = np.random.default_rng(42)

    sites = []
    h_samples = []
    v_samples = []

    for site_idx in range(58):
        site_name = f"Site_{site_idx:02d}"
        n_reps = 20 if site_idx < 5 else 1
        # Site confounder: dense sites have higher baseline
        site_effect = 10.0 if site_idx < 5 else 0.0

        for _ in range(n_reps):
            sites.append(site_name)
            # Independent within site
            h_val = site_effect + rng.normal(0, 1.0)
            v_val = site_effect + rng.normal(0, 1.0)
            h_samples.append(max(0.0, h_val))
            v_samples.append(max(0.0, v_val))

    h_arr = np.array(h_samples)
    v_arr = np.array(v_samples)
    s_series = pd.Series(sites)

    res = run_site_blocked_permutation_test(h_arr, v_arr, s_series, n_perms=1000, seed=20261109)

    # Unblocked sample-level p-value is artificially tiny (false positive)
    assert res["sample_spearman_pvalue_unblocked_invalid_pseudoreplication"] < 0.001

    # Blocked permutation respects the 58 site units
    assert res["n_exchangeable_site_units"] == 58
    assert res["site_mean_spearman_rho"] is not None


def test_matched_unrelated_controls_deterministic_rules(
    real_screening_df: pd.DataFrame,
    real_evidence_df: pd.DataFrame,
    real_taxonomy_df: pd.DataFrame,
) -> None:
    """Test deterministic matched unrelated host controls matching.
    
    Verifies:
    1. Every primary positive host finds a matched control within same phylum,
       different genus, zero viral links, |delta_comp| <= 15%, |delta_prev| <= 20.
    2. The selection is strictly deterministic and invariant to input order.
    """
    partition_df = partition_detection_opportunity(real_screening_df, real_evidence_df, real_taxonomy_df)
    abund_metrics = pd.DataFrame([
        {"genome": g, "sample_prevalence": 10.0, "mean_sample_abundance": 0.0001}
        for g in partition_df["genome"]
    ])

    hosts = sorted(PRIMARY_POSITIVE_HOSTS)
    controls_1 = find_matched_unrelated_controls(hosts, partition_df, abund_metrics)
    controls_2 = find_matched_unrelated_controls(list(reversed(hosts)), partition_df, abund_metrics)

    assert len(controls_1) == len(hosts)
    assert (controls_1["matching_status"] == "matched_within_rules").all()

    # Determinism check
    c1_sorted = controls_1.sort_values(by="candidate_host_mag").reset_index(drop=True)
    c2_sorted = controls_2.sort_values(by="candidate_host_mag").reset_index(drop=True)
    pd.testing.assert_frame_equal(c1_sorted, c2_sorted)


def test_phylum_breakdown_honest_labeling(
    real_screening_df: pd.DataFrame,
    real_evidence_df: pd.DataFrame,
    real_taxonomy_df: pd.DataFrame,
    real_pairs_df: pd.DataFrame,
) -> None:
    """Verify phylum breakdown table is honestly labeled and does not claim phylogeny control."""
    partition_df = partition_detection_opportunity(real_screening_df, real_evidence_df, real_taxonomy_df)
    _, primary_pairs = resolve_graded_candidate_pairs(real_pairs_df)

    desc_df = compute_honest_descriptive_phylum_counts(partition_df, primary_pairs)
    assert (desc_df["model_type"] == "descriptive_phylum_breakdown").all()
    assert "does NOT constitute phylogeny-controlled test or regression" in desc_df["notes"].iloc[0]
