"""Package completeness, receipt provenance, and independent 47-site replay tests for stats_v4."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(BASE / "code"))
sys.path.insert(0, str(BASE))

from stage4_statistical_audit_v4 import (  # noqa: E402
    bh_adjust,
    build_family_assay_maps,
    spearman_permutation_result,
    stratified_permutation_indices,
    validate_assay_provenance,
)
from verify_stats_v4 import verify_all  # noqa: E402


def test_packaged_stage1_and_stage1b_receipts_and_gate_false_retention() -> None:
    full_path = BASE / "inputs/task_a_stage1b_all_candidate_member_loci_verified.tsv.gz"
    stage1 = BASE / "inputs/stage1_summary.json"
    stage1b = BASE / "inputs/stage1b_summary.json"
    full = pd.read_csv(full_path, sep="\t")
    provenance = validate_assay_provenance(full, stage1, stage1b)
    assert provenance["n_screened_genomes"] == 237
    assert provenance["n_original_proteins_in_stage1_scope"] == 777177
    assert provenance["n_candidate_member_locus_rows_complete_membership"] == 41494
    assert provenance["n_exact_family_ids"] == 8776
    assert provenance["n_final_promotion_gate_false_rows_retained_as_membership_evidence"] == 3184

    thermo237 = pd.read_csv(BASE / "inputs/taxonomy_reconciliation_thermo237.tsv", sep="\t")
    assert set(full["genome"].astype(str)) == set(thermo237["genome"].astype(str))

    failed_members = full.loc[~full["final_promotion_gate_pass"].astype(bool), ["family_id", "genome"]]
    assert len(failed_members) == 3184
    members = pd.read_csv(BASE / "inputs/exp4_47_member_evidence.tsv", sep="\t")
    _, member_maps, _, _ = build_family_assay_maps(members, full, provenance)
    frozen_failed = failed_members.loc[failed_members["family_id"].isin(member_maps)]
    assert len(frozen_failed) > 0
    assert all(row.genome in member_maps[row.family_id] for row in frozen_failed.itertuples(index=False))


def test_independent_audit_selected47_raw_observations_recompute_all_pairs() -> None:
    raw47 = pd.read_csv(BASE / "independent_audit/stage4_v4_independent_audit_selected47_raw_site_observations.tsv", sep="\t")
    wide47 = pd.read_csv(BASE / "independent_audit/stage4_v4_independent_audit_selected47_wide_matrix.tsv", sep="\t")
    prod_abund = pd.read_csv(BASE / "results/stage4_v4_environment_matched_abundance.tsv", sep="\t")
    prod_diff = pd.read_csv(BASE / "results/stage4_v4_candidate_control_correlation_differences.tsv", sep="\t")
    receipt = json.loads((BASE / "independent_audit/stage4_v4_independent_audit_extraction_receipt.json").read_text())

    assert len(raw47) == 282  # 6 candidate pairs x 47 sites
    assert len(wide47) == 47
    assert raw47["campaign_trip"].eq("trip5").all()
    assert raw47["compartment"].eq("PlantRoot").all()
    assert raw47["compartment_mag"].eq("RootAdjacent").all()
    assert raw47["mag_abundance_source_sha256"].eq("4915c95b5a72dbef5490df1db8b9e88f4cf0c9390de081f82386addf2a1a635c").all()
    assert raw47["viral_tpm_source_sha256"].eq("eaec54b7506ca9e3fc41376e566193f94d3bd2814c3814c6b2546c798276685b").all()

    field_counts = wide47["dune_field"].value_counts().to_dict()
    assert field_counts == {"uruq": 28, "seif": 11, "megabarchans": 6, "sigmoidal and hooked": 2}

    # Recompute Spearman and 9,999 within-field permutations directly from the 47-site raw observations
    perms = stratified_permutation_indices(wide47.sort_values("site")["dune_field"].to_numpy(dtype=str), 9999, 20261109)
    primary_cand_pvals = []
    for _, group in raw47.groupby("pair_order", sort=True):
        g = group.sort_values("site").reset_index(drop=True)
        host = g["candidate_host_mag"].iloc[0]
        ctrl = g["matched_alternative_mag"].iloc[0]
        target = g["target_id"].iloc[0]
        is_prim = bool(g["is_primary_detected_link_pair"].iloc[0] == 1)
        has_vir = g["viral_catalog_status"].iloc[0] == "observed_in_viral_tpm_catalog"

        if has_vir:
            cand_res = spearman_permutation_result(
                g["candidate_host_site_mean_abundance"].to_numpy(dtype=float),
                g["viral_votu_site_mean_tpm"].to_numpy(dtype=float),
                perms,
            )
            ctrl_res = spearman_permutation_result(
                g["matched_alternative_site_mean_abundance"].to_numpy(dtype=float),
                g["viral_votu_site_mean_tpm"].to_numpy(dtype=float),
                perms,
            )
        else:
            cand_res = {"observed_rho": np.nan, "pvalue": np.nan}
            ctrl_res = {"observed_rho": np.nan, "pvalue": np.nan}

        if is_prim:
            primary_cand_pvals.append(cand_res["pvalue"])
            prod_row = prod_abund.loc[
                prod_abund["analysis_set"].eq("primary")
                & prod_abund["pair_evaluation_role"].eq("candidate")
                & prod_abund["host_mag"].eq(host)
                & prod_abund["target_id"].eq(target)
            ].iloc[0]
            if np.isfinite(cand_res["observed_rho"]):
                assert np.isclose(cand_res["observed_rho"], prod_row["site_mean_spearman_rho"], atol=1e-12)
                assert np.isclose(cand_res["pvalue"], prod_row["field_stratified_permutation_pvalue"], atol=1e-12)
            else:
                assert pd.isna(prod_row["site_mean_spearman_rho"])
                assert pd.isna(prod_row["field_stratified_permutation_pvalue"])

            if has_vir:
                ctrl_row = prod_abund.loc[
                    prod_abund["analysis_set"].eq("primary")
                    & prod_abund["pair_evaluation_role"].eq("matched_alternative_MAG")
                    & prod_abund["control_for_host"].eq(host)
                    & prod_abund["target_id"].eq(target)
                ].iloc[0]
                if np.isfinite(ctrl_res["observed_rho"]):
                    assert np.isclose(ctrl_res["observed_rho"], ctrl_row["site_mean_spearman_rho"], atol=1e-12)
                    assert np.isclose(ctrl_res["pvalue"], ctrl_row["field_stratified_permutation_pvalue"], atol=1e-12)
                else:
                    assert pd.isna(ctrl_row["site_mean_spearman_rho"])

    qvals = bh_adjust(np.array(primary_cand_pvals, dtype=float), family_size=5)
    finite_q = qvals[np.isfinite(qvals)]
    assert np.allclose(np.sort(finite_q), np.array([0.4155, 0.423]), atol=1e-12)

    # Verify descriptive difference for 41PRr2__single_SemiBin_97 -> 9Sr2::9Sr2_284237
    diff_row = prod_diff.loc[
        prod_diff["analysis_set"].eq("primary") & prod_diff["target_id"].eq("9Sr2::9Sr2_284237")
    ].iloc[0]
    assert np.isclose(diff_row["observed_candidate_rho"], 0.3347580880891416, atol=1e-12)
    assert np.isclose(diff_row["observed_alternative_rho"], 0.48310328782590084, atol=1e-12)
    assert np.isclose(diff_row["observed_rho_difference_candidate_minus_alternative"], -0.14834519973675925, atol=1e-12)
    assert receipt["exact_replay_verification"]["environment_matched_abundance_exact_match"] is True


def test_full_8_gate_verifier() -> None:
    result = verify_all(BASE)
    assert result["status"] == "PASS" and result["gates_passed"] == 8
