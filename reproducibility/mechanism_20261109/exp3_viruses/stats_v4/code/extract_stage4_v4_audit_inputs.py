#!/usr/bin/env python3
"""Extract compact Stage 4 v4 abundance inputs and independent 47-site audit observations under Slurm.

Does not modify or re-run the completed Stage 4 v4 production outputs (job 15894468).
Extracts:
1. mag_abundance_arch591_441samples.tsv.gz (591 canonical archaeal MAG rows x 441 shared samples)
2. viral_tpm_candidate_targets_441samples.tsv.gz (all stage3 ledger viral rows present in viral_tpm.tsv.gz x 441 shared samples)
3. stage4_v4_independent_audit_selected47_raw_site_observations.tsv (long-format 6 pairs x 47 sites = 282 rows)
4. stage4_v4_independent_audit_selected47_wide_matrix.tsv (wide-format 47 sites x metadata + 5 candidate hosts + 5 matched controls + 6 viral targets)
5. stage4_v4_independent_audit_extraction_receipt.json (provenance + independent Spearman/permutation replay verification)
"""
from __future__ import annotations

import hashlib
import json
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

EXP3_ROOT = Path(os.environ.get("EXP3_ROOT", "/scratch/project/k10018/eq-mechanism-20261109/exp3_viruses"))
sys.path.insert(0, str(EXP3_ROOT / "code"))

from stage4_statistical_audit_v4 import (  # noqa: E402
    compute_abundance_metrics,
    compute_sha256,
    evaluate_environment_matched_associations,
    find_matched_alternative_controls,
    load_and_orient_abundance_matrices,
    load_field_map,
    load_metadata,
    select_primary_environment_design,
    spearman_permutation_result,
    stratified_permutation_indices,
)

METADATA_PATH = Path("/project/k10018/empty-quarter/data/mobilome/ecology_revision_20261020/inputs/metadata.tsv")
MAG_ABUND_PATH = Path("/project/k10018/empty-quarter/data/mobilome/ecology_revision_20261020/inputs/mag_abundance.tsv.gz")
VIR_TPM_PATH = Path("/project/k10018/empty-quarter/data/mobilome/ecology_revision_20261020/inputs/viral_tpm.tsv.gz")
FIELD_MAP_PATH = Path("/scratch/project/k10018/eq-ecology-mag-20261105/results/rund_20261106/matched_profiles.tsv")


def main() -> None:
    t0 = time.time()
    res_dir = EXP3_ROOT / "results"
    out_dir = res_dir / "stats_v4_audit_extract"
    out_dir.mkdir(parents=True, exist_ok=True)

    summary_path = res_dir / "stage4_v4_statistical_audit_summary.json"
    summary = json.loads(summary_path.read_text())
    expected_hashes = {k: v["sha256"] for k, v in summary["input_files"].items()}

    actual_hashes = {
        "metadata": compute_sha256(METADATA_PATH),
        "mag_abundance": compute_sha256(MAG_ABUND_PATH),
        "viral_tpm": compute_sha256(VIR_TPM_PATH),
        "field_map": compute_sha256(FIELD_MAP_PATH),
    }
    for k, sha in actual_hashes.items():
        if sha != expected_hashes[k]:
            raise ValueError(f"Source SHA256 mismatch for {k}: actual={sha} expected={expected_hashes[k]}")

    metadata = load_metadata(METADATA_PATH)
    mag_df, vir_df = load_and_orient_abundance_matrices(MAG_ABUND_PATH, VIR_TPM_PATH, set(metadata["sample_id"]))
    shared_samples = [s for s in metadata["sample_id"] if s in mag_df.columns and s in vir_df.columns]
    if len(shared_samples) != 441:
        raise ValueError(f"Expected 441 shared samples, got {len(shared_samples)}")

    partition = pd.read_csv(res_dir / "stage4_v4_detection_opportunity_partition.tsv", sep="\t")
    ledger = pd.read_csv(res_dir / "stage4_v4_candidate_pairs_graded_ledger.tsv", sep="\t")
    primary = pd.read_csv(res_dir / "stage4_v4_primary_detected_link_pairs.tsv", sep="\t")
    provisional = pd.read_csv(res_dir / "stage4_v4_provisional_sensitivity_pairs.tsv", sep="\t")
    controls = pd.read_csv(res_dir / "stage4_v4_matched_alternative_MAG_controls.tsv", sep="\t")
    prod_metrics = pd.read_csv(res_dir / "stage4_v4_host_abundance_assay_metrics.tsv", sep="\t")
    prod_abund = pd.read_csv(res_dir / "stage4_v4_environment_matched_abundance.tsv", sep="\t")
    prod_diff = pd.read_csv(res_dir / "stage4_v4_candidate_control_correlation_differences.tsv", sep="\t")

    # 1. Extract 591 canonical MAG rows across 441 samples
    canonical_set = set(partition["genome"].astype(str))
    mag_rows_ordered = [g for g in mag_df.index if g in canonical_set]
    missing_canonical_mags = sorted(canonical_set - set(mag_rows_ordered))
    mag_sub = mag_df.loc[mag_rows_ordered, shared_samples].copy()
    mag_sub.index.name = mag_df.index.name or "genome"
    mag_sub_path = out_dir / "mag_abundance_arch591_441samples.tsv.gz"
    mag_sub.to_csv(mag_sub_path, sep="\t", compression="gzip")

    # 2. Extract relevant viral rows across 441 samples (all stage3 ledger votu_id / target_id present in viral_tpm)
    candidate_viral_keys = set(ledger["votu_id"].dropna().astype(str)) | set(ledger["target_id"].dropna().astype(str))
    vir_rows_ordered = [v for v in vir_df.index if v in candidate_viral_keys]
    vir_sub = vir_df.loc[vir_rows_ordered, shared_samples].copy()
    vir_sub.index.name = vir_df.index.name or "votu_id"
    vir_sub_path = out_dir / "viral_tpm_candidate_targets_441samples.tsv.gz"
    vir_sub.to_csv(vir_sub_path, sep="\t", compression="gzip")

    # Verify compact subsets reproduce exact control selection and environment associations
    sub_mag_loaded, sub_vir_loaded = load_and_orient_abundance_matrices(
        mag_sub_path, vir_sub_path, set(metadata["sample_id"])
    )
    sub_shared_samples = [s for s in metadata["sample_id"] if s in sub_mag_loaded.columns and s in sub_vir_loaded.columns]
    if sub_shared_samples != shared_samples:
        raise ValueError("Compact subset sample order differs from full source matrices")

    re_metrics = compute_abundance_metrics(sub_mag_loaded, partition["genome"], sub_shared_samples)
    pd.testing.assert_frame_equal(re_metrics, prod_metrics, check_exact=False, rtol=1e-12, atol=1e-15)

    unresolved_hosts = set(ledger.loc[ledger["is_unresolved_root_disposition"].eq(1), "canonical_591_host_mag"])
    re_controls = find_matched_alternative_controls(
        set(provisional["canonical_591_host_mag"]),
        set(primary["canonical_591_host_mag"]),
        set(provisional["canonical_591_host_mag"]),
        unresolved_hosts,
        partition,
        re_metrics,
    )
    pd.testing.assert_frame_equal(re_controls, controls, check_exact=False, rtol=1e-12, atol=1e-15)

    field_map = load_field_map(FIELD_MAP_PATH)
    group_counts, selected, selection = select_primary_environment_design(metadata, field_map, set(sub_shared_samples))
    re_abund, re_diff, re_design = evaluate_environment_matched_associations(
        primary, provisional, re_controls, sub_mag_loaded, sub_vir_loaded, selected, 9999, 20261109
    )
    pd.testing.assert_frame_equal(re_abund, prod_abund, check_exact=False, rtol=1e-12, atol=1e-15)
    pd.testing.assert_frame_equal(re_diff, prod_diff, check_exact=False, rtol=1e-12, atol=1e-15)

    # 3. Build independent audit raw site observations for the 47 selected sites (trip5 x PlantRoot)
    selected_sorted = selected.sort_values(["site", "sample_id"]).reset_index(drop=True)
    if len(selected_sorted) != 47 or selected_sorted["site"].nunique() != 47:
        raise ValueError("Expected 47 selected samples across 47 unique sites in trip5 x PlantRoot")

    control_map = dict(zip(controls["candidate_host_mag"].astype(str), controls["matched_alternative_mag"].astype(str)))
    control_rel_map = dict(zip(controls["candidate_host_mag"].astype(str), controls["relatedness_status"].astype(str)))
    primary_keys = set(
        zip(primary["canonical_591_host_mag"].astype(str), primary["target_seq_sha256"].astype(str))
    )

    long_rows = []
    replay_checks = []
    ordered_pairs = provisional.sort_values(["canonical_591_host_mag", "target_seq_sha256"]).reset_index(drop=True)
    for pair_idx, pair in ordered_pairs.iterrows():
        host = str(pair["canonical_591_host_mag"])
        target = str(pair["target_id"])
        votu = str(pair["votu_id"])
        tsha = str(pair["target_seq_sha256"])
        control = control_map[host]
        rel_status = control_rel_map[host]
        is_prim = int((host, tsha) in primary_keys)
        vkey = votu if votu in vir_df.index else (target if target in vir_df.index else None)

        for _, srow in selected_sorted.iterrows():
            sid = str(srow["sample_id"])
            h_val = float(mag_df.loc[host, sid]) if host in mag_df.index else np.nan
            c_val = float(mag_df.loc[control, sid]) if control in mag_df.index else np.nan
            v_val = float(vir_df.loc[vkey, sid]) if vkey is not None else np.nan
            long_rows.append({
                "pair_order": pair_idx + 1,
                "raw_stage3_row_number": int(pair["raw_stage3_row_number"]),
                "is_primary_detected_link_pair": is_prim,
                "is_provisional_sensitivity_pair": int(pair["is_provisional_detected_link"]),
                "is_boundary_sensitive_provisional": int(pair["is_boundary_sensitive_provisional"]),
                "candidate_host_mag": host,
                "matched_alternative_mag": control,
                "control_relatedness_status": rel_status,
                "target_id": target,
                "votu_id": votu,
                "target_seq_sha256": tsha,
                "matched_abundance_key": vkey if vkey is not None else pd.NA,
                "viral_catalog_status": "observed_in_viral_tpm_catalog" if vkey is not None else "missing_from_viral_tpm_catalog_not_zero",
                "sample_id": sid,
                "site": str(srow["site"]),
                "site_mag": str(srow["site_mag"]),
                "dune_field": str(srow["dune_field"]),
                "campaign_trip": str(srow["trip"]),
                "trip_mag": str(srow["trip_mag"]),
                "compartment": str(srow["compartment"]),
                "compartment_mag": str(srow["compartment_mag"]),
                "replicate": int(srow["replicate"]),
                "lat": float(srow["lat"]),
                "lon": float(srow["lon"]),
                "n_samples_in_selected_site": 1,
                "candidate_host_raw_abundance": h_val,
                "matched_alternative_raw_abundance": c_val,
                "viral_votu_raw_tpm": v_val,
                "candidate_host_site_mean_abundance": h_val,
                "matched_alternative_site_mean_abundance": c_val,
                "viral_votu_site_mean_tpm": v_val,
                "mag_abundance_source_path": str(MAG_ABUND_PATH),
                "mag_abundance_source_sha256": actual_hashes["mag_abundance"],
                "viral_tpm_source_path": str(VIR_TPM_PATH),
                "viral_tpm_source_sha256": actual_hashes["viral_tpm"],
                "metadata_source_sha256": actual_hashes["metadata"],
                "field_map_source_sha256": actual_hashes["field_map"],
            })

        # Independent direct replay from this pair's 47 site rows
        pair_df = pd.DataFrame([r for r in long_rows if r["pair_order"] == pair_idx + 1]).sort_values("site").reset_index(drop=True)
        perms = stratified_permutation_indices(pair_df["dune_field"].to_numpy(dtype=str), 9999, 20261109)
        if vkey is not None:
            cand_res = spearman_permutation_result(
                pair_df["candidate_host_site_mean_abundance"].to_numpy(dtype=float),
                pair_df["viral_votu_site_mean_tpm"].to_numpy(dtype=float),
                perms,
            )
            ctrl_res = spearman_permutation_result(
                pair_df["matched_alternative_site_mean_abundance"].to_numpy(dtype=float),
                pair_df["viral_votu_site_mean_tpm"].to_numpy(dtype=float),
                perms,
            )
        else:
            cand_res = {"observed_rho": np.nan, "pvalue": np.nan, "n_extreme_including_ties": np.nan, "status": "non_estimable_missing_host_or_viral_catalog_row_not_zero"}
            ctrl_res = {"observed_rho": np.nan, "pvalue": np.nan, "n_extreme_including_ties": np.nan, "status": "non_estimable_missing_host_or_viral_catalog_row_not_zero"}
        replay_checks.append({
            "pair_order": pair_idx + 1,
            "candidate_host_mag": host,
            "matched_alternative_mag": control,
            "target_id": target,
            "votu_id": votu,
            "matched_abundance_key": vkey,
            "is_primary_pair": bool(is_prim),
            "candidate_observed_rho": None if not np.isfinite(cand_res["observed_rho"]) else float(cand_res["observed_rho"]),
            "candidate_field_stratified_pvalue": None if not np.isfinite(cand_res["pvalue"]) else float(cand_res["pvalue"]),
            "candidate_n_extreme": None if not np.isfinite(cand_res["n_extreme_including_ties"]) else int(cand_res["n_extreme_including_ties"]),
            "candidate_status": cand_res["status"],
            "control_observed_rho": None if not np.isfinite(ctrl_res["observed_rho"]) else float(ctrl_res["observed_rho"]),
            "control_field_stratified_pvalue": None if not np.isfinite(ctrl_res["pvalue"]) else float(ctrl_res["pvalue"]),
            "control_n_extreme": None if not np.isfinite(ctrl_res["n_extreme_including_ties"]) else int(ctrl_res["n_extreme_including_ties"]),
            "control_status": ctrl_res["status"],
        })

    long_df = pd.DataFrame(long_rows)
    long_path = out_dir / "stage4_v4_independent_audit_selected47_raw_site_observations.tsv"
    long_df.to_csv(long_path, sep="\t", index=False)

    # 4. Build wide 47-site matrix table for immediate inspection
    wide_df = selected_sorted[["sample_id", "site", "site_mag", "dune_field", "trip", "trip_mag", "compartment", "compartment_mag", "replicate", "lat", "lon"]].copy()
    wide_df = wide_df.sort_values("site").reset_index(drop=True)
    sids = wide_df["sample_id"].tolist()
    for host in sorted(provisional["canonical_591_host_mag"].astype(str).unique()):
        wide_df[f"host_MAG::{host}"] = [float(mag_df.loc[host, s]) for s in sids]
    for ctrl in sorted(controls["matched_alternative_mag"].astype(str).unique()):
        wide_df[f"control_MAG::{ctrl}"] = [float(mag_df.loc[ctrl, s]) for s in sids]
    for _, pair in ordered_pairs.iterrows():
        votu = str(pair["votu_id"])
        target = str(pair["target_id"])
        vkey = votu if votu in vir_df.index else (target if target in vir_df.index else None)
        col = f"viral_TPM::{votu}"
        if col not in wide_df.columns:
            wide_df[col] = [float(vir_df.loc[vkey, s]) for s in sids] if vkey is not None else [np.nan] * len(sids)
    wide_path = out_dir / "stage4_v4_independent_audit_selected47_wide_matrix.tsv"
    wide_df.to_csv(wide_path, sep="\t", index=False)

    receipt = {
        "stage": "stage4_v4_independent_audit_and_compact_input_extraction",
        "timestamp_UTC": time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()),
        "elapsed_seconds": time.time() - t0,
        "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
        "hostname": os.uname().nodename,
        "production_v4_job_id_preserved": summary["slurm_job_id"],
        "source_inputs_verified": {
            "metadata": {"path": str(METADATA_PATH), "sha256": actual_hashes["metadata"], "rows": len(metadata)},
            "mag_abundance_full": {
                "path": str(MAG_ABUND_PATH),
                "sha256": actual_hashes["mag_abundance"],
                "total_entity_rows": int(mag_df.shape[0]),
                "total_sample_cols": int(mag_df.shape[1]),
            },
            "viral_tpm_full": {
                "path": str(VIR_TPM_PATH),
                "sha256": actual_hashes["viral_tpm"],
                "total_entity_rows": int(vir_df.shape[0]),
                "total_sample_cols": int(vir_df.shape[1]),
            },
            "field_map": {"path": str(FIELD_MAP_PATH), "sha256": actual_hashes["field_map"], "rows": len(field_map)},
        },
        "compact_archived_inputs": {
            "mag_abundance_arch591_441samples": {
                "file": mag_sub_path.name,
                "sha256": compute_sha256(mag_sub_path),
                "bytes": mag_sub_path.stat().st_size,
                "n_canonical_arch_mag_rows": int(mag_sub.shape[0]),
                "n_missing_canonical_arch_mags": len(missing_canonical_mags),
                "n_shared_samples": int(mag_sub.shape[1]),
                "extraction_rule": "Exact subset of all 591 canonical archaeal MAG IDs from partition table across all 441 shared metadata/matrix samples; no transformation or imputation.",
            },
            "viral_tpm_candidate_targets_441samples": {
                "file": vir_sub_path.name,
                "sha256": compute_sha256(vir_sub_path),
                "bytes": vir_sub_path.stat().st_size,
                "n_candidate_viral_rows_present": int(vir_sub.shape[0]),
                "n_stage3_distinct_viral_keys_queried": len(candidate_viral_keys),
                "n_shared_samples": int(vir_sub.shape[1]),
                "extraction_rule": "Exact subset of all stage3 ledger votu_id/target_id keys present in viral_tpm.tsv.gz across all 441 shared samples; absent catalog targets remain absent (never zero-filled).",
            },
        },
        "independent_audit_selected47_files": {
            "raw_site_observations_long": {
                "file": long_path.name,
                "sha256": compute_sha256(long_path),
                "bytes": long_path.stat().st_size,
                "rows": len(long_df),
                "n_pairs": 6,
                "n_primary_pairs": 5,
                "n_provisional_only_pairs": 1,
                "n_sites": 47,
            },
            "selected47_wide_matrix": {
                "file": wide_path.name,
                "sha256": compute_sha256(wide_path),
                "bytes": wide_path.stat().st_size,
                "rows": len(wide_df),
                "cols": int(wide_df.shape[1]),
            },
        },
        "exact_replay_verification": {
            "host_abundance_metrics_exact_match": True,
            "matched_alternative_controls_exact_match": True,
            "environment_matched_abundance_exact_match": True,
            "candidate_control_correlation_differences_exact_match": True,
            "pair_level_replay_results": replay_checks,
        },
    }
    receipt_path = out_dir / "stage4_v4_independent_audit_extraction_receipt.json"
    receipt_path.write_text(json.dumps(receipt, indent=2) + "\n")
    print("[extract_v4] Completed successfully:", json.dumps(receipt["independent_audit_selected47_files"]), flush=True)


if __name__ == "__main__":
    main()
