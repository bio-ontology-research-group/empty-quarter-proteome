#!/usr/bin/env python3
"""Post-hoc collector repair (`v2p1`) for `exp1_targeted_candidate_fva_244_strata_summary.tsv`
(`mechanism-20261109/EXP1-STRATA-POSTHOC-REPAIR.txt`).

Scope and Provenance Rules:
1. Does NOT modify or overwrite `counterfactual_closure_engine.py`
   (`a399da3875be2d1d7f4671d8d2799112376ce340507218c30547dabeefc442e4`) or
   `run_targeted_candidate_fva_closure.py`
   (`180bdcee71cdd3445fb4f02db5394923be418ac3a1f0c2dc8c46f536c9fa80f9`) that were executed
   in Slurm pilot `15901595`, 8-shard array `15901928`, and initial collector `15901929`.
2. Does NOT re-run any LP/FVA or touch the 8 completed shards or the 5 Root-accepted output tables:
   - `exp1_targeted_candidate_fva_306_cases_union_qc.tsv` (`e20ddc61260e05e5bdfbeade566b867bcfccd48f5cabc81194ab57548781e228`)
   - `exp1_targeted_candidate_fva_13750_instances.tsv.gz` (`afa12bfa8b9c188e6298692ee6dc15604dd1786ecd4f8d153ddbd8c1586a69c1`)
   - `exp1_targeted_candidate_fva_63_reactions_summary.tsv` (`e562777802d0a5f5623e3c5538780eaea2f74b5207009d642906b0f4edfb746e`)
   - `exp1_counterfactual_complete_3872_case_ledger.tsv` (`5db32419b2b803346022d093d3b2a94067e620179bdd989b3216bcd2137d44fe`)
   - `exp1_full484_counterfactual_wkf_corrected_v2.tsv.gz` (`2a3b203a2567bca7ce5c61f0166ebbde8ba76ac4ae566721a17561d809bbabdd`)
3. Repairs ONLY the derived 244-strata summary table (`exp1_targeted_candidate_fva_244_strata_summary.tsv`)
   and summary metadata (`exp1_counterfactual_closure_summary.json`) by partitioning the 13,750 instances
   using the full 4-column stratum key:
   `(scenario_id, rxn_id, biomass_sensitivity_classification, followup_selection_tier)`
   instead of `(scenario_id, rxn_id)` alone (which double-counted `R_rxn02322_c0` across its two
   `biomass_sensitivity_classification` / `followup_selection_tier` strata in S3, S4, S5, S6, yielding 14,056
   instead of 13,750).
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path
from typing import Any

import pandas as pd

from counterfactual_closure_engine import compute_file_sha256

THIS_DIR = Path(__file__).resolve().parent

EXPECTED_EXECUTED_V2_CODE_HASHES: dict[str, str] = {
    "counterfactual_closure_engine.py": "a399da3875be2d1d7f4671d8d2799112376ce340507218c30547dabeefc442e4",
    "run_targeted_candidate_fva_closure.py": "180bdcee71cdd3445fb4f02db5394923be418ac3a1f0c2dc8c46f536c9fa80f9",
}

EXPECTED_UNCHANGED_OUTPUT_HASHES: dict[str, str] = {
    "exp1_targeted_candidate_fva_306_cases_union_qc.tsv": "e20ddc61260e05e5bdfbeade566b867bcfccd48f5cabc81194ab57548781e228",
    "exp1_targeted_candidate_fva_13750_instances.tsv.gz": "afa12bfa8b9c188e6298692ee6dc15604dd1786ecd4f8d153ddbd8c1586a69c1",
    "exp1_targeted_candidate_fva_63_reactions_summary.tsv": "e562777802d0a5f5623e3c5538780eaea2f74b5207009d642906b0f4edfb746e",
    "exp1_counterfactual_complete_3872_case_ledger.tsv": "5db32419b2b803346022d093d3b2a94067e620179bdd989b3216bcd2137d44fe",
    "exp1_full484_counterfactual_wkf_corrected_v2.tsv.gz": "2a3b203a2567bca7ce5c61f0166ebbde8ba76ac4ae566721a17561d809bbabdd",
}

INITIAL_15901929_2KEY_STRATA_SHA256 = "9be0e6bed967a313497316c53c4870c3d900f20f4801536937e84e27e711a010"

ORIGINAL_CANDIDATE_SUMMARY_14_COLS: list[str] = [
    "scenario_id",
    "rxn_id",
    "base_id",
    "rxn_name",
    "is_transport",
    "biomass_sensitivity_classification",
    "followup_selection_tier",
    "epistemic_status",
    "n_hq_eq_genomes_requiring",
    "n_hq_eq_families_requiring",
    "median_freq_B2",
    "median_candidate_weight",
    "example_genomes_csv",
    "families_csv",
]

STRATUM_4KEY_COLS: list[str] = [
    "scenario_id",
    "rxn_id",
    "biomass_sensitivity_classification",
    "followup_selection_tier",
]


def repair_244_strata_summary(
    out_dir: Path,
    orig_summary_path: Path | None = None,
) -> dict[str, Any]:
    t0 = time.time()

    # 1. Verify that executed v2 code files were NOT overwritten
    for fname, exp_hash in EXPECTED_EXECUTED_V2_CODE_HASHES.items():
        actual_hash = compute_file_sha256(THIS_DIR / fname)
        if actual_hash != exp_hash:
            raise RuntimeError(
                f"Executed v2 script {fname} hash mismatch: expected {exp_hash}, got {actual_hash}"
            )

    # 2. Verify that all 5 Root-accepted non-strata outputs are untouched
    for fname, exp_hash in EXPECTED_UNCHANGED_OUTPUT_HASHES.items():
        fpath = out_dir / fname
        actual_hash = compute_file_sha256(fpath)
        if actual_hash != exp_hash:
            raise RuntimeError(
                f"Root-accepted artifact {fname} hash mismatch: expected {exp_hash}, got {actual_hash}"
            )

    # 3. Preserve initial 15901929 2-key strata summary and summary JSON if not already archived
    archive_dir = out_dir / "v2_initial_collector_15901929_record"
    archive_dir.mkdir(parents=True, exist_ok=True)
    strata_out_path = out_dir / "exp1_targeted_candidate_fva_244_strata_summary.tsv"
    summary_json_path = out_dir / "exp1_counterfactual_closure_summary.json"

    archived_strata_path = archive_dir / "exp1_targeted_candidate_fva_244_strata_summary_v2_2key_15901929.tsv"
    archived_summary_path = archive_dir / "exp1_counterfactual_closure_summary_v2_15901929.json"

    if strata_out_path.is_file() and compute_file_sha256(strata_out_path) == INITIAL_15901929_2KEY_STRATA_SHA256:
        if not archived_strata_path.is_file():
            archived_strata_path.write_bytes(strata_out_path.read_bytes())
        if summary_json_path.is_file() and not archived_summary_path.is_file():
            archived_summary_path.write_bytes(summary_json_path.read_bytes())

    if not archived_strata_path.is_file():
        raise FileNotFoundError(
            f"Missing preserved initial 15901929 2-key strata file at {archived_strata_path}"
        )
    if compute_file_sha256(archived_strata_path) != INITIAL_15901929_2KEY_STRATA_SHA256:
        raise RuntimeError(
            f"Preserved initial 15901929 2-key strata file hash mismatch at {archived_strata_path}"
        )

    # 4. Load the 14 original base columns of the 244 strata summary
    if orig_summary_path is not None and orig_summary_path.is_file():
        orig_sum_df = pd.read_csv(orig_summary_path, sep="\t")[ORIGINAL_CANDIDATE_SUMMARY_14_COLS].copy()
    else:
        orig_sum_df = pd.read_csv(archived_strata_path, sep="\t")[ORIGINAL_CANDIDATE_SUMMARY_14_COLS].copy()

    if len(orig_sum_df) != 244:
        raise RuntimeError(f"Expected 244 rows in original candidate summary, got {len(orig_sum_df)}")
    if int(orig_sum_df["n_hq_eq_genomes_requiring"].sum()) != 13750:
        raise RuntimeError(
            f"Expected sum(n_hq_eq_genomes_requiring) == 13750, got {orig_sum_df['n_hq_eq_genomes_requiring'].sum()}"
        )
    if len(orig_sum_df.drop_duplicates(STRATUM_4KEY_COLS)) != 244:
        raise RuntimeError("4-key (scenario_id, rxn_id, biomass_sensitivity_classification, followup_selection_tier) is not unique across the 244 strata!")

    # 5. Load the Root-accepted 13,750 instances table and partition by the full 4-key
    cand_inst_df = pd.read_csv(out_dir / "exp1_targeted_candidate_fva_13750_instances.tsv.gz", sep="\t")
    if len(cand_inst_df) != 13750:
        raise RuntimeError(f"Expected 13,750 rows in instances table, got {len(cand_inst_df)}")

    inst_by_stratum: dict[tuple[str, str, str, str], pd.DataFrame] = {
        (str(s), str(r), str(bsc), str(tier)): grp
        for (s, r, bsc, tier), grp in cand_inst_df.groupby(STRATUM_4KEY_COLS, sort=False)
    }
    if len(inst_by_stratum) != 244:
        raise RuntimeError(f"Expected 244 unique 4-key strata in instances table, got {len(inst_by_stratum)}")

    seen_indices: set[int] = set()
    strata_rows: list[dict[str, Any]] = []
    rxn02322_strata_records: list[dict[str, Any]] = []

    for _, srow in orig_sum_df.iterrows():
        key = (
            str(srow["scenario_id"]),
            str(srow["rxn_id"]),
            str(srow["biomass_sensitivity_classification"]),
            str(srow["followup_selection_tier"]),
        )
        if key not in inst_by_stratum:
            raise KeyError(f"Stratum 4-key {key} missing from 13,750 instances table")
        grp = inst_by_stratum[key]
        grp_indices = set(grp.index.tolist())
        if seen_indices & grp_indices:
            raise RuntimeError(f"Overlapping instance indices in stratum {key}")
        seen_indices |= grp_indices

        n_tot = len(grp)
        expected_n = int(srow["n_hq_eq_genomes_requiring"])
        if n_tot != expected_n:
            raise RuntimeError(
                f"Stratum {key} instance count {n_tot} != n_hq_eq_genomes_requiring {expected_n}"
            )

        n_orig = int((grp["was_in_original_top30_fva"] == True).sum())
        n_opt = int((grp["fva_both_optimal_and_finite"] == True).sum())
        n_unk = int((grp["fva_both_optimal_and_finite"] == False).sum())
        n_str_1e6 = int((grp["fva_strictly_nonzero_1e-6"] == True).sum())
        n_str_1e12 = int((grp["fva_strictly_nonzero_1e-12"] == True).sum())
        n_trace = int((grp["fva_interval_class"] == "strictly_positive_trace_1e-12_to_1e-6").sum())
        n_spans0 = int((grp["fva_interval_class"] == "spans_zero_contains_zero_flux").sum())
        n_uni_cyc_pass = int(
            (
                grp["union_network_qc_class"]
                == "mathematical_union_passed_declared_numerical_and_closed_uptake_egc_checks"
            ).sum()
        )
        n_uni_cyc_fail = int(
            (grp["union_network_qc_class"] == "mathematical_union_failed_closed_uptake_atp_cycle_check").sum()
        )

        d = srow.to_dict()
        d.update(
            {
                "n_instances_in_stratum": n_tot,
                "n_assayed_in_original_top30": n_orig,
                "n_newly_evaluated_in_closure": n_tot - n_orig,
                "n_fva_both_optimal_and_ordered": n_opt,
                "n_fva_unknown_or_failed": n_unk,
                "n_fva_strictly_nonzero_1e-6": n_str_1e6,
                "n_fva_strictly_nonzero_1e-12": n_str_1e12,
                "n_fva_trace_positive_1e-12_to_1e-6": n_trace,
                "n_fva_spans_zero": n_spans0,
                "n_instances_in_egc_passing_mathematical_union": n_uni_cyc_pass,
                "n_instances_in_egc_failing_mathematical_union": n_uni_cyc_fail,
                "n_instances_in_egc_unknown_mathematical_union": n_tot - n_uni_cyc_pass - n_uni_cyc_fail,
                "median_fva_min_flux": float(grp["fva_min_flux_at_95pct_opt"].median()),
                "median_fva_max_flux": float(grp["fva_max_flux_at_95pct_opt"].median()),
                "median_fva_span": float(grp["fva_span"].median()),
            }
        )
        strata_rows.append(d)

        if str(srow["rxn_id"]) == "R_rxn02322_c0":
            rxn02322_strata_records.append(
                {
                    "scenario_id": str(srow["scenario_id"]),
                    "rxn_id": str(srow["rxn_id"]),
                    "biomass_sensitivity_classification": str(srow["biomass_sensitivity_classification"]),
                    "followup_selection_tier": str(srow["followup_selection_tier"]),
                    "n_hq_eq_genomes_requiring": expected_n,
                    "n_instances_in_stratum": n_tot,
                    "n_fva_both_optimal_and_ordered": n_opt,
                    "n_fva_strictly_nonzero_1e-6": n_str_1e6,
                }
            )

    if len(seen_indices) != 13750:
        raise RuntimeError(f"Expected disjoint partition of 13,750 instances, got {len(seen_indices)}")

    strata_244_df = pd.DataFrame(strata_rows)
    if int(strata_244_df["n_instances_in_stratum"].sum()) != 13750:
        raise RuntimeError("Sum of n_instances_in_stratum across 244 strata != 13,750")
    if int(strata_244_df["n_fva_both_optimal_and_ordered"].sum()) != 13750:
        raise RuntimeError("Sum of n_fva_both_optimal_and_ordered across 244 strata != 13,750")
    if int(strata_244_df["n_assayed_in_original_top30"].sum()) != 3369:
        raise RuntimeError("Sum of n_assayed_in_original_top30 across 244 strata != 3,369")
    if int(strata_244_df["n_newly_evaluated_in_closure"].sum()) != 10381:
        raise RuntimeError("Sum of n_newly_evaluated_in_closure across 244 strata != 10,381")
    if int(strata_244_df["n_fva_strictly_nonzero_1e-6"].sum()) != 13733:
        raise RuntimeError("Sum of n_fva_strictly_nonzero_1e-6 across 244 strata != 13,733")
    if int(strata_244_df["n_fva_trace_positive_1e-12_to_1e-6"].sum()) != 17:
        raise RuntimeError("Sum of n_fva_trace_positive_1e-12_to_1e-6 across 244 strata != 17")

    strata_244_df.to_csv(strata_out_path, sep="\t", index=False)
    repaired_strata_sha256 = compute_file_sha256(strata_out_path)

    # 6. Verify the 17 trace-positive instances (all R_rxn01603_c0)
    trace_df = cand_inst_df[cand_inst_df["fva_interval_class"] == "strictly_positive_trace_1e-12_to_1e-6"]
    trace_rxns = sorted(trace_df["rxn_id"].unique().tolist())
    if len(trace_df) != 17 or trace_rxns != ["R_rxn01603_c0"]:
        raise RuntimeError(f"Unexpected trace instances: count={len(trace_df)}, rxns={trace_rxns}")

    # 7. Update exp1_counterfactual_closure_summary.json with honest v2 + posthoc v2p1 provenance
    base_summary = json.loads(archived_summary_path.read_text(encoding="utf-8"))
    this_script_sha256 = compute_file_sha256(Path(__file__).resolve())

    base_summary["closure_version"] = "v2_highs_thread1_worker_init__with_v2p1_posthoc_244_strata_4key_repair"
    base_summary["posthoc_244_strata_repair_receipt"] = {
        "repair_timestamp_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "repair_elapsed_s": round(time.time() - t0, 3),
        "repair_script_name": Path(__file__).name,
        "repair_script_sha256": this_script_sha256,
        "preserved_v2_pilot_15901595_array_15901928_collector_15901929_code_sha256": EXPECTED_EXECUTED_V2_CODE_HASHES,
        "preserved_initial_15901929_2key_strata_summary_sha256": INITIAL_15901929_2KEY_STRATA_SHA256,
        "repaired_4key_strata_summary_sha256": repaired_strata_sha256,
        "defect_fixed": (
            "Initial collector 15901929 grouped cand_inst_df by 2-key (scenario_id, rxn_id) instead of full 4-key "
            "(scenario_id, rxn_id, biomass_sensitivity_classification, followup_selection_tier), duplicating the 306 "
            "instances of R_rxn02322_c0 across its 2 biomass_sensitivity_classification strata in S3, S4, S5, S6 "
            "(sum=14,056 instead of 13,750). Post-hoc repair groups by the full 4-key without re-running any LPs."
        ),
        "stratum_4key_columns": STRATUM_4KEY_COLS,
        "n_strata": 244,
        "sum_n_hq_eq_genomes_requiring": int(strata_244_df["n_hq_eq_genomes_requiring"].sum()),
        "sum_n_instances_in_stratum": int(strata_244_df["n_instances_in_stratum"].sum()),
        "exact_rowwise_equality_with_n_hq_eq_genomes_requiring": bool(
            (strata_244_df["n_instances_in_stratum"] == strata_244_df["n_hq_eq_genomes_requiring"]).all()
        ),
        "sum_n_assayed_in_original_top30": int(strata_244_df["n_assayed_in_original_top30"].sum()),
        "sum_n_newly_evaluated_in_closure": int(strata_244_df["n_newly_evaluated_in_closure"].sum()),
        "sum_n_fva_both_optimal_and_ordered": int(strata_244_df["n_fva_both_optimal_and_ordered"].sum()),
        "sum_n_fva_strictly_nonzero_1e-6": int(strata_244_df["n_fva_strictly_nonzero_1e-6"].sum()),
        "sum_n_fva_trace_positive_1e-12_to_1e-6": int(strata_244_df["n_fva_trace_positive_1e-12_to_1e-6"].sum()),
        "rxn02322_split_strata_records": rxn02322_strata_records,
        "trace_17_instances_summary": {
            "n_trace_instances": int(len(trace_df)),
            "rxn_id": "R_rxn01603_c0",
            "rxn_name": str(trace_df["rxn_name"].iloc[0]),
            "followup_selection_tier": str(trace_df["followup_selection_tier"].iloc[0]),
            "min_fva_min_flux_at_95pct_opt": float(trace_df["fva_min_flux_at_95pct_opt"].min()),
            "max_fva_min_flux_at_95pct_opt": float(trace_df["fva_min_flux_at_95pct_opt"].max()),
            "min_fva_max_flux_at_95pct_opt": float(trace_df["fva_max_flux_at_95pct_opt"].min()),
            "max_fva_max_flux_at_95pct_opt": float(trace_df["fva_max_flux_at_95pct_opt"].max()),
            "scientific_interpretation": (
                "Selected reaction requirements within relaxed ensemble-union networks at 95% of primary-model "
                "optimum (min flux 8.925e-7 to 9.946e-7, just below the 1e-6 threshold); not novel biochemical "
                "function or physiological essentiality."
            ),
        },
    }
    base_summary["output_artifacts_sha256"]["exp1_targeted_candidate_fva_244_strata_summary.tsv"] = repaired_strata_sha256

    summary_json_path.write_text(json.dumps(base_summary, indent=2) + "\n", encoding="utf-8")
    return base_summary["posthoc_244_strata_repair_receipt"]


def main() -> None:
    parser = argparse.ArgumentParser(description="Post-hoc v2p1 repair of 244-strata summary table.")
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=THIS_DIR,
        help="Directory containing the v2 closure outputs (default: package directory).",
    )
    parser.add_argument(
        "--orig-summary-tsv",
        type=Path,
        default=None,
        help="Optional path to original exp1_stage4_supported_candidate_deficiencies_summary.tsv.",
    )
    args = parser.parse_args()
    receipt = repair_244_strata_summary(out_dir=args.out_dir, orig_summary_path=args.orig_summary_tsv)
    print(
        f"[Post-hoc 244-Strata Repair v2p1 Complete] "
        f"repaired_sha256={receipt['repaired_4key_strata_summary_sha256']} "
        f"sum_instances={receipt['sum_n_instances_in_stratum']}"
    )


if __name__ == "__main__":
    main()
