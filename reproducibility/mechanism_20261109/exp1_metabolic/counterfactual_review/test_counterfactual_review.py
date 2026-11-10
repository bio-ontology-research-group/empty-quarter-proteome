#!/usr/bin/env python3
"""Pytest suite for Experiment 1 Counterfactual Closure (`counterfactual_review`).

Covers:
1. Frozen `code_v2` SHA-256 integrity using portable package-relative `frozen_code_v2_snapshot/`.
2. HiGHS process-global thread scheduler call-order defect reproduction and spawned-worker fix:
   - Reproduces `default -> threads=1` failing with HiGHS `status = 4` (`"(HiGHS Status 0: Not Set)"`)
     in an uninitialized spawned subprocess.
   - Verifies that calling `initialize_highs_single_thread_scheduler()` (`threads=1`) first in a spawned
     worker allows subsequent default-option frozen `code_v2` calls (`detect_closed_uptake_egc`,
     `solve_fba_highs`) and subsequent `threads=1` FVA solves (`solve_verified_linear_objective_highs`)
     to all succeed with `status == "optimal"`.
3. Unperformed primary-KO and addback fields in `reaction_knockout_and_fva` (`NaN`/`None` with
   `"not_performed_in_executed_code_v2"` status).
4. Corrected withdrawal classification (`nutrient_withdrawal_addback`):
   - Separates partial union growth (`gw_prim < 1e-4` and `1e-4 <= gw_uni < 0.005`, 205 rows) from stimulatory uptake.
   - Separates non-required addback (`gw_prim >= 1e-4`, 10,733 rows) from threshold restoration and failed restoration.
   - Labels `g_addback >= 0.99 * min_growth` as threshold restoration (`restored_primary_growth_above_declared_threshold_0.99_min_growth`)
     and reports continuous `addback_to_baseline_growth_ratio` separately (`EXP1-COUNTERFACTUAL-PRELAUNCH-REVIEW.txt` Item 2).
5. Fail-closed primary reconstruction gate (`EXP1-COUNTERFACTUAL-PRELAUNCH-REVIEW.txt` Item 1):
   - Failed primary solve (`g0 = NaN`) never enters `max(lb, 0.95*g0)` to use `lb=0`; produces explicit `unknown_*` FVA status.
   - Finite `g0` that mismatches `ledger_g0` or `wkf_baseline_g0` (`|g0 - ledger_g0| > 1e-6`) fails closed and produces explicit `unknown_*` FVA status.
6. Mathematical union QC classification (`EXP1-COUNTERFACTUAL-PRELAUNCH-REVIEW.txt` Item 3):
   - Never claims union is physically established; distinguishes unknown EGC from positive ATP-cycle failure.
7. Exact float64 sparse flux witness round-trip (`EXP1-COUNTERFACTUAL-PRELAUNCH-REVIEW.txt` Item 4).
8. `build_corrected_wkf_v2_table` truthfulness for top-30 overlap rows (`EXP1-FVA-THREAD-ORDER-DIAGNOSIS.txt` Section 4.4):
   - Never claims `"verified"` if the targeted rerun failed; preserves original valid top-30 FVA values untouched.
9. Fetched Shaheen v2 closure artifacts (`3,872`-case ledger, `306`-case union QC table, `13,750`-instance targeted FVA table,
   `63`-reaction summary, `244`-strata summary, `exp1_counterfactual_closure_summary.json`, and preserved `v1_failed_thread_order_record/`).
"""
from __future__ import annotations

import gzip
import json
import math
import multiprocessing as mp
import tempfile
from pathlib import Path
from typing import Any, Dict

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from counterfactual_closure_engine import (
    classify_fva_interval_v2,
    classify_reaction_knockout_counterfactual_v2,
    classify_union_network_qc,
    classify_withdrawal_counterfactual_v2,
    make_unverified_primary_fva_result,
    reproduce_highs_thread_order_defect_in_spawned_subprocesses,
    solve_verified_linear_objective_highs,
    verify_code_v2_integrity,
    verify_primary_reconstruction_gate,
)
from run_targeted_candidate_fva_closure import (
    _init_worker,
    build_corrected_wkf_v2_table,
    run_synthetic_toy_checks,
)

THIS_DIR = Path(__file__).resolve().parent
FROZEN_CODE_V2_DIR = THIS_DIR / "frozen_code_v2_snapshot"


def _spawned_worker_production_call_sequence_entry(code_v2_dir_str: str, out_q: Any) -> None:
    """Run the exact production worker initialization and LP call sequence in a fresh spawned process."""
    import pickle
    import run_targeted_candidate_fva_closure as runner

    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)
        pool_pkl = tmp / "pool.pkl"
        with open(pool_pkl, "wb") as fh:
            pickle.dump({"reactions": {}, "species": {}}, fh)

        runner._init_worker(code_v2_dir_str, str(pool_pkl), str(tmp), str(tmp))
        eng = runner._WORKER_CTX["eng"]
        init_res = runner._WORKER_CTX["highs_thread1_init_receipt"]

        mod = eng.SBMLModelData(
            model_id="toy_model",
            model_name="toy_model",
            notes_text="",
            species={
                "M_A_c0": eng.SpeciesRecord(
                    id="M_A_c0", cpd_id="cpd00027", name="A", compartment="c0", charge=0, formula="C6H12O6"
                ),
                "M_B_c0": eng.SpeciesRecord(
                    id="M_B_c0", cpd_id="cpd00029", name="B", compartment="c0", charge=0, formula="C6H12O6"
                ),
                "M_cpd00002_c0": eng.SpeciesRecord(
                    id="M_cpd00002_c0", cpd_id="cpd00002", name="ATP", compartment="c0", charge=-4, formula="C10H12N5O13P3"
                ),
                "M_cpd00001_c0": eng.SpeciesRecord(
                    id="M_cpd00001_c0", cpd_id="cpd00001", name="H2O", compartment="c0", charge=0, formula="H2O"
                ),
                "M_cpd00008_c0": eng.SpeciesRecord(
                    id="M_cpd00008_c0", cpd_id="cpd00008", name="ADP", compartment="c0", charge=-3, formula="C10H12N5O10P2"
                ),
                "M_cpd00009_c0": eng.SpeciesRecord(
                    id="M_cpd00009_c0", cpd_id="cpd00009", name="Pi", compartment="c0", charge=-2, formula="HO4P"
                ),
                "M_cpd00067_c0": eng.SpeciesRecord(
                    id="M_cpd00067_c0", cpd_id="cpd00067", name="H+", compartment="c0", charge=1, formula="H"
                ),
            },
            reactions={
                "EX_cpd00027_e0": eng.ReactionRecord(
                    id="EX_cpd00027_e0",
                    base_id="EX_cpd00027_e0",
                    name="Exchange A",
                    reversible=True,
                    lb=-10.0,
                    ub=10.0,
                    stoichiometry={"M_A_c0": -1.0},
                    is_exchange=True,
                ),
                "rxn00001_c0": eng.ReactionRecord(
                    id="rxn00001_c0",
                    base_id="rxn00001",
                    name="Mid",
                    reversible=False,
                    lb=0.0,
                    ub=10.0,
                    stoichiometry={"M_A_c0": -1.0, "M_B_c0": 1.0},
                ),
                "bio1": eng.ReactionRecord(
                    id="bio1",
                    base_id="bio1",
                    name="Biomass",
                    reversible=False,
                    lb=0.0,
                    ub=10.0,
                    stoichiometry={"M_B_c0": -1.0},
                    is_biomass=True,
                ),
            },
            objective_rxn_id="bio1",
            gene_products={},
        )

        # 1. Default-option frozen code_v2 calls (detect_closed_uptake_egc + solve_fba_highs)
        egc_res = eng.detect_closed_uptake_egc(mod)
        sol0 = eng.solve_fba_highs(mod)
        sp_ids, rxn_ids, S, lb, ub, obj_j = eng.build_sparse_stoichiometric_matrix(mod)
        target_j = rxn_ids.index("rxn00001_c0")
        min_req = 0.95 * float(sol0["growth_rate"])
        lb_fva = lb.copy()
        lb_fva[obj_j] = max(float(lb_fva[obj_j]), min_req)

        # 2. Subsequent single-thread FVA calls (options={"threads": 1})
        fva_min = solve_verified_linear_objective_highs(
            sp_ids,
            rxn_ids,
            S,
            lb_fva.copy(),
            ub.copy(),
            obj_j=obj_j,
            target_j=target_j,
            sense="min",
            min_Required_growth=min_req,
        )
        fva_max = solve_verified_linear_objective_highs(
            sp_ids,
            rxn_ids,
            S,
            lb_fva.copy(),
            ub.copy(),
            obj_j=obj_j,
            target_j=target_j,
            sense="max",
            min_Required_growth=min_req,
        )

        out_q.put(
            {
                "init_res": init_res,
                "egc_status": egc_res.get("egc_check_status"),
                "sol0_status": sol0.get("status"),
                "sol0_obj": float(sol0.get("growth_rate", math.nan)),
                "fva_min_status": fva_min["status"],
                "fva_min_raw_status": fva_min["raw_status"],
                "fva_min_flux": fva_min["target_flux"],
                "fva_max_status": fva_max["status"],
                "fva_max_raw_status": fva_max["raw_status"],
                "fva_max_flux": fva_max["target_flux"],
            }
        )


def test_frozen_code_v2_sha256_integrity() -> None:
    hashes = verify_code_v2_integrity(FROZEN_CODE_V2_DIR)
    assert len(hashes) == 3
    assert hashes["stage2_pilot24_ensemble.py"] == "bdf6e6eb1d75a1ae774c18edc40523d189602b550d40bb32aaf7dc0f718b56ec"


def test_highs_thread_order_reproduction_and_spawned_worker_fix() -> None:
    # 1. Verify direct subprocess reproduction of uninitialized default->threads=1 vs initialized threads=1->default->threads=1
    diag = reproduce_highs_thread_order_defect_in_spawned_subprocesses()
    assert diag["defect_reproduced_status4_when_default_first"] is True
    assert diag["default_first_then_threads1_probe"]["step1_default"]["status"] == 0
    assert diag["default_first_then_threads1_probe"]["step2_threads1"]["status"] == 4
    assert diag["fix_verified_all_optimal_when_threads1_first"] is True
    assert diag["threads1_first_then_default_then_threads1_probe"]["step1_threads1_init"]["status"] == 0
    assert diag["threads1_first_then_default_then_threads1_probe"]["step2_default"]["status"] == 0
    assert diag["threads1_first_then_default_then_threads1_probe"]["step3_threads1"]["status"] == 0

    # 2. Verify actual production worker sequence (_init_worker -> frozen code_v2 default LP -> threads=1 FVA)
    ctx = mp.get_context("spawn")
    q = ctx.Queue()
    p = ctx.Process(target=_spawned_worker_production_call_sequence_entry, args=(str(FROZEN_CODE_V2_DIR), q))
    p.start()
    res = q.get(timeout=30)
    p.join(timeout=10)
    assert p.exitcode == 0
    assert res["init_res"]["status"] == 0
    assert res["sol0_status"] == "optimal"
    assert res["fva_min_status"] == "optimal"
    assert res["fva_max_status"] == "optimal"
    assert abs(res["fva_min_flux"] - 9.5) < 1e-6
    assert abs(res["fva_max_flux"] - 10.0) < 1e-6


def test_reaction_knockout_unperformed_fields_are_nan_not_zero_or_true() -> None:
    res = classify_reaction_knockout_counterfactual_v2(baseline_growth=0.042, ko_growth_union=0.0)
    assert math.isnan(res["perturbed_growth_primary_v2"])
    assert math.isnan(res["addback_rescue_growth_v2"])
    assert res["addback_rescues_growth_v2"] is None
    assert res["primary_perturbation_measurement_status"] == "not_performed_in_executed_code_v2"
    assert res["addback_measurement_status"] == "not_performed_in_executed_code_v2"
    assert res["addback_rescue_disposition"] == "not_performed_in_executed_code_v2"
    assert res["corrected_classification"] == "essential_in_sampled_ensemble_union_lt_1e-4"


def test_withdrawal_partial_union_rescue_and_non_required_addback_and_threshold_label() -> None:
    # 1. Partial union growth (gw_prim < 1e-4 and 1e-4 <= gw_uni < 0.005) -> never stimulatory uptake
    r_gap = classify_withdrawal_counterfactual_v2(
        baseline_growth=0.050,
        withdrawn_growth_primary=0.0,
        withdrawn_growth_union=0.0025,
        addback_rescue_growth=0.010,
        min_growth=0.01,
    )
    assert r_gap["primary_growth_lost_lt_1e-4"] is True
    assert r_gap["corrected_classification"] == "primary_growth_lost_lt_1e-4__partial_subthreshold_union_growth_1e-4_to_0.005"
    assert "stimulatory" not in r_gap["corrected_classification"]
    # Addback tests threshold restoration (0.010 >= 0.99 * 0.01), while continuous ratio is 0.010 / 0.050 = 0.20
    assert r_gap["addback_rescues_growth_v2"] is True
    assert r_gap["addback_rescue_disposition"] == "restored_primary_growth_above_declared_threshold_0.99_min_growth"
    assert abs(r_gap["addback_to_baseline_growth_ratio"] - 0.20) <= 1e-12

    # 2. Primary growth retained (gw_prim >= 1e-4) -> addback is not applicable (None), never failed rescue
    r_ret = classify_withdrawal_counterfactual_v2(
        baseline_growth=0.050,
        withdrawn_growth_primary=0.035,
        withdrawn_growth_union=0.050,
        addback_rescue_growth=0.050,
        min_growth=0.01,
    )
    assert r_ret["primary_growth_lost_lt_1e-4"] is False
    assert r_ret["addback_rescues_growth_v2"] is None
    assert r_ret["addback_rescue_disposition"] == "not_applicable_primary_growth_retained_ge_1e-4"
    assert r_ret["corrected_classification"] == "primary_growth_retained_above_model_min_growth_ge_0.01_non_essential_or_stimulatory"

    # 3. Failed solver (NaN) -> explicit unknown, never False or stimulatory
    r_nan = classify_withdrawal_counterfactual_v2(
        baseline_growth=0.050,
        withdrawn_growth_primary=math.nan,
        withdrawn_growth_union=0.0,
        addback_rescue_growth=0.050,
    )
    assert r_nan["primary_growth_lost_lt_1e-4"] is None
    assert r_nan["addback_rescues_growth_v2"] is None
    assert r_nan["corrected_classification"] == "unknown_solver_or_invalid_baseline"


def test_fail_closed_primary_g0_nan_and_finite_mismatch_guards() -> None:
    sp_ids = ["M_A", "M_B"]
    rxn_ids = ["R_in", "R_mid", "R_bio"]
    S = sparse.csc_matrix(np.array([[1.0, -1.0, 0.0], [0.0, 1.0, -1.0]], dtype=np.float64))
    lb = np.array([0.0, 0.0, 0.0], dtype=np.float64)
    ub = np.array([10.0, 10.0, 10.0], dtype=np.float64)

    # Case A: Primary solve failed (g0 = NaN). Must fail closed and NEVER use lb=0 growth floor!
    ok_nan, reason_nan = verify_primary_reconstruction_gate(
        sol0_status="infeasible",
        prim_cyc_ok=True,
        g0=math.nan,
        ledger_g0=0.05,
        wkf_baseline_g0=0.05,
        max_residual=0.0,
        max_bound_viol=0.0,
    )
    assert ok_nan is False
    assert reason_nan == "unverified_primary_fba_solver_infeasible"

    res_nan = solve_verified_linear_objective_highs(
        sp_ids, rxn_ids, S, lb, ub, obj_j=2, target_j=1, sense="min", min_Required_growth=math.nan
    )
    assert res_nan["status"] == "unknown_unverified_primary_growth_floor"
    assert math.isnan(res_nan["target_flux"])

    res_zero_floor = solve_verified_linear_objective_highs(
        sp_ids, rxn_ids, S, lb, ub, obj_j=2, target_j=1, sense="min", min_Required_growth=0.0
    )
    assert res_zero_floor["status"] == "unknown_unverified_primary_growth_floor"
    assert math.isnan(res_zero_floor["target_flux"])

    # Case B: Finite primary g0 that mismatches ledger_g0 (|g0 - ledger_g0| > 1e-6). Must fail closed!
    ok_mis, reason_mis = verify_primary_reconstruction_gate(
        sol0_status="optimal",
        prim_cyc_ok=True,
        g0=0.050000,
        ledger_g0=0.050010,
        wkf_baseline_g0=0.050010,
        max_residual=0.0,
        max_bound_viol=0.0,
    )
    assert ok_mis is False
    assert reason_mis == "unverified_primary_g0_mismatch_vs_solution_ledger"
    unv_res = make_unverified_primary_fva_result(reason_mis)
    unv_cls = classify_fva_interval_v2(
        unv_res["target_flux"], unv_res["target_flux"], unv_res["status"], unv_res["status"]
    )
    assert unv_cls["fva_both_optimal_and_finite"] is False
    assert unv_cls["fva_strictly_nonzero_1e-6"] is None
    assert unv_cls["fva_interval_class"].startswith("unknown_fva_solver__")


def test_union_network_qc_distinguishes_unknown_egc_from_positive_atp_cycle() -> None:
    c_pass = classify_union_network_qc(
        True, "optimal", {"has_spurious_atp_cycle": False, "egc_check_status": "evaluated_optimal"}, True
    )
    c_fail = classify_union_network_qc(
        True, "optimal", {"has_spurious_atp_cycle": True, "egc_check_status": "evaluated_optimal"}, False
    )
    c_unk = classify_union_network_qc(
        True, "optimal", {"has_spurious_atp_cycle": None, "egc_check_status": "unknown_timeout"}, False
    )
    assert c_pass == "mathematical_union_passed_declared_numerical_and_closed_uptake_egc_checks"
    assert c_fail == "mathematical_union_failed_closed_uptake_atp_cycle_check"
    assert c_unk == "mathematical_union_unknown_closed_uptake_egc_check__unknown_timeout"
    assert "physically_established" not in c_pass


def test_exact_float64_sparse_flux_witness_roundtrip() -> None:
    sp_ids = ["M_A", "M_B"]
    rxn_ids = ["R_in", "R_mid", "R_bio"]
    S = sparse.csc_matrix(np.array([[1.0, -1.0, 0.0], [0.0, 1.0, -1.0]], dtype=np.float64))
    irrational_g = 1.0 / 3.0
    lb = np.array([0.0, 0.0, irrational_g], dtype=np.float64)
    ub = np.array([10.0, 10.0, 10.0], dtype=np.float64)

    res = solve_verified_linear_objective_highs(
        sp_ids, rxn_ids, S, lb, ub, obj_j=2, target_j=1, sense="min", min_Required_growth=irrational_g
    )
    assert res["status"] == "optimal"
    dumped = json.dumps(res["sparse_fluxes"], sort_keys=True, separators=(",", ":"))
    loaded = json.loads(dumped)
    assert loaded["R_mid"] == res["target_flux"]
    assert abs(loaded["R_mid"] - (1.0 / 3.0)) < 1e-15


def test_wkf_v2_top30_preserves_original_valid_fva_if_rerun_fails() -> None:
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)
        wkf_out = tmp / "out_wkf_v2.tsv.gz"

        orig_df = pd.DataFrame(
            [
                {
                    "genome": "G1",
                    "scenario_id": "S1",
                    "biomass_recipe_id": "B2",
                    "record_type": "reaction_knockout_and_fva",
                    "target_id": "R_OK",
                    "baseline_growth": 0.05,
                    "perturbed_growth_primary": 0.0,
                    "perturbed_growth_ensemble_union": 0.0,
                    "addback_rescue_growth": 0.05,
                    "essential_in_primary": True,
                    "rescued_by_ensemble_union": False,
                    "addback_rescues_growth": True,
                    "fva_min_flux": 0.12,
                    "fva_max_flux": 0.34,
                    "classification": "essential_even_in_ensemble_union",
                },
                {
                    "genome": "G1",
                    "scenario_id": "S1",
                    "biomass_recipe_id": "B2",
                    "record_type": "reaction_knockout_and_fva",
                    "target_id": "R_FAIL_RERUN",
                    "baseline_growth": 0.05,
                    "perturbed_growth_primary": 0.0,
                    "perturbed_growth_ensemble_union": 0.0,
                    "addback_rescue_growth": 0.05,
                    "essential_in_primary": True,
                    "rescued_by_ensemble_union": False,
                    "addback_rescues_growth": True,
                    "fva_min_flux": 0.50,
                    "fva_max_flux": 0.80,
                    "classification": "essential_even_in_ensemble_union",
                },
            ]
        )

        inst_df = pd.DataFrame(
            [
                {
                    "genome": "G1",
                    "scenario_id": "S1",
                    "biomass_recipe_id": "B2",
                    "rxn_id": "R_OK",
                    "fva_min_flux_at_95pct_opt": 0.12,
                    "fva_max_flux_at_95pct_opt": 0.34,
                    "fva_both_optimal_and_finite": True,
                    "fva_interval_ordered": True,
                    "fva_min_status": "optimal",
                    "fva_max_status": "optimal",
                    "fva_interval_class": "strictly_positive_flux_gt_1e-6",
                },
                {
                    "genome": "G1",
                    "scenario_id": "S1",
                    "biomass_recipe_id": "B2",
                    "rxn_id": "R_FAIL_RERUN",
                    "fva_min_flux_at_95pct_opt": math.nan,
                    "fva_max_flux_at_95pct_opt": math.nan,
                    "fva_both_optimal_and_finite": False,
                    "fva_interval_ordered": False,
                    "fva_min_status": "numerical_error",
                    "fva_max_status": "numerical_error",
                    "fva_interval_class": "unknown_fva_solver__min_numerical_error__max_numerical_error",
                },
            ]
        )

        stats = build_corrected_wkf_v2_table(
            orig_df, inst_df, wkf_out, expected_wkf_rows=2, expected_cand_matched=2
        )
        src_counts = stats["knockout_fva_measurement_source_counts"]
        assert src_counts.get("original_top30_and_verified_targeted_candidate_fva_closure") == 1
        assert src_counts.get("original_top30_valid__targeted_rerun_failed_numerical_error_numerical_error") == 1

        out_df = pd.read_csv(wkf_out, sep="\t")
        r_ok = out_df[out_df["target_id"] == "R_OK"].iloc[0]
        assert r_ok["fva_measurement_source_v2"] == "original_top30_and_verified_targeted_candidate_fva_closure"
        assert abs(float(r_ok["fva_min_flux_v2"]) - 0.12) < 1e-12

        r_fail = out_df[out_df["target_id"] == "R_FAIL_RERUN"].iloc[0]
        assert "verified" not in str(r_fail["fva_measurement_source_v2"])
        assert str(r_fail["fva_measurement_source_v2"]).startswith("original_top30_valid__targeted_rerun_failed_")
        # Original valid top-30 FVA must be preserved, not overwritten with NaN
        assert abs(float(r_fail["fva_min_flux_v2"]) - 0.50) < 1e-12
        assert abs(float(r_fail["fva_max_flux_v2"]) - 0.80) < 1e-12


def test_synthetic_toy_suite_passes() -> None:
    out = run_synthetic_toy_checks()
    assert out["all_toy_assertions_passed"] is True


def test_fetched_v2_closure_artifacts_when_present() -> None:
    summary_path = THIS_DIR / "exp1_counterfactual_closure_summary.json"
    if not summary_path.exists():
        pytest.skip("v2 closure summary not yet fetched to local package directory")

    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    assert summary["closure_version"] == "v2_highs_thread1_worker_init__with_v2p1_posthoc_244_strata_4key_repair"
    assert summary["executed_closure_code_sha256"] == {
        "counterfactual_closure_engine.py": "a399da3875be2d1d7f4671d8d2799112376ce340507218c30547dabeefc442e4",
        "run_targeted_candidate_fva_closure.py": "180bdcee71cdd3445fb4f02db5394923be418ac3a1f0c2dc8c46f536c9fa80f9",
    }
    assert "posthoc_244_strata_repair_receipt" in summary
    repair_rec = summary["posthoc_244_strata_repair_receipt"]
    assert repair_rec["sum_n_instances_in_stratum"] == 13750
    assert repair_rec["sum_n_hq_eq_genomes_requiring"] == 13750
    assert repair_rec["exact_rowwise_equality_with_n_hq_eq_genomes_requiring"] is True
    assert (
        repair_rec["preserved_initial_15901929_2key_strata_summary_sha256"]
        == "9be0e6bed967a313497316c53c4870c3d900f20f4801536937e84e27e711a010"
    )

    assert summary["complete_3872_case_ledger"]["n_scheduled_cases_total"] == 3872
    assert summary["complete_3872_case_ledger"]["n_evaluated_feasible_cases"] == 3863
    assert summary["complete_3872_case_ledger"]["n_non_evaluable_cases"] == 9

    tcf = summary["targeted_candidate_fva_closure"]
    assert tcf["n_unique_b2_cases_accounted"] == 306
    assert tcf["n_primary_reconstructions_verified"] == 306
    assert tcf["n_candidate_instances_accounted"] == 13750
    assert tcf["n_fva_both_optimal_and_ordered"] == 13750
    assert tcf["n_fva_failed_or_unknown"] == 0
    assert tcf["n_previously_assayed_in_top30"] == 3369
    assert tcf["n_newly_evaluated_beyond_top30"] == 10381
    assert tcf["max_top30_abs_diff_vs_orig_wkf"] <= 1e-6

    inst_df = pd.read_csv(THIS_DIR / "exp1_targeted_candidate_fva_13750_instances.tsv.gz", sep="\t")
    assert len(inst_df) == 13750
    assert "fva_min_message" in inst_df.columns
    assert "fva_max_message" in inst_df.columns
    assert (inst_df["fva_min_status"] == "optimal").all()
    assert (inst_df["fva_max_status"] == "optimal").all()

    # Check the 17 trace-positive instances (all R_rxn01603_c0, min flux 8.92497e-7..9.94603e-7)
    trace_df = inst_df[inst_df["fva_interval_class"] == "strictly_positive_trace_1e-12_to_1e-6"]
    assert len(trace_df) == 17
    assert set(trace_df["rxn_id"].unique()) == {"R_rxn01603_c0"}
    assert float(trace_df["fva_min_flux_at_95pct_opt"].min()) > 8.9e-7
    assert float(trace_df["fva_min_flux_at_95pct_opt"].max()) < 1.0e-6

    qc_df = pd.read_csv(THIS_DIR / "exp1_targeted_candidate_fva_306_cases_union_qc.tsv", sep="\t")
    assert len(qc_df) == 306
    assert (qc_df["worker_highs_thread1_init_status"] == 0).all()

    rxn_df = pd.read_csv(THIS_DIR / "exp1_targeted_candidate_fva_63_reactions_summary.tsv", sep="\t")
    assert len(rxn_df) == 63
    assert int(rxn_df["n_candidate_instances_total"].sum()) == 13750
    assert int(rxn_df["n_fva_both_optimal_and_ordered"].sum()) == 13750

    strata_df = pd.read_csv(THIS_DIR / "exp1_targeted_candidate_fva_244_strata_summary.tsv", sep="\t")
    assert len(strata_df) == 244
    assert (
        len(
            strata_df.drop_duplicates(
                ["scenario_id", "rxn_id", "biomass_sensitivity_classification", "followup_selection_tier"]
            )
        )
        == 244
    )
    assert (strata_df["n_instances_in_stratum"] == strata_df["n_hq_eq_genomes_requiring"]).all()
    assert int(strata_df["n_instances_in_stratum"].sum()) == 13750
    assert int(strata_df["n_fva_both_optimal_and_ordered"].sum()) == 13750
    assert int(strata_df["n_assayed_in_original_top30"].sum()) == 3369
    assert int(strata_df["n_newly_evaluated_in_closure"].sum()) == 10381
    assert int(strata_df["n_fva_strictly_nonzero_1e-6"].sum()) == 13733
    assert int(strata_df["n_fva_trace_positive_1e-12_to_1e-6"].sum()) == 17

    rxn2322_strata = strata_df[strata_df["rxn_id"] == "R_rxn02322_c0"]
    assert len(rxn2322_strata) == 8
    assert int(rxn2322_strata["n_instances_in_stratum"].sum()) == 306

    ledger_df = pd.read_csv(THIS_DIR / "exp1_counterfactual_complete_3872_case_ledger.tsv", sep="\t")
    assert len(ledger_df) == 3872

    # Verify that v1_failed_thread_order_record and v2_initial_collector_15901929_record are preserved untouched
    v1_dir = THIS_DIR / "v1_failed_thread_order_record"
    assert (v1_dir / "exp1_counterfactual_closure_summary.json").exists()
    v1_summary = json.loads((v1_dir / "exp1_counterfactual_closure_summary.json").read_text(encoding="utf-8"))
    assert v1_summary["targeted_candidate_fva_closure"]["n_fva_failed_or_unknown"] == 13750

    v2_init_dir = THIS_DIR / "v2_initial_collector_15901929_record"
    assert (v2_init_dir / "exp1_targeted_candidate_fva_244_strata_summary_v2_2key_15901929.tsv").exists()
