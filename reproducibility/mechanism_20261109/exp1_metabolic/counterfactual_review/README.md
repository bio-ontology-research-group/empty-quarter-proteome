# Experiment 1 Counterfactual Closure & Targeted Candidate FVA Package (`counterfactual_review`)

## 1. Scope and Ownership Boundaries

This self-contained package (`reproducibility/mechanism_20261109/exp1_metabolic/counterfactual_review/`) implements the bounded Experiment 1 counterfactual closure and targeted candidate Flux Variability Analysis (FVA) specified in:
- `mechanism-20261109/EXP1-COUNTERFACTUAL-INDEPENDENT-AUDIT.txt`
- `mechanism-20261109/EXP1-COUNTERFACTUAL-CLOSURE-IMPLEMENTATION.txt`
- `mechanism-20261109/EXP1-COUNTERFACTUAL-ROOT-SERIALIZATION-FINDING.txt`
- `mechanism-20261109/EXP1-COUNTERFACTUAL-PRELAUNCH-REVIEW.txt`
- `mechanism-20261109/EXP1-FVA-THREAD-ORDER-DIAGNOSIS.txt`
- `mechanism-20261109/EXP1-STRATA-POSTHOC-REPAIR.txt`
- `mechanism-20261109/EXP1-COUNTERFACTUAL-FINAL-DOC-REVIEW.txt`

Strict ownership boundaries were maintained throughout:
- Frozen `code_v2/` (`exp1_sbml_highs_engine.py`, `stage2_pilot24_ensemble.py`, `stage3_full484_worker.py`) and all original Stage 3/4 outputs in `/scratch/project/k10018/eq-mechanism-20261109/exp1_metabolic/results/` were never modified.
- All closure computations used targeted HiGHS linear programs (`scipy.optimize.linprog(..., method="highs")`) on the exact frozen gapseq draft models (`*.tar.gz`), medium scenarios, biomass recipes, and accepted MILP solution sets (`sol_0 .. sol_{k-1}`) recorded in `exp1_full484_per_solution_ledger.tsv`. No full MILP reruns were performed.
- **Distinction between Counterfactual Protocol Universe (`3,872` cases) and Targeted Candidate FVA Universe (`306` cases)**:
  - The full counterfactual protocol universe in `exp1_counterfactual_complete_3872_case_ledger.tsv` consists of `484` unique genomes × `4` rich nutrient scenarios (`S3_anoxic_desert_oa_aa`, `S4_oxic_desert_oa_aa`, `S5_anoxic_defined_aa_precursors`, `S6_oxic_defined_aa_precursors`) × `2` biomass recipes (`B0_methanosarcina_template_baseline`, `B2_perturbation_omit_7_non_f420_template_cofactors_keep_f420`) = `3,872` scheduled counterfactual cases (`3,863` evaluated feasible baseline cases + `9` non-evaluable cases).
  - The targeted candidate FVA closure universe in `exp1_targeted_candidate_fva_306_cases_union_qc.tsv` and `exp1_targeted_candidate_fva_13750_instances.tsv.gz` consists of the `306` feasible `(genome, scenario_id, B2_perturbation_omit_7_non_f420_template_cofactors_keep_f420)` cases across the `81` candidate EQ MAGs in `{S3, S4, S5, S6}` that have supported candidate deficiencies (`13,750` total `(genome, scenario_id, B2, rxn_id)` instances).
- **Candidate EQ MAG Quality Qualification (`review-20261109/exp1_candidate_selection_quality_qualification.json`)**:
  - The `81` candidate EQ MAGs were selected using CheckM2 numeric thresholds `completeness >= 85%` and `contamination <= 5%` (`exp1_metabolic/code/stage3_collect_and_summarize_v4.py:634-652`; observed completeness range `85.39%..100.00%`, contamination range `0.00%..4.94%`), together with candidate reaction and model requirements.
  - Of these `81` genomes, `59` also satisfy `completeness >= 90%` and `contamination <= 5%`, while the remaining `22` have `85% <= completeness < 90%`. Preserved legacy `n_hq_eq_*` column names are retained for schema compatibility and do **not** establish formal MIMAG high-quality status.

---

## 2. Execution & Provenance History (`v1` -> `v2` -> `v2p1`)

### 2.1 Preserved `v1` Failed Thread-Order Run (`v1_failed_thread_order_record/`)
- **Slurm jobs**: `15897611` (single-process pilot), `15897849` (8-shard array), `15897850` & `15898403` (failed collector filename attempts), `15898605` (completed collector).
- **Root cause (`EXP1-FVA-THREAD-ORDER-DIAGNOSIS.txt`)**: In `15897849`, all `306` primary reconstructions verified (`max |g0 - ledger_g0| = 9.89e-17`), but all `13,750 / 13,750` FVA min/max LP pairs returned `fva_min_status = "numerical_error"` and `fva_max_status = "numerical_error"` (`raw_status = 4`). Note that the original `v1` instance table did not persist the raw solver `res.message` string per row; the subsequent thread-order diagnostic (`EXP1-FVA-THREAD-ORDER-DIAGNOSIS.txt`) reproduced the exact `raw_status = 4` failure and its solver message `"(HiGHS Status 0: Not Set)"`, showing that each spawned worker process first executed frozen `code_v2` `detect_closed_uptake_egc` / `solve_fba_highs` with default options (initializing HiGHS's process-global scheduler to hardware threads) and subsequently called `solve_verified_linear_objective_highs` with `options={"threads": 1}`.
- **Preservation**: All `v1` scripts, outputs, and Slurm logs are preserved untouched in `v1_failed_thread_order_record/` locally and `/scratch/project/k10018/eq-mechanism-20261109/exp1_metabolic/results/counterfactual_closure/` on Shaheen.

### 2.2 Executed `v2` HiGHS Single-Thread Worker-Init Run (`15901595`, `15901928`, `15901929`)
- **Executed scripts & SHA-256 digests (preserved untouched)**:
  - `counterfactual_closure_engine.py`: `a399da3875be2d1d7f4671d8d2799112376ce340507218c30547dabeefc442e4`
  - `run_targeted_candidate_fva_closure.py`: `180bdcee71cdd3445fb4f02db5394923be418ac3a1f0c2dc8c46f536c9fa80f9`
- **Fix**: `_init_worker()` calls `initialize_highs_single_thread_scheduler()` (`options={"threads": 1}`) as the very first HiGHS call in every spawned worker process before importing or calling any frozen `code_v2` function, and `v2` persists `fva_min_message` and `fva_max_message` on every instance row.
- **Slurm Job `15901595` (`exp1_cf2_pilot`, `COMPLETED 0:0`, `00:00:55`)**:
  - Replicated the uninitialized `default -> threads=1` failure (`status = 4`, `"(HiGHS Status 0: Not Set)"`) and verified the `threads=1 -> default -> threads=1` fix (`status = 0`) in fresh spawned subprocesses (`pilot_verification_receipt.json`).
  - Executed a 2-genome spawned worker pool (`1PRr2__self_SemiBin_81` and `25PRr1__self_SemiBin_1693`, `7` cases, `307` candidate instances), achieving `307 / 307` (`100%`) optimal and ordered FVA intervals (`max_top30_abs_diff_vs_orig_wkf = 9.67e-17`).
- **Slurm Job `15901928_0..7` (`exp1_cf2_fva_arr`, `COMPLETED 0:0`, `00:00:49..00:00:50`) & Initial Collector `15901929` (`exp1_cf2_collect`, `COMPLETED 0:0`, `00:00:46`)**:
  - Verified all `306 / 306` unique `(genome, scenario_id, B2)` cases and solved all `13,750 / 13,750` (`100%`) candidate FVA intervals to optimality (`0` failed/unknown).
  - Root independently audited and accepted the recorded numerical diagnostics and accounting scope in `review-20261109/counterfactual-closure-v2-root/root_independent_counterfactual_v2_audit.json` (which notes that recorded numerical diagnostics and hashes were checked without independently replaying all full sparse-flux witnesses).

### 2.3 Post-Hoc `v2p1` 244-Strata 4-Key Collector Repair (`repair_244_strata_summary_posthoc_v2p1.py`)
- **Issue diagnosed in `EXP1-STRATA-POSTHOC-REPAIR.txt`**: In initial collector `15901929`, `build_candidate_fva_summaries()` indexed `cand_inst_df` by the 2-column key `(scenario_id, rxn_id)` instead of the full 4-column stratum key `(scenario_id, rxn_id, biomass_sensitivity_classification, followup_selection_tier)`. Because `R_rxn02322_c0` appears in two `biomass_sensitivity_classification` groups (`biomass_assumption_sensitive__vanishes_when_specialized_lipids_or_homospermidine_omitted` and `biomass_assumption_sensitive__vanishes_in_core_skeleton_perturbation`) across `S3`, `S4`, `S5`, and `S6` (`8` rows totaling `306` instances: `S3` = `75 + 5 = 80`; `S4` = `62 + 2 = 64`; `S5` = `77 + 4 = 81`; `S6` = `79 + 2 = 81`), the 2-key lookup repeated the full reaction-condition counts across those 8 rows (`sum(n_instances_in_stratum) = 14,056` instead of `13,750`).
- **Post-hoc repair (`repair_244_strata_summary_posthoc_v2p1.py`, SHA-256 `906fbe3c1e678582a9c2e1b47994811748d127d7e5373890f72f1d4b44579b68`)**:
  - Preserves `counterfactual_closure_engine.py` (`a399da...`) and `run_targeted_candidate_fva_closure.py` (`180bdc...`) untouched so the exact code executed on `15901595`, `15901928`, and `15901929` is never overwritten.
  - Archives the initial `15901929` 2-key strata file (`9be0e6bed967a313497316c53c4870c3d900f20f4801536937e84e27e711a010`) and initial summary JSON (`40483c84b63ee3cb6a278d9bcb3d6b2fbfd7108ddb760abb2110f9640a2cebf9`) in `v2_initial_collector_15901929_record/`.
  - Leaves all 5 other Root-accepted output files byte-for-byte untouched and regenerates `exp1_targeted_candidate_fva_244_strata_summary.tsv` (`60192ac693d2b5b0d576005b9ca6a79e08503013f99309cc595ec304a6031382`) and `exp1_counterfactual_closure_summary.json` (`c4c9b9ad6e6c5b69b9a0d08c2726747aef43d77a4e4935ddbf9157dc43e4d44c`) using the full 4-key partition (`sum(n_instances_in_stratum) == sum(n_hq_eq_genomes_requiring) == 13,750`, with exact row-wise equality across all `244` rows).

---

## 3. Summary of Verified Results

1. **Complete `3,872`-Case Counterfactual Protocol Accounting (`exp1_counterfactual_complete_3872_case_ledger.tsv`)**:
   - `3,872` scheduled `(genome, scenario_id, biomass_recipe_id)` cases across `484` unique genomes × `4` rich nutrient scenarios (`S3`, `S4`, `S5`, `S6`) × `2` biomass recipes (`B0`, `B2`).
   - `3,863` evaluated feasible baseline cases (`692,922` rows in `exp1_full484_withdrawal_knockout_fva.tsv.gz`).
   - `9` non-evaluable cases explicitly accounted for (`8` `fba_verification_failed` + `1` `timeout_no_incumbent`).
2. **Targeted Candidate FVA Closure (`13,750` Instances across `306` Unique `(genome, scenario_id, B2)` Cases)**:
   - Primary baseline reconstruction verified on all `306 / 306` cases (`max |g0 - ledger_g0| = 9.89e-17`, `max |g0 - wkf_g0| = 9.89e-17`).
   - Mathematical union network passed declared numerical and closed-uptake EGC checks on all `306 / 306` cases (`union_egc_check_status = "evaluated_optimal"`, `has_spurious_atp_cycle = False`).
   - All `13,750 / 13,750` (`100%`) candidate instances achieved optimal, finite, ordered FVA intervals `[fva_min_flux_at_95pct_opt, fva_max_flux_at_95pct_opt]` on the mathematical union network at `>= 0.95 * g0`:
     - `3,369` instances overlapped the original alphabetical top-30 FVA cap and matched the original recorded FVA values to machine precision (`max_top30_abs_diff_vs_orig_wkf = 1.725e-13`).
     - `10,381` instances beyond rank 30 were newly evaluated in closure.
     - Interval breakdown across the `13,750` instances:
       - `strictly_positive_gt_1e-6`: `12,786` (`92.989%`)
       - `strictly_negative_lt_-1e-6`: `947` (`6.887%`)
       - `strictly_positive_trace_1e-12_to_1e-6`: `17` (`0.124%` — all `17` are `R_rxn01603_c0`, `7,8-dihydropteroate:L-glutamate ligase (ADP-forming)`, with `fva_min_flux_at_95pct_opt` in `8.92497e-7 .. 9.94603e-7`)
       - `spans_zero_contains_zero_flux`: `0` (`0.000%`)
       - `locked_zero_within_1e-12`: `0` (`0.000%`)
     - Thus `13,733 / 13,750` (`99.876%`) are strictly non-zero at `1e-6`, and `13,750 / 13,750` (`100.000%`) are strictly non-zero at `1e-12`.
3. **Corrected Counterfactual Table (`exp1_full484_counterfactual_wkf_corrected_v2.tsv.gz`, `692,922` rows)**:
   - Preserves original raw columns (`raw_perturbed_growth_primary`, `raw_addback_rescue_growth`, `raw_addback_rescues_growth`, `raw_classification`, `raw_fva_min_flux`, `raw_fva_max_flux`).
   - **Knockout rows (`598,310`)**: Unperformed primary-KO and single-reaction addback fields are set to `NaN`/`None` with explicit status `"not_performed_in_executed_code_v2"` (`0` false zeros or false booleans). FVA measurement sources are partitioned into `not_assayed_outside_top30_non_candidate_rxn` (`472,039`), `original_executed_top30_fva` (`112,521`), `verified_targeted_candidate_fva_closure_beyond_top30` (`10,381`), and `original_top30_and_verified_targeted_candidate_fva_closure` (`3,369`).
   - **Withdrawal rows (`94,612`)**: Separates partial sub-threshold union growth (`205` rows with `gw_prim < 1e-4` and `1e-4 <= gw_uni < 0.005`) from stimulatory uptake, and separates retained primary growth where addback is not applicable (`10,733` rows with `gw_prim >= 1e-4`) from threshold restoration (`83,879` rows with `gw_prim < 1e-4` and `g_addback >= 0.99 * min_growth`).

> **Scientific Scope Limit**: All `13,750` candidate FVA intervals measure selected reaction flux requirements within relaxed ensemble-union networks constrained to `>= 95%` of the primary-model optimum (`0.95 * g0`). They do not establish physiological essentiality in vivo or novel biochemical function.

---

## 4. Verification Suite

Run the portable verification suite via `uv`:

```bash
uv run pytest test_counterfactual_review.py -v
```
