# Hypothesis 4 (Structure) Post-Processing Repairs & Recovery Plan

**Target Manuscript**: Nature Communications (Dark-Matter MAGs Proteome)  
**Location**: `reviews/darkmatter-2026-10-13/round7/h4_postprocess_repairs/`  
**Date**: 2026-10-14  
**Status**: Repaired candidates prepared, synthetically verified, ready for Root review and promotion.

---

## 1. Summary of Identified Defects and Root Cause Analysis

### Defect 1: Definite NameError in `audit_hit_cap_saturation`
- **Location**: `v7_crop_and_search_regions.py`, line 631
- **Root Cause**: In Python 3, list comprehensions have their own local scope. The list comprehension:
  ```python
  prob_k = [
      qid for qid in saturated
      if a["total_hits_raw"] > 0 and a["rem_counts"][k] == 0
  ]
  ```
  referenced `a`, which was defined in the preceding comprehension `[qid for qid, a in acc.items() if ...]`. Because comprehension variables do not leak into outer scope in Python 3, executing this line raises `NameError: name 'a' is not defined`.
- **Impact**: In job `53738540`, both massive Foldseek searches (whole-chain on 5,851 queries and cropped confident regions on 7,773 queries against PDB) run cleanly, but the script terminates with error code 1 as soon as it enters `audit_hit_cap_saturation` at Step 5.
- **Repair**: Explicitly index `a = acc[qid]`. Tested and verified in `test_repairs.py`.

---

### Defect 2: Hit-Cap Saturation Scope & Per-Query Cap Awareness
- **Root Cause**: The original code only audited `rem_counts[k] == 0` (`POST_EXCLUSION_EMPTY`). However, queries that are saturated at 1,000 hits and have `LOW_SCORE` ($qTM < 0.5$) may have had their true homologs truncated by the 1,000-hit prefilter cap. Furthermore, queries evaluated after an initial rerun at 5,000 hits must be evaluated against their *actual* per-query cap (5,000), not a hardcoded 1,000.
- **Repair**:
  1. Saturation is evaluated against `a.get("hit_cap", MAX_SEQS)`.
  2. The audit categorizes saturated queries into `problematic_exclusion_empty` (`rem_counts[k] == 0`), `problematic_low_score` (`rem_counts[k] > 0` and $qTM < 0.5$), and `problematic_all_unmatched` (union of both).
  3. Queries lacking `MATCHED` at any tier while saturated are flagged for adaptive rerun at `--max-seqs 5000`.

---

### Defect 3: Rerun Collision, Evidence Overwriting, and Score Monotonicity
- **Root Cause**:
  1. The original rerun function used hardcoded identical paths for whole-chain and crop reruns (`results/rerun_hits.raw.tsv` and `rerun_queries_tmp/`), meaning a subsequent crop rerun would overwrite whole-chain raw evidence.
  2. The original code replaced accumulators whenever `total_hits_raw` increased, which could silently discard a higher-scoring original match if Foldseek's alignment order or heuristics shifted under `--max-seqs 5000`.
- **Repair**:
  1. Strictly separated paths: `v7_wc_rerun_queries/`, `v7_wc_rerun_hits.raw.tsv`, `v7_wc_rerun_fs_tmp/` vs `v7_crop_rerun_queries/`, `v7_crop_rerun_hits.raw.tsv`, `v7_crop_rerun_fs_tmp/`.
  2. Merges original and rerun hit rows with explicit target/coordinate deduplication `(target, qstart, qend, tstart, tend)`.
  3. Enforces score monotonicity: asserts and guarantees `final_top_qtm >= initial_top_qtm - 1e-9`.
  4. Stores `initial_cap`, `final_cap`, `initial_statuses`, `final_statuses`, and `score_monotonicity_verified` for every rerun query in the audit log.

---

### Defect 4: Finite-Population Census vs Sampling-Error CIs
- **Root Cause**: The study evaluates the complete, enumerated finite cohort ($N=5,851$ models, $R=7,773$ regions). Presenting Wilson binomial confidence intervals as "uncertainty" in primary census and figure tables is statistically misleading because there is no sampling error for the enumerated cohort.
- **Repair**:
  1. Primary output `v7_stratum_census.tsv` and figure tidy data `v7_figure_data.tsv` contain exact finite-population counts, denominators, and proportions without sampling CIs.
  2. Undefined conditioning domains (where denominator $n \le 0$, e.g., whole-chain unmatched in $Q \setminus K$ PDB or matched in $K_{dual\_none}$) return `None` (`null` in JSON, empty string in TSV), eliminating degenerate `[0.0, 0.0]` intervals.
  3. Optional Wilson intervals and Fisher/$\chi^2$ p-values are placed under a separate JSON key `model_based_superpopulation_inferential_statistics` with an explicit disclaimer that they model a hypothetical infinite superpopulation under independence assumptions, which may not hold across related protein families.

---

### Defect 5: Checksum Manifest Self-Inclusion and Fail-Open Verifier
- **Root Cause**:
  1. Shell redirection `sha256sum v7_*.tsv v7_*.json > v7_results.sha256.tsv` causes the output file to match the wildcard pattern on reruns, including itself in the manifest.
  2. The closure verifier used `if os.path.exists(fpath):`, silently skipping missing required artifacts (failing open).
- **Repair**:
  1. `generate_sha256_manifest()` lists all `v7_*.tsv` and `v7_*.json` in `results/`, explicitly excluding `v7_results.sha256.tsv`, writes to a temporary file, and performs an atomic rename.
  2. `v7_verify_closure.py` fails closed: every required artifact must exist and match `v7_results.sha256.tsv`. Missing files or checksum mismatches raise `AssertionError`.
  3. `v7_verify_closure.py` writes a verification receipt `v7_closure_verification.json` that is included in the final manifest, avoiding circularity.

---

## 2. Inventory of Candidate Files

All candidate files have been placed exclusively in:
`reviews/darkmatter-2026-10-13/round7/h4_postprocess_repairs/`

| File | Purpose | Key Improvements |
| :--- | :--- | :--- |
| `v7_crop_and_search_regions.py` | Complete pipeline script | Fixed NameError, actual per-query hit cap, adaptive reruns with separate paths and deduplication, exact finite census, clean manifest generation |
| `v7_recovery_postprocess.py` | Dedicated post-crash recovery script | Reuses validated Foldseek search hits from job `53738540`, pins raw search SHA-256 hashes on IBEX, executes hit streaming, repaired saturation audit, adaptive reruns, and closure verification |
| `v7_verify_closure.py` | Independent closure verification | Fail-closed checks on all required artifacts, validates N=5,851 and R=7,773 closure, 100% CIGAR $M = qspan + tspan - alnlen$ audit, outputs `v7_closure_verification.json` |
| `v7_recovery.sbatch` | Slurm recovery job script | Allocates 16 CPUs, 64 GB RAM, 4h walltime on `batch` partition to execute `v7_recovery_postprocess.py` on IBEX |
| `test_repairs.py` | Synthetic unit test suite | 6 comprehensive unit tests verifying NameError fix, saturation detection, rerun deduplication, score monotonicity, undefined domain nulls, and fail-closed verifier |
| `REPAIRS.md` | Audit & promotion documentation | Root cause analysis, repair specifications, test receipts, and promotion instructions |

---

## 3. Synthetic Verification Results

The test suite was executed locally via `uv run python reviews/darkmatter-2026-10-13/round7/h4_postprocess_repairs/test_repairs.py`:

```
.Generated SHA-256 manifest with 7 entries at v7_results.sha256.tsv
=== Running Independent V7 Code-Data Closure Verification in /tmp/tmpobx13vx3 ===
  [EXISTS] v7_full5851_cropped_regions_manifest.tsv (59 bytes)
=== Running Independent V7 Code-Data Closure Verification in /tmp/tmpobx13vx3 ===
  [EXISTS] v7_full5851_cropped_regions_manifest.tsv (59 bytes)
  [EXISTS] v7_region_evaluations.tsv (18 bytes)
  [EXISTS] v7_model_summaries.tsv (41 bytes)
  [EXISTS] v7_stratum_census.tsv (40 bytes)
  [EXISTS] v7_figure_data.tsv (37 bytes)
  [EXISTS] v7_recovery_summary.json (43 bytes)
  [EXISTS] v7_audited_example_hits.tsv (46 bytes)
  [PASS] SHA-256 verified: v7_full5851_cropped_regions_manifest.tsv (249395878abfdbfb...)
.....
----------------------------------------------------------------------
Ran 6 tests in 0.002s

OK
```

All 6 test cases passed:
1. `test_audit_hit_cap_saturation_no_nameerror_and_categories`: Saturated vs uncapped queries, `POST_EXCLUSION_EMPTY` vs `LOW_SCORE` vs `MATCHED` tiers, no `NameError`.
2. `test_rerun_target_deduplication_and_score_monotonicity`: Deduplication of overlapping search rows; guarantees scores never worsen when higher-cap search is merged.
3. `test_undefined_domains_return_none`: Evaluates $n=0$ conditioning domains returning `None`, and boundary Wilson intervals ($k=0$ and $k=n$) avoiding $[0.0, 0.0]$.
4. `test_cigar_aligned_pairs_formula`: Validates exact operation accounting ($M = qspan + tspan - alnlen$) across matches, insertions, and deletions.
5. `test_checksum_manifest_and_verifier_fail_closed`: Verifies that `v7_results.sha256.tsv` never includes itself, and that `verify_closure` immediately raises `AssertionError` if an artifact is missing or altered.
6. `test_tier_monotonicity_validation`: Verifies that tier monotonicity assertions catch violations ($T_{k+1} > T_k$).

---

## 4. Promotion and Recovery Instructions for Root

Once active Slurm job `53738540` terminates on IBEX:

1. **Verify Raw Evidence on IBEX**:
   ```bash
   ssh ibex "ls -lh /ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p7_analyses_clean/h4_structure/results/v7_*_pdb_hits.raw.tsv"
   ```
   Both `v7_wc_pdb_hits.raw.tsv` and `v7_crop_pdb_hits.raw.tsv` will be present and populated.

2. **Promote Repaired Candidate Files to Remote Scripts Area**:
   ```bash
   rsync -avz reviews/darkmatter-2026-10-13/round7/h4_postprocess_repairs/v7_*.py \
     ibex:/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p7_analyses_clean/h4_structure/scripts/
   rsync -avz reviews/darkmatter-2026-10-13/round7/h4_postprocess_repairs/v7_recovery.sbatch \
     ibex:/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p7_analyses_clean/h4_structure/scripts/
   ```

3. **Launch Recovery Job on IBEX**:
   ```bash
   ssh ibex "sbatch /ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p7_analyses_clean/h4_structure/scripts/v7_recovery.sbatch"
   ```
   The recovery job will:
   - Verify previous search completion in `h4_v7_53738540.out`.
   - Pin raw search SHA-256 hashes to `v7_raw_search_hashes.sha256.tsv`.
   - Stream hits and validate 100% of CIGAR strings.
   - Run the repaired saturation audit and adaptive reruns (if any).
   - Write all exact census tables and figure data.
   - Run independent closure verification and generate `v7_results.sha256.tsv`.
