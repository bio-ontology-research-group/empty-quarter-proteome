# Hypothesis 4 (Structure) Post-Processing Repairs & Recovery Plan (Root-Audited Revision)

**Target Manuscript**: Nature Communications (Dark-Matter MAGs Proteome)  
**Location**: `reviews/darkmatter-2026-10-13/round7/h4_postprocess_repairs/`  
**Date**: 2026-10-14  
**Current Run Status**: Slurm Job `53738540` is currently **RUNNING** on IBEX (`cn509-04-l`, 32 CPUs, 128 GB RAM). It has NOT yet failed or succeeded. Candidate repairs have been prepared in isolation to address all known defects upon search completion without disrupting active compute. Initial candidates are preserved in `initial_candidates_v1/` for full auditability.

---

## 1. Root Audit Findings & Implemented Repairs

### Item 1: Strict Log Anchoring, Terminal State Guard, & Baseline Reference Contract
- **Problem**: Recovery prerequisites previously failed open if the log was missing, and broad substring matching (`"Completed"` + filename) would mistakenly accept a completed whole-chain search while the cropped search was merely started or incomplete.
- **Implemented Fix**:
  1. `validate_slurm_log_state()` strictly requires an existing, non-empty log file.
  2. Uses anchored regex for each raw output:
     - Whole-chain: `^Completed <full wc_hits_path> in \d+\.\d+s \(\d+\.\d+m\)$`
     - Cropped regions: `^Completed <full crop_hits_path> in \d+\.\d+s \(\d+\.\d+m\)$`
  3. Uses anchored regex for the full row validation line:
     - `^Validated \d+ whole-chain hit rows and \d+ cropped hit rows\.$`
  4. Anchored check for Step 5 header:
     - `^=== Step 5: Auditing Hit-Cap Saturation & Checking for Exclusion-Empty Capped Queries ===$`
  5. Terminal State Guard:
     - Checks terminal execution state: explicitly requires EITHER the expected `NameError: name 'a' is not defined` traceback OR explicit success `"All V7 outputs successfully written and verified in"`.
     - Refuses to execute if the job is still actively running or if it terminated in an unknown/unexpected error state (e.g., OOM, timeout, segfault).
  6. Emits structured `v7_recovery_job_receipt.json` recording log path, log SHA-256 digest, terminal state, pinned raw search evidence hashes, and verification timestamp.
  7. Validates baseline reference identity: checks `per_rep_structures.tsv` (1,436,582 bytes, `9d2aa6f8...`), `per_control_structures.tsv` (250,513 bytes, `33869ed2...`), and `current_reference_control_manifest_20260802.tsv` (1,502,037 bytes, `a65ba1b0...`) against pinned byte sizes and SHA-256 hashes matching `INPUTS.tsv`.

---

### Item 2: Atomic Foldseek Completion Signatures & Reuse Validation
- **Problem**: `run_foldseek_search(..., reuse=True)` previously checked only `os.path.exists() and os.path.getsize() > 0`, blindly trusting arbitrary non-empty or interrupted files.
- **Implemented Fix**:
  1. Implemented `run_foldseek_search_atomic()`: writes an atomic companion file `{out_tsv}.signature.json` containing:
     - Subprocess returncode (`0`) and completion status (`"COMPLETED"`).
     - Raw TSV SHA-256 hash and exact byte size.
     - Query count and exact SHA-256 hashes for every query PDB in `query_dir`.
     - DB component metadata (`foldseek_pdb.version` hash).
     - Search settings (`--alignment-type 1`, `-e 10`, `-s 9.5`, `--max-seqs`, `--threads`).
  2. When `reuse=True`: validates that `{out_tsv}.signature.json` exists, verifies raw TSV SHA-256 against the signature, verifies settings match, verifies returncode 0, and validates query hashes. Refuses unsigned, partial, or mismatched existing results.
  3. Validates that every requested query source file exists before creating rerun symlinks.
  4. Never performs destructive overwriting of raw evidence.

---

### Item 3: Full-Row Deduplication, Score Monotonicity, & Standalone Rerun Saturation
- **Problem**: Deduplicating by `(target, qstart, qend, tstart, tend)` would discard a higher-scoring alignment from the rerun if the coordinate span was identical. Additionally, saturation was evaluated on the union count rather than the standalone rerun search count.
- **Implemented Fix**:
  1. Full-row exact duplicate removal: two rows are deduplicated ONLY if identical across all 21 Foldseek TSV columns (`row_key = tuple(p[:21])`). Alignments with identical coordinates but distinct scores/CIGAR strings are deterministically retained and evaluated.
  2. Preserves best scores independently at every exclusion tier:
     `merged_top_qtm[k]` and `merged_top_alntm[k]` select the best alignment across original and rerun hits.
     Score monotonicity is verified: asserts `final_top_qtm >= initial_top_qtm - 1e-9`.
  3. Standalone Rerun Saturation: evaluates saturation after rerun using the standalone rerun hit count (`standalone_rerun_hits_raw = len(rerun_lines)`), checking `standalone_rerun_hits_raw >= higher_max_seqs` (not union count).
  4. Final Saturation Audit: records both initial saturation audit and final post-rerun saturation audit in `summary_json`.

---

### Item 4: Uniform Return Dictionary Keys
- **Problem**: When `prob_qids` was empty, `rerun_problematic_queries` returned a dictionary lacking `queries_updated` and `queries_rescued_to_matched`, causing a `KeyError` when only one side (whole-chain or crop) triggered a rerun.
- **Implemented Fix**:
  1. Both branches (rerun executed vs no rerun needed) return an identical dictionary structure:
     `{"rerun_executed": bool, "query_type": str, "problematic_queries_count": int, "queries_updated": int, "queries_rescued_to_matched": int, "records": list}`.
  2. Verified by unit test `test_uniform_rerun_return_dict_keys`.

---

### Item 5: Two-Phase Manifest & Strict Closure Verification
- **Problem**: `v7_verify_closure.py` previously passed when `v7_results.sha256.tsv` was missing. Generating the manifest while verifying caused circularity with `v7_closure_verification.json`.
- **Implemented Fix**:
  1. Standalone Execution Mode: defaults to `mode="strict"`. Fails closed immediately (`AssertionError`) if `v7_results.sha256.tsv` is missing, or if any required file (including `v7_closure_verification.json` and `v7_raw_search_hashes.sha256.tsv`) is missing or has a mismatched checksum.
  2. Two-Phase Pipeline Execution:
     - Phase 1: Generates provisional manifest `v7_results.sha256.tsv` over data files (excluding receipt).
     - Phase 2: Runs `verify_closure(OUT_DIR, mode="provisional")`: verifies data tables against provisional manifest, runs scientific closure assertions, and writes `v7_closure_verification.json`.
     - Phase 3: Generates final manifest `v7_results.sha256.tsv` including all results, `v7_closure_verification.json`, and `v7_raw_search_hashes.sha256.tsv` (explicitly excluding itself).
     - Phase 4: Runs `verify_closure(OUT_DIR, mode="strict")`: strictly validates all entries in final manifest against disk without rewriting receipt.

---

### Item 6: Real Failure-Path Synthetic Tests
- **Test Suite**: `test_repairs.py` includes 7 comprehensive test fixtures testing:
  1. Log validation failure paths: missing log, WC-only completed / crop started, missing validation line, job still running, unexpected terminal state, expected NameError crash, and explicit success.
  2. Baseline reference contract: exact byte sizes and SHA-256 digests.
  3. Full-row deduplication: ensures alignments with same coordinates but different scores retain the higher score; exact duplicate rows removed; score monotonicity preserved.
  4. Standalone rerun saturation evaluation: 5,000 hits -> saturated at 5,000; 3,200 hits -> uncapped at 5,000.
  5. Uniform return keys: zero `KeyError` when only one side reruns.
  6. Strict vs provisional manifest verification: fails closed on missing manifest, missing receipt, or checksum mismatch.
  7. Undefined domains (`None`/`null` instead of `[0.0, 0.0]`) and exact CIGAR $M = qspan + tspan - alnlen$ accounting.
- **Verification Output**: All 7 tests passed in 0.004s via `uv run python test_repairs.py`.

---

## 2. Inventory of Candidate Files

All files are located in `reviews/darkmatter-2026-10-13/round7/h4_postprocess_repairs/`:

| File | Purpose | Key Enhancements |
| :--- | :--- | :--- |
| `v7_crop_and_search_regions.py` | Complete analysis pipeline | Atomic search signatures, full-row deduplication, standalone rerun saturation, baseline reference checks, two-phase manifest protocol |
| `v7_recovery_postprocess.py` | Dedicated post-crash recovery wrapper | Anchored log validation, terminal state guard, pins raw evidence SHA-256, baseline verification, full hit streaming, repaired saturation audit & reruns |
| `v7_verify_closure.py` | Independent closure verifier | Default strict mode fails closed on missing manifest/receipt; provisional mode writes receipt; verifies N=5,851 / R=7,773 / 161 no-reg closure |
| `v7_recovery.sbatch` | Slurm batch recovery script | Allocates 16 CPUs, 64 GB RAM, 4h walltime on `batch` partition to execute recovery post-processing and strict verification on IBEX |
| `test_repairs.py` | Synthetic unit test suite | 7 unit tests covering real failure paths, deduplication, score monotonicity, log validation, and strict manifest verification |
| `REPAIRS.md` | Audit & promotion documentation | Comprehensive record of defects, fixes, test outputs, and execution protocol |
| `initial_candidates_v1/` | Initial candidate snapshot | Preserves pre-review candidate versions for complete diffability and auditability |

---

## 3. Promotion and Execution Protocol for Root

Once active Slurm job `53738540` finishes the search phase:

1. **Inspect Search Log**:
   ```bash
   ssh ibex "tail -n 25 /ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p7_analyses_clean/h4_structure/logs/h4_v7_53738540.out"
   ```
   Confirm that both Foldseek searches finished cleanly and the job terminated with the expected `NameError: name 'a' is not defined`.

2. **Promote Repaired Scripts**:
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
   - Verify anchored search completion lines and terminal state in `h4_v7_53738540.out`.
   - Validate baseline reference files against pinned hashes and sizes.
   - Pin raw search SHA-256 hashes to `v7_raw_search_hashes.sha256.tsv`.
   - Stream hits (100% CIGAR validation) and run repaired saturation audit.
   - Run adaptive reruns with atomic signatures and full-row deduplication if needed.
   - Write exact census tables, figure data, and recovery summary.
   - Run two-phase closure verification and seal `v7_results.sha256.tsv`.
