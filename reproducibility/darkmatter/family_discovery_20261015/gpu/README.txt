================================================================================
ProstT5 AA2fold GPU Embedding Runner: Software Validation and Runtime Manifest
================================================================================
Date: 2026-10-15
Author / Agent: Gemini Flash 3.8 High (Antigravity)
Scope: Standalone reusable ProstT5 AA2fold GPU embedding infrastructure
Local Directory:  reproducibility/darkmatter/family_discovery_20261015/gpu/
ORIX Directory:   /mnt/data/u/hohndor/eq-family-discovery-20261015/gpu/

1. PURPOSE & BOUNDARIES
--------------------------------------------------------------------------------
This runner delivers a standalone, deterministic, audited embedding extraction
pipeline for caller-provided protein FASTA sequences into pooled ProstT5 AA2fold
encoder representations.

Boundaries:
- Zero manuscript modifications.
- Zero subagents or LLM model calls during execution.
- Owns strictly local reproducibility/darkmatter/family_discovery_20261015/gpu/ and
  remote ORIX /mnt/data/u/hohndor/eq-family-discovery-20261015/gpu/.
- No mutation of f12 or earlier H1-H4 structure search analyses.
- Does not compute on unpromoted EQ datasets; tested strictly on public benign
  reference enzymes (software validation only).

2. CLUSTER ASSETS & RUNTIME ENVIRONMENT
--------------------------------------------------------------------------------
- Cluster: ORIX GPU cluster (orix-login.kaust.edu.sa)
- Partition: freecycle-h200 (QoS: freecycle, requeue supported)
- GPU Allocation: --gpus=nvidia_h200:1 (NVIDIA H200 SXM5, 141 GB VRAM)
- Software Environment: Reused existing compatible pixi environment:
    /home/hohndor/pixi-envs/gspa-glm/.pixi/envs/default/bin/python
    Python: 3.10.20
    PyTorch: 2.10.0 (CUDA 12.9)
    Transformers: 4.57.6
    SentencePiece: 0.2.1
    Safetensors: 0.7.0
    NumPy: 1.26.4
- Staged Model Directory: /mnt/data/u/hohndor/eq-family-discovery-20261015/gpu/model
    Symlinked to existing verified ORIX HuggingFace cache:
    /mnt/data/u/hohndor/hf_cache/hub/models--Rostlab--ProstT5/
    - config.json (SHA256: 69c1915585d20b1d3c1d44a566968a27c7d9cf784f223569ffc5dc2a9e0f892c)
    - model.safetensors (SHA256: e067227fca56e6e60423dabeb5d911de2161d5c078afeb2bce563d9a53a18dde)
    - spiece.model (SHA256: 74da7b4afcde53faa570114b530c726135bdfcdb813dec3abfb27f9d44db7324)
    - tokenizer_config.json, special_tokens_map.json, added_tokens.json
  *Note: Model weights SHA256 was computed on the worker node inside the Slurm
   allocation as a short provenance step (per compliance policy, avoiding login parsing).

3. RUNNER IMPLEMENTATION DETAILS (prostt5_embed.py)
--------------------------------------------------------------------------------
- Input & IDs: Preserves exact input FASTA sequence IDs and original order.
- Hashes & Lengths: Computes per-sequence SHA256 and residue lengths.
- Ambiguity Handling:
  * Standard 20 amino acids + X permitted.
  * UZOB ambiguity codes handled per Rostlab ProstT5 documentation:
    U (Sec) -> X, Z (Glx) -> X, O (Pyl) -> X, B (Asx) -> X.
  * Every UZOB substitution is logged per-sequence with exact residue position,
    original code, and replacement in sequence_manifest.tsv.
  * Trailing stop codon ('*') is stripped and logged (trailing_stop_stripped=True).
  * Non-standard characters or internal stops fail closed with informative errors.
- Long-Sequence Policy: Explicit --max-seq-len limit (default 2048). Sequences
  exceeding the limit raise an error. No silent truncation is permitted.
- Tokenization: Prepends "<AA2fold> " and separates residues with spaces.
- Mean Pooling with Special-Token Exclusion:
  * Token 0: <AA2fold> (prefix)
  * Tokens 1 to 1 + L: amino acid residues (length L)
  * Token 1 + L: </s> (EOS)
  * Tokens > 1 + L: <pad> (padding)
  * Pooling averages ONLY residue hidden states across dimensions 1 to 1 + L.
  * Special tokens (<AA2fold>, </s>, <pad>) are strictly excluded.
- Length-Aware Batching: Sorts sequences by length into token-budgeted batches
  (--batch-tokens, default 4096), then restores exact input ID order.
- Deterministic Inference: eval() mode, torch.inference_mode(), fixed seeds (42),
  cudnn.deterministic = True.
- Idempotence / Resume: Computes SHA256 run signature of (fasta_sha256, seq_count,
  model_sha256, precision, max_seq_len, batch_tokens). Automatically skips
  recomputation if matching provenance and valid checksums already exist.
- Atomic Publication: Writes to temporary workspace (.work_<ts>_<uuid>), verifies
  checksums, then atomically publishes to output destination.
- Emitted Artifacts:
  * embeddings.npz: compressed archive with embeddings, ids, lengths, sequence_sha256
  * embeddings.npy: raw float32 numpy matrix (N, 1024)
  * sequence_manifest.tsv: sequence metadata, UZOB transformations, lengths, SHAs
  * provenance.json: hardware, software versions, peak VRAM, runtime, signature
  * checksums.sha256: SHA256 checksums of all published files

4. WORKER EXECUTION GUARD (require_slurm_worker.sh)
--------------------------------------------------------------------------------
Copied from existing verified repository guard:
- Automatically inspects hostname, SLURM_JOB_ID, and scontrol job record.
- Refuses execution on login hosts (exit code 64).
- Requires allocation to be in RUNNING state on an assigned compute worker.
- Both smoke_test.sbatch and run_embed.sbatch enforce this guard at entry.

5. SMOKE TEST VALIDATION RECEIPTS
--------------------------------------------------------------------------------
Original Smoke Test (Immutable):
- Slurm Job ID: 316765 (Node: orix-worker-h200-7)
- Idempotency Test Job ID: 316782 (exited in 3s)
- Preserved untouched in: smoke_results/

Smoke Test v2 (Immutable):
- Slurm Job ID: 316783 (Node: orix-worker-h200-6)
- Idempotency Test Job ID: 316785 (exited in 3s)
- Preserved untouched in: smoke_results_v2/

Production-Ready Smoke Test v3 (Argon Preflight Validated):
- Slurm Job ID: 316812 (Node: orix-worker-h200-6)
- Runner-Level Signature Re-check Job ID: 316813 (Node: orix-worker-h200-7, verified & skipped)
- Target Directory: smoke_results_v3/
- Runner Code SHA256: d5e019d7e1fc2a22ce55668b8c36f73932607082f6fa53b2cb385f9c29729a81
- Comprehensive Run Signature: d1bfe7b6b3306af619664566a6c66fc37013903bd7df7dde30c52973096f3398
- Key Fixes Delivered:
  1. Unified Raw-Sequence Hash Semantics:
     * raw_sha256: computed on exact raw FASTA sequence (pre stop-strip).
     * clean_sha256: computed on sequence after stripping trailing '*'.
     * length: clean residue length (137).
     * Verified identical between runner and independent verifier.
  2. Real Trailing-Stop Input Tested:
     * Added O52251_stop (sequence ending in '*').
     * Runner strips trailing stop, logs trailing_stop_stripped=True.
     * Duplicate check confirms O52251 vs O52251_stop embedding equivalence: max_abs_diff = 0.0.
  3. Removed sbatch COMPLETE Short-Circuit:
     * run_embed.sbatch and smoke_test_v3.sbatch always invoke prostt5_embed.py directly.
     * The runner evaluates the comprehensive signature (code, model, tokenizer, libraries, input).
     * Job 316813 verified that the runner directly evaluates signature and skips in seconds.
  4. Preserved Accurate Hardware & Metrics:
     * Hardware: Genuine NVIDIA H200 (141 GB SXM5, driver 570.211.01, CUDA 12.9).
     * Peak VRAM: 3.26 GB.
     * Min cosine similarity (single vs batched): 0.99999976.
     * Max absolute difference (single vs batched): 4.8828e-04 (well within 5e-03 tolerance).
     * Duplicate equivalence (exact and trailing-stop): 0.0 (strictly identical).

Test cohort (7 public benign characterised enzyme reference sequences):
1. O52251: Halomonas elongata L-ectoine synthase (EctC), 137 aa
2. Q9ZEU6: Chromohalobacter israelensis L-ectoine synthase (EctC), 130 aa
3. P19406: Bacillus subtilis alkaline phosphatase 4 (PhoA), 461 aa
4. P42251: Bacillus subtilis alkaline phosphatase D (PhoD), 583 aa
5. O52251_dup: Exact duplicate of O52251 (verifies duplicate embedding equivalence)
6. O52251_uzob: Controlled variant of O52251 with U, Z, O, B substitutions
7. O52251_stop: Exact copy of O52251 with trailing stop codon '*' (verifies stop strip & hash)

Validation Results v3 (verify_embeddings.py):
- Files present & checksums verified: PASSED (100% match)
- npz and npy identical: PASSED (np.array_equal == True)
- Output dimensions: [7, 1024] (float32): PASSED
- Finiteness: PASSED (all finite, 0 NaNs, 0 Infs)
- ID order preserved against FASTA: PASSED (all 7 IDs match in exact order)
- Raw and clean sequence SHAs match: PASSED
- Sequence manifest integrity: PASSED
- Single-vs-batched consistency:
  * Max absolute difference: 4.8828e-04 (well within half precision tolerance 5e-03)
  * Mean absolute difference: 2.2704e-05
  * Min cosine similarity: 0.99999976 (> 0.9999 threshold)
  * Status: PASSED
- Duplicate sequence equivalence:
  * O52251 vs O52251_dup max abs diff: 0.0 (strictly identical)
  * O52251 vs O52251_stop max abs diff: 0.0 (strictly identical)
  * Status: PASSED
- UZOB ambiguity transformation:
  * 4 substitutions (pos 9 U->X, pos 19 Z->X, pos 29 O->X, pos 39 B->X)
  * Transformed and recorded in sequence_manifest.tsv: PASSED
- GPU utilisation and memory:
  * Device: NVIDIA H200 (141 GB SXM5)
  * Peak VRAM allocated: 3.2556 GB
  * Continuous 1-second utilization logged to gpu_utilization.csv
  * Status: PASSED

6. HOW ROOT CAN INVOKE FOR PROMOTED DATASETS
--------------------------------------------------------------------------------
Option A: Via Slurm sbatch (Recommended)
  sbatch \
    --export=ALL,INPUT_FASTA=/path/to/promoted.faa,OUTPUT_DIR=/path/to/output_dir \
    /mnt/data/u/hohndor/eq-family-discovery-20261015/gpu/run_embed.sbatch

Option B: Inside an active Slurm GPU allocation
  export ROOT=/mnt/data/u/hohndor/eq-family-discovery-20261015/gpu
  $ROOT/require_slurm_worker.sh \
    /home/hohndor/pixi-envs/gspa-glm/.pixi/envs/default/bin/python \
    $ROOT/prostt5_embed.py \
    --fasta /path/to/promoted.faa \
    --output-dir /path/to/output_dir \
    --model-dir $ROOT/model \
    --precision half \
    --batch-tokens 4096

Audit & Verification Command:
  /home/hohndor/pixi-envs/gspa-glm/.pixi/envs/default/bin/python \
    /mnt/data/u/hohndor/eq-family-discovery-20261015/gpu/verify_embeddings.py \
    --output-dir /path/to/output_dir \
    --fasta /path/to/promoted.faa \
    --report /path/to/output_dir/verification_report.json

7. CPU FAMILY EMBEDDING ANALYSIS & DIAGNOSTIC PLOTTING
--------------------------------------------------------------------------------
Script: family_embedding_analysis.py
Slurm Batch Runner: run_analysis_cpu.sbatch
IBEX Verification Job: 53816470 (Worker cn511-18, partition batch)
Python Environment on IBEX: /ibex/scratch/projects/c2014/rob/envs/gspa-neural/bin/python
  Libraries: numpy 2.4.3, scipy 1.17.1, scikit-learn 1.8.0, matplotlib 3.10.9

Key Rules & Capabilities:
- CLI Inputs:
    --embeddings <path.npz> : Output from ProstT5 GPU runner (with embeddings, ids, lengths, raw/clean SHAs)
    --metadata <path.tsv>   : Metadata table with unique IDs, source, family_class, length, phylum, genus
    --output-dir <path>     : Output directory for coordinates, distances, diagnostic plots, and provenance
- Exact ID Closure & SHA Verification:
    Strict 1-to-1 bijection enforced between --embeddings and --metadata.
    Verifies clean and raw sequence SHA256 matches and sequence lengths match.
- Unrelated Folds/Classes Strictly Separated:
    Every family_class is analyzed and plotted independently (e.g. EctC vs Phosphatase; F5beta vs F5gamma).
    Unrelated folds are NEVER mixed into a shared manifold projection.
- Multi-Seed t-SNE & PCA Baseline:
    Evaluates 3 random seeds (42, 123, 20261015) using PCA initialization with seed.
    No cherry-picked 'prettiest' seed: all 3 seeds are plotted side-by-side with PCA.
    Valid perplexities evaluated subject to N. Exact method used for N < 500 to avoid quad-tree stalls.
    Graph trustworthiness (sklearn.manifold.trustworthiness) recorded per seed.
- High-Dimensional Cosine Nearest-Reference Distance:
    Computed directly in native 1024-d ProstT5 space.
    Queries: minimum cosine distance to curated reference in same family.
    References: leave-one-out minimum cosine distance to other curated reference in same family.
- Distance vs Length Disparity Diagnostic:
    Plots cosine distance vs |length(seq) - length(nearest_ref)| to detect fragment/truncation artifacts.
- Exact Duplicate Sequence Handling:
    Exact duplicate clusters identified by sequence SHA; duplicate counts recorded and point sizes scaled.
- Scientific Caveat Policy:
    Prominently annotated on all figures and provenance:
    "Descriptive cosine distances do NOT prove functional novelty or classification.
     Bounded by reference database density, length disparities, and taxonomy sampling."

Execution Command on IBEX:
  sbatch --chdir=/home/hohndor/eq-family-discovery-20261015/cpu run_analysis_cpu.sbatch

8. LOCAL FILE INVENTORY
--------------------------------------------------------------------------------
reproducibility/darkmatter/family_discovery_20261015/gpu/
├── README.txt                     # This document
├── STATUS.json                    # Execution status and verification metadata
├── environment_manifest.json      # Cluster, model, and library environment specification
├── prostt5_embed.py               # Standalone ProstT5 GPU runner (v3)
├── verify_embeddings.py           # Independent audit and verification script (v3)
├── require_slurm_worker.sh        # Compute allocation and login host guard
├── run_embed.sbatch               # Reusable Slurm sbatch script for arbitrary FASTAs
├── family_embedding_analysis.py   # Standalone reproducible CPU analysis and diagnostic plotting script
├── run_analysis_cpu.sbatch        # IBEX Slurm CPU runner for family embedding analysis
├── smoke_test_metadata.tsv        # Smoke test metadata with sources, family classes, taxonomy, and SHAs
├── smoke_test.sbatch              # Slurm sbatch script for original smoke test
├── smoke_test_v2.sbatch           # Slurm sbatch script for smoke test v2
├── smoke_test_v3.sbatch           # Slurm sbatch script for smoke test v3
├── smoke_references.fasta         # Benign enzyme test sequences (incl. trailing stop)
├── smoke_references_manifest.tsv  # Accessions, entry names, and sequence SHA256 table
├── logs/                          # Slurm stdout and stderr receipts (jobs 316765, 316782, 316783, 316785, 316812, 316813)
├── smoke_results/                 # Original verified smoke test outputs (immutable)
├── smoke_results_v2/              # V2 verified smoke outputs (immutable)
├── smoke_results_v3/              # Production-ready v3 verified smoke outputs and provenance
│   ├── embeddings.npz             # Pooled embeddings (7, 1024) with IDs, lengths, raw/clean SHAs
│   ├── embeddings.npy             # Raw float32 embeddings matrix
│   ├── sequence_manifest.tsv      # Per-sequence metadata, UZOB transformations, and raw/clean SHAs
│   ├── provenance.json            # Machine-readable hardware/software/runtime provenance (v3)
│   ├── verification_report.json   # Independent audit report (ALL CHECKS PASSED)
│   ├── gpu_info.csv               # GPU device name, UUID, driver version
│   ├── gpu_utilization.csv        # Time-series utilization and memory consumption
│   ├── model_cache.sha256         # Model cache checksums computed on worker
│   ├── slurm_execution.txt        # Job ID, worker hostname, timestamp
│   ├── checksums.sha256           # SHA256 of published artifacts
│   ├── final_bundle.sha256        # Comprehensive package bundle checksum
│   ├── COMPLETE                   # Idempotence completion marker
│   └── STARTED                    # Execution start marker
└── cpu_smoke_results/             # Verified IBEX CPU smoke test receipts (Job 53816470, cn511-18)
    ├── analysis_cpu_53816470.out  # Slurm worker stdout log
    ├── analysis_cpu_53816470.err  # Slurm worker stderr log
    ├── analysis_provenance.json   # Full software/library provenance, parameters, seeds, trustworthiness
    ├── checksums.sha256           # SHA256 checksums of emitted tables and figures
    ├── COMPLETE                   # Completion marker
    ├── coordinates.tsv            # PCA and multi-seed t-SNE coordinates table
    ├── distances.tsv              # High-dimensional cosine nearest-reference distances and duplicate metadata
    ├── slurm_execution.txt        # Slurm job ID, worker hostname, timestamp, Python path
    ├── tsne_multiseed_EctC_family_perp2.0.pdf / .png          # Multi-seed manifold comparison plot (EctC)
    ├── distance_diagnostics_EctC_family.pdf / .png            # Distance distribution & length delta diagnostic (EctC)
    ├── tsne_multiseed_Phosphatase_family_perp1.0.pdf / .png   # Multi-seed manifold comparison plot (Phosphatase)
    └── distance_diagnostics_Phosphatase_family.pdf / .png     # Distance distribution & length delta diagnostic (Phosphatase)
================================================================================
