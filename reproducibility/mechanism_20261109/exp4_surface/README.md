# Conserved genomic neighbourhoods of archaeal dark proteins

## 1. Overview & Experimental Design

Experiment 4 evaluates whether uncharacterized ("dark") exported protein families identified in uncultivated archaeal metagenome-assembled genomes (MAGs) from the Rub' al Khali (Empty Quarter) desert exhibit conserved genomic associations with predicted surface assemblies, nutrient acquisition machinery, or protective coats.

The computational workflow consists of five stages:
1. **Curated Reference Panel**: 17 curated sequence references and 11 active single-chain PDB structures across six functional machinery guilds. Crucially, the structural panel includes deliberately competing viral and non-cellular negative controls (`REF_VIR_3J9X_SIRV2_gp26_Coat`, `REF_CTRL_8CE8_Entity1_LabelE_Autha_CcmA`, and `REF_CTRL_9UES_ChainA_Propeptide_NegativeCtrl`), rather than 17 verified cellular active structures.
2. **GSPA Evaluation**: Retained as an **exploratory and uncalibrated** heuristic. Exact checkpoint training family/order/graph leakage was not verified, and evaluation on a 3-genome pilot with 5 known labels (18 rows, 0 dark shifts) does not establish fold-level discrimination or calibrated functional separation.
3. **High-Performance Cluster Shard Production**: Parallel execution across 17 serialized shard bundles on the Shaheen III supercomputer (AMD EPYC Genoa compute nodes).
4. **Exact Circular Permutation Null & Genus Sensitivity**: Full-universe exact circular convolution null model evaluated across all $1,460 \times 6 = 8,760$ family-guild hypothesis tests, combined with a 1-MAG-per-genus-block dereplicated sensitivity test.
5. **Targeted Followup Evidence Audit & Dynamic Dossiers**: Independent member-polymer replay, byte-level MSA/tree identity validation, and authoritative contig coordinate audit for 47 candidate families.

---

## 2. Fixed Cohort & Model Availability Partition

The evaluation universe consists of **1,460 protein families** partitioned into five predeclared cohort groups:
- `candidate_route_a1_gapfill_dark_eligible`: 28 families
- `candidate_route_a2_exported_dark_eligible`: 1,190 families
- `control_known_function_exported_background_panel`: 120 families
- `control_non_exported_dark_verified_soluble_background_panel`: 120 families
- `prior_heat_stress_candidate_example`: 2 families

### Representative Model Availability Partition
The global 565 / 802 / 93 partition represents an audit of the **representative query model** per family across all 1,460 families:
- **Exact sequence hash match**: 565 families (38.7%)
- **Missing prediction**: 802 families (54.9%) — *missing precomputed models cannot be interpreted as structural absence*
- **Sequence mismatch / invalid**: 93 families (6.4%)
- **Sum**: Exactly 1,460 families across all five cohort groups.

Full member-level model audits were conducted for **428 individual members across the 47 deep-dive candidate families**, yielding:
- 326 verified complete interior members
- 32 members with exact sequence hash matches
- 27 verified complete interior members with exact sequence matches

### Unequal Model Coverage Across Cohorts
Model availability is unequal across backgrounds:
- Exported dark candidates: 42.4% exact models (504/1190)
- Known exported controls: 7.5% exact models (9/120)
- Soluble dark controls: 22.5% exact models (27/120)

Because model availability differs markedly between cohorts and backgrounds are unmatched, cohort comparisons are strictly descriptive baseline summaries; no calibrated functional negative or enrichment $p$-value is claimed.

---

## 3. Statistical Testing & Genus Sensitivity

Six machinery guilds were tested for genomic colocalization within $\pm 5$ locus windows:
1. `nutrient_acquisition_peptidase`
2. `substrate_binding_or_uptake`
3. `surface_assembly_pilin_archaellin`
4. `cell_coat_or_envelope`
5. `host_defense_antiviral`
6. `viral_or_mobile_element_origin`

### Primary Joint FDR vs Genus Dereplication
- **Primary Joint BH-FDR**: Evaluated across all 8,760 tests. Two dark candidate families achieve statistical significance under primary joint FDR ($q \le 0.05$):
  - `22PRr1::22PRr1_2531955_3` ($p = 4.71 \times 10^{-7}$, $q = 0.00206$) for `surface_assembly_pilin_archaellin`
  - `V60PRr2::V60PRr2_1314979_4` ($p = 4.00 \times 10^{-8}$, $q = 0.00035$) for `cell_coat_or_envelope`
- **Genus-Dereplicated Sensitivity Test**: Corrects for clonal expansion by restricting observations to 1 MAG per genus block. Across all 8,760 tests, the minimum genus-dereplicated joint FDR is $q = 1.0000$ (0 families significant at $q \le 0.10$).
- **Claim Boundary**: The genus sensitivity test is a post-screen sensitivity check, **not an independent replication**. Lower statistical power (due to 1-MAG-per-genus-block reduction) and taxonomic dependence are possible explanations; neither lead survives this genus sensitivity test.

---

## 4. Supported Evidence: Conserved Local Gene Modules

Authoritative contig coordinates from verified candidate loci (`inputs/task_a_stage1b_all_candidate_member_loci_verified.tsv.gz`, SHA256: `263673dc9a7213997d5ad0e97353b4b238cc3a4152ab0da68a5783fe571a6dca`) were joined to flanking evidence:
- **Lead 1 (`22PRr1::22PRr1_2531955_3`)**: 13/13 predicted role context pairs recovered across 6 genus blocks. Ten partners belong to `F37Dr2::F37Dr2_752038_3` and three to `S20PRr1::S20PRr1_1478121_8` (linking PF04964 to PilJ-like or archaellum-related reference proteins). All 13 pairs are **immediately adjacent** (relative offset $\pm 1$), **co-oriented on the same strand**, with signed intergenic gaps from **-4 to 37 bp**.
- **Lead 2 (`V60PRr2::V60PRr2_1314979_4`)**: 11/11 predicted role context pairs recovered across 5 genus blocks. Partners belong to three annotated families (6, 3, and 2 loci respectively). All 11 pairs are **immediately adjacent** (relative offset $\pm 1$), **co-oriented on the same strand**, with signed intergenic gaps from **-4 to 6 bp**.

### Scientific Interpretation & Boundaries
This tight adjacency and co-orientation supports a **conserved local gene-module hypothesis**. However:
1. **Co-transcription and operon status are unproven**: RNA expression or experimental polycistronic data are absent.
2. **No function assigned**: Neither lead has an exact structural model for its complete interior representative (`no_model_available`). Precomputed models exist only for contig-edge seed loci, which cannot support positive or negative functional conclusions for complete proteins.
3. **Whole-genome status undetermined**: Core archaeal pilin machinery (FlaI/ArlI, FlaJ/ArlJ, PibD) is unannotated within the surveyed $\pm 5$ gene neighbourhood; genome-wide status outside this local window remains undetermined.
4. **Dispositions**: Both leads are classified honestly as `uncertain_insufficient_evidence`. Exactly **0 families** are assigned confirmed biochemical, enzymatic, or motility functions.

---

## 5. Directory Layout & Manifests

```
exp4_surface/
├── README.md                                         # This documentation
├── verify_frozen_outputs.py                          # 10-gate verifier: FDR, models, dossiers, coordinates, Slurm receipts
├── code/                                             # Production pipeline, helpers, and Slurm scripts
│   ├── build_curated_reference_panel.py
│   ├── common_science.py
│   ├── exp4_followup_targeted_audit.py
│   ├── exp4_independent_polymer_and_msa_replay.py
│   ├── exp4_surface_helpers.py
│   ├── exp4_surface_pipeline.py
│   ├── sbatch_exp4_array.sh
│   ├── sbatch_exp4_collect.sh
│   ├── sbatch_exp4_followup_audit.sh
│   ├── sbatch_exp4_gspa_eval.sh
│   ├── sbatch_exp4_pilot.sh
│   ├── sbatch_exp4_stage_preflight.sh
│   └── stage3_helpers.py
├── tests/                                            # Pytest test suite (26 passed, 1 skipped)
│   └── test_exp4_surface.py
├── compact_evidence/                                 # Audited compact tables, summaries, MSAs & trees (5.3 MB)
│   ├── 22PRr1__22PRr1_2531955_3.aln.faa
│   ├── 22PRr1__22PRr1_2531955_3.iqtree.treefile
│   ├── V60PRr2__V60PRr2_1314979_4.aln.faa
│   ├── V60PRr2__V60PRr2_1314979_4.iqtree.treefile
│   ├── curated_reference_panel_exclusion_manifest.tsv
│   ├── curated_reference_panel_source_identity_audit.tsv
│   ├── curated_role_reference_panel.faa
│   ├── curated_role_reference_panel.tsv
│   ├── exp4_428_independent_polymer_replay.tsv
│   ├── exp4_47_candidate_dossiers.tsv
│   ├── exp4_47_flank_evidence.tsv
│   ├── exp4_47_member_evidence.tsv
│   ├── exp4_cohort_architecture_and_role_comparison_v2_audited.tsv
│   ├── exp4_convergent_candidates_v2.tsv
│   ├── exp4_independent_polymer_and_msa_replay_audit.json
│   ├── exp4_independent_tree_leaf_sequence_audit.json
│   ├── exp4_lead_flank_coordinates_independent.json
│   ├── exp4_lead_flank_coordinates_independent.tsv
│   ├── exp4_replay_member_identity_join.json
│   ├── exp4_root_reviewed_lead_interpretations.json
│   ├── exp4_surface_family_evaluation_v2_audited.tsv # 1,460 family evaluation table (2.8 MB)
│   ├── exp4_surface_final_summary_v2_audited.json
│   ├── exp4_tree_leaf_name_crosswalk.tsv
│   └── frozen_family_universe.tsv                    # Complete frozen 1,460 family universe & cohort definitions
└── manifests/
    ├── institutional_fetch_instructions.txt          # Research storage provenance (Shaheen / IBEX)
    ├── raw_slurm_execution_records.tsv              # 42 raw Slurm sacct records for hohndor protein jobs
    ├── remote_slurm_log_manifest.tsv                # 100 owned Shaheen log hashes & observation timestamps
    ├── sha256_manifest.tsv                          # 43 exact SHA256 checksums
    ├── slurm_execution_receipts.tsv                 # 11 authoritative Slurm job receipts (15872818..15892964)
    └── upstream_shard_bundles.sha256                # Exact SHA256 of all 17 shard bundles
```

---

## 6. How to Verify & Reproduce

### A. Lightweight Frozen-Output Verification (Local, ~2 seconds)
Does not require cluster access, GPU, or bulk databases:
```bash
cd /home/leechuck/Public/software/empty-quarter-proteome/reproducibility/mechanism_20261109/exp4_surface
uv run python verify_frozen_outputs.py
uv run pytest tests/test_exp4_surface.py
```

#### Test Suite Status & Conditional Skip Conditions:
- `uv run python verify_frozen_outputs.py`: **All 10 verification gates pass cleanly**.
- `uv run pytest tests/test_exp4_surface.py`: **26 passed, 1 skipped** (in environments with `diamond` installed, e.g. `/usr/bin/diamond`).
  - **Single actual skip condition**: `test_followup_audited_neighbourhood_panels_coverage` skips when the uncommitted, bulk neighbourhood panel (`exp4_surface_neighbourhood_panels_v2_audited.tsv`, 5,626 rows) is absent from the local lightweight checkout and remains on cluster storage.
  - **Conditional DIAMOND test**: `test_diamond_positive_fixture_execution` checks for a local `diamond` executable via `shutil.which("diamond")`. When present, it runs and passes cleanly; when absent, it skips with `"diamond binary not installed in test environment"` (yielding 25 passed, 2 skipped in minimal container environments without DIAMOND). No mock requirements are invented.

### B. Full Upstream Reanalysis (Shaheen III HPC)
Requires access to the `k10018` Shaheen project allocation, Singularity, and prerequisite input stages.

#### Historical Root Immutability Contract & Fresh-Root Reanalysis:
- Historical results under `/scratch/project/k10018/eq-mechanism-20261109/exp4_surface` are immutable, frozen evidence.
- Running reanalysis inside the historical root is strictly refused by entrypoint wrappers to protect accepted evidence and avoid false `.done` skips.
- Reanalysis is conducted from frozen shard inputs and actual upstream dependencies; exact staging is necessary.
- Entrypoint wrappers refuse pre-existing outputs and refuse silent skipping on `.done` markers.
- Path comparisons are canonicalized using `realpath -m`, rejecting trailing slashes, dot-segments, and symlinks to the historical root.

#### Step 1: Stage Package from Local Environment to Shaheen Scratch
Transfer the package code, tests, and manifests from your local environment to a dedicated, unique scratch root on Shaheen using `rsync` via the KAUST proxy jump host. Note that Slurm's `#SBATCH --output=logs/...` directive requires the `logs/` directory to exist *before* `sbatch` is invoked; explicitly create the remote scratch root and `logs/` directory prior to `rsync` and batch submission:
```bash
# Executed on local machine (lc-dell):
# Use exact fresh path with timestamp/unique suffix to prevent collisions across multiple runs:
export EQ_REMOTE_ROOT="/scratch/project/k10018/hohndor/reanalysis_$(date +%Y%m%d_%H%M%S)_exp4_surface"

# Explicitly create remote root and logs directory BEFORE rsync / sbatch:
ssh -J ws-ts shaheen "mkdir -p '${EQ_REMOTE_ROOT}/logs'"

# Stage package code, tests, and manifests:
rsync -avz -e 'ssh -J ws-ts' \
  /home/leechuck/Public/software/empty-quarter-proteome/reproducibility/mechanism_20261109/exp4_surface/ \
  shaheen:${EQ_REMOTE_ROOT}/
```

#### Step 2: Slurm Preflight (Staging & SHA256 Integrity Verification)
Shaheen policy prohibits heavy CPU and I/O operations (such as copying multi-gigabyte shard bundles or full-manifest SHA256 checks) directly on login nodes. Execute staging and hash preflight via Slurm compute allocation. Navigate to the exact fresh scratch root (avoiding wildcard `cd reanalysis_*` which fails if multiple runs exist):
```bash
ssh shaheen
# Use the exact fresh path created above:
cd "${EQ_REMOTE_ROOT}" # Or replace with exact timestamped directory path
export EQ_EXP4_ROOT="$PWD"

# Option A: Submit preflight batch job (stages inputs from historical root and verifies SHA256)
STAGE_JOB=$(sbatch --export=ALL,EQ_EXP4_ROOT="$EQ_EXP4_ROOT" code/sbatch_exp4_stage_preflight.sh | awk '{print $4}')

# Option B: Run interactively in a compute allocation
# srun -A k10018 -p workq -N 1 -n 1 -c 16 --time=00:30:00 bash code/sbatch_exp4_stage_preflight.sh
```
`sbatch_exp4_stage_preflight.sh` automatically stages the required immutable inputs from the historical root:
- `inputs/task_a_stage1b_all_candidate_member_loci_verified.tsv.gz`
- `inputs/curated_role_reference_panel.*`
- `inputs/pdb_chains/*.pdb`
- `results/shard_inputs/shard_*.pkl` (Shards 0..16)
- `provenance/frozen_family_universe.tsv`
and verifies the 17 shard bundles against the staged relative manifest:
```bash
sha256sum -c "$EQ_EXP4_ROOT/manifests/upstream_shard_bundles.sha256"
```

#### Step 3: Dependency-Chained Slurm Submission
Submit the full reanalysis pipeline chained to execute after staging preflight completes:
```bash
cd "$EQ_EXP4_ROOT"

# 1. Launch 17-task fullscreen array job across all shards (dependent on stage job)
ARRAY_JOB=$(sbatch --dependency=afterok:${STAGE_JOB} --export=ALL,EQ_EXP4_ROOT="$EQ_EXP4_ROOT" code/sbatch_exp4_array.sh | awk '{print $4}')

# 2. Launch collector stage dependent on array completion
COLLECT_JOB=$(sbatch --dependency=afterok:${ARRAY_JOB} --export=ALL,EQ_EXP4_ROOT="$EQ_EXP4_ROOT" code/sbatch_exp4_collect.sh | awk '{print $4}')

# 3. Launch targeted followup evidence audit dependent on collector completion
sbatch --dependency=afterok:${COLLECT_JOB} --export=ALL,EQ_EXP4_ROOT="$EQ_EXP4_ROOT" code/sbatch_exp4_followup_audit.sh
```
