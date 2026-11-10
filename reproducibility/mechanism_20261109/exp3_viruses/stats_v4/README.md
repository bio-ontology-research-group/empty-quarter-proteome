# Stage 4 v4 Viral Statistical Audit & Independent 47-Site Abundance Replay (`exp3_viruses/stats_v4`)

## 1. Scope & Ownership Boundaries

This subdirectory (`reproducibility/mechanism_20261109/exp3_viruses/stats_v4`) packages the completed **Stage 4 v4 viral statistical audit** (Shaheen III Slurm job `15894468`) and its **compact source-abundance & 47-site independent audit extraction** (Shaheen III Slurm job `15894666`).

- **Strict Subsection Ownership**: This subsection owns **only** the Stage 4 v4 statistical audit (`code/stage4_statistical_audit_v4.py`, `code/sbatch_stage4_statistical_audit_v4.sh`, `tests/test_stage4_statistical_audit_v4.py`, `results/stage4_v4_*`, compact abundance subsets, and 47-site raw independent audit tables). Competitive read mapping, combined reference index validation, and sequence hallmark followups are owned separately by the sibling viral mapping workflow.
- **Preservation of Stage 4 v3**: Earlier Stage 4 v3 outputs remain untouched in the working branch as a documented invalid-inference baseline (exploratory descriptives only; superseded by v4).
- **Exploratory Scientific Status**: All v4 statistical evaluations are exploratory, conditional on CRISPR spacer detection and protein-family assay coverage, and pending manuscript integration.

---

## 2. Root-Reviewed Candidate Pair Ledger & Cohort Partition

### A. Graded Candidate Pair Ledger (`results/stage4_v4_candidate_pairs_graded_ledger.tsv`)
From the 2,121 raw rows in `stage3_archaeal_virus_host_pairs_all_tiers_v2.tsv` (SHA256: `43cca6029736b7334d0755c5cdf97ca969994d064d570208cf1c5c6f52b5fbca`), 8 raw rows carry `Tier1_verified_direct`. Without applying any arbitrary global CheckV host-gene threshold, exact `(canonical_591_host_mag, target_seq_sha256)` deduplication and individual root dispositions partition these 8 rows as follows:
- **1 Alias Duplicate Collapsed**: `S49Dr1__self_SemiBin_1106 -> coasm_48_S::k127_1210458` (exact sequence SHA256 `ca78a3b6...` identical to representative `48Sr2::48Sr2_190170`).
- **1 Individual Unresolved Root Disposition Excluded**: `V45Dr1__self_SemiBin_1192_sub -> 41Sr3::41Sr3_56779` (`results/stage4_v4_root_disposition_41Sr3.json`, SHA256: `c36c9e8f...`; 196,115 bp contig with 179 ORFs, 70 host genes, 1 viral gene, sole Pin recombinase hallmark, no delimited viral interval; retained as unresolved and excluded from both primary positives and negative controls).
- **1 Boundary-Sensitive Provisional Candidate Excluded from Primary**: `43PRr3__self_SemiBin_469 -> 44Dr2::44Dr2_4002037` (`results/stage4_v4_provisional_disposition_second_thermo.json`, SHA256: `6b5580ae...`; canonical host taxonomy `p__Thermoplasmatota / o__JACQPN01 / g__` blank, completeness `91.62%`, contamination `4.85%`; exact 32-bp spacer `archspc_001539` matches target positions `37,670..37,701`, outside the historical catalog viral core `1..34,828` but inside the fresh CheckV viral interval `1..43,323` on the 67,924-bp contig; evaluated in `provisional_sensitivity` only).
- **5 Primary Detected-Link Pairs across 4 Canonical Archaeal Hosts** (`results/stage4_v4_primary_detected_link_pairs.tsv`):
  1. `41PRr2__single_SemiBin_236` (`p__Thermoplasmatota`, `o__JACPAO01`, `g__` unassigned) $\to$ `V60PRr2::V60PRr2_715976` (SHA256: `6b19afac...`)
  2. `41PRr2__single_SemiBin_97` (`p__Halobacteriota`, `g__Halobacterium`) $\to$ `55Dr2::55Dr2_75333` (SHA256: `2ce44dd7...`)
  3. `41PRr2__single_SemiBin_97` (`p__Halobacteriota`, `g__Halobacterium`) $\to$ `9Sr2::9Sr2_284237` (SHA256: `96ff8fd3...`)
  4. `S49Dr1__self_SemiBin_1106` (`p__Halobacteriota`, `g__CAZOFG01`) $\to$ `48Sr2::48Sr2_190170` (SHA256: `ca78a3b6...`)
  5. `V48Dr1__self_SemiBin_3215` (`p__Halobacteriota`, `g__Halorussus`) $\to$ `coasm_48_PR::k127_4172870` (SHA256: `539fdb93...`)
- **6 Provisional Sensitivity Pairs across 5 Canonical Archaeal Hosts** (`results/stage4_v4_provisional_sensitivity_pairs.tsv`): the 5 primary pairs plus `43PRr3__self_SemiBin_469 -> 44Dr2::44Dr2_4002037`.

### B. Detection Opportunity Partition (`results/stage4_v4_detection_opportunity_partition.tsv`)
Across all **591 canonical archaeal MAGs**:
- **447 Zero-Spacer Untested Genomes** (`n_valid_spacers == 0`): Classified as `zero_spacers_untested_not_negative` with `primary_DETECTEDLINK = NA` and `provisional_DETECTEDLINK = NA`. Zero recovered spacers reflect absence of CRISPR-array detection opportunity, never biological absence of viral infection.
- **144 Screen-Positive Eligible Genomes** (`minced_ok & cctyper_ok & n_valid_spacers > 0`): Form the conditional detection cohort (`4` primary detected-link hosts, `5` provisional detected-link hosts).

---

## 3. Complete 237-Genome Protein-Family Assignment Scope & Design Matrix Audit

### A. Independently Verified Complete Family Assignment Scope
- Validated against `inputs/stage1_summary.json` (SHA256: `9eafab23...`) and `inputs/stage1b_summary.json` (SHA256: `4f69536f...`):
  - **237 Empty Quarter Thermoplasmatota genomes** (`777,177` proteins in Stage 1 scope), matching `inputs/taxonomy_reconciliation_thermo237.tsv` 1:1 in exact genome IDs.
  - **41,494 candidate member locus rows** across **8,776 exact family IDs** in `inputs/task_a_stage1b_all_candidate_member_loci_verified.tsv.gz` (SHA256: `263673dc...`).
  - All **3,184 rows with `final_promotion_gate_pass == False`** are retained as membership evidence so that non-promoted loci are never misclassified as protein-family absences.
- Across the **47 frozen candidate protein families** (`inputs/exp4_47_member_evidence.tsv`, SHA256: `5a1fc798...`):
  - **44 families** have exact family IDs in the complete 237-genome assignment table, yielding **58 assay-covered screen-positive eligible Thermoplasmatota genomes** (containing **1 primary detected-link host** and **2 provisional detected-link hosts**).
  - **3 families** (`22PRr3::22PRr3_756560_4`, `43Sr2::43Sr2_645388_1`, `coasm_24_R::k127_1412270_2`) are dossier-only without a verified family-specific absence denominator and remain `NA` (`non_estimable_no_verified_absence_denominator`).
  - Benjamini-Hochberg adjustment uses the fixed hypothesis family size **$m = 47$** (`bh_fdr_qvalue_fixed47`). No family reaches significance (`min q = 1.0`).

### B. Defense/Surface Design Matrix Audit (`results/stage4_v4_design_matrix_audit.json`)
- **Eligible Primary Cohort (`n = 144` raw rows)**:
  - Complete-case rows without completeness or membership imputation: **58 genomes** (all within `p__Thermoplasmatota`), containing **1 primary `DETECTEDLINK` event** (`EPV = 0.125` across 8 slope predictors).
  - Algebraic rank: **5 of 9 parameters** (`rank_deficient = True`; 4 non-reference phylum indicator columns are constant zeros on the 58-genome complete-case subset).
  - Linear programming separation diagnostic: **Quasi-complete separation witness confirmed** (`status = "quasi_complete_separation_LP_witness"`).
  - Both design leads (`22PRr1::22PRr1_2531955_3` and `V60PRr2::V60PRr2_1314979_4`) exhibit `quasi_complete_separation` (`present_positive = 0`).
  - No overparameterized logistic model is fit (`model_fit_performed = False`).

---

## 4. Environment-Matched Abundance Evaluation & Independent 47-Site Audit Data

### A. Metadata-Only Environment Selection (`results/stage4_v4_environment_design.json`)
- Selected **before inspecting any abundance outcomes** by choosing the `(trip, compartment)` stratum with the largest number of distinct sites among the 441 joint metadata/matrix-available samples joined to the current audited Rund dune-field map (`inputs/matched_profiles.tsv`, SHA256: `338d1530...`, where `RootAdjacent` is the audited synonym of `PlantRoot`):
  - **Selected Stratum**: `trip5` $\times$ `PlantRoot` (**47 samples across 47 distinct sites**, 1 sample per site).
  - **Audited Dune-Field Strata (4 permutable strata, 47 permutable sites, 0 fixed singletons)**:
    - `uruq`: 28 sites
    - `seif`: 11 sites
    - `megabarchans`: 6 sites
    - `sigmoidal and hooked`: 2 sites
  - **Permutation Space**: $\log_{10}(\prod_f n_f!) = 40.244$ ($>10^{40}$ within-field site permutations); evaluated with **9,999 Monte Carlo within-field viral permutations** (`seed = 20261109`, resolution $10^{-4}$).

### B. Candidate & Matched Alternative MAG Results (`results/stage4_v4_environment_matched_abundance.tsv`)
- **Focal Thermoplasmatota Pair (`41PRr2__single_SemiBin_236 -> V60PRr2::V60PRr2_715976`)**:
  - Observed across 47 sites (`1` host-present site, `32` virus-present sites, `1` co-detected site): $\rho = 0.2542$, within-field permutation $p = 0.1692$, fixed-family BH $q = 0.4230$ (primary $m=5$) / $q = 0.5076$ (provisional $m=6$). Unblocked sensitivity $p = 0.0222$ ($q = 0.0555$) is hypothesis-generating only.
  - Matched alternative MAG `60PRr1__single_SemiBin_79` (`p__Thermoplasmatota`, `g__` unassigned; `unknown_genus_relatedness_not_established`) has 0 present sites in `trip5` $\times$ `PlantRoot` (`non_estimable_zero_variance_or_insufficient_sites`).
- **Halobacteriota Pair (`41PRr2__single_SemiBin_97 -> 9Sr2::9Sr2_284237`)**:
  - Candidate host (`g__Halobacterium`): $\rho = 0.3348$, within-field $p = 0.0831$, BH $q = 0.4155$ (primary $m=5$).
  - Matched alternative MAG `54PRr2__single_SemiBin_131_sub` (`g__Natronococcus`, `known_different_genus`): $\rho = 0.4831$, within-field $p = 0.0114$, BH $q = 0.0570$.
  - Descriptive correlation difference (`results/stage4_v4_candidate_control_correlation_differences.tsv`): $\Delta\rho = \rho_{\text{cand}} - \rho_{\text{alt}} = -0.1483$ (the matched alternative MAG correlates more strongly descriptively; no p-value or host-preference claim is made because viral-vector permutation tests joint independence from both hosts, not equality of correlations).
- **Remaining Pairs**:
  - `S49Dr1__self_SemiBin_1106 -> 48Sr2::48Sr2_190170` and its control `V48Dr1__self_SemiBin_2681`: 0 host-present sites in `trip5` $\times$ `PlantRoot` (`non_estimable_zero_variance_or_insufficient_sites`).
  - `41PRr2__single_SemiBin_97 -> 55Dr2::55Dr2_75333`, `V48Dr1__self_SemiBin_3215 -> coasm_48_PR::k127_4172870`, and provisional `43PRr3__self_SemiBin_469 -> 44Dr2::44Dr2_4002037`: target absent from historical `viral_tpm.tsv.gz` catalog (`non_estimable_missing_host_or_viral_catalog_row_not_zero`; never zero-filled).

### C. Compact Independent Audit Tables (`independent_audit/`)
Extracted under Slurm job `15894666` directly from the full Shaheen source matrices:
- `independent_audit/stage4_v4_independent_audit_selected47_raw_site_observations.tsv` (SHA256: `d67c2a7ad9769b5dfa099f91e80e3536c51adeac2b8163b221a2a6b2033b5c76`, 282 rows = 6 pairs $\times$ 47 sites): Tidy long-format table with exact host abundance, matched alternative MAG abundance, viral TPM, sample/site/dune_field/campaign/compartment metadata, and source SHA256s.
- `independent_audit/stage4_v4_independent_audit_selected47_wide_matrix.tsv` (SHA256: `be01781e19b701d6a3fde1e65ea9492b13e6b115b95df484cb45f1dd9bd36ed5`, 47 rows $\times$ 27 columns): Wide-format 47-site matrix for direct inspection and independent Spearman/permutation replay.
- `independent_audit/stage4_v4_independent_audit_extraction_receipt.json`: Full Slurm extraction receipt and zero-drift replay verification.

---

## 5. Verification Commands

```bash
cd /home/leechuck/Public/software/empty-quarter-proteome/reproducibility/mechanism_20261109/exp3_viruses/stats_v4
uv run python verify_stats_v4.py
uv run pytest -q tests/test_stage4_statistical_audit_v4.py tests/test_stats_v4_package_and_independent_audit.py
```
