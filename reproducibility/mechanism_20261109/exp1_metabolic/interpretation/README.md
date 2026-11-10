# Experiment 1 Metabolic Modelling — Final Descriptive Interpretation (`484` Genomes)

## 1. Scope and Design Authority

This module implements the final descriptive interpretation of the completed 484-genome metabolic gap-filling ensemble (`237 EQ_desert` MAGs and `247 external_comparator` genomes: `218 GTDB_r232_rep` + `29 MGnify`) and the 7-genome same-medium GPR masking controls (`Slurm job 15894236`), following:

- `mechanism-20261109/EXP1-DESCRIPTIVE-INTERPRETATION-IMPLEMENTATION.txt`
- `mechanism-20261109/EXP1-FULL-COHORT-INTERPRETATION-DESIGN.txt`
- `mechanism-20261109/EXP1-INTERPRETATION-FIGURE-REVIEW.txt`

This subdirectory (`exp1_metabolic/interpretation/`) and its review mirror (`mechanism-20261109/review-20261109/interpretation/`) are self-contained and read-only with respect to upstream model execution, collector, and Stage 4 deficiency search code.

---

## 2. Inputs and Portability

All default paths in `build_exp1_descriptive_interpretation.py` are resolved relative to the package tree (`SCRIPT_DIR = Path(__file__).resolve().parent` and `PACKAGE_DIR = SCRIPT_DIR.parent`), so a reviewer does not need the author's home directory layout:

1. **Full-cohort case summaries (`12,580` rows across `484` genomes)**:
   - Portable combined summary included in package: `exp1_full484_combined_summary_12580.tsv`
   - Or upstream 8 shard summaries via `--shards-dir`: `exp1_shard_00_summary.tsv` .. `exp1_shard_07_summary.tsv`
     - `exp1_shard_00_summary.tsv`: `5182219833aac046f33beccfc7483e355e6073c8b6cd4c815cfc0cd5e3ac1b3e`
     - `exp1_shard_01_summary.tsv`: `b8bc464b2c9c72d83c8a5ecaa7b735614a7e15cada379c738aa845d442000740`
     - `exp1_shard_02_summary.tsv`: `15886a3eed4ddd7790626a54b0b08b626c9c2ed59bc5e36f03ee678faacb00c3`
     - `exp1_shard_03_summary.tsv`: `3b28adffb53c6ca35380178f86e873425f420a3ef9f7747b4bfe18e33032904e`
     - `exp1_shard_04_summary.tsv`: `94372f5d9a8899557cfe06af5ee251e4a1f1337c6beba08fa19a624342c9d1bd`
     - `exp1_shard_05_summary.tsv`: `820a9572bea37e748aa837407199ee302c755891a73d25fc528916823aa21c8e`
     - `exp1_shard_06_summary.tsv`: `3ee7456574b6fe5dc9f0d413f29f175e1607aef7d32542d92408c987d768cc21`
     - `exp1_shard_07_summary.tsv`: `0e30d8d6b2503b4c944913c2ff17cad6bf202929c3e0cd16a5e5be8e727e438b`
2. **Corrected 484-genome covariate supplement**:
   - `../covariates/exp1_484_covariate_metadata_supplement.tsv` (`sha256: c3492b38be9978d3b81fd269ae81ab7a0baa41e7d8554a2772d312e5083a016c`)
3. **Immutable 484-genome model inventory**:
   - `../compact_evidence/exp1_model_inventory_and_qc_484.tsv` (`sha256: 044014df6aa1f7bf9b5d4fcaea3a2daea0ef62f3dfcb3c79ef08e54b0347eb2b`)
4. **Same-medium 7-genome masking controls (`job 15894236`) and paired audit reference**:
   - `inputs/exp1_masking_controls_7hq_v4.tsv` (`sha256: 8268f1f5372c0f60190e5b21b6cb25b8b72e0f342fe65bf52f70b8ce5013a53f`, 63 rows)
   - `inputs/exp1_masking_controls_7hq_v4_baselines.tsv` (`sha256: 80fbbe6761b2c1c4bd0d5a2093cb70c0770c55b84f803c9fca2f5b8e132ab528`, 7 rows)
   - `inputs/exp1_full484_paired_biomass_sensitivity.json` (`sha256: abb636e05de0dd1986371f2c31903c319d58a64b76c7a63da090142fde0db4e8`)

Any explicitly supplied CLI path (`--shards-dir`, `--combined-summary`, `--covariate-supplement`, `--inventory`, `--masking-cases`, `--masking-baselines`, `--root-paired-audit-json`) is validated fail-closed and never silently falls back to a default path if missing.

---

## 3. Reproduction Commands

```bash
# 1. Regenerate all tables, figures (PDF + PNG), JSON report, and SHA-256 manifest
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  uv run python build_exp1_descriptive_interpretation.py

# Or explicitly from the portable combined 12,580-row summary table:
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  uv run python build_exp1_descriptive_interpretation.py \
  --combined-summary ./exp1_full484_combined_summary_12580.tsv

# 2. Run unit, mutation, and acceptance test suite (without using /tmp)
mkdir -p ./.pytest_tmp
TMPDIR=./.pytest_tmp OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  uv run pytest -v test_descriptive_interpretation.py
rm -rf ./.pytest_tmp
```

---

## 4. Output Artifacts and Summary of Findings

### A. Case Universe & Solver Status (`exp1_case_universe_and_solver_status_summary.tsv`)
- **Total genomes**: `484` (`237 EQ_desert`, `247 external_comparator`: `218 GTDB_r232_rep` + `29 MGnify`).
- **Core cases**: `24` per genome (`11,616` core cases total, zero missing or duplicate).
- **Optional `S0_gapseq_native_predicted` cases**: `964` evaluated (`482 × 2`); absent (`4` cases) only for `GCA_016180755.1` and `S22PRr1__self_SemiBin_29`, both independently verified on Slurm to have empty native/desert media.
- **Total evaluated cases (`12,580`)**:
  - `10,578` (`84.09%`) valid feasible ensembles (`heuristic_incumbent_bounded`, `n_distinct_valid_solutions > 0`, max residual `4.51e-10 <= 1e-6`, max bound violation `2.84e-12 <= 1e-6`, min verified growth `0.0100000028 >= 0.01`).
  - `1,936` (`15.39%`) `infeasible_under_fixed_medium` (all `484` genomes across the `4` minimal `H2+CO2+NH4` conditions `S1/S2 × B0/B2`).
  - `66` (`0.52%`) numerically unknown cases: `56` `fba_verification_failed` and `10` `timeout_no_incumbent`.

### B. Within-Genome Paired Biomass & Redox Sensitivity (`exp1_paired_biomass_and_redox_sensitivity.tsv`, `fig1_paired_b0_vs_b2_distributions.pdf`)
- **Paired `B0` minus `B2` (`omit 7 non-F420 cofactors, keep F420`)**:
  - `S0_gapseq_native_predicted`: `EQ_desert` median reduction `73.0` reactions (IQR `70.0–76.0`, usable `235/236` present, `1` optional-medium unavailable, `1` numerical unknown; `211` reductions, `0` ties, `24` increases); `external_comparator` median `73.0` (IQR `69.0–77.0`, usable `242/246` present, `1` optional-medium unavailable, `4` numerical unknown; `216` reductions, `0` ties, `26` increases).
  - `S3_anoxic_desert_oa_aa`: `EQ_desert` median `73.0` (IQR `70.0–76.0`, usable `236/237`, `1` numerical unknown); `external_comparator` median `73.0` (IQR `70.0–76.0`, usable `246/247`, `1` numerical unknown).
  - `S4_oxic_desert_oa_aa`: `EQ_desert` median `73.0` (IQR `69.5–75.0`, usable `231/237`, `6` numerical unknown); `external_comparator` median `73.0` (IQR `70.0–75.75`, usable `246/247`, `1` numerical unknown).
  - `S5_anoxic_defined_aa_precursors`: `EQ_desert` median `80.0` (IQR `78.0–82.0`, usable `237/237`); `external_comparator` median `79.0` (IQR `78.0–81.0`, usable `247/247`).
  - `S6_oxic_defined_aa_precursors`: `EQ_desert` median `79.0` (IQR `78.0–82.0`, usable `237/237`); `external_comparator` median `79.0` (IQR `78.0–81.0`, usable `247/247`).
  - `S1/S2` minimal `H2+CO2+NH4`: `0` usable pairs in both cohorts (`237/237` EQ and `247/247` external are `fixedmedium_infeasible_pairs`, `0` `numerical_unknown_pairs`).
- **Paired `B1` (`omit all 8 cofactors incl. F420`) minus `B2` (`keep F420`)**:
  - Median difference `-11.0` to `-10.0` reactions across `S3..S6` in both cohorts (`B2` requires ~10–11 more added reactions than `B1` when F420 is retained in the biomass objective).
- **Broader recipe sensitivity (`B0` vs `B1`, `B3`, `B4`)**:
  - Median `B0 - B1` is `84.0–90.0` reactions, `B0 - B3` is `103.0–109.0` reactions, and `B0 - B4` is `121.0–128.0` reactions across `S3..S6`, with similar median sensitivity in both cohorts.

### C. Descriptive Cohort Comparisons & Quality/Lineage Stratification (`exp1_order_and_quality_composition_overlap.tsv`, `exp1_checkm2_stratified_cohort_comparisons.tsv`, `exp1_mgnify29_reported_quality_sensitivity.tsv`, `exp1_unadjusted_all484_descriptive_summary.tsv`, `fig2_quality_order_composition_and_overlap.pdf`)
- **Lineage & Quality Method Denominators**:
  - `455/484` genomes (`237 EQ_MAG` + `218 GTDB_r232_rep`) have CheckM2 method-family quality estimates (`method_family_checkm2_*`).
  - `293/455` CheckM2-compatible genomes (`234 EQ_MAG` + `59 GTDB_r232_rep`) belong to the 3 shared orders (`o__JACQPN01`, `o__JACPAO01`, `o__SW-10-69-26`).
  - `29/484` `MGnify` genomes have `reported_quality_method = "unknown_mgnify_catalogue_pipeline"` (`8` in shared orders, `21` in non-shared orders) and are kept strictly separate in `exp1_mgnify29_reported_quality_sensitivity.tsv`.
- **Within-Stratum (`canonical_order × CheckM2 completeness bin × CheckM2 contamination bin`) Contrasts**:
  - Across `35` occupied order × quality cells, only `9` CheckM2 cells have valid models in both cohorts, and `5` are non-sparse (`n >= 5` in both cohorts):
    1. `o__JACPAO01 | [70,90) | [0,5)` (`EQ n=19`, `GTDB n=25`): in `S3`, `B0` median is `182.0` vs `186.0` (`EQ - GTDB = -4.0`); `B2` median is `106.0` vs `109.0` (`-3.0`).
    2. `o__JACPAO01 | [90,101] | [0,5)` (`EQ n=15`, `GTDB n=10`): in `S3`, `B0` median is `176.0` vs `177.5` (`-1.5`); `B2` median is `106.0` vs `101.5` (`+4.5`).
    3. `o__JACQPN01 | [70,90) | [0,5)` (`EQ n=25`, `GTDB n=9`): in `S3`, `B0` median is `173.0` vs `164.0` (`+9.0`); `B2` median is `100.0` vs `97.0` (`+3.0`).
    4. `o__JACQPN01 | [90,101] | [0,5)` (`EQ n=35`, `GTDB n=5`): in `S3`, `B0` median is `158.0` vs `161.0` (`-3.0`); `B2` median is `88.5` (`n=34`) vs `88.0` (`+0.5`).
    5. `o__SW-10-69-26 | [70,90) | [0,5)` (`EQ n=14`, `GTDB n=5`): in `S3`, `B0` median is `169.5` vs `172.0` (`-2.5`); `B2` median is `97.5` vs `100.0` (`-2.5`).

### D. Same-Medium Masking Controls Across 7 Reference Genomes (`exp1_masking_7genomes_descriptive_summary.tsv`, `fig3_masking_recovery_by_7genomes.pdf`)
- **Primary unit of analysis**: `7` high-completeness, low-contamination MGnify reference genomes (`reported_quality_completeness >= 90%`, `reported_quality_contamination <= 5%`; formal MIMAG status unknown because rRNA/tRNA criteria are not recorded in the metabolic inventory), each evaluated across `9` nested masking conditions (`3` mask fractions `10%, 20%, 30%` × `3` deterministic seeds = `63` evaluations under matched `S4_oxic_desert_oa_aa` + `B2`).
- **Per-genome mean recovery rate** (`mean` of 7 genome means = `0.168216`, `median` of 7 genome means = `0.151613`, range `0.099285` in `MGYG000519748` to `0.248099` in `MGYG000520046`).
- **Per-genome mean alternative candidate fraction** (`mean` of 7 genome means = `0.767300`, `median` = `0.791368`, range `0.634500` to `0.815900`).

---

## 5. Scientific Limitations

1. **Template sensitivity vs physiological trajectory**: Similar cohort medians under cofactor removal indicate similar sensitivity to the declared biomass template in both cohorts, not identical reaction-count distributions or a measured physiological trajectory.
2. **Heuristic bounds vs global minima**: Reported `min_added_rxns_count` values are minima among found valid time-limited MILP solutions (`heuristic_incumbent_bounded`), not proven global optima. Weighted reaction objective costs mean that relaxing a biomass recipe can occasionally increase unweighted reaction counts in individual models.
3. **Infeasibility vs numerical failure**: Fixed-medium infeasibility (`1,936` cases in `S1/S2`) is a property of the declared candidate network and minimal `H2+CO2+NH4` medium, distinct from numerical unknown outcomes (`66` cases) and from biological non-growth across all possible pathways.
4. **Confounding and quality provenance**: Sharing the CheckM2 method family does not guarantee identical software versions or training models (`CheckM2 v1.1.0` for EQ MAGs vs GTDB R232 metadata). Order-level stratification does not eliminate within-order phylogenetic structure, and descriptive differences are not interpreted as environmental adaptation.
5. **Masking controls and dark-protein function**: In the 7 tested reference genomes, alternative reactions account for `76.7%` of newly induced gap-fill additions while only `16.8%` of masked GPR reactions are recovered. Gap-filling alone does not establish specific enzymatic functions for dark proteins, and these observations apply to the 7 selected reference genomes rather than all archaea.
