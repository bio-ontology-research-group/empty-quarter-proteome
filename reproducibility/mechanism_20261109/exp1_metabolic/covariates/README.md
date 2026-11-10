# Exp1 484-Genome Cohort: Post-Hoc Covariate Metadata Supplement

This directory provides the companion covariate metadata table, provenance audit, and mutation/fixture test suite for all 484 genomes in the metabolic modeling cohort (`exp1_model_inventory_and_qc_484.tsv`).

> [!IMPORTANT]
> **Data Integrity and Immutability**: The primary model inventory table (`compact_evidence/exp1_model_inventory_and_qc_484.tsv`) remains strictly immutable. This supplement is an independent companion join table designed to support downstream covariate adjustment, lineage stratification, and sensitivity analyses without altering frozen production inputs.

---

## 1. Scientific Background & Methodological Principles

### A. CheckM v1 vs. CheckM2 Separation and Method-Family Compatibility
1. **Tool Discrepancy**:
   - CheckM v1 (Parks et al., 2015) uses lineage-specific marker gene sets identified by HMM search.
   - CheckM2 (Chklovski et al., 2023) uses machine learning models (Gradient Boost and Neural Networks) trained on universal and lineage features.
2. **Empirical Differences on Identical Genomes**:
   - Across the 218 GTDB r232 species representatives, CheckM v1 and CheckM2 yield substantially different estimates on identical assemblies:
     - CheckM v1 completeness: mean 88.10% (std 8.70%, median 90.63%, min 57.07%, max 99.62%).
     - CheckM2 completeness: mean 86.00% (std 9.02%, median 88.05%, min 56.41%, max 100.00%).
     - Individual genome discrepancies reach up to 14.5% in completeness.
     - Based strictly on **numeric completeness/contamination thresholds (`completeness >= 90%` and `contamination < 5%`)**, 108/218 (49.5%) meet the threshold under CheckM v1 vs. 82/218 (37.6%) under CheckM2.
3. **Method-Family Compatibility vs. Full Harmonization**:
   - The 237 Empty Quarter (EQ) desert MAGs were quality-assessed using **CheckM2 v1.1.0** (`evidence/mag-catalogue-20261024/source/checkm2_refined_quality_report.tsv`, models reported as `"Gradient Boost (General Model)"` or `"Neural Network (Specific Model)"`).
   - Official GTDB r232 metadata reports CheckM2 estimates (`checkm2_completeness`, `checkm2_contamination`) and model tier (`"General"` or `"Specific"`), but does not record the exact CheckM2 software version string in `ar53_metadata_r232.tsv.gz` (`checkm2_version = "unknown_gtdb_r232_release"`).
   - Direct CheckM2 provenance across EQ and GTDB r232 therefore makes these **455 genomes** (237 EQ + 218 GTDB) **method-family-compatible** (`checkm2_method_family_compatible = True`), rather than identical-version harmonized. Version and model labels are retained explicitly (`checkm2_version`, `checkm2_model`, `quality_method_sensitivity_flag`).
   - Unknown or unverified quality methods (such as the 29 MGnify entries) are preserved in `reported_quality_*` columns only; their `checkm_*`, `checkm2_*`, and `method_family_checkm2_*` columns are left empty (`""`), never relabeled as CheckM2 even with a `False` compatibility flag.

### B. Formal MIMAG vs. Numeric Quality Threshold Tiers
- Formal MIMAG High Quality (Bowers et al., 2017) requires the presence of 5S, 16S, and 23S rRNA genes plus at least 18 tRNAs in addition to `>= 90%` completeness and `< 5%` contamination.
- Official GTDB r232 release metadata records formal MIMAG flags for the 218 species representatives: **24 High Quality (11.0%)**, **194 Medium Quality (89.0%)**, and **0 Low Quality**.
- Because rRNA and tRNA counts are not present in the metabolic inventory or MGnify comparator inputs, formal MIMAG flags (`mimag_high_quality`, `mimag_medium_quality`, `mimag_low_quality`) are recorded as `"unknown"` for all 237 `EQ_MAG` and 29 `MGnify` entries (never inferred from `>=90%`/`<5%` or set to `False`).
- Numeric completeness/contamination threshold tiers (`comp90_cont5_high`, `comp50_cont10_medium`, `comp_below50_or_cont10_low`, `unknown`) are reported in separate columns (`reported_completeness_contamination_tier`, `checkm_completeness_contamination_tier`, `checkm2_completeness_contamination_tier`).

### C. Lineage Stratification & Confounding Caution
The taxonomic order distributions of the two cohorts differ substantially:
- **EQ Desert Cohort (n=237)**:
  - `o__JACQPN01`: 153 genomes (64.6%)
  - `o__JACPAO01`: 51 genomes (21.5%)
  - `o__SW-10-69-26`: 30 genomes (12.7%)
  - Other/unassigned (`o__RBG-16-68-12`, `o__`): 3 genomes (1.3%)
- **External Comparator Cohort (n=247)**:
  - `o__DHVEG-1`: 131 genomes (53.0%)
  - `o__UBA202`: 46 genomes (18.6%)
  - `o__JACPAO01`: 45 genomes (18.2%) [42 GTDB + 3 MGnify]
  - `o__JACQPN01`: 16 genomes (6.5%) [11 GTDB + 5 MGnify]
  - `o__SW-10-69-26`: 6 genomes (2.4%) [6 GTDB + 0 MGnify]
  - Other (`o__UBA9212`, `o__GCA-002838935`): 3 genomes (1.2%)

Over **71.6%** of the external comparator cohort belongs to orders absent from the desert cohort (`o__DHVEG-1`, `o__UBA202`). Raw, unadjusted whole-cohort comparisons of metabolic gapfill burden or cofactor requirements **could be confounded by** deep phylogenetic differences and cannot by themselves establish desert-specific adaptation.

The supplement flags genomes belonging to the shared order stratum (`shared_order_stratum = True`, covering **301 genomes**: 234 EQ and 67 external; of which **293 genomes** — 234 EQ and 59 GTDB — also have CheckM2 method-family-compatible quality estimates).

---

## 2. Data Provenance, Taxonomy Comparison Types & Usable Counts

| Cohort | Source Label | Total | Taxonomy Comparison Type | Quality Provenance | CheckM v1 Usable | CheckM2 Method-Family Compatible | Formal MIMAG Known |
|---|---|---|---|---|---|---|---|
| EQ Desert | `EQ_MAG` | 237 | `inventory_only_no_independent_source` | `CheckM2_v1.1.0` (`checkm2_refined_quality_report.tsv`) | 0 | 237 | 0 (`unknown`) |
| External | `GTDB_r232_rep` | 218 | `official_gtdb_r232_release_metadata` | Official GTDB r232 (`ar53_metadata_r232.tsv.gz`) | 218 | 218 | 218 (24 HQ, 194 MQ) |
| External | `MGnify` | 29 | `gtdbtk_v2.7.2_r232_classify_wf_rerun` | `unknown_mgnify_catalogue_pipeline` (`reported_quality_*` only) | 0 | 0 | 0 (`unknown`) |
| **Total** | **All** | **484** | — | — | **218** | **455** | **218** |

- **GTDB Release 232 Join Audit (218 genomes)**:
  - Source URL: `https://data.gtdb.ecogenomic.org/releases/release232/232.0/ar53_metadata_r232.tsv.gz`
  - Published/Verified MD5: `0f51fe6f01b4ccbf0d9bab6c07df87af` (SHA-256: `d6ddc0f95d269a6e31af8329e0c3c2af0d90d5dad48f5b3e9d820f340951d8c3`)
  - Accession Join: **218 / 218 exact versioned matches** (`GC[AF]_*.version`) after stripping `GB_`/`RS_` prefix.
  - Taxonomy Verification: **218 / 218 exact full-string and order matches** (0 mismatches) against official `gtdb_taxonomy`.
- **EQ Desert MAG Audit (237 genomes)**:
  - Quality values verified 1:1 against `checkm2_refined_quality_report.tsv` (0 completeness/contamination mismatches; missing or non-finite values fail closed and never default to 0.0).
  - Taxonomy in `exp1_model_inventory_and_qc_484.tsv` derives from the canonical catalogue; because `checkm2_refined_quality_report.tsv` carries no independent taxonomy rerun, `taxonomy_comparison_type` is recorded as `inventory_only_no_independent_source` and `order_match_inventory_vs_independent_r232` is `"unknown"`.
- **MGnify Provenance & NCBI Accession Audit (29 genomes)**:
  - Taxonomy: Independently re-classified with GTDB-Tk v2.7.2 (`classify_wf`, release 232) in `v2_thermo_relatives_gtdbtk_mgnify.tsv`. All **29 / 29 match at Phylum, Class, Order, and Family**; 21/29 match across the full 7-rank string and 8/29 differ at genus or species level.
  - NCBI Accession Candidates: In `v2_thermo_relatives_comparators.tsv`, 9 MGnify rows have no `ncbi_accession` (`no_ncbi_accession`) and 20 carry an **unversioned** assembly stem (e.g. `GCA_013329505`). Because these 20 accessions are unversioned and MGnify FASTAs have `md5_ok = no_published_md5` (no verified assembly identity against NCBI/GTDB), the 7 stems present in `ar53_metadata_r232.tsv.gz` are recorded strictly as `unversioned_candidate_only` (and the other 13 as `unversioned_no_gtdb_r232_match`). None of the 29 MGnify genomes are assigned GTDB `checkm_*` or `checkm2_*` values.

---

## 3. Supplement Schema (`exp1_484_covariate_metadata_supplement.tsv`)

| Column Name | Type | Description |
|---|---|---|
| `genome` | string | Unique genome identifier (primary key, 1:1 join with inventory) |
| `cohort` | string | `EQ_desert` or `external_comparator` |
| `source` | string | `EQ_MAG`, `GTDB_r232_rep`, or `MGnify` |
| `canonical_order` | string | Canonical GTDB order from model inventory |
| `canonical_family` | string | Canonical GTDB family from model inventory |
| `canonical_genus` | string | Canonical GTDB genus from model inventory |
| `canonical_classification` | string | Full canonical 7-rank taxonomy string from model inventory |
| `lineage_stratification_order` | string | Order rank used for lineage-matched sensitivity analysis |
| `shared_order_stratum` | string (`True`/`False`) | `True` if order is shared between EQ and external cohorts (`o__JACQPN01`, `o__JACPAO01`, `o__SW-10-69-26`) |
| `taxonomy_comparison_type` | string | `official_gtdb_r232_release_metadata`, `gtdbtk_v2.7.2_r232_classify_wf_rerun`, or `inventory_only_no_independent_source` |
| `independent_r232_classification` | string / empty | Independent r232 taxonomy string (populated for GTDB and MGnify; empty for EQ_MAG) |
| `order_match_inventory_vs_independent_r232` | string | `True`, `False`, or `unknown` (when no independent source is joined) |
| `full_taxonomy_match_inventory_vs_independent_r232` | string | `True`, `False`, or `unknown` |
| `reported_quality_completeness` | float / empty | Completeness reported in source inventory/comparator table (empty for GTDB_r232_rep) |
| `reported_quality_contamination` | float / empty | Contamination reported in source inventory/comparator table (empty for GTDB_r232_rep) |
| `reported_quality_method` | string | `CheckM2_v1.1.0`, `none_in_frozen_inventory`, or `unknown_mgnify_catalogue_pipeline` |
| `reported_completeness_contamination_tier` | string | Numeric threshold tier (`comp90_cont5_high`, `comp50_cont10_medium`, `comp_below50_or_cont10_low`, `unknown`) |
| `checkm_completeness` | float / empty | Official CheckM v1 completeness (populated for 218 GTDB reps; empty for EQ and MGnify) |
| `checkm_contamination` | float / empty | Official CheckM v1 contamination (populated for 218 GTDB reps; empty for EQ and MGnify) |
| `checkm_marker_count` | int / empty | CheckM v1 marker gene count |
| `checkm_marker_lineage` | string / empty | CheckM v1 marker lineage |
| `checkm_marker_set_count` | int / empty | CheckM v1 marker set count |
| `checkm_strain_heterogeneity` | float / empty | CheckM v1 strain heterogeneity |
| `checkm_completeness_contamination_tier` | string | Numeric threshold tier from CheckM v1 (`unknown` when CheckM v1 absent) |
| `checkm2_completeness` | float / empty | Verified CheckM2 completeness (populated for 237 EQ + 218 GTDB = 455; empty for 29 MGnify) |
| `checkm2_contamination` | float / empty | Verified CheckM2 contamination (populated for 455; empty for 29 MGnify) |
| `checkm2_model` | string / empty | CheckM2 model label (`Gradient Boost (General Model)`, `Neural Network (Specific Model)`, `General`, `Specific`, or empty) |
| `checkm2_version` | string / empty | `1.1.0` (EQ), `unknown_gtdb_r232_release` (GTDB), or empty (MGnify) |
| `checkm2_completeness_contamination_tier` | string | Numeric threshold tier from CheckM2 (`unknown` for MGnify) |
| `method_family_checkm2_completeness` | float / empty | Method-family-compatible CheckM2 completeness (455 genomes; empty for 29 MGnify) |
| `method_family_checkm2_contamination` | float / empty | Method-family-compatible CheckM2 contamination (455 genomes; empty for 29 MGnify) |
| `checkm2_method_family_compatible` | string (`True`/`False`) | `True` for 455 genomes with verified CheckM2 provenance; `False` for 29 MGnify |
| `quality_method_sensitivity_flag` | string | Detailed provenance/sensitivity stratum flag |
| `mimag_high_quality` | string | `True`/`False` from official GTDB r232 metadata (218 genomes); `unknown` for EQ and MGnify |
| `mimag_medium_quality` | string | `True`/`False` from official GTDB r232 metadata (218 genomes); `unknown` for EQ and MGnify |
| `mimag_low_quality` | string | `True`/`False` from official GTDB r232 metadata (218 genomes); `unknown` for EQ and MGnify |
| `mimag_assessment_source` | string | `GTDB_r232_official_release` or `unknown_rrna_trna_not_in_*` |
| `ncbi_assembly_level` | string | Official NCBI assembly level for GTDB reps; `unknown` for EQ and MGnify |
| `mgnify_catalogue` | string / empty | Source MGnify catalogue ID |
| `mgnify_biome` | string / empty | Source MGnify biome path |
| `mgnify_reported_ncbi_accession` | string / empty | Raw `ncbi_accession` reported in `v2_thermo_relatives_comparators.tsv` |
| `mgnify_ncbi_accession_match_status` | string | `not_applicable`, `no_ncbi_accession`, `unversioned_candidate_only`, `unversioned_no_gtdb_r232_match`, etc. |
| `candidate_gtdb_r232_accession` | string / empty | Semicolon-delimited candidate GTDB r232 accession(s) for unversioned/version-mismatched MGnify stems |

---

## 4. How to Reproduce and Test

To regenerate the supplement table and report (including SHA-256 and MD5 verification of all 5 canonical input files):

```bash
uv run python covariates/build_exp1_covariate_metadata_supplement.py
```

To run the unit, acceptance, and negative mutation/fixture test suite:

```bash
uv run pytest covariates/test_covariate_metadata_supplement.py -v
```
