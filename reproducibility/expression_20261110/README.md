# External Metaproteomic Expression Validation (`reproducibility/expression_20261110`)

This module packages the compact reproducible code, frozen Stage 1 v2 and Stage 2 audit tables, and unit/regression test suite for the external soil metaproteome reanalysis of **Karaoz et al. (2026)** (*Nature Communications* 17:10740; DOI: `10.1038/s41467-026-78137-0`; Blodgett Forest whole-soil warming experiment, 3 paired plots x 2 treatments x 5 soil depths = 30 samples collected June 2018).

## Directory Layout

- `code/`
  - `expression_validation_core.py`: Core library implementing translation-table-consistent ORF validation, MS-GF+ `.mzTab` metadata-verified Q-value parsing (`MS-GF:QValue <= 0.01`, `MS-GF:PepQValue <= 0.01`), streaming I/L-normalized Aho-Corasick peptide specificity, unique-spectrum set-union quantification (`aggregate_unique_spectra_for_entity`), zero-pseudocount-free paired differences, and exact whole-plot sign-flip permutation (`n=3` plots, 8 configurations, minimum achievable two-sided `p = 0.25`).
  - `stage1_v2_karaoz_homology_eligibility.py` & `.sbatch`: Stage 1 v2 homology and search-database eligibility audit across 2,321 dereplicated MAGs (113 Archaea, 38 Thermoplasmatota) and 122,370 detected EMSL50964 proteins.
  - `stage2_karaoz_peptide_fdr_and_warming_expression.py` & `.sbatch`: Stage 2 `.mzTab` peptide-FDR, whole-background I/L specificity, and paired warming sensitivity pipeline.
  - `download_massive_mztab.py`: Checksum-verified downloader for the 30 `.mzTab` files in MassIVE `MSV000097826`.
- `results/`
  - `mztab_download_manifest.tsv`: Per-file byte size and SHA-256 manifest for all 30 `.mzTab` files.
  - `stage1_v2_homology_eligibility/`: Frozen Stage 1 v2 summary JSON and eligibility/audit ledgers.
  - `stage2_expression_validation/`: Frozen Stage 2 summary JSON, 30-sample `.mzTab` audit, 1,514-family expression & non-estimability ledger, 333-protein archaeal PSM confidence audit, 2,305-peptide specificity ledger, and 30-sample paired proteome + June 2018 DNA table.
- `manifests/`
  - `bulk_input_and_missing_search_db_manifest.tsv`: External bulk FASTA/mzTab inventory and documentation of the missing 15-split composite search FASTAs.
  - `sha256_manifest.tsv`: Byte sizes and SHA-256 digests of all files packaged in this module.
- `tests/test_expression_validation.py`: Self-contained unit and frozen-output regression tests.
- `verify_expression_20261110.py`: Standalone verifier script for this module.

## Key Scientific & Methodological Boundaries

1. **Non-Estimable Primary Lead 1, Lead 2, and Family/Clade Warming Claims**:
   - **Lead 1 (`22PRr1::22PRr1_2531955_3`, PilJ-profile-neighbour context)**: Zero strict full-length homologs in the 113 Karaoz archaeal MAGs (`1` moderate bidirectional hit in `metabat.3300060323.474` whose protein was not included in `EMSL50964_Proteins.faa`) -> `nonestimable_search_database_inclusion_unverified`.
   - **Lead 2 (`V60PRr2::V60PRr2_1314979_4`, cell-envelope-profile-neighbour context)**: Zero strict or moderate homologs and zero HMMER domain hits (`E <= 1e-5`) in the 113 Karaoz archaeal MAGs or `EMSL50964_Proteins.faa` -> `nonestimable_no_homolog_in_study_genomes`.
   - **Missing 15-Split Search Database FASTAs**: Exhaustive audit of all 604 deposited files in MassIVE `MSV000097826` and ESS-DIVE `2588706`/`3365729` confirmed that the 15 split composite search FASTAs used for two-stage 5%-FDR database filtering were not deposited (`full_background_verified = False`). Consequently, primary family-level and clade-level warming responses are strictly `nonestimable_full_search_background_unavailable`, and all paired differences are reported solely as `observed_subset_provisional_sensitivity` diagnostics (`n = 3` paired plots, minimum two-sided `p = 0.25`).
2. **Conservative Observed-Subset Thermoplasmatota Detection**:
   - Across 78 `EMSL50964` proteins linked to the 38 Thermoplasmatota MAGs: `2` have `>=2` unambiguous peptides at `Q <= 0.01` and `PepQ <= 0.01`, `15` have `1` unambiguous peptide, `4` have `1` clade/family-shared peptide, `34` have positive spectral counts in the published CSV but zero qualifying `Q <= 0.01` unambiguous peptides (`rejected_spectral_count_without_qualifying_q01_unambiguous_peptide`), and `23` are `undetected_in_soluble_metaproteome_not_proof_of_absence`.

## Verification Commands

```bash
uv run pytest reproducibility/expression_20261110/tests/test_expression_validation.py -v
uv run python reproducibility/expression_20261110/verify_expression_20261110.py
```
