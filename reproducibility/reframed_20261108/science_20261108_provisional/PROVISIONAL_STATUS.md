# Code and Manifest Snapshot: science-20261108 (Protein Stages 1, 1b, 2, 2b, 3)

This directory contains the snapshot of the protein follow-up analysis code, Slurm batch scripts, unit/adversarial test suites, bundled canonical reconciliation outputs, and compact comparator/genome manifests from `science-20261108` (preserving the frozen directory name `science_20261108_provisional/` for manifest compatibility).

**Status (2026-11-09):**
- **Completed and integrated (`../inputs/compute_results/science_20261108_completed/`):**
  - Stage 1 / Stage 1b within-desert lineage- and length-adjusted SignalP (`sp_any`), TMbed, and combined export sensitivity models under the canonical foundational GTDB release 232 taxonomy and the frozen legacy analytical taxonomy (`panel_fig3c.tsv`, Supplementary Table `tab:export_adj_completed`).
  - Stage 2 / Stage 2b direct search-harmonised UniRef50 sequence-darkness comparisons across 237 Rub' al-Khali and 247 deduplicated external Thermoplasmatota genomes under the canonical foundational GTDB release 232 taxonomy and the frozen legacy analytical taxonomy (Supplementary Table `tab:ext_uniref50_completed`).
  - Stage 3 canonical candidate-family verification and functional audit across 19 evaluated families (`panel_fig4c.tsv`, Supplementary Table `tab:cand_audit_completed`).
- **Remaining pending placeholders in `manuscript/main.tex`:**
  - `[PENDING:EXT-COMP-HARMON]`: SignalP-dependent external Thermoplasmatota export comparisons.
  - `[PENDING:MECH-FAMILIES-C]`: Residue-, architecture-, and pathway-level biochemical corroboration of surviving candidate families.
  - `[PENDING:TAXONOMY-SNAPSHOT-SENSITIVITY]`: Whole-catalogue recurrent dark-family phylum attribution sensitivity under the canonical foundational taxonomy.
  - `[PENDING:BULK-MAG-ZENODO-DOI]` / `DATA-BULK-ARCHIVE-DOI`: Foundational bulk MAG/genome Zenodo archive release.
  - `[PENDING:PROTEIN-STRUCTURE-BULK-RELEASE]`: Assembly-wide protein catalogue, predicted structures, and bulk analysis archive release.
  - `[PENDING:CODE-SNAPSHOT-RELEASE]`: Versioned public archive of the revised reproducibility package.
- **Test status (bundled package inputs):**
  - From repository root: `EQ_SCIENCE_OUT=reproducibility/reframed_20261108/science_20261108_provisional EQ_CANONICAL_SCIENCE_OUT=reproducibility/reframed_20261108/science_20261108_provisional uv run pytest -q reproducibility/reframed_20261108/science_20261108_provisional/tests` -> `40 passed`.
  - Default without those environment variables: `39 passed, 1 skipped`.
- No private correspondence is included.
