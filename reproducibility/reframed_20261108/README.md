# Reframed Thermoplasmatota Dark-Proteome Package (`reframed_20261108`)

This directory contains the self-contained figure generator, compact source inputs, exported per-panel data tables, SHA-256 source/panel manifest, provenance documentation, and verification tests for the reframed protein manuscript (`../../manuscript/main.tex` and `../../manuscript/supplementary_notes.tex`).

## Entry points

```bash
# Regenerate the four reframed main figures + 14 panel TSVs, run tests, and build manuscript/main.pdf
bash reproducibility/reframed_20261108/run_reframed.sh

# Or run individual steps (using uv or python3):
uv run python reproducibility/reframed_20261108/plot_proteome_reframed.py \
  --out manuscript/figures \
  --panel-dir reproducibility/reframed_20261108/panel_data
uv run pytest -q reproducibility/reframed_20261108/test_proteome_reframed_figures.py
bash reproducibility/build_manuscript.sh
```

Input roots can also be overridden via `--comp-res`, `--repro-dm`, `--arch-res`, or environment variables (`EQ_COMP_RES`, `EQ_REPRO_DM`, `EQ_ARCH_RES`, `EQ_PANEL_DIR`). By default, `plot_proteome_reframed.py` and `test_proteome_reframed_figures.py` resolve all inputs from within this repository (`reproducibility/reframed_20261108/inputs/compute_results`, `reproducibility/darkmatter`, and `reproducibility/integration_20261108/results`).

## Contents

- `plot_proteome_reframed.py`: Generates the four reframed main figures (`fig_archaeal_ecology_reframed.{pdf,png}`, `fig_dark_credibility_reframed.{pdf,png}`, `fig_export_accessory_reframed.{pdf,png}`, `fig_structural_mechanistic_reframed.{pdf,png}`) and exports 14 per-panel TSVs (`panel_fig1a.tsv` through `panel_fig4c.tsv` plus `panel_fig4_c17_struct.tsv`) and `panel_sources_sha256.json`.
- `test_proteome_reframed_figures.py`: Independent verification of all 17 source hashes, 14 panel TSV hashes and row counts, competitive and legacy depth-contrast values, `v2_thermo_metabolism_summary.tsv.gz` module counts, `per_genome_reconciliation.tsv` labelled-call denominators (`dark_fraction == n_dark / n_labelled`, `20,216,302` labelled calls out of `20,451,986` Pyrodigal genes across `5,438` representatives), dN/dS and pN/pS medians, OrthoFinder core/shell/cloud counts (`237,251` dark vs `529,405` non-dark genes), SignalP export rates (`32.5%` vs `11.6%`), canonical Stage 1b adjusted export odds ratios (`panel_fig3c.tsv`), HH-suite/InterProScan all-family summaries, `c17` all-family ESMFold (`19,620` of `105,867`) and ColabFold/AlphaFold2 (`3,453` of `9,644`) structural counts, and the 19-family Stage 3 candidate-audit summary (`panel_fig4c.tsv`).
- `test_three_paper_consistency.py`: Cross-manuscript consistency test suite (10 tests) across the foundational data paper (`empty-quarter-data-current/sn-article.tex`), the ecology companion (`empty-quarter-ecology-current/alt_stats.tex`), and the protein manuscript (`empty-quarter-proteome/main.tex`), including source-to-manuscript verification of completed Stage 1b/2b/3 tables and helper functions, plus a separate submission-readiness gate (`--require-submission-ready`) that fails while `[PENDING:...]` or `DATA-BULK-ARCHIVE-DOI` placeholders remain.
- `rep_ko_v2_provenance.json`: Exact provenance linking Figure 1c (`v2_thermo_metabolism_summary.tsv.gz` via `v2_thermo_metabolism.py`) to `04_functional/repko_clean/results/rep_ko_v2.tsv` and documenting its relationship to `EQ-PROT-BAKTA-v1.0.0` (`19,186,792` proteins) and `EQ-PROT-PYRODIGAL-REP-v1.0.0` (`20,451,986` proteins).
- `inputs/compute_results/`: Checksum-verified compact source tables for Figures 2, 3, and 4, including `science_20261108_completed/` containing the completed Stage 1/1b adjusted export sensitivity tables, Stage 2/2b direct UniRef50 external comparison tables, and Stage 3 candidate-family functional audit tables.
- `panel_data/`: Exported per-panel TSV tables (`panel_fig1a.tsv`, `panel_fig1b.tsv`, `panel_fig1c.tsv`, `panel_fig2a.tsv`, `panel_fig2b.tsv`, `panel_fig2c.tsv`, `panel_fig2d.tsv`, `panel_fig3a.tsv`, `panel_fig3b.tsv`, `panel_fig3c.tsv`, `panel_fig4a.tsv`, `panel_fig4b.tsv`, `panel_fig4c.tsv`, and `panel_fig4_c17_struct.tsv`) and `panel_sources_sha256.json`.
- `science_20261108_provisional/`: Snapshot of the `science-20261108` protein follow-up scripts (Stages 1, 1b, 2, 2b, 3), Slurm batch wrappers, unit/adversarial tests, bundled canonical reconciliation outputs, and compact comparator manifests (`PROVISIONAL_STATUS.md`, `SHA256SUMS`). From the repository root, running `EQ_SCIENCE_OUT=reproducibility/reframed_20261108/science_20261108_provisional EQ_CANONICAL_SCIENCE_OUT=reproducibility/reframed_20261108/science_20261108_provisional uv run pytest -q reproducibility/reframed_20261108/science_20261108_provisional/tests` passes all 40 tests using only bundled package inputs (`39 passed, 1 skipped` without those environment variables). Stage 1/1b adjusted export, Stage 2/2b direct UniRef50 external comparisons, and Stage 3 candidate-family audits are completed and integrated; remaining pending items in `manuscript/main.tex` are `[PENDING:EXT-COMP-HARMON]`, `[PENDING:MECH-FAMILIES-C]`, `[PENDING:TAXONOMY-SNAPSHOT-SENSITIVITY]`, `[PENDING:BULK-MAG-ZENODO-DOI]`, `[PENDING:PROTEIN-STRUCTURE-BULK-RELEASE]`, and `[PENDING:CODE-SNAPSHOT-RELEASE]`.

## Relationship to historical atlas verifiers

In the November 2026 reframing, the general 51.9-million-protein catalogue inventory (`Table S1`, `\label{tab:cohorts}`), external desert read-mapping benchmark, culture/metagenome gene-caller audit, and preliminary 17,221-family whole-chain-versus-crop structural screen (`Figures S1--S5`) were moved from `main.tex` to `supplementary_notes.tex`. Accordingly:
- `python3 reproducibility/darkmatter/check_round7_registered.py` continues to pass all 152 pinned H1--H4 numeric claims.
- All numeric table assertions in `reproducibility/verify.py`, `reproducibility/darkmatter/current_revision_checks.py`, `round4_checks.py`, and `round5_extensions/checks.py` remain intact on the retained tables, while the hardcoded string-presence guards in `verify.py`, `current_revision_checks.py`, and `round4_checks.py` that searched `manuscript/main.tex` for October atlas benchmark sentences now correspond to prose located in `manuscript/supplementary_notes.tex`.
