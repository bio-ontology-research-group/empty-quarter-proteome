# Rub’ al-Khali proteome: reproducibility package

This paper asks which cellular roles could explain the unknown protein repertoire of widespread desert archaea. The reframed manuscript (`manuscript/main.tex` and `manuscript/supplementary_notes.tex`, working title *Widespread desert archaea harbour an export-enriched unknown protein repertoire*) centres on the disproportionate sequence-dark repertoire of widespread shallow-subsurface Thermoplasmatota, its coding and evolutionary support, its association with accessory orthogroups and predicted export, and structural and genomic-context hypotheses. Full catalogue inventory tables, technical benchmarks and historical cohort screens are retained in the Supplementary Information (`manuscript/supplementary_notes.tex`).

## Reproduce the reframed manuscript and figures

Requirements: Python 3.12 (tested) with `pip install -r requirements.txt` (or `uv`), `latexmk`, and a TeX distribution including `natbib` and `authblk`.

```bash
# Regenerate the four reframed main figures + 14 exported panel TSVs, run figure/source tests, and build manuscript/main.pdf
bash reproducibility/reframed_20261108/run_reframed.sh

# Or run the steps individually:
uv run python reproducibility/reframed_20261108/plot_proteome_reframed.py \
  --out manuscript/figures \
  --panel-dir reproducibility/reframed_20261108/panel_data
uv run pytest -q reproducibility/reframed_20261108/test_proteome_reframed_figures.py
bash reproducibility/build_manuscript.sh
```

See [`reproducibility/reframed_20261108/README.md`](reproducibility/reframed_20261108/README.md) for the full manifest of compact inputs, panel data tables (`panel_fig1a.tsv`--`panel_fig4c.tsv`, `panel_fig4_c17_struct.tsv`), `rep_ko_v2.tsv` provenance (`rep_ko_v2_provenance.json`), three-paper consistency and submission-readiness tests (`test_three_paper_consistency.py`), and the snapshot of `science-20261108` follow-up scripts (`science_20261108_provisional/`).

## Retain and run historical catalogue verifiers

```bash
IBEX_USER=your_kaust_id bash reproducibility/fetch_large_files.sh
bash reproducibility/run_all.sh
```

Large derived inputs remain on IBEX. The fetch command requires SSH/rsync and authorised project c2014 access, and verifies SHA-256 checksums. The historical harness regenerates historical figures, verifies frozen numerical evidence and builds `manuscript/main.pdf`. For checks requiring only committed compact inputs without historical `main.tex` prose guards, run `python3 reproducibility/darkmatter/check_round7_registered.py` (which checks all 152 H1–H4 claims). Note that in the November 2026 reframing, the general 51.9-million-protein catalogue inventory and benchmark prose were moved from `manuscript/main.tex` to `manuscript/supplementary_notes.tex`; consequently, hardcoded string guards in the historical October scripts (`reproducibility/verify.py`, `reproducibility/darkmatter/current_revision_checks.py`, and `reproducibility/darkmatter/round4_checks.py`) that search `manuscript/main.tex` for October atlas benchmark sentences now correspond to prose in `manuscript/supplementary_notes.tex`.

The visit-matched Thermoplasmatota analysis and historical atlas figures are in
[`reproducibility/darkmatter/natcomms_20261020`](reproducibility/darkmatter/natcomms_20261020), and the November 2026 competitive-mapping integration scripts and tables are in [`reproducibility/integration_20261108`](reproducibility/integration_20261108).

## Contents and scope

- `manuscript/`: Active reframed LaTeX (`main.tex`, `supplementary_notes.tex`), full author list (`authors.tex`), bibliography (`references.bib`), rendered `main.pdf`, and figures (`figures/`).
- `reproducibility/reframed_20261108/`: Self-contained reframed figure generator, compact inputs, exported per-panel TSVs, SHA-256 manifest, `rep_ko_v2.tsv` provenance, provisional `science-20261108` code snapshot, and active pytest suites.
- `reproducibility/`: Entry points, checksum manifests (`manuscript_snapshot.sha256`), historical claim registry, retained analysis scripts and compact data.
- `code/`, `results/`, `RQ/`: Analysis sources, derived evidence and historical supporting analyses.
- `reviews/`: Review reports, responses and verification evidence.

## Data and repository status

- ENA umbrella BioProject: `PRJEB104209`
- ENA 16S rRNA V3--V4 amplicon BioProject: `PRJEB106069` (`1,186` deduplicated analytical profiles)
- ENA shotgun metagenome and MAG BioProject: `PRJEB125153` (`ERP204894`; `554` runs across `551` samples; `5,437` MAGs `ERZ29258846`--`ERZ29264282` plus boundary representative `48Dr2__self_SemiBin_25` at `10.0%` contamination under `ERZ29258845`; public release scheduled for `25 November 2026`)
- Released Knowledge Graph v3.0.1 Zenodo archive: `10.5281/zenodo.23134168`
- Bulk MAG/genome Zenodo archive: pending (`[PENDING:BULK-MAG-ZENODO-DOI]`)
- Bulk protein/structure archive (`51,929,855` assembly proteins, `28,159,496` clusters, `106,042` dark-family sequences/alignments/models): pending (`[PENDING:PROTEIN-STRUCTURE-BULK-RELEASE]`)

## Integrated evidence, 10 November 2026

The manuscript incorporates the full microTrait profiles, metabolic-model sensitivity and masking controls, archaeal virus–host candidates, and conserved surface-associated gene pairs. The strongest family-level evidence identifies conserved modules and fold-level hypotheses; no new biochemical function is established. Independent soil proteogenomic data lack adequate coverage to validate the two focal leads. Harmonised external export comparisons remain pending.

Reproduce the new evidence with `reproducibility/mechanism_20261109/run_tests.sh`, `reproducibility/expression_20261110/verify_expression_20261110.py`, and the commands in `reproducibility/integration_20261110/surface_figure/README.md`. Full microTrait matrices, provenance and the ecology sensitivity are maintained in the companion data foundation package. Bulk inputs remain on institutional storage.
