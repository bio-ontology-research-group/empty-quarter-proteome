# CLAUDE.md — Empty Quarter MAGs paper

Guidance for coding agents and human collaborators.

## Scope and repositories

This is the independent proteome reproducibility package. The manuscript editing checkout is
`~/Documents/papers/empty-quarter-proteome/`, Overleaf project `6a40af16a879c0f96d4c1761`.
The package is public at `bio-ontology-research-group/empty-quarter-proteome`. The separate Life on Air manuscript is edited at
`~/Documents/papers/empty-quarter-lifeonair/`, Overleaf project `6ad0d1b4338ca7cd7d291188`;
Its separate public analysis package is `bio-ontology-research-group/empty-quarter-lifeonair`,
checked out at `~/Public/software/empty-quarter-lifeonair/`.

The paper answers biological questions. Put scientific findings in Results and headings,
procedural controls and qualifications in Methods, and limitations in Discussion. Avoid
census framing and contrastive headings such as "lineage, not genome quality".

The amplicon paper is in `../empty-quarter-amplicon/`; the
defensome, mobilome, phage, plasmid, CRISPR, and retron work is in
`../empty-quarter-phages/`. Exchange data between repositories only through versioned,
checksum-pinned interface tables.

Heavy computation and raw data live in the shared workbench at
`/home/leechuck/Public/software/empty-quarter-workbench` and on IBEX/ORIX. Keep manuscript source,
compact evidence tables, figure code, manifests, and reproducibility checks here.

## Current scientific boundary

- Current taxonomy: 14,426 refined bins; 5,438 ANI95 HQ/MQ representatives; GTDB-Tk 2.7.2/R232;
  51 named phyla; no first gap at domain/phylum/class; 11 first gaps at order, 142 at family,
  1,343 at genus, and 3,755 at species.
- Ten order-gap representatives were tree-adjudicated and none supports an automatic novel-rank
  proposal. One low-MSA order-gap candidate remains unadjudicated.
- The 15,469-MAG/51.9-million-protein functional snapshot is historical and is not interchangeable
  with the current 14,426-bin catalogue.
- Structural no-hits are reference- and metric-relative candidates. Report query-normalised TM,
  retrieval decomposition, confidence/length gates, and known-fold controls. Foldseek 10
  alignment-normalised TM has an audited normalisation defect; retain its raw values as software
  diagnostics and exclude its threshold subdivisions from biological interpretation. The active
  H4 endpoint pools 1,031 initial dual-reference qTM-low models; see the round7 disposition.
- Historical dark-orthogroup dN/dS is descriptive and conditional. Do not infer a catalogue-wide
  selection regime; the 148-family greater-than-one tail is hypothesis generation only.
- Catalogue specificity means no detected homolog in explicitly named comparisons; never convert
  it into universal absence or geographic endemism.
- Capacity is not activity or flux; analogue/ligand predictions are not validated functions.

## Reproducibility

The public manuscript snapshot is `manuscript/main.tex`; the dedicated editing checkout is named above.
A clean clone must reproduce without the submodule:

```bash
bash reproducibility/fetch_large_files.sh
bash reproducibility/run_all.sh
```

Every active quantitative claim belongs in `reproducibility/verify.py` and
`reproducibility/claims.tsv`. Compact frozen inputs go in
`reproducibility/current_manuscript_artifacts.tsv`; heavy historical values require an exact entry
in `reproducibility/historical_provenance_ledger.tsv`. Never label a TABLE or PROVENANCE claim as
locally recomputed. Inspect the final Git diff and rerun the full harness before committing.

## Working records and precedence

- Open items: `notes/OPEN_ITEMS.md`. Results structure: `notes/PAPER_STRUCTURE.md`. Per-analysis
  records: `notes/dossiers/`. Superseded planning files: `notes/archive/`.
- Recomputed numbers are registered in `notes/CORRECTIONS.md` (old → new, basis, date, source).
- Recomputations run in the workbench, one self-contained folder per analysis:
  `/ibex/scratch/projects/c2014/rob/metagenomics/<area>/<theme>_clean/` — each carries a `README.md`
  with verified numbers, input paths and job IDs, plus its own `scripts/` and `results/`. Analysis
  outputs stay there; this repo holds the manuscript and final cited tables and figures.
- Precedence: a value registered in `notes/CORRECTIONS.md` supersedes older values elsewhere in the
  repo. Superseded result directories are banner-marked rather than rewritten, and entries in
  `reproducibility/claims.tsv` are changed only by coordination, so that registry may lag a correction.

## Compute policy

- Never execute analysis, database extraction, compilation, bulk parsing, or bulk checksumming
  on IBEX login nodes, including small helper steps delegated to agy. A `.sbatch` filename
  does not create an allocation: submit it with `sbatch`. Before executing compute, use
  `reproducibility/darkmatter/round7_analyses/require_slurm_worker.sh COMMAND ...`
  (mirrored at the remote analysis root). It rejects login hosts and requires the current
  host to belong to a running Slurm allocation. Keep agents themselves on the workstation.
- Submit CPU/high-memory workflows to IBEX and actual GPU workloads to ORIX through Slurm.
- Do not reserve a GPU for download, setup, polling, or CPU-only work.
- Do not run heavy jobs on the local workstation.
- Record database versions, parameters, input/output checksums, job IDs, and closure counts in the
  workbench changelog at promotion time.
- Treat `plans/` as historical planning records, not proof that an analysis was delivered.

## Current blockers

The current-catalogue PICRUSt2 replacement is **no longer gated** (resolved 2026-08-16): CoverM profiles,
the corrected MAG-KO matrix (`rep_ko_v2`) and the checksum-pinned PICRUSt2 predictions are all in place, and
§2.1 is computed and written on the 548/448-sample cohorts. Only the PI's own fail-closed 573/573 validation
chain remains stalled on four failed CoverM array tasks; it is not a dependency for this paper. MQ/LQ recovery and targeted
reassembly/coassembly remain prospective routes to higher-order novelty; do not state a new phylum
or class until read support, quality control, dereplication, R232 classification, and phylogenetic
adjudication all pass.
