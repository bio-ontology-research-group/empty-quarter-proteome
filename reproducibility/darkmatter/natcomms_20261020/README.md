# Genome-resolved atlas revision

This directory supports the resource narrative and visit-matched Thermoplasmatota
analysis in the manuscript. It contains compact evidence and executable source.
Bulk data remain on IBEX under
`/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/natcomms_20261020`.

## Ecological analysis

`scripts/paired_ecology.py` reanalyses retained CoverM abundance estimates from the
complete 551-sample matrix. It excludes the wholly unmapped F5Dr1 profile, averages
replicates within site × visit × compartment, pairs compartments within visits,
and averages paired differences within each site. Sites contribute equally.

The two compartment contrasts are deep minus surface and root-adjacent minus deep.
Both all-read and mapped-read denominators are reported. Two-sided signed-rank
tests use site differences. Holm adjustment covers the four primary
contrast-by-denominator tests. Median-difference intervals use 10,000 site
bootstrap draws, seed 20261020. Each visit is omitted in turn as a sensitivity.

The first completed Slurm run, job 54055936, produced the retained ecological
estimates, preserved on IBEX in `archive_run54055936/`. `ecology_run.json` records the worker and numerical software versions;
`ecology_input_manifest.tsv` pins the three input files. The input abundance
matrix was previously assembled from full CoverM profiles, whose source checksums
are retained in `../heat_stress_20261016/results/input_manifest.tsv`.

The new analysis operates on existing mapping estimates. Raw-read mapping was
not repeated for this revision. The primary depth result is based on 37 sites
and 69 matched visits; 36 sites have a positive subsurface difference.

## Compact verification and figures

From the repository root:

```bash
python3 reproducibility/darkmatter/natcomms_20261020/scripts/verify_compact.py
python3 reproducibility/darkmatter/natcomms_20261020/scripts/plot_atlas.py \
  --out manuscript/figures
```

Verification independently aggregates the retained compartment cells, checks the
cohort and paired outcomes, calculates the four-test Holm block, and writes
`results/ecology_primary_verified.tsv`. The figure script combines the new
ecological result with previously retained phylum-composition, metabolic-marker
and matched-reference results. Each display's denominator is specified in its
caption. `outputs.sha256.tsv` pins the compact analysis results.

## Gene-call inventory

`scripts/catalogue_inventory.py` enumerates the retained proteomes of all 5,438
species representatives, preserving partial-call and length information. It also
extracts complete, internally stop-free proteins of at least 100 residues from
the 126 full-MIMAG high-quality genomes into an IBEX-only FASTA.

These representative gene calls form a separate prediction set from the original
assembly catalogue. Joining darkness by gene identifier alone would mix callers
and prediction coordinates. Sequence-and-locus reconciliation is required before
using this inventory for a length- and completeness-stratified dark-fraction
comparison.

## Resource validation still required

The manuscript's quantified external comparison covers the recurrent dark-family
cohort and its control queries. A whole-catalogue increment at multiple sequence
identity thresholds requires a separately defined, quality-controlled query
universe and external reference panel. The inventory supplies a candidate query
set; it is not itself a reference-search result.

Independent utility should be evaluated using desert datasets whose study
identities are excluded from both the construction data and reference additions.
Compare a fixed baseline with baseline plus Rub’ al-Khali sequences and an
equal-size, length-matched soil-reference addition. Report query/sample-level
recovery, taxonomic resolution, source-study overlap checks, and site- or
study-level uncertainty. Within-catalogue recruitment has a different role: it
describes occurrence in the sampled environment.

Catalogue-wide GUNC screening, sequence/locus reconciliation for the full
representative set, external benchmark execution, public MAG/protein deposition,
and an archival release remain outstanding. The public ENA query on 2026-10-20 is
retained in `results/ena_public_filereport_20261020.tsv`.
