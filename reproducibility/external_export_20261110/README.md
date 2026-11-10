# Harmonised export comparison across Thermoplasmatota genomes (10 November 2026)

This package holds the analysis behind the external export comparison in the
manuscript (Results, Methods and the Supplementary Table "Harmonised export
comparison across Thermoplasmatota genomes").

## Question

Do sequence-dark genes carry predicted signal peptides more often than other
genes in Thermoplasmatota from other habitats, when every genome is processed
with the same gene caller, reference search and signal-peptide predictor?

## Inputs

The inputs are the verified tables of the external-validation run on Shaheen
(`/scratch/project/k10018/eq-protein-external-20261110`, handoff status
`COMPLETE_VERIFIED`):

- gene calls: Pyrodigal v3.7.1, metagenomic mode, 1,957 genomes and 3,617,011 genes
- UniRef50 2026_01: MMseqs2 `search -s 7.5 -e 1e-3 --max-accept 10`, one row per gene
- SignalP 6.0h, fast mode, organism group `other`, one row per gene
- genome manifest with cohort, GTDB release 232 taxonomy and quality
- 95% ANI groups across all resources (skani)

Cohorts: 237 Rub' al-Khali representatives, 247 GTDB/MGnify references,
35 Ca. Penumbrarchaeia genomes (Maeke et al. 2025), 354 SMAG and 1,084 SPIRE
Thermoplasmatota genomes. Cross-resource duplicates and SPIRE entries without
genome assets are excluded.

## Method

`code/export_comparison.py` joins the three per-gene tables with exact
identifier closure and computes, per cohort, the Mantel-Haenszel odds ratio of
a predicted signal peptide in genes without versus with a UniRef50 match,
stratified by genome and seven protein-length categories. Intervals are
block-bootstrap percentiles (2,000 draws, seed 20261110); a block is a named
genus, or the 95% ANI group of a genus-unassigned genome. Subsets: all genomes,
genomes with at least 70% completeness and at most 5% contamination, and the
four named orders that contain the desert genomes. The analysis is repeated
for complete interior ORFs.

Run on Shaheen as Slurm job 15940343 with `code/sbatch_export_comparison.sh`.

## Results

- `results/export_by_cohort.tsv`: estimates per gene set, subset and cohort
- `results/export_or_ratio.tsv`: ratio of the desert odds ratio to each external odds ratio
- `results/per_genome_all_genes.tsv`, `results/per_genome_complete_interior.tsv`: per-genome counts
- `results/receipt.json`: seed, draws and totals
