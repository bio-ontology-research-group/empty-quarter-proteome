# RQ-C2 — Why does taxonomy-based gene prediction fail? Are these pangenomes open?

**Status:** historical exploratory output; not promoted to the active manuscript

## Answer
An archived pangenome workflow (IBEX job 47843640, completed 2026-06-28) produced a strong accessory-enrichment
pattern, but the product has not been closed against the current 14,426-bin/5,438-representative catalogue or the
active claim/provenance ledger. It is therefore retained as a historical exploratory result, not as a mechanism
for the PICRUSt2 comparison in RQ-C1.

Across **345 non-redundant Empty Quarter species** (dRep ≥95% ANI clusters with ≥10 genomes; **8,635 genomes**, one
per sample per species to avoid binning-duplicate inflation), MMseqs2 gene-family pangenomes show:

| Metric | Value |
|--------|------:|
| Species with **open** pangenome (Heaps' α < 1) | **337 / 345 = 97.7%** |
| Median Heaps' α | **0.555** (≪1 = strongly open) |
| Median core / pan ratio | **0.065** (only ~6.5% of a species' genes are core) |
| Dark genes: accessory fraction | 8.24% |
| Dark genes: core fraction | 4.79% |
| **Dark enrichment in accessory vs core** | **OR = 1.79, Fisher exact p < 1e-300** |

Within that historical workflow, UniRef50-no-hit genes were approximately 1.8-fold over-represented in its
accessory class relative to its core class. The association does not by itself show that accessory status causes
reference-prediction error, nor that the result transfers to the current catalogue.

## Evidence / reproduce
- IBEX pipeline: `workbench .../scripts/ibex/140_pangenome_openness.{py,sbatch}` (copies in this folder); job 47843640
  (`--account=pi-hohndor`, 32 cpu/256 G). Raw outputs in `$M/06_pangenome/out/` (`family_table.tsv` 221 MB stays on IBEX).
- Pulled here: `enrichment_summary.json`, `species_openness.tsv`. Regenerate the figure + stats:
  `python3 RQ/C2/finalize.py` (reads the pulled summaries; Fisher via scipy).
- Figure: `fig_pangenome_openness.{pdf,png}` (α distribution + dark-in-accessory enrichment).

## Caveats
- The exact input-membership and intermediate-family artifacts are not part of the active checksum-closed
  manuscript chain. Until they are promoted and verified, these statistics must not be presented as current
  manuscript evidence.
- Openness is computed per dRep species cluster with ≥10 genomes; smaller species are excluded (underpowered for
  saturation). De-duplication to one genome per (species, sample) prevents the two binning passes from inflating the core.
- "dark" = UniRef50 no-hit (sequence-dark); the accessory/core split uses ≥95% genome occupancy for core.

**References:** tettelin2005pangenome; mcinerney2017openpangenome (keys in `manuscript/references.bib`).
