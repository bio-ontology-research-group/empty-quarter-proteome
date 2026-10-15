# Are gas-oxidation operons enriched for structurally-novel proteins? — No (clean negative)

Generated 2026-06-16 (`scripts/ibex/82e_operon_darkmatter.sh`, job 47557701). Tests the
"novel-fold-in-gas-operon" hypothesis with the **actual** structurally-novel (reallydark) protein
catalog, replacing the misleading Bakta "hypothetical protein" proxy used in `atmo_operon_context.tsv`.

## Method
Dark-matter proteins are prodigal-called (ids `mag::sample_contig_gene`); gas genes are Bakta locus_tags —
two gene-callings on the same contigs. Bridged by **C-terminal sequence suffix** (prodigal & pyrodigal agree
on stop codons; ~71% of dark proteins recover their Bakta CDS this way). "Dark" = membership in a really-dark
**orthogroup** (reps in `reallydark_ids.txt`, expanded to all members via `orthogroups_cluster.tsv` so we
catch dark proteins in non-representative MAGs, not just the rep's MAG). For each confirmed coxL/hydB gene we
scan ±5 genes and flag dark neighbours; recurrence is tallied by canonical orthogroup.

## Result — the money shot does NOT hold
- **177 / 11,923 gas operons (1.5%)** have ≥1 structurally-novel neighbour (143 MAGs; coxL 89 / hydB 88).
- Structurally-novel proteins exist in only **579 / 6,010 guild MAGs** at all (reallydark proteins span just
  2,213 of 15,469 MAGs) — so for ~90% of guild MAGs the question is moot.
- **No dark orthogroup recurs near gas genes in more than 7 MAGs** (only 2 dark OGs reach 7, both near hydB);
  the rest are essentially MAG-specific. No recurrent novel family travels with the trace-gas machinery.

## Why 82d said 76%
`atmo_operon_context.tsv` counted Bakta **"hypothetical protein"** neighbours (76% of operons). That label
reflects **annotation gaps**, not structural novelty — most such neighbours have UniRef/AFDB homology and are
not in the dark-matter catalog. **Use this file's 1.5%, not 82d's 76%, for any novelty claim.**

## Interpretation
Gas-oxidation operons in EQ MAGs are built predominantly from **characterizable** genes (the enzymes plus
their known maturation/accessory machinery). This is a useful guard-rail: the paper's strength is the
**guild itself** (community-wide atmospheric trace-gas metabolism in taxonomically ordinary organisms), not a
novel-protein-in-operon story. The structural dark matter is a *separate* axis, not specifically coupled to
the gas machinery.

## Files
- `operon_darkmatter_context.tsv` — per gas gene: n_neighbors, n_dark_neighbors, dark-neighbour OG reps.
- `dark_gas_families.tsv` — dark OGs near gas genes ranked by recurrence (n_distinct_mags).
- `dark_gas_neighbors.faa` — the 202 dark neighbour sequences (available for ESMFold/FoldSeek if a fold-level
  recurrence check is ever wanted; note the small absolute count caps the upside).

## Remaining avenue (not pursued)
Recurrence here is by **sequence** orthogroup; a **fold-level** check (map the 202 neighbours onto the tier-2
FoldSeek 3Di clusters) could in principle rescue a structurally-recurrent-but-sequence-divergent family. Given
only 177 operons are involved, the upside is small — flagged, not run.
