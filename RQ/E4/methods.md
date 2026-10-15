# RQ-E4 — Methods

**Question.** Are dark protein families taxon-restricted or broadly distributed across
the tree?

## Input
`results/darkmatter_ecology/fam_taxa_summary.tsv` — 23,382 dark families. Columns:
`family_id, n_members, n_mags, n_classified_members, n_phyla, dom_phylum, n_genera,
dom_genus, n_species, top_lineage, top_lineage_n`. Taxonomy is GTDB-Tk r226 assigned to
the MAGs carrying each family.

## Procedure (`run.py`, stdlib + matplotlib)
1. **Breadth distributions.** Median / mean / max of `n_phyla`, `n_genera`, `n_species`,
   `n_mags` per family.
2. **Restricted vs broad classification.**
   - single-phylum = `n_phyla ≤ 1`; broad = `n_phyla ≥ 2`.
   - single-genus = `n_genera == 1`; multi-genus = `n_genera ≥ 2`. (`n_genera == 0` means
     no member was classified to genus level, kept as a separate "≤1 genus" bucket.)
3. **Carrier phyla.** Counts of `dom_phylum` (the dominant phylum of each family).
4. **Figure.** `fig_taxbreadth.pdf/png`: log-scaled histograms of `n_phyla` and `n_genera`
   per family, plus a bar of the top carrier phyla.

## Caveats baked into the method
- GTDB classification is sparse at the **species** level for these uncultivated MAGs, so
  `n_species` is ~0 for almost all families and is **not** a reliable breadth axis; phylum
  and genus are the informative axes.
- "Taxon-restricted" here means *restricted within this MAG set*; a family confined to one
  genus in the EQ catalogue could still occur elsewhere globally. Breadth is measured
  relative to the 15k-MAG catalogue, not the whole tree of life.
- A family being single-genus is partly a property of protein-family granularity (fine
  clustering → narrow families); the biologically meaningful read is *which* taxa carry the
  dark fraction (carrier phyla), not just the per-family span.

## Outputs
`breadth_summary.tsv`, `carrier_phyla.tsv`, `nphyla_hist.tsv`,
`fig_taxbreadth.pdf/png`, `_summary.txt`.
