# Historical dark-family dN/dS screen

This directory preserves an exploratory historical analysis. It is not a catalogue-wide formal
test of positive selection.

`positive_families.tsv` contains 148 dark families selected because their mean pairwise NG86
dN/dS exceeded one; 123 also pass the recorded `mean_dS > 0.05` and `n_pairs >= 3` heuristic.
Across the 148, the median dN/dS is about 1.34, the maximum is 3.65, the median family contains
three members, and 14 phyla are represented. Eleven records carry an environmental-gradient flag.
These quantities describe an ascertained screen tail and do not establish selection or its
prevalence.

The best-EC fields are low-confidence decoding hints rather than annotations. The per-family
values are suitable for prioritising an independently confirmed codon-model analysis, with
explicit treatment of recombination, saturation, phylogenetic dependence, and small family size.
No claim of adaptation is made from NG86 alone.

The exact historical scripts and result-table hashes are enforced by
`reproducibility/dnds_legacy_provenance.tsv`. The current bounded interpretation and reproducible
summary are in `RQ/E2/`.
