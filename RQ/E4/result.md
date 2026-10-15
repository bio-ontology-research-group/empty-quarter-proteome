# RQ-E4 — Result

## ANSWER (short)
Dark protein families are **taxonomically narrow per family but carried by mainstream,
abundant phyla.** At the family level the dark fraction is overwhelmingly **taxon-restricted**:
**97.5%** of families are confined to a single phylum (median `n_phyla` = 1) and **54.4%**
to a single genus (median `n_genera` = 1); only **2.5%** span ≥2 phyla. So in the strict
sense, dark families are *not* broadly distributed across the tree — each is a lineage-local
innovation. **But** this is fully consistent with the paper headline that dark proteins are
*not genomically/taxonomically isolated*: the families are not hiding in exotic deep-branching
lineages — their carriers are the **dominant arid-soil phyla** (Actinomycetota 35.5%,
Pseudomonadota 10.4%, plus Thermoplasmatota, Myxococcota, Chloroflexota). The dark matter
lives in ordinary, well-classified, abundant organisms; it is narrow *per family* simply
because protein families are fine-grained, not because it sits on isolated branches.

## Key numbers (n = 23,382 families)
| metric | median | mean | max |
|---|---|---|---|
| n_phyla / family | 1 | 1.03 | 7 |
| n_genera / family | 1 | 1.04 | 161 |
| n_species / family | 0 | 0.01 | 51 |
| n_mags / family | 4 | 5.6 | 663 |

- **Single-phylum (n_phyla ≤ 1): 22,787 (97.5%)**; **broad (≥2 phyla): 595 (2.5%)**.
- **Single-genus (n_genera == 1): 12,711 (54.4%)**; ≤1 genus (incl. unclassified):
  19,120 (81.8%); **multi-genus (≥2): 4,262 (18.2%)**.
- n_phyla histogram: {1: 22,787, 2: 546, 3: 41, 4: 4, 5: 1, 6: 2, 7: 1}.

### Carrier phyla (dom_phylum, top 8)
| phylum | families | % |
|---|---|---|
| Actinomycetota | 8,293 | 35.5 |
| Thermoplasmatota | 3,699 | 15.8 |
| Pseudomonadota | 2,436 | 10.4 |
| Myxococcota | 1,635 | 7.0 |
| Chloroflexota | 1,491 | 6.4 |
| Acidobacteriota | 1,248 | 5.3 |
| Gemmatimonadota | 1,182 | 5.1 |
| Planctomycetota | 803 | 3.4 |

These are exactly the mainstream, abundant phyla of the EQ MAG catalogue (Actinomycetota
and Pseudomonadota are the two largest), confirming the dark fraction is broadly *spread
across* the community's major lineages even though each individual family is narrow.

## Interpretation vs the headline
- "**Not isolated**" (commit 82g) refers to genomic/lineage context: dark proteins sit in
  normal genomes of dominant taxa, not on exotic isolated branches. **E4 supports this** —
  carriers are the top arid-soil phyla.
- "**Taxon-restricted per family**" (97.5% single-phylum) is a *different* axis and is also
  true. The two are not in conflict: a fine-grained protein family is naturally
  lineage-local, while the *aggregate* dark reservoir is distributed across all major phyla.

## Caveats
- Breadth is measured **within the EQ MAG catalogue**, not the global tree; a single-genus
  family here may occur elsewhere worldwide.
- `n_species` is ≈0 for almost all families (GTDB species-level sparsity for uncultivated
  MAGs); phylum and genus are the reliable breadth axes — species breadth is uninformative.
- Single-genus dominance is partly an artefact of protein-family granularity; the
  carrier-phylum view is the more biologically meaningful distribution.
- A handful of very broad families exist (max n_genera = 161, n_phyla = 7) — these are the
  rare cross-phylum dark families and are the better candidates for genuinely widespread,
  possibly horizontally-shared, novel proteins.
