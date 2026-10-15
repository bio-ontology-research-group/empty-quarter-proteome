# RQ-B5 — Result

**Question.** How is the dark proteome partitioned by depth/compartment
(Surface / Deep / PlantRoot-rhizosphere)?

## ANSWER

The dark proteome is **weakly but significantly compartment-structured**:
roughly **one-third of dark-protein families (7,156 / 23,382 = 30.6%) show
significant compartment specificity** (Kruskal-Wallis across Surface/Deep/
PlantRoot, BH-FDR q < 0.05), while the **majority (16,226 = 69.4%) are broadly
distributed** across all three compartments. The strongest partition is along
**depth and the rhizosphere axis**: the rhizosphere (PlantRoot) and Deep
compartments each carry ~2.9k enriched families, more than twice the
surface count.

This is fully consistent with the community-level signal: compartment is a
*significant but minor* axis of variation (PERMANOVA **R² = 0.020, p = 0.001**;
`phase_b_summary.txt`), below space (R² = 0.038) and climate (R² = 0.047). So
the dark proteome carries a real depth/rhizosphere imprint, but most dark
families are cosmopolitan across compartments rather than compartment-private —
the "known taxa, unknown biology" reservoir is largely **shared** across the
desert profile.

## Numbers

| Compartment | # samples | # dark families enriched | % of all families |
|---|---|---|---|
| Surface | 116 | 1,287 | 5.5% |
| Deep | 165 | 2,919 | 12.5% |
| PlantRoot (rhizosphere) | 266 | 2,950 | 12.6% |
| **Any (q<0.05, compartment-specific)** | 547 | **7,156** | **30.6%** |
| **Broadly distributed (not significant)** | — | **16,226** | **69.4%** |

- Families tested (finite p): 22,430 of 23,382.
- PERMANOVA compartment R² = 0.020 (p = 0.001); unique R² | space = 0.020.
- Per-compartment survival-marker context (`uv_survival/compartment.tsv`):
  Surface MAGs are the most photolyase-rich (22.4% vs 13.1% Deep) and
  carotenoid-rich (67.6% vs 55.1% Deep), i.e. the compartment partition tracks
  expected UV/desiccation gradients even where dark-family enrichment is modest.

## Interpretation / ecological context

The rhizosphere and deep compartments host the largest sets of compartment-
enriched dark families. This matches the dryland **"fertile islands" /
rhizosphere-effect** paradigm — plant roots concentrate carbon, water, and
microbial activity in otherwise oligotrophic desert soil, building distinctive
sub-communities (Garcia-Pichel et al. and the dryland-microbiome / desertification
synthesis, `resources/literature/dryland-ecology/2024_dryland_microbiomes_desertification_ISME_wrae056.pdf`).
The pre-computed rhizosphere indicators
(`indicators/plantroot_vs_other_top.tsv`) provide candidate decodable members of
this rhizosphere-enriched dark set. Surface yields the fewest enriched families,
in part real (thin, UV/desiccation-stressed habitat) and in part reduced
statistical power (fewest samples).

## Caveats

- Unequal sampling depth (Surface 116 < Deep 165 < PlantRoot 266) inflates power
  for the larger compartments; the smaller-compartment counts are conservative.
- CoverM abundance = community membership/coverage, not expression — this is
  *where dark families reside*, not in-situ activity.
- "enriched_in" is assigned by maximum mean abundance and can be marginal for
  families with near-equal compartment means.

## Files
- `run.py` — reproducible pipeline (streams the 92 MB matrix).
- `compartment_partition.tsv` — per-family means, prevalences, KW H, p, q, enriched_in.
- `compartment_summary.tsv` — headline counts.
- `fig_compartment.{pdf,png}` — enriched-family counts per compartment + example splits.
