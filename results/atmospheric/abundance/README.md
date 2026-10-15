# 16S abundance weighting of the atmospheric trace-gas guild

Generated 2026-06-16 by `code/make_16s_weighting.py`. Translates the genome-resolved guild (MAGs with
high-affinity hydB 1h/1l or form-I coxL) into a **community-abundance** estimate using the companion 16S
amplicon dataset (1,227 samples; `empty-quarter-amplicon/repository/cache/`). Answers: *what fraction of the
desert community can live on air?*

## Method
Amplicon taxonomy is GTDB-style, so MAG and ASV taxonomy link directly at **genus** (precise) and **phylum**
(full-coverage fallback). Per taxon we compute guild **prevalence** = guild MAGs / total MAGs (denominator =
all 15,469 GTDB-classified MAGs). Per sample, community guild fraction = Σ_taxon relative_abundance(taxon) ×
prevalence(taxon). Reported for the **full** community and the **alive** (non-relic-DNA) subset, by compartment.
- `prev_any/H2/CO` = prevalence-weighted (primary, defensible estimate).
- `any_presence` = presence-based upper bound (any taxon with ≥1 guild MAG counts fully).
- `covered_fraction` = share of community abundance in taxa that have ≥1 MAG (the rest is not assessable → the
  prevalence-weighted number treats it as 0, i.e. it is a **lower bound** at genus level).

## Files
| File | What |
|---|---|
| `guild_abundance_summary.tsv` | mean guild fraction by table × rank × metric × scope (overall + compartment) |
| `guild_abundance_by_sample.tsv` | per-sample guild fraction + covered fraction |
| `per_genus_prevalence.tsv`, `per_phylum_prevalence.tsv` | guild prevalence per taxon (n_total, n_H2, n_CO, n_any) |

## Headline (full community, prevalence-weighted)
| Rank | any guild | CO (form I) | H₂ (1h/1l) | presence upper-bound | covered |
|---|---|---|---|---|---|
| **genus** (conservative) | **14.3%** | 13.2% | 2.4% | 29.9% | 64% |
| **phylum** (full coverage) | **28.7%** | 19.5% | 15.7% | 73.9% | 98% |

**So ~14–29% of the EQ soil community (by 16S) belongs to taxa with the genomic capacity for atmospheric
trace-gas oxidation — up to ~74% as an upper bound.** CO-oxidation capacity is more widespread than H₂ across
the community. The signal is **uniform across compartments** (Surface 14.3% / Deep 14.8% / PlantRoot 14.0% at
genus) — a pervasive baseline survival strategy, not a niche specialisation.

## Secondary observation — the 'alive' (non-relic) community
Guild fraction is markedly **lower** in the alive subset (genus 0.2%, phylum 6.7%) than in the full community.
This is intriguing — consistent with trace-gas scavenging being a **persistence/maintenance strategy of dormant
or low-activity cells** rather than the most active live populations, or with guild taxa being enriched in the
relic-DNA pool. **Interpret cautiously**: the alive ASV set is small (863 ASVs) and depends on the amplicon
paper's relic-correction. Flagged for follow-up, not a headline claim.

## Caveats
- Capacity (gene presence), not in-situ flux (see `../CLIMATE_IMPLICATIONS.md`; resolve with GC incubations + metatranscriptomics).
- Prevalence-weighting assumes within-genus guild prevalence estimated from MAGs applies to amplicon members of that genus.
- Genus-level underestimates (36% of community in genera lacking any MAG → counted as 0); phylum-level overestimates within-phylum uniformity. The truth lies between; we report both.
