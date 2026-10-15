# EQ-MAGs: proposed structure, narrative and display set

Status: PROPOSAL (2026-09-16, revised 2026-09-17). Replaces `PAPER_STRUCTURE.md` only when Alejandra
confirms. This file holds the STRUCTURE only. Open items live in `TODO_RESTRUCTURE.md`, current
values in `CATALOGUE_NUMBERS.md`, what changed and why in `CORRECTIONS.md`, method detail in each
workbench `_clean/README.md`, and how each analysis works in `dossiers/`.

**Writing method.** The manuscript is rewritten from scratch: `paper/main.tex` becomes `main_v1.tex`,
a new `main.tex` keeps the preamble, macros and bibliography, and no prose is reused. Drafting reads
only the dossiers, the registry and `READING_NOTES.md`. Nothing enters from an agent's memory or
from an exploratory folder without a registered rerun.

**Lead decision, changed openly.** The 2026-09-09 decision made structure the headline. This
proposal moves structure to one controlled layer and leads with the controlled census, because no
armour item has run, the alignment-normalised result is now a null, and both Nature-tier dark-matter
precedents carried wet-lab validation. If the armour lands with the candidate rate clearly above the
floor, the structure section moves into the title and nothing else changes.

## 1. The story in one paragraph

Earth's largest contiguous sand sea was sampled five times over 31 months at 60 sites and three soil
positions, with geochemistry. The organisms are ordinary at every rank above order, and the companion
16S paper already describes their ecology. What 16S cannot see is the gene content, and a tenth of the
genome-resolved proteome has no homolog in any reference. This paper measures that unknown against a
control at every layer and asks three questions the amplicon paper could not: is it real, who carries
it, and why is it unknown. It is real: the same families recur across all five expeditions, they
evolve under the purifying selection of annotated genes, and long-read genomes of species cultured
from the same soils carry dark protein at the same rate within measurement error. It has owners:
lineage, not assembly quality, sets how dark a genome is, and the largest single owner is an unnamed
archaeal radiation. Two hundred and thirty-seven Thermoplasmatota species with no named genus, in a
GTDB class holding twenty-one genomes worldwide, are a third unknown against four to seven per cent
for archaea and bacteria matched on abundance, prevalence and genome quality. They occur in nine
samples in ten, dominate the archaea, increase with depth, and are aerobic heterotrophs rather than
methanogens or rhodopsin halophiles. It is unknown because no reference has sampled here: what does
match is architecturally rearranged known parts, and biome-matched sequence buys more coverage per
residue than generic database depth. The annotated fraction encodes the expected desert toolkit, and
salt, not aridity, filters both the known and the unknown gene content. Known phyla do not imply
known biology; the gap is in the reference set.

## 2. Title and abstract

Title candidates:
- Known phyla, unknown proteome: a controlled census of Rub' al-Khali soil (73 chars)
- An unnamed archaeal radiation carries the unknown proteome of a sand sea (72)
- Lineage sets the unknown fraction of the soil proteome (54; safest if the archaea are cut)

Abstract (Nature Microbiology, 150 words, unreferenced):

> Most microbial proteins match no characterised family, and whether that unknown fraction is biology
> or artefact cannot be settled in pooled catalogues. We shotgun-sequenced 551 soil samples from 60
> sites across the Rub' al-Khali on five expeditions over 31 months, recovering 14,426 genomes in
> 5,438 species, all within named phyla. Of the proteins encoded on binned contigs, 9.89 % have no
> reference homolog. They form 106,042 families; 79.4 % recur in all five expeditions, they evolve
> under the purifying selection of annotated genes, and they reappear at the same rate in long-read
> genomes of species cultured from the same soils. Lineage, not assembly, sets how dark a genome is,
> and one owner dominates: 237 Thermoplasmatota species with no named genus are a third unknown
> against matched controls, and supply a fifth of all dark families from 4 % of the species.
> Biome-matched references recover more of this proteome than databases five times deeper. Salinity,
> not aridity, filters known and unknown gene content alike.

Basis rule: 9.89 % is the 37,949-bin union, not the 5,438 representatives. Name the basis wherever a
dark rate appears; four different rates circulate.

## 3. Results sections

Heading style is per section, not uniform: a section earns a declarative heading only where a control
in it licenses the claim. R1 has no control, so it takes a topic label, not a claim. No section opens
with "We asked whether".

⚠ The abstract drafted above quotes archaeal figures (237 species, 0.310, a fifth of dark families)
that are NOT in the registry. They cannot go into prose until registered, whatever the tree shows.

| # | section | carrying evidence | control | figure |
|---|---|---|---|---|
| R1 | A five-expedition genome catalogue of the Rub' al-Khali is taxonomically ordinary | 551 samples, 60 sites, 5 expeditions over 31 months, 3 positions, 21.6 Gbp median depth; 14,426 bins, 5,438 species, 51 phyla; first empty rank at order (11), family (142), genus (1,343), species (3,755); 10 order gaps adjudicated, 0 proposals | none, and say so. GUNC passed on the 10 candidates | Fig 1 |
| R2 | A tenth of the encoded proteome has no reference homolog, and it is not an assembly artefact | 9.89 % genome-resolved; dark tracks length not completeness; 106,042 families, median prevalence 98/551, 79.4 % in all five expeditions; dN/dS 0.425 vs 0.388; cultured same-species genomes at the same rate within error, residual confined to single-strain pairs | length/size-matched annotated families; cultured long-read genomes | Fig 2 |
| R3 | Lineage, not genome quality, sets how much of a genome is unknown | quality R² 0.115, +phylum 0.768, +genus 0.835; variance unique to lineage 0.720 vs 0.020; novelty gradient 0.217 to 0.017; darkness falls with genome size | quality covariates, novelty-matched stratum, size normalisation | Fig 3 |
| R4 | An unnamed archaeal radiation carries a fifth of the unknown proteome | 237 species, 0 named genera, class holds 21 genomes worldwide; 0.310 vs 0.048 against matched archaea (6.5x) and 0.069 against matched bacteria (4.5x); 19.8 % of dark families from 4.4 % of species; 481/549 samples, deeper than surface at 40/45 sites; aerobic, no mcrABG, no amoABC, no rhodopsin | matched null on completeness, contamination, size, N50, abundance, prevalence; answers Maeke 2025 | Fig 3 |
| R5 | Structural novelty is architectural, not local, against a blinded floor | 2,838 confident models, 1,031 with no fold match; 6.6x the control floor standardised on confidence, 8.1x on confidence and length; 96.7 % locally familiar, 63.7 % covered end to end | 3,013 blinded known-fold controls co-searched | Fig 4 |
| R6 | The gap is in the reference set: biome-matched content beats depth | 28.74 % absent from 1.26 G sequences vs 0.41 % control; desert catalogues leave 42.6 % unexplained at a fifth of MGnify90's size; depth resolves 34 % of control absences but 6 % of dark; 55.6 % of genes carry no KO; unpredicted-but-present genes are zero in 420 of 448 samples | matched annotated control in the same pass | Fig 5 |
| R7 | The annotated fraction encodes the expected desert toolkit | trace gas, repair, osmolytes, BGC, resistome, metals tracking geology rather than human input | clade placement, biosafety census | Extended Data |
| R8 | Salt, not aridity, filters known and unknown gene content alike | geochemistry R² 0.132, 0.089 beyond space; partition on the XRF cohort gives geochemistry twice the unique variance of space; S 54,227 families, 91.7 % negative; light sorts species with depth, desiccation uniform, salt replaces the community | vanadium null, within-phylum and within-species controls, space partialling | Fig 6 |

Order: organisms, then how much is unknown, who owns it, the archaea, what it looks like, why it is
unknown, what the known part does, how the environment sorts both.

Word budget (3,500 for Introduction, Results and Discussion): Intro 400 · R1 250 · R2 450 · R3 300 ·
R4 450 · R5 350 · R6 400 · R7 200 · R8 300 · Discussion 450. Trim R7 first, then R8; if still long,
R3 folds into R4 as its opening paragraph. Methods are uncounted and re-cut; the frozen M.1–M.14
numbering does not carry over.

Two hazards while writing: 55.6 % is both the no-KO ceiling and the NHEJ prevalence, so never put
them unqualified in one paragraph; and four dark rates circulate on four bases, so each carries its
denominator in its own sentence.

## 4. Introduction, four paragraphs

1. Drylands are 45 % of land and expanding; hot hyperarid soil is the low-biomass case the 2025
   contamination consensus names; this interior had no genome-resolved survey. The companion 16S
   paper mapped the ecology and listed what 16S cannot do: genes without primer bias, function to
   genome, strains, archaea.
2. The unknown proteome: the fraction depends on reference and rule; fold space is near-saturated and
   residual novelty is architectural; unknown families are reported as reference absence, with
   lineage but without biome persistence, environment, genomes behind every protein, cultured
   complements or a measured error floor.
3. The two design limits of pooled catalogues: context stripped, error rate unmeasured.
4. Here we: one biome, five expeditions, genomes, cultured strains, geochemistry, a control at every
   layer. Close on the gap sentence.

## 5. Discussion

Comparisons live here only, never in Results. Source: `DISCUSSION_NOTES.md`.

- What the intersection adds over the taxon-axis, fold-space and global-catalogue papers: persistence,
  a cultivation control, and the biome-matched reference gain.
- Content beats depth is an actionable recommendation for soil reference building.
- Trace gas here is poorer than in cold deserts and CO-weighted; not a chemosynthesis paper.
- Aridity is not null for composition; the nulls are trait-level. Salt agrees with the amplicon mosaic.
- Limits: capacity not activity; half the representatives below 90 % complete; XRF is total elements
  with no pH; the cultivation control reaches only cultivable bacteria; contamination decision pending.
- The companion annotation study, stated as ceilings rather than a promise: the enzyme route reaches
  8.4 % of the dark set and recovery falls from 0.45 to 0.05 with target identity.
- Outlook: species already in hand as strains; Thermoplasmatota as the culture and fold target; strain
  persistence across expeditions as the named next study.

## 6. Display set

Six main figures, the Nature Microbiology allowance.

| item | content | source | state |
|---|---|---|---|
| Fig 1 | map, timeline, positions, funnel, first-empty-rank ladder | `fig_study_design.png` | redraw in vector, add the ladder, fix the caption basis |
| Fig 2 (R2) | dark ladder; length vs completeness; expedition persistence; dN/dS vs control; cultured vs MAG | `fig_darkfam_reality` (missing) + PoC 05 | rebuild with new panels |
| Fig 3 (R3+R4) | variance decomposition; residual dark by phylum; novelty gradient; archaeal dark distribution, depth and physiology | PoC 01 | new; two sections share one figure |
| Fig 4 (R5) | structural calibration, four panels | `fig_structural_calibration.pdf` | exists; update to the standardised floor |
| Fig 5 (R6) | absence vs control; content vs depth; two ceilings | `fig_specificity_controls` (rebuilt 2026-09-17) + `fig_refbias_ceiling` | compose |
| Fig 6 (R8) | PCoA and PERMANOVA blocks; per-axis families; site-paired stress modes | `fig_darkmatter_pcoa/permanova` + `fig_stress_modes` | compose |
| ED 1–5 | trace-gas guild; BGC and resistome; unbinned reservoir; genome size vs KO richness; contamination diagnostics | existing + PoC | as is |

Cut: the endemism set, `fig_picrust2_bias`, `fig_lifeonair_by_phylum` (dead bases). Drop Box 1 and
Box 2; negatives move into the sections that need them, each with its detection limit.

## 7. What enters from the dark_function project

| result (722 cultured genomes) | placement |
|---|---|
| definition overlap of the darkness rule | supplementary table |
| coding-evidence gate: real in aggregate, unattributable per protein | one Discussion sentence and a supplementary table |
| pN/pS 0.355 vs 0.112 | supplementary, same sign as the catalogue result |
| culture-to-MAG synteny, 96.5 % of 737 pairs | supplementary; 20 representatives is too thin for main text |
| mobile-element and defence enrichment, and all the negatives | companion chapter, not this paper |
| annotator reach and identity collapse | one Discussion sentence, as a ceiling |
| binning loss | companion chapter; its own verdict is an upper bound |

## 8. Target journal

**Nature Microbiology, Article**: 3,500 words, 6 display items, 150-word abstract, headed Results
required. Its Resource type (4,000 words, 6 items) is the re-entry if the Article is desk-rejected.
**Fallback: Nature Communications**, 5,000 words and 10 items, at effectively the same time to
acceptance (242 days against 241), so the fallback costs no schedule.
Nature and Science are not realistic without experimental validation of a dark function or a
cross-biome comparator; both Nature dark-matter precedents carried wet-lab work. ISME J is
deprioritised despite topical fit: its scope page names MAG-only descriptive work as out of scope,
Methods precede Results and count toward the limit, and it has no Resource vehicle.
The October deadline buys a submission, not a decision.
