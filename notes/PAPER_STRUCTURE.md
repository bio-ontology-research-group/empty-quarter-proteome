# EQ-MAGs: structure, narrative and display set

Status: structure in force for the v2 draft (2026-09-20), written from the 2026-09-17 proposal after an
independent evaluation; Alejandra's sign-off pending on the two changes marked ★. This file holds the
STRUCTURE only. Open items: `OPEN_ITEMS.md`. Current values: `CATALOGUE_NUMBERS.md`. What changed and
why: `CORRECTIONS.md`. Method detail: each workbench `_clean/README.md`. How each analysis works and
what it licenses: `dossiers/`. Archived predecessors: `archive/`.

**Writing method.** Rewritten from scratch: `paper/main_v1.tex` is the archive, `paper/main.tex` the
new manuscript, no prose reused. Drafting reads only the dossiers, the registry and
`resources/literature/READING_NOTES.md`. Nothing enters from an exploratory folder without a
registered rerun in a workbench `_clean/` folder; until then it is marked `\unreg{}` in the source.

**Two decisions already made (Alejandra):** the archaea get their own Results section but not the
title; the journal is not fixed, so the draft is written to the tightest standard (3,500 words for
Introduction, Results and Discussion; six display items; 150-word abstract) so the top targets stay open.

## 1. The story in one paragraph

Earth's largest contiguous sand sea was sampled five times over 31 months at 60 sites and three soil
positions, with geochemistry. The organisms are ordinary at every rank above order, and the companion
16S paper already describes their ecology. What 16S cannot see is the gene content, and a tenth of the
genome-resolved protein clusters have no homolog in any reference. This paper measures that unknown
against a control at every layer and asks three questions a pooled catalogue cannot: is it real, who
carries it, and why is it unknown. It is real: the same families recur across all five expeditions,
they carry the purifying signature of coding genes, and long-read genomes of species cultured from the
same soils carry dark protein at the same rate once the species is sampled on both sides. It has
owners: lineage, not assembly quality, sets how dark a genome is, and the largest single owner is two
unnamed archaeal classes. It is unknown because no reference has sampled here: what does match at the
structural level is known parts in unmatched arrangements, and biome-matched sequence buys more
coverage than generic database depth. Geochemistry sorts the known and the unknown gene content alike,
beyond what distance explains. Known phyla do not imply known biology; the gap is in the reference set.

## 2. Title and abstract

Working title: **Known phyla, unknown proteome: a controlled census of Rub' al-Khali soil**.
Alternatives if the archaea are cut: "Lineage sets the unknown fraction of the soil proteome".
Archaea stay out of the title (decision).

Abstract: written last, 150 words, unreferenced. Sentence 1 states a concrete fact about the object;
gap within two sentences; "Here we" carries the design and the headline count with its denominator;
one sentence per Results section; close on the reference-set sentence. Every archaeal figure it quotes
must be registered first (R3/R4 rows in the registry).

Basis rule: 9.89 % is a rate over 90/80 cluster representatives with a member on a contig of the
37,949-bin union; it is not a rate over proteins or over the 5,438 representatives. Four dark rates
circulate on four bases; each carries its denominator in its own sentence.

## 3. Results sections (seven) ★

Heading style is per section: a declarative heading only where a control in the section licenses the
magnitude it asserts; a label for construction and for uncontrolled censuses. No section opens with a
question or "We asked whether".

| # | heading | carrying evidence (registered unless marked) | control | figure |
|---|---|---|---|---|
| R1 | A five-expedition genome catalogue of the Rub' al-Khali (label) | 551 samples, 60 sites, 1,043 km, five expeditions over 31 months, three positions, median 21.58 Gbp; 14,426 bins, 5,438 species, 51 phyla; first empty rank order 11 / family 142 / genus 1,343 / species 3,755; 10 order gaps adjudicated, 0 proposals | none, and say so; GUNC clean on the 10 candidates | Fig 1 |
| R2 | A tenth of the encoded protein clusters have no reference homolog, and they behave as coding sequence | 9.89 % (4,292,446 of 43,423,322); 106,042 families; median prevalence 98/551; 79.4 % in all five expeditions (28.9 % re-assembled in all five); dN/dS 0.425 vs 0.388, direction stated; cultured same-species genomes, paired n = 23 (unreg.) | length- and size-matched UniRef50-matched families; 23 cultured long-read bacterial genomes (reach stated: bacteria only) | Fig 2 |
| R3 | Lineage, not genome quality, sets how much of a genome is unknown | quality R² 0.115, +phylum 0.768, +genus 0.835; unique lineage 0.720 vs quality 0.020; novelty gradient 0.217 → 0.017; phylum term keeps 52 % after novelty; darkness falls with genome size (all unreg.) | quality covariates; novelty rank; size | Fig 3 |
| R4 | Two poorly sampled archaeal classes carry a fifth of the unknown proteome | 237 species in c__SW-10-69-26 (185) and c__E2 (51), 0 named genera; GTDB r232 holds 21 and 197 species representatives in those classes; median dark 0.310 vs 0.048 matched archaea (6.5x, 218 pairs) and 0.069 matched bacteria (4.5x, 237 pairs); 19.8 % of single-phylum dark families from 4.4 % of species; occurrence and depth pattern (all unreg.) | 1:1 matched null on completeness, contamination, size, N50, abundance, prevalence, within domain and across | Fig 3 |
| R5 | Structural novelty is architectural, not local, against a blinded floor | 2,838 confident models; 1,031 (36.3 %) dual no-hit; floor 3.35 %; 6.6x standardised on confidence, 8.1x on confidence × length; alignment-normalised 1.38x with CI crossing 1 = null; 96.7 % locally familiar, 63.7 % covered end to end | 3,013 blinded known-fold controls co-searched | Fig 4 |
| R6 | The gap is in the reference set: biome-matched content beats depth | 55.6 % of protein-coding genes carry no KO; 22 % of encoded KOs outside the predictor's set; 26.83 % of dark families absent from 1.387 G sequences (SMAG included, A11 2026-09-27) vs 0.36 % control; desert catalogues +24.33 and SMAG +27.69 points beyond size; depth resolves 35.8 % of control absences vs 5.4 % of dark | matched annotated control in the same DIAMOND pass; rarefaction self-check 0.15 pt | Fig 5 |
| R7 | What the annotated fraction encodes (label) | trace-gas guild moved to P3 (`lifeonair/`, 10 Oct 2026); repair floors; osmolytes 79.5 %; BGC 97.5 %, 85.3 % contig-edge; resistome 37.4 %, arsenic 17.3 % tracking Ca and precipitation | none: prevalence only, no enrichment claim; pathogen census withheld | Extended Data table |
| R8 | Geochemistry structures known and unknown gene content beyond distance ★ | geochemistry R² 0.132, 0.089 beyond space; unique variance geochem 0.106 vs space 0.055 (n = 251); S 54,227 families (91.7 % negative), aridity 23,031 (95.9 %); vanadium null; photolyase 24.3 vs 12.1 % by depth, species sorting; Halobacteriota vs S ρ +0.451 | vanadium block, space partialling, within-phylum and within-species controls, site-paired design | Fig 6 |

Change 1 (2026-09-20, revised the same day on Alejandra's decision): the functional censuses (trace
gas, repair and osmolyte floors, biosynthesis, resistome) stay in the main text as R7 under a LABEL
heading, written as description with every floor named. A control is required for a declarative
claim, not for a census to be reported. The pathogen census is withheld until the E. coli
representative and the blanks are resolved.
★ Change 2: the former heading "Salt, not aridity, filters known and unknown gene content alike"
is withdrawn: the dark families do track aridity (23,031 families, 95.9 % negative, the highest
retention under space partialling of any axis). Geochemistry explains more variance than climate;
aridity is not null. The section is R8 with the heading above.

Order: organisms · how much is unknown and is it real · who owns it · the archaea · what the novelty
looks like · why the references miss it · what the known part encodes · how chemistry sorts known and
unknown alike.

Word budget (3,500 for Introduction, Results, Discussion; 3,700 is the working ceiling before the
true-up trim): Intro 400 · R1 250 · R2 480 · R3 230 · R4 340 · R5 350 · R6 340 · R7 200 · R8 430 ·
Discussion 550. Trim the Discussion and R8 first. Methods uncounted, re-cut from `main_v1.tex`
M.1–M.14; the frozen numbering does not carry over.

Pending decision: `eqmags_poc/04_differential_function` (module-level compartment, trip and chemistry
contrasts, site-paired, phylum sign control) is presented to Alejandra before any promotion; if
promoted it becomes one paragraph in R7.

Hazards while writing: 55.6 % is both the no-KO ceiling and the complete-NHEJ prevalence, never
unqualified in one paragraph; four dark rates on four bases; the cultured control reaches 23
bacterial species and no archaea, so it is never cited from R4; "radiation" is not licensed (two
classes, no phylogenomic claim); "genomes worldwide" is wrong, write "species representatives in
GTDB r232"; never write anthrax; never write "hydrogenase-rich"; comparisons with other studies live
in the Discussion only.

## 4. Introduction, four paragraphs (400 words)

1. Drylands, the hot hyperarid case, and the absence of a genome-resolved survey of this interior;
   what the companion 16S paper mapped and what it listed as out of reach (genes without primer
   bias, function to genome, archaea).
2. The unknown proteome: the fraction depends on reference and rule; the global family catalogues
   report reference absence with lineage but without biome persistence, environment, cultured
   complements or a measured error floor.
3. The two design limits of pooled catalogues: context stripped, error rate unmeasured.
4. Here we: one biome, five expeditions, genomes, cultured strains, geochemistry, a control at every
   layer. Close on the gap sentence.

## 5. Discussion (500 words)

First paragraph states a generalisation the Results did not. Then one paragraph per competing
explanation, each at full strength and then bounded with a number: reference bias, assembly
artefact, rarity or quality of the archaeal genomes, contamination, lineage composition as a
confound. Then what the evaluation does not establish, then the specific next tests. Comparisons
live here only. Source: `DISCUSSION_NOTES.md`.

- What the intersection adds over the taxon-axis, fold-space and global-catalogue papers:
  persistence, a cultivation control, biome-matched reference gain.
- Content beats depth as a recommendation for how soil references should grow.
- The dissolved censuses, one sentence each: CO-weighted where cold deserts are H2-weighted;
  biosynthesis as a target list not a mineable set; a pre-anthropogenic resistome.
- Limits with magnitude or mechanism: capacity not activity; half the representatives below 90 %
  complete; XRF is total elements with no pH; the cultivation control reaches only cultivable
  bacteria; the selection test has no non-coding arm; contamination decision pending.
- Outlook: species in hand as strains; the two archaeal classes as culture and fold targets; strain
  persistence across expeditions as the named next study.

## 7. The dark_function project

Not in this paper (Alejandra, 2026-09-20): no result, ceiling or sentence from `alelopezv/dark_function/`
enters the manuscript. The paper states only that it makes no functional claim about any dark family.

## 8. Target journal

## 6. Display set (six main, five Extended Data)

| item | content | source | state |
|---|---|---|---|
| Fig 1 (R1) | map, timeline, positions, funnel, first-empty-rank ladder | `fig_study_design` | redraw in vector, add ladder |
| Fig 2 (R2) | dark ladder with bases; length vs completeness; expedition persistence; dN/dS vs control; cultured vs MAG paired | `fig_darkfam_reality` (missing) + cultured comparison | rebuild |
| Fig 3 (R3+R4) | variance decomposition; residual dark by phylum; novelty gradient; matched null; archaeal occurrence and depth | lineage darkness `_clean` | new |
| Fig 4 (R5) | structural calibration, four panels, standardised floor | `fig_structural_calibration.pdf` | update |
| Fig 5 (R6) | absence vs control; content vs depth; two ceilings | `fig_specificity_controls` + `fig_refbias_ceiling` | compose |
| Fig 6 (R7) | PERMANOVA blocks and partition; per-axis families; site-paired stress modes | `fig_darkmatter_permanova` + `fig_stress_modes` | compose |
| ED 1–5 | trace-gas guild, BGC and resistome censuses; unbinned reservoir; genome size vs KO richness; contamination diagnostics | existing + `_clean` reruns | as is |

Cut: the endemism set, `fig_picrust2_bias`, `fig_lifeonair_by_phylum` (dead bases), Box 1 and
Box 2 (negatives move into the sections that need them, each with its detection limit), the
energetics model (Supplementary Note 1 stays as a supplementary note only if the trace-gas ED item
needs it).


Not fixed. Draft to Nature Microbiology Article limits (3,500 words, 6 items, 150-word abstract,
headed Results), which keeps Nature and Science formats reachable by trimming rather than rewriting.
Fallback Nature Communications (5,000 words, 10 items). ISME J deprioritised (MAG-only descriptive
work named out of scope; Methods count toward the limit).
