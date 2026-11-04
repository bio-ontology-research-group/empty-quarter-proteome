# EQ-MAGs: structure, narrative and display set

Status: structure in force = the PI's five Results sections (decided 2026-11-04 by Alejandra); the title is the PI's.
Pointers: open items `OPEN_ITEMS.md` · values `CATALOGUE_NUMBERS.md` · changes `CORRECTIONS.md` · text = Overleaf main (`paper/` submodule of empty-quarter-mags); this repo holds the record.
Rule: nothing unregistered enters the paper; a computed but unregistered value stays `\unreg{}` in the source.

## 1. Story in one paragraph

A genome-resolved protein atlas from 551 Rub' al-Khali soil metagenomes links sequence-dark proteins to the
genomes, sites and soil positions that carry them. A tenth of genome-linked clusters lack a UniRef50 match, and
106,042 families recur across samples, sites and expeditions. Thermoplasmatota hold a disproportionate share:
4.4 % of species, 19.8 % of single-phylum dark families. These archaea are widespread and enriched below the
surface. External soil references recover much of the repertoire, more than a quarter stays unmatched, and
structural crops plus recurring genomic neighbourhoods nominate candidates for functional characterisation.

## 2. Title and abstract

Title in force (PI, 20 Oct): **A genome-resolved protein atlas reveals the unknown repertoire of widespread desert archaea**.
Abstract carries: 551 metagenomes; 5,438 species representatives; Thermoplasmatota 4.4 % of species vs 19.8 % of
single-phylum dark families; widespread, enriched below the surface; more than a quarter of recurrent families
unmatched after expanded searches; confident-region crops recover structural matches for most whole-chain misses.

## 3. Results (five sections)

| # | heading (main.tex) | carrying numbers (registry section) | control or caveat → where | fig |
|---|---|---|---|---|
| 1 | A genome-resolved atlas links desert protein diversity to recurrent families | 551 samples, 61 locations, five expeditions; 5,438 reps (126 HQ / 5,312 MQ), 51 phyla, 187 with species; 9.89 % (4,292,446 of 43,423,322); 106,042 families; median 9 sites / 98 samples; 28.9 % assembled, 79.4 % detected in all five expeditions; dN/dS < 1; culture overlap 47 families (§1.1, §1.2, §1.4, PI October) | visit-bundle permutation, culture pairs (22, *E. coli* excluded), AntiFam/shadow → Note S1; dark call basis → Methods | 1, Table 1 |
| 2 | Thermoplasmatota contribute a disproportionate share of unknown families | median dark 31.0 % vs 4.6 % other archaea, 6.2 % bacteria; 236/237 above bacterial P95; matched-null medians 6.5× (218 pairs) / 4.5× (237 pairs); 19.8 % of 77,069 from 4.4 % of species, 4.6-fold; relatives ANI ≥ 95 % 0/237, AAI median 64.8 % vs 75.4 % congeneric, 73.1 % nearest desert rep; PGLS −0.1375 (R4, PI October) | quality-matched null → Methods; PGLS → Note S2, Table S3 | 2, 3a |
| 3 | The unknown archaeal repertoire occurs in widespread subsurface communities | 482 of 550 samples; median 21 species, 0.38 % of reads; 36 of 37 sites (69 pairs) higher below surface, +0.46 pp (0.33–0.86); glycolysis, glyoxylate, peptide transport; precedents Thermoprofundales, Lutacidiplasmatales, Atacama halite type genome (R4) | Holm over four tests, mapped-read and leave-one-expedition sensitivities → Note S3, Table S5; capacity not activity | 3b,c |
| 4 | External soil references recover much of the unknown sequence repertoire | text: 30,479 unmatched on 1.258 G, 27.08 % after MGnify90 2026_07 (registry headline 26.83 % on 1.387 G with SMAG, see B2); occupancy 37 vs 38 of 61; desert +24.33 pp beyond size; SMAG 28.35 vs GEM 9.81 % = +18.54 pp, control +3.10 (§3.1, PI October) | matched UniRef50 control, rarefaction, GEM draws, length strata → Note S2, Table S4, Supp Fig S2, Suppl. Methods | 4 |
| 5 | Structural and genomic context connect unknown sequences to biological hypotheses | 17,221 modelled; 2,838 qualified; 1,807 matched initially; 964 whole-chain unmatched, crops recover 714 (74.07 %); flanking recurrence 353 families, contrast unresolved; GrpE 16/18, DnaJ/GrpE 14/18 (§1.3, PI October) | 3,013 known-fold controls, paired PDB run → Note S2; genus aggregation q = 1 → Note S4 | 5 |

## 4. Discussion (six paragraphs)

1. The atlas ties an unresolved repertoire to widespread subsurface Thermoplasmatota; evidence-explicit candidate selection.
2. Complement to global catalogues; Maeke et al. attribute unknown content in another Thermoplasmatota class to rarity and reference composition; here low abundance with wide occurrence.
3. Soil references recover more dark families at matched depth; habitat, taxonomy and recovery remain confounded; unmatched families as restricted-range candidates.
4. Structural crops shift the task to divergent proteins; neighbourhoods add context but are taxonomically concentrated; chaperone families as assay targets.
5. Limits: coding potential and DNA occurrence only; site independence; partial genes; chimerism screen; expression and proteomics; independent benchmark.
6. Path from sequence discovery to biology: cultivation, targeted expression, biochemistry.

## 5. Display set

| item | content (file) |
|---|---|
| Fig 1 | sampling and genome quality (`fig_study_current.pdf`) |
| Fig 2 | lineage dark fractions, matched null (`fig_lineage_current.pdf`) |
| Fig 3 | family shares, depth contrast, marker sets (`fig_archaeal_ecology.pdf`) |
| Fig 4 | SMAG–GEM, rarefaction, occupancy (`fig_reference_atlas.pdf`) |
| Fig 5 | TM-scores, whole-chain vs crop outcomes (`fig_structure_current.pdf`) |
| Table 1 | analytical units (`tab:cohorts`) |
| Supp Figs S1–S4 | dN/dS and culture sharing · reference recovery and PICRUSt2 · protein maintenance across visits · heliorhodopsins |
| Supp Tables S1–S5 | taxonomy depths · culture gene-call support · PGLS coefficients · matched references · depth contrasts |
| Notes S1–S5 | catalogue quality and coding evidence · reference proximity, context, search controls · ecology and functional profiles · protein-maintenance families · light-response and osmoprotection candidates; plus Supplementary Methods |

## 6. Pending that touches structure

| item | effect |
|---|---|
| B2 SMAG in the headline union | R4 numbers and Fig 4b (26.83 % / 1.387 G vs the text's 30,479 / 27.08) |
| B4 structural excess as a claim | crops-only in main; 6.6× / 8.1× standardised only in Note S2 |
| Contributions statement (B5) | missing section |
| Manuscript snapshot in this repo | blocked on the PI's phrase checkers (read only `manuscript/main.tex`) |

## 7. Superseded

- Eight-section plan (R1–R8): replaced by the five sections above.
- Title "Known phyla, unknown proteome": replaced by the PI's atlas title.
- Census section R7: now Note S3 and the Life-on-Air paper.
- Geochemistry section R8: now Note S3, exploratory.
- dark_function project: still excluded; no result or sentence from it enters the paper.

Hazards while writing: four dark rates on four bases, each with its denominator; the cultured control reaches
bacteria only, never cited for the archaea; "radiation" is not licensed; "species representatives in GTDB R232",
never "genomes worldwide"; never anthrax; never "hydrogenase-rich"; comparisons with other studies in the
Discussion only; 60 = design count of sites, 61 = sampled locations.
