# EQ-MAGs restructure: what to do next

Written 2026-09-17 from the seven dossiers in `notes/dossiers/` and the blocker table in
`PAPER_STRUCTURE_PROPOSED.md`. Numbers in brackets are that table's blocker ids.

## Where corrections live

No separate ledger. Three files, each with one job, as the repo already worked:
`CATALOGUE_NUMBERS.md` = what the paper must say (corrected values are written straight in, with the
superseded one named in the same cell) · `CORRECTIONS.md` = dated old → new and why ·
each workbench `_clean/README.md` = method detail for that analysis. Proof-of-concept numbers stay
in their own README and never enter the registry until you decide they are in the paper.

**PoC promotion, pending your decision.** `alelopezv/eqmags_poc/01_dark_by_lineage` (lineage
darkness + the archaeal matched null) and `05_rakpr_overlap` (cultured comparison) now hold
paper-grade results but sit outside the workbench `_clean` convention. They are not moved yet
because 19 scripts and 19 note references carry those paths, and the results may not be used. The
moment you decide either is in the paper, promote it to
`06_darkmatter/{lineage_darkness_clean,cultured_comparison_clean}` and repoint the paths in one
pass. Until then their numbers stay out of the registry, so nothing is inconsistent.

## The sequencing question, answered first

**Neither "fix everything then write" nor "write then polish". Do a short fix sprint on the items that
change a number or a claim, then write while the long compute runs in parallel.**

Reason: prose written on numbers that are about to move is exactly what produced the current draft, and
rewriting it costs more than the fixes do. But most defects found are disclosure items that cost nothing
to handle while drafting, and the two big computations only decide two section titles, so they must not
block the other six sections.

| phase | what | how long | blocks what |
|---|---|---|---|
| 0 | List A, the number-changing reruns | about 1 day | writing R2, R4, R5, R6, R8 |
| 1 | Write the manuscript from the dossiers, list B handled inline | 5 days | nothing |
| 1 (parallel) | List C long compute | 1 to 2 weeks | only R5's title and R4's depth |
| 2 | True-up: registry rows, claim registry, figures | 2 days | submission |

---

## Other open items, not in lists A to C

Carried over when the blocker table was retired from `PAPER_STRUCTURE_PROPOSED.md`.

- Shotgun blanks SNTC and SPositive exist under trip4 controls and were never profiled against the
  5,438. Hours of compute, and the contamination argument leans on them.
- `figures/fig_darkfam_reality.pdf` is referenced in the old manuscript and does not exist.
- M.1 field facts are still owed: XRF workflow and units, aridity source, samples collected.
- `env_matrix.tsv` cannot be rebuilt, its sources are gone. The salinity proxy is an undocumented sum
  of Na, Cl and S, the aridity index is de Martonne rather than UNEP, and there is no pH anywhere.
- Dark-on-mobile-element sign is unresolved between the PI's pooled odds ratio of 0.62 and the
  cultured 2.3x. Print neither until reconciled.
- The selection result has no upper arm: no shuffled or intergenic set near neutrality, so the
  distance of 0.425 from 1.0 is uncalibrated. Cheapest upgrade to R2's strongest claim.
- The dark-family control was admitted on member count alone while the dark arm cleared three gates;
  and family admission uses the unrefined bins while phylum attribution uses the refined catalogue.
  Both need one Methods clause.
- R8's unanswered attack: 87 to 96 % of associations are negative, which is what a presence rule
  produces if biomass falls toward the salty end. No depth or biomass covariate is in any model.
- Literature to ingest into READING_NOTES and the bibliography: AFESM v2, Schaeffer 2026, SMAG
  (Ma 2023), Maeke 2025, Diamond 2022, Rinke 2019.

## List A. Mechanical, my call, no approval needed

All are reruns of existing scripts on existing inputs, or arithmetic. None invents a method. Each changes a
number that prose would otherwise be written on. I report every result before it enters any text.

1. Rerun the geochemistry block with sulfur, chlorine and sodium included, and drop the variance-partitioning
   claim that the code never computed. About 4 minutes. [17, 18]
2. Rerun the occupancy join on the corrected union so the sites result leaves the superseded basis.
   Minutes. [36]
3. Standardise the structural control floor to the candidate confidence and length distribution, and report
   the stratified multiplier in place of 10.8-fold. One pass over two files on disk. [29]
4. Recompute the cultured comparison over the 23 distinct genomes instead of the 132 pseudo-replicated
   pairs. Minutes. [21]
5. Rebuild the archaeal contrast as a proper matched null on domain, completeness, size, abundance and
   prevalence. Minutes, inputs all on disk. [16, 23]
6. Rerun the darkness-versus-genome-size test on the native key, since it is the last analysis still on the
   superseded rename basis. Under a minute. [dossier 07]
7. Fix the two arithmetic slips: the ceiling is 58.3 not 58.4 per cent, and the predictor does not emit a
   constant number of families in every sample. [39]
8. Register the facts that exist but are unrecorded: sequencing depth per sample, and the two named samples
   behind the drop from 551 to 549. [10, 35]
9. Correct the stale records: the pipeline guide's provisional structural header, its wrong presence-rule
   citation, and the flat claim that chimerism screening never ran when it did run on the ten order-gap
   candidates. [34, 39]

## List B. Needs your approval, because it changes a claim or a registered number

Each is a one-line decision. None needs compute.

1. **The selection sentence.** State plainly that dark families are less constrained than the matched
   control. The section currently prints both medians and lets the reader infer, while the abstract already
   says it. [26]
2. **The control's admission rule.** Disclose in Methods that control families were admitted on member count
   alone while dark families also had to clear three samples and two sites, or rebuild the control under the
   identical rule. Disclosure is honest and cheap; rebuilding is cleaner and costs hours. [24]
3. **The matched-control wording.** The shipped query files are not length-matched even though the
   construction was. One Methods clause. The error runs in your favour. [38]
4. **The MIMAG sentence.** Either run ribosomal and transfer RNA detection on the 5,438, or say that only
   completeness and contamination were applied. The data currently do not support the claim as written. [32]
5. **The carbon-monoxide figure.** 22 per cent is a floor, with 60 per cent of proteins left unresolved.
   Say so, improve the calling, or drop the CO-weighted framing. [41]
6. **The aridity and dark-supply negatives** in the trace-gas section are on the dead 15,469 basis and are
   cited in the present tense. Recompute or mark them historical. [42]
7. **The ectoine triangulation.** Methods say three methods; two is accurate. [43]
8. **The undisclosed structural selection step**, where models already built by the PI were preferred.
   Disclose in Methods. [30]
9. **The E. coli representative and the pathogen census.** Keep with a contamination annotation or remove,
   which would move the locked 5,438; and fix the name matching so the count is 5 rather than 2. Never write
   anthrax. [1, 2]
10. **The sample cohort**: pin at 548 and say why, or rerun at 551 now that the failed profiling tasks were
    repaired. [37]

## List C. Complex decisions, yours

1. **Does the archaeal section lead the paper?** It is the strongest unreported finding and the cheapest to
   secure, but the contrast needs its matched null first (A5) and the cultured control cannot support it.
2. **Structural armour: how far?** Confidence standardisation (A3) is mandatory and cheap. Pfam on the
   candidates is installed and half a day. Domain-level calling needs a tool installed first, and the
   structure-level database search is not runnable here at all. Decide how much armour before submission, and
   accept that the section title depends on the answer.
3. **Target journal.** Nature Microbiology as an Article, with Nature Communications as the fallback at
   effectively the same time to acceptance. Nature and Science are not realistic without experimental
   validation. The heading and length rules differ, so this decides the word budget before writing.
4. **Strain persistence as a second paper.** Feasible and unclaimed, about 1,900 species qualify, but it needs
   a full re-map of 551 samples. Confirmed today that no measured coverage table exists to shortcut it.
5. **Authorship and the byline**, still unresolved in the manuscript.
6. **What to do about the claim registry**, which has no rows for two sections and holds withdrawn values for
   others. This is the PI's to rebind, but someone has to raise it. [8, 40]
7. **The biosynthesis and resistome censuses have no control and no size normalisation.** Normalise, or
   present them as description only. [44]

## Not doing

- Nothing from the exploratory folders enters the paper without a registered rerun.
- Nothing from the companion annotation project enters the paper beyond two sentences stating its ceilings.
- The environment matrix cannot be rebuilt, its sources are gone. Work with what is registered and say so.
