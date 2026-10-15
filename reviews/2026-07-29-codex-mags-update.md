# Empty Quarter MAG manuscript corrections — 2026-07-29

This log records each incremental manuscript or claims-ledger edit and the
source verification that preceded it. Numbers that do not reproduce are not
written into the manuscript.

- Log initialization: the requested path
  `/home/leechuck/Public/software/empty-quarter-workbench/reviews/2026-07-29-codex-mags-update.md`
  is outside the writable workspace and was rejected by the sandbox. This
  workspace copy is being maintained incrementally so no work is buffered.

- Source verification for the fold correction: direct IBEX access failed
  because DNS/network access is unavailable. A same-day concatenation of the
  source `plddt_shard*.tsv` files survived from the preceding run at
  `/tmp/claude-1000/-home-leechuck-Public-software-empty-quarter/cc872bc3-f629-4f01-b2ac-37717694cc04/scratchpad/plddt_all.tsv`
  (SHA-256
  `45bc076e22db314548fc30be11c769b4c54ecf50c47f912f871fe5e9b2861827`).
  It contains 63,654 `ok` models and 281 `too_long` records with `NA`
  pLDDT. Joining the 63,654 models to `RQ/A3/best_afdb_tm.tsv` and
  `best_pdb_tm.tsv`, with no hit treated as novel, independently reproduced
  all six supplied threshold counts and rates.

- Edit 1 — `paper/main.tex`, B1: removed the incorrect 10,007
  dual-database count as a “validated” conservative floor. B1 now introduces
  the explicit-coordinate result as confidence- and length-dependent and
  points to the threshold sensitivity in B1a.

- Edit 2 — `paper/main.tex`, B1a: replaced the single-floor/high-confidence
  interpretation with the reproduced pLDDT sensitivity series:
  15,248/63,654 (24.0%) at ≥0.00; 9,258/50,154 (18.5%) at ≥0.40;
  5,514/33,524 (16.4%) at ≥0.50; 2,836/17,196 (16.5%) at ≥0.60;
  1,164/5,709 (20.4%) at ≥0.70; and 15/187 (8.0%) at ≥0.80. The
  paragraph now explains that the rate does not collapse across thresholds
  with substantial support, whereas the count follows the retained
  population. The same join reproduced median pLDDT 0.440 vs 0.520 and
  median length 182 vs 153 aa for dual-novel vs matched representatives.
  Within pLDDT ≥0.60 it reproduced the five supplied length strata
  (n=9,660/1,839/3,421/1,831/445; dual-novel
  26.3%/5.4%/2.7%/4.8%/4.9%) and exactly 300 dual-novel
  representatives at least 100 aa. The text states both unseparated
  short-protein explanations and the missing length- and pLDDT-matched
  known-fold false-novelty control.

- Edit 3 — `paper/main.tex`, Table `tab:ladder`: removed the 10,007
  “validated novel folds” row and replaced it with all six reproduced
  pLDDT-threshold numerators, denominators, and rates. The caption now labels
  dual-novel representatives as candidates rather than a validated-fold
  count.

- Edit 4 — `RQ/A3/fold_validation_summary.tsv`: replaced the stale
  10,007-floor-era summary with the reproduced confidence and length
  sensitivity values, cohort accounting (63,654 scored and 281 too-long
  records), medians, the 300 confident representatives at least 100 aa, and
  an explicit `NOT_RUN` record for the matched known-fold false-novelty
  control.

- Edit 5 — `reproducibility/claims.tsv`, B1a ledger: removed
  `b1a.novelfloor=10007` and added ledger entries for the scored cohort, all
  six pLDDT thresholds, confidence/length medians, five length strata, the
  300 confidently modelled representatives at least 100 aa, and the unrun
  matched control. These entries point to the updated RQ/A3 summary.

- Edit 6 — `reproducibility/verify.py`, claims registry: mirrored the B1a
  ledger replacement so `verify.py --manifest` will preserve rather than
  overwrite the corrected claims ledger.

- Edit 7 — `reproducibility/verify.py`: removed the obsolete
  `c_a3_novel_floor()` calculation, which operated on the union of TM-hit
  tables and encoded the superseded 10,007 “floor” rather than the
  pLDDT-scored sensitivity cohort.

- Edit 8 — `RQ/A3/result.md`: rewrote the result around the reproduced
  pLDDT-threshold and length-stratified sensitivity analyses. It now explains
  why the stable 16.4–24.0% rate through pLDDT ≥0.70 argues against a purely
  low-confidence artefact, why the 187-model ≥0.80 subset cannot define a
  corrected count, and why the unrun matched known-fold control prevents
  calling the candidates validated novel folds.

- Taxonomy verification: an independent join of
  `gtdbtk.bac120.summary.tsv`, `gtdbtk.ar53.summary.tsv`, and
  `checkm2_all_samples.tsv` reproduced 15,469 classified MAGs; novelty at
  domain/order/family/genus/species of 17/12/297/3,677/10,815 (HQ
  0/3/80/996/2,585); zero at phylum and class; and 3,664 HQ novel
  SeqCode-eligible candidates. It also reproduced that all 17 domain-unplaced
  bins are MQ, with completeness 51.38–65.96%, contamination 2.71–9.74%,
  all PlantRoot (`PR`), and three assigned to the ar53 marker set.

- Non-reproducing supplied detail: the 17-bin source join gives **four**, not
  three, bins from `62PRr2`; `F9PRr2` and `23PRr2` each have three. I will
  not write “three each” into the paper. The same-day workbench `CHANGELOG.md`
  further reports that six of the 17 resemble eukaryotic sequence and four
  gene-dense bins resemble large mobile elements, leaving seven
  uncharacterised. I will therefore label these as unplaced candidate bins,
  not assert that all 17 are cellular deep-branching lineages.

- Edit 9 — `RQ/A3/run.py`: replaced the obsolete printer for the 10,007
  floor-era keys with a printer for the corrected cohort accounting,
  six-threshold sensitivity, medians, five length strata, and the status of
  the matched known-fold control.

- Length-bin boundary check: the supplied `n=445` final stratum reproduces
  only for length **≥400 aa**, not `>400 aa`; eight confident models are
  exactly 400 aa. The corrected RQ/A3 table uses `≥400`, and the manuscript
  wording will be kept consistent with that verified boundary.

- Edit 10 — `paper/main.tex`, B1a: changed the final length-bin label from
  “above 400 aa” to “at least 400 aa”, matching the boundary that reproduces
  `n=445`.

- Edit 11 — `RQ/A3/finalize_tm.py`: replaced the union-of-hit-tables
  10,007-floor calculation and figure with a reproducible join over the
  pLDDT-scored model universe. The script now accepts shard/concatenated
  pLDDT TSVs, treats missing database hits as novel, generates all confidence
  and length sensitivities plus medians, records the unrun matched control,
  and draws confidence/length sensitivity panels.

- Edit 12 — regenerated `RQ/A3/fold_validation_summary.tsv` and
  `fig_fold_validation.{pdf,png}` with the corrected script and the same-day
  concatenated source shards. The executable output reproduced all six
  supplied threshold rows, both medians, all five verified length bins, and
  the 300-model ≥100 aa count; `RQ/A3/run.py` reads the regenerated summary
  successfully.

- Edit 13 — `RQ/A3/methods.md`: replaced the “high-confidence subset”
  method with the scored-universe join, explicit no-hit handling,
  dual-database definition, six inclusive pLDDT thresholds, non-overlapping
  verified length bins, medians, and the missing matched-control limitation.
  Added the current reproduction command for `finalize_tm.py`.

- Edit 14 — `paper/main.tex`, abstract: changed “2.53M ... adopt no known
  fold” to the supported observation that they lacked an AFDB match in the
  initial predicted-3Di screen. This avoids converting a screening no-hit
  into an observed fold.

- Edit 15 — `paper/main.tex`, abstract taxonomy claim: replaced “44 phyla, 0
  novel phyla—novelty is not phylogenetic” with the reproduced rank census
  (12 novel orders, 297 families, 3,677 genera, 10,815 species) and 3,664
  high-quality SeqCode-eligible novel-taxon candidates, while retaining the
  true zero novel-phylum/class result.

- Edit 16 — `paper/main.tex`, A0: replaced “discovery is not taxonomic” with
  the full reproduced shallowest-unclassified-rank census, including all HQ
  counts and the 14,752 bacterial/717 archaeal split. The interpretation now
  states explicitly that zero novel phyla reflects placement in GTDB’s
  high-level scaffold, not absence of taxonomic novelty.

- Edit 17 — `paper/main.tex`, new A0a candidate paragraph: added the
  independently reproduced 17 domain-unplaced bins, all MQ and rhizosphere,
  their rounded completeness/contamination ranges, and the three ar53 bins.
  The paragraph explicitly rejects a novel-phylum claim, notes the
  non-cellular alternative supported by the same-day follow-up audit, and
  requires targeted rhizosphere co-assembly plus cellular-marker recovery
  before interpreting any remaining bin as a deep-branching lineage.

- Edit 18 — `paper/main.tex`, Table `tab:scale`: changed the caption from a
  “taxonomically conventional” interpretation to the GTDB-Tk novelty census
  and added the verified domain/order/family/genus/species and SeqCode rows,
  including HQ counts. The table states that rank categories are mutually
  exclusive.

- Edit 19 — `reproducibility/claims.tsv`, A0 ledger: replaced the
  provenance-only `0 / 0` novelty entry with table-backed entries for the
  bacterial/archaeal split, every reproduced novelty rank and HQ count,
  SeqCode eligibility, and the domain-unplaced bins’ MQ/rhizosphere/ar53
  status and quality ranges.

- Edit 20 — `reproducibility/verify.py`, A0 registry: mirrored the expanded
  taxonomy ledger so manifest regeneration preserves the corrected
  rank-specific and candidate-bin claims.

- Edit 21 — `paper/main.tex`, working title and alternate title: changed
  “Known organisms”/“recognise every microbe” to the supported contrast
  between known phyla/high-level placement and extensive lower-rank plus
  functional novelty.

- Edit 22 — `paper/main.tex`, abstract interpretation: removed the claim that
  novelty lies inside “ordinary taxa” and now states the supported inference:
  high-level placement neither excludes lower-rank novelty nor predicts gene
  content, candidate structural novelty, or community metabolism.

- Edit 23 — `paper/main.tex`, Introduction frame: changed “known organisms,
  unknown biology” to “known phyla do not imply known organisms or
  predictable biology.”

- Edit 24 — `paper/main.tex`, Introduction study summary: replaced the
  “known organisms defeat references” wording with high-level placement not
  making gene content predictable. Immediate inspection found that this patch
  duplicated the opening clause; the next logged edit removes that
  duplication.

- Edit 25 — `paper/main.tex`, Introduction study summary: removed the
  duplicated clause introduced by Edit 24; the paragraph now reads once and
  cleanly.

- Edit 26 — `paper/main.tex`, Discussion: replaced “discovery need not be
  phylogenetic/ordinary taxa” with the supported coexistence of high-level
  placement, extensive lower-rank taxonomy, and distinct functional and
  candidate-structural unknowns.

- Edit 27 — `paper/main.tex`, Honest Negatives box: replaced “conventional
  taxonomy (no novel phyla)” with zero novel phyla/classes alongside
  substantial lower-rank novelty and the explicitly candidate 17 MQ
  domain-unplaced rhizosphere bins.

- PICRUSt2 verification: `per_sample_concordance.tsv` independently gives
  medians NSTI 0.1591, 10,543 predicted KOs, 4,482.5 encoded KOs (reported as
  4,482), and 0.9327 presence-recall over 436 samples.
  `dark_fraction.tsv` gives 31,413,266 of 55,380,378 proteins without a KO
  (56.7227%). The rebuilt 440-sample `threshold_sweep.tsv` gives a
  PICRUSt2/MAG KO-presence richness ratio of 2.39619 at threshold zero and
  1.1148 at a shared $10^{-5}$ threshold.

- Edit 28 — `paper/main.tex`, abstract PICRUSt2 claim: removed 2.4× as a
  biological enrichment, retained the verified 93% presence-recall and 56.7%
  unrepresented-protein bias measurements, and states that 2.4× is a
  zero-threshold artefact falling to 1.11× at $10^{-5}$.

- Edit 29 — `paper/main.tex`, C1: recast the 10,543/4,482 comparison as a
  threshold-sensitive richness measurement rather than enrichment, added the
  verified 93% presence-recall, 2.40→1.11 threshold result and reversal above
  $10^{-5}$, explicitly drops the enrichment claim, and retains the 56.7%
  KO-space ceiling plus signature-function abundance decoupling. Removed the
  unneeded marquee false-positive example from the prose.

- Edit 30 — `paper/main.tex`, C1 heading: replaced the absolute “cannot
  predict their genes” claim with the supported result that 16S prediction
  recovers the KO core but has a reference-defined ceiling.

- Edit 31 — `paper/main.tex`, Table `tab:picrust`: removed “2.4×” from the
  median KO-count row, added the 2.40-at-zero/1.11-at-$10^{-5}$/reversal
  sensitivity row with the correct 440-sample rebuild denominator, and
  removed the marquee false-positive row. The 436-sample NSTI, KO medians,
  recall, abundance decoupling, and no-KO fraction remain.

- Edit 32 — `paper/main.tex`, PICRUSt2 figure caption: replaced the
  “2.4× KO over-prediction” panel claim with threshold sensitivity and the
  retained abundance-decoupling/56.7%-outside-KO-space results.

- Edit 33 — `code/make_picrust2_bias_figure.py`: rebuilt the PICRUSt2 figure
  specification around the corrected evidence. Panel a now plots the verified
  2.40-at-zero to 1.11-at-$10^{-5}$ threshold sensitivity (and the reversal at
  stricter thresholds), panel b retains abundance decoupling without the
  unsupported marquee false-positive framing, and panel d pairs NSTI with
  presence-recall rather than implying a threshold-independent enrichment.

- Edit 34 — `figures/fig_picrust2_bias.pdf` and
  `figures/fig_picrust2_bias.png`: regenerated both figure outputs from the
  corrected script and current source tables (`matched=436`,
  `shared_KOs=7,220`).

- Edit 35 — `paper/figures/fig_picrust2_bias.pdf`: synchronized the
  manuscript’s embedded figure with the regenerated, threshold-sensitive
  PICRUSt2 figure.

- Edit 36 — `reproducibility/claims.tsv`: added separate computed claims for
  the verified PICRUSt2/MAG KO-richness ratios at the zero threshold (2.40)
  and at $10^{-5}$ (1.11), both sourced to the 440-sample threshold sweep.

- Edit 37 — `reproducibility/verify.py`: added source-table computations and
  registry entries for both PICRUSt2 threshold claims, so the ledger values
  are checked directly rather than retained as prose-only corrections.

- Edit 38 — `paper/main.tex`, novelty comparison: replaced the categorical
  “GTDB succeeds/PICRUSt2 fails” contrast with the supported paired result:
  substantial below-phylum novelty, recovery of the KO core, and a
  reference-defined ceiling for the no-KO protein majority.

- Edit 39 — `paper/main.tex`, A1: retained the atmospheric-hydrogenase
  census but explicitly withdrew the aridity-gradient interpretation. Added
  both independent negatives: no within-Empty-Quarter association between
  group-1l share and precipitation ($\rho=-0.1716$, permutation
  $p=0.19610$, 58 sites) and no monotonic cross-biome ordering.

- Edit 40 — `reproducibility/claims.tsv`: added computed ledger entries for
  the within-biome aridity test ($\rho$, permutation $p$, and site count) and
  a table-backed entry recording the negative cross-biome gradient result.

- Edit 41 — `reproducibility/verify.py`: added direct extraction and
  verification of the group-1l/precipitation statistic, permutation
  $p$-value, and 58-site denominator, plus the cross-biome negative in the
  manifest registry.

- Edit 42 — `paper/main.tex`, conceptual-figure plan: changed the proposed
  opening panel from a generic GTDB tree to a GTDB rank-novelty census, so
  the visual cannot imply that “no novel phyla” means no taxonomic novelty.

- Edit 43 — `paper/main.tex`, Honest Negatives box: added both negative
  hydrogenase-gradient tests and explicitly marked the 2.4× PICRUSt2
  enrichment as withdrawn, while retaining the 93% KO-core recall and 56.7%
  no-KO ceiling.

- Verification exception — refined-bin quality census: the requested raw
  source `$M/05_refine/_checkm2_refined/quality_report.tsv` is not present in
  the local workspace (`$M` is unset), and the remote data mount could not be
  reached in this run because network name resolution failed. Workbench
  documentation repeats 14,426 total / 4,225 HQ / 8,209 MQ / 1,992 LQ and
  304 LQ bins at 30--50% completeness and at most 10% contamination, but
  those are derived prose, not the requested per-bin source. In accordance
  with “verify every number against the source files,” these counts have
  **not** been written into the manuscript or claims ledger.

- Edit 44 — `reproducibility/claims.tsv`: corrected the terminal confident-
  model length-bin label from `>400 aa` to `>=400 aa`. The supplied count
  $n=445$ includes eight exactly-400-aa models and therefore does not
  reproduce under a strict greater-than boundary.

- Edit 45 — `reproducibility/verify.py`: mirrored the corrected `>=400 aa`
  boundary in the claims registry.

- Edit 46 — `results/picrust2_bias/SUMMARY.md`: removed the stale 2.4×
  enrichment and “generic reference metabolism” interpretation, added the
  440-sample 2.40→1.11→reversal sensitivity result, and rewrote the bottom
  line around the retained KO-core recall, abundance decoupling, and no-KO
  target-space ceiling.

- Edit 47 — `results/picrust2_bias/README.md`: aligned the analysis README
  with the correction by withdrawing enrichment, reporting threshold
  sensitivity, and retaining only the supported abundance and target-space
  limitations.

- Edit 48 — `code/picrust2_bias_analysis.py`: updated the SUMMARY generator
  so a future analysis rerun reads the threshold sweep and cannot regenerate
  the withdrawn 2.4× enrichment, marquee false-positive framing, or
  “hallucinated metabolism” conclusion.

- Edit 49 — `code/picrust2_bias_analysis.py`: simplified the corrected
  generator’s threshold-sweep denominator label to the directly verified
  440 samples.

- Validation checkpoint: all five modified Python programs compile without
  syntax errors. `python3 reproducibility/verify.py` reports **42 computed
  claims passed, 0 failed**, with 35 table/provenance claims documented. The
  newly added aridity and PICRUSt2 threshold checks all pass exactly at the
  manuscript precision.

- Edit 50 — `reproducibility/claims.tsv`: normalized the zero-threshold
  ratio’s numeric spelling from `2.40` to `2.4` so the hand-maintained ledger
  is byte-for-byte identical to the verifier-generated manifest; the
  manuscript keeps 2.40 to show two-decimal reporting precision.

- Validation checkpoint: `reproducibility/claims.tsv` now has no diff from
  `python3 reproducibility/verify.py --manifest`.

- Edit 51 — `paper/main.tex`, B1 ladder prose: corrected the provenance of
  the 2.53M screen-level count to the initial ProstT5-predicted-3Di/FoldSeek
  AFDB search, removing the implication that all 2.53M underwent explicit
  ESMFold-coordinate validation.

- Edit 52 — `paper/main.tex`, B1 ladder prose: changed the 2.53M screen-level
  label from “structurally novel” to “initial structural no-matches,”
  reserving candidate-fold language for the explicit-coordinate sensitivity
  analysis.

- Edit 53 — `paper/main.tex`, Table `tab:ladder`: relabelled the 2.53M row
  as an initial predicted-3Di no-AFDB-match screen rather than a validated
  no-fold result.

- Edit 54 — `reproducibility/claims.tsv`: applied the same screen-level
  2.53M label in the claims ledger.

- Edit 55 — `reproducibility/verify.py`: mirrored the corrected 2.53M
  provenance label in the manifest registry.

- Edit 56 — `paper/main.tex`, Table `tab:ladder`: corrected the 68,808
  candidate-family screen label from “no AFDB sequence hit” to “no initial
  AFDB predicted-3Di match,” matching the actual selection stage.

- Edit 57 — `reproducibility/claims.tsv`: mirrored the corrected 68,808
  selection-stage label in the ledger.

- Edit 58 — `reproducibility/verify.py`: mirrored the corrected 68,808
  selection-stage label in the manifest registry.

- Build checkpoint — `paper/main.pdf` and LaTeX auxiliaries: rebuilt the
  manuscript with `latexmk -pdf -interaction=nonstopmode -halt-on-error`.
  The build completed successfully (14 pages) with resolved citations and
  cross-references. Existing layout warnings remain, including an overfull
  box in the widened novelty-ladder table; there are no LaTeX errors.

- Edit 59 — `paper/main.tex`, Table `tab:ladder`: constrained the stage
  column to a wrapping width so the corrected, explicit selection-stage
  labels fit the page instead of extending into the margin.

- Build checkpoint — rebuilt `paper/main.pdf` after the table-width fix.
  LaTeX again completed successfully (14 pages). The post-build warning
  filter command itself could not find `rg` in that nested-shell invocation;
  this did not affect the successful manuscript build.

- Build verification: checking `paper/main.log` from the repository shell
  confirms that the novelty-ladder overfull-box warning is gone and that
  there are no undefined references/citations or LaTeX errors. Three
  pre-existing overfull paragraphs remain elsewhere in the manuscript.

- Additional PICRUSt2 source check: `marquee_functions.tsv` gives abundance
  Spearman $\rho=-0.123$ to 0.417 (median 0.086) across the ten signature
  functions present in at least half of MAG samples. This does not support
  the blanket “all decoupled / approximately 0 to -0.12” wording.

- Edit 60 — `paper/main.tex`, C1: replaced that blanket abundance claim with
  the verified heterogeneous range ($-0.12$ to 0.42; median 0.09) and states
  plainly that it is not a uniform failure.

- Edit 61 — `paper/main.tex`, Table `tab:picrust`: replaced the stale
  “approximately 0 to -0.12” abundance row with the verified $-0.12$ to
  0.42 range and 0.09 median.

- Edit 62 — `paper/main.tex`, PICRUSt2 figure caption: changed
  “signature-function abundance is decoupled” to the supported description,
  “abundance agreement is heterogeneous.”

- Edit 63 — `reproducibility/claims.tsv`: added computed minimum, median,
  and maximum signature-function abundance correlations to cover the newly
  corrected manuscript range.

- Edit 64 — `reproducibility/verify.py`: added direct computations of those
  three correlations over the ten functions with MAG prevalence at least
  0.5 and registered them in the verifier manifest.

- Edit 65 — `results/picrust2_bias/SUMMARY.md`: replaced uniform abundance-
  decoupling language with the verified function-dependent range and
  softened the bottom line accordingly.

- Edit 66 — `results/picrust2_bias/README.md`: made the same correction in
  the analysis overview, retaining the 56.7% no-KO result as the hard
  reference-space limit.

- Edit 67 — `code/picrust2_bias_analysis.py`: updated the generated summary
  to calculate and report the signature-correlation minimum, median, and
  maximum and to avoid regenerating a uniform-decoupling claim.

- Edit 68 — `code/picrust2_bias_analysis.py`: aligned the nearby code comment
  with the corrected “abundance agreement” interpretation.

- Edit 69 — `code/make_picrust2_bias_figure.py`: changed panel b’s title
  from uniform decoupling to heterogeneous agreement.

- Edit 70 — `figures/fig_picrust2_bias.pdf` and `.png`: regenerated the
  corrected figure after the abundance-agreement wording fix
  (`matched=436`, `shared_KOs=7,220`).

- Edit 71 — `paper/figures/fig_picrust2_bias.pdf`: synchronized the
  manuscript copy with the regenerated figure.

- Edit 72 — `code/picrust2_bias_analysis.py`: removed the last stale
  “abundance decoupling” header comment.

- Edit 73 — `paper/main.tex`, Part C heading: replaced the residual “known
  organisms defeat references” framing with “high-level placement does not
  guarantee reference coverage.”

- Edit 74 — `paper/main.tex`, Discussion: replaced categorical reference-
  database “failure” wording with the supported conclusion that reference-
  defined ceilings persist despite high-level placement.

- Edit 75 — `RQ/A3/run.py`: added the verified count of 300 dual-novel
  representatives at pLDDT at least 0.60 and length at least 100 aa to the
  printed summary, matching the regenerated TSV and manuscript.

- Validation checkpoint: `python3 RQ/A3/run.py` prints the full six-threshold
  series, both medians, all five length bins, the 300-model count, and
  `NOT_RUN` for the matched known-fold control without error.

- Edit 76 — `code/make_picrust2_bias_figure.py`: removed the obsolete
  full-matrix loads used only by the withdrawn overprediction panels. The
  corrected figure now depends only on its plotted summary/source tables,
  avoiding a slow and unnecessary recomputation and external amplicon path.

- Edit 77 — `figures/fig_picrust2_bias.pdf` and `.png`: regenerated from the
  streamlined script in 1.5 seconds; it still reports all 436 matched
  samples.

- Edit 78 — `paper/figures/fig_picrust2_bias.pdf`: synchronized the
  manuscript figure after the streamlined regeneration.

- Validation checkpoint: the source and manuscript PICRUSt2 PDFs have the
  identical SHA-256
  `a21b45e53a7817a8f3c4559df2c8d35a6efea575e4866a1d36ab761a50559e4e`.

- Edit 79 — `reproducibility/claims.tsv`: added a computed 440-sample
  threshold-rebuild denominator, sourced to the 441-line (header plus 440
  samples) per-sample rebuild table.

- Edit 80 — `reproducibility/verify.py`: added and registered the direct
  per-sample rebuild row count.

- Edit 81 — `code/picrust2_bias_analysis.py`: replaced its hard-coded
  “440-sample” summary label with a row count from
  `rebuild/per_sample_rebuilt.tsv`.

- Validation checkpoint: all modified Python files compile; the claims
  ledger again has no diff from `verify.py --manifest`; and the final
  verifier run reports **46 computed claims passed, 0 failed**, plus 35
  documented table/provenance claims.

- Final manuscript build: `latexmk -pdf -interaction=nonstopmode
  -halt-on-error main.tex` completed successfully and produced the updated
  14-page `paper/main.pdf`.

- Final QC: neither repository reports whitespace errors under
  `git diff --check`; the final LaTeX log has no undefined citations,
  undefined references, or errors. Four pre-existing overfull-box warnings
  remain outside the corrected novelty-ladder table.

## Handoff

Completed: novel-fold threshold/length correction; lower-rank taxonomy and
domain-unplaced candidate framing; PICRUSt2 threshold correction; both
negative group-1l aridity results; manuscript/figure/source-summary alignment;
and claims-ledger verification.

Intentionally omitted: the 14,426-bin refined CheckM2 census, because its raw
per-bin source was unavailable; and the supplied “three bins from 62PRr2,”
because the source join gives four. The requested public-workbench log target
remains unwritable in this sandbox; this workspace log is the complete copy
to transfer there.
