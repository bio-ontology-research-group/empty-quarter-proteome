# Root fact audit of Argon H1/H2 draft

The major finding of erroneous H1_COMPLETENESS.md numerical transcription is valid. Root replaced that report using root_checks/generate_h1_summary.py, reading authoritative h1_final TSV/JSON values directly. No analysis data changed. h2_initial retains the faulty original for traceability.

The review itself contains factual errors that must be corrected before acceptance. This is not a request for a more favorable verdict.

1. Current manuscript main.tex reports551 samples,5438 species-level genomes,51 phyla (see Introduction/Results and CLAUDE.md). Draft opening's522/4768/1622/32/17.5M is unsupported. Read the current full paper and use exact source values, or omit ancillary counts.
2. K is the low-whole-query-TM group and Q\\K the fold-matched comparator; draft second paragraph reverses their interpretations as known-fold K / novel-fold Q\\K. Neither group warrants proven novel-fold labels.
3. The draft blockquote of H1_COMPLETENESS.md is not verbatim; the actual table has four columns, and several quoted values differ from the erroneous table itself. Quote actual text or paraphrase explicitly. The requested correction to authoritative coefficients is nonetheless valid.
4. S6 authoritative coefficients are AAI+.001568,p8.387793e-05;RBH-.145879,p7.485088e-13. S7:AAI+.001550,p7.245799e-04;RBH-.104756,p1.293680e-04. Draft gives different unsupported values. S7 is complete-gene-proxy sensitivity, not harmonized-denominator sensitivity (S3 is denominator sensitivity).
5. Five missing CheckM2FAA genomes in v3_closure_verification.json are48PRr1__SemiBin_112,141,510,539,849; draft020,062,095,111,112 incorrect.
6. H2 tier2 MH raw recurrence difference=.043859,p=.232277 in pilot_conservation_summary.json; draft.0258,p.5067 incorrect. Primary and tier2 excess contrasts remain.035358,p.406859 and.018391,p.561844.
7. H2 dependence summary reports877evaluable loci,869unique bin-contig pairs,8shared contigs,617bins,568species,145sharedbins,153sharedspecies. Draft846pairs/30sharedcontigs incorrect.
8. Species-disjoint subset is112families,51K/61Q,meanexcess.261855/.198430,difference+.063425. Draft76families31K/45Q,difference-.0156 incorrect.
9. Do not claim conserved “non-mobilome” neighborhoods based on limited KO keyword annotations, or that unknown genes form a disjoint “remaining accessory proteome”. These overinterpret the audited bounds.
10. Lambda~.916 is not percent variance explained or proof vertical lineage history dominates all causes. Describe residual phylogenetic covariance under the PGLS model.
11. Current Data/Code Availability text and line numbers must be checked against actual manuscript, not stale placeholder strings/counts. Do not imply external legal/editorial requirements were verified if no journal policy source was consulted.
12. Claimed exact independent re-execution should cite saved scripts/results/logs. Save the actual audit evidence in review directory; distinguish arithmetic recomputation, read-only code review and inherited prior verification. Avoid certifying all estimators from a few spot checks.

Recheck every number and quoted statement, beyond these examples. Preserve original review and verdicts separately; produce corrected review with exact evidence pointers. Keep any remaining major/critical scientific issues and do not promote a positive outcome merely to close review.
