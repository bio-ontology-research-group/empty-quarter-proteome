# Root response to Gemini 4 Argon

The fourth review ran through `agy --model gemini-4-argon-high --effort high
--dangerously-skip-permissions`. The CLI listed that exact model before launch;
`argon_result.json` records SUCCESS and 1,034.96 seconds. The receipt does not
include a separate provider-resolved model field. The unedited report and audit
scripts/outputs are retained. Root independently checked each actionable finding.
Nature Communications remains the selected stretch target, not an assured fit.

## Reviewer reliability and rejected assertions

The report contains substantial inaccuracies despite useful computational audits.
Its opening gives 17,677 bins, 132.4 million calls, 7.5% darkness, 235 Poseidoniales,
and multi-database absence as the definition of darkness. These are not the current
paper: it reports 14,426 refined bins, 43,423,322 genome-linked 90%-identity clusters,
4,292,446 dark representatives (9.89%), and a UniRef50 no-hit definition. Its claims
of purifying selection established after adjustment and multidomain structure as
an established cause overstate the manuscript. Controls are PDB-derived; the review
incorrectly says they come from desert MAGs. Neither the claimed MIMAG abstract
passage nor the claimed Tables 2/3 and Figure 4d exists. Several quoted passages and
line numbers differ materially from the frozen source. Raw-file blockers are 78
missing paths and 302 permission failures, not 23 and 28. These assertions were
rejected and were not introduced into the manuscript.

## Main comments

1. **PICRUSt saturation: accepted with corrected interpretation and numbers.**
   The exact retained table has 419/448 samples with all 10,543 KOs emitted and
   420 with every within-reference MAG KO predicted. Added these counts, median
   Jaccard 0.634 and median 3,077.5 predictor-only KOs, and clarified reference-set
   overlap in Results and Figure 5. Rejected the review's “false positives” and
   “specificity is zero” conclusions: incomplete MAG recovery is not a complete
   negative truth set. Recall equals reference overlap in the 420 stated samples,
   not automatically every sample. The review's minimum 9,683 and 438 samples above
   10,500 are wrong; the table gives 1,948 and 421. Spearman 0.839 remains a separate
   abundance-rank comparison.

2. **Metabolic matching: accepted; several proposed methods rejected.**
   Separated 237-genome prevalences from the 181 matched Thermoplasmatota versus
   Halobacteriota pairs. Reported glyoxylate 139/181 versus 63/181; carbon-skeleton
   markers 3/181 versus 165/181; catalase 49/181 versus 153/181; sugar transport
   4/181 versus 132/181. Archaellum-associated markers occur in 236/237, and the
   wording now identifies marker calls. Methods give the one-point caliper,
   181/216/237 pair counts, module and marker rules. Actual module threshold is
   two thirds, not 0.75; other archaea include 354 genomes, not 216 before matching;
   Halobacteriota total 269, not 229. Oxoglutarate markers use any of K00174/K00175,
   and CRISPR any of seven KOs, not the invented multi-KO rules in the review.
   Cytochrome oxidase explicitly requires both K02274/K02275, distinguishing it
   from the separate single-subunit screen. Removed stronger metabolic wording
   unsupported by marker presence alone.

3. **Depth reporting: accepted.**
   Replaced unpaired medians with paired-site medians 0.72% subsurface and 0.27%
   surface across 45 sites, 40 higher below, p=4.177e-7. Normalised comparison is
   38/45, p=7.936e-7. Removed “unchanged” and “indifferent.” The review's quoted
   299/548 detection passage is absent from the frozen manuscript; the existing
   481/549 result is preserved.

4. **Ecology bugs: accepted and recomputed, rather than only disclosed.**
   IBEX job 53695688 restricted the binary soil-position contrast to 547 samples
   with recorded positions and propagated membership to downstream phylum and
   detection summaries. It now tests 101,897 families, with 20,136 enriched and
   41,046 depleted. Updated all dependent percentages in Note S1 and corrected
   Myxococcota_A versus Myxococcota. The excluded samples are 62PRr2, 64Dr3,
   64PRr2, 64Sr2, not Camel1/HotSpring/OilSpill/WaterWell; minimum presence is 25,
   not the review's 10. Root retained the existing 94.1% geography-conditioning
   statistic on its documented 547-sample cohort. The review's phylum numerator
   examples also confuse all tested families with phylum-attributed denominators.
   Corrected the two logged PERMANOVA core-model sample counts to 251/456 without
   altering estimates or rerunning permutations. Corrected partial-Spearman
   residual degrees of freedom and reran 224 tests: 57 pass 5% FDR versus 58
   previously; sulfur/pigment occurrence changes q=0.04981 to 0.05117. No effect
   estimates change. Existing Methods already identify the 548-sample stress
   cohort; the reviewer’s claimed switch to 550-sample UV results is unsupported.
   Original outputs remain frozen alongside the corrected versions.

5. **SPIRE anchor: accepted as a provenance footnote.**
   Methods and Figure 5 now disclose the 271,000,000-sequence expectation anchor
   versus the actual 269,496,580-sequence index. The review's 74.87% observed absence
   and 13.7-point residual are wrong: the retained table gives 39.90% observed,
   61.20% expected and +21.30 points. We preserve this explicitly approximate
   analytical anchor; no new random subsampling analysis is implied.

6. **Production: partially accepted.**
   Enlarged Figure 2 numeric labels, added headroom to Figures 3c/4c, inspected the
   renders, changed bare numeric citealp calls to bracketed citations, repaired
   embedded “others and” author lists and the Soding umlaut, removed internal
   bibliography notes, and protected verified title capitalisation. Added NAR
   volume 54, issue 18, article gkag938 for Aplakidou from the retained DOI-verified
   primary XML. Rejected nonexistent Supplementary Table S3/Extended Data Table 1
   references: neither string occurs in the submitted manuscript, and the claimed
   quoted availability paragraph is absent. Several allegedly missing reference
   keys are also absent; Makhalanyane already had volume and pages. No references
   or claims were added to satisfy these fabricated findings.

## Minor comments

1. Removed the orphaned cross-culture collection search procedure; retained the
   executed paired-species and reciprocal protein-hit methods used in Results.
2. Explicitly stated that the genome CoverM invocation sets no read-identity or
   aligned-read-fraction threshold; did not infer new filter values from a newer
   software version.
3. Added non-partial rRNA criteria and 80/900/1,900 nt floors. Audited all 5,438
   original GFFs on IBEX: preferring a qualifying complete hit changes zero calls.
   Preserved a corrected parser and regression fixture for the real edge case of
   a longer partial hit masking a complete hit. MIMAG counts remain 126/5,312.
4. Documented the extreme source Fe values and site-44 median without accepting
   the review's unverified percent units or proposed column-swap explanation.
   Raw laboratory metadata remain necessary; no values were silently corrected.
5. Removed the stray space before the supplementary sentence-ending period.

## Verification and boundaries

The extended compact driver includes 49 prior checks plus 40 round-4 checks.
The independent root key-level audit checks culture/locus/family population joins
and structural zero-row identities. The input ledger now pins 77 retained inputs
plus 162 historical upstream source files (239 total); it preserves older hashes.
The rRNA GFF/source/input hashes from the new remote analysis are retained separately.
The historical verifier still has 31 unreconciled manuscript gates, with its 134
numerical checks passing; no whole-paper pass is claimed. Final compact and remote
verification receipts are retained alongside this response.

Raw-file access/location, duplicate-run merge provenance, laboratory XRF metadata,
historical registry migration, public ENA/reviewer access and MAG/code-data archive
release remain submission prerequisites. No submission, release, push or external
message was made. Root's revisions after this fourth report have been tested and
audited; they have not undergone an additional external review.
