# Current status after round 4

Four external reviews are complete: three Opus 5.5 rounds and the requested
Gemini 4 Argon High review through Antigravity. Nature Communications remains the
selected stretch target. Root independently checked the fourth report because it
contains both useful findings and inaccurate quotations, counts and methods.
The full report and point-by-point dispositions are in `round4/`.

The current manuscript is `manuscript/main.pdf`, 28 pages. Round 4 corrected:

- PICRUSt interpretation: all 10,543 reference orthologues are predicted in
  419/448 samples, so presence-recall largely describes reference overlap.
- Metabolic comparisons: full-lineage prevalences and the 181 matched pairs are
  reported separately, with explicit marker rules and denominators.
- Depth comparisons: paired-site medians and the distinct mapped-fraction result.
- Soil-position analysis: 547 samples with recorded positions; 20,136 enriched
  and 41,046 depleted families, with updated phylum and detection summaries.
- Partial-correlation inference: corrected degrees of freedom; 57/224 tests pass
  5% FDR, one fewer than before. Correlation coefficients are unchanged.
- SPIRE anchor provenance, figure legibility, citations and bibliography rendering.

IBEX job 53695688 completed the new analyses in 40 seconds. Auditing all 5,438
rRNA GFFs found no changed calls under the corrected complete-hit preference;
126 high-quality MIMAG representatives remain supported. Original upstream tables
and code are preserved alongside corrections and input/source checksums.

Validation: the full compact driver passes 49 prior plus 40 round-4 checks,
verifies 239 retained input/source hashes, regenerates five figures and builds the
PDF cleanly. Root also checked key-level culture/locus/family joins and structural
zero-hit classifications. The final PDF has no tiny-text flags; its three white
heatmap counts are legitimate. The historical verifier still has 31 unreconciled
manuscript gates, with 134 numerical checks passing. No whole-paper pass is claimed.

The literature library remains 77 distinct DOI-verified fulltexts. This round did
not add a literature-search campaign. The earlier 1,494 large-file hash audit
(292.8 GB) remains retained; it was not redundantly rerun.

The curated working package is mirrored under:
`/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p4_review_20261013/review_package/`.
It reproduces the stated compact scope, not the full raw-to-final pipelines.

Before submission: resolve 78 missing raw-file paths and 302 permission failures,
duplicate-run merge provenance, laboratory XRF metadata, historical claim-registry
migration, ENA release/reviewer access, and MAG/code-data deposition with an archive
DOI. Companion manuscripts in preparation must accompany the editorial package.
No submission, public release, Git push or external message was made. Root’s final
round-4 revisions were verified, but have not had a fifth external review.

Final IBEX verification completed successfully in job 53696382 (four seconds):
239 input/source hashes, all 89 compact checks, the independent cross-table audit,
and all 52 final manuscript/figure/script hashes passed. Receipt:
`round4/ibex_package_verification.log`.
