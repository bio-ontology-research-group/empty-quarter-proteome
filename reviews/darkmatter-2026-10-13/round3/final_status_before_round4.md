# Final status: Nature Communications revision

All three requested independent reviews completed using `claude-opus-5-5`.
Root revised the paper after each report. Antigravity/Gemini Flash guided the
literature expansion; 77 distinct DOI-verified full-text papers are retained,
with inaccessible or invalid leads recorded separately.

The final draft is `manuscript/main.pdf` (28 pages; abstract 148 words). It presents
the desert genome resource alongside evidence that gene calling, reference
retrieval and comparison populations affect conclusions about dark proteins.
The Nature Communications target remains a stretch; the draft is not ready for
submission until the remaining prerequisites below are resolved.

The most consequential changes are:

- Full MIMAG assessment gives 126 high-quality genomes, distinct from 2,259 genomes
  meeting the CheckM2-only near-complete thresholds.
- Most artifact-flagged families also have dN/dS below one. The paper reports that
  calibration and avoids treating the statistic alone as evidence of coding status.
- Among DNA-supported calls, 46.8% of dark calls and 0.90% of annotated calls lack
  a Bakta CDS covering half the aligned region. Only 65 culture-comparison dark
  calls map into 47 admitted recurrent families; these populations are separated.
- Structural retrieval failures are distinguished from below-threshold alignments.
  The paper no longer treats query-score failures as evidence of novel folds.
- Environmental results are explicitly exploratory supplementary context.

Validation: 49 current revision claims, 228 input/source checksums, regenerated
Figures 1–5, a byte-identical environmental-table rebuild, and a clean LaTeX build.
The final PDF was inspected, including the corrected Figure 2 layout. The only
hidden-text flags were legitimate white sample-count labels on a blue heatmap.
The historical whole-paper verifier still fails stale manuscript gates and is
not represented as passing.

IBEX job 53680819 additionally hashed 1,494 unique large analysis files (292.8 GB)
with zero errors. The completed clean-analysis inventory has 447,858 checksum
records and 14 broken temporary links; duplicate link references are included.
The separate upstream inventory remains partly metadata-only.

A curated working package is on IBEX at:
`/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p4_review_20261013/review_package/`.
Original large inputs remain at the paths in the provenance inventories. The
compact package is reproducible within its stated scope; it is not a complete
raw-to-final archive or a public deposition.

Before submission, resolve the 78 missing raw-file paths and 302 permission
failures, duplicate-run merge provenance, laboratory XRF units/detection limits
and zero coding, historical claim-registry migration, ENA release or reviewer
access, and MAG/code-data deposition with an archived DOI. The two companion
manuscripts in preparation also need to accompany the editorial package.

No submission, public release, Git push or external message was made. The reports,
root responses and frozen submitted versions are retained in the three round
folders. Final edits were made after round 3; no fourth review is implied.

IBEX package verification completed successfully (job 53681298, three seconds):
228 input/source checksums, all 49 revision checks, and the final manuscript,
figure and analysis-script hashes match the local package. The receipt is
`round3/ibex_package_verification.log`.
