# Current status after round 5

Five external reviews are complete: three Opus 5.5 rounds and two Gemini 4 Argon
High reviews through agy. Nature Communications remains the selected stretch
target. The fifth review followed Flash scouting and three new IBEX analyses.
Root checked the recommendations, rejected unsupported statements and revised the
manuscript. The report, response, model receipts and ranked extension plan are in
round5/. Argon reviewed the frozen manuscript plus new results; the final prose
has not had another external review.

Current manuscript: manuscript/main.pdf, 30 pages. New results:

- Recalling the same 22 culture assemblies with the original MAG caller reduces
  dark calls without a CDS overlap from 344/735 to 235/735. Gains are concentrated
  across pairs; the admitted-family subset remains 18/42. This establishes caller
  sensitivity at fixed DNA-supported loci, while coding status remains unresolved.
- Five-expedition detection is 79.38%, against an individual-label null mean of
  82.74%. A reviewer-driven whole-visit sensitivity gives 80.54%, with a central
  95% range of 78.95–82.11% and lower-tail p = 0.11. Neither supports excess
  detection; the apparent deficit depends on the permutation unit. Both are reported.
- Paired functional profiles have a correlation advantage of 0.0185 over the
  expedition/position null, falling to 0.0023 within site and visit (146 samples).
  The full comparisons are in Supplementary Note S2. At the highest abundance
  threshold, 446/448 MAG profiles retain no KOs and have undefined recall.
- Other revisions disclose structural-confidence imbalance, compare the two dark
  Pfam subsets, clarify contaminant flags and a limited exclusion sensitivity,
  repair two dangling table references and culture-group counts, and complete
  seven citations from retained DOI-matched sources.

The further hypotheses are ranked in round5/EXTENSIONS.md. Root recommends testing
Thermoplasmatota reference distance, shared-gene coverage and lineage-specific gene
content with quality and ancestry controls, followed by targeted conserved gene
neighbourhoods. Matched-reference searches and domain-level structural calibration
are additional directions. Argon's claims of established accessory expansion,
136 expressed targets, an 8.2-fold soil-reference gain, or equivalence to the PDB
control baseline are unsupported and were rejected.

IBEX jobs 53697241 and 53697281 completed the initial analyses; job 53698522 completed
the whole-visit sensitivity. Scripts, seeds, outputs and input hashes are retained.
The compact driver verifies 239 original input/source hashes and passes 49 prior,
40 round-4 and 159 round-5 checks, regenerates figures and builds the PDF. The final
integrity scan flags only legitimate heatmap counts 69, 90 and 61. Changed pages
were visually inspected. Historical verification remains separate: 62 COMPUTE and
72 FROZEN checks pass, six TABLE sources are present, 25 PROVENANCE entries are
pinned, and 31 manuscript gates remain unresolved. No raw-to-paper pass is claimed.

The 77 DOI-verified fulltexts and earlier 1,494-file checksum audit (292.8 GB) remain
preserved. The curated package is mirrored at:
/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p4_review_20261013/review_package/.
Large matrices, assemblies and new culture GFFs remain at the workbench paths
recorded in round5_extensions/README.md and its input manifests.

Before submission: resolve 78 missing raw-file paths, 302 permission failures,
duplicate-merge provenance, XRF laboratory metadata, historical registry migration,
ENA/reviewer access and MAG/code-data archival release. Abundance-weighted
contaminant effects remain an uncompleted sensitivity. No submission, public
release, push or external message was made. The round-4 status is preserved in
round5/prior_status.md.

Final IBEX verification passed in job 53699101 (eight seconds): 239 original
input/source hashes, all 248 compact checks, the independent cross-table audit,
61 final artifact hashes and 101 remote extension input/generated-GFF hashes.
The pending job was moved from batch to the debug partition for this short CPU
verification. Receipts are in round5/ibex_package_verification.log and
round5/ibex_verification_job.txt.
