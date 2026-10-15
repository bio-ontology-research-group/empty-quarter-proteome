# Response to corrected Argon H1/H2 review

- H1-MAJ-1: H1_COMPLETENESS.md is generated directly from authoritative h1_final TSV/JSON sources with root_checks/generate_h1_summary.py. Argon already confirmed this repair in its corrected review. No H1 numerical outputs changed.
- H1-MIN-1: Generated H1 summary explicitly states finite-cohort/model conditioning,32unique comparators,dual callers and five missingCheckM2FAA IDs.
- H2-MAJ-1: RESULTS.md now explicitly states75.84%primary-family attrition,arm-specific primary andTier2 eligibility,and thatTier3 lacks positional/adjusted analysis and cannot establish a group effect. Replaced equivalence/generalization wording with failure to establish a difference within eligible families.
- H2-MAJ-2: RESULTS.md explicitly gives877evaluable loci/869bin-contig pairs/617bins/568species and shared counts;112species-disjoint families;four coarse phylum groups including14pooled phyla;85missing annotation bins;and limits of the mobilome keyword screen.
- H2-MIN-1: Moved matching-summary calculations outside the10,000sign-flip loop. IBEX replay job53740945 reran the analysis against copied inputs; all14result files,including the13data outputs and existing manifest,remain byte-identical. Evidence:round7/root_checks/h2_review_replay_verification.json. Original run/results and frozen snapshots are retained.
- H2-MIN-2: Corrected species-resampling comment to568full-run species and its descriptive fixed-outcome nature; RESULTS.md explains that legacypilot_* filenames contain full1,461-family results andrun_full.sbatch is the workflow.

No new positive biological claim was introduced. Manuscript integration remains a separate step after completed H3/H4 review.
