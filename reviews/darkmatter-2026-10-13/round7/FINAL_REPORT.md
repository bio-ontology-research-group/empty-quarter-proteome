# Completed H1–H4 audit, integration and two review/refine rounds

The revised paper is `manuscript/main.tex` and `manuscript/main.pdf` (35 pages).
Nature Communications remains the target. Supported results are integrated into
Results, Methods, Discussion, the 150-word abstract, Supplementary Note S1 and
Figure 4. Supplementary Notes are in first-citation order.

## Scientific judgement

- **H1:** phylogeny-conditioned associations across 237 Thermoplasmatota representatives are supported. Shared-protein coverage and nearest-comparator identity remain conditional associations, with six sensitivity specifications and explicit comparator/tree limitations.
- **H2:** callable multi-species families show recurrence above their positional background. The adjusted between-group contrast remains inconclusive (difference 0.035358; p=0.406859). Differential window eligibility and dependence remain explicit.
- **H3:** equal-depth, exact-length-matched reference searches recover 28.35% of dark queries from SMAG and 9.81% from GEM (18.54 percentage-point difference). Fixed SMAG sampling, GEM draw sensitivity, caller/rematching effects and unmatched taxonomic composition bound the interpretation.
- **H4:** pooled qTM endpoints support crop recovery for 714/964 whole-chain-unmatched dark models (74.07%); GPU sensitivity gives 731/961. Cropping jointly changes sequence content, length-dependent scoring and search opportunity. Defective alignment-TM subdivisions and heuristic rescaling are excluded from biological interpretation.

## Review trail and refinements

Round1: `postintegration_round1/review.txt`, `verdict.json`, executed verification,
and `root_response.txt`. Argon conversation `600688d7-e052-4ef5-b6fd-d7b58a64ccb4`.

Round2: `postintegration_round2/review.txt`, `verdict.json`, executed verification,
`root_response.txt`, and **`closure.txt` / `closure_verdict.json`**. Argon conversation
`944f3f83-0316-46fc-8b8b-970e5aea569e`. The closure supersedes identified errors in the
original review narrative. Argon withdrew invented catalogue numbers/definitions,
incorrect figure and XRF descriptions, causal overstatement and overstated submission
readiness. Original reports remain preserved. Final verdict: **ACCEPT_WITH_BOUNDS**
for H1–H4, with no remaining material issue identified for those bounded additions.
Root independently agrees with this limited judgement.

Refinements added exact H4 transitions/control counts, CPU row-exclusion and GPU
sensitivity results, H3 pre-/post-quota overlap accounting, H1 comparator-pool and tree
details, and H2 locus ranking, opportunity definition, weights and permutation seeds.
The paper's scientific values were checked against source tables rather than model
narratives. Flash's erroneous registry narrative is retained as superseded.

## Verification and retained evidence

- H3 raw audit **IBEX 53796328 PASS**: matched FASTAs/index hashes, identifiers,
  raw alignments and all summaries/length strata independently checked.
- H4 GPU raw audit **IBEX 53796103 PASS**: 6,081,552 rows; hashes, IDs, lengths,
  coordinates, CIGAR aggregate closure and qTM outcomes. Searches ran on ORIX
  freecycle H200 array315233; raw GPU outputs are retained on IBEX.
- Root CPU census: **118,857 assertions PASS**. Original strict alignment-TM
  verifier failures remain documented; acceptance uses the audited qTM endpoint.
- Current and full local drivers **PASS**. The global registry checks **62 COMPUTE
  and 216 FROZEN claims**, including the **152 H1–H4 claims** backed by **27 pinned
  inputs**. These overlapping totals must not be added together.
- All 25 remote provenance paths and two local root-audit source paths exist.
  A falsely attributed, unused preclosure receipt was removed from the active
  manifest and replaced by the executed root census audit.
- PDF build: resolved references, 35 pages. Integrity scan flags only legitimate
  white plot labels 69, 90 and 61 on page32. Final diff passes whitespace checks.
- Retained literature includes full text for SMAG, its correction and Pyrodigal;
  XML serves as the full-text archive where publisher PDF retrieval failed.

The IBEX review mirror is
`/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p7_analyses_clean/manuscript_review_20261014/`.
Its scope and file inventory are in `manuscript/IBEX_MIRROR_SCOPE.md` and
`manuscript/ibex_mirror.sha256`; the Slurm verification receipt is recorded separately
in `manuscript/ibex_mirror_verification.json`. Heavy analysis data remain at the
production paths mapped in `reproducibility/darkmatter/round7_analyses/README.md`.

## Remaining whole-paper prerequisites

This is a supported revision, not a complete raw-to-final or submission-ready
certification. Open items include 78 missing raw-file paths, 302 permission failures,
unverified duplicate-run merge commands, incomplete historical claim migration,
public/reviewer access to ENA PRJEB125153, MAG deposition with MIMAG metadata and an
archived release DOI. The laboratory XRF numerical table is retained and the derived
environmental table rebuilds exactly; units, instrument/detection-limit metadata,
original climate API responses/model selection and relevant confounding remain open.
Final author/contribution information and journal submission packaging also require
completion. No submission or public push was performed.
