# pathway_check_clean — genome-level verification layer (§2.3/§2.4 censuses vs gapsmith)

Owner: velazqal. Created 2026-08-15 (session 03-bgc-resistome).

## STATUS (2026-08-15) — DONE, armor role only
Three independent methods call the same capacities on the current 5,438 reps:

| comparison | carriers | agreement | Jaccard |
|---|---|---|---|
| ectoine: KO census vs antiSMASH | 895 vs 899 (860 shared) | 98.6% | **0.921** |
| trehalose: KO census vs gapsmith draft | 4,099 vs 4,075 (3,620 shared) | 82.8% | **0.795** |
| ectoine: KO census vs gapsmith binary | 895 vs 15 | 83.8% | 0.017 |

- The ectoine KO×antiSMASH concordance is the §2.3/§2.4 verification armor: protein-KO
  panel and nucleotide rule-based detection are fully independent pipelines.
- gapsmith's **binary** ectoine call (P101-PWY `Prediction`) is insensitive, not
  contradictory: it is gated on the key reaction, which the draft org-db rarely maps for
  ectoine synthase. The **continuous** completeness tracks the KO calls: median 60 in
  KO+ carriers vs 40 in KO−; 647/895 carriers ≥60%.
- Basis gates verified: 5,438 gapsmith model dirs = the current dRep reps (2026-08-13);
  `all-Pathways.tbl` is written by `gapsmith find` **before** gap-filling (log step
  [1/5] vs [4/5]; checked 2026-08-15), so completeness is annotation-derived and safe
  for capacity claims. Gap-filled models are used for NO capacity claim;
  `filled-added.tsv` is the audit trail.
- Role: verification armor for the §2.3 osmolyte census and §2.4 ectoine BGC number.
  NOT a new Results subsection (PI decision on genome-level metabolism pending —
  coordination item 1).

## Pipeline
- scripts/check_1_triangulate.py → results/triangulation_per_rep.tsv,
  results/triangulation_summary.tsv
- Inputs: osmolyte_clean/results/rep_traits.tsv (KO calls), bgc_clean/results/
  bgc_per_mag.tsv (antiSMASH ectoine), 06_darkmatter/gapsmith_models/<rep>/
  <rep>-all-Pathways.tbl (P101-PWY; trehalose biosynthesis I–VII MetaCyc IDs).
