# RQ-A5 — Result

## ANSWER (short)
The decoded dark fraction is **enzymatically conventional, not chemically exotic.** When
the dark families that survive a top-quartile decode-confidence cut (C_confidence ≥ 0.369,
n = 5,863; 5,418 of them resolving to a fully-specified 4-level EC) are projected onto EC
space, they fall squarely into the seven canonical enzyme classes, dominated by
oxidoreductases and transferases, and they map overwhelmingly onto **well-characterised
central-metabolism pathways** (fatty-acid biosynthesis, TCA, amino-acid metabolism,
tRNA charging). We find **no signal that confidently-decoded dark families encode
fundamentally novel reaction chemistry.** The novelty of these families is at the
**sequence/fold level** (they have no UniRef50 ortholog) while the *predicted reaction*
is a known one — i.e. **structurally novel proteins performing canonical chemistry**, not
new EC numbers. This is the honest, defensible reading; a true "absent-from-all-genomes"
or pathway-hole-completion claim cannot be made from this table alone (see caveats).

## Key numbers
- Total decoded dark families: **23,382**; with any `best_ec`: 23,083.
- Confidence: min 0.045 / median 0.276 / Q3 **0.369** / max 0.806. Cutoff = top quartile.
- Confident decodes (C ≥ 0.369): **5,863 (25.1%)**; of those with a fully-specified EC:
  **5,418 (23.2% of all families)** — the analysis set.
- core_secondary among the analysis set: **secondary 5,195 / core 223** (the decode field
  is dominated by the "secondary" tier across the whole catalogue, 22,886 / 23,382).

### EC-class landscape (families counted once per class; % of 5,418 confident-specific)
| EC class | families | % |
|---|---|---|
| 1.x Oxidoreductases | 1,707 | 31.5 |
| 2.x Transferases | 1,664 | 30.7 |
| 3.x Hydrolases | 1,340 | 24.7 |
| 4.x Lyases | 944 | 17.4 |
| 6.x Ligases | 425 | 7.8 |
| 5.x Isomerases | 421 | 7.8 |
| 7.x Translocases | 9 | 0.2 |

(Sum > 100% because a family may carry alternative ECs in different classes.)

### Top decoded pathways (446 distinct pathways hit)
PWY-5994 (fatty-acid biosynthesis, 398) · TCA / TCA-1 (255 / 142) · PWY-8173/8174/8175
(172–173 each) · VALDEG-PWY valine degradation (168) · PWY0-1337 (157) · PWY-8279 (138) ·
PWY-5913 (137) — i.e. the busiest decoded routes are core carbon/amino-acid metabolism.

## Candidate novel-chemistry shortlist
Full ranked list: `candidate_novel_chemistry.tsv` (3,833 confident + named-reaction +
pathway-mapped families). The point of the shortlist is *pathway-completion candidacy* and
high-confidence named chemistry on otherwise-dark sequences — not a claim of new EC numbers.
Top 10:

| family_id | EC | reaction | pathway | C |
|---|---|---|---|---|
| V38PRr3…463465_4 | 4.1.99.17 | phosphomethylpyrimidine synthase | PWY-6890 | 0.806 |
| V30PRr1…825193_17 | 3.1.2.14 | palmitoyl-[acp] thioesterase | PWY-5994 | 0.776 |
| 28Dr1…1292762_3 | 2.3.1.8 | phosphate acetyltransferase | PWY-1281 | 0.770 |
| F41PRr1…611041_35 | 2.6.1.1 | aspartate aminotransferase | TCA-1 | 0.750 |
| F32PRr1…605508_222 | 2.6.1.1 | aspartate aminotransferase | TCA-1 | 0.743 |
| 62PRr2…504022_43 | 2.5.1.54 | DAHP synthase | PWY-6164 | 0.740 |
| F1PRr3…1280418_2 | 6.3.4.18 | N5-CAIR synthase | PWY-6123 | 0.740 |
| V24Dr1…1034307_7 | 3.1.2.14 | palmitoyl-[acp] thioesterase | PWY-5989;PWY-5994 | 0.730 |
| F15Dr3…21868_29 | 2.1.1.14 | cobalamin-indep. methionine synthase | PWY-5041 | 0.715 |
| F1Dr2…36471_5 | 2.3.1.30 | serine acetyltransferase | PWY-7870 | 0.700 |

Reading: these are **dark (no canonical ortholog) sequences confidently decoded to known,
named reactions on annotated pathways** — the most defensible "interesting" cases are dark
proteins that could be the *unrecognised* catalysts for these steps in EQ genomes
(candidate gap-fillers), e.g. cobalamin-independent methionine synthase (B12-independent
metabolism — ecologically relevant in cobalamin-poor desert soil) and palmitoyl-ACP
thioesterase (fatty-acid biosynthesis termination).

## Honest caveats
- **No new EC numbers.** Every confident decode lands in a known EC; "novel chemistry" in
  the strict reaction sense is **not supported** by this evidence. The robust claim is
  *sequence/fold novelty carrying canonical function*.
- **Cannot prove genome-wide absence.** No local per-MAG KO/EC matrix is available here, so
  "reaction absent from all genomes" and definitive broken-pathway completion are **not
  testable locally**. The shortlist gives *candidates*, ranked by decode confidence, for
  follow-up against MAG KO coverage.
- **Decode confidence is a model score.** C_confidence orders candidates; it is not
  experimental validation. The cutoff (top quartile) is defensible but arbitrary at the
  margin; lowering it adds families but not new EC classes.
- **`core_secondary` is a decode tier, not a metabolism class** — the catalogue is ~98%
  "secondary", so this field is reported for completeness rather than as a biological
  secondary-metabolism enrichment signal.
