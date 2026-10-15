# Review track — items to resolve

Reviewing and cleaning the manuscript, alongside `PLAN.md`'s section roadmap. One conversation each.
Numbers below were measured on disk 2026-08-18.

## R1 · Citations
13 open claim-support decisions in `resources/literature/CITATION_AUDIT.md` §4c, plus D-16 in §9.
Sharpest is D-1: §3.1's "MGnify90 is dominated by aquatic and engineered biomes" is contradicted by the
paper it cites.

## R2 · Scaffold vs Results — ✅ done 2026-08-18
Compared claim by claim. A0-A5, B1, B2-B5, C1-C3 fenced (`\iffalse`); B1a left open for the §1.3 session;
§1.1 gained the six-single-gaps / two-candidate-pairs detail; Box 1's retracted negatives trimmed.

## R3 · Numbers in prose vs the registry
**Partly done 2026-08-19:** the abstract no longer carries 573, 12.1%, 51.9M or the one-bin gap (see CORRECTIONS). Still open: `tab:scale`'s 573/551/550 and 550/527/126/304/184 rows, and L1241.
Every digit against `notes/CATALOGUE_NUMBERS.md`, both directions. Known: the abstract says **573
paired-read metagenomes** against the 551 headline, and still leads with the historical R226 block
(51.9M → 23,416,070 → 2,839,030 / 12.1%).

## R4 · Captions and `\ref` integrity
Six floats are now **uncited in live text** — their only references were in the scaffold R2 fenced:
`fig:concept`, `fig:dnds`, `fig:lifeonair`, `fig:picrust`, `tab:decode`, `tab:picrust`.
`tab:atmo`'s only live mention sits inside a `\todo`.

## R5 · Supplementary alignment
S1-S7 numbering and contents in `paper/supplementary.tex` against the plain-text "Supplementary Table~S*"
strings in `main.tex` (L154, L182, L198, L208).

## R6 · Methods — ✅ done 2026-08-19
**7,028 → 5,001 words (−29%)**, all fourteen subsections. Numbering unchanged.

| sub | before | after | sub | before | after | sub | before | after |
|---|---|---|---|---|---|---|---|---|
| M.1 | 716 | **475** | M.6 | 627 | **469** | M.11 | 405 | **317** |
| M.2 | 441 | **242** | M.7 | 872 | **664** | M.12 | 587 | **509** |
| M.3 | 503 | **147** | M.8 | 405 | **321** | M.13 | 448 | **235** |
| M.4 | 638 | **294** | M.9 | 314 | **297** | M.14 | 185 | **224** |
| M.5 | 418 | **361** | M.10 | 469 | **446** | | | |

Rules applied are in memory as `feedback_methods_trim_rules`. Errors corrected (all pre-existing):
M.4's "9.89% of the catalogue" → of the 43,423,322 reps in an HQ/MQ bin · M.3's dRep "default weighting
except contamination ×5" (×5 **is** the default) · M.4's pyrodigal "translation table 11" (never passed).

**Structural:** M.5 now defines both recruitment runs, both denominators and both presence rules once ·
M.13 split into M.13 *Statistics* + M.14 *Data, code and reproducibility* · section title
"Methods (pointer-level)" → "Methods" · 573/552/569 dropped so 551 and 549 are the only sample counts.

**Results-side edits made here:** §2.3 gained the ectoine cross-validation clause (Jaccard 0.92), which
had been sitting unused in M.11 · three stale `M.7` pointers repointed to M.5 (Results L160,
`fig:specificity` caption, M.12) · M.1's transect length restored to 1,043 km, which the first trim had
dropped while Results still cited it.

**Dropped and recoverable only from the workbench:** M.11's trehalose Jaccard 0.80 and the explanation
that the reconstruction's binary ectoine call is insensitive rather than discordant (15 calls in 5,438;
completeness 60% vs 40%).

**Open for R3:** 573 orphaned in the abstract (L39), `tab:scale` and L1241; the one-bin dRep gap
(12,433 of 12,434) still in the abstract at L41.

## R7 · Flow and redundancy — Introduction done 2026-08-20
**Intro done (session R-intro):** 1,023 printing words -> 776; PI skeleton + delta block fenced, not deleted.
Still open: Results is 4,299 words against the 2-3k narrative target. Repeated gradient statements across
§2.3/§2.4/§2.5; bridges between subsections.

## R8 · Voice — Introduction done 2026-08-20
4 live `---` in the Introduction cleared (1 in the heading, 3 in the now-fenced bullets). Remaining em
dashes are in Results. Against the writing guides.

## Decisions still mine
- **D-16** — supplement citations to `\citep` (needs a bibliography attached) or leave as author-year prose.
- **`tab:picrust` + `fig:picrust`** — dead 15,469/436 basis, captioned superseded; cut or keep.
- **Tawfiq companion citation** — "in preparation" or drop.

## Stale text noticed
- `CORRECTIONS.md` open-block table: says 16 citation decisions (it is 14), and shows the structural search
  as RUNNING.
- `PAPER_STRUCTURE.md` §2.2 holds the current and the superseded numbers in one block.
