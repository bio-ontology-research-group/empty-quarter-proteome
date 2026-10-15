# Ladder rewrite — splitting Table 3 into two catalogue-specific tables (2026-07-29)

Scope: `paper/main.tex` (abstract, B1, B1a, B2, B3, B5, C3, Tables 3/4/6, four figure captions),
`reproducibility/verify.py` (CLAIMS registry) and its generated `reproducibility/claims.tsv`.
Nothing else was touched. Nothing was committed.

Driver: A5 in `workbench/reviews/2026-07-29-opus-number-verification.md` §2 and Check 4 in
`reviews/2026-07-29-number-verification.md`. Both audits agree that every individual number in the
old ladder is correct and that the *chain* is not: it presented three clusterings across two protein
catalogues as one nested descending sequence.

Shorthand: `$M = /ibex/scratch/projects/c2014/rob/metagenomics`,
`$W = /home/leechuck/Public/software/empty-quarter-workbench` (= `workbench/` here).

---

## Before

One table (`tab:ladder`, old Table 3), 12 rows, one implied chain, no denominators:

```
Proteins (total)                                    ~51.9M
UniRef50 no-hit reps (sequence-dark)                 2.84M
  initial predicted-3Di screen: no AFDB match        2.53M    89.2%
Novel multi-member families (30% clustering)       334,810
Candidate novel-fold families (>=2 members, ...)    68,808
  explicit models with pLDDT score                  63,654
  dual-novel, pLDDT >= 0.00 ... 0.80              15,248 ... 15
```

Prose (B1, old `:164-166`) read it as a chain:
`51.9M -> 2.84M -> 2.53M -> 334,810 -> 68,808`.

Three defects: (i) the 23,416,070-representative denominator never appeared, so 2.84M looked like
5.5% of 51.9M rather than 12.1% of the representatives; (ii) the arrow into 334,810 reversed the
pipeline (those families come from the full 2,839,030 sequence-dark set, not from the 2,533,446
structural no-matches); (iii) 68,808 was presented as a subset of 334,810 although it comes from the
tier1 cohort of the *753,305,533-protein assembly catalogue* — a different catalogue entirely.
The abstract additionally attributed selection/structuring/endemism (properties of a *third* object,
the 23,382 dark orthogroups) to the 334,810 families.

## After

Two tables, no arrows, every row carrying its own aggregation level and denominator, and no relation
asserted across the catalogue boundary.

**Table 3 `tab:ladder` — "Dark matter in the MAG-Bakta catalogue"** (Chain A + Object B)

```
Bakta protein calls, 15,469 MAGs            protein          51,920,778   ---
mmseqs2 linclust 90% id / 80% cov           representative   23,416,070   of 51,920,778 proteins
 [Sequence screen vs UniRef50 — denominator 23,416,070 representatives]
   UniRef50 hit                             representative   20,577,040   of 23,416,070
   UniRef50 no-hit (sequence-dark)          representative    2,839,030   12.1%
 [Predicted-3Di screen vs AFDB — denominator 2,839,030 sequence-dark representatives]
   no AFDB match (ProstT5+FoldSeek)         representative    2,533,446   89.2%
   AFDB hit                                 representative      305,584   10.8%
 [Family clustering at 30% id of the SAME 2,839,030 sequence-dark representatives]
   clusters formed                          family            1,875,540   of 2,839,030
   multi-member families (>=2)              family              334,810   of 1,875,540 clusters
     representatives in them                representative    1,298,300   of 2,839,030
   singletons                               family            1,540,730   of 1,875,540 clusters
```

**Table 4 `tab:tiers` — "Structural novelty in the assembly catalogue"** (tier1/tier2 + Object C +
the B1a ESMFold sensitivity rows, which belong to this catalogue and were moved out of Table 3)

```
Assembly protein catalogue                  protein         753,305,533   ---
  clusters at 90% id                        cluster         449,265,306   of 753,305,533
  UniRef50-dark clusters                    cluster          68,646,640   15.3% of 449,265,306
 [Tier1 — MAG-associated dark clusters]
   tier1 input                              representative   11,418,079   of 68.8M dark reps
   AFDB structural hit                      representative      883,688   7.74%
   no AFDB match (raw; fragment-inflated)   representative   10,534,391   92.26%
   no AFDB match, >=100 aa (length-gated)   representative    2,339,154   ~74% at >=100 aa
 [Tier2 — non-MAG, >=100 aa; for contrast, not used elsewhere]
   tier2 (non-MAG)                          representative   57,379,961   of 68.8M dark reps
   tier2 structural input (>=100 aa)        representative   13,088,646   of 57,379,961
   no AFDB match                            representative   10,113,496   77.27%
   novel-fold families                      family              156,743   of 11,705,404 tier2 clusters
 [Structural clustering of tier1 (single-step FoldSeek)]
   clusters formed                          family           10,898,521   of 11,418,079
   multi-member (>=2)                       family              220,177   of 10,898,521
   candidate novel-fold families            family               68,808   of 220,177
     proteins in them                       protein             169,651   of 11,418,079
 [Explicit ESMFold validation of the 68,808 family reps (B1a) — unchanged values]
   models built / scored                    model         68,527 / 63,654
   dual-novel pLDDT >=0.00 ... >=0.80       model         15,248 ... 15   24.0% ... 8.0%
```

The Table 4 caption carries both tier1 caveats: the fragment inflation (median no-hit 56 aa; 77.8%
below 100 aa) and the fact that the cohort is the *pre-correction* partition (re-tiering against
DAS_Tool-refined bins demoted 9,515,846 reps to tier2, leaving 1,902,233).

---

## What changed, item by item

1. **Abstract** (`:43-51`). Now states the 23,416,070 denominator, calls the 2,839,030 and 2,533,446
   objects *representatives* rather than "proteins", says the 334,810 come from re-clustering the
   same 2,839,030 sequence-dark representatives at 30% identity, and re-attaches purifying
   selection / environmental structuring / 21.4% endemism to the **23,382 dark orthogroups**
   (50% identity over the 5,438 dereplicated MAGs) instead of to the 334,810 families.
2. **B1** rewritten into two labelled blocks (MAG-Bakta catalogue; assembly catalogue) with no
   arrows. States explicitly that the 334,810 families are parallel to, not downstream of, the
   predicted-3Di screen, and that the 68,808 families and the B1a work do not descend from the
   MAG-catalogue ladder. Carries the fragment-inflation caveat and the superseded-partition caveat.
3. **B1a** first sentence now identifies the 68,527/68,808 cohort as tier1 assembly-catalogue family
   representatives carrying tier1's fragment inflation. All B1a numbers unchanged.
4. **B2** now opens by defining the 23,382-orthogroup object and stating that B2/B3/B5/C3 all use
   23,382 as their denominator; adds the evaluable-alignment counts n=22,767 (dark) and n=21,579
   (control), which were previously unstated.
5. **B3 / B3d** attributed to the 23,382 orthogroups; occupancy now printed as 10,669/23,382 (45.6%)
   and 543 (2.3%); the abundance statistics labelled as the 23,382 x 551 matrix.
6. **B5** decode attributed to the same 23,382 orthogroups (the "10 well-supported EC" figure itself
   untouched — see "not verified" below).
7. **C3** attributed to the 23,382 orthogroups; 21.4% now printed as 5,011/23,382. Table 6
   (`tab:specificity`) caption now says "orthogroups (50% id over 5,438 dereplicated MAGs — not the
   334,810 families of Table 3)".
8. **Figure captions** for B2, B3 (x2) and C3 now name the 23,382 dark orthogroups.
9. **`reproducibility/verify.py` CLAIMS + regenerated `claims.tsv`**. The six rows that cited a
   repository or a run rather than a file now cite files:

   | claim_id | before | after |
   |---|---|---|
   | `scale.proteins` | `workbench (Bakta; IBEX_SOIL_PROGRESS)` | `workbench/CHANGELOG.md:2039` |
   | `b1.dark` | `workbench/results/darkmatter (mmseqs/UniRef50)` | `$M/04_structure/seq_vs_uniref50/nohit_ids.txt` |
   | `b1.structnovel` | `workbench (ProstT5+FoldSeek vs AFDB)` | `$M/04_structure/foldseek_results/struct_nohit_ids.txt` |
   | `b1.families` | `workbench novel-fold clustering` | `$M/04_structure/novelty_analysis/novel_family_summary.tsv` |
   | `b1.folds` | `workbench is_novel_family` | `$M/04_structure/tier1_cluster/novel_fold_analysis/family_summary.tsv` |
   | `b1a.structures` | `ORIX ESMFold run` | `$M/06_esmfold_novelfold/pdb/` |

   (paths written in full in the ledger; `$M` abbreviated here only). 19 rows were added for the
   numbers the two new tables now print — every one of them citing a file — and the B2/B3/B5/C3 row
   descriptions were amended to name the 23,382-orthogroup denominator. Ledger grew 82 -> 102 rows.

## Verification performed

Every count written into the manuscript was checked against its source before it was written:

| value | how checked |
|---|---|
| 23,416,070 / 2,839,030 / 20,577,040 | `wc -l` on `$M/04_structure/seq_vs_uniref50/{all_rep_ids,nohit_ids,hit_ids}.txt` (IBEX); partition closes |
| 2,533,446 / 305,584 | `wc -l` on `struct_nohit_ids.txt` and `merged_afdb_hit_ids.txt`; partition closes to 2,839,030 |
| 1,875,540 / 334,810 / 1,540,730 / 1,298,300 | `awk` over `$M/04_structure/novelty_analysis/novel_family_summary.tsv`; 1,540,730 + 1,298,300 = 2,839,030 |
| 68,808 / 169,651 | `awk '$8==1'` over `$M/04_structure/tier1_cluster/novel_fold_analysis/family_summary.tsv`; matches `wc -l $M/06_esmfold_novelfold/novelfold_reps.txt` = 68,808 and `$W/results/tier2_structural_novelty/headline_tier1.txt` |
| 10,898,521 / 220,177 | same file + `headline_tier1.txt` |
| 753,305,533 / 449,265,306 / 68,646,640 | `wc -l assembly_proteins.mmseqsdb.index`; `$M/04_protein_catalog/catalog_ladder/catalog_ladder_summary.tsv` |
| 11,418,079 / 883,688 / 10,534,391 / 92.26% / 56 aa / 77.8% / 2,339,154 / ~74% | `$W/results/tier1_structural_novelty/SUMMARY.md` |
| 57,379,961 / 9,515,846 / 1,902,233 | `$W/CHANGELOG.md:1306-1307` |
| 13,088,646 / 10,113,496 / 77.27% / 156,743 / 11,705,404 | `$W/results/tier2_structural_novelty/{README.md,headline.txt}` |
| 51,920,778 | `$W/CHANGELOG.md:2039` (job 46528255) |
| 23,382 | `grep -c '^>' $M/06_darkmatter/darkfam_ecology/endemism/query.faa` (IBEX) = 23,382; definition from `docs/dark-matter-pipeline.md:25` |
| 22,767 / 21,579 | recount of the `set` column of `results/darkmatter_ecology/dnds.tsv` |
| 10,669 / 543 / 45.6% / 2.3% | `results/darkmatter_ecology/ENDEMISM.md` |
| 5,011 / 23,382 = 21.4% | `results/global_specificity/README.md:52` |
| 5,438 dRep MAGs, 50% id | `docs/dark-matter-pipeline.md:15-25` |

No number was invented, altered, or re-rounded. Two derived percentages that would have been new
(the UniRef50-hit share of 23,416,070, and the representative-to-protein reduction) were **not**
printed; those cells give the denominator as a count instead.

## Could not verify / disagreements recorded rather than fixed

1. **The tier denominators do not close to a single figure.** tier1 11,418,079 + tier2 57,379,961 =
   68,798,040 ("68.8M" in `$W/CHANGELOG.md:1302`), but `catalog_ladder_summary.tsv` gives
   `clusters_uniref50_dark = 68,646,640` — a 151,400 difference between two products of the same
   catalogue. I did not reconcile or adjust either number. The tables therefore say "of 68.8M dark
   reps" for the tier rows and quote 68,646,640 only in its own row, with its own denominator
   (449,265,306). This needs an owner decision, not an edit.
2. **51,920,778 rests on the CHANGELOG entry alone.** Unlike every other Chain A count there is no
   id-list file to `wc -l`; job 46528255's record is the only artefact. The downstream 23,416,070 in
   the same line *is* independently confirmed by `all_rep_ids.txt`.
3. **Tier1's fragment statistics and tier2's numbers were taken from the workbench SUMMARY/README
   files**, not recomputed from the FASTAs. They are internally consistent
   (883,688 + 10,534,391 = 11,418,079; 13,088,646 − 2,975,150 = 10,113,496) but I did not re-derive
   the length distribution.
4. **`b5.enzymes` = 10 and `P(correct) ≈ 0.34–0.48` remain unverifiable** (audit B-3:
   `family_ec.tsv` has 23,382 rows and no "well-supported" column). I re-attached the object to the
   23,382 orthogroups but left the 10 untouched.
5. **`claims.tsv:7` (`scale.samples`, "547 / 60x5") still cites "sampling design metadata"**, the
   one remaining non-file source. It is not a ladder row, and the audits' D-2 shows the sample count
   is itself unsettled (547 / 550 / 551 / 573 across products), so choosing a file here would have
   meant choosing a number. Left for the D-2 fix.
6. `claims.tsv:5` (`scale.nphyla`) cites a directory rather than a file. Out of scope here; noted.
7. The 4,592-model gap between "68,527 built" and the 63,935 pLDDT records (audit §3) is untouched;
   the table prints 68,527 and 63,654 with their denominators, as before.

## Checks

- `python3 reproducibility/verify.py` → **COMPUTE: 46 pass, 0 fail | TABLE/PROVENANCE documented: 58**,
  exit code 0. No MISS/FAIL/ERROR rows.
- `pdflatex -interaction=nonstopmode -draftmode main.tex` (run twice for references) → exit 0,
  **zero lines beginning with `!`**, no undefined references. Tables render as Table 3 and Table 4;
  downstream tables renumber to 5/6/7. Three `Overfull \hbox` warnings remain, all pre-existing and
  outside this edit (B4 prose line 249; Methods line 323). Both new tables fit the text block.
