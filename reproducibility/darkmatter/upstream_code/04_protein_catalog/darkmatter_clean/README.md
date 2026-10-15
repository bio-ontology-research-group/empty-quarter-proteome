# darkmatter_clean

Sequence-dark basis for the unified 90/80 protein catalogue. **Dark = no UniRef50 hit.**
Assembly namespace (`sample::sample_contig_gene`) only. Owner: velazqal.

## Inputs
- clusters: `../all_linclust90/all_cluster90.tsv` (449,265,306 reps / 753,305,533 proteins)
- rep seqs: `../all_linclust90/all_reps90.faa.gz`
- delta UniRef50 search: `../delta_vs_uniref50/{hit_ids,nohit_ids}.txt`
- UniRef50 DB: `../../04_structure/uniref50/uniref50`

## Steps (`scripts/`)
1. `darkclean_1_uncovered_ids.sh` → `01_uncovered_rep_ids.txt` — reps with no delta status (~37.8M)
2. `darkclean_2_extract.sh` → `02_uncovered_reps.faa.gz`
3. `darkclean_3_uniref.sh` → `03_uncovered_vs_uniref50/` — array search
4. `darkclean_4_combine.sh` → `04_catalog_dark/` — delta + uncovered → dark set + summary
5. `darkclean_5_mag.sh` → `05_mag_dark/` — dark ∩ MAG-mapping

## Results
| metric | value |
|---|---|
| catalogue clusters | 449,265,306 |
| UniRef50 hit | 379,326,678 (84.4%) |
| **UniRef50-dark (no hit)** | **69,938,628 (15.6%)** |
| still uncovered | 0 |
| MAG-mapped reps | 43,423,322 (of 449.3M) |
| MAG-mapped **dark** reps | **4,292,446 (9.89% of mapped)** |

## Job log
- darkclean_1 (uncovered ids) ✅ 37,828,539 uncovered (matches catalog_ladder) — job 49580861 (~20 min)
- darkclean_2 (extract seqs)  ✅ 37,828,539 — job 49582072
- darkclean_3a (chunk 120)     ✅ — job 49582600
- darkclean_3b (UniRef50 array)✅ 120/120 COMPLETED — job 49583033
- darkclean_3c (merge)         ✅ uncovered: 36,536,551 hit / 1,291,988 dark — job 49852641 (8:33)
- darkclean_4 (combine → dark set) ✅ 69,938,628 dark / 449,265,306 (15.6%); still_uncovered=0 — job 49856452 (28:41)
- darkclean_5 (MAG dark, member-level rewrite) ✅ 2026-08-03 job 49969014 (22 min, 39G)
    → genome-resolved: 88,247,222 proteins / 43,423,322 reps; dark reps 4,292,446 (9.89% of mapped reps);
      dark proteins 6,151,650 (6.97%). per_bin_dark.tsv.gz = 37,949 rows. ↔ old single-sample 51.9M/23.4M/2.84M.
