# cultured_comparison_clean: RAK-PR cultured genomes vs the 5,438 soil MAG species reps

Owner: velazqal. Promoted 2026-09-20 from the exploratory folder
`alelopezv/eqmags_poc/05_rakpr_overlap` (step one and step three).
This folder is authoritative; the exploratory copy is superseded.

## STATUS (2026-09-21): COMPLETE on all 722 catalogue genomes

Step three previously ran on the 651 genomes annotated under `E/results/annotation/bakta`. The
other 71 (soil batch, thermophilic) are annotated in a second tree,
`EmptyQuarter_Data/cultures/annotation/thermophilic/bakta`, which `step3_00_prepare.py` now reads
as a second Bakta root. All 722 catalogue genomes are in; `results/step3_missing_bakta.txt` is
empty and the script fails if it is not.

| step | script | how it ran | wall |
|---|---|---|---|
| inputs + path lists | `00_prepare_inputs.py` | login | seconds |
| skani sketch + search + dist | `01_skani.sbatch` | Slurm **52131417**, 32 CPUs, 128G | 00:00:35 |
| binning + shared species | `02_classify.py` | login | seconds |
| ANI histogram | `03_plot.py` | login | seconds |
| step-three FASTAs + chunks, 722 | `step3_00_prepare.py` | login | ~4 min |
| DIAMOND database, 722 | `step3_01_makedb.sbatch` | Slurm **52196464**, 8 CPUs, 32G | 00:00:07 |
| dark reps vs cultured (array 0-20) | `step3_02_dark_vs_cultured.sbatch` | Slurm **52196465**, 32 CPUs, 128G each | 00:02:03 to 00:04:08 per task |
| cultured vs UniRef50, chunks 00-11 | `step3_03_cultured_vs_uniref50.sbatch` | Slurm **52131486**, 16 CPUs, 256G each | 07:49 to 11:36 per chunk |
| cultured vs UniRef50, chunks 12-13 | `step3_03_cultured_vs_uniref50.sbatch --array=12-13` | Slurm **52196463**, same resources | 05:39 and 06:53 |
| summary tables | `step3_04_summarise.py` | login | ~4 min |
| length comparison | `step3_05_lengths.py` | login | ~2 min |

Chunk boundaries are now fixed per Bakta root rather than recomputed over the whole file: root 0
keeps its 12 chunks and root 1 gets chunks 12-13 at the same records per chunk. Chunks 00-11 were
verified byte-identical before and after the rebuild (md5, 12 of 12; the 20 dark chunks are
unchanged too), so the 12 UniRef50 searches of job 52131486 did not need rerunning.

**Nothing is carried over any more.** `results/uniref50_hits_job51696940/` (the 12 hit tables of
job 51696940, copied from the exploratory folder) was deleted; `step3_04_summarise.py` reads
`tmp/uniref_hits/` for all 14 chunks. The two searches were compared first, on chunks 00-11:
2,664,556 proteins hit in both, 30 hit only in the rerun, 20 only in job 51696940, so 25 of 651
genomes differ by 1 to 8 proteins and the dark total moves from 46,337 to 46,327 (1.7093 % to
1.7089 %). They are **not** bit-identical: the MMseqs2 prefilter is not reproducible to the hit
across runs. The 50 proteins involved are 0.0018 % of the 2.71 M searched, and job 51696940's
tables are still on disk at
`alelopezv/eqmags_poc/05_rakpr_overlap/tmp/uniref_hits/` if the comparison is ever needed again.

Python `/home/velazqal/miniconda3/envs/python-utils/bin/python` (scipy 1.11.3). skani 0.3.1
(GTDB-Tk 2.5.2 env), DIAMOND 2.0.15 (`eggnog` env), MMseqs2 18.8cc5c (`mmseqs2` env).

## STATUS (2026-09-22): step 3.6, anatomy of the 23 pairs. EXPLORATORY, NOT REGISTERED

`step3_06_pair_anatomy.py` / `.sbatch` (Slurm **52269563**, 32 CPUs, 128G, 00:03:04) takes the
23 same-species pairs apart protein by protein. MAG side on the catalogue key: the pyrodigal
proteins of `04_protein_catalog/assembly_proteins/<sample>` on the representative's contigs
(`contig_to_bin_refined`), with the dark flag from a fresh parallel scan of `all_cluster90.tsv`
against `dark_rep_ids.txt.gz`. That scan returns 87,381 proteins and 1,513 dark, per
representative identical to `v2_rep_dark.tsv`. Cultured side from `tmp/cultured_all.faa` with
the step3_04 dark and hypothetical id lists and the Bakta `.tsv` products. Overlap = DIAMOND
2.0.15 blastp both ways (`--very-sensitive -e 1e-5 -k 5`), reciprocal best hit at >= 50 %
identity and >= 50 % coverage of both; the sensitivity rule is 30 % identity / 70 % coverage.
Outputs in `results/pair_anatomy/`. **These numbers are not registered** and nothing in the
paper rests on them.

| pooled over the 23 pairs | value |
|---|---|
| MAG proteins, cultured proteins, reciprocal pairs | 87,381 / 85,949 / 69,207 (67,054 under the sensitivity rule) |
| annotated vs dark MAG proteins with a partner | 80.3 % vs 15.7 % |
| the same on complete calls of >= 100 aa only | 88.9 % (n = 65,450) vs 15.7 % (n = 166) |
| annotated vs dark cultured proteins with a partner | 81.2 % vs 23.6 % |
| MAG-unique remainder, n = 18,174: partial, contig-edge, dark | 45.4 %, 50.8 %, 7.0 % |
| culture-unique remainder, n = 16,742: hypothetical, dark | 23.6 %, 4.5 % |

Annotated proteins pair up more often than dark ones in 23 of 23 pairs on both sides. Per-pair
medians: MAG-unique 137 aa against 268 aa for shared proteins, culture-unique 210 aa; shared
pairs are 99.7 % identical, 97.7 % of them at >= 95 %. The culture-unique products are mobile
and accessory (lipoprotein, phage protein, transposase, IS3 transposase, MFS transporter).

## Inputs

| input | count |
|---|---|
| RAK-PR catalogue `E/catalogue.tsv` (`E = alelopezv/eqmgc`) | 722 genomes, 722/722 FASTAs located |
| batches | enrichment 520, isolate 131, soil (thermophilic) 71 |
| Bakta proteomes, root 0 `E/results/annotation/bakta` | 651 genomes (enrichment + isolate) |
| Bakta proteomes, root 1 `EmptyQuarter_Data/cultures/annotation/thermophilic/bakta` | 71 genomes; 11 further Bakta runs there are not catalogue genomes and are skipped |
| cultured proteins | 3,024,784, of them 288,520 hypotheticals |
| MAG species reps `W/05_refine/_drep/out/dereplicated_genomes` | 5,438, all with a GTDB r232 call |
| abundance `W/10_thermal/abundance/rep_relative_abundance.tsv` | 5,438 x 550 samples |
| MAG dark reps `W/06_darkmatter/darkfam_basis_clean/results/dark_reps.faa` | 4,292,446 |
| matched control `W/06_darkmatter/specificity_clean/results/query_all.faa` | 106,042 dark family reps + 106,042 matched annotated family reps |
| MAG rep dark fractions | `../lineage_darkness_clean/results/v2_rep_dark.tsv`, 5,438 rows |

## Headline numbers

**Step one** (`results/summary.tsv`, `results/shared_species.tsv`) is unchanged: it always ran on
all 722. 601 pairs at the 75 % screen; shared = ANI >= 95 and aligned fraction >= 30.

| set | n | shared | related | none |
|---|---|---|---|---|
| all cultured | 722 | 189 (26.2 %) | 114 | 419 |
| enrichment | 520 | 100 (19.2 %) | 83 | 337 |
| isolate | 131 | 32 (24.4 %) | 19 | 80 |
| soil batch, thermophilic | 71 | 57 (80.3 %) | 12 | 2 |
| MAG species reps, reverse | 5,438 | 25 (0.46 %) | 31 | 5,382 |

203 shared-species pairs, 189 distinct cultured genomes, 25 distinct MAG species. All 189 now have
a proteome (was 132 of 189). On best ANI the 189 resolve to 23 MAG representatives, the same 23 as
before: the 57 genomes added here fall on representatives that were already paired.

**Step three, MAG-dark reps with a cultured homolog** (`results/step3_dark_with_cultured_homolog.tsv`,
denominator 4,292,446), 722 vs the earlier 651:

| metric | 651 | 722 |
|---|---|---|
| reciprocal coverage >= 50 % | 3,988 (0.093 %) | **4,226** (0.098 %) |
| any hit at e-5 | 25,365 | **25,403** |
| best identity >= 90 % | 891 | **1,042** |
| homolog in a shared-species genome | 1,275 | **1,300** |
| of those at >= 90 % | 648 | **679** |
| only in an unshared genome | 2,713 | **2,926** |
| >= 100 aa, of 1,466,943 | 1,520 (0.104 %) | **1,552** (0.106 %) |
| control, annotated family reps | 49,291/106,042 = 46.5 % | **49,170/106,042 = 46.4 %** |
| control, dark family reps | 170/106,042 = 0.16 % | **187/106,042 = 0.18 %** |

The control gap is 263-fold (was 290-fold). The annotated control fell by 121 although the target
database only grew: DIAMOND keeps 500 hits per query, so for a few queries the new cultured
proteins displace lower-scoring hits that had passed the reciprocal 50 % rule.

**Step three, cultured dark rate** (`results/step3_cultured_dark_rate_summary.tsv`):

| set | 651 | 722 |
|---|---|---|
| all cultured genomes | 46,337/2,710,913 = 1.71 % | **48,741/3,024,784 = 1.61 %** |
| shared-species genomes | 132 genomes, 5,194/525,525 = 0.99 % | **189 genomes, 6,713/777,916 = 0.86 %** |
| no shared species | 519 genomes, 1.88 % | **533 genomes, 1.87 %** |
| Bakta hypotheticals | 45,804/269,871 = 17.0 % | **48,128/288,520 = 16.7 %** |
| distinct matched MAG reps | 23 reps, 1,513/87,381 = 1.73 %, median 1.29 % | **unchanged, 23 reps** |
| cultured dark proteins recovered as MAG-dark reps | 1,961/5,194 = 37.8 % | **2,746/6,713 = 40.9 %** |

By batch: soil 2,414/313,871 = 0.77 %, enrichment 34,346/2,125,475 = 1.62 %, isolate
11,981/585,438 = 2.05 %. The soil batch is the least dark of the three, which is what the
shared-species split predicts: 57 of its 71 genomes share a species with a MAG (0.60 % dark),
against 1.45 % for the 14 that do not.

**Paired by species** (`results/step3_paired_by_species.tsv`, n = 23 both times, the unit is the
MAG representative):

| statistic | 651 | 722 |
|---|---|---|
| cultured lower | 16 of 23 | **18 of 23** |
| median difference | -0.003258 (-0.33 pp) | **-0.003258 (-0.33 pp)** |
| Wilcoxon two-sided | p = 0.0415 | **p = 0.0214** |
| sign test | p = 0.0931 | **p = 0.0106** |
| >= 10 strains | n = 4, median +0.07 pp | **n = 5, median -0.19 pp** |
| single strain | n = 10, median -0.49 pp | **n = 9, median -0.62 pp** |

The pairs are the same 23 species, with more cultured strains aggregated into each: the deeply
cultured stratum gains one species (`V53PRr2__self_SemiBin_82`, 6 strains to 29) and its median
flips sign, and the single-strain stratum loses one (`V35PRr2__self_SemiBin_99`, now 2 strains). Both p-values and both strata are now written by `step3_04_summarise.py`; on the old table
the same code reproduces 0.0415, 0.0931, +0.07 pp and -0.49 pp exactly.

**Lengths** (`results/step3_length_summary.tsv`): MAG dark reps unchanged, n = 4,292,446, median
70 aa, 65.83 % under 100 aa, 42.41 % under 60 aa, 56.40 % incomplete gene calls. Cultured dark
proteins n = 48,741 (was 46,337), median 59 aa (was 60), 78.49 % under 100 aa (was 77.60 %),
50.38 % under 60 aa (was 49.28 %).

## Reconciliation against the exploratory README

The step-one files are still byte-identical to the exploratory copies: `cultured_meta.tsv`,
`mag_meta.tsv`, `skani_search.tsv`, `cultured_best_hit.tsv`, `mag_best_cultured_hit.tsv`,
`shared_species.tsv`, `summary.tsv` (and `skani_dist.tsv` once sorted, same 601 pairs). Every
step-three file now supersedes its exploratory counterpart, which was built on 651 genomes.

## Notes

- `tmp/mag_dark_len.tsv.gz` (MAG dark rep id, codon count, prodigal partial flag) is derived by
  `step3_05_lengths.py` from the prodigal headers of `dark_reps.faa` when it is absent.
- `LIN` in `step3_04_summarise.py` points at `../lineage_darkness_clean/results`, so the MAG side
  of every cultured-vs-MAG comparison comes from the promoted native-key scan.
- Adding a further annotation tree means appending it to `BAKTA_ROOTS` and running
  `sbatch --array=<new chunks>` for step3_03; finished chunks are skipped.

## Caveats (carried forward)

The cultured genomes come from saline pans, oases and campground sites 41, 50, 52, 59 and 60, so
this is a site-matched and cultivation-biased comparison; only the same-species contrast is
licensed and the all-722 rate of 1.61 % is not comparable to the catalogue-wide MAG rate. Gene
calling differs (Bakta on cultured genomes, pyrodigal on assemblies), and the units differ (MAG
darkness is a property of a 90 % linclust representative, cultured darkness is per protein). The
722 are long-read assemblies but not all closed: median 2 contigs, 49 % single contig (the 651
alone had median 1). 28 of the 189 shared genomes are MQ. Four shared species are human-associated
taxa (Escherichia coli in 316 of 550 samples, Acinetobacter baumannii, Enterococcus_B lactis, the
Bacillus_A anthracis group); contamination remains a hypothesis to test, not a finding. MR59-69
carries 9,840 proteins at 2.49 % contamination, high for an isolate.
