# bgc_clean — antiSMASH BGC census on the current 5,438 catalogue (§2.4 recompute)

Owner: velazqal. Created 2026-08-12 (session 01-basis-and-launch).

## STATUS (2026-08-15, session 03-bgc-resistome) — DONE + SHIPPED (GitHub dbe367c / Overleaf 2fef715)
- Array **50413466** complete: 5,438/5,438 JSONs (antiSMASH 8.0.4 --minimal, prodigal),
  272/272 logs clean, no truncated files. Aggregation job **50533513**.
- Census (denominator 5,438 species): **97.48% carry ≥1 BGC** (5,301); 44,002 regions;
  **85.29% of regions on a contig edge** (measured; replaces the inherited "90%").
  terpene 3,379 (62.14%) / terpene-precursor 4,556 (83.78%, kept SEPARATE — decided
  2026-08-15) / T3PKS 41.49% / RiPP-like 36.52% / NRPS-like 31.28% / NRPS 25.06% /
  T1PKS 22.36% / **ectoine 899 (16.53%)** / NI-siderophore 14.90%. 81 product types.
- Community carriage (548 canonical samples, dual denominators as §2.2/§2.3):
  any_BGC median 48.2% whole-community / 99.4% of mapped.
- Ectoine cross-check: antiSMASH 899 vs §2.3 KO panel 895, 860 shared, **Jaccard
  0.921** → verification armor (see 06_metabolism/pathway_check_clean/).
- Historical census SUPERSEDED (terpene 64.5% / NRPS 25.3% / ectoine 17.0% of
  15,469 redundant bins; 58% of current reps were never scanned by it).
- Written: §2.4 para 1 + Methods M.10 + fig_bgc_product rebuilt (repo
  code/make_bgc_amr_figures.py; make_figures.py §2 superseded, harmonisation
  session pending). Registered: CORRECTIONS 2026-08-15, CATALOGUE_NUMBERS §2.4.
- SETTLED / do-not-re-litigate: terpene-precursor stays separate from terpene
  (category rollup would move the headline to 93% and break comparability);
  novelty % never computed (deferred list, not just un-headlined).
- Open: none for this folder.

## Why
The historical antiSMASH census (06_metabolism/antismash/, terpene 64.5% / NRPS 25.3% /
ectoine 17%) ran on the OLD 15,469 single-binning MAG set; 58% of the current reps
(`__self_` bins) were never scanned. CORRECTIONS row A4 = OLD → recompute on 5,438.
Census only; BGC novelty % is an upper bound (contig-fragmented) — do NOT headline.

## Inputs
- results/rep_index.tsv — 5,438 rows (mag \t fna), from 05_refine/_drep/out/dereplicated_genomes
- antiSMASH 8.0.4 (IBEX module), --minimal, --genefinding-tool prodigal (same design as
  historical 86_antismash_array.sh)

## Pipeline
- scripts/bgc_antismash_array.sh — sbatch --array=1-272%60 (BATCH=20, 8 CPU/24G/10h) →
  results/antismash/<mag>/<mag>.json (+gbk); idempotent, version logged per task
- scripts/bgc_1_census.py — JSON walk (areas[].products, old counting semantics;
  protoclusters[].category rollup; contig_edge from region features) →
  results/bgc_per_mag.tsv, bgc_product_summary.tsv, bgc_category_summary.tsv,
  bgc_census_summary.tsv. Taxonomy join = 05_refine/_gtdbtk_refined (tax_class naming).
- scripts/bgc_2_community.py — CoverM carriage, dual denominators →
  results/community_by_sample.tsv, community_summary.tsv

## Job log
- antismash array: **50413466** (1-272%60, submitted 2026-08-13) — COMPLETE.
- census aggregation: **50533513** (16 CPU, 2026-08-15) — COMPLETE, 0 areas-vs-gbk
  mismatches.
