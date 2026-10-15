# C1 — Reference bias: PICRUSt2 (16S) vs MAG-encoded gene content

The methodological pillar (Subversion 2): GTDB *classifies* the Empty Quarter taxa, but is their gene content
predictable from 16S? We compare PICRUSt2 functional predictions (already computed for the amplicon companion)
against the gene content of genome-resolved MAGs from the **same samples**.

## Design
- **PICRUSt2 side** (amplicon paper, `empty-quarter/data/processed/functional/picrust2/merged/`): KO metagenome
  (`ko_metagenome_unstrat.tsv`, 10,544 KO × 1,270 samples), weighted NSTI, sample metadata. PICRUSt2 v2.4.1.
- **MAG side** (`mag/`, built on IBEX by `scripts/ibex/130_mag_ko_metagenome.py`): per-sample KO metagenome from
  EggNOG `KEGG_ko` over all 15,469 HQ+MQ single-sample MAGs — CoverM-relab-weighted (`mag_ko_relab_by_sample.tsv`,
  PICRUSt2-comparable), unweighted pangenome copies (`mag_ko_count_by_sample.tsv`), per-sample MAG presence
  (`mag_ko_presence_by_sample.tsv`), and per-MAG annotation stats (`mag_annotation_stats.tsv`).
- **Join**: MAG sample = `mag_id.split('__')[0]` == 16S `sample_id` (identical scheme). **436 matched samples**.
- Analysis: `code/picrust2_bias_analysis.py` → `per_sample_concordance.tsv`, `marquee_functions.tsv`,
  `dark_fraction.tsv`, `SUMMARY.md`. Dark-fraction denominator = CheckM2 `Total_Coding_Sequences`.

## Headline (see SUMMARY.md)
NSTI is typical-soil (median 0.16) — **placement is not the problem**. PICRUSt2 recovers the housekeeping core
(median 93% KO presence-recall). Under the non-zero rule PICRUSt2 has exactly 10,543 predicted KOs in every
matched sample (a constant, not a median); MAG richness has median 4,482. The rebuilt 440-sample sweep exposes
two incompatible denominator rules. Thresholding PICRUSt2 against fixed MAG any-presence gives 2.40 at zero and
1.11 at \(10^{-5}\); applying the same relative-abundance threshold to both sides gives 31.41 and 14.63.
Neither supports a threshold-independent biological enrichment or reversal claim. Across ten prevalent
signature functions, abundance agreement is heterogeneous (Spearman **−0.12 to
0.42**, median **0.09**) rather than uniformly absent. The supported hard limit is that **56.7% of all MAG
proteins carry no KO** (19.1% no eggNOG ortholog), a ceiling on what any KO-based 16S predictor can represent,
plus 3,376 encoded KOs PICRUSt2 never predicts. The methodological pillar reinforces the dark-matter pillar.

## Provenance and denominator rules
- MAG-KO metagenome: IBEX job 47804253 (`05_picrust_validation/`), 532 WGS samples, 15,469 MAGs, 12,410 KOs.
- PICRUSt2: amplicon-companion runs (Trips 1–4 job, Trip 5 job 45925397), v2.4.1, NSTI<2 cutoff.
- `threshold_sweep.tsv:ratio_picrust_over_mag_relab` applies the listed relative-abundance threshold to both
  PICRUSt2 and CoverM-weighted MAG profiles.
- `threshold_sweep.tsv:ratio_picrust_over_mag_presence` applies the threshold only to PICRUSt2 and divides by
  the fixed count of KOs encoded by any recovered MAG; it is an asymmetric diagnostic retained for provenance.
