# C1 — Reference bias: PICRUSt2 (16S) vs MAG-encoded gene content

Matched samples (both 16S-PICRUSt2 and genome-resolved MAGs): **436**. GTDB classifies these taxa, so this tests *gene-content* predictability, not placement.

## Phylogenetic distance is NOT the failure mode
Weighted NSTI median **0.159** (p95 0.322) — within the typical-soil range and well under PICRUSt2's NSTI<2 cutoff. ASVs place fine; the failure is downstream, in gene content.

## What PICRUSt2 gets right — and the remaining limits
PICRUSt2 *does* recover the housekeeping core: median **93.3%** presence-recall of MAG-encoded KOs.

**(1) KO-richness differences are threshold-sensitive, not evidence of 2.4-fold enrichment.**
The original presence comparison counted every non-zero value. PICRUSt2 has exactly **10,543 predicted KOs in
every matched sample** (a constant, not a median), whereas MAG-encoded richness has median **4,482**. In the
rebuilt 440-sample sweep, thresholding PICRUSt2 against fixed MAG any-presence gives ratios of **2.40** at zero
and **1.11** at \(10^{-5}\). Applying the same relative-abundance threshold to both sides gives **31.41** and
**14.63**, respectively. These incompatible denominator rules do not support a biological enrichment or
reversal claim.

**(2) Abundance agreement is function-dependent.** Across ten prevalent EQ
signature functions (atmospheric CO `coxLMS`, hydrogenase `hyaB`, osmolyte
`otsAB`/`ect`, and UV-repair `phr`), PICRUSt2-vs-MAG across-sample Spearman
correlations range from **−0.12 to 0.42** (median **0.09**). This is
heterogeneous agreement, not a uniform failure.

**(3) The KO-space ceiling — the dominant failure.**
PICRUSt2 operates only in KEGG-KO space. Across 15,469 MAGs, **56.7% of all MAG-encoded proteins have no KO** (median per-MAG 56.1%); 19.1% have no eggNOG ortholog at all. This is a hard upper bound: even a perfect KO predictor represents none of that fraction, and the 3,376 KOs the MAGs encode that PICRUSt2 never predicts in any sample are an additional reference-database blind spot.

## Bottom line — the methodological pillar reinforces the dark-matter pillar
An amplicon-only survey of this biome recovers most of the KO core but cannot represent the no-KO majority of
the proteome, including by construction the structural dark matter central to this paper. It also does not
uniformly recover the sample abundance patterns of the atmospheric/osmolyte/UV signature functions. The
KO-richness contrast must be reported with threshold sensitivity and is not retained as a biological enrichment.

_Outputs: per_sample_concordance.tsv, marquee_functions.tsv, dark_fraction.tsv,
rebuild/threshold_sweep.tsv._
