# RQ-B3 — Methods: gene–environment association map for dark protein families

## Question
Which specific geochemical/climatic variables select which dark protein families? We
synthesise the per-variable indicator scans into one unified gene–environment map and
quantify, per axis, how many dark families are recruited (positive) vs filtered (negative),
how robust the associations are to geography, and how interpretable the recruited families are.

## Data
- **Dark-family abundance matrix** `results/darkmatter_ecology/darkfam_abundance.tsv`
  (23,382 multi-member dark protein families × 551 samples; abundance = read-recruitment
  signal per family per sample).
- **Environment matrix** `results/darkmatter_ecology/env_matrix.tsv` (547 samples; site,
  compartment, trip, lat, lon, temp, precip, aridity_index, salinity_proxy + 40 XRF
  elements; `has_xrf=1` for the 251 samples with XRF geochemistry).
- **Annotation joins**: `family_ec.tsv` (C·V·T decoded best EC + reaction/pathway, calibrated
  confidence 0.22–0.43), `dnds.tsv` (median dN/dS), `darkfam_meta.tsv` (membership/endemism).

## Axes tested (15)
XRF elements **V, Mo, Ni, Fe, Cu, Co, Mn, Zn, S, Cl, Na, Ca, K** and chemistry proxies
**salinity_proxy** (both gated on `has_xrf=1`, n=251) and **aridity_index** (climate, ungated,
n=456). Plant-host / compartment contrasts are handled in the sibling indicator tables and are
out of scope for this geochemical map.

## Association statistic (per `code/phase_e_indicators.py`)
For each axis we take the subset of samples carrying the variable, keep only families present
(abundance > 0) in **≥ 25** subset samples (`MIN_PRESENT`), and compute a vectorised
**Spearman rank correlation** of each family's abundance against the axis value. Two-sided
p-values from the t-approximation `t = r·√((n−2)/(1−r²))`. Per axis we apply
**Benjamini–Hochberg FDR** across all tested families and report counts at **q < 0.05** and
**q < 0.10**, split by sign of r (positive = recruited with the stressor, negative = filtered).

`run.py` re-implements this from the raw matrix rather than reading
`indicators/<VAR>_top.tsv`, because those files store only the **top 40 hits** per axis and so
cannot yield the full positive/negative split or the q<0.10 tier. The recomputation reproduces
the documented `INDICATORS.md` q<0.05 counts exactly (V=0, Mo=177, Fe=1644, S=3964, Cl=3594,
Ca=4066, Na=705, K=955, salinity=3179, aridity=491), which validates the pipeline.

## Space-partial control (per `code/phase_e_partial_osmolyte.py`)
For the five desert master axes **salinity_proxy, S, Cl, Ca, aridity_index** we additionally
compute a **partial Spearman** controlling for ranked latitude and longitude: rank-transform
abundance, axis and lat/lon, residualise abundance and axis on the design `[1, rank(lat),
rank(lon)]` via the hat matrix, and correlate the residuals (df = n−4). This separates genuine
environmental chemistry from spatial autocorrelation. We report the marginal q<0.05 count, the
partial q<0.05 count, and the **retained %** = fraction of marginally-significant families that
remain significant under the partial test (intersection-based; matches `PARTIAL.md`).

## Interpretability score
For each axis's significant set we report `frac_with_ec` (fraction carrying any decoded best
EC) and `n_osmolyte_ec` (families whose EC is in a curated osmolyte set: glucosylglycerol-P /
mannosyl-3-P-glycerate / trehalose-6-P phosphatases and related solute-synthesis ECs).

## Outputs
`gene_env_summary.tsv` (one row per axis), `fig_gene_env.pdf/.png` (a: signed #significant
families per axis; b: strongest dark family per key axis). Reproduce with
`python3 RQ/B3/run.py` from the repo root (numpy + scipy + matplotlib).

## Caveats
- EC labels are **low-confidence C·V·T decoder calls** (calibrated P 0.22–0.43): plausible,
  not certain. ~99% of families carry one, so "EC-less" is negligible — interpretability is
  limited by decoder *confidence*, not label presence.
- Abundance is relative read-recruitment; correlations are associative, not causal.
- The chemistry axes (Ca, Cl, S, salinity_proxy, Na) are strongly collinear, so they recruit
  overlapping family sets — per-axis counts are not independent.
- XRF axes use only 251 samples; aridity 456. Families below the 25-sample prevalence floor
  are untested.
