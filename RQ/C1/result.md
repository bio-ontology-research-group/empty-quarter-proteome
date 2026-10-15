# RQ-C1 — Does 16S/PICRUSt2 predict the gene content of these classified taxa?

**Status:** ✅ answered

## Answer
**Only partly.** Weighted NSTI is 0.16 and median KO presence-recall is 93%, yet 56.7% of MAG proteins lack a KO and signature-function abundance agreement is heterogeneous. The historical 2.40× richness ratio uses an asymmetric denominator; the shared relative-abundance rule gives 31.41×, so neither is retained as biological enrichment evidence.

**Evidence / reproduce:** `results/picrust2_bias/`; `code/picrust2_bias_analysis.py`; `reproducibility/verify.py` (`c1.*`). Fig `figures/fig_picrust2_bias.pdf`.

**References:** douglas2020picrust2; matchado2024limits; sun2020inference (keys in `manuscript/references.bib`).
