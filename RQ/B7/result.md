# RQ-B7 — Result

**Question.** Dispersal limitation vs environmental filtering — how much of dark-gene
biogeographic structure is pure-space vs pure-environment vs shared?

## Answer

**Both processes act, and environmental filtering dominates the explainable fraction.**
The dark proteome shows a significant positive **distance-decay** (Mantel community vs
geography **r = +0.161, p = 0.001**) — the signature of **dispersal limitation** — *and*
environmental gradients (climate, geochemistry) explain community variance **beyond** what
space alone captures, the signature of **environmental filtering**. In a clean Borcard
space-vs-environment variance partition, the **pure-environment fraction is larger than the
pure-space fraction**, and both are significant. The shared (spatially-structured
environment) fraction is small, so the two processes are largely separable here.

## Numbers

### Borcard variance partition (recomputed on the real 547-sample matrix)

| environment block | n | pure space (dispersal) | pure env (filtering) | shared | residual |
|---|---:|---:|---:|---:|---:|
| SPACE vs CLIMATE  | 456 | **0.025** (p=0.002) | **0.033** (p=0.002) | 0.014 | 0.928 |
| SPACE vs GEOCHEM  | 251 | **0.034** (p=0.002) | **0.052** (p=0.002) | 0.014 | 0.901 |

In both partitions pure environment > pure space; both unique fractions are significant
(p=0.002); the shared fraction is modest (≈0.014).

### Supporting per-block PERMANOVA (from `phase_b_summary.txt`)

| block | n | R² (alone) | R² (\|space) | p |
|---|---:|---:|---:|---:|
| space (lat/lon) | 547 | 0.038 | — | 0.001 |
| compartment | 547 | 0.020 | 0.020 (unique) | 0.001/0.002 |
| climate | 456 | 0.047 | **0.033** | 0.001 |
| geochem | 251 | 0.066 | **0.051** | 0.001 |
| vanadium (alone) | 251 | 0.004 | — | 0.325 (n.s.) |
| geochem \| space [core fams] | 547 | — | 0.056 | 0.001 |
| climate \| space [core fams] | 547 | — | 0.036 | 0.001 |

Climate and geochemistry retain most of their explanatory power after conditioning on space
(climate 0.047→0.033; geochem 0.066→0.051), confirming filtering operates **on top of**
dispersal limitation. Vanadium alone is not significant (p=0.325).

### Distance-decay

Mantel community Bray–Curtis vs geographic distance: **r = +0.161, p = 0.001** — significant
spatial turnover → dispersal limitation is real but weak-to-moderate.

## Interpretation

Dark-gene biogeography is jointly governed by the two canonical drivers of microbial
community assembly: **dispersal limitation** (spatial turnover; significant distance-decay
and a significant pure-space fraction) and **environmental filtering / species sorting**
(climate and geochemistry structure communities beyond space; the larger pure-environment
fraction). The dark proteome is therefore not spatially inert — it is ecologically sorted by
the abiotic Empty-Quarter gradient, while distance still imposes additional structure.

## Honest caveat

All R² values are **modest** (pure fractions 2.5–5.2 %; residual ≈ 90–93 %). This is **typical
of soil-microbiome biogeography**, where Bray–Curtis variance is dominated by sample-specific
and stochastic/unmeasured variation, and absolute PERMANOVA R² is routinely small even when
highly significant. The *ranking* (environment > space, both > 0, shared small) is robust and
significant; the absolute magnitudes should not be over-interpreted. The partition is computed
on environment-complete subsets (climate n=456, geochem n=251), and the geochem subset is
geographically uneven, so the geochem partition carries more uncertainty than the climate one.

## References (microbial biogeography processes)

- Borcard D., Legendre P., Drapeau P. (1992) *Partialling out the spatial component of
  ecological variation.* Ecology 73:1045–1055. (variance partitioning method)
- Martiny J.B.H. et al. (2006) *Microbial biogeography: putting microorganisms on the map.*
  Nat. Rev. Microbiol. 4:102–112. (dispersal limitation vs environmental filtering framework)
- Hanson C.A., Fuhrman J.A., Horner-Devine M.C., Martiny J.B.H. (2012) *Beyond biogeographic
  patterns: processes shaping the microbial landscape.* Nat. Rev. Microbiol. 10:497–506.
- Anderson M.J. (2001) *A new method for non-parametric multivariate analysis of variance.*
  Austral Ecology 26:32–46. (PERMANOVA / adonis)

## Files

- `run.py` — reproduces tables + figure (`python3 RQ/B7/run.py` from repo root).
- `b7_blocks.tsv` — per-block PERMANOVA table.
- `b7_varpart.tsv` — Borcard space/env/shared/residual partition.
- `fig_varpart.pdf` / `.png` — stacked-bar partition (pure space / shared / pure env).
- Source: `results/darkmatter_ecology/phase_b_summary.txt`, `code/phase_b_ecology.py`,
  `code/eco_stats.py`.
