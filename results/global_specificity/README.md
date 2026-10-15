# C3 — EQ-specificity vs undersampling: dark families vs global catalogues + other deserts

Measures database-relative EQ specificity and its sensitivity to comparator choice by
searching the **same 23,382 UniRef50-dark family rep proteins** used for the MGnify90 test against additional
global catalogues and other deserts. **All searches use the identical hit definition** (DIAMOND blastp
`--sensitive --evalue 1e-3`, ≥1 hit = "detected in that catalogue"; no hit = "catalogue no-hit") so no-hit
fractions are directly comparable to the existing MGnify result. They do not establish geographic absence.

## Catalogues
| Catalogue | What | Source | Method |
|-----------|------|--------|--------|
| MGnify90 | 718M env protein clusters (all biomes) | EBI peptide DB | DIAMOND (done: 36.6% no-hit) |
| GEM | 52,515 metagenome MAGs (Nayfach 2021) | portal.nersc.gov/GEM `faa.tar` (26.7 GB, AA) | DIAMOND blastp |
| SPIRE | dereplicated representative genomes (Schmidt 2024) | black.embl.de rep genomes (86.8 GB, nt → prodigal) | DIAMOND blastp |
| Other deserts | Atacama, Mojave, arid, desert trace-gas autotrophy, Antarctic cold-desert | SPIRE per-study genecall FAA (predicted proteins) | DIAMOND blastp, per-desert + combined |

Desert studies (SPIRE per-study predicted proteins, all open):
- **atacama** — Uritskiy_2019_soil_Atacama
- **mojave** — JGI_soil_desert_Mojave
- **arid** — JGI_soil_arid
- **desert_autotroph** — Bay_2021_autotrophy_desert (trace-gas chemosynthesis desert; on-theme with the atmospheric pillar)
- **antarctic_colddesert** — JGI_soil_Antarctica

## Scripts / jobs
- `workbench/scripts/ibex/131_global_specificity.sh` — GEM + SPIRE (job 47810126, largemem).
- `workbench/scripts/ibex/131b_desert_specificity.sh` — 5 deserts + combined (job 47810171, batch).
- Query: `06_darkmatter/darkfam_ecology/endemism/query.faa` (23,382 reps; same as MGnify test).
- Outputs (pulled here): `global_specificity_summary.tsv` (GEM/SPIRE), `desert_specificity_summary.tsv` (deserts).

## Results (2026-06-26)
% of the 23,382 dark families without a detected homolog in each named catalogue — see `MULTICATALOG.md` +
`fig_darkmatter_global_specificity`:

| Catalogue | depth | absent |
|-----------|------:|-------:|
| MGnify90 (all biomes) | 718M | 36.6% |
| GEM (52k MAGs) | 123M | **47.4%** |
| arid soils | 31M | 41.7% |
| Mojave | 54M | 51.9% |
| Antarctic cold-desert | 49M | 55.9% |
| trace-gas-autotrophy desert | 7M | 62.8% |
| Atacama | 6M | 89.1% |
| **ALL deserts (union)** | 147M | **31.9%** |
| SPIRE (rep. genomes) | — | finishing (job 47810126) |

**Interpretation:** 47.4% of the EQ dark families have no detected homolog in GEM, 42–89% have no detected
homolog in each individual desert catalogue, and **31.9% have no detected homolog in the union of all five
deserts**. The no-hit fraction remains 32–47% in the deepest named comparisons (MGnify 718M=37%, GEM 123M=47%,
combined-desert 147M=32%). This persistence prioritises families for broader searches; it cannot distinguish
biogeographic restriction from incomplete reference sampling. Shallow catalogues (Atacama 6M, autotroph 7M)
inflate the fraction, so the deeper comparisons are the descriptive anchors.

**Cross-catalogue intersection:** pooling the no-hit sets, **5,011 / 23,382 = 21.4%** of EQ dark families lack a
detected homolog in the named MGnify90, GEM and five-desert searches (a hit occurs in at least one named
catalogue for 18,371). This database-relative intersection can only decrease as additional references are
searched. The historical no-hit/hit dN/dS contrast is descriptive and does not establish adaptation or
geographic restriction.
