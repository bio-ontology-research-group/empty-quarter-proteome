# RQ-D6 — Result

## ANSWER
**Rhodopsin-based light-energy harvesting is a rare, niche-agnostic trait in the Empty Quarter — it
does NOT track a light niche.** Only **337 of 15,469 MAGs (2.18%)** encode a genuine microbial
type-1 rhodopsin (Bac_rhodopsin / PF01036), 391 opsin proteins in total. Carriage is essentially
**flat across compartments**: Surface 2.50% (56/2,241), Deep 2.34% (91/3,890), PlantRoot 2.04%
(187/9,170). Surface-vs-Deep enrichment is **not significant** (OR = 1.07, Fisher exact **p = 0.73**).

This is a clean **negative / contrast** result. The light-dependent **photoprotection** machinery in
the same MAG set partitions sharply by sun exposure — photolyase 22.4% Surface vs 13.1% Deep
(`results/uv_survival/compartment.tsv`) — whereas **light harvesting via rhodopsins does not**.
Empty Quarter soil bacteria invest heavily in *surviving* light (UV/desiccation defence) but only a
small minority attempt to *harvest* it, and those few are scattered across depth rather than
concentrated at the illuminated surface.

## Key numbers
| metric | value |
|---|---|
| MAGs with ≥1 type-1 rhodopsin | 337 / 15,469 = **2.18%** |
| rhodopsin proteins | 391 |
| Surface carriage | 2.50% (56/2,241) |
| Deep carriage | 2.34% (91/3,890) |
| PlantRoot carriage | 2.04% (187/9,170) |
| Surface vs Deep | OR 1.07, Fisher **p = 0.73** (n.s.) |
| photolyase reference (light-tracking) | 22.4% Surface vs 13.1% Deep |

### Sub-types (coarse; KO/name-resolved only)
| sub-type | proteins | MAGs |
|---|---|---|
| type-1 rhodopsin (subtype unresolved) | 299 | 280 |
| sensory-rhodopsin-type (K04642/3) | 57 | 51 |
| bacteriorhodopsin-type (K04641) | 35 | 35 |

Most opsins (299/391) carry only the generic Bac_rhodopsin Pfam with no informative KO, so EggNOG
cannot separate light-driven proton-pump proteorhodopsins from sensory rhodopsins at this resolution.
Notably, a substantial minority (~57 proteins) carry **sensory**-rhodopsin KOs — i.e. light *sensing*
rather than energy capture — further weakening any "life on light" energy-harvesting interpretation.

## Interpretation
Microbial rhodopsins are light-driven retinal proton pumps that can supplement ATP synthesis and are
common in sunlit aquatic and surface environments (Béjà et al. 2000; Rozenberg et al. 2021,
`rozenberg2021rhodopsins`). In these arid Empty Quarter soils they are present but **uncommon and not
surface-enriched**, consistent with rhodopsin phototrophy being a marginal supplementary strategy here
rather than a dominant, niche-structured energy source — in contrast to the genome-wide UV/desiccation
hardening (RQ on UV survival) and the atmospheric H₂/CO chemolithotrophy pillar.

## Honest caveats
- **Annotation source is homology-based EggNOG/Pfam only.** No HMMER-confirmed retinal-binding-lysine
  check and no phylogenetic placement; subtype (proteorhodopsin vs sensory) calls are coarse and very
  divergent soil opsins could be missed. Prevalence is robust (PF01036 is specific), but the subtype
  split and any proton-pump claim need a targeted Pfam/HMMER + retinal-residue validation.
- **Curated 7 eukaryotic GPCR "rhodopsin-family" false positives** (Pfam 7tm_1/Orexin_rec2/adh_short,
  all from 4 PlantRoot MAGs) were removed before counting — they are animal/plant opsins
  (contamination), not microbial type-1 rhodopsins.
- The compartment denominators use **all** MAGs per compartment (rhodopsin is rare); the null result
  is not a denominator artefact (Surface 2.50% vs Deep 2.34% are near-identical raw fractions).

## Files
- `rhodopsin_census.tsv` — prevalence, sub-types, compartment split, Fisher test.
- `rhodopsin_hits.tsv` — per-protein census pulled from IBEX (398 raw rows; 391 retained after QC).
- `fig_rhodopsin.{pdf,png}` — prevalence (left) + compartment split vs photolyase reference (right).
- `run.py` — documents the IBEX census commands and reproduces the analysis.
- IBEX: `scripts/ibex/D6_rho_census.sh` (job 47843530); `06_metabolism/rhodopsin/rhodopsin_hits.tsv`.

**Status: RQ-D6 COMPLETE.** Real EggNOG annotations located and used; no fabricated counts.
