# Proteome review and integration — 8 November 2026

The defensible story is a conserved unknown repertoire in widespread subsurface archaea, placed in the context of functional differentiation across Rub’ al-Khali soil compartments. These results strengthen an ecological and evolutionary paper; they do not yet establish a new biochemical mechanism or justify a Nature-level mechanistic claim.

## Evidence and disposition

| Result | Evidence | Integration |
|---|---|---|
| Surface betaine uptake | HMM mapped-community excess 7.55 percentage points, spatial-block 95% interval 5.11–10.16; independent KO projection agrees | Main results and new Figure 1 |
| Surface peroxide defence | AhpCF excess 6.95 points, interval 4.04–9.65; KO projection agrees | Main results and new Figure 1 |
| Community replacement | Between shared genera contributes 6.12 and 5.28 points respectively; within-genus contributions small | Interpret gene-content differences primarily through community composition, not induced gene expression |
| Subsurface Thermoplasmatota | Competitive reference includes 5,438 EQ and 301 external genomes; 193 profiles, 68 visits, 36 sites. At 20% breadth, 35/36 sites favour subsurface; mean excess 0.528 points, interval 0.258–0.903 | Main results, methods and Supplementary Note S6 |
| Breadth sensitivity | At 50% breadth, mean excess 0.472 points, interval 0.210–0.816 | Retained as robustness evidence |
| Trehalose synthesis | HMM enrichment disagrees with KO projection | Not a headline claim; disagreement documented in S6 |
| Protein-maintenance panel | Direction and support depend on annotation and denominator | No claim that this panel explains the subsurface niche |
| Dark-protein candidates | 7,381 unique genome–reaction–gene proposals; 3,766 pass priority filters, including 3,387 dark proposals | Candidate evidence only, not validated enzyme assignments |
| Context benchmark | Mean top-five recovery improvement over within-contig null is 4.25 percentage points; genome-bootstrap interval 3.69–4.83 | Modest ranking evidence conditional on known pathway context, not calibrated dark-protein accuracy |

The HMM and KO routes are distinct annotation procedures, not interchangeable implementations. Functional contrasts also retain positive intervals with the original-read denominator, higher-quality genomes and the amplicon-overlap cohort. Taxonomic decomposition is exact, but 27.5% of mean mapped abundance lacks genus assignments and is represented as individual genome groups, whose within-group component is zero by construction.

## Remaining limitations and submission work

1. **Function remains the main scientific gap.** Conserved neighbourhoods, structure and gap-filling prioritise candidates but do not determine biochemical activity. A small number of convergent, independently supported computational case studies would strengthen the story; this review does not prescribe new wet-lab experiments.
2. **Spatial replication is incomplete.** Rund’s full INLA input bundle, fitted models and environment remain unavailable. Current spline-based sensitivities are exploratory; the manuscript explicitly distinguishes them from INLA/SPDE replication.
3. **Mapping claims remain aggregate.** The MAPQ20 raw mean-coverage sensitivity has a wide block interval crossing zero and is not scaled to the original read denominator. It cannot be compared numerically with the relative-abundance effect. Stronger species-level claims require further mapping and denominator checks.
4. **ANI components are provisional.** The reciprocal 95% ANI / 50% aligned-fraction graph has 825 components for 892 genomes, 66 nonsingletons, and no mixed EQ/external component. Connected components can chain and are not a validated species delineation.
5. **Different cohorts must remain explicit.** Original archaeal occurrence analysis has 550 usable profiles and 37 paired sites; the audited core competitive analysis has 193 profiles and 36 paired sites. The new functional analysis uses 546 nonzero core profiles, with 442 in the amplicon intersection. These are not contradictory totals.
6. **Existing submission blockers persist.** Data/reviewer access, MAG deposition, the versioned code/results release and historical raw-input provenance gaps need resolution. Existing orange numbers mark computed results not yet registered in the prior evidence system; these markings were preserved. A full re-audit of every older analysis was outside this integration.
7. **Paper boundaries.** Rund’s ecology manuscript was not edited. This proteome paper uses ecological context to explain the importance of unresolved archaeal genes; it does not claim ownership of the entire community ecology story or force a defensome merger.

## Changes and validation

Reframed title, abstract, introduction and discussion; added functional contrasts, competitive mapping, reproducible methods, Figure 1 and Supplementary Note S6; shortened the dense coding-validation paragraph while retaining numerical detail in the supplement. Repaired four overflowing supplementary tables and an oversized supplementary figure. The microTrait citation was checked against the publisher article (DOI 10.3389/fbinf.2022.918853).

The standalone manuscript build passes without undefined references/citations or overfull boxes. New figure and representative main/supplementary pages were rendered and visually checked. The final compact-review Slurm job 15852347 and figure job 15852345 completed with exit code 0. This integration did not rerun the raw-data pipeline. Pre-edit manuscript files and the earlier package snapshot are backed up in this directory; no commit, push or external distribution was performed.

## Evidence locations

- `../proteome-integration-20261106/results/`: paired HMM/KO tests, turnover decomposition, spatial sensitivities, candidate audit and masked benchmark.
- `results/mapping_group_review.tsv`: competitive mapping aggregate and breadth/MAPQ sensitivities.
- `results/ani_component_review.tsv`: provisional component composition.
- `../proteome-integration-20261106/code/`: analysis and plotting scripts, including `review_mapping.py` and `plot_integration.py`.
- Remote integration root: `/scratch/project/k10018/eq-proteome-integration-20261106`.
- Remote mapping root: `/scratch/project/k10018/eq-archaeal-function-20261106`.
