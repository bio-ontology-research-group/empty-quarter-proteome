# A1 — Phylogenetic placement of atmospheric-metabolism markers (reviewer-armor)

Confirms the EQ hydrogenase / CO-dehydrogenase hits are the **bona-fide high-affinity atmospheric** enzyme
groups, not "any old" homologs — pre-empting the "group 1h is a wastebasket" objection. EQ marker reps were
classified against the **HydDB** reference set (Søndergaard 2016; group-labelled `[NiFe]_Group_*`) and a form-I
CO-dehydrogenase reference set, and placed on phylogenetic trees (`hydB.nwk`, `coxL.nwk` = EQ reps + references).

## [NiFe]-hydrogenase (hydB) — HydDB group breakdown (proteins / MAGs)
| Group | proteins | MAGs | note |
|-------|---------:|-----:|------|
| **1l** | 5,292 | 2,609 | high-affinity, atmospheric H₂ (Actinobacteria-type) |
| **1h** | 986 | 893 | canonical high-affinity atmospheric H₂ uptake |
| 2a | 246 | 235 | (cyanobacterial/uptake) |
| 1d,1c,1f,2b,1a,1i,1e,2c,1b | <60 total | — | minor |
→ **Atmospheric groups 1l+1h dominate** (6,278 of 6,602 hydB; 3,386 MAGs); the assignment is phylogenetic, not
just HMM-hit.

## CO-dehydrogenase (coxL) — form classification
| Form | proteins | MAGs | note |
|------|---------:|-----:|------|
| **form I** | 5,645 | 3,828 | aerobic, atmospheric CO oxidation (AYRGAGR motif) |
| ambiguous | 8,003 | 4,284 | form-I/II boundary — not counted as confirmed |
| form II/other | 27 | 27 | non-atmospheric |

## Guild core
Both confirmed atmospheric H₂ (group 1h/1l) AND form-I CO: **1,204 MAGs (~8%)** (`guild_core_mags.txt`).

## Status
Trees + classification DONE (workbench `06_atmospheric/classify/`; pulled here). Remaining: render the
publication tree figure (EQ reps coloured by HydDB group, reference clades labelled) — see figure task.
**Gap:** V/Fe-nitrogenase (vnf/anf) phylogeny not yet built (`07_nitrogen` empty) — optional A4 add.
