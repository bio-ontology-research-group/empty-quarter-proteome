# Bounded F3 Genomic-Context Test: Synthesis & Hypothesis Evaluation (FULL)
**Date:** 2026-10-15 07:34:59 UTC  
**Total Evaluated Rhodopsin Loci:** 205 (Type-1: 140, Heliorhodopsin: 65)  
**Contig Flank Completeness:** Full-flank: 91 (44.4%), Edge-truncated: 114 (55.6%)  
**Matched Controls:** Same-contig: 115 (56.1%), Same-MAG surrogate: 90 (43.9%)  

---

## 1. Executive Summary & Biological Adjudication

The F3 discovery hypothesis asks whether Empty Quarter (EQ) rhodopsins represent divergent photosensory
systems with conserved sensory/signaling genomic contexts (distinguishable from ion pumps), as opposed to
simply additional sequence hits. Following Bulzu et al. 2021 (*mSphere*) and microbial sensory transducer
precedents, we tested a same-contig +/- 5kb genomic window for source-verified Pfam signaling domains
(histidine kinases, response regulators, MCP chemoreceptor transducers, HAMP/CheW linkers, PAS/GAF sensor
domains, and cyclic di-GMP enzymes) using curated gathering thresholds (`--cut_ga`) against matched background controls.

### Key Findings:
#### Type-1 Microbial Rhodopsin (n = 140):
- **Any Sensory/Signaling Marker:** Rhodopsin neighborhood = 66/140 (47.14%) vs Matched Control = 15/140 (10.71%)  
  * Fisher's Exact p = 1.28e-11, Odds Ratio = 7.432 (95% CI: 3.958 - 13.956)  
  * Paired McNemar p = 4.29e-13 (discordant: rho-only = 54, ctrl-only = 3)  
- **Two-Component Kinase (HisKA/HATPase):** Rho = 9 (6.43%) vs Ctrl = 6 (4.29%) (p = 0.5971)  
- **Two-Component Regulator (Response_reg/Trans_reg_C):** Rho = 16 (11.43%) vs Ctrl = 7 (5.0%) (p = 0.07975)  
- **Transducer/Linker (MCP/HAMP/CheW):** Rho = 42 (30.0%) vs Ctrl = 3 (2.14%) (p = 3.634e-11)  

#### Heliorhodopsin Type-3 (n = 65):
- **Any Sensory/Signaling Marker:** Rhodopsin neighborhood = 8/65 (12.31%) vs Matched Control = 14/65 (21.54%)  
  * Fisher's Exact p = 0.2418, Odds Ratio = 0.511 (95% CI: 0.198 - 1.318)  
  * Paired McNemar p = 0.2632 (discordant: rho-only = 7, ctrl-only = 13)  
- **Two-Component Kinase (HisKA/HATPase):** Rho = 2 (3.08%) vs Ctrl = 7 (10.77%) (p = 0.1643)  
- **Two-Component Regulator (Response_reg/Trans_reg_C):** Rho = 3 (4.62%) vs Ctrl = 9 (13.85%) (p = 0.1271)  
- **Transducer/Linker (MCP/HAMP/CheW):** Rho = 0 (0.0%) vs Ctrl = 0 (0.0%) (p = 1)  

#### All Complete EQ Rhodopsins (n = 205):
- **Any Sensory/Signaling Marker:** Rhodopsin neighborhood = 74/205 (36.1%) vs Matched Control = 29/205 (14.15%)  
  * Fisher's Exact p = 3.955e-07, Odds Ratio = 3.428 (95% CI: 2.11 - 5.57)  
  * Paired McNemar p = 2.418e-07 (discordant: rho-only = 61, ctrl-only = 16)  
- **Two-Component Kinase (HisKA/HATPase):** Rho = 11 (5.37%) vs Ctrl = 13 (6.34%) (p = 0.8339)  
- **Two-Component Regulator (Response_reg/Trans_reg_C):** Rho = 19 (9.27%) vs Ctrl = 16 (7.8%) (p = 0.7242)  
- **Transducer/Linker (MCP/HAMP/CheW):** Rho = 42 (20.49%) vs Ctrl = 3 (1.46%) (p = 1.128e-10)  

---

## 2. Comprehensive Results Table

| Stratum | Subset | Marker | Rhodopsin Hits | Control Hits | Odds Ratio (95% CI) | Fisher p | McNemar p |
|---|---|---|---|---|---|---|---|
| All_Complete_EQ_Rhodopsins | all_loci | Any Sensory / Signaling | 74/205 (36.1%) | 29/205 (14.15%) | 3.428 (2.11-5.57) | 3.955e-07 | 2.418e-07 |
| All_Complete_EQ_Rhodopsins | all_loci | TCS Histidine Kinase | 11/205 (5.37%) | 13/205 (6.34%) | 0.837 (0.366-1.915) | 0.8339 | 0.8318 |
| All_Complete_EQ_Rhodopsins | all_loci | TCS Response Regulator | 19/205 (9.27%) | 16/205 (7.8%) | 1.207 (0.602-2.418) | 0.7242 | 0.7201 |
| All_Complete_EQ_Rhodopsins | all_loci | Transducer / Linker (MCP/HAMP/CheW) | 42/205 (20.49%) | 3/205 (1.46%) | 17.35 (5.282-56.994) | 1.128e-10 | 3.638e-12 |
| All_Complete_EQ_Rhodopsins | all_loci | Sensor Domain (PAS/GAF) | 23/205 (11.22%) | 16/205 (7.8%) | 1.493 (0.764-2.917) | 0.3125 | 0.324 |
| All_Complete_EQ_Rhodopsins | all_loci | Second Messenger (GGDEF/EAL) | 6/205 (2.93%) | 2/205 (0.98%) | 3.06 (0.61-15.344) | 0.2842 | 0.2891 |
| All_Complete_EQ_Rhodopsins | full_flank_only | Any Sensory / Signaling | 35/91 (38.46%) | 15/91 (16.48%) | 3.167 (1.578-6.354) | 0.00146 | 0.0008214 |
| All_Complete_EQ_Rhodopsins | full_flank_only | TCS Histidine Kinase | 6/91 (6.59%) | 5/91 (5.49%) | 1.214 (0.357-4.129) | 1 | 1 |
| All_Complete_EQ_Rhodopsins | full_flank_only | TCS Response Regulator | 10/91 (10.99%) | 8/91 (8.79%) | 1.281 (0.481-3.409) | 0.8046 | 0.7905 |
| All_Complete_EQ_Rhodopsins | full_flank_only | Transducer / Linker (MCP/HAMP/CheW) | 18/91 (19.78%) | 2/91 (2.2%) | 10.973 (2.465-48.848) | 0.0001904 | 3.052e-05 |
| All_Complete_EQ_Rhodopsins | full_flank_only | Sensor Domain (PAS/GAF) | 12/91 (13.19%) | 10/91 (10.99%) | 1.23 (0.503-3.01) | 0.8207 | 0.8238 |
| All_Complete_EQ_Rhodopsins | full_flank_only | Second Messenger (GGDEF/EAL) | 5/91 (5.49%) | 1/91 (1.1%) | 5.233 (0.599-45.706) | 0.2109 | 0.2188 |
| All_Complete_EQ_Rhodopsins | edge_truncated_only | Any Sensory / Signaling | 39/114 (34.21%) | 14/114 (12.28%) | 3.714 (1.882-7.332) | 0.0001368 | 0.0001702 |
| All_Complete_EQ_Rhodopsins | edge_truncated_only | TCS Histidine Kinase | 5/114 (4.39%) | 8/114 (7.02%) | 0.608 (0.193-1.917) | 0.5696 | 0.5488 |
| All_Complete_EQ_Rhodopsins | edge_truncated_only | TCS Response Regulator | 9/114 (7.89%) | 8/114 (7.02%) | 1.136 (0.422-3.056) | 1 | 1 |
| All_Complete_EQ_Rhodopsins | edge_truncated_only | Transducer / Linker (MCP/HAMP/CheW) | 24/114 (21.05%) | 1/114 (0.88%) | 30.133 (3.999-227.039) | 4.442e-07 | 2.384e-07 |
| All_Complete_EQ_Rhodopsins | edge_truncated_only | Sensor Domain (PAS/GAF) | 11/114 (9.65%) | 6/114 (5.26%) | 1.922 (0.686-5.388) | 0.3134 | 0.3323 |
| All_Complete_EQ_Rhodopsins | edge_truncated_only | Second Messenger (GGDEF/EAL) | 1/114 (0.88%) | 1/114 (0.88%) | 1.0 (0.062-16.185) | 1 | 1 |
| All_Complete_EQ_Rhodopsins | mag_dereplicated | Any Sensory / Signaling | 55/170 (32.35%) | 23/170 (13.53%) | 3.057 (1.773-5.269) | 5.346e-05 | 5.778e-05 |
| All_Complete_EQ_Rhodopsins | mag_dereplicated | TCS Histidine Kinase | 10/170 (5.88%) | 10/170 (5.88%) | 1.0 (0.405-2.468) | 1 | 1 |
| All_Complete_EQ_Rhodopsins | mag_dereplicated | TCS Response Regulator | 16/170 (9.41%) | 14/170 (8.24%) | 1.158 (0.546-2.453) | 0.8487 | 0.8506 |
| All_Complete_EQ_Rhodopsins | mag_dereplicated | Transducer / Linker (MCP/HAMP/CheW) | 26/170 (15.29%) | 1/170 (0.59%) | 30.514 (4.09-227.659) | 1.597e-07 | 5.96e-08 |
| All_Complete_EQ_Rhodopsins | mag_dereplicated | Sensor Domain (PAS/GAF) | 21/170 (12.35%) | 12/170 (7.06%) | 1.856 (0.882-3.904) | 0.1418 | 0.1496 |
| All_Complete_EQ_Rhodopsins | mag_dereplicated | Second Messenger (GGDEF/EAL) | 6/170 (3.53%) | 2/170 (1.18%) | 3.073 (0.611-15.447) | 0.2832 | 0.2891 |
| Type-1_Microbial_Rhodopsin | all_loci | Any Sensory / Signaling | 66/140 (47.14%) | 15/140 (10.71%) | 7.432 (3.958-13.956) | 1.28e-11 | 4.29e-13 |
| Type-1_Microbial_Rhodopsin | all_loci | TCS Histidine Kinase | 9/140 (6.43%) | 6/140 (4.29%) | 1.534 (0.531-4.432) | 0.5971 | 0.6072 |
| Type-1_Microbial_Rhodopsin | all_loci | TCS Response Regulator | 16/140 (11.43%) | 7/140 (5.0%) | 2.452 (0.976-6.16) | 0.07975 | 0.06357 |
| Type-1_Microbial_Rhodopsin | all_loci | Transducer / Linker (MCP/HAMP/CheW) | 42/140 (30.0%) | 3/140 (2.14%) | 19.571 (5.897-64.959) | 3.634e-11 | 3.638e-12 |
| Type-1_Microbial_Rhodopsin | all_loci | Sensor Domain (PAS/GAF) | 19/140 (13.57%) | 11/140 (7.86%) | 1.841 (0.842-4.029) | 0.1754 | 0.1849 |
| Type-1_Microbial_Rhodopsin | all_loci | Second Messenger (GGDEF/EAL) | 3/140 (2.14%) | 1/140 (0.71%) | 3.044 (0.313-29.624) | 0.6223 | 0.625 |
| Type-1_Microbial_Rhodopsin | full_flank_only | Any Sensory / Signaling | 31/68 (45.59%) | 9/68 (13.24%) | 5.492 (2.351-12.829) | 5.738e-05 | 2.98e-06 |
| Type-1_Microbial_Rhodopsin | full_flank_only | TCS Histidine Kinase | 6/68 (8.82%) | 3/68 (4.41%) | 2.097 (0.502-8.752) | 0.4927 | 0.5078 |
| Type-1_Microbial_Rhodopsin | full_flank_only | TCS Response Regulator | 8/68 (11.76%) | 3/68 (4.41%) | 2.889 (0.732-11.398) | 0.207 | 0.125 |
| Type-1_Microbial_Rhodopsin | full_flank_only | Transducer / Linker (MCP/HAMP/CheW) | 18/68 (26.47%) | 2/68 (2.94%) | 11.88 (2.634-53.585) | 0.0001416 | 3.052e-05 |
| Type-1_Microbial_Rhodopsin | full_flank_only | Sensor Domain (PAS/GAF) | 10/68 (14.71%) | 7/68 (10.29%) | 1.502 (0.536-4.211) | 0.6052 | 0.6072 |
| Type-1_Microbial_Rhodopsin | full_flank_only | Second Messenger (GGDEF/EAL) | 3/68 (4.41%) | 1/68 (1.47%) | 3.092 (0.314-30.5) | 0.6194 | 0.625 |
| Type-1_Microbial_Rhodopsin | edge_truncated_only | Any Sensory / Signaling | 35/72 (48.61%) | 6/72 (8.33%) | 10.405 (4.004-27.041) | 8.825e-08 | 1.309e-07 |
| Type-1_Microbial_Rhodopsin | edge_truncated_only | TCS Histidine Kinase | 3/72 (4.17%) | 3/72 (4.17%) | 1.0 (0.195-5.128) | 1 | 1 |
| Type-1_Microbial_Rhodopsin | edge_truncated_only | TCS Response Regulator | 8/72 (11.11%) | 4/72 (5.56%) | 2.125 (0.61-7.401) | 0.3666 | 0.3877 |
| Type-1_Microbial_Rhodopsin | edge_truncated_only | Transducer / Linker (MCP/HAMP/CheW) | 24/72 (33.33%) | 1/72 (1.39%) | 35.5 (4.645-271.293) | 1.831e-07 | 2.384e-07 |
| Type-1_Microbial_Rhodopsin | edge_truncated_only | Sensor Domain (PAS/GAF) | 9/72 (12.5%) | 4/72 (5.56%) | 2.429 (0.712-8.281) | 0.2442 | 0.2668 |
| Type-1_Microbial_Rhodopsin | edge_truncated_only | Second Messenger (GGDEF/EAL) | 0/72 (0.0%) | 0/72 (0.0%) | 1.0 (0.02-51.085) | 1 | 1 |
| Type-1_Microbial_Rhodopsin | mag_dereplicated | Any Sensory / Signaling | 47/111 (42.34%) | 10/111 (9.01%) | 7.417 (3.5-15.717) | 1.246e-08 | 3.021e-09 |
| Type-1_Microbial_Rhodopsin | mag_dereplicated | TCS Histidine Kinase | 8/111 (7.21%) | 4/111 (3.6%) | 2.078 (0.607-7.111) | 0.3742 | 0.3877 |
| Type-1_Microbial_Rhodopsin | mag_dereplicated | TCS Response Regulator | 13/111 (11.71%) | 5/111 (4.5%) | 2.812 (0.967-8.177) | 0.08281 | 0.07681 |
| Type-1_Microbial_Rhodopsin | mag_dereplicated | Transducer / Linker (MCP/HAMP/CheW) | 26/111 (23.42%) | 1/111 (0.9%) | 33.647 (4.476-252.957) | 8.923e-08 | 5.96e-08 |
| Type-1_Microbial_Rhodopsin | mag_dereplicated | Sensor Domain (PAS/GAF) | 17/111 (15.32%) | 7/111 (6.31%) | 2.687 (1.067-6.765) | 0.04991 | 0.05248 |
| Type-1_Microbial_Rhodopsin | mag_dereplicated | Second Messenger (GGDEF/EAL) | 3/111 (2.7%) | 1/111 (0.9%) | 3.056 (0.313-29.835) | 0.6216 | 0.625 |
| Heliorhodopsin_Type-3 | all_loci | Any Sensory / Signaling | 8/65 (12.31%) | 14/65 (21.54%) | 0.511 (0.198-1.318) | 0.2418 | 0.2632 |
| Heliorhodopsin_Type-3 | all_loci | TCS Histidine Kinase | 2/65 (3.08%) | 7/65 (10.77%) | 0.263 (0.053-1.318) | 0.1643 | 0.125 |
| Heliorhodopsin_Type-3 | all_loci | TCS Response Regulator | 3/65 (4.62%) | 9/65 (13.85%) | 0.301 (0.078-1.168) | 0.1271 | 0.146 |
| Heliorhodopsin_Type-3 | all_loci | Transducer / Linker (MCP/HAMP/CheW) | 0/65 (0.0%) | 0/65 (0.0%) | 1.0 (0.02-51.159) | 1 | 1 |
| Heliorhodopsin_Type-3 | all_loci | Sensor Domain (PAS/GAF) | 4/65 (6.15%) | 5/65 (7.69%) | 0.787 (0.202-3.073) | 1 | 1 |
| Heliorhodopsin_Type-3 | all_loci | Second Messenger (GGDEF/EAL) | 3/65 (4.62%) | 1/65 (1.54%) | 3.097 (0.314-30.581) | 0.6191 | 0.625 |
| Heliorhodopsin_Type-3 | full_flank_only | Any Sensory / Signaling | 4/23 (17.39%) | 6/23 (26.09%) | 0.596 (0.144-2.479) | 0.7222 | 0.7539 |
| Heliorhodopsin_Type-3 | full_flank_only | TCS Histidine Kinase | 0/23 (0.0%) | 2/23 (8.7%) | 0.183 (0.008-4.03) | 0.4889 | 0.5 |
| Heliorhodopsin_Type-3 | full_flank_only | TCS Response Regulator | 2/23 (8.7%) | 5/23 (21.74%) | 0.343 (0.059-1.986) | 0.414 | 0.4531 |
| Heliorhodopsin_Type-3 | full_flank_only | Transducer / Linker (MCP/HAMP/CheW) | 0/23 (0.0%) | 0/23 (0.0%) | 1.0 (0.019-52.535) | 1 | 1 |
| Heliorhodopsin_Type-3 | full_flank_only | Sensor Domain (PAS/GAF) | 2/23 (8.7%) | 3/23 (13.04%) | 0.635 (0.096-4.207) | 1 | 1 |
| Heliorhodopsin_Type-3 | full_flank_only | Second Messenger (GGDEF/EAL) | 2/23 (8.7%) | 0/23 (0.0%) | 5.465 (0.248-120.377) | 0.4889 | 0.5 |
| Heliorhodopsin_Type-3 | edge_truncated_only | Any Sensory / Signaling | 4/42 (9.52%) | 8/42 (19.05%) | 0.447 (0.124-1.619) | 0.3502 | 0.3438 |
| Heliorhodopsin_Type-3 | edge_truncated_only | TCS Histidine Kinase | 2/42 (4.76%) | 5/42 (11.9%) | 0.37 (0.068-2.025) | 0.4326 | 0.375 |
| Heliorhodopsin_Type-3 | edge_truncated_only | TCS Response Regulator | 1/42 (2.38%) | 4/42 (9.52%) | 0.232 (0.025-2.166) | 0.3597 | 0.375 |
| Heliorhodopsin_Type-3 | edge_truncated_only | Transducer / Linker (MCP/HAMP/CheW) | 0/42 (0.0%) | 0/42 (0.0%) | 1.0 (0.019-51.573) | 1 | 1 |
| Heliorhodopsin_Type-3 | edge_truncated_only | Sensor Domain (PAS/GAF) | 2/42 (4.76%) | 2/42 (4.76%) | 1.0 (0.134-7.451) | 1 | 1 |
| Heliorhodopsin_Type-3 | edge_truncated_only | Second Messenger (GGDEF/EAL) | 1/42 (2.38%) | 1/42 (2.38%) | 1.0 (0.06-16.535) | 1 | 1 |
| Heliorhodopsin_Type-3 | mag_dereplicated | Any Sensory / Signaling | 8/62 (12.9%) | 13/62 (20.97%) | 0.558 (0.213-1.461) | 0.3384 | 0.3593 |
| Heliorhodopsin_Type-3 | mag_dereplicated | TCS Histidine Kinase | 2/62 (3.23%) | 6/62 (9.68%) | 0.311 (0.06-1.606) | 0.2727 | 0.2188 |
| Heliorhodopsin_Type-3 | mag_dereplicated | TCS Response Regulator | 3/62 (4.84%) | 9/62 (14.52%) | 0.299 (0.077-1.165) | 0.1261 | 0.146 |
| Heliorhodopsin_Type-3 | mag_dereplicated | Transducer / Linker (MCP/HAMP/CheW) | 0/62 (0.0%) | 0/62 (0.0%) | 1.0 (0.02-51.195) | 1 | 1 |
| Heliorhodopsin_Type-3 | mag_dereplicated | Sensor Domain (PAS/GAF) | 4/62 (6.45%) | 5/62 (8.06%) | 0.786 (0.201-3.078) | 1 | 1 |
| Heliorhodopsin_Type-3 | mag_dereplicated | Second Messenger (GGDEF/EAL) | 3/62 (4.84%) | 1/62 (1.61%) | 3.102 (0.314-30.67) | 0.6188 | 0.625 |

---

## 3. Methodological Audit & Pre-registration Criteria

1. **Pre-registered Thresholds:** All searches used curated Pfam Gathering Thresholds (`--cut_ga`). No ad-hoc e-value or score tuning was applied.
2. **Self-Exclusion:** The candidate rhodopsin gene was strictly excluded from its own neighborhood.
3. **Tandem Rhodopsin Accounting:** Tandem rhodopsin duplicates were tracked and distinguished from signaling transducers.
4. **Contig-Edge Limitations:** Every contig boundary was tracked. Full-flank and edge-truncated loci were analyzed separately; the absence of markers in edge-truncated contigs was explicitly recorded as constrained/untested rather than assumed negative.
5. **Matched Controls:** Background control loci were sampled on the same contig (outside the +/- 5kb rhodopsin window) with matching flank parameters, falling back to same-MAG surrogate controls only when contig span was insufficient.

---

## 4. Scientific Disposition & Conclusion

Genomic context alone cannot definitively prove biochemical ion pumping or light sensing, but it establishes
whether the candidate lineages recurrently associate with classical sensory transducers or two-component systems.
Based on these bounded empirical data, this report provides the exact evidentiary basis to either support or
appropriately downgrade the sensory rhodopsin hypothesis in the Empty Quarter manuscript without overclaiming.