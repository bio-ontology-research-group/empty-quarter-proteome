# Bounded F3 Genomic-Context Test: Synthesis & Hypothesis Evaluation (PILOT)
**Date:** 2026-10-15 07:34:26 UTC  
**Total Evaluated Rhodopsin Loci:** 20 (Type-1: 10, Heliorhodopsin: 10)  
**Contig Flank Completeness:** Full-flank: 7 (35.0%), Edge-truncated: 13 (65.0%)  
**Matched Controls:** Same-contig: 10 (50.0%), Same-MAG surrogate: 10 (50.0%)  

---

## 1. Executive Summary & Biological Adjudication

The F3 discovery hypothesis asks whether Empty Quarter (EQ) rhodopsins represent divergent photosensory
systems with conserved sensory/signaling genomic contexts (distinguishable from ion pumps), as opposed to
simply additional sequence hits. Following Bulzu et al. 2021 (*mSphere*) and microbial sensory transducer
precedents, we tested a same-contig +/- 5kb genomic window for source-verified Pfam signaling domains
(histidine kinases, response regulators, MCP chemoreceptor transducers, HAMP/CheW linkers, PAS/GAF sensor
domains, and cyclic di-GMP enzymes) using curated gathering thresholds (`--cut_ga`) against matched background controls.

### Key Findings:
#### Type-1 Microbial Rhodopsin (n = 10):
- **Any Sensory/Signaling Marker:** Rhodopsin neighborhood = 5/10 (50.0%) vs Matched Control = 1/10 (10.0%)  
  * Fisher's Exact p = 0.1409, Odds Ratio = 9.0 (95% CI: 0.809 - 100.143)  
  * Paired McNemar p = 0.125 (discordant: rho-only = 4, ctrl-only = 0)  
- **Two-Component Kinase (HisKA/HATPase):** Rho = 1 (10.0%) vs Ctrl = 0 (0.0%) (p = 1)  
- **Two-Component Regulator (Response_reg/Trans_reg_C):** Rho = 2 (20.0%) vs Ctrl = 1 (10.0%) (p = 1)  
- **Transducer/Linker (MCP/HAMP/CheW):** Rho = 3 (30.0%) vs Ctrl = 0 (0.0%) (p = 0.2105)  

#### Heliorhodopsin Type-3 (n = 10):
- **Any Sensory/Signaling Marker:** Rhodopsin neighborhood = 1/10 (10.0%) vs Matched Control = 2/10 (20.0%)  
  * Fisher's Exact p = 1, Odds Ratio = 0.444 (95% CI: 0.034 - 5.88)  
  * Paired McNemar p = 1 (discordant: rho-only = 1, ctrl-only = 2)  
- **Two-Component Kinase (HisKA/HATPase):** Rho = 1 (10.0%) vs Ctrl = 0 (0.0%) (p = 1)  
- **Two-Component Regulator (Response_reg/Trans_reg_C):** Rho = 0 (0.0%) vs Ctrl = 1 (10.0%) (p = 1)  
- **Transducer/Linker (MCP/HAMP/CheW):** Rho = 0 (0.0%) vs Ctrl = 0 (0.0%) (p = 1)  

#### All Complete EQ Rhodopsins (n = 20):
- **Any Sensory/Signaling Marker:** Rhodopsin neighborhood = 6/20 (30.0%) vs Matched Control = 3/20 (15.0%)  
  * Fisher's Exact p = 0.4506, Odds Ratio = 2.429 (95% CI: 0.512 - 11.511)  
  * Paired McNemar p = 0.4531 (discordant: rho-only = 5, ctrl-only = 2)  
- **Two-Component Kinase (HisKA/HATPase):** Rho = 2 (10.0%) vs Ctrl = 0 (0.0%) (p = 0.4872)  
- **Two-Component Regulator (Response_reg/Trans_reg_C):** Rho = 2 (10.0%) vs Ctrl = 2 (10.0%) (p = 1)  
- **Transducer/Linker (MCP/HAMP/CheW):** Rho = 3 (15.0%) vs Ctrl = 0 (0.0%) (p = 0.2308)  

---

## 2. Comprehensive Results Table

| Stratum | Subset | Marker | Rhodopsin Hits | Control Hits | Odds Ratio (95% CI) | Fisher p | McNemar p |
|---|---|---|---|---|---|---|---|
| All_Complete_EQ_Rhodopsins | all_loci | Any Sensory / Signaling | 6/20 (30.0%) | 3/20 (15.0%) | 2.429 (0.512-11.511) | 0.4506 | 0.4531 |
| All_Complete_EQ_Rhodopsins | all_loci | TCS Histidine Kinase | 2/20 (10.0%) | 0/20 (0.0%) | 5.541 (0.249-123.086) | 0.4872 | 0.5 |
| All_Complete_EQ_Rhodopsins | all_loci | TCS Response Regulator | 2/20 (10.0%) | 2/20 (10.0%) | 1.0 (0.127-7.893) | 1 | 1 |
| All_Complete_EQ_Rhodopsins | all_loci | Transducer / Linker (MCP/HAMP/CheW) | 3/20 (15.0%) | 0/20 (0.0%) | 8.2 (0.396-169.908) | 0.2308 | 0.25 |
| All_Complete_EQ_Rhodopsins | all_loci | Sensor Domain (PAS/GAF) | 1/20 (5.0%) | 1/20 (5.0%) | 1.0 (0.058-17.182) | 1 | 1 |
| All_Complete_EQ_Rhodopsins | all_loci | Second Messenger (GGDEF/EAL) | 0/20 (0.0%) | 1/20 (5.0%) | 0.317 (0.012-8.261) | 1 | 1 |
| All_Complete_EQ_Rhodopsins | full_flank_only | Any Sensory / Signaling | 0/7 (0.0%) | 0/7 (0.0%) | 1.0 (0.017-57.315) | 1 | 1 |
| All_Complete_EQ_Rhodopsins | full_flank_only | TCS Histidine Kinase | 0/7 (0.0%) | 0/7 (0.0%) | 1.0 (0.017-57.315) | 1 | 1 |
| All_Complete_EQ_Rhodopsins | full_flank_only | TCS Response Regulator | 0/7 (0.0%) | 0/7 (0.0%) | 1.0 (0.017-57.315) | 1 | 1 |
| All_Complete_EQ_Rhodopsins | full_flank_only | Transducer / Linker (MCP/HAMP/CheW) | 0/7 (0.0%) | 0/7 (0.0%) | 1.0 (0.017-57.315) | 1 | 1 |
| All_Complete_EQ_Rhodopsins | full_flank_only | Sensor Domain (PAS/GAF) | 0/7 (0.0%) | 0/7 (0.0%) | 1.0 (0.017-57.315) | 1 | 1 |
| All_Complete_EQ_Rhodopsins | full_flank_only | Second Messenger (GGDEF/EAL) | 0/7 (0.0%) | 0/7 (0.0%) | 1.0 (0.017-57.315) | 1 | 1 |
| All_Complete_EQ_Rhodopsins | edge_truncated_only | Any Sensory / Signaling | 6/13 (46.15%) | 3/13 (23.08%) | 2.857 (0.528-15.474) | 0.411 | 0.4531 |
| All_Complete_EQ_Rhodopsins | edge_truncated_only | TCS Histidine Kinase | 2/13 (15.38%) | 0/13 (0.0%) | 5.87 (0.255-135.154) | 0.48 | 0.5 |
| All_Complete_EQ_Rhodopsins | edge_truncated_only | TCS Response Regulator | 2/13 (15.38%) | 2/13 (15.38%) | 1.0 (0.119-8.421) | 1 | 1 |
| All_Complete_EQ_Rhodopsins | edge_truncated_only | Transducer / Linker (MCP/HAMP/CheW) | 3/13 (23.08%) | 0/13 (0.0%) | 9.0 (0.417-194.077) | 0.22 | 0.25 |
| All_Complete_EQ_Rhodopsins | edge_truncated_only | Sensor Domain (PAS/GAF) | 1/13 (7.69%) | 1/13 (7.69%) | 1.0 (0.056-17.904) | 1 | 1 |
| All_Complete_EQ_Rhodopsins | edge_truncated_only | Second Messenger (GGDEF/EAL) | 0/13 (0.0%) | 1/13 (7.69%) | 0.309 (0.011-8.3) | 1 | 1 |
| All_Complete_EQ_Rhodopsins | mag_dereplicated | Any Sensory / Signaling | 6/20 (30.0%) | 3/20 (15.0%) | 2.429 (0.512-11.511) | 0.4506 | 0.4531 |
| All_Complete_EQ_Rhodopsins | mag_dereplicated | TCS Histidine Kinase | 2/20 (10.0%) | 0/20 (0.0%) | 5.541 (0.249-123.086) | 0.4872 | 0.5 |
| All_Complete_EQ_Rhodopsins | mag_dereplicated | TCS Response Regulator | 2/20 (10.0%) | 2/20 (10.0%) | 1.0 (0.127-7.893) | 1 | 1 |
| All_Complete_EQ_Rhodopsins | mag_dereplicated | Transducer / Linker (MCP/HAMP/CheW) | 3/20 (15.0%) | 0/20 (0.0%) | 8.2 (0.396-169.908) | 0.2308 | 0.25 |
| All_Complete_EQ_Rhodopsins | mag_dereplicated | Sensor Domain (PAS/GAF) | 1/20 (5.0%) | 1/20 (5.0%) | 1.0 (0.058-17.182) | 1 | 1 |
| All_Complete_EQ_Rhodopsins | mag_dereplicated | Second Messenger (GGDEF/EAL) | 0/20 (0.0%) | 1/20 (5.0%) | 0.317 (0.012-8.261) | 1 | 1 |
| Type-1_Microbial_Rhodopsin | all_loci | Any Sensory / Signaling | 5/10 (50.0%) | 1/10 (10.0%) | 9.0 (0.809-100.143) | 0.1409 | 0.125 |
| Type-1_Microbial_Rhodopsin | all_loci | TCS Histidine Kinase | 1/10 (10.0%) | 0/10 (0.0%) | 3.316 (0.12-91.607) | 1 | 1 |
| Type-1_Microbial_Rhodopsin | all_loci | TCS Response Regulator | 2/10 (20.0%) | 1/10 (10.0%) | 2.25 (0.17-29.769) | 1 | 1 |
| Type-1_Microbial_Rhodopsin | all_loci | Transducer / Linker (MCP/HAMP/CheW) | 3/10 (30.0%) | 0/10 (0.0%) | 9.8 (0.438-219.26) | 0.2105 | 0.25 |
| Type-1_Microbial_Rhodopsin | all_loci | Sensor Domain (PAS/GAF) | 1/10 (10.0%) | 1/10 (10.0%) | 1.0 (0.054-18.575) | 1 | 1 |
| Type-1_Microbial_Rhodopsin | all_loci | Second Messenger (GGDEF/EAL) | 0/10 (0.0%) | 0/10 (0.0%) | 1.0 (0.018-55.271) | 1 | 1 |
| Type-1_Microbial_Rhodopsin | full_flank_only | Any Sensory / Signaling | 0/2 (0.0%) | 0/2 (0.0%) | 1.0 (0.014-73.27) | 1 | 1 |
| Type-1_Microbial_Rhodopsin | full_flank_only | TCS Histidine Kinase | 0/2 (0.0%) | 0/2 (0.0%) | 1.0 (0.014-73.27) | 1 | 1 |
| Type-1_Microbial_Rhodopsin | full_flank_only | TCS Response Regulator | 0/2 (0.0%) | 0/2 (0.0%) | 1.0 (0.014-73.27) | 1 | 1 |
| Type-1_Microbial_Rhodopsin | full_flank_only | Transducer / Linker (MCP/HAMP/CheW) | 0/2 (0.0%) | 0/2 (0.0%) | 1.0 (0.014-73.27) | 1 | 1 |
| Type-1_Microbial_Rhodopsin | full_flank_only | Sensor Domain (PAS/GAF) | 0/2 (0.0%) | 0/2 (0.0%) | 1.0 (0.014-73.27) | 1 | 1 |
| Type-1_Microbial_Rhodopsin | full_flank_only | Second Messenger (GGDEF/EAL) | 0/2 (0.0%) | 0/2 (0.0%) | 1.0 (0.014-73.27) | 1 | 1 |
| Type-1_Microbial_Rhodopsin | edge_truncated_only | Any Sensory / Signaling | 5/8 (62.5%) | 1/8 (12.5%) | 11.667 (0.922-147.57) | 0.1189 | 0.125 |
| Type-1_Microbial_Rhodopsin | edge_truncated_only | TCS Histidine Kinase | 1/8 (12.5%) | 0/8 (0.0%) | 3.4 (0.12-96.706) | 1 | 1 |
| Type-1_Microbial_Rhodopsin | edge_truncated_only | TCS Response Regulator | 2/8 (25.0%) | 1/8 (12.5%) | 2.333 (0.167-32.585) | 1 | 1 |
| Type-1_Microbial_Rhodopsin | edge_truncated_only | Transducer / Linker (MCP/HAMP/CheW) | 3/8 (37.5%) | 0/8 (0.0%) | 10.818 (0.463-252.804) | 0.2 | 0.25 |
| Type-1_Microbial_Rhodopsin | edge_truncated_only | Sensor Domain (PAS/GAF) | 1/8 (12.5%) | 1/8 (12.5%) | 1.0 (0.052-19.361) | 1 | 1 |
| Type-1_Microbial_Rhodopsin | edge_truncated_only | Second Messenger (GGDEF/EAL) | 0/8 (0.0%) | 0/8 (0.0%) | 1.0 (0.018-56.466) | 1 | 1 |
| Type-1_Microbial_Rhodopsin | mag_dereplicated | Any Sensory / Signaling | 5/10 (50.0%) | 1/10 (10.0%) | 9.0 (0.809-100.143) | 0.1409 | 0.125 |
| Type-1_Microbial_Rhodopsin | mag_dereplicated | TCS Histidine Kinase | 1/10 (10.0%) | 0/10 (0.0%) | 3.316 (0.12-91.607) | 1 | 1 |
| Type-1_Microbial_Rhodopsin | mag_dereplicated | TCS Response Regulator | 2/10 (20.0%) | 1/10 (10.0%) | 2.25 (0.17-29.769) | 1 | 1 |
| Type-1_Microbial_Rhodopsin | mag_dereplicated | Transducer / Linker (MCP/HAMP/CheW) | 3/10 (30.0%) | 0/10 (0.0%) | 9.8 (0.438-219.26) | 0.2105 | 0.25 |
| Type-1_Microbial_Rhodopsin | mag_dereplicated | Sensor Domain (PAS/GAF) | 1/10 (10.0%) | 1/10 (10.0%) | 1.0 (0.054-18.575) | 1 | 1 |
| Type-1_Microbial_Rhodopsin | mag_dereplicated | Second Messenger (GGDEF/EAL) | 0/10 (0.0%) | 0/10 (0.0%) | 1.0 (0.018-55.271) | 1 | 1 |
| Heliorhodopsin_Type-3 | all_loci | Any Sensory / Signaling | 1/10 (10.0%) | 2/10 (20.0%) | 0.444 (0.034-5.88) | 1 | 1 |
| Heliorhodopsin_Type-3 | all_loci | TCS Histidine Kinase | 1/10 (10.0%) | 0/10 (0.0%) | 3.316 (0.12-91.607) | 1 | 1 |
| Heliorhodopsin_Type-3 | all_loci | TCS Response Regulator | 0/10 (0.0%) | 1/10 (10.0%) | 0.302 (0.011-8.332) | 1 | 1 |
| Heliorhodopsin_Type-3 | all_loci | Transducer / Linker (MCP/HAMP/CheW) | 0/10 (0.0%) | 0/10 (0.0%) | 1.0 (0.018-55.271) | 1 | 1 |
| Heliorhodopsin_Type-3 | all_loci | Sensor Domain (PAS/GAF) | 0/10 (0.0%) | 0/10 (0.0%) | 1.0 (0.018-55.271) | 1 | 1 |
| Heliorhodopsin_Type-3 | all_loci | Second Messenger (GGDEF/EAL) | 0/10 (0.0%) | 1/10 (10.0%) | 0.302 (0.011-8.332) | 1 | 1 |
| Heliorhodopsin_Type-3 | full_flank_only | Any Sensory / Signaling | 0/5 (0.0%) | 0/5 (0.0%) | 1.0 (0.017-59.998) | 1 | 1 |
| Heliorhodopsin_Type-3 | full_flank_only | TCS Histidine Kinase | 0/5 (0.0%) | 0/5 (0.0%) | 1.0 (0.017-59.998) | 1 | 1 |
| Heliorhodopsin_Type-3 | full_flank_only | TCS Response Regulator | 0/5 (0.0%) | 0/5 (0.0%) | 1.0 (0.017-59.998) | 1 | 1 |
| Heliorhodopsin_Type-3 | full_flank_only | Transducer / Linker (MCP/HAMP/CheW) | 0/5 (0.0%) | 0/5 (0.0%) | 1.0 (0.017-59.998) | 1 | 1 |
| Heliorhodopsin_Type-3 | full_flank_only | Sensor Domain (PAS/GAF) | 0/5 (0.0%) | 0/5 (0.0%) | 1.0 (0.017-59.998) | 1 | 1 |
| Heliorhodopsin_Type-3 | full_flank_only | Second Messenger (GGDEF/EAL) | 0/5 (0.0%) | 0/5 (0.0%) | 1.0 (0.017-59.998) | 1 | 1 |
| Heliorhodopsin_Type-3 | edge_truncated_only | Any Sensory / Signaling | 1/5 (20.0%) | 2/5 (40.0%) | 0.375 (0.022-6.348) | 1 | 1 |
| Heliorhodopsin_Type-3 | edge_truncated_only | TCS Histidine Kinase | 1/5 (20.0%) | 0/5 (0.0%) | 3.667 (0.118-113.736) | 1 | 1 |
| Heliorhodopsin_Type-3 | edge_truncated_only | TCS Response Regulator | 0/5 (0.0%) | 1/5 (20.0%) | 0.273 (0.009-8.46) | 1 | 1 |
| Heliorhodopsin_Type-3 | edge_truncated_only | Transducer / Linker (MCP/HAMP/CheW) | 0/5 (0.0%) | 0/5 (0.0%) | 1.0 (0.017-59.998) | 1 | 1 |
| Heliorhodopsin_Type-3 | edge_truncated_only | Sensor Domain (PAS/GAF) | 0/5 (0.0%) | 0/5 (0.0%) | 1.0 (0.017-59.998) | 1 | 1 |
| Heliorhodopsin_Type-3 | edge_truncated_only | Second Messenger (GGDEF/EAL) | 0/5 (0.0%) | 1/5 (20.0%) | 0.273 (0.009-8.46) | 1 | 1 |
| Heliorhodopsin_Type-3 | mag_dereplicated | Any Sensory / Signaling | 1/10 (10.0%) | 2/10 (20.0%) | 0.444 (0.034-5.88) | 1 | 1 |
| Heliorhodopsin_Type-3 | mag_dereplicated | TCS Histidine Kinase | 1/10 (10.0%) | 0/10 (0.0%) | 3.316 (0.12-91.607) | 1 | 1 |
| Heliorhodopsin_Type-3 | mag_dereplicated | TCS Response Regulator | 0/10 (0.0%) | 1/10 (10.0%) | 0.302 (0.011-8.332) | 1 | 1 |
| Heliorhodopsin_Type-3 | mag_dereplicated | Transducer / Linker (MCP/HAMP/CheW) | 0/10 (0.0%) | 0/10 (0.0%) | 1.0 (0.018-55.271) | 1 | 1 |
| Heliorhodopsin_Type-3 | mag_dereplicated | Sensor Domain (PAS/GAF) | 0/10 (0.0%) | 0/10 (0.0%) | 1.0 (0.018-55.271) | 1 | 1 |
| Heliorhodopsin_Type-3 | mag_dereplicated | Second Messenger (GGDEF/EAL) | 0/10 (0.0%) | 1/10 (10.0%) | 0.302 (0.011-8.332) | 1 | 1 |

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