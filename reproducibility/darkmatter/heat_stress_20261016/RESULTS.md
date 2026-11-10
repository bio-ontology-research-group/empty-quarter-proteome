# Heat-stress hypotheses: results and paper options

Analysis scope: 5,438 representative soil genomes, eight protein-maintenance families, 550 usable metagenomic profiles, and five visits. Results refer to encoded genomic capacity, not stress-induced expression. The February visit is the sole winter visit. Root executed and audited the analysis after Argon was unavailable; no Argon endorsement is claimed.

## HS1 — Do hot-season communities encode more protein-maintenance capacity?

**Not supported by these comparisons.** None of the primary hot-visit family comparisons passes BH q<0.05; minimum q=0.833. Minimum q values also remain above 0.05 under strict annotation (0.723), completeness adjustment (1.000), and copies-per-Mbp normalisation (0.625).

The August/February comparison contains 9 surface, 13 deep and 32 root-adjacent paired sites. July/February has only 3, 5 and 6 respectively. The July surface comparisons are descriptive. Across the 32 root-adjacent August/February pairs, the mean HSP20 difference is +0.0625 abundance-weighted copies per represented genome (signed-rank p=0.949; q=1). The corresponding HSP70 difference is +0.0262 (p=0.678; q=1). Within-Thermoplasmatota comparisons also do not pass correction (minimum q≈0.888).

This does not demonstrate absence of heat adaptation: transcriptional regulation, constitutive protection, untested pathways and fine-scale strain differences remain outside this test. It does rule out presenting these analyses as evidence for a general summer enrichment of the eight tested families.

## HS4 — Does depth or measured temperature distinguish a heat-associated response?

**Not supported by the available paired contrasts and environmental models.** Only five sites have both depths in August and February; the minimum adjusted q for the depth interaction is 0.5. July has only three fully paired sites.

The temperature/humidity models use 391 site/visit/compartment cells across 59 sites; none of the 32 declared coefficients passes correction (minimum q≈0.0756). The strongest temperature coefficients are negative, including HSP90, rather than a general positive heat-response pattern. Rainfall sensitivities also have no corrected association (minimum q≈0.0696). Temperature-only models retaining July measurements remain nonsignificant after correction (minimum q≈0.196).

Environmental measurements are collection-time conditions. The limited temporal replication and covariance of temperature, moisture, collection hour and visit prevent a heat-specific causal interpretation. The additional models do not justify calling any candidate a thermal-adaptation mechanism.

See `figures/HS1_HS4_effects.pdf`, its source tables and caption; confidence intervals in the figure are pointwise rather than multiplicity-adjusted.

## HS2 — Do unknown archaeal families recur beside protein-maintenance genes?

**Supported as a genomic-context observation; not established as heat-stress function or seasonal adaptation.** After replacing the invalid identifier-only join with exact sequence links, 2,630 recurrent unknown families were tested across 237 Thermoplasmatota representatives. Five pass the exact within-contig rotation null. Three of those five also pass a gene-length-matched comparison:

| Family identifier | Near an anchor / eligible genomes | Anchor | Exact-rotation q | Length-matched q |
|---|---:|---|---:|---:|
| F41Dr3::F41Dr3_256695_46 | 16 / 18 | GrpE | 3.47e-9 | 4.90e-6 |
| 49Dr1::49Dr1_177214_2 | 14 / 18 | DnaJ or GrpE | 1.19e-5 | 5.48e-4 |
| 24Dr1::24Dr1_1575644_3 | 5 / 29 | HSP20 | 0.00115 | 0.0464 |

The other two primary discoveries fail the length-matched threshold. One additional family passes only the length-matched sensitivity, not the primary rotation null; it is not promoted to a concordant finding. All tested results are retained.

No existing full-Pfam hit was found among the loci belonging to the five primary candidate families. This observation holds after the sequence correction; the incompatible annotations in the initial result were consequences of incorrect locus links. This does not prove absence of all homology or establish a new biochemical function.

All five candidates are concentrated within an archaeal order. None passes the sparse sensitivity selecting one representative per genus-label stratum. In the strongest family, 13 of 18 genomes lack an assigned genus; the remaining five belong to two named genus labels. The four strata in that sensitivity therefore must not be described as four identified genera. Shared ancestry can preserve neighbourhoods, and the data do not establish repeated independent evolution of the association.

Candidate-carrier abundances show no corrected hot-visit association (minimum exploratory q≈0.933). The supported observation is recurrent genomic proximity, strongest for the GrpE- and DnaJ-associated families. A role in protein maintenance remains a hypothesis.

## HS3 — Do persistent populations change gene dosage or alleles between visits?

**No convincing heat-associated change in the completed bacterial comparison; the focal archaeal pilot is undercovered.** The highest shared-abundance Thermoplasmatota pilot fails the prespecified coverage/breadth gate. Its outcome is inconclusive, not a negative result for archaea generally.

The Prauserella_A representative at site 14 in root-adjacent soil passed the pilot and was analysed using all reads from both February and August. Median background-gene depth is 259.5× and 319.25× respectively, with median breadth 1.0. All 14 protein-maintenance loci pass the full coverage gate; 56 of 60 background loci do so.

No maintenance-gene position among 12,873 jointly callable positions shows an allele-frequency total-variation difference of at least 0.20. Three such positions occur among 57,799 background positions. Gene-depth ratios provide no general summer increase: the median marker log2(August/February), after background-depth normalisation, is −0.246. The HSP20 ratio is approximately 0.776. These ratios describe relative read coverage; they do not establish physical copy-number loss or selection. Background loci also show substantial coverage-ratio variation, and local mapping, strain composition and genome-replication effects remain alternatives.

This is one matched site for one bacterial representative. It is technically well covered, but it cannot establish a replicated seasonal mechanism. Other matched sites for the same representative have much lower abundance; the next-highest minimum abundance is 0.116%, versus 5.846% at the analysed site. The pilot does not justify an unrestricted claim about all species.

## Root judgement and incorporation options

1. **Recommended: a short functional-context addition centred on the two strongest HS2 families.** Present their repeated association with GrpE/DnaJ as a concrete biological observation, with neighbourhood diagrams and full evidence in the supplement. Keep seasonal tests as supporting analyses. Do not call these proteins new chaperones or claim a heat mechanism. This would modestly strengthen the paper's biological interpretation of unknown Thermoplasmatota proteins.
2. **Conservative option: keep the entire heat-stress follow-up in the supplement/reproducibility package.** This preserves a clear main narrative while documenting a tested explanation and identifying candidates for functional work.
3. **A larger heat-adaptation section is not justified by these results.** That would require evidence connecting the candidate proteins to function or repeated heat-related activity, beyond adjacency and one deeply sequenced population pair.

The new work supplies useful candidates, not a demonstrated thermal-adaptation mechanism or a major increase in the case for Nature Communications. Robert approved the supplementary treatment with one or two main Results sentences. Supplementary Note S5 now reports HS1–HS4, the candidate table and methods; Supplementary Figure S1 shows seasonal and temperature effects. Two main Results sentences describe the two strongest genomic-context findings. Argon has not reviewed these results; this judgement is root's.
