# Root response to the final review

All three requested Opus 5.5 reviews completed successfully. This response records root's own decisions and corrections after round 3; it does not claim a fourth external review or reviewer endorsement of the final edits.

| Finding | Disposition |
|---|---|
| R3-01 Culture interpretation | Strengthened the gene-call dependence result. Reported the union of protein/DNA support as a union, not a DNA-only rate. The review's suggested “identical DNA” was rejected: the search threshold is95%identity. Its “DNA present” denominators752/70,354 also include protein-supported calls failing the nucleotide threshold; those denominators cannot describe DNA-only detection. Extended the locus audit to every DNA-supported call and traced eligible dark calls through cluster90 to admitted families; results are recorded below. |
| R3-02 Flagged-family dN/dS | Added the flagged distributions, including AntiFam47evaluablefamilies,median0.8307,76.6%below1, and the91-familyunion,median0.7622,76.9%below1. These make the limitation of below-one dN/dS explicit. Rejected “necessary but insufficient”: a coding gene need not have dN/dS below one. Abstract describes the annotated-versus-shifted frame contrast directly. |
| R3-03 Overlap power | Stated the original Prodigal overlap design and observed195ntmaximum. The screen is conditional on both calls surviving gene selection and has low power for suppressed alternatives. Added primary Hyatt/Nelson/Yooseph sources. |
| R3-04 Persistence | Removed the Discussion-lead percentage; qualified supplementary detection rates as prevalence-dependent; corrected the Introduction's inaccurate matched-recurrence-control description. |
| R3-05 Rarefaction | Abstract names the combined desert catalogue and SPIRE. Discussion distinguishes GEM's4.27-point gain, SPIRE21.30 and Atacama−4.67. Figure caption gives Atacama/control residuals and identifies the analytical expectation. No real subset search is implied. |
| R3-06 Selection wording | Harmonised below0.1, named dark/control medians, identified93,220/95,541common-support families used for weighted means before equal-count sampling, cited Fig2b. |
| R3-07 Figures | Fixed colliding labels, renamed “Exact” to “Base-level,” cited missing panels, made Fig5c start at zero, replaced “small in every stratum” with explicit44vs50example. Map remains optional; no new map implied. |
| R3-08 Enrichment | **Rejected with exact source evidence.** `v2_dark_composition_by_phylum.tsv` gives0.198381/0.043582 and an explicitly stored enrichment4.551874, correctly rounding to4.6. The reviewer recomputed from rounded19.8%, causing the proposed4.5. Frozen exact source and added verification. |
| R3-09 Minor text/cohorts | Changed25to24shared species after contaminant exclusion; removed unsupported paralog-set explanation; trimmed supplementary-only Discussion; clarified548-sample stress and549-sample Thermoplasmatota snapshots and their44/45paired-site counts. They are distinct retained cohorts, not a numerical contradiction. |
| R3-10 Environment | Further removed directional language; explicitly stated that more site-aggregate sulfur associations can reflect between-site structure. Aggregation is not claimed to remove spatial confounding. |
| R3-11 Data/verification | Replaced TODOs with explicit unresolved availability statements; no claim that local revisions are already published. Separated culture diagnostic production from the checker; checker is read-only unless explicitly asked to write its fixed ledger. External raw-access, merge, chemistry and deposition gaps remain open. Historical verifier was run and its failures retained, not described as passing. |
| R3-12 Citations | Added verified overlapping-gene and gene-prediction primary sources. Checked Yeo's exact title against current records: retained bioRxiv designation; no verified journal version found. Additional free fulltexts retained and identity-audited. |

## Framing and remaining decision points

Nature Communications remains the selected stretch target. Root recommends a genome-resource paper with a main methodological finding: evolutionary, culture and structural interpretations of dark proteins depend on gene calling, reference retrieval and comparison populations. The author was asked whether to prefer that emphasis or lead almost entirely with the resource. No manuscript submission, public release, Git push or external message has occurred.

The substantive optional extensions in the review are documented, not treated as mandatory new rounds: harmonised gene calling, exhaustive alternative-frame homology, domain segmentation, dark fractions in close relatives, and actual random reference subsets. The existing relatives ANI/AAI analysis is preserved; it does not replace a relative-specific dark-fraction comparison.

Required external work before submission remains raw-file location/access, duplicate-run merge provenance, XRF metadata/non-detect semantics, ENA public release or reviewer access, MAG deposition with MIMAG fields, archived release DOI, and provision of the companion manuscripts to the editor. Historical registry rebinding remains explicitly separate from the passing compact revision checks.

## Completed final culture audit (IBEX job 53680886)

All DNA-supported calls were classified against Bakta, including the protein-shared group. The corrected denominators are735dark and67,802annotated calls. No major CDS overlap occurs in344/735dark calls(46.8%) and607/67,802annotated calls(0.90%). There is one such call in each protein-shared group, confirming that protein support cannot simply be assumed to imply the same aligned locus. The original unshared-only counts343/606remain correct and are preserved.

The cluster90-to-family join maps65of1,341eligible dark calls into47admitted families. The broader culture comparison is therefore explicitly distinguished from the recurrent-family analysis. Full per-call outputs and input/source provenance are retained under `culture_locus_all/` and `culture_family_membership/`.
