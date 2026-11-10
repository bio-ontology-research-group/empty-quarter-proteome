# Reproducibility and interpretation audit

This is root's audit of root-executed analyses. It is not independent Argon review. Argon was unavailable because of its recorded quota failure; Robert subsequently authorised direct execution.

## Accepted ecology analysis

The 5,438 representative-genome records and all eight-family HMM counts were joined by unique genome identifier. The newly parsed HMM locus counts agree with the existing uniformly generated per-family count matrix. Strict counts require complete ORFs and HMM coverage of at least 0.70. The analysis uses 550 nonempty sample profiles, averaged to 414 site/visit/compartment cells. The all-unmapped profile F5Dr1 is excluded explicitly.

Paired comparisons use sites as units. Forty of the 48 planned primary hot-visit comparisons have enough paired sites for inference; the eight July surface comparisons have only three sites and remain descriptive. BH corrections were independently recomputed with SciPy's `false_discovery_control`, rather than the implementation used by the analysis. Paired differences and signed-rank p values were recomputed from the retained site-level table. The between/within/unshared-lineage decomposition sums to the original effect within numerical tolerance.

Encoded copy counts are abundance-weighted within the represented genome catalogue. They do not measure transcription, protein abundance, chaperone activity, or unrepresented organisms. A nonsignificant result is not an equivalence test or proof that heat does not matter.

## Environmental data

The July field sheet has nine headers but ten values per nonblank row, with an empty pressure column. The parser explicitly restores the missing header and asserts row widths. October's trailing blank rows are ignored. Humidity 194 at site 40 is preserved in the raw audit but set missing for analysis; no guessed decimal correction is used. These fixes invalidate the initial environmental model, retained under `ecology_initial` for provenance.

Temperature/humidity models use 391 complete cells across 59 sites. Temperature-only models retain 410 cells; rainfall models use the pre-existing environmental cohort and explicitly report their smaller denominators. Models use site-clustered uncertainty, site and compartment indicators, collection hour, and separate versions with visit indicators. Temperature at collection is not continuous soil-temperature exposure. A single February visit does not provide replicated winter seasons or separate all visit/year effects. Site 52 has a coordinate discrepancy and is excluded in sensitivity analyses. No temperature or moisture coefficient passes its declared multiplicity threshold.

## Input integrity and compute

The four complete-read bacterial FASTQs and the complete competitive mapping index were independently rehashed. All five SHA256 values match prior provenance. Input HMM files, abundance profiles, field sheets and supplementary annotations have source manifests. All resource-intensive work runs in Slurm workers; no login-node analysis was used. Full mapping requests 128 GB after the 96 GB pilot approached its memory allocation. Bulk reads, indexes, alignments and matrices remain on IBEX.

A public script-and-manifest package is not equivalent to public deposition of raw data. Re-execution currently requires institutional access to the listed inputs and environment. No external review is claimed.

## HS2 identifier failure and correction

The initial catalogue link used original-looking gene identifiers. Direct amino-acid comparison found only 100 of 136candidate loci matched their purported original sequences. Those results, including their statistical findings, are rejected and retained under `hs2_identifier_only_REJECTED`. The replacement uses exact amino-acid identity on the same contig and a unique original locus, allowing changed gene numbers. Unmatched and ambiguous sequences are excluded rather than assigned by approximate identifiers. All family tests must be rerun after this correction.

## Accepted HS2 rerun

The exact join retains 656,648 unique loci, with 120,069 unmatched and 460 ambiguous; 86,642 links use changed identifiers. It maps 83,115 loci to admitted unknown families. The accepted test universe contains 2,630 families and 28,253 eligible family/genome loci. The exact marginal null is a Poisson-binomial convolution, because each family contributes one locus per genome. The small Monte Carlo p-value floor is not used as a discovery limit. Primary and sensitivity tables remain separate; only three candidates pass both the exact rotation and length-matched tests. One further candidate passes only the length-matched sensitivity. No candidate is claimed to represent an independently evolved association across genera or a demonstrated biochemical function.

The statistical audit passes 2,925 table/identity checks. These include one-locus-per-genome checks and independent numerical recomputations, not independent biological validation. Existing Pfam results were audited for the primary candidate loci; zero hits remain after correcting sequence links. The original identifier-based findings and their annotation conflicts are rejected, not silently overwritten.

## Accepted HS3 comparison

All 14 marker loci pass the full coverage criterion in both samples. Four background loci fail the joint criterion; results retain their coverage flags. Normalisation uses the median across 60 background CDSs. Restricting to the 56 jointly eligible background loci changes the background normalisers from 261.13/319.43 to 263.41/321.13 mean-depth units (February/August), less than 1% in either sample. The large negative or positive gene-depth ratios are not interpreted as integer copy-number changes. Base qualities, alignment identity, unique competitive mapping and overlapping mates are handled explicitly in `hs3_map.py`. Allele comparisons require at least 20 observations in both samples. No read-level p value is presented as a biological replicate test.

The archaeal pilot is a feasibility failure, not evidence against adaptation. The bacterial pair is a deeply covered single-site observation, not a replicated seasonal experiment. No manuscript changes were made.

## Final acceptance checks

The final audit passes all 14 additional checks, including exact sequence-link use, all 14 marker coverage gates, allele-count identities, and single-part competitive reference mapping. The mapping logs account for exactly 189,755,536 February input reads and 207,238,870 August input reads, matching prior source records. Four background loci fail the joint coverage gate; zero-depth ratios there are not interpreted as biological copy-number changes. The artifact manifest retains checksums for remote-only outputs as well as compact evidence.
