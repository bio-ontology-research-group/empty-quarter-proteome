# Round 1: author assessment and response

13 October 2026. Reviewer: Claude `claude-opus-5-5`, auto mode with permission prompts disabled as requested. See `prompt.txt`, `opus_result.json` and `review.txt`. This is the first of three planned rounds. Rounds 2 and 3 have not started. The manuscript revision is provisional and is not submission-ready.

## Confirmed revisions and outstanding analyses

| Finding | Disposition in this revision | Remaining work for the next round |
| --- | --- | --- |
| R1-01 Coding-evidence denominator | Abstract and Results restrict selection evidence to recurrent families; report the admission funnel and 12.5% of dark representatives. | Test coding evidence in singleton/doubleton classes or maintain the restricted claim. |
| R1-02 Length and partial-gene cohorts | Correct genome-resolved dark values to 65.8% short and 56.4% partial; distinguish assembly-wide comparisons. | Derive the matching genome-resolved annotated distribution and render the figure. |
| R1-03 DNA recurrence is not translation evidence | Remove universal coding/accessory conclusions; distinguish recurrence from frame-specific constraint. | Other-frame homology, AntiFam and overlapping-ORF checks; evaluate tractable additional coding controls. |
| R1-04 Divergence-confounded dN/dS controls | Independently reproduce divergence strata and remove the claim that the control is conservative. | Match or stratify jointly on divergence, length, GC and lineage; document residual imbalance. |
| R1-05 Shifted-frame null and attrition | Correct input/evaluable/paired counts and rounding. Add independently audited stop counts. | Trace family-level attrition and rerun sensitivity where justified. The widespread internal-stop inference is not supported; see below. |
| R1-06 Cultured comparison | Exclude confirmed E. coli contaminant; recompute 22-pair rates and sharing; correct actual reciprocal-hit methods; describe protein-call sharing rather than proven accessory content. | Harmonise gene callers and test unmatched genes against partner nucleotide assemblies. Culture-derived assemblies are not all closed isolates. |
| R1-07 Taxonomic novelty definitions | Separate placeholder rank from first empty rank. | Ensure figure and registry use these definitions consistently. |
| R1-08 Thermoplasmatota contrast | Add independently reproduced same-placeholder-rank contrast (2.48-fold); delimit the family-share denominator; remove causal interpretation of matching. | Search close GTDB/MGnify relatives and quantify reference proximity. |
| R1-09 Structural-control calibration | Describe positive-control retrieval rather than a universal error floor; qualify whole-query versus local metrics. | Decompose no-row/threshold failures, inspect near-threshold scores, and assess comparable control search difficulty. |
| R1-10 Architectural novelty | Demote architecture claims to hypotheses. | Domain segmentation and domain-aware comparison if architecture remains central. |
| R1-11 Persistence/mapping | Coding claims no longer rest on persistence alone. | Mapping-specificity sensitivity and site/expedition-aware null remain open. |
| R1-12 Sequencing/deposition | Reconcile 554 runs to 551 samples and private successful ENA receipts; correct NovaSeq X/6000 methods. State that public access was unavailable at audit time. | Resolve duplicate-run merge provenance, permissions and missing raw folder; public release remains outstanding. |
| R1-13 Reference coverage/rarefaction | Narrow desert-specific interpretation; describe full-catalogue endpoint agreement as consistency. | Actual subsampling and comparison with non-desert soil where required by the retained claim. |
| R1-14 Environmental analyses | Recorded as unresolved; no new causal environmental claim introduced. | Establish XRF zero semantics from laboratory records, model repeated sites/batches and compositional effects, document climate provenance. |
| R1-15 Reproducibility/figures | Add current compact inputs, independent arithmetic checks, sensitivity code and IBEX inventories; revised LaTeX builds. | Replace stale manuscript-era gates with evidence-based current checks; build missing figures and current supplement; freeze environments and database versions. |
| R1-16 Contaminants/controls | Confirmed E. coli exclusion applied to culture comparison; single-control limitations stated. | Catalogue-wide contaminant sensitivity and control-mapping sensitivity remain open. |
| R1-17 Phylogenetic dependence | Lineage associations described without causal attribution. | Evaluate phylogenetically informed modelling if quantitative lineage explanation remains central. |
| R1-18 Metabolic reconstruction | Replace inferred physiology/auxotrophy with provisional annotation patterns; qualify missing markers and pathway completeness. | Audit module rules and lineage-sensitive marker searches before strengthening these claims. |
| R1-19 Literature | Update Aplakidou to verified NAR 2026 DOI; correct interpretation of Maeke; preserve source-identity audit. | Verify outstanding full texts, add targeted references as analyses and journal framing settle. |

## Reviewer statements requiring correction or qualification

* An IBEX audit (job 53676424) of all 10,000 frame-test input families found 50,898 sequences, 43,953 terminal stops and only 21 sequences with internal stops (16 families). Terminal stops dominate. This does not establish an extraction failure; it does not validate the entire null either. Code and hashed input list are in `reproducibility/darkmatter/frame_stops/` and `audit_frame_stops.py`.
* The Fesenko download titled *The Cryptic Bacterial Microproteome* is a corresponding preprint, not an unrelated paper. Record version differences when citing the final article.
* Thermoplasmatota's 19.8% share refers to the 77,069 single-phylum-attributable families. Its share of all 106,042 families is approximately 14.4%, not the review's approximate 14.6%.
* XRF zeros may encode non-detection, but laboratory metadata must establish their meaning before treating them as censored measurements.
* Failure to show excess no-hits with one local structural metric does not establish absence of local novelty. The revision describes metric-dependent evidence.

## Evidence and verification

`reproducibility/darkmatter/check_current_tables.py` verifies compressed-input checksums and recomputes selected current quantities. `round1_sensitivity.py` reproduces the contaminant exclusion, divergence strata and placeholder-rank contrast. These pass; they are table-level checks, not complete upstream reruns. The public manuscript builds to 26 pages with resolved bibliography. Both manuscript source copies and bibliography copies agree. The older verification harness had 30 baseline gate failures, chiefly obsolete manuscript assertions; those remain a separate migration task. Figures 2 and 3 remain placeholders and Figure 1 needs updating.

Flash supplied five literature directions. The current manifest has 67 DOI leads: 45 DOI-verified retained full texts, 20 metadata-only records and two unresolved original leads. Incorrectly identified downloads are retained in quarantine, excluded from evidence. Search output is not itself a verified source. Full-text acquisition remains incomplete.

The first IBEX sweep contains 447,872 records (446,361 small-file hashes, 1,497 large-file metadata records and 14 broken-link errors). The second sweep adds 443,134 upstream metadata-only records; these inventories overlap in scope and must not be summed as unique data files. Raw reads have not been fully checksummed. 302 raw-file paths are inaccessible and 78 paths (39 runs) are missing. No claim of full reproducibility is made.

## Journal decision before round 2

The working recommendation is Nature Communications as an ambitious first target, conditional on resolving the central computational controls and reproducibility gaps. Microbiome and The ISME Journal are credible alternatives depending on the retained biological contribution. The current evidence does not justify choosing Nature Microbiology solely for impact; a substantially stronger biological validation would change that assessment.

Robert chooses the primary journal at this agreed stopping point. The next review packet must contain the original review, this response, the manuscript/code diff, actual rerun outputs, source audit and unresolved items. Round 2 tests the revised argument and new controls; round 3 is a final adversarial evidence/reproducibility check followed by author corrections. Unresolved critical findings remain explicit even after the third round.
