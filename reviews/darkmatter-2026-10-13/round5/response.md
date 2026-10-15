# Root adjudication of Argon's fifth review

The requested agy Gemini 4 Argon High review completed successfully in723.16seconds,
conversation fcda4839-2e2f-4792-ae8f-6031f94600bc. The exact CLI model and prompt,
unaltered report, model receipt and audit scripts are retained. Flash scouting
preceded three root-designed IBEX extensions. The report reviews the frozen round-4
manuscript plus the candidate extension results; root's final prose is subsequent.
Nature Communications remains the selected stretch target, not an acceptance forecast.

## Reliability

Argon's executable audits contain useful results, but its prose again includes
fabricated quotations and substantial numerical/method errors. The opening's62sites,
110isolate genomes, five-criterion/multispecies admission gate and expression-based
validation are not the manuscript's design. Current sampling is59core sites plus
two auxiliary locations; culture source has722assemblies; admission requires at
least three distinct90%-cluster representatives across three samples and two sites.
No136expressed-target cohort or metatranscriptomic validation is established here.
Several section titles and supposedly exact quotations are absent from the frozen
source. Recommendations were checked against sources rather than accepted wholesale.

## Main comments

1. **Structural controls: partially accepted.** Added the independently recomputed
   post-gate confidence imbalance (mean pLDDT0.775dark,0.867controls). The old
   standardized-floor table exists, but its all-model failure ratios combine retrieval
   failures with score failures. It also uses PDB-derived positive controls whose
   chains overlap searched references. Confidence/length weighting cannot repair
   that ascertainment or make them exchangeable remote-homology controls. We retain
   the already explicit retrieval decomposition and reject restoring the1.38ratio
   as an inferential novelty comparison or treating an interval spanning1 as
   equivalence. The review's length bins and within-cell-bootstrap description
   differ from05e_stratified_floor.py, which uses100–149,150–249,250–399,>=400aa
   and resamples whole arms. Domain-resolved, homology-excluded calibration is a
   genuine further test in EXTENSIONS.md; the old reweighting is not new work.

2. **Pfam comparison: accepted in bounded form.** Independently aggregated
   per_protein_pfam.tsv. Added one-representative comparisons93/1,031(9.02%) below
   qTM0.5 versus101/1,807(5.59%) above it, avoiding member-versus-representative
   replication mismatch. The124/427DUF-only observation is present in the source
   but is not proof of known function for the remaining named domains. Pfam absence,
   long uncovered intervals and low whole-query TM do not prove domain architecture
   or structural novelty. Those causal claims are rejected.

3. **Length and temporal detection: accepted, with an additional experiment.**
   Added the retained representative-length-stratified detection percentages and
   the read-length/breadth-threshold qualification. The original new null reports
   all551samples, not550; it has84,174all-five families, not84,180. Its10–24occupancy
   bin is31.15%observed versus41.43%null, not2.51%versus12.78%; the25–49bin is68.73%
   versus74.97%. Expedition2 was July2023, not March2022. Counts of always/never
   detected in499permutations do not prove exact null probabilities0or1.
   The review correctly raises within-visit replication. Root therefore designed
   and ran a second null that permutes whole visits only among equal-size visits
   within site and position. Job53698522 completed in2m1s. Its result is reported
   alongside the original, not substituted selectively: nullmean80.54%, central95%range
   78.95–82.11%, lower-tail p0.11. The observed79.38%lies within that range. Visit/batch exchangeability
   remains an assumption; the review's attribution of the deficit to replicate
   clustering or spatial occupancy alone is not established.

4. **Caller substitution: accepted with correct populations and transitions.**
   Integrated H1, complete>=100-aa and recurrent-family subsets, pair-level direction
   and concentration (44/118gains in one pair), and117complete calls among118gains.
   Actual nucleotide thresholds are>=90%coverage and>=95%identity, not80%/70%.
   Overlap is assessed before strand/phase classification; no-CDS does not mean
   no same-strand CDS. Annotated baseline is607/67,802, not592; complete-long dark
   baseline is25/59, not19/59. Of118gains,116have compatible first-base phase, not
   all118. The review's proposed explanation of every residual as strain divergence
   or dynamic programming is unsupported. Shared caller settings are a controlled
   substitution at fixed loci, not independent coding validation or a full reciprocal
   reannotation experiment. We report descriptively without a gene-level test.

5. **Contamination: partially accepted.** Named both SNTC read-profile representatives
   (E.coli and C.acnes) and identified the separate vector-contig screen, including
   Mixta alhagi. Root independently reproduced12unique flagged representatives,
   42,831Pyrodigal instances and1,792dark instances; these are explicitly different
   denominators from Bakta annotation counts and genome-linked clusters. Removing
   the12changes median per-genome darkness6.31%->6.34%. Counts alone do not bound
   abundance-weighted ecological effects, so the remaining sensitivity requirement
   is retained rather than falsely closed. Vector sequence in a bin does not prove
   the entire bin is contamination. The review's quoted Methods/Discussion passages
   are not exact passages in the baseline.

6. **Cross-references and culture counts: accepted.** Both dangling references are
   real: Supplementary TableS3 and line-wrapped Extended Data Table1. Round4's root
   response incorrectly denied them; this is corrected in notes/CORRECTIONS.md.
   Added the exact13-KO solute panel and logical rules to Methods and removed the
   unsupported generic Extended Data sentence. Discussion now uses22retained
   representative–culture groups, consistent with Results.

## Minor comments and completed-extension placement

- No8,857shared-KO count appears in the baseline text. The new estimand is explicitly
  the8,860identifier intersection including zero ties; no invented count discrepancy
  was patched. The review's14blocks is15; its paired median0.8442 differs from its
  audit0.8510. Its grand0.8215matrix mean includes diagonals, despite being described
  as off-diagonal. The review's housekeeping-copy-number mechanism is an inference,
  not a result. H3 now has one concise main-text sentence and the full paired-null
  and threshold results in Supplementary NoteS2, including the distinct146-sample
  within-visit comparison. The degenerate1e-3row is disclosed with446undefined recall
  denominators, not silently omitted after inspecting results.
- Completed seven bibliography entries from retained DOI-matched primary XML.
  Argon's tRNAscan pages, OLGenie end page and Orr journal are wrong; verified values
  are9077–9096,2440–2449 and Bioinformatics Advances6(1):vbag160. Yooseph is5(3):e16.
  Existing valid BibTeX `and others` renders as et al.; no fabricated author expansion.
- The partial-flag parsing concern is valid for downstream readers. Root's summary
  uses original CSV strings and its independent audit explicitly sets string dtype.
  Argon's first audit parsed them numerically and incorrectly reported0%complete;
  its second audit corrected that to117/118. Its claim that culture_caller.py line161
  calls zfill is invented. Added explicit reading guidance; original data preserved.
- PDF alerts are the legitimate heatmap labels69/90/61, with no instruction injection.
  The review's figure/page attributions are inconsistent; root inspected the renders.

## Further hypotheses

EXTENSIONS.md ranks four corrected directions. We do not accept the review's claim
that extra SW-10-69-26 genes prove accessory expansion or that null AAI correlation
refutes reference bias; ancestry, quality, gene-calling basis and ratio dependence
need explicit tests. SMAG40.47%is absence, not hits; its hit rate is59.53%versus
GEM35.62%, not4.93%. The proposed8.2-fold contrast is false. The descriptive union
absence28.74%->26.83%is real, but fair matched-reference tests remain to be done.
The invented IBEXproject/c2188 path and136expressed targets are excluded from the
neighbourhood proposal. No new phage, amplicon or life-on-air analysis is claimed.

## Validation scope

The final compact checks verify retained evidence and manuscript-facing values;
remote jobs retain input hashes, commands, seeds and outputs on IBEX. The original
large-file audit and77DOI-verified fulltexts remain preserved. No new literature
campaign was needed to support these control analyses. The historical verifier
still has31unreconciled manuscript gates, while62COMPUTE and72FROZENchecks pass.
Raw-file location/permission, merge provenance, XRF metadata, ENA/reviewer access,
MAG/code-data deposition and historical registry migration remain submission items.
No public release, submission, push or external message was made.
