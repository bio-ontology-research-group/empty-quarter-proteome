# Exploratory extension results for round 5

All three designs were recorded in DESIGN.md before these runs. These post hoc
analyses use the existing data; results are conditional on the original catalogue,
gene-calling, mapping and reference choices. Code and full compact tables accompany
this summary. IBEX jobs 53697241 (profiles, 1m57s) and 53697281 (caller, 1m36s)
completed with exit 0. Job 53697242 was cancelled while pending to replace the
approximate nucleotide-length filter with the original exact amino-acid length.

## H1: controlled caller substitution

The same 22 culture assemblies were called with Pyrodigal 3.7.1 meta mode, matching
the original MAG-caller configuration. Existing >=90%-coverage, >=95%-identity DNA
alignments were held fixed. At >=50% aligned-target-span CDS overlap, dark calls
without a culture CDS decreased from 344/735 (46.8%) with Bakta to 235/735 (32.0%)
with Pyrodigal. Of the original 344, 118 gained an overlapping CDS (116 compatible
with orientation and phase at the first aligned base); nine other calls lost overlap.
Annotated calls changed from 607/67,802 to 498/67,802; 182 gained overlap and 73 lost it.
This demonstrates caller dependence at fixed DNA-supported loci. Phase here checks
only the first aligned base, not full downstream codon correspondence; neither
caller provides independent evidence of translation or function.

For complete dark calls >=100 aa, no-CDS counts changed from 25/59 to 17/59 (nine
gains, one loss); for annotated calls, 220/53,969 to 198/53,969. The DNA-supported
admitted-family subset has only 42 dark calls from 10 pairs: no-CDS totals remain
18/42 with each caller (one gain and one loss). At 90% span overlap the overall dark
counts are 400/735 versus 281/735, annotated 1,700/67,802 versus 1,397/67,802.
The family subset again remains 18/42. Thus the aggregate improvement cannot be
transferred to the recurrent-family set. Genome-pair summaries are retained;
genes are not treated as independent biological replicates for a significance test.

## H2: five-expedition detection conditioned on sampling

Among all 106,042 admitted families and 551 mapping samples, observed all-five
frequency is 0.7937798. Shuffling expedition labels within site x soil position
(499 global permutations, seed 20261015) yields a mean 0.8273839 and central 95%
null range 0.8111640–0.8419018: observed minus null = -3.3604 percentage points.
The lower-tail permutation p is 0.002, upper-tail p 1.0 (exploratory, unadjusted).
There are 159 strata, 121 with multiple expedition labels, covering 509 samples.
All family sample occupancies and co-occurrence patterns, and within-stratum
expedition sample counts, are fixed. This does not support excess five-expedition
detection relative to this exchangeability model. It does not identify turnover,
activity or its cause: trip-associated sequencing/mapping differences and persistent
DNA remain alternatives. The six fixed occupancy-bin summaries are retained.

## H3: paired functional-profile null and threshold sensitivity

448 paired samples were rebuilt using the original amplicon duplicate-collapse
rule. Over the fixed common universe of 8,860 KOs, with zero ties included, mean
within-sample rank correlation is 0.8430955. Permuting PICRUSt sample labels within
expedition x soil position (4,999 permutations, seed 20261016; 15 blocks, all448
exchangeable) gives mean 0.8245498 (central95% null range0.8235045–0.8256062), an
advantage0.0185457, upper-tail p0.0002. With site also held fixed, only146samples
in66blocks are exchangeable: mean paired0.8430753 versus null0.8407527, advantage
0.0023226 (p0.0002). Each null uses its own exchangeable-sample cohort. Small
site-conditioned effects limit the interpretation despite low permutation tails.
This estimand differs from the manuscript's median shared-positive-only rho0.839.

Applying the fixed within-sample KO-fraction grid to both profiles gives median
Jaccard0.63386 at>0,0.69710 at>1e-6,0.69609 at>1e-5,0.68315 at>1e-4, and0 at>1e-3.
At the final threshold median MAG KO count is0, so it is a degenerate sensitivity,
not evidence of an informative threshold. MAG recall at that threshold is undefined
for zero-denominator samples; pandas' summary median excludes those NaNs. The full
per-sample table retains them. The two profiles' distinct constructions mean an
equal relative-abundance cutoff need not have equal biological detection meaning.
Predictor-only KOs cannot be classified as false positives from incomplete MAGs.
