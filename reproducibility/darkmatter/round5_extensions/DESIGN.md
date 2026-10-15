# Exploratory extension design, fixed before computation

These post hoc hypotheses were selected after four rounds of review. Results will
be reported together, including negative outcomes; they are not preregistered or
confirmatory. Root is inspecting feasibility while Flash independently scouts
extensions. Argon will receive all outputs and assess their value and limitations.

H1: Does replacing Bakta with the original MAG-caller configuration recover CDSs
at DNA-supported culture loci? Hold the existing nucleotide alignments and eligible
MAG calls fixed. Call each of the 22 paired culture assemblies with Pyrodigal meta
mode, compare >=50% aligned-span overlap with Bakta, and assess same-orientation
first-base phase. Report dark/annotated arms, genome-pair differences, complete
calls >=100 aa, and membership in the 47 recurrent families. This tests caller
sensitivity, not coding truth; it conditions on nucleotide support and strain pairs.
A 90% span-overlap sensitivity will assess the loose primary overlap criterion.

H2: Does detection in all five expeditions exceed the expectation given family
prevalence and sampling structure? Permute expedition labels within site x soil
position, using the same permutation for every family to preserve family
co-occurrence. Hold the family x sample detection matrix fixed; every family's
within-stratum prevalence and every stratum's expedition sample counts remain
unchanged. Use 499 permutations, seed20261015; primary statistic is the proportion
of all106042admitted families detected in all5expeditions. Also report fixed
occupancy bins. Interpret results relative to label exchangeability; sequencing
batch/platform, mapping ambiguity and DNA persistence remain possible causes.

H3: Does PICRUSt carry sample-specific information beyond shared reference-wide
KO ranks? Rebuild the448paired profiles with the original duplicate-collapse rule.
Compute Spearman correlations over a fixed common KO universe (including zero ties)
for every MAG/PICRUSt pairing. Compare correctly paired mean correlation with4999
permutations of PICRUSt sample identities within expedition x soil position; a
sensitivity further restricts exchanges to site. Report numbers of exchangeable
samples and block sizes. Separately describe presence overlap at within-sample
KO-fraction thresholds0,1e-6,1e-5,1e-4,1e-3, applied to both normalised profiles.
MAG absence is not a negative truth set and these overlap metrics are not accuracy.
The common-universe rank estimand differs from the existing shared-positive-only rho.

All heavy work runs under Slurm on IBEX; input hashes, scripts, seeds, full compact
outputs and job receipts will be retained. No database searches or raw pipelines
are silently substituted. No neighbouring phage/defensome or life-on-air content
is incorporated.
