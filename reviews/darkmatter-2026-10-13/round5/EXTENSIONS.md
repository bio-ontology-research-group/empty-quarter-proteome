# Further testable hypotheses after Argon's fifth review

Root's ranked plan, after checking Argon's recommendations against the actual files.
These are exploratory follow-ups, not established mechanisms or prerequisites for
calling the three completed tests complete. Nature Communications remains the
stretch target. The best route to a stronger biological contribution is the first
question, followed by targeted gene context; merely adding more catalogue counts
is less likely to sharpen the paper.

## 1. What explains Thermoplasmatota darkness within the lineage?

**Question:** Does the dark fraction track external sequence divergence, shared-gene
coverage, or lineage-specific gene content after accounting for genome quality and
ancestry? Compare SW-10-69-26 and E2 without treating two classes as independent
replicates of an evolutionary mechanism.

**Available inputs:** the verified 237-row `v2_thermo_relatives_pairs.tsv` (retained
under round5/reviewer_inputs/06_darkmatter/lineage_darkness_clean/results/), and
`reproducibility/darkmatter/inputs/06_darkmatter/lineage_darkness_clean/results/v2_reps.tsv.gz`.
Full RBH outputs and identically single-genome-called comparator proteins are on
IBEX under `06_darkmatter/lineage_darkness_clean/`; the preserved
`v2_thermo_relatives_aai.py` documents their derivation. Family-presence and original
protein-call mappings would be needed to distinguish accessory content directly.

**Test:** use species representatives, with genus/family blocks and the retained
phylogeny where support permits; estimate within-lineage darkness-versus-AAI and
shared-gene-coverage associations controlling completeness, contamination, N50,
length and GC. Separate dark gene counts from ratios. Reconcile the gene-calling
basis of RBH coverage and the dark denominator before joining them. Compare
orthologous gene-family retention explicitly rather than labelling extra CDSs
"accessory expansion". Use matched callable complete/long proteins as sensitivity.

**Interpretation:** Argon's exploratory external-AAI correlation is near zero, a
useful lead, not proof that reference distance is irrelevant. A robust association
would quantify a contribution of representation; a null would restrict that simple
explanation within the sampled range. Larger genome/CDS counts alone cannot prove
evolutionary expansion, equivalent core divergence, or direction of gene gain.

**Feasibility:** existing compact data support an initial quality/ancestry audit;
family-content work requires an IBEX CPU analysis. No new wet lab. Highest priority.

## 2. Are selected dark families embedded in conserved cellular gene neighbourhoods?

**Question:** Do low-query-TM families repeatedly occur beside the same annotated
neighbours across distinct species, compared with fold-matched sequence-dark families?
This can generate specific functional hypotheses for a bounded candidate set.

**Inputs:** original coordinate-verified Pyrodigal headers used by the completed
shadow-overlap audit, admitted-member/90%-cluster mappings, reference contigs and
Bakta/KO annotations on IBEX. The 1,031 versus 1,807 structural sets are identified
in the retained per-model and Pfam tables. Coordinate joins across callers and
coverage of full neighbourhoods must be verified before a test is declared feasible.
Argon's `/ibex/project/c2188/empty-quarter` path is not this project's workbench;
use `/ibex/scratch/projects/c2014/rob/metagenomics`. No expressed-136-family cohort
has been established in this paper.

**Test:** analyse one family-level outcome, with members in at least two distinct
species and callable flanks. Predefine an annotated-neighbour window (for example
within two CDSs), strand and spacing rules; match length, contig completeness,
lineage, member count and annotation availability. Compare conserved flanking KO
identities against within-lineage shuffled neighbourhoods; cluster uncertainty by
family/lineage. Record failures at contig edges separately.

**Interpretation:** repeated conserved context suggests a pathway association; it
does not establish operon transcription or biochemical function. A negative result
may reflect missing annotations or fragmented flanks, not mobile-island biology.
Exclude viral/mobile analyses assigned to the phage companion. Restrict any overlap
with life-on-air content before promoting specific metabolic examples.

**Feasibility:** an IBEX CPU prototype on a small, explicitly chosen set, then scale
only if coordinate and annotation coverage are adequate. No RNA claim or wet-lab
validation can be inferred from the available DNA evidence.

## 3. Does soil-reference composition improve discovery after technical matching?

**Question:** Does a soil-focused reference recover more dark families than a
similarly sized broad reference after controlling sequence length, redundancy,
gene-caller basis and overlap between the databases?

**Available inputs:** preserved SMAG search outputs and catalogue manifests under
`06_darkmatter/specificity_clean/`; compact `smag_content_vs_depth.tsv`,
`content_vs_depth.tsv` and `union_specificity_with_smag_summary.tsv` are pinned in
round5/reviewer_inputs. Existing SMAG absence is 40.47% (hit rate59.53%), while GEM
absence is64.38% (hit rate35.62%). Argon incorrectly treated absence as hits and
invented a4.93% GEM hit rate. The existing union-absence decrement28.74% ->26.83%
is descriptive; it is not itself a new experiment.

**Test:** use families as paired units, audit reference overlap and versions, and
repeat matched-size sequence sampling/searches under identical settings. Stratify
family length and prevalence; retain annotated controls and evaluate incremental
hits beyond the named union. Match or model gene-caller and redundancy differences.
Keep catalogue deduplication and missing-hit definitions fixed across contrasts.

**Interpretation:** a remaining contrast supports useful reference-composition
information conditional on those controls. Disappearance after matching is an
informative negative result. Neither result establishes geographic endemism or a
causal biome effect. No unverified claim about SMAG lacking hyperarid habitats.

**Feasibility:** inputs exist, but the fair search comparison needs IBEX CPU time,
input-identity checks and verification of the catalogue paper before new citation.

## 4. Which structural regions remain unmatched under domain-level calibration?

**Question:** Among confident models with low whole-query similarity but substantial
local hits, are unmatched regions well-supported domains, uncertain extensions,
or failures of reference retrieval?

**Inputs:** existing models, raw Foldseek rows, per-residue confidence and Pfam
coordinates on IBEX; compact model gates and Pfam per-protein tables are retained.
The historical confidence/length reweighting is already computed and is not a new
hypothesis. Its PDB-derived controls overlap the searched references and cannot
supply an exchangeable biological novelty baseline.

**Test:** define domain boundaries independently of the outcome, compare exact
matched targets and coverage, and use length/confidence-matched controls with
homologous reference chains excluded where feasible. Separate no retrieval from
low scores. Treat a model/family as the unit, preserve multiple-domain dependence,
and check boundary sensitivity. Reuse existing models; GPU work is needed only if
new predictions are demonstrably necessary.

**Interpretation:** a validated unmatched domain is a candidate for further study.
A local match or uncertain boundary narrows that candidate set. Non-significance
against a biased positive-control set cannot establish equivalence or "no novel
folds". Avoid the ratio of independently maximised TM scores as a coverage proxy.

**Feasibility:** domain and reference-audit work is needed, so this is more involved
than Argon's suggestion to copy the old standardized-floor table into the paper.

## Completed in this round

The caller substitution, expedition null and paired-profile tests are already done.
A reviewer-driven whole-visit temporal sensitivity was added after the review to
preserve replicate groups. Their outcomes and negative findings are documented in
`reproducibility/darkmatter/round5_extensions/`; they are not future-work promises.
The next scientific choice is which new biological direction to prioritise. My
recommendation is question1, with question2 as a targeted follow-up if callable
neighbourhoods are sufficiently complete.
