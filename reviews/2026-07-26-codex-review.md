# Critical pre-submission review — 2026-07-26

## Recommendation

**Reject in its present form.** This is not yet a Nature/Science-family manuscript. It is an
ambitious computational catalogue organized around a premise that the taxonomic results do not
support, with two load-bearing analyses—the novel-fold ladder and the PICRUSt2 comparison—that are
currently misinterpreted. Several other conclusions are associative, are based on incomplete MAGs,
or use invalid units of replication. The project also remains literally an outline: figures are
placeholders, the Results are bullet-like assertions, and the Methods consist of a pointer paragraph
followed by an instruction to expand it (`paper/main.tex:240-248`). That alone prevents meaningful
submission, independently of the scientific objections below.

This review treats the analyses recorded on 2026-07-26 as current, including the 753-million-protein
catalogue census, complete-BGC sensitivity run, dark-on-MGE negative, kit-ome/decontam screens, and
LRO benchmark (`/home/leechuck/Public/software/empty-quarter-workbench/CHANGELOG.md:37-74`). They do
not repair the central problems. The manuscript PDF passed the hidden-text/instruction integrity
screen; that is not a concern here. In citations below, `CHANGELOG.md` and paths beginning
`scripts/ibex/` are relative to
`/home/leechuck/Public/software/empty-quarter-workbench/`; all other relative paths are under the
manuscript-project root.

## 1. The first four attacks a hostile reviewer will make

### 1. “Known organisms” is not established

The paper equates placement in an existing phylum with a known organism:

> “The tree is ordinary: GTDB-Tk recovers a conventional arid-soil profile (44 phyla, 0 novel
> phyla)—novelty is not phylogenetic.” (`paper/main.tex:40-41`)

> “Discovery here is not taxonomic.” (`paper/main.tex:98-101`)

That inference is false. A direct audit of
`/home/leechuck/Public/software/empty-quarter-workbench/results/soil_gtdbtk/gtdbtk.bac120.summary.tsv`
and `gtdbtk.ar53.summary.tsv` gives:

- 651/15,469 MAGs (4.21%) with a non-empty species assignment;
- 11,466/15,469 (74.12%) with a non-empty genus assignment, leaving 25.9% without one;
- 4,771 classifications explicitly labelled “taxonomic novelty determined using RED”;
- 17 MAGs not classified even at domain/phylum level.

The project’s own plotting code calls a genome with genus but no species “novel species” and one
with family but no genus “novel genus” (`code/make_picrust2_bias_figure.py:37-45`). Thus the data
contradict the title, abstract, A0, and the “same taxa” premise of Part C. “No novel phyla” means that
the MAGs can be placed within GTDB’s high-level scaffold; it does not mean the organisms are known,
represented by close genomes, or functionally predictable. GTDB-Tk also necessarily attempts to
place genomes within its reference taxonomy, so absence of a newly erected phylum is not an
independent biological null.

What is needed is a proper novelty profile: ANI and alignment fraction to the closest reference and
type genome, RED novelty category, fraction assigned at each rank, and phylogenetic branch-length
distribution. The central contrast would be defensible only if “known” were defined in advance at a
biologically relevant resolution—probably species or at least genus—and the functional-dark result
were then shown *within that subset*. On the current numbers, a reviewer can invert the story:
“mostly undescribed species and many undescribed genera also contain undescribed proteins.”

### 2. The structural-dark headline confuses “no hit from a predicted representation” with a novel fold

The abstract claims that approximately 2.53 million proteins “adopt no known fold”
(`paper/main.tex:42-44`). The Results then call 10,007 representatives a “high-confidence
novel-fold floor” (`paper/main.tex:141-157`). Neither statement follows from the analysis.

The 2.53 million count is the no-hit output from ProstT5-predicted 3Di searches, not an observation
that 2.53 million proteins adopt new folds. The explicit ESMFold recheck shows why that distinction
matters: mean pLDDT is 0.507, only 5,990/63,935 models reach pLDDT >= 0.70, and 72% of those
trustworthy models already have a PDB analogue (`RQ/A3/result.md:6-29`). The 10,007 count is merely
the intersection of TM < 0.5 against AFDB and PDB among 62,018 comparable models
(`RQ/A3/result.md:31-39`); it is not filtered to high-confidence structures. The same result file
states explicitly that the defensible claim is the “high-confidence, reconfirmed core (<= a few
thousand reps),” not 10,007 (`RQ/A3/result.md:48-55`). The manuscript says the opposite.

Searching two databases does not validate a low-confidence model. It only gives two no-hit calls
for that model. Low-confidence/disordered/short proteins, fragments, coiled coils, low-complexity
sequences, and domain-boundary errors can all produce low TM scores. Nor is recurrence established
for this structural set: only 46/68,808 identifiers join to the family recurrence table because the
two clustering derivations use incompatible IDs (`RQ/A3/result.md:41-46`). The manuscript’s
“at family scale” language therefore also outruns the join.

A defensible analysis needs the exact intersection of high pLDDT, sufficient aligned length and
coverage, low-complexity/disorder filters, dual-database no analogue, and recurrence in independently
assembled samples. It also needs positive controls (known folds withheld from the databases),
negative/decoy controls, threshold sensitivity, and manual/domain-level inspection of a stratified
sample. Until then, “2.53 million fold-dark candidates” and “at most a few thousand
high-confidence candidate novel folds” are the upper limits of the wording.

### 3. The 2.4-fold PICRUSt2 result is an artefact of the presence threshold and an invalid gold standard

Part C is advertised as the conceptual bridge:

> “16S-based prediction fails on the very taxa GTDB classifies: PICRUSt2 over-predicts a generic
> metabolism 2.4x and is blind to 56.7% of the proteome.” (`paper/main.tex:47-48`)

The implementation reads every PICRUSt2 value that is not the literal string `"0"` as present
(`code/picrust2_bias_analysis.py:71-80`) and defines KO richness from those sets
(`code/picrust2_bias_analysis.py:129-146`). In the curated output, **all 436 samples have exactly
10,543 predicted KOs**. Direct inspection of
`/home/leechuck/Public/software/empty-quarter-workbench/data/processed/functional/picrust2/merged/ko_metagenome_unstrat.tsv`
shows why: 98.24% of all KO-by-sample cells are positive, including values as small as approximately
6.6e-243. Even at thresholds of 0.01, 0.1, 1, 10, and 100, the positive-cell fractions remain
89.45%, 87.10%, 83.81%, 78.72%, and 70.01%, respectively. The 10,543-KO “richness” is therefore
the supported output universe of a continuous ancestral-state predictor, not a credible presence
call. The 2.4-fold ratio is thresholded dust divided by a different measurement.

The purported “MAG-encoded reality” is also not reality. It includes only EggNOG KOs found in
HQ+MQ binned genomes, not all genes in the shotgun metagenomes. Medium-quality MAGs dominate the
catalogue, so absent KOs can reflect incomplete genomes, unbinned organisms, missed genes, or
annotation failure. The 56.7% statistic is exactly the fraction of MAG proteins lacking a KO
annotation (`code/picrust2_bias_analysis.py:156-184`); it is an annotation ceiling shared by any
KO-based method, not an observed PICRUSt2 false negative. The marquee comparison is worse:
PICRUSt2 prediction is correlated with unweighted MAG gene-copy counts rather than CoverM-weighted
gene abundance (`code/picrust2_bias_analysis.py:187-200`), even though a relative-abundance matrix
was loaded.

Most importantly, samples are the only link. The code intersects sample identifiers
(`code/picrust2_bias_analysis.py:83-109`); it never maps amplicon ASVs to the MAGs. NSTI = 0.16 does
not show that PICRUSt2 and GTDB operate on “the same taxa,” and the permissive placement of ASVs
does not validate organism-level pairing. The sentence “where GTDB succeeds yet PICRUSt2 fails on
the same taxa” (`paper/main.tex:85-87`) is factually untested.

This analysis needs to be rebuilt around compositional abundance, a declared detection threshold,
and a genuine shotgun reference: all assembled/unbinned predicted genes or direct read-level KO
profiles, with sensitivity to genome completeness and annotation method. To retain a “same taxa”
claim, ASVs must be linked phylogenetically or by exact/near-exact 16S sequence to MAGs, and
prediction error must be assessed within those pairs. Otherwise the paper can only say that
PICRUSt2 profiles and MAG-derived KO profiles differ at the sample level under two non-equivalent
measurement processes.

### 4. The paper reports capacity, prediction, and association as biology, without a decisive biological result

The project has many large tables but no direct evidence that a dark protein is expressed, that a
predicted enzyme performs the stated reaction, that atmospheric gases support in situ activity, or
that the inferred environmental variables select the implicated families. A Nature reviewer will
see a resource paper whose integration is rhetorical:

- “A hot sand-sea running on trace gases” is inferred from marker genes and taxon-weighting
  (`paper/main.tex:103-111`). The atmospheric analysis itself warns that it measures capacity, not
  flux (`results/atmospheric/README.md:38-43`). Only 974/6,010 guild genomes in the mechanistic
  subset carry validated CBB carbon fixation; the majority may use gas oxidation for maintenance
  rather than growth (`results/atmospheric/integrate/MECHANISM.md:8-18,38-43`).
- UV/desiccation claims use universal DNA-repair genes and MAG-level compartment tests without a
  matched non-desert baseline or control for lineage, genome size, completeness, site, or trip
  (`paper/main.tex:113-117`; `code/uv_survival_analysis.py:66-91`).
- The environmental analysis tests 23,382 selected families, not the complete dark proteome,
  using marginal correlations; its own Methods acknowledge relative abundance and collinear
  chemistry (`RQ/B3/methods.md:58-66`).
- The ten decoded enzymes have calibrated P(correct) of only 0.34-0.48
  (`paper/main.tex:181-187`; `RQ/A4/result.md:5-8`). Calling them “well-supported specific-EC
  enzymes” converts hypotheses with expected error into discoveries.
- The cross-feeding result is a 16-community pilot on complete-environment models. It is an upper
  bound on potential, the detailed SMETANA analysis was infeasible, and the proposed atmospheric
  guild mechanism is null (`results/smetana/CROSSFEED.md:1-22`).
- The new dark-on-MGE result is a useful negative—OR 0.619, not enrichment—but it removes another
  proposed integration route rather than supporting one
  (`/home/leechuck/Public/software/empty-quarter-workbench/CHANGELOG.md:49-52`).

The manuscript needs at least one causal or experimentally grounded result around which the
catalogue is organized: expression/protein detection, isotope or gas-flux evidence, or biochemical
validation of a dark enzyme/fold. Without it, the scale does not substitute for mechanism.

## 2. Claims that outrun the evidence

| Claim | Why it is vulnerable | Evidence required |
|---|---|---|
| “60 sites x 5 expeditions” (`paper/main.tex:38-40,258-263`) | This reads as a 300-cell repeated factorial. The environment table has 547 samples from 59 sites; only 7 sites occur in all five trips. Trip sizes are 206, 37, 128, 85, and 91, and trip 2 covers only seven sites (`RQ/B8/methods.md:17-23`). | Describe the design exactly; provide a site-by-trip matrix and model the unbalanced repeated structure. |
| “Novelty is not phylogenetic” (`paper/main.tex:40-44`) | Only phylum novelty was excluded. Species labels exist for 4.2% of MAGs and genus labels for 74.1%; 4,771 were classified by RED novelty. | ANI/AF/RED novelty profile and an analysis restricted to demonstrably represented species/genera. |
| “~2.53M ... adopt no known fold” (`paper/main.tex:42-44`) | These are ProstT5-3Di no-hit candidates. Explicit models are mostly low confidence and recover known analogues when searched properly. | High-confidence, coverage-filtered explicit structures; calibrated false-negative rate; recurrence and domain-quality controls. |
| “334,810 novel families ... under purifying selection, environment-structured, and largely endemic” (`paper/main.tex:43-45`) | Selection and ecology were measured for selected subsets: dN/dS for 22,767 recurrent families and ecology/global specificity for 23,382 families, not all 334,810. The sentence distributes subset properties over the entire set. | Exact denominators in every clause, inclusion-flow diagram, and a test of whether analyzed families are representative of the full set. |
| “A hot sand-sea running on trace gases” (`paper/main.tex:108-109`) | Marker presence and 16S taxon weighting show genetic capacity, not activity, energy contribution, or carbon flux. | Expression/proteomics, gas uptake or isotope tracing, and ideally comparison with non-desert controls. |
| “a genome pre-loaded for radiation + drying, partitioned by stressor” (`paper/main.tex:113-117`) | NER/BER are nearly universal housekeeping. Photolyase contrasts treat MAGs as independent and do not control lineage, completeness, site, trip, or multiple marker tests. | Taxonomically matched desert/non-desert genomes and a phylogenetic or mixed model with sample/site/trip effects. |
| “Pathogen species = 0.05% ... therefore pristine, low-anthropogenic biome (doubles as a contamination control)” (`paper/main.tex:131-136`) | Seven named pathogens—five *E. coli*, two *A. baumannii*—do not measure reagent contamination, relic DNA, or anthropogenic input (`results/pathogens/PATHOGENS.md:1-14`). | Extraction/sequencing blanks, prevalence/frequency/source-tracking analyses, anthropogenic markers, and appropriate environmental comparators. |
| dN/dS “refutes the random/erroneous-ORF null” (`paper/main.tex:159-163`) | The control is annotated families matched only by family size, not simulated/random ORFs. Pairwise NG86 values are not tests that each family has omega < 1; identical pairs are dropped, paralogy/recombination are not handled, and one valid pair is sufficient. | Matched intergenic/frameshift/shuffled ORFs, coding-potential and codon-usage controls, orthology/recombination filtering, and codon-model likelihood tests with independent support from transcription or translation. |
| “Abundant and widespread, not a rare contamination tail” (`paper/main.tex:165-172`) | Presence is any TPM > 0 in downstream analyses, although the pipeline produced a covered-fraction >= 0.5 presence matrix. A recurrent reagent sequence can also be widespread. | Reanalyse the intended coverage-qualified presence matrix, include blanks and competitive mapping, and show robustness to detection thresholds. |
| “active arms-race selection at recognition interfaces” (`paper/main.tex:175-179`) | The 148-family mean-NG86 > 1 list is a screen, not a positive-selection test. Only 12 families had sufficient depth for HyPhy; pocket enrichment is null (p = 0.66), and only eight confident structure-function hits exist (`results/darkmatter_ecology/positive_selection/POSITIVE_SELECTION.md:1-3,49`; `results/darkmatter_ecology/positive_selection/STRUCT_FUNCTION.md:3,22-27`). | Codon-model evidence with recombination control, corrected site tests, confident structures, and a predeclared interface-enrichment test. |
| “10 well-supported specific-EC dark enzymes” (`paper/main.tex:181-187`) | P(correct) = 0.34-0.48 implies that several top calls are expected to be wrong. The new LRO summary reuses the existing calibration tests and compares only with random; it has no homology baseline (`CHANGELOG.md:53-58`; `scripts/ibex/168_decode_benchmark.py:4-31`). | A disjoint benchmark of the deployed model, homology/context-only baselines, dark-like decoys, blind annotation, and biochemical validation. |
| “2.4x over-prediction” and “same taxa” (`paper/main.tex:192-199`) | All PICRUSt2 samples contain all 10,543 KOs under the nonzero rule; the comparator omits unbinned and missing MAG sequence; ASVs are not linked to MAGs. | Thresholded abundance benchmarking against a comprehensive shotgun profile plus direct ASV-MAG taxonomic pairing. |
| Open pangenomes are “the mechanism” of PICRUSt2 failure (`RQ/C2/result.md:5-26`) | The analysis uses HQ+MQ MAGs but defines core as presence in >=95% of genomes. With incomplete MAGs, a truly universal gene is unlikely to be observed in 95% of ten or more bins, depressing the core and inflating openness. No analysis shows housekeeping closed while gas/osmolyte/BGC/dark loci are open, and no held-out prediction connects openness to PICRUSt2 error. | HQ-only/completeness-adjusted pangenomes, marker recovery simulations, class-specific saturation, and prospective held-out gene-content prediction per species. |
| “21.4% have no homolog anywhere sequenced” and “endemism plateaus” (`paper/main.tex:214-219,310-321`) | The analysis means no DIAMOND `--sensitive` hit at e < 1e-3 in selected catalogues. There are no coverage/identity/low-complexity criteria and the databases are non-nested, so three large-catalogue points do not establish saturation. | Sensitive HMM/profile and structural searches, coverage and complexity controls, nested subsampling/rarefaction within catalogues, and explicit database versions. |
| “Endemic families are under stronger purifying selection” (`paper/main.tex:217-219`) | Within-EQ occupancy gives the stated direction, but global MGnify endemics have median dN/dS 0.393 versus 0.390 for shared families—slightly *less*, not more, constrained (`results/darkmatter_ecology/GLOBAL_ENDEMISM.md:1-8`). | Separate the within-EQ and global definitions and report effect sizes; do not claim a global selection result that has the opposite direction. |
| “98.81% novel GCFs” is a robust natural-product novelty result (current 2026-07-26 state) | The complete-BGC rerun removes the edge-fragmentation explanation: 6,298/6,374 complete-BGC GCFs lack a MIBiG member (`CHANGELOG.md:63-68`). But MIBiG is only ~2,502 characterized BGCs, T = 0.4 is a single clustering choice, and “not co-clustered with MIBiG” does not imply a novel molecule or globally novel predicted BGC. The methods themselves identify BiG-FAM as the necessary stricter comparison (`RQ/F2/methods.md:83-94`). | BiG-FAM/global predicted-BGC comparison, threshold sensitivity, leave-one-out MIBiG recovery calibration, and chemistry/expression before claiming a product reservoir. |
| “test with explicit nulls” and Box 2’s nulls are defeated (`paper/main.tex:48-50,402-406`) | The supposed contamination null uses pathogens and an underpowered decontam test whose positive control is also null. The spurious-ORF null lacks random ORFs. The open-pangenome “null” is not a PICRUSt2 validation. The global-catalogue test is not a saturation analysis. | Real negative and positive controls for each pillar, defined before inspecting outcomes. |

## 3. Statistical and methodological problems

### Confounding and design

1. **Trip is batch, and site/trip/compartment are nested and unbalanced.** The temporal analysis
   acknowledges that every trip is also a sequencing batch (`RQ/B8/methods.md:54-63`) but then uses
   unrestricted permutations and compares separate one-factor PERMANOVAs for trip and site
   (`RQ/B8/methods.md:45-48`). The conclusion that “time is by far the dominant structuring factor”
   (`RQ/B8/result.md:48-59`) is not estimable from those models. R-squared values from factors with
   4 versus 58 degrees of freedom are not a contest, and “per degree of freedom” is not a valid
   resolution.

2. **Environmental predictors are proxies for each other and for geography/taxonomy.** Fifteen
   separate Spearman screens are run, while Ca, Cl, S, Na, salinity, climate, and location are
   collinear (`RQ/B3/methods.md:19-31,58-66`). Regressing ranked latitude and longitude does not
   “separate genuine environmental chemistry from spatial autocorrelation,” contrary to
   `RQ/B3/methods.md:39-46`; it removes only a planar spatial trend. There is no spatial covariance,
   site blocking, trip blocking, compartment adjustment, lineage adjustment, or MAG-community
   composition adjustment. “Recruited” and “filtered” should be replaced with “positively” and
   “negatively associated.”

3. **Genome quality and taxonomic composition confound nearly every MAG-level comparison.**
   Completeness, contamination, protein count, contig fragmentation, taxon, sample sequencing
   depth, and binning recoverability affect annotation, BGCs, KO richness, pangenome core size,
   marker prevalence, and dark-protein counts. These covariates are rarely modeled jointly.

4. **Atmospheric abundance is inferred across taxonomic ranks.** The 16S weighting assigns guild
   propensity by genus where available and by phylum otherwise; only about 64% of abundance is
   genus-resolved (`results/atmospheric/abundance/README.md:8-16`). This assumes within-genus or
   within-phylum conservation of guild status and cannot support “community-wide” genomic
   prevalence or flux.

### Multiple testing

1. The environmental scan applies BH separately within each of 15 axes
   (`RQ/B3/methods.md:25-31`). There is no experiment-wide correction or multivariate shrinkage,
   despite reporting thousands of discoveries and then selecting examples from those discoveries.
   Axis-level collinearity makes the repeated lists non-independent.

2. Marker-panel and compartment analyses do not consistently report family-wise correction.
   Testing many repair, osmolyte, gas, resistance, and BGC markers and then highlighting the most
   extreme p-values creates an unreported garden of forking paths.

3. Extremely small p-values largely report sample size, not importance. Examples include
   p < 1e-300 for millions of gene-family observations and p < 1e-67 for occupancy categories.
   Effect sizes, uncertainty intervals, and blocked/resampled estimates are more relevant.

### Pseudoreplication

1. The 109,377 across-site/across-trip pairwise Bray-Curtis values are treated as independent in a
   Mann-Whitney test (`RQ/B8/result.md:31-45`). Every sample occurs in hundreds of pairs. The
   reported p-value is invalid. Resample at site/trip level or fit an appropriate repeated-distance
   model.

2. MAGs from the same sample and closely related MAGs are repeatedly treated as independent in
   repair-marker, kit-ome, MGE, pathogen, and metabolic comparisons. The sampling unit is normally
   sample/site or lineage, not protein or MAG.

3. The kit-ome and MGE Fisher tests use proteins as independent observations. Proteins are clustered
   within genomes, samples, and taxa. The ORs may be descriptive, but the nominal inferential
   precision is meaningless.

4. Pairwise dN/dS values within a family share sequences, and families can share taxa and samples.
   Averaging up to 80 pairs does not create independent evidence. The pipeline also drops identical
   sequence pairs and accepts a family with one valid divergent pair
   (`code/ibex/94_darkfam_dnds.sh:86-122`), which changes the estimand and selects on divergence.

### Circularity and ascertainment

1. Dark families are discovered by UniRef50 no-hit and recurrence, then only recurrent families are
   tested for selection and ecology. Recurrence and apparent constraint cannot be generalized to
   singletons or all 2.84 million representatives. The headline nevertheless assigns these
   properties to all 334,810 families.

2. Structural “validation” starts from candidates selected by a predicted-3Di no-hit. The explicit
   structure analysis usefully shows that most fail, but the remaining no-hit set is still not an
   independent validation unless model confidence and known-fold recovery are calibrated on held-out
   controls.

3. The 2026-07-26 LRO “benchmark” reads the existing calibration table and summarizes stored ranks
   (`scripts/ibex/168_decode_benchmark.py:4-18`). It therefore demonstrates the behavior of the
   calibration harness on its own tests, not independent external validation of the ten dark-enzyme
   calls. A random ordering is an inadequate baseline. The homology-only job was still running at
   the stated cutoff (`CHANGELOG.md:69-74`).

4. The global-specificity analysis begins with UniRef50-dark queries, then presents failure against
   other partly overlapping reference catalogues as a general plateau. Database construction,
   clustering, taxonomic composition, and search sensitivity are part of the ascertainment and must
   be modeled.

5. Candidate examples—osmolyte enzymes, lectins, structural targets—were selected after inspecting
   several correlated evidence streams. Any reported confidence must reflect this selection, and
   final validation must use held-out evidence.

### Database and reference biases

1. **Sequence novelty:** an e-value-only DIAMOND criterion is length dependent. There are no
   minimum query/subject coverage, identity, low-complexity, or profile-HMM criteria
   (`results/global_specificity/README.md:3-7`). “No homolog” should be “no hit under this search
   configuration.”

2. **Structural novelty:** AFDB is predicted, versioned, and confidence-biased; PDB is
   experimentally sparse and taxonomically/functionally biased. TM < 0.5 on a low-confidence model
   is not a biological class.

3. **BGC novelty:** MIBiG represents characterized clusters, not global biosynthetic sequence space.
   The 98.81% result is robust to edge fragmentation but not to reference-set choice or clustering
   threshold.

4. **Taxonomic novelty:** GTDB reference placement and empty species labels measure different things.
   The manuscript uses high-rank placement as evidence that reference genomes “recognise” the
   organism.

5. **Functional novelty:** failure to receive an EggNOG KO is not necessarily novel biology. It
   mixes real novel function, divergent homologues, fragments, annotation sensitivity, and
   proteins outside KO’s scope.

### Missing or failed controls

1. **Contamination:** The new genus-list screen reports depletion of dark proteins in putative
   reagent genera, but reagent-genus membership is a crude proxy and is taxonomically confounded
   (`CHANGELOG.md:38-40`). The frequency screen has only 69 overlapping samples and fails its own
   positive control—0/63 Salter-genus MAGs flagged (`CHANGELOG.md:41-47`). Its implementation drops
   zero-abundance observations and correlates concentration only among samples where each MAG was
   observed (`scripts/ibex/171_decontam_frequency.py:67-92,202-227`), which is not an adequate
   frequency method. Box 2 cannot say the kit-ome null is defeated.

2. **Spurious ORFs:** The size-matched annotated control is deterministic and matched only on log2
   family size (`code/ibex/93_darkfam_member_cds.sh:47-62`). There is no intergenic, shuffled,
   frameshifted, or low-coding-potential negative control.

3. **Ecology:** There is no matched known-family null controlling family size, protein length,
   taxon, mean abundance, mappability, and prevalence. Such a null is necessary to establish that
   dark families are unusually structured rather than merely tracking their host lineages.

4. **PICRUSt2:** There is no positive-control benchmark showing recovery of known sample functions
   under the same thresholds, no mock-community calibration, and no comprehensive shotgun truth.

5. **Structure:** There is no held-out known-fold recovery curve across pLDDT, length, disorder, and
   TM thresholds.

6. **BGCs:** There is no leave-one-out MIBiG positive control showing that the chosen BiG-SLICE
   threshold recognizes known related clusters.

The clean negative results should remain negative. Dark proteins are depleted on MGEs, and the
atmospheric guild is not a preferential donor in the pilot. Recasting either as indirect support
would make the paper less credible.

## 4. Missing analyses ranked by value per compute

The ranking below is by expected improvement to the paper per unit of computation, not by novelty.
`$M` denotes `/ibex/scratch/projects/c2014/rob/metagenomics`.

| Rank | Analysis | Why it matters | Runnable now? Concrete inputs |
|---:|---|---|---|
| 1 | **Audit taxonomic novelty at species/genus resolution and condition the headline analyses on “knownness.”** | This is cheap and decides whether the title survives. Report rank assignment, ANI, AF, RED and reference distance; then compare dark burden within genuinely reference-proximal versus novel MAGs while controlling completeness and genome size. | **Yes; existing data, negligible compute.** GTDB summaries at `/home/leechuck/Public/software/empty-quarter-workbench/results/soil_gtdbtk/gtdbtk.{bac120,ar53}.summary.tsv`, plus `$M/04_drep_self_vs_single/out/data_tables/Cdb.csv` and GTDB-Tk closest-reference fields. |
| 2 | **Rebuild the PICRUSt2 benchmark.** | The current 2.4x claim is invalid. First plot abundance/detection curves across thresholds; then compare compositional KO abundance against all shotgun genes, not only KO-annotated HQ+MQ bins. Link ASVs to MAGs where possible and stratify error by actual reference distance/NSTI. | **Threshold/sensitivity work is runnable now, trivial compute** from `.../picrust2/merged/ko_metagenome_unstrat.tsv`, `results/picrust2_bias/mag/mag_ko_{count,relab}_by_sample.tsv`, and `per_sample_concordance.tsv`. **A comprehensive shotgun truth requires new compute but no new sequencing** from `$M/04_protein_catalog/assembly_proteins/<sample>/<sample>.faa.gz`, assemblies/reads, and the existing functional-annotation outputs. |
| 3 | **Use the intended coverage-qualified dark-family presence calls and fit the sampling hierarchy.** | This tests whether the “abundant, persistent, environment-structured” pillar survives removal of trace mappings and invalid permutations. Use covered fraction >=0.5, competitive mapping, abundance/prevalence sensitivity, blocked permutations or mixed models for site/trip/compartment, PERMDISP, and adjusted multivariate R-squared. | **Yes; existing data, low-to-moderate compute.** `$M/06_darkmatter/darkfam_ecology/darkfam_presence.tsv`, `darkfam_abundance.tsv`, per-sample CoverM TPM/covered-fraction outputs from `code/ibex/91_darkfam_coverm_array.sh`, and `results/darkmatter_ecology/env_matrix.tsv`. |
| 4 | **Correct the pangenome analysis for MAG incompleteness.** | A >=95% core threshold on MQ genomes almost guarantees false absences. Restrict to near-complete genomes, model per-genome detection probability, simulate recovery of universal markers, vary core threshold, and test class-specific openness. Add held-out gene-content prediction rather than declaring a mechanism from openness. | **Yes; existing data, moderate compute.** `$M/04_drep_self_vs_single/out/data_tables/Cdb.csv`, `$M/04_drep_self_vs_single/inputs/*.fa.gz`, `$M/04_protein_catalog/assembly_proteins/`, per-sample CheckM2 `quality/*/quality_report.tsv`, and `$M/06_pangenome/out/{species_openness,family_table,species_members}.tsv`. |
| 5 | **Calibrate the structural ladder and report only a confidence-qualified set.** | This is necessary before any “novel fold” claim. Intersect pLDDT, aligned length/coverage, disorder/complexity, AFDB/PDB TM, and independent recurrence; benchmark the same pipeline on held-out known folds. Manually audit a stratified sample. | **Yes; existing sequences/results, moderate-to-high compute, no new data.** `$M/04_structure/foldseek_results/{all_query_ids,afdb_hit_ids,struct_nohit_ids}.txt`, `$M/06_esmfold_novelfold/{best_afdb_tm.tsv,best_pdb_tm.tsv,pdb/plddt_shard*.tsv}`, and the family/member tables. Some refolding/domain segmentation may be needed. |
| 6 | **Replace the global “plateau” plot with calibrated, nested search sensitivity.** | It would turn a database-size cartoon into a specificity analysis. Apply length/coverage/complexity filters; search profiles/HMMs; subsample one catalogue at nested depths; report database versions and uncertainty. | **Mostly yes; existing queries and downloaded databases, moderate compute.** `results/global_specificity/`, the 23,382 family representatives, and the MGnify/GEM/five-desert databases already used on IBEX. Profile-HMM construction is new computation, not new data. |
| 7 | **Harden the complete-BGC result.** | The complete-only rerun resolves fragmentation but not reference bias. Query complete BGCs against BiG-FAM/global predicted BGCs, vary T, and perform leave-one-out MIBiG recovery. | **Yes if the planned BiG-FAM model/database is available; moderate compute.** `$M/06_metabolism/bgc_novelty/results_complete/`, `$M/06_metabolism/antismash/<MAG>/*.region*.gbk`, and the existing MIBiG 3.1 files. Otherwise the BiG-FAM model must be obtained; no new biological samples are needed. |
| 8 | **Finish a genuinely independent decoder benchmark.** | A random baseline and recycled calibration ranks cannot justify “well-supported enzymes.” Complete the homology/proximity/recurrence baselines, evaluate the deployed scoring version on a disjoint train/test split, use dark-like decoys, and report calibration/reliability, not only top-k. | **Yes; existing data and the already launched job 172, low-to-moderate compute.** `$M/06_darkmatter/calibration_sub/lro_tests.tsv`, Gapsmith reaction/candidate data, and `$M/06_darkmatter/benchmark_i1/`. It still cannot validate the ten biological calls without new evidence. |
| 9 | **Run contamination analysis with the actual controls and all available concentration samples.** | This is required before using abundance as a contamination argument. Include zeros, widen the CoverM matrix from 69 toward the 433 concentration samples, block trip/batch, and map extraction/sequencing controls competitively to the dark catalogue. | **Partly runnable now.** `$M/metadata/dna_conc.tsv`, `$M/04_functional/<sample>/abund.tsv`, raw reads/assemblies for `SNTC` and `SPositive`, and `$M/06_darkmatter/decontam/`. If the controls were not processed with the soil libraries or lack usable reads, **new control sequencing is required**. |
| 10 | **Add experimental evidence around one dark-function mechanism.** | This is the only change likely to move the work from a strong computational resource toward a general-interest mechanistic paper. Priorities: heterologous expression/assay of predeclared decoded enzymes; metatranscriptomics/proteomics of target families; H2/CO uptake and 13CO2 incorporation under dry-soil conditions. | **Needs new work and biological material.** Existing candidates, MAGs, predicted structures, and primer/construct targets are sufficient for design, but the evidence does not exist on IBEX. |

Two additional low-cost analyses should accompany the above: a full inclusion/denominator flow chart
from 753 million proteins to every analyzed subset, and a matched-known-family null for every
dark-family ecology/selection comparison.

## 5. Framing and venue

### The proposed through-line is not supported

“Known organisms, unknown biology” currently fails on both nouns.

- The organisms are not known at organism-level resolution: 95.8% lack a species assignment and
  25.9% lack a genus assignment.
- The unknown is mostly *unknown annotation*, not demonstrated unknown biology. No-hit proteins,
  predicted no-hit folds, and unassigned KOs are informative computational categories, but they do
  not establish novel biochemical activity.
- The claimed bridge between the two is invalid because PICRUSt2 and MAGs are not paired by taxon
  and are compared using non-equivalent presence rules.
- The integration is weakened by its own clean negatives: dark proteins do not supply the trace-gas
  enzymes, are depleted rather than enriched on MGEs, the atmospheric guild is not a preferential
  metabolic donor, and structured dark families show no class-level osmolyte enrichment. Those
  results may be worth reporting, but they mean the three parts are adjacent surveys rather than a
  demonstrated mechanism.

A more accurate current framing would be: **“A genome-resolved catalogue of taxonomic and protein
sequence novelty in Rub al-Khali soils.”** That is less dramatic because it matches the evidence.
The phrase “classification without comprehension” should not be used while the genomes are mostly
unassigned at species level and comprehension is assessed through an invalid PICRUSt2 comparison.

### What would have to be true for Nature or Science

For a Nature/Science-family claim, at least the following would need to hold simultaneously:

1. The paradox survives a rank-appropriate definition of “known”: many MAGs are close to represented
   reference species/genomes, yet their gene content is demonstrably unpredictable.
2. A calibrated, confidence-filtered structural analysis establishes a substantial number of
   genuinely new recurrent folds rather than low-confidence no-hits.
3. The reference-prediction result is measured against a comprehensive shotgun truth in matched
   taxa and remains large after threshold, abundance, completeness, and annotation sensitivity.
4. Ecological organization survives blocked/hierarchical statistics, spatial autocorrelation,
   taxonomic composition, compositionality, and matched known-family controls.
5. At least one dark-function or trace-gas mechanism is supported by expression, flux, isotope
   incorporation, or biochemical assay.
6. The findings replicate beyond this single survey—through an independent Rub al-Khali cohort or a
   predeclared desert/non-desert comparison—rather than relying on many correlated samples from one
   programme.

Without those conditions, the credible target is a specialist microbiome/genome-resource venue.
With rigorous reanalysis but no experimental mechanism, the ceiling is more plausibly an ISME
Journal/Microbiome/Genome Biology or possibly Nature Communications-style resource paper, depending
on how much of the central signal survives. The current project plan itself already recognized that
a single strong result plus wet-lab validation would be required for a Nature-level main story
(`plans/darkmatter-proteins.md:176-193,224-237`).

## 6. Apparent real errors and internal contradictions

1. **The taxonomic premise is numerically wrong.** “Databases recognise every microbe”
   (`paper/main.tex:28-29`) is contradicted by 651 species-level assignments and 17 unclassified
   genomes. “Zero novel classes” should also be reconciled with the actual GTDB rank strings rather
   than inferred from placement.

2. **Representatives are called proteins.** The 2.84 million sequence-dark objects are clustered
   representatives in the ladder (`paper/main.tex:141-144`), but the abstract calls them proteins
   (`paper/main.tex:42-44`). Those denominators are not interchangeable.

3. **The 10,007 structural count is mislabeled high-confidence.** The manuscript calls it a
   “high-confidence novel-fold floor” (`paper/main.tex:150-156`), while the source analysis states
   that only 9.4% of models are high-confidence and that the confidence-filtered, reconfirmed set is
   at most a few thousand (`RQ/A3/result.md:16-20,48-55`).

4. **The ladder table has a category error.** The 68,808 row is labelled “0 AFDB *sequence* hit”
   (`paper/main.tex:281-293`), although it is a structure/3Di-search stage.

5. **PICRUSt2 KO richness is constant by construction.** All 436 curated rows have
   `n_pic_ko = 10543`; “10,543 KOs/sample” is not an empirical richness distribution. It results
   from treating arbitrarily small nonzero predictions as presence
   (`code/picrust2_bias_analysis.py:71-80,129-146`).

6. **“Same taxa” is not implemented.** Only sample identifiers are intersected
   (`code/picrust2_bias_analysis.py:83-109`). No ASV-to-MAG mapping or taxon-paired comparison exists.

7. **The marquee abundance correlation uses the wrong MAG quantity.** It correlates PICRUSt2
   abundance with `mag_count_s`, not the loaded `mag_relab_s`
   (`code/picrust2_bias_analysis.py:187-200`).

8. **The dark-family presence criterion is not the declared one.** The aggregation pipeline creates
   `darkfam_presence.tsv` using covered fraction >= 0.5
   (`code/ibex/92_darkfam_aggregate.sh:17,62-63`), but `code/darkmatter_abundance.py:28-43`,
   `code/phase_b_ecology.py:74-78`, and `RQ/B8/run.py:125-135` use any TPM > 0. Direct re-audit of
   `results/darkmatter_ecology/darkfam_abundance.tsv` changes median prevalence from 52 samples at
   >0 to 31 at TPM >=1, 17 at >=10, and 3 at >=100. The contamination/persistence narrative is
   threshold sensitive.

9. **The sign explanation in the endemism file is wrong.**
   `results/darkmatter_ecology/ENDEMISM.md:13-15` reports Spearman r = +0.176 but immediately
   interprets “negative => ...”. The numerical group comparison supports lower dN/dS for local
   endemics, but the parenthetical sign statement is a typo/error.

10. **Within-EQ and global endemism are conflated.** Global MGnify-endemic families have dN/dS
    0.393 versus 0.390 for shared families (`results/darkmatter_ecology/GLOBAL_ENDEMISM.md:7`).
    The manuscript’s unqualified statement that endemic families are more strongly constrained is
    false for the global comparison.

11. **The positive-selection/interface files contradict each other.** `RQ/E2/result.md` asserts
    localization to interfaces and later says the current evidence cannot support that conclusion;
    related summaries report 23/332 versus 39/520 pocket/background sites. Those denominators need
    one audited source of truth before any result is quoted.

12. **The pangenome result is structurally biased toward “open.”** A core threshold of >=95% across
    incomplete HQ+MQ MAGs (`RQ/C2/methods.md:45-48`) treats missing sequence as biological absence.
    The reported median core/pan ratio of 0.065 (`RQ/C2/result.md:13-20`) is therefore not currently
    interpretable as true gene-content turnover.

13. **The catalogue-wide result does not end the lower-bound caveat as stated.** The changelog says
    it “ends” the caveat while reporting 37,828,539 clusters uncovered
    (`CHANGELOG.md:59-62`). The script explicitly separates uncovered clusters
    (`scripts/ibex/169_catalog_wide_ladder.py:16-21,149-169`), and the submitted wrapper supplies no
    `EXTRA_HIT_IDS`/`EXTRA_NOHIT_IDS`; the Python defaults load only
    `delta_vs_uniref50/{hit_ids,nohit_ids}.txt`
    (`scripts/ibex/169_catalog_wide_ladder.py:36-46`;
    `scripts/ibex/169_catalog_wide_ladder.sh:10-16`). Until those 37.8 million are reconciled, the
    result is 68.65 million dark clusters among the *covered* clusters, not a complete
    753-million-protein dark ladder. It is also a sequence-only layer, not a catalogue-wide
    structure/selection/ecology ladder.

14. **The decontam implementation conditions on detection.** Zero observations are omitted at load
    time and each MAG is correlated only across samples where it was detected
    (`scripts/ibex/171_decontam_frequency.py:67-92,202-227`). This can remove the very
    low-concentration/high-prevalence pattern the method is meant to detect. The failed positive
    control confirms that the null result is not validating.

15. **The new decode “benchmark” is not new held-out validation.** It summarizes already persisted
    LRO calibration ranks and lacks candidate lists for alternative baselines
    (`scripts/ibex/168_decode_benchmark.py:4-31`). “5,623,008 held-out tests” should not be described
    as independent validation of the final decoder or the ten dark candidates.

16. **Some supposedly reproducible local scripts contain expired temporary paths.**
    `code/uv_survival_analysis.py:18` and `code/phase_c_struct_integrate.py:25` hard-code a
    `/tmp/claude-1000/.../scratchpad` path that is absent. These scripts cannot be rerun from the
    repository as written.

17. **The reproducibility accounting is stale and overstates what is verified.**
    `handover.md:109-118` says 51 claims (37 COMPUTE, 6 TABLE, 8 PROVENANCE), while
    `reproducibility/claims.tsv` currently contains 52 (38, 6, 8). The 38 computational checks pass,
    but TABLE and PROVENANCE entries mainly verify file existence or recorded summaries, not the
    scientific validity of heavy analyses (`reproducibility/README.md:33-46`).

18. **Project status is internally stale.** `README.md:35-40` says 27/40 questions are answered and
    C2 is running; `RQ.md:15-22` says 29/40 and C2 is done; the manuscript still labels BGC novelty
    and pangenome openness as planned (`paper/main.tex:125-129,201-206`), while both have results.
    `RQ/INDEX.md` still lists the dark-on-MGE, complete-BGC, and formal decontam work as pending or
    deferred despite the 2026-07-26 runs. These discrepancies are not cosmetic when claims,
    denominators, and caveats changed.

19. **The BGC result’s old and new caveats coexist.** `RQ/F2/result.md:23-34` calls 99.29% a
    fragmentation-inflated ceiling pending a complete-only run. The 2026-07-26 complete-only result
    gives 98.81%, but the RQ and manuscript have not been updated. The residual BiG-FAM/reference
    caveat remains even though the fragmentation caveat can now be retired.

20. **The submission has no full Methods.** Versions, thresholds, sample exclusions, statistical
    units, multiple-testing families, permutation constraints, and provenance are scattered across
    RQ files and scripts. A one-paragraph tool list (`paper/main.tex:240-248`) is not a Methods
    section and makes the headline results unauditable within the manuscript.

The core editorial decision should be made before polishing prose: either rebuild the paper around
a narrower, valid catalogue claim, or obtain the organism-level reference comparison,
confidence-qualified structural result, corrected PICRUSt2 benchmark, hierarchical ecology, and
experimental validation needed for the present framing. Cosmetic revision cannot solve these
problems.
