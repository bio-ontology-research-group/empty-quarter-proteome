# Independent design checks

These are checks against early drafts, not conclusions about the final dossiers.
Resolve them against executed code and output before accepting any claim.

## H1

- A common denominator does not automatically create a negative correlation:
  with fixed numerators both fractions decrease as the denominator increases.
  A Poisson offset also does not automatically remove compositional coupling.
- An independent permutation of dark and RBH counts can violate their bounds
  relative to total proteins. It can illustrate an algebraic effect but is not
  automatically a valid biological null. State the direction/tail correctly and
  use the plus-one correction for Monte Carlo tests.
- Check catalogue versus AAI/CheckM2 caller denominators and identifier joins.
  CheckM2 partial flags cannot silently stand in for the catalogue gene calls.
- The actual 237-genome table has eight named families (one contains 140 genomes)
  plus nine genomes without a family. Genus is missing for 135 genomes. Missing
  taxa are not one biological clade. A few-cluster asymptotic sandwich p-value is
  not convincing ancestry control. Use the available phylogeny with tested tip
  mapping or explicitly descriptive/block-robust estimates.
- Regex extraction of tip/terminal lengths is not a phylogenetic covariance
  model. Failure to reject a class coefficient does not show complete explanation.

## H2

- Require actual same-gene coordinate and annotation joins, flank opportunity,
  family as the unit, and independent species support. No function, operon,
  expression or adaptation inference from proximity alone.

## H3

- SMAG database access is denied; old top-one rows cannot implement a matched
  database experiment. A GEM-only demonstration of top-one bias is useful but
  does not establish the feasibility of the requested SMAG versus GEM contrast.
- Rebuild a justified bounded public SMAG subset and execute a fresh comparison.
  Preserve primary-source download and caller provenance. Quantify any caller,
  overlap, redundancy and length differences left in the bounded estimand.
- Annotated-family sensitivity need not respond proportionally to dark-family
  sensitivity. Do not use controls to assert technical confounding is eliminated.
- Avoid causal biome or evolutionary explanations from an observational contrast.

## H4

- Confidence-threshold runs are confident regions, not validated domains, and
  low confidence alone does not establish disorder.
- Whole-chain hits overlapping a region do not provide region-specific TM scores.
  If testing region recovery, crop independently defined regions and actually
  search/align them. Use query-aligned residue counts, not independently maximized
  qTM/alnTM ratios, as coverage.
- Exclusion of a source PDB identifier is not exclusion of sequence homologues.
  Preserve retrieval truncation and reference-leakage limitations explicitly.
- Test PDB confidence scale, residue mapping, control eligibility and sampling.
  Define both positive and negative outcomes without a new-fold interpretation.

### H4 first dossier (unaccepted)

The first DOSSIER.md labels proximity-derived regions bona fide novel domains,
claims linker reorientation and disordered-tail mechanisms without testing them,
and even sets a “Confirmed Novel Structural Domains” go criterion. These are
unjustified. It also misrepresents the current round5 manuscript as committing
these errors: check actual current text, which already forbids fold novelty and
equivalence inference. The table labels controls “Homology-excluded PDB” although
only source-PDB identifiers were excluded. The full-test 30%-identity exclusion
was not performed in the pilot. A nonsignificant difference is not equivalence.
Require a genuinely cropped-region search and modest reference-conditional
recovery estimand; make no causal/mechanistic/new-fold claim from this pilot.

### H1 first dossier (unaccepted)

The dossier reverses the sign of its own output: `results/decoupling_models_summary.txt`
records the independent-count common-denominator null mean **+0.3017**, CI
[+0.1789,+0.4194], whereas the dossier prints **-0.3017** and a negative interval.
The script's interpretation text also calls a positive null negative. Neither
this artificial null nor raw-count correlations establishes an exclusively
mathematical explanation of the observed negative ratio association. Reject the
“100% artifact”, “refuted”, “vanishes completely”, “definitive debunking”, and
causal fragmentation claims. A nonsignificant coefficient does not establish zero
or absence. Audit all numerical prose against saved outputs, not model memory.
The original requested hypothesis concerned continuous within-lineage associations;
a post-hoc shift to asserting quality dominates needs evidence and honest scope.
Run actual phylogenetic covariance/tip-mapping feasibility and a defensible bounded
association analysis. Characterize supported ranges/uncertainty and leave
macroevolutionary gain/loss claims outside this experiment. Do not label solitary
family instances as fragmented gene calls without checking those calls.

### H2 early code

`03_analyse_conservation.py` describes within-lineage/contig permutation but
actually collects all nonempty flank KOs across all families into a single pool
and shuffles globally. This destroys host/contig-specific annotation composition.
Implement the specified constrained null or narrow and label the analysis honestly;
a host-preserving random focal-position/whole-window null is preferable. A family
with fewer than two callable species is not an eligible negative. Code currently
counts a species as callable when even one of four flanks exists, whereas DESIGN
requires nontruncated windows. Match/condition family/species breadth, annotation
opportunity and host lineage for the arm comparison. Preserve cluster/member
mapping and cross-caller identity verification; same species repeated MAGs must
not add independent evidence. Remove candidate novel-fold/operon mechanisms from
claims based only on neighbourhood recurrence.

### H1 v2 code checks (before final dossier)

`02_pgls_ancestry_models.py` constructs `y_tilde` as shape (n,) and
`ones_tilde * mu_null` as (n,1); their subtraction broadcasts to (n,n). Its
reported R-squared denominator is therefore wrong. Flatten consistently and test
GLS against a trusted implementation for fixed covariance. ML `RSS/n` belongs in
the profile likelihood; standard fixed-covariance coefficient t intervals require
`RSS/(n-p)`. Report profiled-lambda intervals as conditional unless lambda
uncertainty is propagated. Verify covariance eigenspectrum/regularization and
state midpoint rooting plus correlation normalization assumptions.

### H4 v2 code checks (before final dossier)

Crops are now actually searched. However, whole-chain values are read from a
historical search while crops use `-e 10 -s 9.5 --max-seqs 50`; demonstrate identical
settings/reference or execute same-run whole-chain controls. Use the SAME exclusion
tier for both baseline and crops. Query span `qend-qstart+1` is not the number of
aligned residues when internal gaps occur. Record aligned query characters from
an alignment field or call it span coverage explicitly. Distinguish no eligible
region from zero retrieval, and raw retrieval from no remaining hit after exclusion.
A family-level conditional recovery proportion is estimable; don't confuse the
length-stratified pilot's unweighted estimate with the full-cohort proportion.

H1 v2 dossier also asserts `N_dark + N_RBH <= N_prot` as a disjoint partition.
This is not established: sequence-dark annotation status and RBH to an external
GTDB genome need not be mutually exclusive, and the callers differ. Verify actual
intersection or remove bounded-composition proof. AAI comparators are GTDB genomes;
do not call them cultivated/isolate references without metadata proof. The new
post-hoc positive-AAI hypothesis must be labelled exploratory, not a prospective
one-sided confirmatory test. Corrected ancestry modelling can establish a viable
association experiment without these causal or partition assumptions.

H2 critical identifier bug: in `02_extract_neighbourhoods.py`,
`bakta_by_coord[(start, stop, strand)]` drops contig identity. Coordinates are
only unique WITHIN a contig, so this can overwrite or misassign annotations across
contigs. The correct key must include proven source-to-Bakta contig mapping and
bin/sample, then start/end/strand. Audit duplicates and ambiguous joins fail closed.
Its `exact_coordinate_matches=554` counter only measures discovery of the focal
protein in original FAA headers; it does NOT measure Pyrodigal-to-Bakta coordinate
matches. Report the latter separately among eligible binned neighbours, separating
missing source files, no cross-caller match, mapped/noKO and KO-present. Also never
add blank species as a second species when aggregating KO recurrence.

Independent small-table eigen-whitening check is saved in
`root_checks/audit_h1_gls.py` and `root_checks/h1_v2_gls_audit.json`.
At rounded lambda0.9158, correct same-covariance intercept R² is0.46905,
not0.99916; coefficients agree. Unbiased fixed-covariance sigma² is0.00105894
versus ML0.00102767. This is an audit of frozen v2, not a replacement analysis.

H2 admission also uses `bin_to_species.get(b, b)` (unknown bins become fictitious
species) while emitted members use `.get(b, '')`. Reconcile admission and analysis
missingness. Explicitly define these as structural clusters and state whether
`04_cluster/clusters.tsv` contains only 90%-sequence-cluster representatives or
all gene instances. Demonstrate the membership basis from upstream provenance;
either expand correctly or condition the population on the selected members.

### H3 v2 code checks (before final dossier)

`06_match_gem_and_search.py` samples whole GEM proteomes until reaching the SMAG
protein count, allowing an overshoot. It then only DESCRIBES sequence lengths,
exact duplicate counts and overlap; it does not actually match length strata or
deduplicate. Do not label this length/redundancy matched. A concrete repair is
canonicalized unique sequences, explicit cross-release overlap policy, exact
length-stratum quotas/common support, deterministic/repeated samples, and fresh
searches of each constructed database. Distinguish a feasibility/convenience
subset from inference for either entire catalogue. Keep caller and species
composition restrictions visible; one pair of reference subsets is not random
replication of whole catalogues.

## Reproducibility checks across dossiers

`INPUTS.tsv` entries marked DYNAMIC are discovery templates, not pinned inputs.
The final executed package must enumerate and checksum the actual per-sample FAA,
Bakta/eggNOG/contig maps, PDB models and reference files read by each run (or point
to an already verified immutable manifest with matching checksums). Keep absence
and denied access explicit. A hash of a summary does not pin raw inputs. Run these
I/O audits under Slurm when substantial; retain every failed job and version.

H3 v2 archive-integrity flaw: `05_extract_and_genecall_smag.py` runs partial
`tar -xz ... 2>/dev/null || true`, then declares a FASTA complete if it starts `>`
and ends in newline. A truncated archive member can end on a newline while losing
other contigs/sequence. This does NOT establish complete genome extraction.
Use tar header declared member size and read exactly that many bytes from the
verified ordered part stream; reject incomplete members even if FASTA is syntactically
valid. Retain tar-member order, declared/observed bytes, complete-member status and
SHA-256s, plus the actual extracted source genomes used (no cleanup deletion).
The current code deletes extracted genomes and translated proteins without source
content hashes; rebuild provenance and rerun if any incomplete member was included.

### H2 v2 code checks (during repair; verify final version)

The revised mapping pairs raw-bin and Bakta contigs by FASTA order and length only.
Identical lengths/order are not proof of nucleotide identity. Verify actual sequence
hashes and strand orientation, including reverse-complement handling or fail-closed
exclusion. Preserve one-to-one mapping/ambiguity records and actual input hashes.
The revised positional null draws `randint(0,n_cds-1)`, admitting edge-truncated
windows, while the observed primary population requires complete four-flank windows.
Use identical callable-window/annotation-source gates for observed and null positions,
otherwise the permutation creates an artificial recurrence deficit. Retain actual
eligible position counts and don't classify unmappable annotations as observed noKO.

H2 v2 callability remains inconsistent in `03_analyse_conservation.py`: source-file
unavailable and coordinate-mismatch flanks set `has_any_callable=True` and leave
`is_complete=True`. Thus missing annotation evidence can still become a negative.
The script adds KO observations to `ko_by_sp_primary` before verifying that the
entire locus qualifies; partial loci can contribute primary recurrence once two
other complete species make the family eligible. Apply locus eligibility FIRST,
then aggregate only eligible species/loci to primary numerator AND denominator.
The dossier calls contig identity 'formally proven' from lengths and repeats a
'candidate novel folds' label despite claiming removal. Verify actual code/data,
not those assurances. A nonconverged adjusted logistic model is not valid evidence
of covariate adjustment; its failed status must govern interpretation.
