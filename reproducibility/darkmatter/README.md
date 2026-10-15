> Current state (14 October 2026): H1–H4 integrated, independently audited, and
> refined through two post-integration Gemini 4 Argon High/root review rounds.
> Both local drivers pass; 152 current extension claims use 27 pinned inputs.
> See `../../reviews/darkmatter-2026-10-13/round7/FINAL_REPORT.md`. Dated sections
> below preserve earlier states; the round7 record supersedes conflicting summaries.

# Current P4 data and code audit

Started 13 October 2026 for the approved three-round review. This directory
documents the current September/October manuscript. The full registry combines
historical claims with current H1–H4 claims; the current-paper driver additionally
checks the intervening review refinements and rebuilds the current figures.

## Evidence levels

`sample_input_manifest.tsv` inventories the paths in the 551-row IBEX sample
metadata table. It records existence, byte counts and resolved symlink targets,
not content verification. Its 2,204 entries cover assemblies and processed paired
and orphan reads. Raw sequencer files and public accessions require a separate
join; processed reads must not be described as raw reads.

`inventory_remote.py` runs under Slurm and inventories current analysis trees,
hashing files up to 16 MiB and explicitly labelling larger files metadata-only.
IBEX job 53675321, work directory:
`/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p4_review_20261013/`.
This initial sweep is not a complete raw-data checksum audit.

`source_snapshot/` is a local, untracked evidence cache copied from
`/ibex/scratch/projects/c2014/rob/metagenomics/`. It retains original scripts,
README files and selected text/tabular results. The text transfer limits each file
to 16 MiB; absent files are not evidence that the original does not exist.
The broad BGC text transfer was interrupted after it entered many intermediate
files; that subtree is partial. Original source files have not been changed.
Scripts retain their original paths and environment requirements; a copied script
is not yet a portable reproduction workflow. No computation is claimed solely
because a result table or README was copied.

## Claim-to-analysis map

| Current claim | IBEX directory under the root above | Key dependencies |
| --- | --- | --- |
| Assemblies, refined bins and species representatives | `05_refine`, metadata; taxonomy R232 | Processed reads, assemblies, CheckM2, dRep membership, GTDB-Tk |
| Genome-resolved sequence-dark census | `04_protein_catalog/darkmatter_clean` | 552 assembly-file universe, Pyrodigal, cluster-90 membership, UniRef50, unrefined HQ/MQ bin membership |
| Dark families, recurrence, selection and shifted-frame contrast | `06_darkmatter/darkfam_basis_clean` | Genome-resolved dark representatives, family membership, CDS, member-aware CoverM, metadata |
| Cultured/MAG same-species comparisons | `06_darkmatter/cultured_comparison_clean` | 722 culture assemblies/proteomes, skani, UniRef50 and reciprocal protein hits |
| Lineage darkness, Thermoplasmatota and matched comparisons | `06_darkmatter/lineage_darkness_clean` | Native contig keys, cluster darkness, R232 taxonomy, abundance, quality and annotation tables |
| Structure-screen calibration | `04_structure/00_structure_clean` | Model coordinates, control identities, AFDB/PDB versions and coordinate-search hits |
| Sequence-domain support | `04_structure/pfam_candidates_clean` | Candidate members, other confident models, controls, Pfam/HMMER |
| Reference coverage and rarefaction | `06_darkmatter/specificity_clean` | Identical dark/control queries, catalogue versions, hit sets and rarefaction scripts |
| PICRUSt2 and annotation ceilings | `06_reference_bias/refbias_clean`, `04_functional/repko_clean` | Rep annotations, CoverM, pinned amplicon predictions |
| Dark-family environment associations | `06_darkmatter/darkfam_ecology_clean` | Family TPM/presence, chemistry/climate matrix, repeated-site design |
| Controls and contaminant evidence | `06_darkmatter/contamination_clean` | Control reads, mapping settings, DNA concentrations |
| Annotated environmental traits | `06_metabolism/{differential_function,uv_survival,osmolyte,bgc,resistome}_clean` | Rep annotations, modules/panels, abundance, metadata |

## Outstanding reproducibility work

Resolve raw-to-filtered provenance and accession mapping; validate missing or empty
sample inputs; identify all external database snapshots and environments; pin
claim-generating scripts and data; distinguish superseded from active results;
build the missing figures; replace manuscript-era mismatches in the claim registry
with checks of current evidence. Full remote re-execution and publicly accessible
deposition are separate from local table verification.

## Raw-read and ENA checks (13 October 2026)

`reconcile_ena.py` reconciles all 554 runs (551 samples) to successful private
submission/modify receipts in `ena_evidence/`; `ena_run_map.tsv` provides accessions.
A fresh public ENA filereport returns only its header (zero public runs).
Submission manifests confirm 333 NovaSeq X and 221 NovaSeq 6000 runs; these are
run counts, not independently sampled soils.

`raw_input_manifest.tsv` checks the 1,108 paired raw-file paths in the submission
table: 728 match the recorded byte size, 302 fail traversal with PermissionError,
and 78 are absent. The absent paths belong to the 39-run December-30 sequencing
folder; its ENA staging symlinks also target the absent location. The permission
failures concern the 2024 lane directories. Robert has been asked for the moved
folder location and read/traverse access for hohndor. No permission changes or
raw-data rewrites were attempted. File sizes do not substitute for checksums.

## Round-1 completion and additional inventories

The first inventory job completed successfully: 447,872 records, including
446,361 SHA-256 records, 1,497 metadata-only large files and 14 broken temporary
symlinks. `ibex_files.tsv.gz` preserves the inventory. A second job, 53676339,
recorded 443,134 upstream metadata-only entries under protein catalogues, refined
genomes, databases, abundance and functional annotation. Its output is
`ibex_upstream_files.tsv.gz`; the sweeps overlap and are not a unique-file census.

Frame-stop audit job 53676424 completed successfully. `frame_stops/` contains
per-family results and hashes of the 10,000 input FASTA files. Of 50,898 sequences,
43,953 have terminal stops and 21 have internal stops. The review's suggestion
of widespread internal-stop contamination is unsupported by this audit.

Run `python3 reproducibility/darkmatter/check_current_tables.py` and
`python3 reproducibility/darkmatter/round1_sensitivity.py` from the repository root
to verify the frozen compact inputs and reproduce the round-1 table calculations.
These do not rerun the upstream searches. See
`reviews/darkmatter-2026-10-13/round1/response.md` for the full disposition and
remaining analyses. Work is at the agreed journal-selection checkpoint before
round 2; all three reviews are not yet complete.

## Nature Communications revision after round 1

Robert selected Nature Communications and authorised rounds 2 and 3. New IBEX
jobs 53678755 and 53678768 completed successfully: family FASTA covariates and
culture/MAG nucleotide comparisons. Job scripts, input checksums and results are
in this directory. `selection_covariates/family_covariates.tsv.gz` permits the
matching calculation without an IBEX connection. `culture_nucleotide/` retains
per-gene results and the input checksum ledger; alignment files remain on IBEX.

`bash reproducibility/darkmatter/run_current.sh` reruns the current compact
arithmetic, matching sensitivity and Figures 1--5, then builds the manuscript.
Install the versions in `requirements.txt` in a separate Python environment;
LaTeX/latexmk/BibTeX are also required. This driver uses the public manuscript
and works without the private paper submodule. It does not run the historical
figure driver or claim that every upstream analysis is reproduced. On 14 October 2026, both this driver and the full
`reproducibility/run_all.sh` passed after stale manuscript gates were updated.
The full driver checks 62 COMPUTE and 216 FROZEN claims; the current driver also
checks 152 H1–H4 registered claims from 27 checksum-pinned files. These are local
checks of available evidence, not a complete raw-to-final rerun. The raw access,
deposition and historical provenance gaps remain open.

## Round-2 revision and final-review package

The preceding checkpoint text is historical. Nature Communications was selected,
round 2 is complete, and round 3 is the final authorised review. The public source
and optional private submodule are kept identical; submission is not authorised.

New current outputs:

| Evidence | Script/output | IBEX job |
|---|---|---|
| Explicit nucleotide base alignments | `culture_nucleotide_check.py`, `culture_nucleotide_cigar/` | 53679350 |
| Culture Bakta locus diagnostic | `culture_locus_audit.py`, `culture_locus/` | 53680202 |
| AntiFam v8.0 screen | `antifam_screen.sbatch`, `antifam/` | 53678944 |
| Opposite-strand overlap | `shadow_overlap.py`, `shadow_overlap/` | 53680096 |
| Stratum bootstrap and raw structure audit | `selection_estimands.py`, `structural_raw_audit.py` | 53679366 |
| Site aggregate sensitivity | `site_level_sensitivity.py` | 53679245 |
| Observed environment manifests | `capture_environments.sbatch`, `environments/` | 53679246 |

**Superseded nucleotide run:** `culture_nucleotide/` and round2 submitted copies
used approximate PAF identities without `-c`. Use `culture_nucleotide_cigar/` for
all current numbers. Original results are retained for audit. Culture locus phase
is a first-aligned-base diagnostic, not full codon correspondence. The overlap
audit only compares existing calls, so its low flag rate does not validate all ORFs.
The first GFF-based overlap attempt failed on truncated gzip; the successful retry
used original Pyrodigal protein-header coordinates and verified all536,437queries.

`run_current.sh` now reproduces matching, bootstrap intervals, overlap sensitivity,
35 current revision claims, the environmental matrix and Figures1–5, then builds
the manuscript. `review_inputs.sha256.tsv` pins59retainedinputs;
`upstream_code_manifest.tsv` pins162originalsourcefiles. These scripts include
absolute historical paths and are preserved evidence, not a portable upstream
workflow. `requirements.txt` pins the local compact environment. Remote conda
manifests describe audit-time installations; the inaccessible python-utils
installation is explicitly recorded as PermissionError. Current environment
capture does not prove which build was used historically.

The environmental matrix rebuilds from `environment_inputs/` to
`environment_rebuild/env_matrix.tsv`, SHA256
`581d602009f36a03a8409903a0e51802b9f1728e65bb387a6b222b2c186be9eb`.
`climate_retrieval_source.py` is the original companion source for provenance,
not automatically invoked: re-querying the provider would not reproduce the
unarchived raw response/model selection. Laboratory XRF units, instrument and
non-detect semantics still need confirmation. Handheld field XRF is a separate
source, not included in this matrix.

`current_claims.tsv` is a scoped revision ledger generated from retained inputs.
It does not replace the165historicalclaims or silently relax their gates.
Full registry rebinding and raw-to-final upstream reproduction remain open;
`notes/CORRECTIONS.md` records the current values and reasons.

## Final state after round 3

All three Opus5.5reviews completed and root implemented the final corrections.
The current ledger contains49revisionclaims, with66pinnedreviewinputs and162
preserved upstream source files. Current full-locus outputs are in
`culture_locus_all/`; the earlier unshared-only audit remains in `culture_locus/`.
Job53680886also produced `culture_family_membership/`:65of1,341eligible dark calls
belong to47admitted recurrentfamilies. These populations are kept separate in the paper.

`produce_culture_diagnostics.py` generates per-pair and post hoc summaries;
`current_revision_checks.py` verifies without writing unless its explicit
`--write-ledger` curator option is supplied. `overlap_sensitivity.py` now also
reports the flagged families' own distributions, showing why dN/dS below one
cannot by itself establish coding status.

A curated copy of this package, the manuscript, review reports and retained
literature is placed on IBEX at
`/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p4_review_20261013/review_package/`.
It excludes the18GBlocal source cache and raw uncompressed inventory duplicates;
its preserved compact inputs are checksum verified. The original large upstream
data remain at the paths documented in the inventories. This copy is a research
working package, not a public archived release with a DOI.

## Completed large-analysis checksum follow-up

IBEX job53680819 completed in12m27s with zero errors: SHA256 for1,494unique
large analysis paths,292,813,193,914bytes. Their1,497inventory references are
filled in `ibex_files_complete.tsv.gz`. This combined clean-analysis inventory
has447,858checksum-bearing records and14broken temporary links; link duplicates
are present. The original inventory is preserved. The separate5.62TBupstream
inventory remains metadata-only except for specifically pinned analysis inputs;
raw-read provenance/access gaps are unchanged. Hashing existing outputs does
not demonstrate that their historical computation was correct.

Final remote verification job53681298 completed successfully:228input/source
checksums, all49revisionchecks, and final manuscript/figure/script hashes match
the local manifest. See `reviews/darkmatter-2026-10-13/round3/ibex_package_verification.log`.

## Round 5: exploratory extension package

`round5_extensions/` retains three new post hoc analyses, executed on IBEX before
Gemini 4 Argon's fifth review: controlled gene-caller substitution in 22 culture
assemblies, an expedition-label null over the full family-presence matrix, and
paired functional-profile nulls with a fixed abundance-threshold grid. Designs,
input/output hashes, Slurm scripts, detailed results and scope are documented there.
The compact driver now includes its 159 checks, alongside 49 prior and 40 round-4
checks and 239 original input/source hashes. Sixteen new manuscript-facing values
have a separate reviewed claim ledger; this does not resolve historical gates.
The remote scripts regenerate heavy analyses; the local checks validate compact
outputs. Review receipts and a ranked further-hypothesis plan are under
`reviews/darkmatter-2026-10-13/round5/`.

## H1–H4 integration and two review rounds (14 October 2026)

The extension execution map is `round7_analyses/README.md`; compact pinned inputs
and extraction helpers are in `round7_final/`. `check_round7_registered.py` writes
`round7_registered_claims.tsv`. The new Figure 4 uses the pooled qTM endpoint and
is regenerated after historical figure scripts so that withdrawn alignment-TM
subdivisions cannot overwrite it. Review receipts, root responses, independent
raw audit code and completed scheduler records are retained in
`reviews/darkmatter-2026-10-13/round7/`.
