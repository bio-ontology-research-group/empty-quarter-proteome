# Family follow-up reproducibility

FAMILY FOLLOW-UP — 15 October 2026

Scope
Root implemented and ran F3, then F1/F4, with independent Gemini 4 Argon High interpretation and manuscript reviews. F5 remains archived after its previous negative novelty assessment. H1–H4 and the previous family-discovery bundle remain separate, retained inputs.

Locations
Authoritative IBEX root (N):
/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/family_followup_20261015_clean
Previous input root (D):
/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/family_discovery_20261015_clean
Workbench (W): /ibex/scratch/projects/c2014/rob/metagenomics
Local evidence: reproducibility/darkmatter/family_followup_20261015
Reviews: reviews/darkmatter-family-followup-20261015
GPU outputs: /mnt/data/u/hohndor/eq-family-followup-20261015 on ORIX.
The original ProstT5 runner/cache remains in /mnt/data/u/hohndor/eq-family-discovery-20261015/gpu on ORIX.

Sources and denominators
The previous bundle documents all 5,438 representative genomes and 20,451,986 Prodigal translations, raw profile searches, sampled external protein panels and source/coordinate joins. The 19,186,853 Bakta CDS denominator elsewhere in the paper is a different caller-defined count.

SMAG taxonomy comes from the original Extended Data Table 2 at Caiyulu-818/SMAG. inputs/smag_extended_table2.xlsx and sheet1.tsv retain the complete 40,039-record parse. GEM taxonomy comes from portal.nersc.gov/GEM/genomes/genome_metadata.tsv. Its column labelled ecosystem contains GTDB taxonomy; the code checks the d__/p__ syntax. Explicit phylum-name aliases reconcile nomenclature, without claiming reclassification against a common GTDB release.

The expanded HeR panel comes from Bulzu et al. (2021), DOI 10.1128/mSphere.00661-21, dataset 10.6084/m9.figshare.13286486.v3. All nine deposit files are retained in inputs/bulzu2021, with publisher MD5/size and local SHA256 checks. The ungapped phylogeny alignment contains 1,699 records, including type-1 outgroups; 1,566 pass the HeR profile and aligned-lysine gates. Alignment-derived sequences may lack terminal residues, and their gene completeness is unavailable. The original 25-record curated panel contains annotated homologues as well as experimental anchors. This is an explicit published panel, not an exhaustive current database.

The structure reference is RCSB 6SU3, chain X. Its original PDB and metadata are retained. The associated paper's title, authors and DOI were verified through RCSB and PubMed. Automated full-text retrieval from PMC, Europe PMC and the institutional repository was blocked; literature/kovalev2020_unavailable.json records the attempts. The failed responses were not accepted as papers. The other family literature and all Bulzu deposit files are retained locally and on IBEX.

ESMFold used checkpoint 75a3841ee059df2bf4d56688166c8fb459ddd97a. The immutable checkpoint identifier, input identities, executed runner and GPU provenance are retained; ProstT5 additionally retains its model-cache file hashes. The 175 unique predictions restore 178 source records. Protein translation, PDB residue identity and confidence-scale checks are explicit.

Ecology uses the canonical 551-sample metadata and current full CoverM profiles with exact 5,438-MAG closure. F5Dr1 has 100% unmapped reads and undefined catalogue fractions, leaving 550 usable profiles. Pairing gives 69 surface/deep visits at 37 sites. profile_inputs.tsv records every input path/hash.

Reads come from original source-assembly sorted BAMs under /ibex/scratch/projects/c2014/EmptyQuarter_Data/soil/wgs_assemblies/binning/<sample>/. alignment_sources.tsv records paths, sizes, mtimes and header/index hashes. Fifty regional BAMs and indices are retained on IBEX with SHA256 identities; their filtered depth output exactly reproduces the original whole-BAM output. Whole-BAM hashes were not recomputed. These are original assembly mappings, not independent remapping, assembly or expression evidence.

Accepted executions
IBEX 53825139: metadata and input preparation; supersedes an XML namespace parser failure.
IBEX 53825414: taxonomy/core comparison and HMM-match-column FastTree.
IBEX 53825474: finer-rank comparison.
IBEX 53825413: complete A/B neighbours, actual contig lengths and translations; supersedes strand parsing and incomplete all-source site-table failures.
IBEX 53825793: exact ecology decomposition; supersedes correctly rejected non-finite values in the zero-mapped sample.
ORIX 317013: four-sequence ESMFold pilot. ORIX 317016: all 175 unique sequences, four recycles, chunk size 64, H200.
IBEX 53825794: TM-align and structure/translation checks; supersedes an output-label parser mismatch.
IBEX 53825833: download and verify all Bulzu deposit files.
IBEX 53825969: expanded published homology search.
IBEX 53826156: expanded full/core reference clustering.
IBEX 53826077: read depth and fusion junctions, including regional BAM replay.
ORIX 317043: expanded ProstT5, 1,744 source records by 1,024 dimensions, duplicate-equivalent encodings; independent verification PASSED.
IBEX 53826408: initial expanded-panel analysis/plots.
IBEX 53826635: final projection sensitivity with random initialisations. PCA initialisation alone made the earlier seed variations identical; the primary PCA-initialised projection and all numerical biological results are unchanged.
IBEX 53826474: independent depth replay, junction recount, algebra/count and analytical rarefaction audit. Any later closure job is identified in results/audit/summary.json.
IBEX 53825970: earlier small-panel figures, retained as historical sensitivities.
IBEX 53826483 and 53826636: failed paper-download attempts; no scientific result depends on their success.

Re-running
Use a fresh output directory and change N in scripts/common.py to that directory. Preserve D/W paths or restore their frozen input trees. Submit each .sbatch through Slurm. CPU jobs use account pi-hohndor and batch/debug with the gspa-neural Python environment. No analysis runs on a login node.

Dependency order:
prepare -> f3_comparison -> f3_finer;
f4_integrity -> ecology;
published_helio -> published_clusters -> expanded embeddings -> expanded_analysis;
read_support and structure_audit must finish before expanded_analysis;
audit last.
The source commands, versions and outputs are retained in scripts, summaries and tool logs. Existing inputs are reused at their frozen paths, not downloaded silently from a new release.

The retained GPU run_embed.sbatch accepts INPUT_FASTA and OUTPUT_DIR. ProstT5 uses AA2fold, the legacy tokenizer, residue-only mean pooling, UZOB-to-X conversion and no truncation. Identical cleaned sequences are computed once, then restored to all IDs. expanded_embeddings contains model/cache and output hashes, the sequence manifest, provenance and independent verification. ESMFold inputs, predictions, confidence arrays and provenance are in inputs and structures.

Clustering uses fixed joint MMseqs partitions at 50%/80% identity with 80% bidirectional coverage. Exact sequences are deduplicated within source before equal-size draws. There are 1,000 draws with seed 20261015+j. Their 2.5–97.5 percentiles are conditional sampling ranges, not population confidence intervals. Independent hypergeometric expectations check the Monte Carlo estimates. Core trimming removes residues outside the first/last occupied HMM match positions and preserves internal insertions.

Ecology reports two separate estimands. The original contrast uses replicate and site medians. The exact B+W+R decomposition uses replicate means, then means over visits within site and equal weights across sites. Genus-unclassified MAGs are separate groups. Bootstrap intervals use 2,000 site resamples with seed 20261015. All six displayed median contrasts receive Holm correction. Rank-specific components are descriptive; fixed representative traits cannot measure within-strain adaptation.

Local manuscript reproduction
verify_family.py registers explicit manuscript values and extracts them from pinned compact evidence. Its sources are included in current_manuscript_artifacts.tsv. plot_family.py rebuilds the manuscript figure from retained numerical coordinates and tables, without rerunning t-SNE, embeddings or predictions. The heavy rerun remains the Slurm workflow above. The repository's run_all.sh runs the compact claim checks, figure regeneration and manuscript build.

Interpretation limits
The supported claims concern sequence coverage within known families, predicted fold compatibility and support in original read mappings. They do not establish biochemical activity, expression, desert endemism, a new fold, universal reference absence or absence from ProstT5 training data. F3 is not an assertion that these recognizable family members belong to the UniRef50-dark cohort. F1 is an operational CPF panel, not an exhaustive photolyase census. F4 neighbourhoods are not experimentally established operons. This bounded follow-up does not repair the manuscript's previously documented raw-data access and deposition prerequisites.

Review-prompted checks
IBEX 53826937 independently re-clustered only the 34 environmental Saccharimonadia sequences, excluding every curated and published reference. Per ten sequences, mean core cluster counts are 8.04 EQ versus 4 GEM at 50% identity and 9.605 versus 10 at 80%. This is a separate partition, not the earlier partition with reference clusters merely unsubtracted.
The same job inspected the seven EQ and one GEM loci whose A/B assignments share an ORF. Their A/B domain matches do not overlap and each covers at least 70% of its profile. Corrected marker-order labels are retained in results/review_checks; the original tie-sorted A/B labels do not establish physical gene order.
IBEX 53827015 passed twelve audit groups after these checks and the random-initialisation projection rerun. Its manifest closed 1,124 retained files. The final snapshot closure, if later, is identified by results/audit/summary.json.
