# ENA deposit — RAK-WGS (PRJEB125153): reads, primary assemblies, bins, MAGs

Scripts `code/ena_deposit/` (README there, one line per script). Work dir, not in git:
`/ibex/scratch/projects/c2014/alelopezv/ena_wgs/` — `xml/` + `receipts/` (what was sent, what ENA answered), `logs/*.tsv`
(status, md5s, accessions), `cache/samples/` (soil BioSamples before the cleanup MODIFY), `assemblies/` and `bins/`
(staging: gz, manifests, coverage, logs), `webin-cli-10.0.0.jar`. Credential only `~/.webin_pass`. Test service
(wwwdev) = daily mirror of production; a MODIFY of production records cannot be rehearsed there.

## Studies and release
| study | accession | hold | note |
|---|---|---|---|
| RAK-16S (amplicon) | PRJEB106069 | 2026-11-25 | owns 445 of the 551 soil BioSamples |
| RAK-PR (cultures) | PRJEB111735 | 2026-11-25 | |
| **RAK-WGS (this deposit)** | **PRJEB125153** | **2026-11-25** | title "Rub' al-Khali Soil Shotgun Metagenomes" |
| RAK umbrella | not yet | | PRJEB104209 cannot be converted; new umbrella "Rub' al-Khali Desert Microbiome" goes LAST, PI |

## Status (what is on ENA)
| object | count | state | accessions |
|---|---|---|---|
| soil BioSamples | 551 (445 shared with RAK-16S, 106 new) | registered, ERC000022 | `wgs_biosample_map.tsv` |
| reads | 554 runs, 6.60 TB | archived | `wgs_runs.tsv` (ERX/ERR) |
| primary metagenome assemblies | 549 | submitted, private until release | `wgs_assemblies.tsv`, ERZ30552864–ERZ30555284 (final) |
| binned metagenomes | 12,434 samples + assemblies | submitted, private until release | `wgs_bins.tsv`, ERZ30573912–ERZ30593965 (final) |
| MAGs | 5,437 samples + assemblies | submitted, private until release; GCA assigned for 5,388 (49 awaiting ENA's next accession email) | `wgs_bins.tsv`: ERZ30595485–ERZ30606248, GCA_990079656–GCA_990090661 + contig / chromosome accessions (`gca_mag`, `sequence_acc_mag`) |
| not deposited | 1,992 LQ bins; QC-rerun reads; protein catalogue, dark families, ESMFold models → Zenodo (later) | | |

## Decisions
### Soil samples and reads
| item | decision |
|---|---|
| study naming | no short names inside titles; "Rub' al-Khali" hyphenated |
| BioSample | checklist ERC000022; title = bare id; description "Rub' al-Khali desert soil, site N, <deep\|surface\|root-adjacent> compartment, expedition N" (62 = site road (location 1), 64 = site campground); project name "Rub' al-Khali Desert Microbiome"; local context = field-sheet `feature`; coordinates 4 decimals; "nucleic acid amplification" dropped, target gene + primers kept; site 52 coordinates untouched (ENA correct, ecology-data wrong: packet below) |
| reads | RAW (Rob); public names `<sample>_wgs_R{1,2}.fastq.gz`, dual-batch `_wgs_2024` / `_2026`; 48Dr2 = processed pair (raw R1 truncated); instrument from header serial (A01018 NovaSeq 6000: 221, LL00134 NovaSeq X: 333); kit from index codes (QIAseq FX UDI 403 + 70 assumed, NEBNext Ultra II 81); insert size measured on the binning BAMs (547); description "Shotgun metagenome of soil sample ID" |

### Primary assemblies
| item | decision |
|---|---|
| objects | 549 = 551 − 51Sr1 (no contigs) − F5Dr1 (empty); MEGAHIT 1.2.9, contigs ≥ 500 bp, deposited unfiltered |
| source | organized tree `EmptyQuarter_Data/processed/sequencing-results/sites/<trip>/site<N>/<s>/wgs/assembly/<s>_contigs.fa` (link to `soil/wgs_assemblies/assemblies/<s>/`) |
| manifest | `ASSEMBLY_TYPE primary metagenome`; SAMPLE = soil BioSample; RUN_REF = all runs of the sample; ASSEMBLYNAME `<s>_assembly`; PROGRAM "MEGAHIT v1.2.9"; PLATFORM ILLUMINA; MOLECULETYPE genomic DNA; DESCRIPTION "Primary metagenome assembly of Rub' al-Khali desert soil sample <s>, contigs >= 500 bp" |
| COVERAGE | length-weighted mean depth over all contigs, `samtools coverage` 1.16.1 on the SemiBin2 self-supervised BAM (BWA-MEM 0.7.17, own reads on own assembly) |

### Bins and MAGs
| item | decision |
|---|---|
| sets | 12,434 HQ+MQ refined bins (CheckM2 completeness ≥ 50, contamination ≤ 10) all as `binned metagenome`, THEN the 5,437 dRep species representatives with contamination < 10 again as `Metagenome-Assembled Genome (MAG)` (ENA's two-level path; `48Dr2__self_SemiBin_25` at exactly 10 stays binned only: checklist "if contamination ≥ 10 % submit as a binned assembly"); 1,992 LQ not deposited |
| names | internal bin names kept (traceability to the evidence tables): `<sample>__single|self_SemiBin_N`, `single` = SemiBin2 pretrained global-model run, `self` = self-supervised run (both single-sample), `_sub` = DAS Tool contig subset; MAG record = same name + `_MAG` |
| taxonomy | GTDB-Tk 2.7.2 / r232: reps `M/05_refine/_gtdbtk_r232/`, non-reps `M/05_refine/_gtdbtk_r232_nonreps/`. NCBI taxon from the GTDB lineage (script 15, every name checked submittable at ENA): species → "Genus sp." → family/order/class/phylum as "<GTDB name> bacterium|archaeon" or the NCBI organism name of that rank's genomes (r232 metadata) → "uncultured bacterium|archaeon"; step recorded in `ncbi_level`. No taxon registration: MAGs take organism-level environmental names |
| rRNA/tRNA | barrnap 0.9 + tRNAscan-SE 2.0.12, `M/mimag_rrna_trna_clean/` (reps + non-reps); MIMAG tier from completeness, contamination, 5S/16S/23S and ≥ 18 tRNAs |
| binned sample (ERC000050) | alias = bin name; taxon as above; description "This sample represents a metagenomic bin from the metagenomic sample ERS<soil>. GTDB release 232: <lineage>."; project name, metagenomic source "soil metagenome", sample derived from = soil ERS, isolation_source "Rub' al-Khali desert soil", collection date / coordinates / country / ENVO contexts copied from the soil sample; binning software "SemiBin2 2.2.1; DAS Tool 1.1.7" + parameters; assembly software "MEGAHIT;1.2.9;<params>"; completeness / contamination (CheckM2 1.1.0; "100" never "100.0"); taxonomic identity marker "multi-marker approach (GTDB-Tk 2.7.2, 120 bacterial / 53 archaeal marker genes, GTDB release 232)"; assembly quality (exact checklist strings); 16s recovered + barrnap; number of standard tRNAs extracted (total) + tRNAscan-SE. No depth, elevation, investigation type |
| MAG sample (ERC000047) | identical attributes; sample derived from = the binned ERS; description "This sample represents a MAG derived from the metagenomic sample ERS<soil> through binned sample ERS<bin>. GTDB release 232: <lineage>."; alias `<bin>_MAG` |
| manifests | as primaries with the tier's type and sample; coverage = length-weighted mean depth over the bin's contigs; single-contig MAG (`11PRr3__self_SemiBin_25287`, 973 kb, 97.5 %) at chromosome level with a chromosome list (name "1", Linear-Chromosome, circularity not asserted) |

## ENA constraints the deposit follows
| topic | rule |
|---|---|
| routes | assemblies only via Webin-CLI genome context, one manifest + gz per object, Webin-CLI uploads the file itself (nothing pre-uploaded); samples by XML to the drop-box (`ena_common.py`), batches of 200; reads with pre-uploaded files by XML (08), because Webin-CLI `-submit` re-uploads |
| accessions | ERZ is final for primary and binned metagenomes; MAGs and isolates get GCA + contig accessions; cite PRJEB + SAMEA (+ GCA), never ERZ |
| accession authority | ERZ: ENA's live analysis list (`report/analyses`), folded by `10/18 --fold`; the Webin-CLI status logs are secondary. GCA + sequence accessions: ENA's "genome assembly has been accessioned" emails, saved verbatim as `ena_wgs/receipts/emails/ena_email_<date>.txt` and parsed by `18 --tier mag --fold` |
| concurrency | ≤ 12 Webin-CLI submissions at once (ENA FTP session limit); arrays of 25 manifests per task (user limit 2,000 jobs); retries as a targeted array over the manifests without an accession |
| test service | daily mirror; used for sample registration and validation rehearsal, not for MODIFY of production records |

## Pending decisions
- Soil BioSample taxon: 256318 "metagenome" (current) vs 410658 "soil metagenome" (ENA's preference); one MODIFY pass over the 551, shared with RAK-16S — hers.
- BioSample depth 40 mm vs "5–10 cm" in the Methods; Methods also to state two instruments, the 2×101 Nov-2025 batch, kit per batch, measured insert size — hers.
- Umbrella study and aligned release of the three studies — PI.
- Paper Data Availability text (PRJEB125153: reads, assemblies, bins, MAGs; release 2026-11-25) and the supplementary bin table (name, sample, run key, quality, ERS, ERZ/GCA) — hers.
- RAK-PR §5.2 cites PRJEB125153.

## Housekeeping after ENA processing
Delete the gz staging (`ena_wgs/assemblies/fasta/` 165 GB, `ena_wgs/bins/fasta/` 15 GB) and the root copies of the 46 Webin-CLI read pairs in the upload area (~0.5 TB) once the Analyses Report shows everything processed and the MAG GCA emails have arrived.

## Inputs
Raw reads: four sequencing batches under `EmptyQuarter_Data/soil/wgs_reads/raw/` (batch, flowcell, instrument, kit per run in
`wgs_runs.tsv`; 22 non-cohort tokens — controls, QC reruns, typo ids — in `ena_wgs/logs/raw_non_cohort.tsv`). Processed
reads: `wgs_assemblies/processed_reads/` (R1/R2/SE). Assemblies and bins: see the decisions tables.

## Commands
```
cd /ibex/scratch/projects/c2014/alelopezv/repos/empty-quarter-proteome/code/ena_deposit
# reads
python3 00_build_runs_table.py --symlinks && python3 04_generate_manifests.py --study PRJEB125153
MODE=validate sbatch --array=1-554%16 05_webin_cli.sbatch;  python3 08_submit_runs_xml.py;  python3 09_fold_accessions.py
# primary assemblies
python3 10_build_assembly_table.py --prepare;  sbatch --array=1-549%64 11_gzip_assemblies.sbatch;  sbatch --array=1-549%64 12_coverage_bam.sbatch
python3 10_build_assembly_table.py --study PRJEB125153
MODE=validate sbatch --array=1-549%32 13_webin_cli_genome.sbatch;  MODE=submit sbatch --array=1-549%12 13_webin_cli_genome.sbatch
python3 10_build_assembly_table.py --fold
# bins (tier = binned, then mag)
srun -c 2 --mem 16G python3 14_build_bin_table.py;  srun -c 2 --mem 16G python3 15_resolve_bin_taxa.py;  sbatch --array=1-25 17_gzip_bins.sbatch
python3 16_register_bin_samples.py --tier <tier> --service test --only <bins> --submit     # rehearsal
python3 16_register_bin_samples.py --tier <tier> --service prod --submit
python3 18_build_bin_manifests.py --tier <tier> --study PRJEB125153
W=/ibex/scratch/projects/c2014/alelopezv/ena_wgs/bins
CHUNK=25 MANIFEST_DIR=$W/manifests_<tier> LOGDIR=$W/logs MODE=validate sbatch -t 02:00:00 --array=1-N%64 13_webin_cli_genome.sbatch
CHUNK=25 MANIFEST_DIR=$W/manifests_<tier> LOGDIR=$W/logs MODE=submit   sbatch -t 04:00:00 --array=1-N%12 13_webin_cli_genome.sbatch
python3 18_build_bin_manifests.py --tier <tier> --fold          # reads ENA's live list (logs/ena_live_analyses.tsv) first
```
Submit steps are irreversible and run on her word. Retries: a targeted array over the manifests without an OK record, never the full range.
