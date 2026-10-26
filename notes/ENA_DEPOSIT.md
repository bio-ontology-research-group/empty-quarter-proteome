# ENA deposit — 551 soil shotgun read sets (RAK-WGS = PRJEB125153)

Scripts `code/ena_deposit/`; work dir (not in git) `/ibex/scratch/projects/c2014/alelopezv/ena_wgs/`: `xml/` + `receipts/`
(what was sent, what ENA answered), `logs/*.tsv` (status, md5s, accessions), `cache/samples/` (BioSamples before the
cleanup MODIFY), `assemblies/` (step 2 staging), `webin-cli-10.0.0.jar`. Reads-phase symlink farm, reads manifests and
per-task Slurm logs deleted 2026-10-25 (regenerable by 00/04). Password only `~/.webin_pass`.
Test service (wwwdev) = daily mirror; MODIFY of production samples cannot be rehearsed there.

## Roadmap under PRJEB125153 (decided 2026-09-03; steps 2–3 revised 2026-10-25)
1. reads (554 runs) — DONE · 2. primary metagenome assemblies (549) — IN PROGRESS, scripts 10–13 · 3. bins: ALL
12,434 HQ+MQ refined bins as `binned metagenome` (ERC000050 samples derived from the soil BioSample), THEN the 5,438 dRep
species representatives again as `Metagenome-Assembled Genome (MAG)` (ERC000047 GSC MIMAGS samples derived from their
binned sample) — ENA's recommended two-level path; LQ (1,992) not deposited · 4. non-ENA datasets (protein catalogue,
dark families, ESMFold models, tables) → Zenodo DOIs · 5. new umbrella (PI) + aligned release dates. Study hold moved to
2026-11-25 (live, Webin report API 2026-10-25).

### Step 2 — primary metagenome assemblies (2026-10-25)
| item | decision / fact |
|---|---|
| objects | 549 = 551 cohort − 51Sr1 (no contigs) − F5Dr1 (empty); MEGAHIT v1.2.9, contigs ≥ 500 bp, headers `<sample>_<N>` |
| source | organized tree `EmptyQuarter_Data/processed/sequencing-results/sites/<trip>/site<N>/<s>/wgs/assembly/<s>_contigs.fa` (symlink to `soil/wgs_assemblies/assemblies/<s>/`); trip/site from workbench `metadata/sample_metadata.tsv` |
| route | Webin-CLI 10.0.0 `-context genome`, `ASSEMBLY_TYPE primary metagenome`, SAMPLE = soil BioSample, RUN_REF = all runs of the sample (3 dual-batch samples carry two) |
| COVERAGE | length-weighted mean depth over all contigs, Σ(depth·len)/Σlen, from `samtools coverage` 1.16.1 on the SemiBin2 self-supervised BAM (BWA-MEM 0.7.17, own reads on own assembly, `soil/wgs_assemblies/binning/<s>/`), the only mapping present for all 549; the pretrained-run minimap2 BAM in the tree covers 545 and gives ~10 % different values, not used |
| name / text | ASSEMBLYNAME `<sample>_assembly`; DESCRIPTION "Primary metagenome assembly of Rub' al-Khali desert soil sample <s>, contigs >= 500 bp" (optional field; RAK-PR used none; carries the contig cut-off, the one fact stored nowhere else); PLATFORM ILLUMINA; MOLECULETYPE genomic DNA |
| 1 M-sequence cap | stated on the genome + MAG doc pages only; dry runs: 25Sr3 (1,004,824 contigs, 9.0.3) and F5PRr3 (5,254,851 contigs, 10.0.0) both validated OK → does not apply to primary metagenomes; all 549 deposited unfiltered |
| staging | gzipped copies in `ena_wgs/assemblies/fasta/` (~165 GB) + `stats/ coverage/ manifests/ logs/ reports/`; kept (her decision 2026-10-25) |
| table | `code/ena_deposit/wgs_assemblies.tsv` (sample, name, BioSample, run_refs, n_contigs, total_bp, coverage, source, ERZ) |
| accession | ERZ is final for primary metagenomes; cite PRJEB125153 + SAMEA, never ERZ |

### Step 3 — bins (12,434 binned metagenome, then 5,438 MAG) — compute started 2026-10-25
| item | decision / fact |
|---|---|
| sets | 12,434 HQ+MQ refined bins = 5,438 dRep reps + 6,995 non-reps in dRep clusters + `25PRr2__self_SemiBin_2_sub` (dropped inside dRep, still HQ+MQ); LQ 1,992 not deposited |
| taxonomy | reps: `M/05_refine/_gtdbtk_r232/`; non-reps: NEW `M/05_refine/_gtdbtk_r232_nonreps/` (same GTDB-Tk 2.7.2 / r232 / options; job 54542503). NCBI taxon per bin from the GTDB lineage: species → NCBI species; else "Genus sp."; else "Family bacterium/archaeon" (all submittable per ENA taxonomy REST) |
| rRNA/tRNA | reps: `M/mimag_rrna_trna_clean/results/`; non-reps: same folder, `SET=nonreps` (array 54542640), domain inherited from the bin's dRep representative |
| bin coverage | length-weighted mean depth over the bin's contigs from the sample's `samtools coverage` table (step 2 output) + DAS_Tool contig2bin |
| samples | binned: ERC000050 by XML (new script 14, helpers of 03), description "This sample represents a metagenomic bin from the metagenomic sample ERSxxx. GTDB r232: <lineage>"; MAG: ERC000047 derived from the binned sample, "This sample represents a MAG derived from the metagenomic sample ERSxxx" |
| order | all 12,434 as `binned metagenome` first, then the 5,438 reps as `Metagenome-Assembled Genome (MAG)` (her decision 2026-10-25; ENA's recommended path); no "representative" clause, the MAG link carries it |
| taxa requests | none for bins; for MAGs only if the test-service dry run rejects the "Genus sp." names |
| taxon resolution, all 12,434 bins (script 15, 2026-10-26) | every name ENA-submittable; step 1 species 352 (276 ENA name + 76 GTDB-metadata majority) · step 2 "Genus sp." 4,373 · step 3 family 7,398 (3,338 GTDB family name + 4,060 NCBI organism name of the family's genomes in the r232 metadata, e.g. "uncultured Candidatus Saccharibacteria bacterium" for Patescibacteriota placeholders), order 289, class 21, phylum 1 · step 4 fallback **0**; `ncbi_level` column records the step (reps alone: 122 / 1,792 / 3,352 / 158 / 13 / 1 / 0) |

## Open items
- **Reads DONE (2026-09-06)** — 554/554 runs registered; next = watch the portal for "archived" on all files (ENA
  emails on errors), then the 549 assemblies (Webin-CLI genome context uploads the FASTAs itself: no pre-upload).
- MAGs paper Data Availability: state PRJEB125153 (reads; assemblies and MAGs to follow under the same study) and
  the release date 2026-11-01 — her edit in `paper/main.tex`.
- History of the reads upload/submission, kept for the assemblies step: **UPLOAD COMPLETE 2026-09-05 06:00** — 1,108/1,108 files OK (first array 51309737: 890 OK; resume pass
  51347541: 95 verified by size, 123 re-sent, 0 refusals). md5 vs the sequencing centre's lists: 342/346 comparable
  files match (48Dr2 = processed pair by design; F2Dr1 pair = list-matching artefact). Validation 51351949 PASSED
  09-05 06:20: 554/554 OK. **Submission 09-06:** Webin-CLI `-submit` (array 51384714) registered **46 runs**
  (`logs/webin_submit_status.tsv`) but RE-UPLOADS every file pair into `webin-cli/reads/<run>/` — it ignores files
  pre-uploaded to the root — so it was cancelled after 46 (≈15 h of duplicate transfer avoided). The remaining
  **508 runs go by XML** (`08_submit_runs_xml.py`: experiment + run per run, files referenced by name + md5 from
  `logs/md5sums.tsv`, aliases `webin-reads-<run>` like Webin-CLI, ENA's live run list checked first, hold
  2026-11-01). Test service OK for 13Dr2_wgs. Log `logs/runs_xml.prod.tsv`. Housekeeping AFTER archiving completes:
  the 46 Webin-CLI runs reference their copies under `webin-cli/reads/<run>/`, so the ROOT copies of those 92 files
  (~0.5 TB, pre-upload) plus partial uploads left by the failed/cancelled Webin-CLI tasks are orphans in the upload
  area — list it and delete what no run references. Do not touch the upload area while ENA is still processing.
  History of the first array (6× 120 s retry on FTP 530): ENA refuses new logins
  ("Access denied: 530") when the account holds too many sessions: a jump 8→16 at 23:04 (09-03) caused 40 min of
  refusals; a stepwise ramp to 16 ran clean for ~30 min, then refusals returned at 00:34 → **12 is the sustainable
  level**; change the throttle in steps. Overnight 09-04: refusal waves ~00:00, 03–04 h, 11 h regardless of level;
  clean hours move 150–300 GB/h; 41 tasks gave up after 6 retries. Script updated 09-04 11:50 for the resume pass:
  12 retries × 300 s, wall time 12 h, and SIZE_MISMATCH "?" files (transfer succeeded, size query refused) are
  re-queried and marked OK instead of re-sent. Resume: resubmit the same sbatch — only OK files in `logs/upload_status.tsv` are skipped, so FAILED /
  SIZE_MISMATCH files re-upload; `logs/md5sums.tsv` keeps the local md5s. If 530 returns for hours: ENA helpdesk.
  Speed: FTP ≈ 9 MB/s per stream (66 GB/h at 4 streams). Aspera (`module load aspera-connect/3.9.6`, ENA's recommended
  transport) FAILS from IBEX: TCP 33001 to webin.ebi.ac.uk is refused (port 22 is not Aspera) → ask IBEX support to
  open TCP+UDP 33001 to webin.ebi.ac.uk; RAK-PR never used Aspera (curl FTP, 86 GB).
  ⚠ Never run two upload arrays at once (same task index = same run): on 09-03 the old and new arrays overlapped on
  13Dr2 / 14Dr2 / 16Sr3 (0-byte remote copies); their OK lines were dropped so they re-upload. FINAL CHECK before
  Webin-CLI: Webin-CLI compares the remote md5 to the local file — any corrupt remote copy fails validation there. + insert-size array on the
  binning BAMs `wgs_assemblies/binning/<sample>/<sample>_metaG.reads.sorted.bam` (first run on the mobilome vOTU BAMs
  was cancelled) running (2026-09-03) → rebuild manifests → validate → submit (her word). After upload: compare `md5sums.tsv` with the
  sequencing centre's lists in `raw/*/Lane*/version_01/meta/md5sum/` (batches A, C, D) as an integrity check.
- Umbrella: PRJEB104209 cannot become an umbrella (ENA refuses type conversion). NEW umbrella "Rub' al-Khali Desert
  Microbiome" (NAME RAK) with children RAK-16S, RAK-PR, RAK-WGS goes LAST, for the PI; script mode not written.
- Ecology-data owner: site 52 fix packet below; live KG also mislabels RAK_1000068/1000070 vs `site_aliases.tsv`;
  43 special aliases (18 EB blanks, 8 negatives, 17 trip-5 T/O/R/RE) untouched, 24 controls carry coral taxid 1242988.
- RAK-PR side: put PRJEB125153 into `paper/sections/05_validation.tex` §5.2 and close `notes/PENDING.md` §3.
- MAGs Methods corrections: two instruments, not one; Nov-2025 batch is 2×101; kit per batch; insert size measured,
  not the 350 bp target; BioSample depth 40 mm vs "5–10 cm" in M.1 (decide, then one MODIFY pass).
- Release coupling: releasing RAK-WGS exposes the shared BioSamples before the amplicon hold (2026-11-25).
- Operational: `ena_wgs/cache/samples/` = pre-cleanup records; `logs/cleanup.prod.tsv` blocks every OK alias
  (a future pass needs `--refresh` + renamed log).

## Decisions (Alejandra, 2026-09-03)
| item | decision |
|---|---|
| studies | RAK = new umbrella (last) · RAK-16S = PRJEB106069 (retitled ✅) · RAK-PR = PRJEB111735 · **RAK-WGS = PRJEB125153**, hold 2026-11-01; no short names inside titles; "Rub' al-Khali" hyphenated |
| BioSample title / description | bare id · "Rub' al-Khali desert soil, site N, <deep\|surface\|root-adjacent> compartment, expedition N" (no date, no replicate; 62 = site road (location 1), 64 = site campground) |
| project name | "Rub' al-Khali Desert Microbiome" on all soil samples |
| checklist / fields | ERC000022 for all; drop "nucleic acid amplification" (on every trip 1–4 experiment); keep target gene + primers; no extraction kit on new samples; local context = field-sheet `feature`; coordinates 4 decimals; site 52 untouched |
| reads | RAW (Rob); never rename on disk; public names `<sample>_wgs_R{1,2}.fastq.gz`; dual-batch samples = `<sample>_wgs_2024` + `_2026`; QC reruns not deposited; 48Dr2 = processed pair (raw R1 truncated, no copy anywhere) |
| run metadata | instrument from header serial: A01018 → NovaSeq 6000 (221), LL00134 → NovaSeq X (333); kit per run from index codes (QIA-UDI → QIAseq FX UDI 403; D7xx-D5xx → NEBNext Ultra II 81; none → QIAseq FX UDI ASSUMED 70, recorded only in `wgs_runs.tsv` `kit_source`); no read length; insert size measured from workbench BAMs; description = "Shotgun metagenome of soil sample ID" |

## Status
| step | state |
|---|---|
| primary assemblies (step 2) | **DONE 2026-10-26 03:30: 549/549 submitted** (Webin-CLI 10.0.0, validate 549/549 OK, submit array 54553238 + retry 54568981 for two transient FTP failures); ERZ30552864–ERZ30555284 in `wgs_assemblies.tsv`, verified against ENA's live analysis list (549, 0 mismatches); private under the study hold 2026-11-25. Next: watch the Webin Portal Analyses Report for processing errors; `fasta/` staging (165 GB) can go once all 549 are processed |
| bins (step 3) | 2026-10-26 21:30: **`wgs_bins.tsv` COMPLETE for her review** — 12,434 rows, all with coverage, GTDB r232 lineage (non-reps: `M/05_refine/_gtdbtk_r232_nonreps/run/`, 6,996 classified, domains = their representatives'), rRNA/tRNA (non-reps: full-MIMAG high 78 / medium 6,918), assembly quality 204 rRNA-tier / 12,230 many-fragments, NCBI taxon + level for every bin (0 fallbacks); scripts 14–16 written and documented; NOTHING registered or submitted for bins — her review first (sample wording in 16, taxon table, bin names) |
| RAK-WGS study | PRJEB125153 registered (2026-09-03) |
| 106 missing BioSamples | registered (`logs/register.prod.tsv`) |
| cleanup of all 1,301 soil BioSamples | done: 1,204 modified, 0 failed, 97 unchanged; verified by re-fetch |
| runs table / symlink farm / manifests | 554 runs, 551 samples, 6.60 TB; kit + instrument per run; INSERT_SIZE measured for 547 samples (job 51307977, binning BAMs; means agree with fastp peaks: 2024 ≈227 bp, Nov-2025 ≈136, Dec-2025 ≈148, Mar-2026 ≈139); none for 51Sr1, F5Dr1, F5Dr2, V32PRr1 (no BAM). Manifests rebuilt 2026-09-03 night; `wgs_runs.tsv` change uncommitted |
| upload · validate · submit | **DONE 2026-09-06: 554/554 experiments + runs registered** (46 by Webin-CLI, 508 by XML; ENA live list confirms 554/554); files "submitted" → "archived" as ENA processes them; accessions in `wgs_runs.tsv` (`ena_experiment`, `ena_run`) |
| title harmonisation | DONE 2026-09-07: the 46 Webin-CLI experiments + runs retitled from "Raw reads: <run>" to "<run>" by a one-off MODIFY pass (fetch-then-patch, script not kept; 46/46 OK, log `logs/titles_modify.prod.tsv`, receipts `receipts/titles/`); all 554 verified by re-fetch. Drop-box view never renders NOMINAL_LENGTH (both routes) |

## Raw batches (all 551 covered)
| batch | dir under `wgs_reads/raw/` | runs | TB | serial / flowcell | model |
|---|---|---|---|---|---|
| A_2024 | `2024/Lane1-4/version_01` | 152 | 1.71 | A01018 / HF3FJDSXC | NovaSeq 6000 |
| B_20251127 | `20251127/merged` | 69 | 0.65 | A01018 / HCV25DSXY | NovaSeq 6000 (2×101) |
| C_20251230 | `20251230_…/Lane1/version_01` | 39 | 0.51 | LL00134 / 235YJ5LT4 | NovaSeq X |
| D_20260308 | `20260308_…/merged_reads` (4 lanes merged) | 294 | 3.74 | LL00134 / 23J2WMLT4 | NovaSeq X |
Excluded tokens (22: controls, QC reruns, typo ids) → `ena_wgs/logs/raw_non_cohort.tsv`. Original md5 lists: `*/Lane*/version_01/meta/md5sum/`.
The 551 processed set is `wgs_assemblies/processed_reads/` (R1/R2/SE), not `wgs_reads/processed/` (150 stale).

## Site 52 fix packet (ecology-data owner)
Correct on all five expeditions: **20.82784 N, 53.57835 E** (= ENA). Trip-1/3 sheets carry site 53's fix
(20.85151 N, 53.75788 E): `metadata/samplesheets/trip1-2023.tsv:53`, `trip3-2024.tsv:53` (+ time; true visit ~16:27) and `:77` (52AD)
→ `environmental_measurements_curated.tsv:53,125,149`, `geodata/trip{1,3}_geodata.tsv`, `paper/env_table.tex:71,72,167`,
`scripts/rdf/generate_site_ontology.groovy:134-153` (RAK_1000052 = midpoint; 52AD snapped into RAK_1000053 with "saline pan").
Proof: GPS tracks `repos/empty-quarter/data/metadata/track/trip1/2023-03-20_07-13_Mon.gpx` (08:30 at site 52) and
`repos/empty-quarter/website/frontend/public/tracks/2024-02-20.gpx`. Regenerate geodata → curated table → OWL → env_table → evidence tables.

## Commands (assemblies phase, 2026-10-25)
```
cd /ibex/scratch/projects/c2014/alelopezv/repos/empty-quarter-proteome/code/ena_deposit
python3 10_build_assembly_table.py --prepare                       # sample_list.tsv: sample, contig FASTA, BAM (549)
sbatch --array=1-549%64 11_gzip_assemblies.sbatch                 # fasta/<s>_assembly.fa.gz + stats/<s>.tsv
sbatch --array=1-549%64 12_coverage_bam.sbatch                    # coverage/<s>.samtools_coverage.tsv
python3 10_build_assembly_table.py --study PRJEB125153            # wgs_assemblies.tsv + 549 manifests
MODE=validate sbatch --array=1-549%16 13_webin_cli_genome.sbatch  # local check, nothing sent
MODE=submit   sbatch --array=1-549%12 13_webin_cli_genome.sbatch  # IRREVERSIBLE, on her word
python3 10_build_assembly_table.py --fold                         # ERZ accessions into wgs_assemblies.tsv
```

## Commands (reads phase)
```
cd /ibex/scratch/projects/c2014/alelopezv/repos/empty-quarter-mags
python3 code/ena_deposit/00_build_runs_table.py --symlinks && python3 code/ena_deposit/04_generate_manifests.py --study PRJEB125153
MODE=validate sbatch --array=1-554%16 code/ena_deposit/05_webin_cli.sbatch   # local check, nothing sent
MODE=submit   sbatch --array=1-554%4  code/ena_deposit/05_webin_cli.sbatch   # IRREVERSIBLE, on her word
```
