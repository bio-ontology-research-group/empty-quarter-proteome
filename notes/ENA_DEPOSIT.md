# ENA deposit — 551 soil shotgun read sets (RAK-WGS = PRJEB125153)

Scripts `code/ena_deposit/`; work dir (not in git) `/ibex/scratch/projects/c2014/alelopezv/ena_wgs/`
(`fastq/` symlink farm, `manifests/`, `xml/`, `receipts/`, `logs/`, `cache/`). Password only `~/.webin_pass`.
Test service (wwwdev) = daily mirror; MODIFY of production samples cannot be rehearsed there.

## Roadmap under PRJEB125153 (decided 2026-09-03)
1. reads (554 runs) — in progress · 2. primary metagenome assemblies (549) — next, reuse eqmgc `05_generate_assembly_manifests.py`
· 3. later: bins + MAGs (proposal: 5,438 species reps as MAG type with MIMAG BioSamples and GTDB→NCBI taxonomy;
6,996 other HQ+MQ bins as binned metagenome; LQ not deposited) — decision with Rob · 4. non-ENA datasets (protein
catalogue, dark families, ESMFold models, tables) → Zenodo DOIs · 5. new umbrella (PI) + aligned release dates.

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

## Commands (reads phase)
```
cd /ibex/scratch/projects/c2014/alelopezv/repos/empty-quarter-mags
python3 code/ena_deposit/00_build_runs_table.py --symlinks && python3 code/ena_deposit/04_generate_manifests.py --study PRJEB125153
MODE=validate sbatch --array=1-554%16 code/ena_deposit/05_webin_cli.sbatch   # local check, nothing sent
MODE=submit   sbatch --array=1-554%4  code/ena_deposit/05_webin_cli.sbatch   # IRREVERSIBLE, on her word
```
