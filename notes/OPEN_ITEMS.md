# EQ-MAGs: open items

The one open-items register (2026-09-20). Merges the former `TODO_RESTRUCTURE.md` lists A–C, the
"Open blocks" table that lived in `CORRECTIONS.md`, the open rows of `REVIEW_TRACK.md`, and the open
questions closing each dossier. Predecessors are in `archive/`. Structure lives in
`PAPER_STRUCTURE.md`; numbers in `CATALOGUE_NUMBERS.md`; what changed in `CORRECTIONS.md`.

Rule: a proof-of-concept number never enters the registry or the manuscript until it has been rerun
in a workbench `_clean/` folder and registered with its basis.

## ⚠ FIRST THING NEXT SESSION (state at 2026-09-20 session end, credits ran out)

1. **Cultured proteomes: all 722 exist, the comparison used 651.** The 71 soil-batch genomes (D<n>_<m>.<k>,
   thermophilic cultures) have Bakta at `EmptyQuarter_Data/cultures/annotation/thermophilic/bakta/<id>/<id>.faa`,
   not in `eqmgc/results/annotation/bakta/` (651 dirs). Fix at the root in
   `M/06_darkmatter/cultured_comparison_clean/scripts/step3_00_prepare.py` (add the second Bakta root), rerun
   step3_00 → 01 → 02 (DIAMOND, minutes) → 03 (UniRef50, only the 71 new proteomes, ~1 h) → 04 → 05; update
   the registry R2 rows and every "651" in `paper/main.tex` (R2, Methods) to 722. The paired test may gain
   pairs (the 71 include the 57 shared thermophilic genomes lacking Bakta).
2. **Two agents were KILLED by the session limit (2026-09-20 ~18:00 Riyadh), state uncertain:**
   - The 722-proteome extension (item 1) did NOT start; nothing in `cultured_comparison_clean/` was changed.
   - The PERMANOVA hardening (task 7) had edited `darkfam_ecology_clean/scripts/` (eco_stats `strata`
     argument, eco_3 within-site p column, a new `eco_0` depth step parsing the 573 fastp JSONs in
     `EmptyQuarter_Data/soil/wgs_assemblies/processed_reads/*_report.json`, since no per-sample depth
     table exists on disk) and was re-running the chain for a determinism check when it died. FIRST verify:
     `git`-less folder, so diff the scripts against the 2026-09-20 job 52135803 state described in the
     README, check `logs/` for the last eco_run job and whether it COMPLETED, confirm every registered
     number (0.132 / 0.089 / 0.092, 0.030, 0.044, 0.079, vanadium 0.003, partition 0.106/0.055/0.057,
     Mantel 0.341, indicators md5) is unchanged, then read the within-site p-values and depth-conditioned
     R² from `permanova_blocks.tsv` and the README STATUS, register, and fix the Methods sentence that
     currently concedes free permutations.
   - Task 8 (archaeal ecology) COMPLETED and is registered + written (R4).
3. **UniRef50 rerun job 52131486** (cultured proteomes, 12 tasks): when done, flip `UNIREF_HITS` (see B2).
   Combine with item 1 so the leg runs once on 722.
4. Alejandra's open pushback (unanswered when credits ran out): (a) she considers the structural-novelty
   statement an important finding and questioned the Pfam-driven rewording of R5's closing sentence; the
   heading and all numbers are unchanged, only "known parts in new arrangements" became "parts that resemble
   known folds structurally but are not recognisable sequence domains". Revisit the wording WITH her.
   (b) She feels the recent analyses (contamination especially) got too much weight in the Discussion;
   candidate trim = the contamination item in the competing-explanations paragraph and the Methods
   blank paragraph. (c) Her third point was cut off.

## A. Decisions for Alejandra

| # | decision | bears on |
|---|---|---|
| A1 | Confirm the seven-section structure: former R7 (trace gas, BGC, resistome, biosafety) dissolved to Extended Data + Discussion; former R8 heading withdrawn ("Salt, not aridity" is false for the dark families) | whole paper |
| A2 | Two reagent contaminants CONFIRMED by the no-template blank (2026-09-20): *E. coli* `V42PRr2__single_SemiBin_33` and *C. acnes* `51Dr3__single_SemiBin_390`. Keep both with a contamination annotation (current draft), or remove and move the locked 5,438. PI call. Never write anthrax | R1, R2, Discussion |
| A2b | Vector sequence in the catalogue: the no-template assembly is a laboratory plasmid (TEM β-lactamase, aminoglycoside phosphotransferase, chloramphenicol acetyltransferase) and matches contigs in 10 representatives. Check how many AMRFinderPlus hits in the resistome census sit on those contigs before R7's β-lactam and aminoglycoside numbers are defended | R7 |
| A3 | Cultured comparison: write it as "23 well-referenced bacterial species", or hold it back. It cannot reach the archaea | R2 |
| A4 | Structural armour before submission: confidence standardisation is done and reported; Pfam 38.2 on the 4,880 K members is half a day; domain calling needs an install; ESM Atlas is not runnable here | R5 |
| A5 | Target journal, once the draft exists (drafted to the tightest limits) | word budget |
| A6 | Strain persistence as a second study (feasible, ~1,900 species, needs a full re-map of 551 samples) | Discussion outlook |
| A7 | Authorship and byline | front matter |
| A8 | Claim registry `reproducibility/claims.tsv`: no rows for §2.1/§3.1, withdrawn values elsewhere. PI rebinds at end-game; someone has to raise it | M.14 |
| A9 | Citation decisions D-1 (MGnify90 biome clause is factually inverted: restate or delete), D-2, D-3 ("effectively absent" → "poorly constrained"), D-12, D-16; `tab:picrust`/`fig:picrust` on the dead basis: cut or keep; Tawfiq companion citation "in preparation" or drop | Intro, R6, Discussion |
| A10 | XRF units and instrument; aridity index named as de Martonne P/(T+10) and its climate source; does pH exist anywhere in the companion data | Methods (sampling) — BLOCKED on her |
| A11 | Data availability: MAG deposition target, archived-release DOI, whether the companion KG is cited | Methods (availability) |
| A12 | Within-population pN/pS as a second selection signal: NOT possible from existing outputs (the inStrain profiles map to the 18,604 vOTUs; 0.1 % of dark member genes touched). New run = bowtie2 + inStrain 1.10.0 on 551 samples against `darkfam_members.fna` + `ctrl_members.fna` together (1.085 M CDS, 0.49 Gbp; both arms symmetric; cost anchor 1,353 CPU-h for the CoverM pass on the same reference → budget 3,000–6,000 CPU-h, 16 CPU / 64 G / 24 h × 551). Recommendation: revision-response item, not needed for submission now that the frame-shifted null calibrates R2 | R2 |

## B. Compute

| # | item | state | cost | unblocks |
|---|---|---|---|---|
| B1 | Promote `eqmags_poc/01_dark_by_lineage` to `M/06_darkmatter/lineage_darkness_clean/`, rerun, register | **DONE 2026-09-20**, byte-identical | | R3, R4 |
| B2 | Promote `eqmags_poc/05_rakpr_overlap` to `M/06_darkmatter/cultured_comparison_clean/`, rerun, register | **DONE 2026-09-20** except the UniRef50 leg: rerun job 52131486 in flight (7–9 h); when it lands, point `UNIREF_HITS` in `step3_04_summarise.py` at `tmp/uniref_hits`, rerun step3_04 and step3_05, confirm, delete `results/uniref50_hits_job51696940/`, update the README STATUS | hours | closes the last copied input |
| B3 | Promote `eqmags_poc/03_genome_size_function` dark-vs-size alongside B1 | **DONE 2026-09-20** | | R3, ED |
| B4 | Non-coding null arm for the selection test (frame-shifted or shuffled CDS through `fam_6c_dnds.sh`) | proposed | one array job | R2's "coding sequence" claim gets a calibrated 1.0 |
| B5 | Profile shotgun blanks SNTC and SPositive (`trip4/controls`) against the 5,438 | to do | hours | contamination argument in R2/R6; A2 |
| B6 | R7 PERMANOVA with permutations restricted within site, plus a per-sample sequencing-depth covariate (depth is now registered) | to do | minutes | answers the pseudoreplication and detection-direction attacks |
| B7 | Thermoplasmatota slice of the existing structural no-hits (~2 h, no new search) | optional | 2 h | may move R4 from annotational to structural |
| B13 | Gene-level check of the compartment contrasts (if a referee asks; decided NOT to run 2026-09-22): annotate the assembly protein catalogue with the ~150 orthologue profiles behind the 23 passing modules and known traits, take per-sample gene coverage from the binning contig depths, recompute the paired root/bulk and surface/subsurface contrasts on all assembled genes, unbinned half included. Confirms or softens "mapped half only"; loses the phylum control. Hours of compute, no new mapping | deferred | hours | one sentence in R7 |
| B14 | Culture–MAG pair anatomy: REGISTERED and written into R2 2026-09-22 | done | | R2 |
| L1 | Literature pass the August search never ran: regional and national soil genome surveys with ecology (Microflora Danica and companions, SMAG / Ma 2023, Diamond 2022, Australian and Chinese soil MAG sets); read for how they frame novelty, ecology and reference coverage; record in READING_NOTES + bib; then revise the Discussion placement paragraph and Intro paragraph 1. Use the `mags-literature` agent. Why missed: the search was organised by paper section with only global and desert comparator shelves | to do | hours | Discussion, Intro |
| B16 | Comparative genomics of the 237 against every known genome of their two classes (21 + 197 GTDB r232 species representatives, plus any NCBI/MAG genomes of c__SW-10-69-26 and c__E2 not in GTDB): ANI/AAI and placement, pangenome of the classes (core/accessory), which KOs and which dark families are specific to the desert genomes vs shared with the known members, whether the 1,031 structural candidates include Thermoplasmatota-specific families, and whether the two classes appear in other catalogues (GEM, SPIRE, MGnify, Microflora Danica, SMAG) at all. Answers "are these genomes and their proteins new" beyond taxonomy | next session | 1–2 days | R4, Discussion |
| L2 | Thermoplasmatota literature (running 2026-09-22 with the mags-literature agent → READING_NOTES + bib): what is known of c__E2 and c__SW-10-69-26 and of Penumbrarchaeia (Maeke 2025); ecology, metabolism and cultivation status of MAG-only Thermoplasmatota classes; soil archaea surveys; anything that pre-empts or supports "aerobic peptide/polymer scavenger, BCAA auxotroph" | running | hours | R4, Discussion |
| B15 | Thermoplasmatota encoded metabolism: DONE, registered, R4 paragraph written 2026-09-22 (`lineage_darkness_clean/results/v2_thermo_metabolism_*`). Note for the true-up: gapseq drafts cover only 24 of the 237, no Bakta annotation exists for the representatives, so transporter/CAZyme counts are KO-based floors | done | | R4 |
| B7b | R4 physiology by clean absences (methanogenesis, ammonia oxidation, rhodopsin markers) and the occurrence/depth pattern of the 237, on the registered KO table; the proposal's "481/549 samples, deeper at 40/45 sites, no mcrABG/amoABC/rhodopsin" has no registered source and is NOT in the draft | to do | minutes | one R4 paragraph |
| B8 | Size-normalise the BGC and AMR censuses, or present them as description only in the ED table | to do | minutes | ED table |
| B9 | MOVED TO P3 (10 Oct 2026; `lifeonair/`): life-on-air text is commented out of the paper; strict-call sensitivity done there (H1 −2.19 vs −2.23). Original item: Identity/coverage floor on the HydDB best-hit (weakest 1l call 37.8 % id / 20.5 % cov) and report sensitivity of 968; state that 2,815 of 4,727 coxL proteins are unresolved so 1,170 is a floor | to do | minutes | ED table, Methods |
| B10 | Rerun §2.1 cohort now that all four CoverM profiles are repaired (551 usable, not 548), or pin the August directory state in Methods | decision + minutes | | R6 |
| B11 | GUNC on the 14,426 refined bins (tool and reference staged) or fence the claim | optional | hours | R1 |
| B12 | rRNA/tRNA detection on the 5,438 for the MIMAG sentence, or say completeness and contamination only | decision + hours | | R1 |

Done in the exploratory folder on 2026-09-17, awaiting B1/B2: paired cultured comparison (16 of 23
lower, median −0.33 pp); matched null (0.310 vs 0.048 archaea 218 pairs; vs 0.069 bacteria 237
pairs); dark-vs-size on the native key (slope −0.0957, R² 0.067); Thermoplasmatota class tree
(Slurm 52004102, 237 + 13 outgroup).

Done and registered 2026-09-17 (uncommitted): geochemistry block with S/Cl/Na (0.132 / 0.089);
variance partition on the XRF cohort; occupancy on the corrected union (37 of 61); structural floor
standardised (6.6x / 8.1x, alnTM null); 58.3 %; 419 of 448; sequencing depth per sample; the two
named lost assemblies.

## C. Disclosures owed in Methods (no compute)

- Control families admitted on member count alone (≥3) while dark families also cleared ≥3 samples
  and ≥2 sites; shipped query sets not length-matched (dark median 105 aa vs control 75 aa);
  direction of both biases is conservative.
- Bin-set switch: admission and control on the 37,949 unrefined HQ+MQ bins, phylum attribution on
  the refined catalogue (74 % attribution coverage).
- Structural model selection preferred family members already folded by the PI (3,451 of 17,221;
  K rate 42.8 % vs 35.3 % on our own folds; headline survives without them).
- Detection = CoverM relative abundance > 0, an inherited 10 % breadth default; family presence =
  any member gene covered fraction ≥ 0.5, length-biased, so occupancy is compared within length strata.
- Ectoine: two independent methods, not three.
- Hydrogenase and CODH typing is DIAMOND best-hit against group-labelled references; trees are
  display only. The five hand-dropped `_sub` representatives (970/1,173 → 968/1,170) are in no script.
- `easy-cluster -s 4` default; aridity index is de Martonne; salinity proxy = Na + Cl + S, unweighted,
  defined nowhere in the paper and spatially confounded (41 % retention), so per-element axes are used.
- MIMAG: only completeness and contamination applied; contamination ceiling ≤10 inclusive (2 bins at 10.0).
- Gene prediction ran on 552 assembly files; canonical assemblies are 549 (`51Sr1`, `F5Dr1` lost).
- Two GTDB releases were used during the analysis (1,898 reclassifications r226 → R232); all live
  numbers are R232.

## D. Records to fix

- `PIPELINE_GUIDE.md`: §8 still headed provisional (stage 05c running); line ~490 cites
  `refbias_2_magko.py:82` for a covered-fraction rule that is not there; r226 warning stale; "GUNC
  never ran" is false for the 10 order-gap candidates.
- `CATALOGUE_NUMBERS.md`: section labels still follow the v1 numbering (§1.1–§3.1) with R2/R3/R4
  tables inserted after §1.4; relabel to R1–R7 at the true-up. Register the adjudication path
  (`M/rebuild_20260729/coassembly_novelty/taxonomic_novelty_evidence/legacy_phylogeny_adjudication_r232_v2/`);
  add R3/R4 tables when B1/B2 land; move 10.8x and "1.5x, p = 0.011" to the superseded map.
- `reproducibility/claims.tsv` and `REPORT.md` pass on withdrawn numbers; `run_all.sh` calls
  dead-basis figure scripts. PI's to rebind (A8).
- Supplementary Note 1: "60–80 Tg yr⁻¹" conflates Ehhalt's soil term (60 +30/−20) with total sinks;
  `ENERGETICS.md` and the script disagree on three values (script authoritative). Only matters if the
  note survives.
- `figures/fig_darkfam_reality.pdf` referenced in v1 does not exist.
- `env_matrix.tsv` cannot be rebuilt (sources gone); availability statement must ship the matrix.
- Literature to ingest into `READING_NOTES.md` and `references.bib`: AFESM v2, Schaeffer 2026, SMAG
  (Ma 2023), Maeke 2025 (bib added), Diamond 2022, Rinke 2019, Hutchinson 2025 (bib added), Reji 2025
  (bib added), Ray 2023 (bib added). **Microflora Danica** (the Danish nationwide soil metagenome and MAG
  catalogue with ecology, Singleton et al.) is NOT in the literature folder or the bib and is the closest
  national-scale soil MAG comparator: download and read before the Discussion is finalised.
- `notes/README.md` updated 2026-09-20 for the new layout.
- Hutchinson et al. is now published (ISME J 20:wrag020, 2026; bib entry added): re-check the 58 %
  (322/554) figure quoted in the Discussion against the published version.
- Figure 6 panel b (`fig_darkmatter_permanova.pdf`) still shows the superseded geochemistry block
  (0.123 / 0.080); rebuild on the eleven-axis block. Figures 2 and 3 are unbuilt placeholders.
- The C2 cross-feeding pilot (SMETANA, 16 samples) is cut from this paper.

## E. Blocked

| item | on |
|---|---|
| Sampling-methods facts: XRF workflow and units, aridity source, samples collected | Alejandra (A10) |
| Availability: reads = PRJEB125153 (554 runs, hold 2026-11-01); genomes and code target, DOI | Alejandra / PI (A11) |
| E. coli representative in or out of the 5,438 | PI (A2) |
| Plant-host genus axes (amplicon cache off-IBEX) | infrastructure |

## Not doing

- Nothing from the exploratory folders enters the paper without a registered rerun.
- Nothing from the `dark_function` annotation project enters this paper (decided 2026-09-20).
- The environment matrix is not rebuilt; work with what is registered and say so.
- No geographic-endemism claim; database-relative absence under a named search rule only.
- No "radiation", no "genomes worldwide", no "hydrogenase-rich", no anthrax, no 10.8x, no 1.5x.
