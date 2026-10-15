# darkfam_basis_clean — strict-rule re-search of the dark-family universe (§1.4/§2.5/§3.1 basis gate)

Owner: velazqal. Created 2026-08-12 (session 01-basis-and-launch).

## STATUS (2026-09-20) dN/dS calibration: frame-shifted null arm

**Why.** The paper reports 0.425 (dark, n=97,200) against 0.388 (length- and size-matched annotated
control, n=95,647). Both arms are coding, so nothing anchors the scale at 1 and the distance of
0.425 from neutrality is uncalibrated.

**Arm** (implemented at the root in `scripts/fam_6c_dnds.sh`, no wrapper script). `ARM=null_frame1`
reads the SAME dark member CDS in the +1 frame: drop the first nucleotide, trim to a codon multiple,
rewrite every codon that translates to a stop under table 11 as `NNN`; the existing pairwise gap/N
stripping then removes those codons from the counts. `ARM=dark_subset` is the same code path with
shift 0 (true frame) and is the positive control for the masking. Every other filter is identical
(>=90 nt, protein >=25 aa, >=3 usable sequences, >=30 clean codons, 0<dS<=2, <=80 thinned pairs).
`ARM=dark` and `ARM=ctrl` are byte-for-byte unchanged (re-verified against the shipped part files)
and `results/selection/dnds.tsv` was not touched.

- **Subset:** `MODE=prep ARM=null_frame1` draws a seeded (20260920) random 10,000 of the 106,042
  evaluable (>=3 CDS) dark families -> `results/selection/family_index_null_frame1.tsv`.
- **Jobs:** null_frame1 **52135782**, dark_subset **52135794** (both `--array=0-49`, BATCH=200,
  same resources as the dark arm: 4 CPU / 16G / 12h). `MODE=merge ARM=<arm>` writes
  `results/selection/dnds_null_frame1.tsv` and `results/selection/dnds_dark_subset.tsv`;
  per-task accounting in `results/selection/dnds_qc/`.

| arm | n evaluable | median dN/dS | IQR | frac <1 | failed >=3 usable seq |
|---|---|---|---|---|---|
| null_frame1 (+1 frame) | 9,148 / 10,000 | **1.087** | 0.901 to 1.326 | 0.377 | 693 |
| dark_subset (true frame) | 9,182 / 10,000 | **0.427** | 0.295 to 0.580 | 0.970 | 635 |
| dark, full arm (shipped) | 97,200 | 0.425 | 0.297 to 0.579 | 0.967 | n/a |

- **Positive control passes:** dark_subset 0.427 vs 0.425 for the full arm, and 94.2% of the subset
  families return exactly the shipped per-family value, so the masking code path does not move the
  true-frame answer. Paired over the 9,029 families evaluable in both arms the null sits 0.655
  higher (median) and exceeds the true frame in 92.4% of families.
- **Stop masking:** median 1.98% of codons per family masked in the +1 frame (mean 2.44%) vs 0.72%
  in the true frame (terminal stop plus occasional internal ones). Both are well under the 4.7% a
  random sequence would give, which is itself a coding signature of these CDS.
- **CAVEAT:** a +1 frame of a coding sequence is NOT random non-coding sequence. It inherits the
  base composition, codon usage and substitution pattern of the real gene, so it calibrates "is the
  annotated frame the selected one", not neutral evolution in general.

## STATUS (2026-08-14, session 02-dark-families) — §1.4 SHIPPED
- **Delivered:** §1.4 Results rewritten + Methods M.9 + Supplementary Table S5 (`tab:darkfam-evidence`)
  + abstract minimally corrected + Discussion deposit. Commits: Overleaf `5504ea7`+`ec3b0a1`,
  GitHub main `bd45b1d`+`cbd46d4`+`55af8d4`.
- **Headline numbers (all re-verified on disk 2026-08-14 against the regenerated matrices):**
  106,042 families / 536,437 member CDS; median 12 samples / 9 sites; dN/dS 0.425 (96.7%<1,
  n=97,200) vs matched control 0.388 (99.6%<1, n=95,647; MWU p<1e-300, double-precision underflow);
  median prevalence 98/550, 4 single-sample, median peak 99 TPM; persistence 79.4% all-5 trips
  (94.8% ≥4), 28.9% re-assembled all-5; phylum 97.9% single-phylum of 78,727 attributable (74%);
  robustness 108,826→130,843 ≥3-member clusters at 50%→30% (cluster stage, NOT vs admitted count).
- **SETTLED — do not re-litigate:** universe definition (darkclean-derived, 50/80, ≥3 members /
  ≥3 samples / ≥2 sites); the matched control as the only defensible dN/dS comparison; B4/FUBAR
  retired (tail = 3.3% vs 0.4% from this universe); presence rule (any member covered_fraction
  ≥0.5, family TPM = summed members).
- **OPEN (updated 2026-08-14 evening, session 1):** ~~F13PRr1~~ COMPLETED 17:17 (25% of 133.8M
  read pairs mapped; fam_8 refresh 50531118 → matrix 551/551, all digits unchanged: 98/551, 79.4%
  — the S5 "550/551" footnote can be DROPPED at the next manuscript touch); ~~fam_9~~ DONE
  (`scripts/fam_9_phylum.sh` materialized); ~~§2.5 rerun~~ DONE 2026-08-16 (session 4,
  `../darkfam_ecology_clean/`, shipped 32e24ae/4cdfcab);
  §3.1 specificity recompute on this universe (session 5); PI briefing (June object banner,
  verify.py b2/b3/c3 rebind).

## Why
The 23,382 dark-family June tables (dN/dS, abundance, specificity) were basis-verified:
genome set = current 5,438 dRep reps (3-way ID match vs rep_ko_v2 + dereplicated_genomes) ✅,
sample set = 551 canonical ✅. BUT their darkness call is the June decode-era flag
(76_candidate_isdark.sh: mmseqs -s 5.7 --max-seqs 5, dark = no hit OR pident<30), NOT the
paper's locked rule (darkclean: -s 7.5 -e 1e-3 --max-seqs 10, dark = no hit). Measured
exposure: 4,701/23,382 family reps (20.1%) are dark only via the <30% clause. Additionally
the June screen only covered gapsmith decode candidates: 667,566 of 985,820 >=3-member
orthogroup reps (67.7%) were never darkness-screened at all.

This folder re-searches, under the locked rule (scope A′ + gap patch):
- chunk_1: all 23,382 family reps → SURVIVING/DROPPED under strict darkness
- chunks 2-9: the 667,566 never-screened >=3-member reps → candidate NEW strict-dark families

Manuscript framing = two-stage screen→confirm search (lossless: higher sensitivity only adds
hits; screen exclusions are guaranteed exclusions under the final rule). Full history stays here
and in notes/CORRECTIONS.md.

## Inputs
- family reps: results/family_rep_ids.txt (23,382; = darkfam_ecology/catalog/darkfam_meta.tsv col1)
- gap reps: results/gap_rep_ids.txt (667,566; = og_recurrence n_members>=3 minus isdark-searched)
- sequences: 06_darkmatter/orthogroups/all_rep_proteins.faa (9.3 GB, June concat of the 5,438 proteomes)
- UniRef50 mmseqs DB: 04_structure/uniref50/uniref50 (same DB as darkclean + June search)

## FINAL design (decided 2026-08-13, Alejandra + assistant; supersedes the re-search chain below)
REBUILD, not patch: §1.4/§2.5/§3.1 families are rebuilt from the paper's own dark call so all of
Part 1 shares ONE object — §1.2's genome-resolved 4,292,446 MAG-mapped dark cluster90 reps
(darkclean_5, 37,949 HQ+MQ bins both binnings; deliberate, documented gate vs unbinned noise).
- Universe: 4.29M dark reps → mmseqs easy-cluster 50%/80% → families; admission = ≥3 diverged
  members AND ≥3 canonical samples AND ≥2 sites (instance spread via all_cluster90.tsv membership,
  551-canonical filter per review fix). CORRECTED FINAL: **106,042 admitted** (paralog-only dropped
  2,784; young recurrent 1–2-variant class 782,616; member CDS catalog 536,437).
- Four §1.4 evidence legs: recurrence (n_samples/n_sites/n_trips) · dN/dS vs rebuilt size-matched
  annotated control · member-aware CoverM abundance (551 samples) · temporal persistence (trips).
  Plus phylum attribution (contig→bin→species join). Specificity = §3.1 launch-ahead, not §1.4.
- No UniRef50 re-search needed: darkness comes from darkclean (complete, strict).
- June 23,382 object: banner-marked superseded for this paper, never deleted; PI item.

## Pipeline (scripts/)
1. fam_1_cluster.sh — extract 4.29M dark rep seqs + easy-cluster 50/80 → results/fam_* + size hist
2. fam_2_recurrence.sh — recurrence/persistence tables + ADMISSION: >=3 diverged members AND >=3 samples AND >=2 sites (paralog-only dropped; recurrent 1-2-member "young dark" reported as a named class, stays in §1.2/§1.3). Phylum attribution + 30%-recluster robustness one-liner: follow-up passes
3. (step 3) member-aware CoverM catalog + 551-sample array; dN/dS dark + control arrays
4. (step 4) specificity (GEM/MGnify/deserts/SPIRE) on family reps — §3.1 launch-ahead

## Job log
- fam_1_cluster: **50421060** (resubmit 2026-08-13 09:1x after 50413899 died on a version-print bug post-extraction; dark_reps.faa 4,292,446 verified, extraction skipped on resubmit)
- fam_2_recurrence: **50430067** (submitted 2026-08-13)
- SUPERSEDED, never run: basis_1_prep / basis_2_search / basis_3_verdict (strict re-search chain;
  kept for the record — the 20.1%-clause + 68%-never-screened findings above motivated the rebuild)

## Verified numbers so far (2026-08-12, pre-search)
- family reps strict-no-hit under June search: 18,681 (79.9%); clause-dark (<30% hit): 4,701 (20.1%)
- June search coverage: 5,954,339 proteins (decode candidates), NOT the full orthogroup rep space
- orthogroup universe: 6,043,585 reps (mmseqs easy-cluster 50%id/80%cov, 74_mmseqs_orthogroups.sh);
  >=3-member reps: 985,820
- erratum: PIPELINE_GUIDE calls families ">=3 dark members"; the builder (90_darkfam_catalog.sh)
  requires a dark REPRESENTATIVE with >=3 members

## Job log (2026-08-13 evening, post-review)
- fam_2 corrected rerun: 50439896 ✅ (106,042 admitted; canonical 551 filter)
- fam_3 CDS: 50435064 ✅ → darkfam_members.fna 536,437 (corrected universe; stale pre-correction intermediates deleted — regenerable via fam_3)
- fam_4 member-aware CoverM: **50440132** (551 tasks, 16CPU/64G/24h, resized per review)
- fam_5 control (rewritten: HQ+MQ gate, full-dark exclusion, length-stratified sampling): **50439953**
- Multi-agent review: 8 confirmed must-fixes applied, 3 refuted, 31 notes (workflow wf_de9e4efa)

## Job log (2026-08-14, wrap day)
- fam_6a control CDS slicing: 50480204 ✅ 548,831/548,831 (100%)
- fam_6b/c dN/dS: dark prep 50480237 ✅ → dark array **50480238** (0-530); ctrl prep 50480447 → ctrl array **50480461**
  (worker = verbatim June algorithm 94_darkfam_dnds.sh, outputs repointed to results/selection/)
- fam_7 robustness 30% recluster: **50480205**
- fam_8 family×sample aggregation: **50480972** (rerun in minutes when the last sample lands; 550/551, F13PRr1 running)
- Phylum attribution ✅ (results/fam_phylum.tsv): 78,727/106,042 families attributable via refined catalogue
  (12,433 bins → 5,438 reps, R232 at 05_refine/_gtdbtk_refined); **97.9% single-phylum** (old E4: 97.5%);
  top carriers Actinomycetota 32% / Thermoplasmatota 19% / Pseudomonadota 9%.
- fam_8 refresh after F13PRr1 landed: 50531118 ✅ — matrix complete 551/551; leg-3/persistence digits
  unchanged (98/551; 79.4%). fam_9_phylum.sh added (materializes the 2026-08-14 ad-hoc attribution;
  closes session-2 open item).
