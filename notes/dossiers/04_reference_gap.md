# Dossier 04 — the reference gap (§2.1 / M.8 · §3.1 / M.12)

One argument in two halves: references cannot **predict** this gene content (§2.1) and do not **contain** it (§3.1). Both closed and shipped. `REG` = present in `notes/CATALOGUE_NUMBERS.md`, which is what registration means there (header line 8; no literal tag exists). `M = /ibex/scratch/projects/c2014/rob/metagenomics`.

## Questions (both)

1. §2.1 — does confident placement make gene content predictable? Read the predictor on its own terms (NSTI, recall, rank agreement), then measure what it cannot reach at any accuracy.
2. §3.1 — are the dark families absent from the named catalogues, and is that absence an artefact of (a) being short and rare, or (b) catalogue size under a fixed E-value threshold?

## Basis

**§2.1** `M/06_reference_bias/refbias_clean/`, 5,438 reps on `rep_ko_v2`. Ceiling denominator = **19,186,853** Bakta protein-coding genes (`>` headers in the `.faa`, *not* the eggNOG table, which omits queries with no seed hit); numerator **8,516,190** with a non-`-` `KEGG_ko`. Pooled protein-weighted, not a mean of per-rep ratios. Alternatives: CheckM2 `Total_Coding_Sequences` **20,445,320** → 58.35 %; HQ-only (2,259 reps) → 53.94 %. Emittable reference = **10,543** KO families = every row of `ko_metagenome_unstrat.tsv`. Cohort = **448** = (548 CoverM-profiled) ∩ (non-control amplicon `sample_id`), exact string join, no normalisation.

**§3.1** `M/06_darkmatter/specificity_clean/`. Query = **212,084** in ONE pass: 106,042 dark family reps + 106,042 matched annotated controls (`query_all.faa`; raref log `unmatched-query rows: 0`). MGnify90 717,738,164 · GEM 122,978,191 · SPIRE reps 269,496,580 · deserts_combined 147,531,178 = **1.258 G**. Five deserts: atacama 6,193,359 · desert_autotroph 7,200,701 · arid 30,569,562 · antarctic_colddesert 49,206,550 · mojave 54,361,006.

## How it ran

| leg | tool + parameters as run |
|---|---|
| §2.1 prediction | PICRUSt2 **v2.4.1**, 10,543 KOs × 1,270 profiles + per-profile weighted NSTI, `M/rebuild_20260729/picrust2_current_20260803/inputs/`; sha256 `c71c2256…` re-verified this session (3/3 OK) but **never enforced in code** |
| §2.1 genome side | rep KO catalogue × CoverM rel-abund; presence = **`v > 0`** (`refbias_2_magko.py:82`), i.e. CoverM's never-passed `--min-covered-fraction` 10 % breadth default |
| §2.1 metrics | recall = \|shared\|/\|MAG\| (directional — predictor-only KOs never penalise) · Jaccard = \|shared\|/\|union\| · Spearman over the **intersection only**, no zero-filling, own impl., average ranks |
| §2.1 duplicates | 4 of the 448 sequenced twice (`F1PRr2`, `F41Dr3`, `F58Sr3`, `T1Dr1`) → sum, then within-sample renormalise (applied to all 448) |
| §3.1 primary | **DIAMOND v2.0.15** `blastp --sensitive --evalue 1e-3 -k 1`; ≥1 hit = occurs, 0 = absent under the rule |
| §3.1 rarefaction | one pass `--sensitive --evalue 1 -k 500 --max-hsps 1` vs MGnify90; curve **analytic**, no subsampled DBs. `--max-hsps 1` is load-bearing: `h(p)` counts distinct *subjects* |
| §3.1 depth | **MMseqs2 v18.8cc5c** `search -s 7.5 -e 1e-3 --max-seqs 10` vs MGnify90 **2026_07**, 1,661,142,385 seqs — PI's recipe; different aligner **and** release, so reported separately, never a `tab:specificity` row, never unioned |

**No NSTI filter, and no filter of any kind, is applied on the amplicon side** — no minimum depth, richness or NSTI cut. Two visibly broken profiles stay in (`F1Dr1` 1,948 predicted KOs, ρ 0.141; `T7PRr1` 2,181, ρ 0.442) and set the `min` column of the summary.

## Controls — CRITICAL

**The matched annotated control** (`darkfam_basis_clean/scripts/fam_5_control.sh`, seed `20260813`) is searched in the *same* DIAMOND pass, one hit table split by the two ID lists (`spec_1_global.sbatch:39`). Construction: (1) pool = non-dark 90 % reps on contigs of the same **37,949 HQ+MQ bins** (`assert len(hqmq)==37949`), darkness excluded against the **full 69,938,628** dark set; (2) **length match at pool level** — the dark 10-aa histogram over all 4,292,446 dark reps reproduced by stratified sampling, log `sampled 3757129 / target 4292446 (shortfall 535317, 12.47%)`, almost all in bins 20-29 aa (409,948) and 30-39 aa (125,366); (3) clustered **identically**, `mmseqs easy-cluster --min-seq-id 0.5 -c 0.8 --cov-mode 0`; (4) **size match** by family-size bin (3–9 exact, 10-19, 20+), one per dark family, shortfall 0.

⚠ **Two asymmetries I verified on disk, neither stated in the paper.** (a) **Admission rules differ**: dark needs `members ≥3 AND samples ≥3 AND sites ≥2` (`fam_2_recurrence.sh:87-97`), control needs only `m >= 3` (`fam_5_control.sh:180`). (b) **The shipped query sets are not length-matched**: measured from the FASTAs, dark reps median **105 aa** (IQR 62–182), control **75 aa** (IQR 47–142). Matching happened before clustering and admission; the admitted dark families ended up longer.

**Judgment (INFERENCE, direction only).** Both asymmetries push the *control* toward more absence — shorter queries score worse at a fixed E-value, less-recurrent families are rarer. A properly matched control would sit **below** 0.41 %, so 28.74 / 0.41 is conservative. **The control is sound as a floor**, but it is not a symmetric twin: "length- and size-matched" is true of the construction, not of the shipped FASTAs, and a referee who measures the two files will find the 105-vs-75 gap.

**The size correction** (`raref_curve.py` docstring 12-26). Retain each subject with probability `p`; `E_sub = p·E_full`; a hit is reportable at the 1e-3 rule iff `E_full ≤ 1e-3/p`; so `h(p) = #{subjects with E_full ≤ 1e-3/p}` and `P(no reportable hit) = (1-p)^h(p)`. Expected absence = the mean of that over families, read at `p = n_target/717,738,164`. Uniform random retention holds composition constant **by construction**, so the residual is what size cannot explain. **Gap = expected − observed; positive = the catalogue leaves LESS unexplained than its size predicts.**

**Judgment (labelled).** Defensible and unusually well controlled: the model self-validates at `p = 1` to **0.15 pt** (52.68 observed vs 52.53 expected dark; 0.01 control), which is the check that licenses reading the rest. Three limits, all disclosed in M.12/S7: the `-k 500` cap truncates `h(p)` and so **inflates the expectation and therefore the gap**, bounded by `(1-p)^500 ≈ 1.3 %` at the smallest anchor; SPIRE's expectation is read at grid anchor 271,000,000 against a true 269,496,580 (~0.1 pt); the grid stops at `p = 1e-3` because the search was at `--evalue 1`. The **+24.33** claim rests on one genuinely untestable assumption — that a real catalogue of that size behaves like a random MGnify90 slice. `atacama` proves it can fail badly (control gap **−29.97**), which is exactly why the gap is restricted to the large, broad catalogues; that restriction is the right call and is made in prose, caption and S7.

## Numbers

**§2.1** — REG (`CATALOGUE_NUMBERS.md` §2.1). I recomputed every row from `ceiling_per_rep.tsv` and `per_sample_comparison.tsv`; all reproduce.

| number | value |
|---|---|
| no-KO ceiling | **55.6 %** (55.61) of 19,186,853; 8,516,190 with KO |
| CheckM2 denominator | **58.35 %** of 20,445,320 — the paper says 58.4, which rounds to **58.3** (defect 3) |
| HQ only | 53.9 % (53.94), 2,259 reps |
| by binning provenance *(not registered)* | single_sample 54.97 % · self_supervised 56.09 % |
| weighted NSTI | median 0.159 (IQR 0.115–0.209; min 0.005, max 0.523) |
| presence-recall | median 86.0 % (IQR 85.0–87.1; **min 0.200**) |
| Spearman ρ, shared KOs | median 0.839 (IQR 0.813–0.858; **min 0.141**) |
| Jaccard | median **0.634** (IQR 0.616–0.646) — computed, named in M.8, never reported |
| MAG KO richness | median 8,666 (IQR 8,169–9,003) |
| KO universe | 11,374 encoded (≥1 of the 448; catalogue total 11,395) · 8,857 inside · **2,517 outside (22.1 %)** · **1,686** reference-only; median 14.0 % per sample |

**§3.1** — REG (`CATALOGUE_NUMBERS.md` §3.1). Dark / matched control, % absent:

| catalogue | seqs | dark | ctrl | size-only expected | gap |
|---|---|---|---|---|---|
| MGnify90 | 717,738,164 | 52.68 | 1.36 | 52.53 | −0.15 (self-check) |
| GEM | 122,978,191 | 64.38 | 2.97 | 68.65 | +4.27 |
| SPIRE reps | 269,496,580 | 39.90 | 0.98 | 61.20 | **+21.30** |
| Five deserts combined | 147,531,178 | 42.60 | 2.07 | 66.93 | **+24.33** |
| *atacama — not interpretable* | 6,193,359 | 95.61 | 43.04 | 90.94 | −4.67 (ctrl **−29.97**) |
| **Union of all four** | **1.258 G** | **28.74 %** (30,479) | **0.41 %** (435) | — | — |

GEM ∪ SPIRE (known *genome* space, 392.5 M) 39.10 / 0.89. **Depth leg**: of the 30,479 dark with no DIAMOND hit, 1,760 gain one at 1.66 G → **28,719 unmatched = 27.08 %**; control 149 of 435 → 286 = 0.27 %. **Depth asymmetry: the deeper reference resolves 34.3 % of the control's residual absences but 5.8 % of the dark set's.** Measured E-value/DB-size confound: 2,544 families (2.40 %) hit only because the desert DBs were searched split rather than concatenated. Occupancy: median **36 of 61** locations for the absent set vs 38 with a hit; strata 44/50 · 38/40 · 33/37 · 31/32 · 25/28 — but see defect 1.

**The PI's parallel rerun (comparison only).** `rebuild_20260729/protein_novelty/mgnify_2026_07/integrated/mgnify_family_summary.json`: 299,927 queries, 155,609 with any hit → **48.1 % no-hit**. It differs because it is a different object (UniRef50-dark family queries: `mag_bakta_nohits` 215,652 + `corrected_refined_tier1` 83,219 + 1,056), different clustering, and **unconditional** over his whole universe, whereas our 27.08 % is conditional on already being absent from four DIAMOND catalogues. Not interchangeable with 28.74 %. His summary states our framing verbatim: *"no hit means not detected in this versioned reference under these thresholds; it does not establish geographic endemicity."*

## Framing locks

`PAPER_STRUCTURE.md:73-75`, verbatim: **"Still LOCKED: database-relative specificity under an explicit search rule, NOT geographic endemism; provenance-level wording only — 'no Arabian or Saharan site in any named desert comparator', never 'absent from all catalogues'."** (§3.1 reframed exceptionalism → reference coverage 2026-08-17, approved, do not re-litigate.) `PAPER_STRUCTURE.md:40` for §2.1: **"Quote the no-KO ceiling, NOT the over-prediction ratio (still withdrawn)."** Gap sign, `CATALOGUE_NUMBERS.md:260`: *"Never write 'a positive gap explains more absence' — that states the opposite of the result."*

## Caveats (quoted)

- M.12: *"a representative with none is absent from it under this rule, which is a statement about the catalogue and the threshold rather than about the organism's distribution."* · *"Raw cross-catalogue ordering is therefore treated as descriptive."*
- M.8: *"…it is reported as a coverage difference, not as prediction error."*
- `specificity_clean/README.md:75-80`: *"⚠ atacama is the exception and must NOT be read naively… State the gap only for the large, broad catalogues; flag atacama's narrowness explicitly."*
- `DISCUSSION_NOTES.md:181-190`: *"A tool default, not a chosen threshold, defines 'absent' across §2.1–2.4"* (CoverM's never-passed 10 % breadth).
- `DISCUSSION_NOTES.md:144-146`: *"the amplicon survey covers 1,217 non-control samples and the metagenome cohort 551; only 448 overlap. Neither paper's cohort is a subset of the other's."*

## Status

Compute **CLOSED and SHIPPED** (§2.1 2026-08-16, §3.1 2026-08-17; `CORRECTIONS.md` rows 43/58). Methods numbering frozen: **§2.1 = M.8, §3.1 = M.12**. §3.1 ships `tab:specificity` + `fig:specificity` + Supp S7; **§2.1 ships no live float** — `fig:refbias` (main.tex:1389) is a session-09 candidate, uncited in prose. Not a defect: eggNOG binary **2.1.13** (M.4) vs reference data **2.1.12** (M.8) is a recorded deliberate choice (`CORRECTIONS.md:400`) — do not "harmonise" it.

⚠ **Six defects found this session, none previously recorded.**

1. **The occupancy panel is pre-SPIRE.** `b6_absence_vs_occupancy.tsv` and `b6_absence_occupancy_by_length.tsv` (written 11:56) split 106,042 as **35,158 absent / 70,884 with a hit** — exactly the `mgnify90+gem+deserts` union, *not* the 30,479 headline; the five strata sum to those two figures exactly. `spire_hit.ids` and `dark_nohit_union_all.ids` were rewritten at 14:54/14:55. So §3.1's "median 36 of the 61 … against 38" and `fig:specificity`b rest on a superseded set. Fix = rerun `b6_occupancy_join.py`, which already reads the corrected `dark_nohit_union_all.ids`.
2. **The 548/448 cohort is stale against its own input.** `refbias_2_magko.py` builds the cohort by listing `profiles_full_recompute/`. The PI repaired all four failed CoverM tasks on 2026-09-06/08 (`F5Dr1`, `SPositive`, `51Sr1`, `F13PRr1`, all now `PROFILE_COMPLETE`). Re-running the script unchanged today yields **551 usable, not 548**, and more comparison samples. Nothing pins the cohort, so 548/448 reproduce only against the August directory state.
3. **58.4 % should be 58.3 %** (`main.tex:772`; computed 0.583465). The workbench README has it right.
4. **"the predictor returns its full 10,543 families in every sample — a constant" is false.** 419 of 448 return exactly 10,543; **29 return fewer**, minimum **1,948**. `comparison_summary.tsv` already records `n_pic_ko … min 1948`. Medians are robust (recall 0.8597 on the 419 vs 0.8596 on all 448); only "every" and the registry's parenthesis need fixing.
5. **`reproducibility/claims.tsv` is stale for both sections.** 166 rows on the old A0–C3 labels, no 2.1/3.1 rows; C1/C3 still hold the superseded 56.7 / 93 / 4,482 / 0.16 and 21.4 / 36.6 / 47.4 / 31.9. M.14's *"Every quantitative claim in this paper is held in a machine-checked registry"* is not true of §2.1 or §3.1 as the file stands.
6. **`PIPELINE_GUIDE.md:490` is wrong.** It cites "`refbias_2_magko.py:82` for orthologue presence", implying the `covf ≥ 0.5` rule; line 82 is `elif v > 0:` and there is no coverage-fraction test anywhere in `refbias_clean/`. The 0.5 rule belongs to the *member-aware* run (§1.4/§2.5). The same line also carries the pre-freeze "M.11" numbering.

Lower-order, one clause each if wanted: the MAG rel-abund column sums over every carrying representative, so it is **not compositional** (fine for Spearman and presence, wrong if ever read as a fraction); only the mapped-community normalisation was produced though M.5 promises both; Jaccard 0.634 is computed but unreported.

## What can and cannot be claimed

**Can.** 55.6 % of protein-coding genes carry no KO — on either denominator (58.3 / 53.9) and in both binning provenances (54.97 / 56.09), so a reference property, not a quality artefact. The predictor works within its remit: NSTI 0.159, recall 86.0 %, ρ 0.839. 22.1 % of encoded KOs (2,517/11,374) sit outside the emittable set — a vintage gap, not error — and in **420 of 448 samples zero** KOs inside the reference go unpredicted (`mag_in_ref_not_predicted` median 0), so the recall shortfall *is* the universe gap. 28.74 % dark vs 0.41 % control across 1.258 G under a named rule. Content beyond depth for the large, broad catalogues: +24.33 (deserts), +21.30 (SPIRE) against control gaps near one point. Depth asymmetry 34.3 % vs 5.8 %.

**Cannot.** Geographic endemism, or "absent from all catalogues" (MGnify90 biome labels carry no geography). Any per-catalogue *ordering* as a claim (2,544 families, 2.40 %, from splitting the DBs). `atacama`'s gap as content. Unioning the MMseqs2 depth leg with the DIAMOND absences. Prediction accuracy outside the shared-KO intersection. Any biological reading of a dark family (that is the companion annotation paper).

## Reviewer attack surface

| # | attack | answer |
|---|---|---|
| 1 | "63.4 % of 608,258 novel families are single-ecosystem (Aplakidou 2025) — 28.74 % is unremarkable." | Conceded in prose as the norm; the contribution is relocated to *which* biome is missing. **D-15**: its "ecosystem" is GOLD's 3/8-level scheme, coarser than "biome" as §3.1 uses it, and the 608,258 were *defined* by a ≥5 %-one-ecosystem rule, so restriction bias sits in its denominator. Both make the concession safer; a one-clause hedge in M.12 is still open. |
| 2 | **D-1, highest priority, open.** §3.1's *"MGnify90 itself is dominated by aquatic and engineered biomes"* is **factually inverted** — Richardson 2023 p.D754 leads with human/host-associated and marine, **soil is itself in the top nine**, and "engineered" never appears. Looks like June's dead, non-recomputable "83.5 % aquatic" re-attributed. | Unfixed. Restate to what the paper says, or delete the clause — the desert-comparator sentence already carries the argument. |
| 3 | Cohort overlap: amplicon **1,217** non-control samples, metagenome **551**, only **448** overlap; neither a subset of the other. Reconciles exactly — the 1,270 PICRUSt2 profiles = 24 controls + 1,246 non-control rows over 1,217 distinct `sample_id`s; 29 ids are duplicated, 4 of them inside the 448. | M.8 gives 448 and the 0.4-pt compartment match; the three-way is only in DISCUSSION_NOTES. Never transfer per-sample claims between the two papers. |
| 4 | **D-2, open.** "largest desert **soil** collections … no Arabian or Saharan site in any of them" cites only `uritskiy2019atacama` — one of five, halite endoliths not soil, and a negative it cannot support; three of the five are unpublished JGI studies. | Drop "soil"; attribute as "five desert metagenome studies assembled in SPIRE"; rest the negative on M.12's accession list. |
| 5 | **D-12, open.** §2.1's `matchado2024limits` is human gut + simulation, **no soil**. | Pair with `sun2020inference`, the one soil source, already in the bib. |
| 6 | "Your control isn't matched." | True at the FASTA level (105 vs 75 aa) and in the admission rule. Both conservative — see Controls. Needs a Methods clause. |
| 7 | "A fixed E-value flatters small catalogues." | Pre-empted and measured (2.40 %); the claim rides on the rarefaction, which self-validates to 0.15 pt. |
| 8 | "You kept samples with broken amplicon profiles." | True — no NSTI or depth filter anywhere; recall min 0.200, ρ min 0.141. Medians are unaffected (0.8597 vs 0.8596), but say so rather than be caught. |
| 9 | "Occupancy from a covered-fraction rule is length-biased." | Conceded and stratified — but on the pre-SPIRE set (defect 1). |

## Open questions for Alejandra

1. Rerun `b6_occupancy_join.py` on the post-SPIRE 30,479 set? (Defect 1 — the only §3.1 number still on a superseded basis, and it sits inside a figure panel.)
2. §2.1's cohort: leave 548/448 as shipped and pin the August directory state in Methods, or re-run now that the PI has repaired all four profiles and quote 551? (Defect 2 — the bigger call of the two.)
3. D-1 — restate the MGnify90 biome clause, or delete it? Same sitting: D-2 and D-12.
4. Fix 58.4 → 58.3 (defect 3) and "in every sample" → "in 419 of the 448" (defect 4)?
5. Add one Methods clause disclosing the control's residual length gap and the recurrence-gate asymmetry, with the direction of bias? (Recommend yes — cheap, and it pre-empts attack 6.)
6. Report Jaccard 0.634 or cut it from M.8? Promote `fig:refbias` to the live §2.1 display set, or leave §2.1 float-free?
7. `claims.tsv` rebinding for §2.1/§3.1 — onto the PI coordination list now, or at end-game as `CATALOGUE_NUMBERS.md` line 9 anticipates?
