# specificity_clean — §3.1 "the dark proteome and the reference gap"

**STATUS 2026-09-27 — B6 with-SMAG leg (not registered).** `b6_occupancy_join.py` gained an optional
`TAG` arg (`with_smag` = `dark_nohit_union_with_smag.ids`, 28,446 absent vs 77,596 with-hit; default
verified unchanged by rerun). Overall median n_sites: absent **37**, has-hit **38** (same as shipped
37/38). By length stratum (absent/has-hit): 44/50 · 38/40 · 34/36 · 31/32 · 27/27. Outputs
`results/b6_*_with_smag.tsv` (5 files); shipped `b6_*.tsv` untouched.

**STATUS 2026-09-24 — SMAG leg ADDED (not registered, not written).** SMAG (40,039 soil MAGs, Ma et al.
2023) searched with the identical rule. Dark **40.47%** / ctrl **1.57%** absent; size-only expectation
68.16 → gap **+27.69 pts** (ctrl +1.63). Union with SMAG: **26.83% / 0.36%** (vs shipped 28.74 / 0.41).
Separate files, shipped tables untouched: see "SMAG leg" below.
Depth check on the with-SMAG absent sets (spec_5 `TAG=with_smag`): dark **1,546/28,446 = 5.43%** gain a
hit, ctrl **138/385 = 35.8%** resolved, **26,900 remain = 25.37%** of the universe (ctrl 247 = 0.23%).
These are the exact subsets of the original run's hits (same 1,546 and 138).

**STATUS 2026-08-17 (session 05) — SECTION CLOSED: all five jobs done, written and shipped.**
Basis = the rebuilt **106,042-family** universe (`../darkfam_basis_clean/`).
**Headline: 28.74% of the dark families have no homolog in MGnify90 ∪ GEM ∪ SPIRE ∪ five deserts
(1.26 G seqs) vs 0.41% of the length/size-matched annotated control**, same DIAMOND pass. 70-fold.
**June's 21.4 / 47.4 / 36.6 / 31.9% are DEAD** (superseded 23,382-orthogroup object) — never quote them.
Numbers registry = repo `notes/CATALOGUE_NUMBERS.md` §3.1 · what changed = `notes/CORRECTIONS.md` 2026-08-17.
**Written:** §3.1 (3 ¶, reframed exceptionalism → reference coverage) + Methods **M.12** (old M.12/M.13 →
M.13/M.14) + M.9 length-disclosure sentence + abstract updated.
**Display:** `tab:specificity` rebuilt · `fig:specificity` → `fig_specificity_controls.pdf` (2 inferential
panels, generator `code/make_specificity_figures.py`) · Supp Table **S7** · **`fig:endemism` CUT**.
**Open, needs a person:** bibkeys `bay2021aridity` and Aplakidou 2025 (`\todo` in `main.tex`); registering
the new generator in the PI's `reproducibility/regenerate_figures.sh` (deliberately untouched).

Method rationale and every "why" lives in the **script docstrings**, not here. This file is numbers +
state + decisions.

## Jobs
| job | what | state |
|---|---|---|
| 50594992 | MGnify90 + GEM | ✅ 2h21m |
| 50594993 | 5 deserts + combined | ✅ 44m |
| 50613399 | SPIRE reps (makedb via FIFO → 269,496,580 seqs, 121.8 GB index) | ✅ search done; job died on a `score()` bash bug, scored manually |
| 50613400 | rarefaction vs MGnify90 (`-e 1 -k 500`) | ✅ 3h52m |
| 50618838 | depth: no-hit set vs MGnify90 2026_07 (1.66 G, MMseqs2, PI's recipe) | ✅ 1h37m |
| 52468882 → 52492434 | SMAG fetch (`spec_6a`): 373 parts md5-verified, stream-extracted, parts deleted | ✅ 3h00m + 1m (first run downloaded and extracted, then stopped at the genome-count check, see SMAG leg) |
| 52492435 [1-81] | SMAG Pyrodigal 3.7.1 `-p meta` (`spec_6b`), 500 genomes/chunk | ✅ 50s–20m per task, 20m total |
| 52498970 | depth (`spec_5`, `TAG=with_smag`): 28,446 dark + 385 ctrl vs MGnify90 2026_07, outputs `*_with_smag` | ✅ 1h45m |
| 52492436 | SMAG makedb + search + score + expectation + union (`spec_6c`); genomes, chunks, FASTA (226 GB) deleted | ✅ 1h27m (search 67m) |

## Depth robustness — the headline barely moves (`depth_robustness_summary.tsv`)
Of the **30,479** dark families with no DIAMOND hit in any of the four catalogues, only **1,760 (5.8%)**
gain a hit against the complete **1,661,142,385**-sequence MGnify90 2026_07 searched with MMseqs2
`-s 7.5 -e 1e-3`. **28,719 remain unmatched = 27.08% of the 106,042 universe** (vs 28.74% before).
Control: 435 queried → 149 hit → 286 still unmatched = **0.27%**. So a 2.3× deeper reference with a more
sensitive tool moves the headline by 1.66 points and the dark/control contrast stays ~100-fold.
**The asymmetry is itself a result, quotable without inference:** the deeper reference resolves
**149/435 = 34.3%** of the control's residual absences but only **1,760/30,479 = 5.8%** of the dark set's.
The control's leftover absences are therefore largely a depth artefact; the dark ones are not.
Surviving set: `dark_nohit_plus_deepmgnify.ids`. **MMseqs2 ≠ DIAMOND — robustness only, never a
`tab:specificity` row, never unioned with the DIAMOND absences.**

## Content vs depth — EARNED, no longer descriptive (`content_vs_depth.tsv`, `rarefaction_curve.tsv`)
The e-value-corrected rarefaction gives, for any size N, the absence a **random N-sequence slice of
MGnify90** would produce. Comparing each real catalogue against its size-matched expectation separates
content from depth. Gap = expected − observed, so **positive = explains more than size alone predicts**.

| catalogue | DB seqs | dark obs | dark exp (size only) | **gap** | ctrl obs | ctrl exp | ctrl gap |
|---|---|---|---|---|---|---|---|
| arid | 30,569,562 | 52.59 | 80.97 | **+28.38** | 3.91 | 6.26 | +2.35 |
| deserts_combined | 147,531,178 | 42.60 | 66.93 | **+24.33** | 2.07 | 3.00 | +0.93 |
| spire | 269,496,580 | 39.90 | 61.20 | **+21.30** | 0.98 | 2.21 | +1.23 |
| desert_autotroph | 7,200,701 | 72.81 | 90.23 | +17.42 | 10.48 | 11.89 | +1.41 |
| mojave | 54,361,006 | 62.54 | 76.14 | +13.60 | 6.16 | 4.84 | −1.32 |
| antarctic_colddesert | 49,206,550 | 67.55 | 77.02 | +9.47 | 8.67 | 5.07 | −3.60 |
| gem | 122,978,191 | 64.38 | 68.65 | +4.27 | 2.97 | 3.29 | +0.32 |
| atacama | 6,193,359 | 95.61 | 90.94 | −4.67 | 43.04 | 13.07 | **−29.97** |
| mgnify90 | 717,738,164 | 52.68 | 52.53 | **−0.15** | 1.36 | 1.35 | −0.01 |

**The model validates itself**: at p=1 the expectation reproduces the observation to **0.15 points**
(dark) and 0.01 (ctrl). That is the check that licenses reading the rest — lead M.12 with it.

Two exactness notes for M.12. (i) The `spire` row reads its expectation at the grid's **271,000,000**
anchor while the true index is **269,496,580** — a ~0.1-point difference, but relabel or footnote it so the
table is internally exact. (ii) The `-k 500` cap's bias **direction flipped with the reframe**: under the
gap framing, truncating `h(p)` inflates the *expectation* and therefore inflates the gap. It is negligible
— a capped query has h ≥ 500 and (1−p)^500 ≈ 1.3% even at the atacama anchor — but state it in one clause
rather than repeating the old "conservative for the absence claim" wording, which was correct only for the
bare-absence framing.

**The claim this supports.** `deserts_combined` explains **24.3 points** more of the dark proteome than a
size-matched random slice of MGnify90 would, and `spire` **21.3 points** more — while their *control*
gaps are +0.93 and +1.23, an order of magnitude smaller. So the biome-matched / genome-resolved advantage
is (a) real after correcting for the e-value/DB-size confound and (b) **specific to dark families**, not a
general property of those catalogues. This replaces the withdrawn "the reversal proves it isn't DB size"
argument with a proper one.

**⚠ atacama is the exception and must NOT be read naively.** Its control gap is **−29.97** — a single
35-sample study is a deeply redundant, narrow slice of protein space, so it is far *worse* than a random
draw even for conserved proteins. The gap metric assumes the comparator is a random draw; one narrow study
violates that badly. Read atacama's −4.67 dark gap **relative to its own −29.97 control penalty**: it is
in fact favourable, not evidence of no relevant content. **State the gap only for the large, broad
catalogues; flag atacama's narrowness explicitly.**

## Numbers
Rule: `diamond blastp --sensitive -e 1e-3 -k 1`, v2.0.15, identical to June. Query = 106,042 dark reps
**+ 106,042 length/size-matched annotated controls** in one pass (`results/query_all.faa`, 212,084 seqs,
0 ID collisions) so every catalogue gets a known-protein floor — the paper's §1.3/§1.4 control idiom,
absent from June.

| catalogue | source | DB seqs | dark %abs | ctrl %abs | June |
|---|---|---|---|---|---|
| atacama | Uritskiy_2019 | 6,193,359 | 95.61 | 43.04 | 89.1 |
| desert_autotroph | Bay_2021 | 7,200,701 | 72.81 | 10.48 | 62.8 |
| arid | JGI_soil_arid | 30,569,562 | 52.59 | 3.91 | 41.7 |
| antarctic_colddesert | JGI_soil_Antarctica | 49,206,550 | 67.55 | 8.67 | 55.9 |
| mojave | JGI_soil_desert_Mojave | 54,361,006 | 62.54 | 6.16 | 51.9 |
| deserts_combined | the five concatenated | 147,531,178 | 42.60 | 2.07 | 31.9 |
| gem | GEM 52,515 MAGs | 122,978,191 | 64.38 | 2.97 | 47.4 |
| spire | SPIRE reps (genome-level) | 269,496,580 | 39.90 | 0.98 | never ran |
| mgnify90 | MGnify90 (assembly-level) | 717,738,164 | 52.68 | 1.36 | 36.6 |

| union | DB seqs | dark %abs | ctrl %abs |
|---|---|---|---|
| gem + spire (known *genome* space) | 392.5M | 39.10 | 0.89 |
| mgnify90+gem+deserts (pre-SPIRE) | 988.2M | 33.15 | 0.54 |
| **mgnify90+gem+deserts+spire** | **1.258 G** | **28.74** | **0.41** |

**28.74% vs 0.41% = the section's central contrast, 70-fold.** SPIRE moved it 33.15 → 28.74.
**`main.tex:52`'s 21.4% AND its "pending SPIRE" hedge are both closed.** Every dark fraction is ~10 pts
above June, direction-stable — a 4.5× larger universe admits more small/rare families. Don't mix bases.

**Comparators are three KINDS, state it in Methods.** Deserts = SPIRE per-**study gene calls**
(`SAMN*.psa_megahit.prodigal.faa.gz`, every gene on every assembled contig) → assembly-level, the *more*
inclusive kind, so their absence result is **conservative**. SPIRE leg = `rep.tar`
(`spire_representative_genomes/specI_v4_*.fa.gz`) → genome-level, blind to anything unbinned. So claims
separate: not in known *assembled sequence* space (MGnify90) · not in known *genome* space (GEM+SPIRE) ·
not in five desert metagenomes *at assembly level*. Union carries the headline.

**e-value/DB-size confound is real and measured**: the 5 desert DBs searched separately find 63,413 dark
families, the same sequences as one 147.5M DB find 60,869 → **2,544 (2.40%) hit only because the DBs were
split**. So cross-catalogue ordering stays **DESCRIPTIVE**; the rarefaction carries the claim.
(A "reversal proves it isn't DB size" argument was **withdrawn** — see `raref_curve.py` docstring.)

**Dark-minus-control gap: computed, less useful than hoped.** Quote the dark/ctrl **pair**, not the
difference — where the ctrl floor is low the gap is just the dark fraction; where it is high it
compresses (atacama ctrl is itself 43.04% absent). Values: `results/union_specificity_summary.tsv`.

## SMAG leg (2026-09-24) — a fourth genome-level comparator
SMAG = Zenodo 8223844 (Ma et al. 2023 *Nat Commun* 14:7318, CC-BY-4.0), 40,039 soil MAGs, nucleotide
only. Gene-called here with the project's Pyrodigal 3.7.1 env and flags (`-p meta --no-stop-codon`),
headers `<MAG>::<gene>`. Catalogue = SD2's `user_genome` column (`results/smag_supplementary_data_2.xlsx`,
per-MAG GTDB-Tk + CheckM2; Strict OOXML, openpyxl reads no sheets). The tarball holds 40,352 files:
the 40,039 MAGs, 312 timestamped duplicate copies (byte-identical to their originals), a few bins not in SD2 and two stray files.
Only the 40,039 are used. **129,550,340 proteins**, 37,862,083,488 residues, all 40,039 MAGs yield proteins.
Rule, DIAMOND binary (v2.0.15, `eggnog` env; not module 2.1.16), queries and `score()` identical to spec_1–3.

| catalogue | DB seqs | dark obs | dark exp (size only) | **gap** | ctrl obs | ctrl exp | ctrl gap |
|---|---|---|---|---|---|---|---|
| smag | 129,550,340 | 40.47 | 68.16 | **+27.69** | 1.57 | 3.20 | +1.63 |

The expectation is read at SMAG's exact size (new anchor, `raref_curve.py ... smag:129550340`) from the
same spec_4 hit table. Every pre-existing grid row reproduces `rarefaction_curve.tsv` exactly (0 rows differ).
SMAG's dark gap is second only to arid (+28.38) and above deserts_combined (+24.33) and spire (+21.30).
At 129.6M proteins SMAG is about GEM's size (123.0M) but leaves 23.9 points fewer dark families absent
(40.47 vs 64.38).

| union | dark %abs | ctrl %abs |
|---|---|---|
| mgnify90+gem+deserts+spire (shipped, recomputed: matches) | 28.74 | 0.41 |
| **mgnify90+gem+deserts+spire+smag** | **26.83** (28,446) | **0.36** (385) |
| gem+spire (shipped, recomputed: matches) | 39.10 | 0.89 |
| gem+spire+smag | 33.53 | 0.66 |

Kept: `/ibex/user/velazqal/eqmags_31/smag/{smag_db.dmnd (54 GB), makedb.log, smag_db_manifest.tsv}`
(dmnd md5 `213ca2a1799ed8dd9d5c3b82c040dd66`; proteins recoverable with `diamond getseq`). Deleted:
genomes (121 GB), parts (39 GB), protein chunks and concatenated FASTA. `tmp/` is gone.
Outputs: `results/smag_{specificity_summary,content_vs_depth}.tsv` · `smag_hits.m8` / `smag_hit.ids` /
`smag.log` · `rarefaction_curve_with_smag.tsv` · `union_specificity_with_smag_summary.tsv` ·
`{dark,ctrl}_nohit_union_with_smag.ids` · `smag_db_manifest.tsv` · `smag_supplementary_data_2.xlsx`.
**Open:** whether SMAG enters the headline union is a decision to make; the numbers are not registered in
`CATALOGUE_NUMBERS.md`. Depth leg rerun on the with-SMAG sets (52498970): `results/depth_robustness_summary_with_smag.tsv`,
`depth_mgnify2026_07_hits_with_smag.tsv`, `depth_hit_with_smag.ids`, `dark_nohit_plus_deepmgnify_with_smag.ids`;
work dir `/ibex/user/velazqal/eqmags_31/depth_with_smag/` (29 MB). `spec_5` gained the `TAG` switch and the
query-FASTA build step, which had been done by hand before.

## Read-mapping audit (2026-08-17, prompted by Alejandra) — nothing to recompute
Three mapping products had been conflated. `M/scripts/` excluded as superseded PI code.

| product | rule | used by | verdict |
|---|---|---|---|
| `fam_presence.tsv` (536,437 dark genes) | `covf ≥ 0.5` any member | §1.4 persistence | **biased** |
| `fam_abundance.tsv` (same run, TPM) | none, length-normalised | §2.5 | clean |
| `profiles_full_recompute/*/abund.tsv` (5,438 rep genomes) | rel-abund > 0 ⇒ CoverM's **unset 10% default** | §2.1–2.4 | clean (tested) |

Bias **demonstrated** (`presence_rule_length_test.py`): swapping the rule reverses the length–occupancy
relation — `covf≥0.5` 146→25 aa, `covf≥0.9` 147→23 aa, `TPM>0` 76→175 aa. So any distribution claim from
that matrix must be **length-stratified**.

**§1.4 is the whole blast radius and the claim SURVIVES** (`s14_persistence_by_length.py`): median
samples 201→107→87→71→**55** and all-5-trip persistence 89.7→81.9→78.4→73.6→**67.1%** across <60aa→≥250aa;
≥4 trips ≥90.8% everywhere; all 4 single-sample families are <60 aa, **none ≥60 aa is single-sample**.
Stronger stated that way. §1.4's universe, recurrence, dN/dS and phylum figures are assembly/annotation
based → untouched.

**Four writing tasks, no fixes:**
1. §1.4 — present persistence length-stratified (98 spans 201–55).
2. §2.5 — wording only: `main.tex:197` "detectable" and `supplementary.tex:285` "present" both describe a
   gate that is `TPM>0` in code (`eco_1:32`, `eco_3:80`); M.9's "non-zero" is correct.
3. **M.5 is a `\todo` stub** (`main.tex:485`) governing all of §2.1–2.4. The **literal argv**, verified at
   `rebuild_20260729/scripts/run_current_coverm_profile.py:193-210` — nothing else is passed:
   `coverm genome --reference <prepared>/reps.mmi --minimap2-reference-is-index --genome-definition
   <prepared>/genome_def.tsv --coupled <R1> <R2> -m relative_abundance -p minimap2-sr -t 16`, CoverM 0.7.0.
   **FOUR consequential parameters are therefore INHERITED, not chosen:** `--min-covered-fraction` (10%,
   the detection floor), `--min-read-percent-identity` (0), `--min-read-aligned-percent` (0) and
   `--contig-end-exclusion` (75 bp).
   **CORRECTED 2026-08-17** (caught by session 7a): an earlier version of this block listed
   `--contig-end-exclusion 75` among the *passed* flags. It is not passed —
   `grep -n 'contig.end.exclusion' run_current_coverm_profile.py` returns nothing. Presenting a CoverM
   default as a deliberate setting would have been a false statement in Methods.
   Reference `coverm_current/prepared/reps.mmi` from `04_functional/refs/reps.fa` (the 5,438 reps).
   **Cohort, reconciled (corrected 2026-08-17 by session 7a — the arithmetic now closes exactly):**
   573 dirs under `profiles_full_recompute/`, of which **4** are abandoned `.tmp` staging dirs —
   `.51Sr1.tmp.50030409.0`, `.F5Dr1.tmp.50063405.0`, `.F13PRr1.tmp.50051922.0`,
   `.SPositive.tmp.50089048.0` — so **573 − 4 = 569 carry `run_provenance.json`**, and **548** are the
   usable canonical soil cohort. **Quote "548 usable of 569 profiled"; 573 is a directory count, not a
   cohort.** (I first found only 2 because I tested for a missing `abund.tsv`, and `SPositive` and `F5Dr1`
   left a non-empty partial one; enumerate `-name '.*tmp*'` instead.)
   **Keep two counts distinct, they are both correct about different things:** **four** failed CoverM array
   tasks in the PI's chain (the four staging dirs) versus **three** failed paper-facing samples
   (`51Sr1`, `F5Dr1`, `F13PRr1`, per `CORRECTIONS.md`) — `SPositive` is a positive control that was never
   in the canonical soil cohort. Do not "fix" one into the other.
   All **569/569 are FRESH runs** — every provenance file reads
   `profile_source = new_prebuilt_index_mapping_full_cohort_after_equivalence_rejection`. The script does
   contain a `validated_legacy_exact_reference_reuse` path (lines 181–191, sha256 fail-closed), but for
   this cohort the equivalence check was **rejected and nothing was reused**, so M.5 may state these as
   fresh recomputes.
4. §2.1–2.4 — hedge **absence** wording: "not detected" = below a 10% breadth floor, not absent
   (e.g. §2.3's "Halobacteriota detected in only 36 of 248"). A `DETECT2 = 0.01` sensitivity variant
   already exists in the community scripts.

Also: a `CORRECTIONS.md` entry for the presence-rule bias — **proposed, not written unilaterally**.
**⚠ Exception to "live code is only in `*_clean`":** this rep-level run is the **PI's**, in
`rebuild_20260729/.../coverm_current/`, with no `*_clean` equivalent. Don't rebuild or tidy it.

## Decisions with Alejandra (do not re-litigate)
1. Scope = core recompute + SPIRE + rarefaction + B6. Reciprocal recruitment and the per-biome density
   baseline **deferred by decision**, not undone work.
2. MGnify depth: the 718M DIAMOND index is the **primary rule**; the 1.66 G 2026_07 build is a separate
   **robustness leg only** (MMseqs2 ≠ DIAMOND — never a `tab:specificity` row, never unioned).
3. Display: rebuild `tab:specificity` and `fig:specificity` as-is; rarefaction → supplementary.
4. **Section REFRAMED** exceptionalism → reference-coverage (approved; recorded in `PAPER_STRUCTURE`).
   Aplakidou 2025: 63.4% of 608,258 novel families are single-ecosystem, so restriction is the NORM.
   **Provenance-level wording only** — "no Arabian or Saharan site in any named desert comparator",
   NEVER "absent from all catalogues" (MGnify90 biome labels carry no geography).

## Open
- ✅ all five jobs harvested and registered (2026-08-17).
- ❓ **Needs Alejandra:** (a) `fig:endemism` occupancy panel — recommend replacing it with the
  length-stratified absence-vs-occupancy comparison (`b6_absence_occupancy_by_length.tsv`); amends
  decision 3, record in `PAPER_STRUCTURE` if approved. (b) shipped §1.4 — recommend no edit now, log and
  raise with the PI.
- ✍ then draft §3.1 + Methods **M.12** (pushes Statistics/Availability to M.13/M.14) + rebuild displays.
- 🧾 **new bibkey `bay2021aridity`** — `desert_autotroph` = SPIRE `Bay_2021_autotrophy_desert` = Bay SK
  et al. *ISME J* 15:3339–3356 (2021), doi:10.1038/s41396-021-01001-0, PRJNA657906. **NOT `bay2021trace`**
  (*Nat Microbiol* 6:246–256) — that would be a mis-citation. Three of five deserts are unpublished JGI
  studies → Methods wording "five desert metagenome studies assembled in SPIRE \citep{schmidt2024spire}".
  Accessions: atacama PRJNA484015 · autotroph PRJNA657906 · mojave PRJNA5676xx/5678xx · arid
  PRJNA442xxx/466xxx/502xxx/518444 · antarctic PRJNA364xxx/366xxx/405xxx/406xxx.
- June's MGnify used `-k 25`, ours `-k 1`: absence fractions fine, but **June's biome-of-hits breakdown
  (83.5% aquatic) is NOT recomputable** from our output.

## B6 occupancy — done, and the occupancy × dN/dS trend is NOT reportable
`b6_occupancy_join.py` → `results/b6_*.tsv`. Definitions differ 7.4×: read-mapping presence puts 1,690
families (1.59%) at ≤3 sites, `fam_recurrence` assembly detection puts 12,523 (11.81%). The
occupancy × dN/dS trend (0.398→0.493) is **dropped**: dN/dS coverage falls 96.0%→59.7% across the same
bins because `fam_6c_dnds.sh` drops short CDS, and the top bin is made of short families. Independently
reinforces the `RQ/B6/result.md` walk-back. **Usable result:** reference-absent families are **not local
rarities**. **Do not call it a null** (same direction every stratum ⇒ a rank test will be significant);
say "marginally less widespread, never confined".

**CORRECTED 2026-09-17.** The b6 outputs had been written before the SPIRE leg landed, so they split the
universe 35,158 / 70,884, the pre-SPIRE union, while `dark_nohit_union_all.ids` already held the
corrected 30,479. Rerun of `b6_occupancy_join.py` on the same inputs (551 columns, 61 sites, 106,042
families, 30,479 absent):

| class | n | median sites | median samples | % ≤3 sites | median dN/dS |
|---|---|---|---|---|---|
| absent from all references | 30,479 | **37** | 95 | 1.59 | 0.4718 |
| has a hit somewhere | 75,563 | 38 | 100 | 1.59 | 0.4232 |

Length strata (absent vs has-hit median sites): 44/50 · 38/40 · 33/36 · 31/32 · 26/27. So the paper's
"36 of 61 against 38" becomes **37 against 38**, the gap narrows, and the conclusion is unchanged and
slightly stronger. Definition check unchanged: read-mapping presence 1,690 families (1.59%) at ≤3 sites
against 12,523 (11.81%) by assembly detection. Not yet registered.

## Deferred by decision 2026-08-17 (not undone work)
- **Reciprocal hot-desert recruitment.** Local desert data is 16S amplicon only (`rob/deserts` = Negev
  PRJNA642232 + Mojave PRJNA544067; `rob/crossdesert/du2025` = PRJNA1137527). Revival target =
  **Leung et al. 2025, PRJNA832417** (911 MAGs incl. Namib/Negev/Kuwait/Sahara — Kuwait is the nearest
  Arabian-Peninsula neighbour in any public genome-resolved set); SPIRE-registered but genecall tar 404s.
- **Per-biome dark-matter density baseline** (EQ vs marine/gut per Mb CDS at matched completeness) —
  ~1 day; the analysis that would have let the OLD heading stand. Deferred with the reframe.

## Files
- `results/query_{106042,ctrl,all}.faa` · `{catalog}_hits.m8` / `_hit.ids`
- `results/{specificity,desert_specificity,union_specificity}_summary.tsv` · `{dark,ctrl}_nohit_union_all.ids`
- `results/rarefaction_curve.tsv` (50613400) · `depth_robustness_summary.tsv` (50618838)
- `results/b6_*.tsv` · `presence_rule_length_test.tsv` · `s14_persistence_by_length.tsv` · `mag_detection_size_test.tsv`
- `scripts/spec_{1..5}_*.sbatch` · `spec_6{a,b,c}_smag_*.sbatch` · `raref_curve.py` (optional `name:N` anchors, 2026-09-24) · `b6_occupancy_join.py` ·
  `presence_rule_length_test.py` · `s14_persistence_by_length.py` · `mag_detection_size_test.py`
- Off-volume (own quota, project volume is 99% full — see `project_c2014_storage_pressure`):
  `/ibex/user/velazqal/eqmags_31/{spire,raref,depth}/`

## Basis
`../darkfam_basis_clean/results/` — `FAMILIES_FINAL.txt` (106,042) · `CONTROL_FAMILIES.txt` (106,042) ·
`fam_rep_seq.fasta` / `ctrl_rep_seq.fasta` · `fam_presence.tsv` (551 samples) · `selection/dnds.tsv`
(dark = 97,200).
