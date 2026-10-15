# Dossier 03 — Structural novelty (§1.3, Methods M.6)

Read-only audit, 2026-09-17. Sources: `paper/main.tex` (§1.3, M.6, `fig:structural`), `notes/CATALOGUE_NUMBERS.md` §1.3, `notes/CORRECTIONS.md`, and workbench `M = /ibex/scratch/projects/c2014/rob/metagenomics`, folder `M/04_structure/00_structure_clean/` (README, every script, `05_validate/`); `M/rebuild_20260729/structural_calibration/` read for comparison only.
Tags: **[V]** re-derived from source tables this session · **[reg]** registered value as is · **[new]** computed here, in no registry, not paper-ready.

## Question
Do the sequence-dark proteins of the genome-resolved catalogue adopt new three-dimensional architectures, or known folds diverged past sequence recognition — and what error rate does the procedure that answers this commit on proteins already known?

## Basis — every filter in order
| # | filter | count | source |
|---|---|---|---|
| 0 | protein representatives (catalogue) | 43,423,322 | §1.2 |
| 1 | sequence-dark AND encoded in ≥1 HQ/MQ bin = genome-resolved dark | **4,292,446** (9.89 %) | `darkmatter_clean/05_mag_dark/mag_dark_rep_ids.txt.gz` |
| 2 | ∩ length ≥100 aa → **the query set** | **1,466,943** [V] | `03a_prep_query.sh` (`06_dark_lengths/dark_rep_len.tsv.gz`, `awk '$2>=100'`) → `03_search/target_ids.txt` |
| 2b | 3Di provenance of the 4.29 M: reused tier1/tier2 vs newly predicted | 2,898,693 / 1,393,753 | `01_prep/` |
| 3 | no hit in **either** AFDB v6 or PDB (triage) | **997,646 = 68.0 %** [V] | `03_search/SUMMARY.txt` |
| 4 | Foldseek structure clustering of those 997,646 | 901,638 clusters | `04_cluster/SUMMARY.txt` |
| 5 | admitted at ≥3 members (each already a 90 %-id rep) | **17,221 families**, 80,205 members | same; excluded: 851,393 singletons + 33,024 pairs, no claim |
| 6 | one model per family, ESMFold-folded or PI model reused | **17,221** = 13,770 + 3,451 [V, PDBs on disk] | `05_validate/SUMMARY_05a.txt` |
| 7 | gate length ≥100 aa | 17,218 of 17,221 [V] — nearly free; confidence is the real gate | `per_rep_structures.tsv` |
| 8 | gate mean pLDDT ≥0.70 → **Q** | **2,838** (16.5 %) [V] | `summary_05c.json` |
| 9 | dual no-hit, query-normalised TM <0.5 vs both refs → **K** | **1,031** [V] | `summary_05c.json` |

**Selection note [new].** At step 6 `05a_prep.sh` prefers any family member that already has a PI ESMFold model over the cluster representative. Those 3,451 reps therefore come from the PI's prior novel-fold cohort — a novelty-enriched pick. Measured: K/Q = 165/386 = 42.8 % on REUSE reps vs 866/2,452 = 35.3 % (CI 33.5–37.2) on our own folds. Dropping REUSE entirely still gives ~10.5× the floor, so the headline survives; the asymmetry is undisclosed in M.6.

## How it ran
| stage | tool + build | parameters as run |
|---|---|---|
| 3Di | ProstT5 `prostt5-f16.gguf` via Foldseek **GPU build `8dc75c74`** | `createdb --prostt5-model … --gpu 1`; no length or batch flags |
| triage search | same GPU build, 215 chunks of ~7 k | `--alignment-type 2 -s 9.5 -e 1e-3 --max-seqs 50 --db-load-mode 2 --split-memory-limit 650G` |
| clustering | Foldseek `cluster`, same build (**not** MMseqs2) | `--single-step-clustering --cluster-mode 2 -s 4`; ≥3-member admission applied in shell |
| folding | ESMFold `facebook/esmfold_v1`, HF `transformers` 4.40.2 / torch 2.2.2 [V from the `esmfold` env] | half-precision trunk, `set_chunk_size(64)`, sequences truncated to 1,022 aa, no explicit `num_recycles` (config `max_recycles: 4` [V]) |
| validation screen | Foldseek **packaged `10-941cd33`** (`module load`) — the only step on this build | `--alignment-type 1 -e 10 -s 9.5 --max-seqs 50 --db-load-mode 2 --split-memory-limit 650G`, format includes `alntmscore,qtmscore,ttmscore` |

References at both stages: **AFDB v6 / UniProt 2025_03, 241,070,489 entries** and **PDB release 250101, 324,204 records**; the triage's `fsdb_pdb` and the validation's `pdb_20260730` are the same snapshot (`20260730` is a freeze date). `afdb_pad` (v4) is used nowhere. Two TM normalisations are recorded per hit: query-normalised (qTM) and alignment-normalised (alnTM). Gates: length ≥100 aa and mean pLDDT ≥0.70, applied identically to both cohorts; a model returning no row scores zero rather than being dropped.
**Citable lineage:** 50342556 + requeue 50542223 (215/215) → 03c 50559990 → cluster 50560204 → prep 50573919 → ESMFold 50573938 + 50589898 (28/28) → validation 50668615 → integrate 50685374; all nine have logs on disk [V]. **Never cite:** 50253325, 50560000/11/36, 50560586, 50574010, 50603505, 50603687, 50603810 (the timed-out single-job screen; only its PDB leg was salvaged).

## Controls — the critical part
**Origin.** `rebuild_20260729/scripts/24_prepare_known_fold_controls.sbatch` → `prepare_pdb_positive_controls.py`. The pool is the frozen PDB snapshot itself via `foldseek convert2fasta`: 324,204 records → 233,144 eligible unique chains after removing non-canonical residues (8,419), length outside 30–4,096 aa (53,922), duplicate sequences (28,719) and sequences identical to a target (0). Controls are matched to **the PI's 6,456 prior dual-novel targets — not to our 17,221** — on length bin, entropy bin and Prodigal partiality class, seed 20260730, cap 4,096; 3,684 matched: 1,777 complete chains + 771 N-terminal + 807 C-terminal + 329 internal crops (= 1,907 crops), reproducing this catalogue's fragmentary gene calls. Crops are deterministic in the parent and the target length.

**"Known-fold" basis.** Membership of the reference itself, quoted from `control_match_summary.json`: *"Every control is a complete sequence or deterministic fragment of a chain present in the frozen experimental PDB structure DB."* Not a similarity call. (A weaker sibling script, `match_known_fold_controls.py`, defines known folds by ≥40 % identity / ≥80 % qcov / e ≤1e-10; it did **not** produce this manifest, which carries `pdb_structure_id`/`crop_start_0based`/`control_mode` — only the PDB-positive script writes those.)

**Blinding.** Renamed `KF000001…KF003684`, mapping held back; selection uses only length, composition and target partiality, with no predicted structure or search outcome in existence at selection time (`"selection_outcome_blind": true`). Folded by `fold_calibration_controls.py` with ESMFold `esmfold_v1`, half-precision trunk, chunk 64, `num_recycles=4` — same predictor and, since our config default is `max_recycles: 4` [V], the same recycle count as our candidates.

**Co-searched in one run under one command — VERIFIED by content.** `05c_screen.sh` builds a single query directory `05_validate/query_pdb/` = 13,770 folds + 3,451 reused + 3,684 `KF*.pdb` = **20,905 files, present on disk** [V]. `05c_afdb.sh` slices that directory `NR % 4`, so controls and candidates interleave across all four AFDB parts under one `module load foldseek/10-941cd33` and one parameter string; the PDB leg is the single banked run of 50603810 before it timed out, i.e. all 20,905 in one `easy-search`. Both banked outputs contain both cohorts [V]: PDB leg 3,570 controls + 16,399 candidates with ≥1 row, AFDB leg 3,522 + 15,846.

**3,684 → 3,013.** The gate applies to the *model*, not the sequence: 671 control models fall below pLDDT 0.70 or 100 aa, exactly as 14,383 of the 17,221 candidate models do. `05_integrate.py` takes control length/pLDDT from the manifest (`fold_length`, `fold_mean_plddt`) and candidate values from our pLDDT shard tables, or re-measured from the model file (CA count, mean CA B-factor) for reused models; both are ESMFold per-residue pLDDT means.

**Identically computed — yes at code level.** `classify()` runs twice, candidates then controls, on the same best-score dictionaries, the same 0.5 threshold and the same dual-no-hit predicate. Floors: 101/3,013 = 3.35 % (Wilson 2.77–4.06) qTM, 68/3,013 = 2.26 % alnTM [V].

**Judgment — SOUND, with one bounded weakness.** Better than anything in the comparator literature: real PDB chains, fragmentary in the same way as the queries, blinded before folding, folded by the same predictor, searched in the same job under the same command. It correctly excludes prefilter misses, TM-threshold failure on genuine homologs, and the no-row artefact (66 of the 101 false calls are no-row [V]). It does **not** exclude confidence-composition bias: gate-passing controls are much more confident than candidates (mean pLDDT 0.867 vs 0.775; mean length 308 vs 246 aa [V]), and the control false-call rate falls steeply with pLDDT — 5.49 % (0.70–0.75), 7.32 %, 4.20 %, 3.30 %, 1.53 % (≥0.90) [new]. Standardising the floor to the candidate pLDDT distribution gives **5.54 %** (→ 6.6×); jointly with length, **4.49 %** (→ 8.1×) [new]. The conclusion is untouched; the multiplier 10.8× is not. The PI's own files ask for exactly this adjustment (quoted below).

## Numbers
| quantity | value | source |
|---|---|---|
| queries ≥100 aa | 1,466,943 | `03_search/target_ids.txt` [V] |
| triage AFDB-any / PDB-any / both / neither | 463,316 (31.6 %) / 232,793 (15.9 %) / 226,812 (15.5 %) / **997,646 (68.0 %)** | `03_search/SUMMARY.txt` [V] |
| admitted fold families | 17,221 (80,205 members) | `04_cluster/SUMMARY.txt` [V] |
| confidence-qualified **Q** | **2,838** of 17,221 | `summary_05c.json` [V] |
| **K**, qTM dual no-hit | **1,031 = 36.3 %**, 95 % CI 34.6–38.1 | `summary_05c.json`; Wilson reproduced [V] |
| control floor, qTM | **101/3,013 = 3.35 %**, CI 2.77–4.06 | `per_control_structures.tsv` [V] |
| excess, significance, floor-subtracted | **10.8×**, two-proportion z p = 1.7 × 10⁻²²³, ≈936 models | reproduced [V] / [reg] |
| alnTM dual no-hit | **95 = 3.35 %** vs floor 68/3,013 = 2.26 %; **1.5×**, p = 0.011 | reproduced [V] |
| locally familiar / covered end-to-end (of Q) | **96.7 % / 63.7 %** — the complements of 3.35 % and 36.3 % | `per_rep_structures.tsv` [V] |
| K vs other confident models: mean length | **300.6 vs 214.7 aa** (medians 265 vs 216) | same [V] |
| K vs other confident models: mean pLDDT | **0.7626 vs 0.7825** | same [V] |
| gate-passing controls: mean length / pLDDT | 308.3 aa / **0.8674** | `per_control_structures.tsv` [V] |
| member proteins in the K families | **4,880** | `fam_reps.tsv` [V] |
| no-row vs retrieved-below-threshold | K: 22 / 1,009 · floor: 66 / 35 | both tables [V] |

## The two normalisations
Foldseek scores one alignment twice. **qTM** divides the TM sum by the *query* length: "does a known fold cover this protein end to end?" — an unmatched half of a two-domain protein pulls the score under 0.5 even when the matched half is textbook. **alnTM** divides by the *aligned* length: "does any aligned region have a known fold?", blind to what was left out. alnTM ≥ qTM on all 20,905 × 2 rows [V], so the 95 alnTM-novel models are a strict subset of the 1,031. The 10.8× counts *architectures without a covering counterpart*; the 1.5× counts *proteins with no locally familiar piece at all*; 936 candidates sit between the two — locally recognisable, globally unmatched. That gap is the paper's strongest defensible claim (novelty is **architectural, not local**, consistent with recombination of known domains rather than new topology) and is why "N novel folds" is a property of the chosen normalisation, not of a proteome.

## Armour not yet run
| item | already on disk | missing | rough cost |
|---|---|---|---|
| Domain calling (Merizo/Chainsaw) on 2,838 + 3,013 | all 20,905 model PDBs in `05_validate/query_pdb/` (+ `models/fold`, `models/reuse`) | no merizo/chainsaw/DPAM/UniDoc anywhere under `rob/`; needs install + weights | install ~1 h, inference GPU-hours. **The one that can move the result**: it turns "architectural novelty" from an inference about length into a measurement |
| Foldseek vs ESM Atlas (third reference) | both Foldseek builds; the query models | **no ESM Atlas / highquality_clust30 structure DB on this system**; `rebuild_20260729/references/` holds only gmsc_v1_90aa, gtdbtk_r232, gunc, pfam_38_2, and `protein_novelty/mgnify_2026_07` is a *sequence* DB | multi-TB download with scratch at 87 % — a >1-week infrastructure item, not a quick run |
| Re-query the "223 published novel folds" | nothing | **223 appears only in two armour rows of `notes/PAPER_STRUCTURE_PROPOSED.md`** — no source file, no citation, no set on disk | undefined until the source is named |
| Pfam / InterPro on the K members | Pfam **38.2** complete (`references/pfam_38_2/`, `PFAM_38_2_REFERENCE_COMPLETE`, `hmmer_version.txt`); membership in `04_cluster/clusters.tsv`; sequences in `darkmatter_clean/02_uncovered_reps.faa.gz` | a 4,880-sequence extraction; no InterProScan | hours. Cheapest, highest value: a Pfam hit on part of a K model is direct evidence for the multi-domain reading |
| AF2/AF3 cross-architecture refold | **a 200-protein AF2-vs-ESMFold pilot already exists**: `04_structure/pilot_af2/` (`pilot_200.fasta`, 202 MSAs from the ColabFold remote MMseqs2 API), `scripts/pilot_af2/01_curate…04_compare.py`, `colabfold-env` | the AF2 leg runs on ORIX, not IBEX; no local MSA DB (uniref30/envdb/BFD) | K = 1,031 is small: MSAs ~1–2 h by remote API, AF2 a few ORIX GPU-days. The obvious revision request |

## Caveats (quoted)
- §1.3, live: *"They are reference-relative candidates under an explicit search rule, not experimentally validated new folds, and they are reported against controls that fail at 3.35\% on known folds. Structure alone does not say what these proteins do."*
- M.6, live: *"A model for which the screen returns no row is scored zero rather than dropped, so an absent result counts as novel rather than shrinking the denominator."*
- Fenced B1a: *"Candidate and control cohorts are not covariate matched for partial status, length or prediction confidence, so their descriptive rate differences are not adjusted effect estimates. AlphaFold-DB community datasets and complexes also lie outside the v6 UniProt monomer reference"*.
- PI, `control_current_structures_20260802/claim_boundary.txt`: *"This is a known-fold false-nohit calibration pool. Apply matched strata before interpreting recurrent-family dual no-hits."*
- PI, same folder `summary.json`: *"Candidate-specific release still requires matching these controls by partial status, length and prediction confidence."*
- PI, `current_structural_novelty_20260730/…summary.json`: *"A dual current-reference no-hit is a calibrated candidate fold, not proof of a new fold. … Experimental structure is required for confirmation."*
- `notes/DISCUSSION_NOTES.md` §1.3: *"ESMFold is single-sequence … one representative per family speaks for its members; the controls were length/partiality-matched to the PI's candidate distribution, not re-matched to ours; the ≥3-member family admission is our choice and the pair-family sensitivity (+33,024 families) was not run."*
- `00_structure_clean/README.md` defect 7: *"224 of the 3,451 reused reps have empty length/mean_plddt in the PI's per-candidate table"* — now measured from the model file instead.

## Status
Computed and shipped. §1.3 is written with all six `\todo` placeholders resolved, M.6 needed no change, the abstract was moved off the historical numbers, the B1a scaffold is fenced. The workbench chain is closed and self-contained (9.8 GB; all 215 chunk hit dirs and both banked screen legs intact, plus `/ibex/user/velazqal/sc5c_parts/`). Two records are stale: `notes/PIPELINE_GUIDE.md` §8 is still headed *"⚠ PROVISIONAL — stage 05c is still running (50603810)"*, and `fig:structural` exists (`paper/figures/fig_structural_calibration.pdf`, generator `code/figures/make_structural_calibration_figure.py`) but is **never `\ref`-ed in live text** [V]. That generator reads the workbench absolute path directly, so it cannot rebuild away from Rob's scratch.

## What can and cannot be claimed
**Can:** 36.3 % of 2,838 confident models have no covering match in AFDB v6 or PDB under qTM, against a 3.35 % floor measured on blinded known folds in the same run — an order-of-magnitude excess, p ≪ 10⁻²⁰⁰; under alnTM the excess is 1.5×. Therefore novelty here is architectural rather than local: 96.7 % of confident models contain a locally familiar region, 63.7 % are covered end to end, and candidates are longer at indistinguishable confidence. A single "N novel folds" figure is a property of the normalisation, not of a proteome — the section's strongest general contribution.
**Cannot:** that these are new folds; that any individual model is correct; that the count is biome-specific; that 10.8× is a covariate-adjusted effect (it is a crude ratio of two pooled rates on cohorts that differ in confidence); that the 68.0 % triage means anything as a fold count; anything at all about function. Never quote 608/1,875, 10,007 or 89.2 % as current: 608/1,875 is the superseded historical cohort, alive only in the explicitly historical table and the fenced B1a paragraph.

## Reviewer attack surface
1. **Controls are more confident than candidates** (pLDDT 0.867 vs 0.775). Standardising the floor to the candidate pLDDT distribution moves it to 5.54 %, or 4.49 % jointly with length [new] — 10.8× becomes 6.6–8.1×. The PI's own claim-boundary files demand this adjustment. Prepare the number before submission; it is one pass over two TSVs already on disk.
2. **Controls were matched to the PI's old cohort, not ours.** Only 275 of his 3,684 matched pairs share a pLDDT bin; `matched_calibration_by_stratum.tsv` shows he built the stratified machinery and we did not use it.
3. **Novelty-enriched rep choice.** 3,451 of 17,221 reps were chosen because they already had a PI model from his prior novel-fold set; K rate 42.8 % vs 35.3 % for ours [new]. Undisclosed in M.6.
4. **Single-sequence ESMFold.** A low-pLDDT model is not evidence of novelty, and the candidate rate falls monotonically with confidence: 45.3 % at pLDDT 0.70–0.75 down to 15.8 % above 0.90 [new]. The AF2 refold is the rebuttal.
5. **Reference staleness.** PDB release 250101 leaves ~20 months of depositions outside the reference at submission; AFDB v6 excludes community datasets and complexes.
6. **One representative speaks for a family**, and the ≥3-member admission is our choice — the 33,024 pair families were never tested.
7. **Truncation.** 290 reps were folded truncated at 1,022 aa; 14 reach Q and 9 reach K [new] — small, but a truncated query cannot be covered end to end by construction.
8. **"324,204 entries"** in M.6 is the chain/record count of the Foldseek PDB database, not PDB entries.

## Open questions for Alejandra
1. Run the pLDDT-stratified floor and report the adjusted multiplier (≈6.6–8.1×) in §1.3, or hold it as reviewer-response armour? Recommendation: report it — cheap, pre-empts the PI's own stated caveat, and the claim survives it.
2. Disclose the REUSE-rep selection and the FOLD-only sensitivity (35.3 %) in M.6?
3. What is the source of the "223 published novel folds"? It is unsourced in the repo and the armour item cannot be planned until it is named.
4. `fig:structural` is uncited in live text — cite it in §1.3, or leave it to the figures session?
5. Pfam 38.2 on the 4,880 K members is hours of work and would turn the multi-domain reading from inference into evidence. Go / no-go — and does it count as touching Rob's tree (it reads his files, writes only into yours)?
6. Is the AF2 pilot on ORIX still live, and do we extend it to the 1,031?
