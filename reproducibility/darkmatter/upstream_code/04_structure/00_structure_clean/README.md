# 00_structure_clean — structural novelty of the genome-resolved dark matter

## STATUS (2026-09-17) — floor standardised on model confidence

Controls are higher-quality models than candidates (pLDDT 0.867 vs 0.775, 308 vs 246 aa) and their
false-call rate falls with confidence, so the pooled floor was measured where errors are rarest.
`scripts/05e_stratified_floor.py` reweights the control rate onto the candidate pLDDT (and length)
distribution → `05_validate/stratified_floor.{tsv,txt}`. 2,000 bootstraps, seed 1.

| normalisation | crude | on pLDDT | on pLDDT × length |
|---|---|---|---|
| query-normalised | 10.84x | **6.56x** (4.97–9.08) | **8.09x** (5.96–11.49) |
| alignment-normalised | 1.48x | 1.74x (1.12–2.99) | **1.38x (0.87–2.41)** |

Crude values reproduce exactly. Control false-call rate 5.49 % at pLDDT 0.70–0.75 → 1.53 % above
0.90, while 70 % of candidates sit in the two lowest bins (15 % of controls).
**Report the standardised qTM multiplier; withdraw the alnTM 1.5x**, whose interval now crosses 1.
The architectural-not-local claim is unaffected. Not registered.

## STATUS — ✅ §1.3 COMPLETE (2026-08-18, session 6)

**Result (paper-facing; registered in `notes/CATALOGUE_NUMBERS.md` §1.3):**
1,466,943 queries (≥100 aa genome-resolved dark) → triage no-hit **997,646 (68.0%)** → **17,221**
admitted fold families (≥3 members) → 17,221 models (13,770 folded + 3,451 reused) → **Q = 2,838**
confidence-qualified (≥100 aa, pLDDT ≥0.70) → **K = 1,031 (36.3%, CI 34.6–38.1)** dual no-hit under
query-normalised TM, against our own measured control floor **101/3,013 = 3.35%** (CI 2.77–4.06):
**10.8× the floor, p < 10⁻²⁰⁰**, ≈936 models after subtracting it.
Sensitivity (alignment-normalised TM): **95 (3.35%)** vs its own floor 2.26% → only **1.5×** (p = 0.011).
Reading: novelty is **architectural, not local** — 96.7% of confident models contain a matching aligned
region, 63.7% are covered end-to-end; candidates are longer (301 vs 215 aa) at equal pLDDT (0.763/0.782).
Triage retrieval: AFDB-any 463,316 (31.6%) · PDB-any 232,793 (15.9%) · dual 226,812 (15.5%) ·
AFDB-only 236,504 · PDB-only 5,981.

**Job lineage (cite these):** search 50342556 + requeue 50542223 (215/215 chunks, 1,466,943 queries
verified exactly) → 03c 50559990 → cluster 50560204 → prep 50573919 → ESMFold 50573938 + 50589898
(28/28 shards, 0 OOM) → validation search 50668615 (4-way split) → integrate 50685374.
Do NOT cite: 50253325, 50560000/11/36, 50560586, 50574010, 50603505, 50603687, 50603810 (see log below).

**References (verified by `.version` md5):** AFDB v6 / UniProt 2025_03 (241,070,489) + PDB release
**250101** — the triage's `fsdb_pdb` and the validation's `pdb_20260730` are the SAME snapshot
(`pdb_20260730` is a freeze date, not a release). `afdb_pad` (v4) is used nowhere here.
Params: triage `--alignment-type 2 -s 9.5 -e 1e-3 --max-seqs 50 --db-load-mode 2 --split-memory-limit 650G`;
validation (PI's 167_postfold, verbatim) `foldseek/10-941cd33 --alignment-type 1 -e 10 -s 9.5 --max-seqs 50`.

**Control apparatus:** the PI's **3,684 blinded known-fold models** were co-searched inside our own
validation run (20,905 queries = 13,770 + 3,451 + 3,684). Our floor reproduces his exactly — the
cross-check that this is the same instrument, and what licenses comparing our result to his cohort.

**Scripts:** `03a/03b_prep_chunks/03b_search/03c_candidates` · `04_cluster` · `05a_prep` · `05b_fold`
(+`05b_fold_v100` variant, unused) · **`05c_afdb` (array) + `05d_integrate` + `05_integrate.py`**.
⚠ `05c_screen.sh` is SUPERSEDED (single job, timed out at 24h) and marked so in its header.

**Problems hit and fixed (so they are not re-hit):**
1. 8 search chunks timed out at ~13h → requeued at 24h; stale foldseek tmp must be deleted first.
2. `03c` read `results/` but `03b` writes `hits/` → would have produced a silent empty union.
3. `04_cluster` bugs: bash `local a=$1 b=${a}` expansion order; **tier `q_t1/q_t2 .lookup` were 0-byte
   broken symlinks** → rebuilt from the header DBs with `foldseek prefixid` (extraction then verified
   exact: 580,244 + 4,712 + 412,690 = 997,646); `concatdbs` ERRORS on lookup/source maxSetId mismatch →
   subsets now drop `.lookup/.source` and stderr is no longer suppressed.
4. GPU: `--partition=gpu` is not permitted; use `--partition=batch --gres=gpu:a100:1`. Small shards
   (500 seqs / 4 h) backfill despite low fairshare; the V100 pool looked idle but its slots turn over
   so slowly that the estimated start was 6 days out — A100 was ~1 h/shard.
5. **`/ibex/scratch/projects/c2014` is 99% full** → the validation died on ENOSPC. Heavy FoldSeek work
   now runs out of `/ibex/user/velazqal` (1.5 T) with only compact outputs copied back.
6. The single-job validation needed ~36 h vs a 24 h wall → split the AFDB leg 4 ways (~9 h each); the
   PDB leg (6 min) was salvaged and not re-run.
7. **224 of the 3,451 reused reps have empty length/mean_plddt in the PI's per-candidate table** →
   `05_integrate.py` now measures both from the model file itself (CA count; mean CA B-factor = pLDDT),
   so all 17,221 reps are measured by one instrument and his table is a last-resort fallback only.

**Paper side:** §1.3 written (all 6 `\todo` placeholders resolved); Methods M.6 needed no change; the
**abstract** and the **spine figure** were still carrying the historical 608/1,875 and were updated; the
live **B1a scaffold paragraph was fenced**. Build 37 pp, 0 undefined refs/citations, 0 overfull boxes.
Open for Alejandra: the historical B1/B1a Table (captioned historical) now sits beside a current-basis
§1.3 — keep or cut is a figures-session decision; and the spine figure's panel 2 still shows the dead
23,382 (§1.4 owns it; 106,042 is current).

Clean, self-contained analysis on the **current** data. It **reuses** the existing ProstT5 3Di (no regeneration)
but touches the messy parent `04_structure/` **exactly once** (a harvest step); after that everything is clean
joins in this folder. Owner: velazqal, 2026-08. (Named `00_` to sort to the top.)

## Goal
For the **dark proteins that map to a MAG** — genome-resolved dark, **4,292,446** rep IDs — determine structural
novelty (known fold vs structurally novel) → **ESMFold-validated** novel folds.

## Established facts (this session)
- **3Di already exist for ~all MAG proteins** (ProstT5), stored as two **disjoint** FoldSeek batches:
  `../tier1_cluster/concat/g0000–g0007` (11,418,079) + `../tier2_cluster/concat/g0000–g0020` (13,088,646) =
  **24,506,725**, overlap 0. So we reuse these — no ProstT5 re-run.
- The old structural *results* are **not** reused. The predicted-3Di triage is reported as a triage
  (68.0% no-hit) and never as a fold count; the paper's number is the calibrated K/Q from step 5.

## Inputs (read-only)
| what | path |
|---|---|
| dark-MAG rep IDs (4,292,446) | `../../04_protein_catalog/darkmatter_clean/05_mag_dark/mag_dark_rep_ids.txt.gz` |
| existing 3Di (24.5M, both tiers) | `../tier1_cluster/concat/`, `../tier2_cluster/concat/` |
| reference structure DBs | AFDB **v6** `../../rebuild_20260729/structural_calibration/references/afdb_uniprot_v6_uniprot_2025_03/database/afdb` + PDB 250101 `../fsdb_pdb` (**not** `../afdb_pad`, which is v4) |
| ESMFold assets | `../../06_esmfold_novelfold/` (reuse structures where proteins overlap; fold only new ones) |

## Steps (`scripts/`)   [all complete — see STATUS]
1. **01_prep** ✅ — split the 4.29M into `reuse_ids` (**2,898,693** already have a 3Di) + `missing_ids`
   (**1,393,753**); extract the missing reps' sequences, chunked (~8×175k) for the A100 array.
2. **02_prostt5** ✅ — ProstT5 3Di + FoldSeek-vs-AFDB for the 1.39M missing (A100, ~8 chunks, ~3h wall).
   Commands copied verbatim from the proven `scripts/06d_prostt5_array.sh`. Keeps the queryDB (3Di) + `afdb_hits.tsv`.
3. **03_combine** — the structurally-novel set over all 4.29M = (reuse: AFDB **no-hit** from the tiers'
   `afdb_hit_ids`) ∪ (missing: AFDB no-hit from step 02). Reuses the tiers' deterministic AFDB result for the 2.9M
   (same AFDB, same params) → **no DB surgery, no re-search**.
4. **04_cluster** ✅ — structural clustering of the no-hit 3Di (62b params) → 17,221 admitted families.
5. **05a/05b/05c_afdb/05d** ✅ — ESMFold on family reps (reusing existing models where they exist), then
   full-coordinate TM validation with the PI's blinded known-fold controls as co-queries → calibrated K/Q.

## Results
| metric | value |
|---|---|
| dark-MAG reps (input) | 4,292,446 |
| have a 3Di → **reuse** | 2,898,693 (67.5%) |
| missing 3Di → **ProstT5** | 1,393,753 (32.5%) |
| queries searched (≥100 aa) | 1,466,943 |
| triage no-hit (candidates) | 997,646 (68.0%) — a triage, not a result |
| admitted fold families (≥3 members) | 17,221 |
| models (folded + reused) | 17,221 (13,770 + 3,451) |
| confidence-qualified **Q** | 2,838 |
| **K** dual no-hit, query-normalised TM | **1,031 (36.3%)** vs 3.35% control floor |
| dual no-hit, alignment-normalised TM | 95 (3.35%) vs 2.26% control floor |

## Rules
- Reuse **assets** (3Di, AFDB/PDB, ESMFold weights); never reuse the old **results** (tiers' fold calls, `foldseek_results`).
- One namespace (assembly `sample::contig_gene`). Nothing is ever written back into `04_structure/`.
