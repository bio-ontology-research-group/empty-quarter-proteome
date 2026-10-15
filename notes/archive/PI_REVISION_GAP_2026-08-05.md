# PI revision gap analysis (2026-08-05)

What the PI changed in the manuscript + repo while Alejandra worked separately. Baseline = paper `3e479b7`
(her WIP: 3-part skeleton + draft 1.1) → PI's `origin/main` (paper `de6164c`; repo main). **Her work is on
branch `ale` and is untouched; his is on `main` — they diverged, nothing overwritten.**

## Bottom line
He **rebuilt** the paper: different catalogue, different GTDB release, reverted section structure, de-hyped
("bounded/calibrated") voice. Almost every headline number she had is replaced. His version passes an
independent reproducibility audit (PASS), but several core analyses are still running (not final for ~3 weeks).

## What changed (old → new)
- **Title/thesis:** "Known **phyla**, unknown biology" (was "organisms"). Retreats from *"the tree is ordinary"* →
  foregrounds **bounded lower-rank novelty** (order/family/genus/species gaps).
- **Catalogue:** merged **37,949 MAGs → 9,253 species / r226 / 56 phyla** → **DAS_Tool-refined 14,426 bins →
  5,438 reps / GTDB R232 / 51 phyla**. Novel ranks now first-empty-rank: **11 order / 142 family / 1,343 genus /
  3,755 species** (10 of 11 order gaps adjudicated → 0 order proposals; 1 unadjudicated → recovery). Old 15,469 /
  51.9M / 44-phyla kept only as a dated "historical R226 snapshot".
  - **Source of his catalogue:** DAS_Tool refines the two binnings per sample (drops/reassigns contigs) → 14,426 →
    dRep → 5,438. Ours skipped DAS_Tool (raw SemiBin bins → dRep → 9,253). His is the more conventional choice.
- **Structure:** reverted to **A / B / C** (dropped the 3-part narrative).
- **Dark ladder:** the single "51.9M→…→10,007 folds" chain **withdrawn**; split into two historical tables —
  MAG-Bakta (51.9M → 23.4M reps → **2.84M UniRef50 no-hit = 12.1%**) and assembly (753M → 449M reps →
  **68.6M dark = 15.3%**). 21.4%-no-hit re-attached to the **23,382 dark-orthogroup** object.
- **Structural novelty:** "10,007 validated folds" **withdrawn**; also drops the "~12 validated" idea. Replaced by
  a current-reference calibration: **608/1,875 (32.4%) qTM dual-no-hit** and **97 (5.2%) alignment-TM**, vs a
  measured known-fold error floor (3.35% / 2.26%). Survivors = "reference-relative candidate folds".
- **dN/dS:** 0.39-vs-0.29 comparator **withdrawn** (28-sample confound + darkness correlated with rate). Keeps only
  dark-only median **0.391 → 0.499 at ≥30-member floor**.
- **"Endemism" → "database-relative specificity"** everywhere. **Arms-race (B4) withdrawn** → "recognition-protein
  hypotheses". **PICRUSt 2.4× reversal withdrawn** (denominator artefact); keeps 56.7% no-KO ceiling. C1a
  open-pangenome subsection dropped.

## Three hard conflicts with our work
1. **Catalogue: 9,253 (ours) vs 5,438 (his).** Different universes. Her Part 1.1 is built on 9,253 → superseded.
2. **Structure framing opposite.** We planned "structurally shallow / functional novelty"; he makes structural dark
   matter a **co-equal pillar** (in the title) = "bounded calibrated candidates". Thesis-level disagreement.
3. **Our 4.29M / 9.89% genome-resolved dark is not in his paper.** He kept assembly-wide + historical MAG-Bakta
   ladders. The "MAG-mapped" recalibration (what `02_prostt5` feeds) wasn't adopted.
   - NB: the assembly-wide dark is catalogue-independent; **our darkclean 69.9M/15.6% is the *corrected* version of
     his 68.6M/15.3%** (we closed the 37.8M-uncovered gap) — a fix he'd want, not a conflict.
   - The genome-resolved (MAG-mapped) dark maps to **all bins pre-dRep** (so 9,253-vs-5,438 dereplication doesn't
     change it) — only **raw (37,949) vs DAS_Tool-refined (14,426) bins** matters.

## Where we agree
89%/10,007 overclaim dead both sides · structure = candidates not proven folds · historical vs current catalogues
separated · retreat from "endemism" · the 23,382 dark-orthogroup object + its 21.4% shared and correctly attributed.

## Still running (per his `takeover-status`) — do NOT treat as final for 3 weeks
573-sample CoverM / MAG-KO rebuild for C1 ("do not retain any historical C1 number"); AFDB-v6 / MGnify90 merges +
Orix folding manifest; COMEBin co-assembly high-rank recovery. A newer UniRef50-2026_02 refresh already exists
(43,409 of 1,902,233 now hit; 1,858,824 remain) — newer than the printed figures.

## The 5 things to raise before he travels
1. **Which catalogue is canonical — refined 5,438 or merged 9,253?** (Everything downstream depends on it.)
2. **A/B/C or the 3-part narrative?** Agree one skeleton + who owns which section for 3 weeks.
3. **Structural dark matter — co-equal pillar or "shallow/functional"?** (Thesis-level.)
4. **What's frozen vs still running** — what's safe to write against for 3 weeks?
5. **Do you want a genome-resolved (MAG-mapped) dark layer at all** — and if so, raw or refined bins? (Decides
   whether `02_prostt5` / the structure re-run is worth running.)

Key files: PI paper `git -C paper show origin/main:main.tex`; baseline `3e479b7:main.tex`; his review docs under
`origin/main:reviews/` (`takeover-status`, `number-verification`, `dnds-comparator-withdrawal`, `ladder-rewrite`,
`current-literature-synthesis`, the two 2026-08-04 opus audits).
