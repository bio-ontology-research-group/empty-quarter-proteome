# Claude Opus final repaired audit

- Date: 2026-08-04
- Reviewer: Claude Code 2.1.221, `claude-opus-5`, high effort
- Mode: read-only plan mode; no web, subagents or writes
- Audited checkout: `/tmp/eq-mags-final2.0YmB6d/repo`
- Public harness log: `/tmp/eq-mags-final2.0YmB6d/run_all.log`
- Process exit status: 0

The report below is the reviewer's verbatim final result. The six non-blocking
notes were subsequently addressed and submitted to a focused Opus follow-up
audit.

## VERDICT

**PASS**

No submission or reproducibility blocker remains. All five audit questions resolve affirmatively.

## BLOCKERS

None.

## NON-BLOCKING NOTES

1. **`reproducibility/verify.py` — `--quiet` is a no-op.** Running `python3 reproducibility/verify.py --quiet` still prints all 165 per-claim rows (exit 0). Correction: gate the per-claim `STATUS/ID` table on `not args.quiet` and keep only the `GATES:`/tier-summary lines. Cosmetic; the README does not advertise `--quiet`, so nothing is mis-documented.
2. **`docs/ibex-workdir-inventory.md:24,69,140` — unhedged "novel fold" wording.** The inventory describes `06_esmfold_novelfold` and `family6_validation` as "novel fold families" / "one novel fold family", which is the exact framing the manuscript withdraws ("reference-relative candidate models, not established new folds"; `manuscript/main.tex:243`). The guard's forbidden list catches `68,808 strict novel` but not this phrasing. Correction: write "candidate novel-fold family (historical directory name)" at those three lines, and add `"novel fold famil"` to the forbidden dict with an allowance for the hedged form. Not a blocker — it is an internal directory inventory, not a claim-bearing output, and no manuscript number depends on it.
3. **`results/global_specificity/MULTICATALOG.md` — 148M vs 147M desert-union depth.** The table row reads `ALL deserts (union) | 148M`, while the same file's closing paragraph and `manuscript/main.tex:549` (Table `tab:specificity`) both say 147M. Correction: set the table cell to 147M. Catalogue depth is a descriptive annotation, not a registered claim, and the no-hit percentages (31.9%, 7,469/23,382) are unaffected and verifier-checked (`c3.deserts` PASS).
4. **Guard is exact-substring only.** `verify_manuscript_boundaries()` (`reproducibility/verify.py:572–690`) rejects 40 literal phrases; a paraphrase of a withdrawn claim would pass. `reproducibility/README.md:8` says the verifier "rejects superseded or unbounded language", which reads stronger than the implementation. Correction: add one clause — "by exact-phrase match over a curated withdrawn-language list; it is a regression guard, not a semantic checker."
5. **RQ-level figures are outside the harness assertion set.** `regenerate_figures.sh` asserts exactly the 11 manuscript figures (matching the README's "all 11"); `RQ/C2/fig_pangenome_openness.pdf` and `RQ/F2/fig_bgc_novelty.pdf` are rebuilt only by their own `finalize.py`/`make_fig_bgc_novelty.py`. This is correct scoping (neither is a manuscript figure), but the README table does not state it.
6. **Abstract elides the one-bin gap.** `manuscript/main.tex:38–39` presents 14,426 bins and 5,438 representatives adjacently; the 5,438 actually derive from the 12,433-member historical dRep run, not the 12,434 current HQ/MQ set. This is disclosed explicitly in Methods (`main.tex:372–374`), in `docs/binning-and-annotation-scope.md:30–35`, and machine-enforced as a *required* phrase. Adequate as-is; a half-clause in the abstract would be belt-and-braces.

## EVIDENCE CHECKED

**Q1 — bounding of the four withdrawn overclaims (consistent in outputs *and* producers).** Confirmed in both directions. Outputs: `ENDEMISM.md` (occupancy labelled "within-EQ", MWU flagged unclustered), `GLOBAL_ENDEMISM.md` ("MGnify90-relative … does not establish geographic absence, desert specificity, or endemism"), `SELECTION_ADAPTIVE.md` ("exploratory unclustered output; not promoted"), `RQ/C2/result.md` ("historical exploratory output; not promoted"; accessory enrichment stated as association, not mechanism for C1), `RQ/F2/result.md` (99.3% explicitly a fragmentation-sensitive upper bound; 107,316/119,090 = 90.1% contig-edge). Producers emit the same bounds: `phase_d_global.py:47,54,63`, `phase_c_adaptive.py:92,100,105`, `phase_d_endemism.py:7`, `ibex/96_darkfam_dnds_aggregate.sh:4,42,70`, `RQ/C2/140_pangenome_openness.py:6–8`, and both figure scripts stamp on-canvas provenance footers (`finalize.py:38`; `make_fig_bgc_novelty.py:67–70`). All eight producers are named individually in the guard's `paths` list, so regenerating them cannot reintroduce banned strings.

**Q2 — can a producer regenerate withdrawn language?** No. `RQ/_gen_answered.py` carries bounded text for B2/B6/D5/E1/H3 and bounded stubs for F2/G2/H1; its `NOTNOW` loop skips any folder where `result.md` already exists, so the richer manual F2/G2/H1 records are preserved rather than reverted to stubs. `RQ/C2` is in neither dict and is never machine-written. `_gen_answered.py` is itself swept by the boundary guard.

**Q3 — dRep reconciliation.** Internally coherent and machine-checked, not merely asserted. `drep_membership_reconciliation.tsv` records 12,434 / 12,433 / 5,438 / current-only 1 (`25PRr2__self_SemiBin_2_sub`) / historical-only 0, each with its upstream path and SHA-256 (manifest `613e060a…`, Cdb `43f0625a…`, Wdb `393238e1…`; set-comparison rows correctly carry `NA`). `sha256sum -c` over all 13 entries of `current_manuscript_artifacts.tsv` returned OK, including the reconciliation file's own digest. Verifier claims `scale.refined_hq_mq`, `scale.drep_inputs`, `scale.drep_current_only`, `scale.drep_historical_only`, `scale.representatives` all PASS, and the non-closure sentence is a *required* manuscript phrase.

**Q4 — public clean-checkout harness.** `/tmp/eq-mags-final2.0YmB6d/run_all.log`: figures rc=0, verify rc=0, manuscript rc=0, "MANUSCRIPT: OK (20 pages; BibTeX resolved; public snapshot)", zero `FAIL`/`MISS`/`Traceback`. `paper/` is absent and the guard tolerates it (`os.path.isfile` filter). `diff status.before status.after` is empty — the harness is idempotent and produces no tracked-file churn. My own `verify.py` run reproduced GATES PASS, 13 frozen artifacts, `COMPUTE: 62 | FROZEN: 72 | TABLE: 6 | PROVENANCE: 25 | FAIL: 0`, exit 0. Tier semantics match implementation: TABLE rows print as `SOURCE … [TABLE] <- path` (presence-only, as README states), PROVENANCE rows as `LEDGER … [PROVENANCE]`; 62+72+6+25 = 165 = `claims.tsv` body rows.

**Q5 — manuscript arithmetic, denominators, labels, caveats.** Recomputed: 4,225+8,209+1,992 = 14,426; 2,259+3,179 = 5,438; 608/1,875 = 32.43%, 97/1,875 = 5.17%, 101/3,013 = 3.35%, 68/3,013 = 2.26%; 2,839,030/23,416,070 = 12.1%; 5,011/23,382 = 21.4%; 10,669/23,382 = 45.6%, 543 = 2.3%; order gaps 11 = 10 staged + 1 unadjudicated. Cross-file consistency holds: `SELECTION_ADAPTIVE.md` 9,425+13,342 = 22,767 = `b2.dark.n`; `GLOBAL_ENDEMISM.md` 14,830+8,552 = 23,382; F2 42,064+301 = 42,365 and 119,090−107,316 = 11,774. The three universes stay separated and labelled throughout — current 14,426/5,438 R232 catalogue, historical 15,469-MAG R226 functional snapshot, and the parallel 23,382-orthogroup object explicitly stated (`main.tex:249–254`) as *neither* the 2,839,030 representatives *nor* a subset of the 334,810 families. Inferential caveats are present at every headline: dN/dS as coding-like constraint not purifying selection, no-hits as database-relative not geographic absence, empty ranks as candidate depths not validated taxa, PICRUSt2 ratios shown under both denominator rules (2.40/1.11 vs 31.41/14.63) with neither promoted as enrichment. Explicitly pending items (G2 biome-scale SMETANA, H1 decontam/relic-DNA, C1 CoverM recomputation, F2 `--complete` rerun) are labelled as such and not used as headline evidence — not counted against this verdict.

**Files read:** `manuscript/main.tex` (abstract, B2–B4, Methods, Tables `tab:scale`/`tab:specificity`/`tab:decode`, pipeline figure); `reproducibility/{README.md,verify.py:570–700,claims.tsv,current_manuscript_artifacts.tsv,drep_membership_reconciliation.tsv,regenerate_figures.sh,REPORT.md,canonical_catalogue_summary.json}`; `docs/{binning-and-annotation-scope.md,ibex-workdir-inventory.md}`; `RQ/{README.md,_gen_answered.py,B6,C2,D5,E1,F2,G2,H1,H2}`; `results/{README.md,darkmatter_ecology/{ENDEMISM,GLOBAL_ENDEMISM,SELECTION_ADAPTIVE}.md,global_specificity/MULTICATALOG.md}`; `code/{phase_c_adaptive.py,phase_d_endemism.py,phase_d_global.py,ibex/96_darkfam_dnds_aggregate.sh}`. Read-only throughout; no web access, no subagents, no writes.

## Post-fix audit trail

The first focused follow-up correctly returned **FAIL** because the manual
147M correction in `MULTICATALOG.md` was overwritten as 148M by
`make_global_specificity_figure.py`. Codex accepted the finding, fixed the
producer, strengthened the multiline novel-fold-family regression guard,
clarified the result/figure regeneration boundary, reran the full harness,
and committed the synchronized private manuscript at `5222c55`.

A second read-only Claude Opus run then returned:

> **PASS. All six cleanup items are resolved, no regression was found, and
> the 147M value survives producer regeneration byte-for-byte.**

It reported no blockers and independently verified:

- `python3 reproducibility/verify.py --quiet` exits zero and emits only the
  two intended summary lines, with 62 COMPUTE, 72 FROZEN, six TABLE and 25
  PROVENANCE claims and zero failures;
- the bounded language regex handles hyphens, whitespace and newline-separated
  “novel fold family” forms, while the documentation accurately describes it
  as a curated regression guard rather than a semantic checker;
- the producer's exact desert-union target count remains 147,531,178 while its
  manuscript-consistent display label is 147M, and a live regeneration leaves
  the Markdown, PDF and PNG byte-identical;
- the README and regeneration script agree on derived result prose, 11 active
  manuscript figures and the separately asserted bounded RQ/E2 products;
- the abstract, Methods, frozen dRep reconciliation and verifier consistently
  distinguish 12,434 current HQ/MQ bins, 12,433 historical dRep members and
  5,438 representatives;
- `git diff --check` is clean, the full report ends with zero failures and a
  clean 20-page BibTeX-resolved manuscript, `manuscript/main.tex` is
  byte-identical to `paper/main.tex`, and the three regenerated private
  figures match their public copies.

The final audit used Claude Code `claude-opus-5` at high effort, made no
writes, used no web access or subagents, and exited successfully.
