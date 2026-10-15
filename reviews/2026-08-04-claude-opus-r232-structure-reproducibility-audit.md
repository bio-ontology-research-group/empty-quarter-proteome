# Independent second-reviewer audit — Empty Quarter MAG / dark-matter manuscript

Read-only audit of the working tree at `/home/leechuck/Documents/papers/empty-quarter-mags`
(branch `main`, uncommitted). No files were modified.

## Verdict: **FAIL** (3 blockers, 9 major, 10 minor)

The bounded-language discipline is genuinely good and most arithmetic checks out. It fails on
three things: an in-tree structural-novelty artifact that contradicts the B1a headline by ~2x in
the opposite direction from what newer references can produce; a reproducibility harness that
cannot run from a committed tree; and an order-gap conclusion that silently drops 1 of 11 candidates.

---

## Independently verified as correct

- All 9 SHA-256 digests in `reproducibility/current_manuscript_artifacts.tsv` match their files.
- `verify.py` runs clean: 96 COMPUTE pass, 0 fail, 34 DOC, exit 0. Matches `REPORT.md`.
- First-empty-rank census closes exactly: 11+142+1,343+3,755+187 = 5,438 representatives;
  LQ salvage 1+5+76+208+14 = 304. HQ splits (6/64/530/1,555/104) also close.
- Ladder partitions close: 2,533,446+305,584 = 2,839,030; 334,810+1,540,730 = 1,875,540.
  Percentages 12.1 / 89.2 / 10.8 / 92.26 / 7.74 / 77.27 / 15.3 all recompute correctly.
- 151,400 non-closure arithmetic correct (11,418,079+57,379,961-68,646,640).
- 608/1,875 = 32.427% -> 32.43%; 101/3,013 = 3.352% -> 3.35%; Wilson CI 2.77-4.06% correct.
- Specificity: 23,382-18,371 = 5,011 = 21.4%; 47.4 / 36.6 / 31.9 match the summary TSVs;
  21.4 <= min(31.9) is internally consistent.
- Endemism 10,669 (45.6%) / 543 (2.3%) and dN/dS 0.332 / 0.420 match `ENDEMISM.md`.
- RED adjudication closes: 6 singleton + 4 RED-incompatible = 10; 0 assignments, 0 proposals.
- A1/A2/A5 values recompute from `results/`. No undefined citations; all 128 bib entries carry
  DOIs; ~20 spot-checked DOIs are correct.

---

## BLOCKERS

### BL-1 — In-tree RQ/A3 contradicts the B1a headline; the TM normalisation change is undisclosed
`RQ/A3/fold_validation_summary.tsv:40` records `plddt_ge_0.60_length_ge100_dual_novel_n = 300`,
and `RQ/A3/result.md` states "Only 300 dual-novel representatives are both pLDDT >=0.60 and at
least 100 aa". `paper/main.tex:210-211` claims **608** dual no-hits at the strictly tighter gate
(pLDDT >= 0.70 **and** >= 100 aa). pLDDT>=0.70 is a nested subset of pLDDT>=0.60, and the B1a scan
adds *newer* references — both directions can only reduce the count. Unfiltered rates diverge the
same way: 15,248/63,654 = 24.0% (RQ/A3) vs 41,127/68,527 = 60.0% (`main.tex:210`).

Cause: `RQ/A3/methods.md:28-30` retains FoldSeek `alntmscore` (alignment-normalised);
`main.tex:209` uses "query-normalised FoldSeek/TMalign TM score". Query-normalisation penalises
partial coverage and inflates no-hits. Nothing in the manuscript, `reproducibility/README.md`, or
RQ/A3 marks the earlier result as superseded or explains the switch.

Why it is load-bearing: under the alignment-normalised rule the confident >=100 aa candidate rate
is 300/7,536 = **4.0%**, against the 3.35% known-fold error floor — the "calibrated candidate
folds" conclusion (`main.tex:216-220`) does not survive. Under query-normalisation it is 32.43%.
An 8x swing driven entirely by a normalisation choice, reported without sensitivity.

Required: state and justify the normalisation, report both, and mark RQ/A3 superseded in-tree.

### BL-2 — The reproducibility harness does not run from a committed tree
Everything the harness depends on is untracked or uncommitted:
- All 9 frozen artifacts (`reproducibility/*.json`, `r232_*.tsv`, `current_manuscript_artifacts.tsv`) — untracked.
- `verify.py:49-54` loads five of those JSONs at *import* time; a clean clone dies before the first check.
- `results/picrust2_bias/rebuild/` (drives `c1.rebuild_n`, `c1.richness_*`) — untracked.
- `results/atmospheric/aridity/` (drives `a1.aridity.*`) — untracked.
- `code/make_darkmatter_dnds_figure.py` — untracked, yet imported at
  `code/make_darkmatter_figures.py:15`, so `regenerate_figures.sh` fails on a clean clone.
- `paper/` submodule: `main.tex` is uncommitted, and the recorded submodule commit `45607b2`
  still asserts the superseded "dual-database-validated novel-fold floor 10,007 (16.1%)".

As deposited, the claim "fail-closed CI gate" (`reproducibility/README.md:22-23`) is unverifiable.

### BL-3 — 1 of 11 order gaps was never adjudicated; the "no new order" conclusion overreaches
`workbench/results/rebuild_20260729/taxonomic_novelty_r232/README.md`: "**Ten of those 11**
canonical genomes passed the marker/MSA staging gate ... One was retained in the audit but
**excluded from staging because its GTDB-Tk MSA percentage was below 50%**."

The manuscript reports only the 10 everywhere — `main.tex:45-46` (abstract), `120-128` (A0a),
`378` (Table 1: "staged / advanced 10 / 0"), `578` (Box 1) — and concludes at `main.tex:127-128`
"the existing catalogue supports no new phylum, class or order". Low marker-MSA recovery is
positively correlated with deep divergence, so the excluded genome is exactly the one a reviewer
will ask about. Disclose: 11 order gaps, 10 staged, **1 excluded (MSA < 50%), unadjudicated**, and
weaken the conclusion accordingly. (Zero domain/phylum/class gaps *is* supported by the primary
census for all 5,742 genomes — that part stands.)

---

## MAJOR

### MJ-1 — C1's "shared relative-abundance threshold" is not shared
`main.tex:279` and Table `main.tex:478`. The quoted column is
`ratio_picrust_over_mag_presence` in `results/picrust2_bias/rebuild/threshold_sweep.tsv`.
Generator `workbench/scripts/analysis/174_picrust2_rebuild.py`: line 167 thresholds PICRUSt2 at
`> t`; line 171 computes the MAG denominator at `v[qi] > 0` — **fixed at zero, never at `t`**.
`mean_mag_kos_presence` is constant 4346.48 across the whole sweep, which proves it.
The genuinely shared-threshold column `ratio_picrust_over_mag_relab` gives **31.4x at 0 and 14.6x
at 1e-5** — i.e. the honest shared-threshold comparison supports *more* over-prediction, not less.
The script's own docstring (line 13) claims the threshold is "applied identically to both sides".
The 2.4x withdrawal may still be right, but the stated reason is false as written.

### MJ-2 — A degenerate constant is reported as a per-sample median
`results/picrust2_bias/per_sample_concordance.tsv` has `n_pic_ko = 10543` for **all 436 rows**
(1 distinct value). `main.tex:277` ("median 10,543 predicted ... KOs per sample") and the Table row
`main.tex:477` present it as a median. The upstream script says so explicitly
(`174_picrust2_rebuild.py:6-9`: "a constant, not an observed richness"). `verify.py:380` takes a
median of a constant and passes. Label it as an artefact of the non-zero presence rule.

### MJ-3 — "~74%" is misdescribed as a survival fraction
`main.tex:196`: "only 2,339,154 no-matches (~74%) survive a >=100 aa gate".
Source `workbench/results/tier1_structural_novelty/SUMMARY.md:32,39`: 74% is the **no-AFDB-match
rate among >=100 aa representatives**, not the fraction surviving the gate. The surviving fraction
is 2,339,154/10,534,391 = **22.2%**. Table row `main.tex:444` repeats the ambiguity. Note the
correct reading also weakens the "fragment-inflated" framing: length-gating moves the rate
92.26% -> 74%, not to something small.

### MJ-4 — The known-fold control may be measuring retrieval failure, not TM discrimination
`known_fold_current_reference_summary.json`: `structure_defined_known_fold_controls` 3684,
`current_structure_union_recovered` 3244. **3684 - 3244 = 440**, which is exactly the ungated false
dual-no-hit count quoted at `main.tex:216`. If a missing search row is scored TM = 0 — as the
sibling analysis explicitly does (`RQ/A3/methods.md:41`) — then the ungated control error is 100%
retrieval failure and 0% TM-threshold error, and the control does not calibrate the TM < 0.5 rule
at all. `current_structural_gate_summary.json`'s `derivation` string never says how missing/empty
`current_best_*_qtm` values are handled. Required: state the missing-row rule, and break the 101
into "retrieved, TM < 0.5" vs "no row returned".

### MJ-5 — ~45 of the 96 "COMPUTE" passes are literal-vs-literal, not recomputations
`verify.py:71-116`: `c_catalogue`, `c_r232`, `c_r232_gap`, `c_r232_named`, `c_phylogeny`,
`c_structure`, `c_structure_gate`, `c_known_fold` each read one scalar out of the same frozen JSON
whose hash is pinned by `verify_artifact_manifest` (`verify.py:465`). Every `scale.*` and `b1a.*`
claim is of this form. `reproducibility/README.md:43-46` describes COMPUTE as "recomputed here from
a local data table ... reproduced from data artifacts". That is overstated — nothing upstream of the
frozen JSON is checked. Split the tier, or rename it.

### MJ-6 — `verify.py` fails **open** on the single most load-bearing claim
`verify.py:83-85`: if `domain`/`phylum`/`class` are absent from
`r232_first_unassigned_census.tsv`, the function returns `0.0`. `scale.gap_domain`,
`scale.gap_phylum`, `scale.gap_class` therefore pass **by construction** — a regenerated census that
dropped or renamed rows would still report zero gaps. Absence of a row is being read as evidence of
zero, which is precisely the conservative-treatment rule this paper otherwise follows. Make it raise
unless an explicit zero row exists.

### MJ-7 — `verify.py` does not verify the manuscript, and has already drifted from it
Every `PAPER` value is a Python literal in `CLAIMS` (`verify.py:527-676`).
`verify_manuscript_boundaries` (`verify.py:488-518`) only substring-checks 8 required and 6
forbidden phrases. 51 named phyla, 12.1%, 2,839,030, 56.7%, 45.6%, 21.4% are never read from
`main.tex`. Confirmed drift: `verify.py:585` records the paper value for terpene BGCs as **64.6**
while `main.tex:156` prints **64.5** (recomputed 9,984/15,469 = 64.54); `tol=0.3` hides it. Parse
the numbers out of `main.tex`.

### MJ-8 — The deposited source-of-truth README still carries the retracted language
`results/global_specificity/README.md` (Results section, ~line 48-52): "no detectable homolog
**anywhere sequenced on Earth**" and "the signature of **genuine endemism, not undersampling**".
These are verbatim the phrases `verify.py:495-497` forbids in `main.tex`. The guard covers only
`paper/main.tex`; a reviewer reading the repro package sees the overclaim the paper retracted.
Extend the guard to `results/**/*.md`.

### MJ-9 — dN/dS regime is not addressed
B2 (`main.tex:225-241`) runs pairwise NG86 between members of the same 50%-identity orthogroup,
which spans within- and between-population divergence. dN/dS is not a valid selection statistic in
the within-population regime; a median of 0.391 is expected there without invoking purifying
selection. The existing caveats cover alignment/saturation filters and the withdrawn comparator,
but not this. Add the caveat and the primary reference (see below).

---

## MINOR

| # | Finding | Location |
|---|---------|----------|
| m1 | 9 of 13 `paper/figures/*.pdf` are stale copies of `figures/*` (e.g. `fig_picrust2_bias.pdf` 2026-07-29 vs 2026-08-04); built `main.pdf` renders the old panels | `paper/figures/` |
| m2 | Pinned artifact internally incoherent: `matched_pairs_same_length_and_plddt_bins` 275 vs `matched_target_dual_novel` 3661 / `matched_target_novel_control_not` 3225; `matched_two_sided_exact_sign_p` = 0.0 exactly (not a valid exact p) | `reproducibility/current_structural_novelty_summary.json` |
| m3 | Docstring says TABLE tier asserts the source "contains the value"; the code only checks existence | `verify.py:11-12` vs `verify.py:732-734` |
| m4 | PROVENANCE is unconditionally `DOC` (`tier == P or exists`), so 51,920,778 / 23,416,070 / 2,839,030 / 753,305,533 / 68,808 get zero checking. Disclosed in README but should be stated in the paper as the reproducibility ceiling | `verify.py:734` |
| m5 | `npres <= 1` conflates absent-everywhere with single-sample; `nsamp` computed and never used, so the 551-column shape is never asserted | `verify.py:331,339` |
| m6 | Dead code pointing at the superseded R226 catalogue: `c_nphyla`, `c_actino`, `c_pseudo` unused | `verify.py:122-149` |
| m7 | `c_mgnify` hardcodes numerator 8552; `c_tier_snapshot_nonclosure` hardcodes all three operands and cites `verify.py` as its own source — it cannot fail | `verify.py:118-119,459-462` |
| m8 | `fig_darkmatter_endemism.pdf` is regenerated and asserted by the harness but never cited in `main.tex` | `regenerate_figures.sh` |
| m9 | 184 named species (`r232_named_taxon_census.tsv`) never reported next to the 187 genomes "assigned through species"; readers will take 187 as a species count | `main.tex:113,377` |
| m10 | Draft markers still live: `\todo` in the abstract, Fig 1 a placeholder `\fbox`, date "Draft outline --- August 2026" | `main.tex:31,57,518-524` |

---

## Missing primary references (DOI-verified, high confidence)

Tools and thresholds that are load-bearing in Methods but have no citation anywhere in
`references.bib`:

| Purpose | Reference | DOI |
|---|---|---|
| DAS_Tool (`main.tex:335`) | Sieber et al. 2018, Nat Microbiol 3:836-843 | 10.1038/s41564-018-0171-1 |
| GUNC (`main.tex:338`) | Orakov et al. 2021, Genome Biol 22:178 | 10.1186/s13059-021-02393-0 |
| **TM = 0.5 threshold** (the whole B1a headline) | Xu & Zhang 2010, Bioinformatics 26:889-895 | 10.1093/bioinformatics/btq066 |
| TM-align | Zhang & Skolnick 2005, NAR 33:2302-2309 | 10.1093/nar/gki524 |
| RCSB PDB (headline reference DB) | Berman et al. 2000, NAR 28:235-242 | 10.1093/nar/28.1.235 |
| MMseqs2 **linclust** (`main.tex:178`) | Steinegger & Soding 2018, Nat Commun 9:2542 | 10.1038/s41467-018-04964-5 |
| dN/dS validity regime (MJ-9) | Kryazhimskiy & Plotkin 2008, PLoS Genet 4:e1000304 | 10.1371/journal.pgen.1000304 |
| Wilson CI (`main.tex:215`) | Wilson 1927, JASA 22:209-212 | 10.1080/01621459.1927.10502953 |
| AFDB current release (paper uses **v6**, cites only the 2022 paper) | Varadi et al. 2024, NAR 52:D368-D375 | 10.1093/nar/gkad1011 |
| smORF context for the 77.8%-under-100-aa fraction | Sberro et al. 2019, Cell 178:1245-1259 | 10.1016/j.cell.2019.07.016 |

Already in `references.bib` but **uncited**, while Methods names the tool/standard:
`chaumeil2022gtdbtk`, `parks2022gtdb`, `parks2018gtdb` (RED), `bowers2017mimag`,
`chklovski2023checkm2`, `olm2017drep`, `pan2023semibin`, `schwengers2021bakta`,
`cantalapiedra2021eggnog`, `steinegger2017mmseqs2`, `suzek2015uniref`, `larralde2022pyrodigal`,
`eddy2011hmmer`, `buchfink2021diamond`, `jain2018fastani`, and the three contamination-control refs
Box 2 leans on (`salter2014contamination`, `eisenhofer2019lowbiomass`, `davis2018decontam`).

Framing gap: the "What is new" head-to-head (`main.tex:85-95`) does not engage the small-ORF /
gene-biogeography literature (`duan2024gmsc`, `coelho2022gmgc`, both in bib, both uncited), even
though 77.8% of tier1 no-hits are under 100 aa and the paper's response is to gate them away.

---

## Disposition — exact file:line

| Ref | Location | Disposition |
|---|---|---|
| BL-1 | `paper/main.tex:209-212`, `RQ/A3/result.md`, `RQ/A3/fold_validation_summary.tsv:40`, `RQ/A3/methods.md:28-30` | Disclose normalisation change; report both; mark RQ/A3 superseded |
| BL-2 | `reproducibility/*.json`, `results/picrust2_bias/rebuild/`, `results/atmospheric/aridity/`, `code/make_darkmatter_dnds_figure.py`, `paper` submodule | Commit all; bump submodule |
| BL-3 | `paper/main.tex:45-46,120-128,378,578` | Report 11 gaps, 10 staged, 1 excluded (MSA<50%) unadjudicated |
| MJ-1 | `paper/main.tex:279,478`; `174_picrust2_rebuild.py:167,171` | Rename to the actual rule or quote `ratio_picrust_over_mag_relab` |
| MJ-2 | `paper/main.tex:277,477` | Label 10,543 as a constant, not a median |
| MJ-3 | `paper/main.tex:196,444` | 22.2% survive the gate; 74% is the rate among >=100 aa |
| MJ-4 | `paper/main.tex:216`; `known_fold_current_reference_summary.json`; `current_structural_gate_summary.json` | State missing-row rule; split the 101 |
| MJ-5 | `verify.py:71-116`; `reproducibility/README.md:43-46` | Retier or rename; stop calling it recomputation |
| MJ-6 | `verify.py:83-85` | Raise instead of returning 0.0 |
| MJ-7 | `verify.py:527-676,488-518,585` vs `paper/main.tex:156` | Parse numbers from `main.tex`; fix 64.6 -> 64.5 |
| MJ-8 | `results/global_specificity/README.md` | Rewrite; extend guard to `results/**/*.md` |
| MJ-9 | `paper/main.tex:225-241` | Add within-population dN/dS caveat + citation |
| m1-m10 | see table above | as noted |
