# §2.1 — Reference and marker-gene methods cannot predict this gene content

Workbench section for paper section 2.1 (block **C1**). `M = /ibex/scratch/projects/c2014/rob/metagenomics`.

## STATUS (2026-08-16, session 5) — **DONE + SHIPPED**

Compute done (jobs 50586647, 50586699), §2.1 written (Results + Methods M.11), registered in
CATALOGUE_NUMBERS / CORRECTIONS / PLAN / PAPER_STRUCTURE, Discussion material deposited.
Shipped: **GitHub main `a5ac6d4`** · **Overleaf `8078aed`** · manuscript mirror synced. Paper compiles
clean (33 pp, 0 errors, 0 undefined citations).

**§2.1 IS NOT BLOCKED. The "BLOCKED on the PI's PICRUSt rebuild" label was inherited and is stale.**

The blocker as written in `paper/main.tex` and `notes/PLAN.md` says the section waits on a rebuild that is
"recomputing all CoverM profiles from raw reads before rebuilding the MAG–KO matrix". Both of those have
since landed, and neither is what actually stalled:

| Leg | State | Evidence |
|---|---|---|
| CoverM profiles on the current catalogue | **DONE 2026-08-06** | `rebuild_20260729/picrust2_current_20260803/coverm_current/profiles_full_recompute/` — 569 sample dirs, 5,438 reps each. §2.2/§2.3/§2.5 have shipped on it since 2026-08-11. |
| MAG–KO table | **DONE 2026-08-11** | `04_functional/repko_clean/results/rep_ko_v2.tsv` — built by us; 0 unannotated, `_sub` contig-subsetted, 8,134 spurious KOs removed. |
| 16S / PICRUSt2 predictions | **DONE, staged + checksummed** | `picrust2_current_20260803/inputs/ko_metagenome_unstrat.tsv`, 10,543 KOs × 1,270 profiles, `INPUTS_COMPLETE`. Plus `weighted_nsti.tsv`. |
| The PI's *comparison* (`mag_ko/`, `comparison/`) | **NEVER RAN** | Jobs `49985773/74/75` are `DependencyNeverSatisfied`: 4 of 573 CoverM array tasks failed in `49985116`. |

What stalled is the PI's **fail-closed 573/573 validation gate**, not any input we need. The 4 failures are
`51Sr1`, `F5Dr1`, `SPositive` (a control) and one timeout (`F13PRr1`). Excluding the control, exactly three
canonical samples lack a profile — which is precisely why the usable cohort is **548**, the basis §2.2/§2.3/§2.5
already ship on. We do not rerun them and we write nothing under `rebuild_20260729/`.

**The amplicon parquet cache is NOT a §2.1 blocker.** That blocker belongs only to §2.2's 16S *abundance*
cross-check. The PICRUSt2 chain reads real `.biom`/`.fasta` inputs that exist on IBEX under
`/ibex/project/c2014/EmptyQuarter_Data/soil/amplicon_16S/`; no `parquet` is referenced anywhere in it.

**Provenance closure:** the staged `ko_metagenome_unstrat.tsv` is sha256 `c71c2256…`, **byte-identical** to the
file the companion amplicon paper pins in its own `input_manifest.tsv`
(`repos/empty-quarter-ecology/analysis/v3/measured_function_summary_results/`). Both papers use the same
predictions; `sample_metadata.tsv` (`4f6b4f67…`) and `weighted_nsti.tsv` (`2a3468f5…`) match too.

## Cohorts

| Cohort | N | Definition |
|---|---|---|
| canonical soil | 551 | `metadata/sample_metadata.tsv` |
| MAG-profiled | 548 | canonical ∩ CoverM-profiled (minus `51Sr1`, `F13PRr1`, `F5Dr1`) |
| **§2.1 comparison** | **448** | MAG-profiled ∩ has a 16S/PICRUSt2 profile |

The 448 preserve the compartment balance of the 548 almost exactly — deep 30.4% vs 30.1%, plant-root 48.4% vs
48.7%, surface 21.2% vs 21.2% — so the 16S-matched subset is not compartment-biased. 4 of the 448 were
sequenced twice on the amplicon side; those columns are collapsed `sum_then_within_sample_normalize`, the rule
the companion paper documents in its own crosswalk.

## RESULT — the no-KO ceiling (headline; step 1 complete)

Recomputed per protein on the current 5,438-representative catalogue. All 5,438 resolved, 0 failures.

| Cohort | denominator | proteins | with a KO | **no KO** |
|---|---|---|---|---|
| **all 5,438** | **Bakta (primary)** | 19,186,853 | 8,516,190 | **55.61%** |
| all 5,438 | CheckM2 (historical convention) | 20,445,320 | 8,516,190 | **58.35%** |
| HQ (2,259) | Bakta | 8,541,971 | 3,934,561 | 53.94% |
| MQ (3,179) | Bakta | 10,644,882 | 4,581,629 | 56.96% |
| `_sub` (1,010) | Bakta | 3,427,370 | 1,504,509 | 56.10% |
| non-`_sub` (4,428) | Bakta | 15,759,483 | 7,011,681 | 55.51% |

**Two denominators, deliberately.** The KO numerator can only come from Bakta gene calls scored by eggNOG, so
the Bakta protein count is the internally consistent denominator and is what we quote. CheckM2's
`Total_Coding_Sequences` uses its own (pyrodigal) calls and runs ~6% higher per bin; it is reported only because
the historical **56.7%** was computed that way, and on that same convention the current catalogue gives
**58.35%** — i.e. the ceiling survives the basis change and is marginally higher.

**Not a completeness artifact:** HQ MAGs (≥90% complete, ≤5% contaminated) still carry 53.9% no-KO proteins.
The ceiling is a property of the reference, not of bin quality — mirroring §1.2's "length, not completeness,
discriminates darkness".

Cross-check: for `10Dr2__single_SemiBin_0` our KO numerator (1,754) reproduces the historical
`results/picrust2_bias/dark_fraction.tsv` row exactly; only the denominator differs, by the known caller gap.

## RESULT — the comparison (steps 2–3 complete, job 50586699)

MAG side: 11,395 distinct KOs across 5,438 reps; median MAG-encoded richness **8,666 KOs/sample**.
PICRUSt2 side: 10,543 KOs, predicted non-zero in **every** sample (a constant, not a median — as the
historical analysis also found).

| Metric (n = 448) | median | IQR |
|---|---|---|
| weighted NSTI | **0.159** | 0.115–0.209 |
| MAG-encoded KO richness | 8,666 | 8,169–9,003 |
| PICRUSt2 KO richness | 10,543 | constant |
| presence-recall of MAG KOs | **0.860** | 0.850–0.871 |
| Jaccard | 0.634 | 0.616–0.646 |
| Spearman ρ on shared KOs | **0.839** | 0.813–0.858 |
| fraction of a sample's MAG KOs outside the PICRUSt2 reference | **0.140** | 0.128–0.150 |

Catalogue-wide KO universe: 11,374 KOs detected in the cohort; 8,857 inside PICRUSt2's reference,
**2,517 (22.1%) outside it entirely**; 1,686 reference KOs never seen in any MAG.

**The three-tier ceiling — this is the section's argument.**
1. Placement is not the problem: NSTI 0.159, matching the historical 0.16 and the companion paper's 0.165.
2. On the shared core, agreement is genuinely good: ρ = 0.839 — we *independently reproduce* the companion
   paper's ρ ≈ 0.74–0.81 on a 3.6× larger matched cohort (448 vs 125 samples). We do not contest it.
3. The ceiling is everything outside that core: **55.6% of MAG-encoded proteins carry no KO at all**, and of
   the KOs that MAGs do encode, **22.1% lie outside PICRUSt2's reference universe** (14.0% per sample).

Headline the no-KO ceiling (55.6%). Never headline the over-prediction ratio — locked decision, and the
asymmetric-denominator problem that sank it is unchanged.

## Framing — the companion paper does NOT contradict this

`empty-quarter-ecology` (`main.tex`, "Predicted functional potential and shotgun comparison") reports PICRUSt2
vs shotgun agreement of **per-sample ρ = 0.743** and **community-mean ρ = 0.806** over **3,471 shared KOs** in
125 samples, concluding it "supports broad functional comparisons while leaving individual pathways as
predictions". Taken at face value that reads as the opposite of §2.1's title.

It is not, and §2.1 must say so explicitly rather than ignore it. That correlation is computed over the KOs the
two methods **share** — it measures accuracy *conditional on a gene family already being inside the KO
reference*. The ceiling measures the size of what is outside: a protein with no KO cannot enter that
intersection at any accuracy. Both are true and they compose into the section's actual claim:

> Placement is fine and the annotated core is predicted well; what cannot be predicted is the majority of the
> gene content, which the reference cannot represent at all.

This is also exactly the locked framing already in the frozen scaffold ("*not* a placement failure (median NSTI
0.16) … the remaining bias is a representational ceiling") and the locked decision to **headline the no-KO
ceiling, never the over-prediction ratio**. Citing the companion result as the setup strengthens §2.1: even at
ρ≈0.8 agreement on the shared core, the ceiling stands.

## Scripts

| Script | Does | State |
|---|---|---|
| `refbias_1_ceiling.py` | per-protein no-KO ceiling on 5,438; Bakta + CheckM2 denominators; `_sub` reps rebuilt with `repko_2_submap.py`'s contig-subset logic verbatim | **DONE** (job 50586647) |
| `refbias_2_magko.py` | MAG-encoded KO profiles per sample: `rep_ko_v2` × `profiles_full_recompute`, the `uv_2`/`osmo_2` join generalised to all KOs | running (50586699) |
| `refbias_3_compare.py` | PICRUSt2 vs MAG: recall, Jaccard, ρ on shared KOs, NSTI, and the KO-universe decomposition | running (50586699) |

Denominator note for `refbias_1`: `.emapper.annotations` contains rows **only** for queries with a seed-ortholog
hit — proteins with no hit are absent from the file, not written as `-`. Counting annotation rows as the
denominator would silently delete them and understate the ceiling, so the denominator always comes from `.faa`.

## Open / for the PI

1. His `178`/`180` chain is dead in the queue (`DependencyNeverSatisfied`); it needs the 4 CoverM tasks repaired
   and resubmission. We did not touch it.
2. When his `comparison/` eventually lands it **will disagree with ours**, for a known and documented reason:
   `178_current_picrust_mag_ko.sbatch` builds from his annotation backfill, which fixes the 168 unannotated reps
   but not the 1,010 `_sub` reps that inherit a parent superset — the defect `rep_ko_v2` corrects. Ours is the
   corrected basis. Reconciliation item, not a conflict.
3. His `--min-samples-with-abundance 573` gate is stricter than our 548 cohort by exactly the 3 failed canonical
   samples; ours matches the shipped §2.2/§2.3/§2.5 basis.
