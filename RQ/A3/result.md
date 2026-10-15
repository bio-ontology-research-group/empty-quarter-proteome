# RQ-A3 — Result

> **Superseded historical sensitivity analysis.** This snapshot used
> alignment-normalised TM scores, older PDB/AFDB references, and a looser
> pLDDT ≥0.60 gate for its 300-model length-qualified count. The current
> manuscript result is the full-coordinate two-normalisation analysis frozen in
> `reproducibility/current_structure_metric_sensitivity_summary.json`: at
> length ≥100 aa and mean pLDDT ≥0.70, candidate no-hits are 608/1,875 by
> query-normalised TM and 97/1,875 by alignment-normalised TM; known-fold false
> no-hits are 101/3,013 and 68/3,013. This file remains only to reproduce the
> historical RQ/A3 figure and must not supply headline counts.

**Question.** Are the candidate novel folds genuine rather than
low-confidence ESMFold artefacts?

## Answer

**The dual-novel rate persists across pLDDT thresholds, but the analysis does
not establish a single validated novel-fold count.** The original 68,808
candidate representatives were selected using ProstT5-predicted 3Di and an
AFDB no-hit criterion. We re-tested explicit ESMFold coordinates against both
AlphaFold-DB and the curated PDB. Of 63,935 records in the pLDDT shards, 63,654
had an `ok` model and numeric mean pLDDT; 281 were marked `too_long` with no
pLDDT and were excluded.

A representative is **dual-novel** here when its best TM score is below 0.5
against both databases; no search hit counts as novel. Joining the 63,654
scored representatives to `best_afdb_tm.tsv` and `best_pdb_tm.tsv` gives:

| Minimum mean pLDDT | Representatives retained | Dual-novel | Rate |
|---:|---:|---:|---:|
| 0.00 | 63,654 | 15,248 | 24.0% |
| 0.40 | 50,154 | 9,258 | 18.5% |
| 0.50 | 33,524 | 5,514 | 16.4% |
| 0.60 | 17,196 | 2,836 | 16.5% |
| 0.70 | 5,709 | 1,164 | 20.4% |
| 0.80 | 187 | 15 | 8.0% |

Across thresholds that retain substantial numbers of models (unfiltered
through pLDDT ≥0.70), the novelty rate remains 16.4–24.0% rather than
collapsing as confidence rises. This is evidence against low-confidence
models alone manufacturing the novelty signal. The absolute count falls
because progressively fewer models remain. The pLDDT ≥0.80 threshold retains
only 187 of 63,654 models (0.3%), so its count is a selected-population result,
not a corrected novel-fold estimate.

Dual-novel representatives have lower median pLDDT than representatives with
a match (0.440 versus 0.520), while their median lengths are 182 versus 153 aa.

## Length sensitivity among confident models

Within models with mean pLDDT ≥0.60, novelty is concentrated among short
proteins:

| Length | Representatives | Dual-novel | Rate |
|---|---:|---:|---:|
| <100 aa | 9,660 | 2,536 | 26.3% |
| 100–150 aa | 1,839 | 99 | 5.4% |
| 150–250 aa | 3,421 | 91 | 2.7% |
| 250–400 aa | 1,831 | 88 | 4.8% |
| ≥400 aa | 445 | 22 | 4.9% |

Only 300 dual-novel representatives are both pLDDT ≥0.60 and at least 100 aa.
The short-protein concentration has two explanations that this analysis does
not separate: TM-align scores are unreliable for short structures, and
structural databases genuinely under-sample small proteins. A length- and
pLDDT-matched false-novelty control on known-fold proteins has not yet been
run. The candidates therefore cannot yet be described as validated novel
folds.

## Historical remaining gap (subsequently resolved)

At the time of this snapshot, per-family recurrence was not available because the novel-fold
representative identifiers and the family-level taxonomy table use different
clustering-derived identifier schemes. This analysis supports a candidate
structural-novelty signal that is not explained solely by low pLDDT, but it
does not determine how many genuinely new folds exist.

## Files

- `best_afdb_tm.tsv`, `best_pdb_tm.tsv` — best structural-analogue TM score
  per representative against AFDB and PDB.
- `fold_validation_summary.tsv` — the threshold, median, and length summaries
  reported above.
- `methods.md`, `run.py`, `finalize_tm.py` — provenance and processing code.
- IBEX source:
  `$M/06_esmfold_novelfold/{best_afdb_tm.tsv,best_pdb_tm.tsv,pdb/plddt_shard*.tsv}`.
