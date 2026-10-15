# RQ-B8 — Result: temporal stability of the dark proteome

## Answer
**Both — a stable persistent core with significant quantitative turnover.** Dark-protein
family *membership* is largely stable across the five expeditions, but the *relative-
abundance composition* of the dark proteome turns over significantly between trips, and
trip (time) structures the community at least as strongly as site (space). The Empty
Quarter dark proteome is a large persistent reservoir whose internal abundance ranking
re-shuffles between expeditions.

## Sampling design
59 sites; 57 sampled in ≥2 trips, 7 in all 5 — so the temporal vs spatial contrast is
real, not confounded by one-off sites. Samples/sites per trip: trip1 206/58, trip2 37/7,
trip3 128/47, trip4 85/51, trip5 91/52 (547 samples, 23,382 dark families).

## 1. Persistence — a large temporally persistent core
Of 23,374 dark families detected in ≥1 sample:

| # trips detected | # families | % of detected |
|---|---|---|
| 1 (trip-restricted) | 606 | 2.6% |
| 2 | 1,440 | 6.2% |
| 3 | 1,929 | 8.3% |
| 4 | 2,971 | 12.7% |
| **5 (persistent core)** | **16,428** | **70.3%** |

**70.3%** of detected dark families are present in **all five** expeditions; only **2.6%**
are confined to a single trip. At the level of *which* novel families exist, the dark
proteome is highly stable over time.

## 2. Turnover — composition re-shuffles, and time beats space
Mean Bray–Curtis dissimilarity (depth-normalised dark-family vectors):

| pair class | n pairs | mean BC |
|---|---|---|
| within-site / within-trip (replicates) | 755 | **0.374** |
| across-site / within-trip (spatial) | 36,819 | **0.438** |
| within-site / across-trip (temporal) | 2,380 | **0.500** |
| across-site / across-trip | 109,377 | **0.514** |

Key reading: a sample resembles a **different site in the same trip (0.438)** more than it
resembles the **same site in a different trip (0.500)**. Returning to a fixed location a
year later yields a *more* dissimilar dark proteome than sampling a neighbouring location
on the same expedition. Repeat-site pairs are only marginally more similar than fully
different pairs (0.500 vs 0.514; one-sided MWU p = 6.1e-5 — significant but a small effect),
i.e. site identity barely persists through time.

## 3. PERMANOVA — trip explains more variance than site
(999 permutations, Bray–Curtis)

| factor | groups | pseudo-F | R² | p |
|---|---|---|---|---|
| **trip** (time) | 5 | **27.3** | **0.168** | 0.001 |
| **site** (space) | 59 | 1.25 | 0.130 | 0.008 |

Trip explains 16.8% of compositional variance with a large pseudo-F (27.3) over just 4 df;
site explains 13.0% but spread over 58 df with pseudo-F ≈ 1.25. Per degree of freedom,
**time is by far the dominant structuring factor** for the dark proteome's abundance
composition. Both factors are statistically significant.

## Synthesis
The dark proteome is **compositionally dynamic on a stable backbone**: ~70% of novel
families persist across all expeditions (true core, not undersampling artefacts), yet their
relative abundances turn over strongly enough between trips that the trip a sample came from
predicts its dark-proteome composition better than where it was collected. This matches the
boom–bust, pulse-driven dynamics typical of dryland microbial communities, where rare/dormant
taxa bloom episodically after rain and depth fluctuates between visits (e.g. Angel & Conrad
2013; Aanderud et al. 2015; Bay et al. 2021 on dormant arid-soil reservoirs).

## Caveats (important)
- **Trip = sequencing batch.** Each expedition is also a library/sequencing batch and
  CoverM abundance is depth-sensitive, so the trip R² (0.168) is an **upper bound** on
  biological temporal turnover; some of it is technical. The *persistence* result (presence
  across trips) is robust to this; the *turnover magnitude* is not fully decomposable here.
- **Unbalanced design:** trip2 has only 37 samples / 7 sites.
- Detection = nonzero coverage, so the persistent core is enriched for abundant families.
- Compartment (Deep/Root/Surface) variance is left in the residual, not modelled.

## Files
- `RQ/B8/run.py` — reproducible analysis (`python3 RQ/B8/run.py` from repo root)
- `RQ/B8/fig_temporal.pdf` / `.png` — turnover boxplot + persistence histogram
- `RQ/B8/turnover_by_class.tsv`, `persistence_histogram.tsv`, `permanova.tsv`, `summary.json`
