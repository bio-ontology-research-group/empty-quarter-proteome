# RQ-B8 — Methods: temporal stability of the dark proteome across the 5 expeditions

## Question
Is the dark proteome temporally stable across the 5 Empty Quarter expeditions ("trips"),
or does it turn over between trips?

## Data
- `results/darkmatter_ecology/env_matrix.tsv` — 547 soil samples with `site`,
  `compartment` (Deep / PlantRoot / Surface), and `trip` / `trip_num` (expedition 1–5).
  Trip is also encoded in the sample-name prefix: trip-1 samples are digit-led
  (e.g. `10Dr2`), trip-2 = `T*`, trip-3 = `F*`, trip-4 = `S*`, trip-5 = `V*`.
- `results/darkmatter_ecology/darkfam_abundance.tsv` — 23,382 dark-protein family ×
  551 per-sample CoverM abundance columns (~92 MB; **streamed** line-by-line in
  `run.py`, never loaded with an editor). 4 abundance columns (sites 62/64) have no
  env-matrix row and were dropped → 547 samples analysed.

## Sampling design (the design that makes the test possible)
The expeditions revisit the same field sites over time: of 59 sites, 57 were sampled
in ≥2 trips, 7 in all 5. This permits a direct **within-site-across-trip** (temporal)
vs **across-site-within-trip** (spatial) contrast rather than a confounded one.

Samples & sites per trip: trip1 206/58, trip2 37/7, trip3 128/47, trip4 85/51,
trip5 91/52. Trip2 is small and spatially restricted (7 sites) — see caveats.

## Analyses (`RQ/B8/run.py`, stdlib + numpy + scipy + matplotlib)

1. **Sample → trip map** built from `env_matrix.tsv`; counts of samples/sites per trip.

2. **Temporal persistence.** A family is "detected in trip *t*" if its CoverM abundance
   is > 0 in any sample of that trip. We count the number of distinct trips (0–5) each
   family is detected in → histogram, persistent core (all 5 trips) vs trip-restricted
   (1 trip).

3. **Community turnover.** Each sample's family-abundance vector is depth-normalised to
   sum 1, and **Bray–Curtis** dissimilarity is computed for all sample pairs
   (`scipy.spatial.distance.pdist`). Pairs are classified into four strata:
   - within-site / within-trip  (field/sequencing replicates; reference floor)
   - within-site / across-trip  (**temporal** turnover at a fixed location)
   - across-site / within-trip  (**spatial** turnover at a fixed time)
   - across-site / across-trip  (both differ)
   The key test — *are repeat sites more similar over time than different sites are?* —
   compares within-site-across-trip vs across-site-across-trip with a one-sided
   Mann–Whitney U test.

4. **PERMANOVA** (one-way, Anderson 2001; implemented from the Bray–Curtis distance
   matrix, 999 free permutations, seed 42) for `trip` and, separately, `site` as the
   grouping factor, reporting pseudo-F, R², and permutational p. This quantifies how
   much compositional variance trip (time) explains relative to site (space).

## Outputs
`turnover_by_class.tsv`, `persistence_histogram.tsv`, `permanova.tsv`, `summary.json`,
and `fig_temporal.pdf/.png` (turnover boxplot + persistence histogram).

## Caveats
- **Batch / depth confound.** Each trip is also a sequencing batch; CoverM abundances
  are depth-sensitive. A "trip" effect therefore partly reflects library prep / depth,
  not only biology. The PERMANOVA trip R² is an **upper bound** on true temporal turnover.
- **Unbalanced design.** Trip2 (n=37, 7 sites) is small; spatial and temporal strata
  have very different pair counts.
- Detection is presence-of-coverage; the persistent core is enriched for abundant
  families that are detectable at any sequencing depth.
- Compartment (Deep/Root/Surface) is not separated out here; it is an additional
  variance source folded into the residual.
