# RQ-B7 — Methods

**Question.** Dispersal limitation vs environmental filtering: how much of dark-gene
biogeographic structure is pure-space vs pure-environment vs shared?

## Data

- **Community matrix.** `results/darkmatter_ecology/darkfam_abundance.tsv` — 23,382 dark
  (UniRef50 no-hit, structurally novel) protein families × 551 soil samples (TPM-like
  read-recruited abundance). After aligning to samples with environmental metadata and
  applying a prevalence filter (family present in ≥5 samples), the working matrix is
  **547 samples × 21,923 families**.
- **Environment.** `results/darkmatter_ecology/env_matrix.tsv` — per-sample geography
  (lat/lon), climate (temp, precip, aridity index) and XRF geochemistry (V, Fe, Ca, K,
  Si, Ti, Mn, salinity proxy).

## Statistics (all self-contained in `code/eco_stats.py` — no scikit-bio/vegan)

1. **Community dissimilarity.** Hellinger transform of abundances → **Bray–Curtis**
   distance matrix.
2. **Distance-based PERMANOVA (adonis, McArdle–Anderson).** Per environmental block, with
   999 permutations: compartment, space (lat/lon), climate, geochemistry, vanadium. Each
   environmental block is run both *alone* (marginal R²) and *conditioned on space*
   (partial R², `block|space`) so that variance explained **beyond** geographic position is
   isolated. Produced by `code/phase_b_ecology.py` → `phase_b_summary.txt`.
3. **Distance-decay.** Mantel test between community Bray–Curtis distance and haversine
   geographic distance (999 permutations). A significant positive r is the classic
   signature of **dispersal limitation**.
4. **Borcard variance partition (this RQ, recomputed here on the real matrix).** A proper
   2-block partition (Borcard, Legendre & Drapeau 1992) of community variance into
   *pure space*, *pure environment*, *shared* and *residual* fractions, computed from
   PERMANOVA R²:
   - `ab  = R²(space alone)`, `bc = R²(env alone)`, `abc = R²(space+env combined)`
   - `pure_space = abc − bc`, `pure_env = abc − ab`,
     `shared = ab + bc − abc`, `residual = 1 − abc`.
   - Significance of each unique fraction = partial PERMANOVA of the block conditioned on
     the other (`block | other`, 499 permutations).
   Two environment blocks are partitioned against space: **CLIMATE** (temp, precip,
   aridity) on the climate-complete subset (n=456) and **GEOCHEM** (XRF elements) on the
   XRF subset (n=251).

## Reproduce

```bash
cd <repo-root>
python3 RQ/B7/run.py
```

Re-extracts the PERMANOVA blocks + Mantel r from `phase_b_summary.txt` into `b7_blocks.tsv`,
recomputes the Borcard space/env/shared partition into `b7_varpart.tsv`, and writes
`fig_varpart.pdf/.png`. Underlying pipeline: `code/phase_b_ecology.py`, `code/eco_stats.py`.

## Caveat on scope

The phase-B `varpart` block in `phase_b_summary.txt` only partitions *compartment vs space*
(climate was not complete for all 547 samples, so it was excluded from that 3-way call).
The clean **space vs environment** decomposition that answers RQ-B7 is the Borcard partition
recomputed here on the climate-complete and XRF subsets.
