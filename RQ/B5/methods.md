# RQ-B5 — Methods

**Question.** How is the dark proteome partitioned by depth/compartment
(Surface / Deep / PlantRoot-rhizosphere)?

## Inputs
- `results/darkmatter_ecology/env_matrix.tsv` — `sample`, `compartment`
  (labels: `Surface`, `Deep`, `PlantRoot`), plus site/space/climate/geochem.
- `results/darkmatter_ecology/darkfam_abundance.tsv` — 23,382 dark-protein
  families × 551 per-sample CoverM relative-abundance columns (~92 MB;
  **streamed line-by-line**, never loaded whole).
- `results/darkmatter_ecology/darkfam_meta.tsv` — family → MAG/locus/size.
- Context: `results/darkmatter_ecology/phase_b_summary.txt` (PERMANOVA,
  compartment R² = 0.020, p = 0.001) and `results/uv_survival/compartment.tsv`
  (per-compartment survival-marker prevalence).

## Procedure (`run.py`, stdlib + scipy + matplotlib)
1. **Sample → compartment.** Parse `env_matrix.tsv`; keep the three compartment
   labels. 547 samples carry a compartment (Surface 116, Deep 165,
   PlantRoot 266). 547/551 abundance columns map to a compartment; the 4
   unmatched columns are dropped.
2. **Per-family statistics.** Streaming the abundance matrix, for each family
   collect its abundance vector split by compartment, and record per-compartment
   **mean abundance** and **presence prevalence** (fraction of samples with
   abundance > 0). Test compartment specificity with a **Kruskal-Wallis** test
   across the three compartment groups (non-parametric; CoverM abundances are
   zero-inflated and non-normal). Families with < 2 non-empty/non-zero groups are
   left untested (NaN).
3. **Multiple testing.** Benjamini-Hochberg FDR across the 22,430 families with a
   finite p-value; monotonicity enforced. Significant = **q < 0.05**.
4. **Compartment assignment.** Each significant family is labelled by the
   compartment with the **highest mean abundance**
   (surface- / deep- / rhizosphere-enriched).
5. **Outputs.** `compartment_partition.tsv` (per-family means, prevalences,
   H, p, q, enriched_in), `compartment_summary.tsv` (headline counts), and
   `fig_compartment.{pdf,png}` (panel A: enriched-family counts per compartment;
   panel B: the strongest example family per compartment).

## Caveats
- **Unequal power.** Surface has the fewest samples (116 vs Deep 165,
  PlantRoot 266), so its lower enriched-family count is partly reduced statistical
  power, not only biology. Counts should be read as a floor for the
  smaller compartments.
- KW tests *distributional* differences; "enriched_in" uses the mean and can
  flip for families with similar means across compartments.
- Abundance is CoverM read-recruitment, i.e. community-level relative abundance,
  not expression; this is partitioning of *who carries* the dark families, not
  in-situ activity.
- Compartment explains little community variance overall (PERMANOVA R² = 0.020),
  so per-family enrichment, while real, sits on top of a weakly
  compartment-structured assemblage dominated by space/climate.

## Reproduce
```bash
python3 RQ/B5/run.py   # from repo root
```
