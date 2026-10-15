# RQ-E2 — Result

## Answer

This analysis is an exploratory follow-up to a historical pairwise-NG86 screen, not a formal
positive-selection test. The screen nominated 148 dark families with mean pairwise dN/dS greater
than one; 123 also pass the recorded minimum-dS and pair-count heuristic. Eight of the 148 have a
confident structural analogue under the deposited Foldseek and model-confidence gates. Manual
reading of those analogue descriptions places five in a broad recognition, binding, surface, or
transport hypothesis group. That small, selected subset generates experimental candidates; it
does not estimate category enrichment in the dark proteome.

## Evidence

- The 148-family tail is shallow (median three members), reaches dN/dS 3.65, and spans 14 phyla.
  Pairwise NG86 is sensitive to sampling, recombination, saturation, and small denominators, so
  dN/dS greater than one is used only as a screening rule. Codon-model confirmation is pending.
- Structural triage classifies 100/148 as having no detected PDB, Swiss-Prot, or AFDB analogue in
  the historical search; 40 have a lower-confidence analogue and eight pass all confidence gates.
  These are database- and threshold-relative categories.
- Five of the eight confident analogue descriptions motivate recognition/binding/surface/transport
  hypotheses; three are uncharacterised or otherwise outside that group. This is a descriptive
  manual categorisation of eight preselected families, not a statistical enrichment test.
- In families with a predicted pocket, 23/332 nominal FUBAR-positive sites fall inside predicted
  pockets. This recomputed tally is descriptive and was not a predeclared interface-enrichment
  test. It cannot localise selected substitutions to binding interfaces.
- Boltz-2 gives a modest GlcNAc-binding hypothesis for the lectin-like candidate
  `58PRr3_401505_1`; the scores support experimental prioritisation, not a validated function.

## Interpretation and limits

The defensible result is a short list of structurally interpretable candidates for codon-model and
experimental follow-up. It does not establish positive selection, interface diversification, or
an arms race. The shallow families that have interpretable analogues are also poorly powered for
site-level inference, while the deeper families tend not to have confident structural analogues.

## Reproducible products

- `table_positive_families_summary.tsv` — screen-tail size and descriptive properties.
- `table_structural_tiers.tsv` — database-relative structural tiers and the 23/332 pocket tally.
- `table_functional_categories.tsv` — eight confident analogues and manual hypothesis categories.
- `table_boltz_hypotheses.tsv` — deposited Boltz-2 hypothesis scores.
- `fig_posel_function.{pdf,png}` — descriptive structural-tier/category figure.

Regenerate all of them from deposited inputs with `python3 RQ/E2/run.py` from the repository root.

## Reference

Murrell B et al. 2013. FUBAR: a fast, unconstrained Bayesian approximation for inferring
selection. *Molecular Biology and Evolution* 30:1196–1205 (`murrell2013fubar`).
