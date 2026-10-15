# RQ-E2 — Methods

**Question.** What experimentally testable hypotheses emerge from the historical dark-family
pairwise-NG86 dN/dS-greater-than-one screen tail?

## Deposited inputs

All inputs are under `results/darkmatter_ecology/positive_selection/`:

- `positive_families.tsv`: 148 screen-tail families and their descriptive NG86, membership,
  taxonomy, low-confidence function, and environmental fields.
- `struct_master.tsv`: deposited ESMFold confidence, Foldseek analogue, P2Rank pocket, and nominal
  FUBAR-site fields.
- `struct_function_confident.tsv`: eight records passing Foldseek probability at least 0.9,
  TM-score at least 0.5, and model confidence at least 0.5.
- `boltz_affinity.tsv`: deposited Boltz-2 ligand-hypothesis scores for two exemplars.

The exact historical darkness and dN/dS provenance is checksum-enforced by
`reproducibility/dnds_legacy_provenance.tsv`. The historical sequence-extraction pipeline predates
the current corrected catalogue; its outputs are retained only for bounded sensitivity and
hypothesis generation.

## Local analysis

Run `python3 RQ/E2/run.py` from the repository root. The script uses the Python standard library
and matplotlib to:

1. summarise the 148-family screen tail and its 123-family recorded heuristic subset;
2. count structural-decodability tiers from `struct_master.tsv`;
3. attach a manual, description-based functional-hypothesis category to each of the eight
   confident analogues;
4. recompute predicted-pocket counts and the nominal-site tally directly from `struct_master.tsv`;
5. pass through the deposited Boltz-2 hypothesis table; and
6. regenerate the E2 tables and figure.

## Statistical boundary

Pairwise NG86 dN/dS greater than one is a screen, not a closed positive-selection test. Formal
claims require an appropriately specified codon model and controls for recombination, saturation,
phylogeny, and family ascertainment. Nominal FUBAR sites are not treated as confirmed selected
sites. The 23/332 in-pocket tally was not a preregistered enrichment test and has no inferential
p-value attached. Structural analogue labels and ligand scores generate hypotheses; they do not
validate molecular function or identify interaction interfaces.

FUBAR reference: Murrell et al. 2013 (`murrell2013fubar`).
