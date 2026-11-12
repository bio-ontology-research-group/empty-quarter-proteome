# Exported sequence-dark repertoire: families, anchoring, structures and predicted complexes

This package holds the scripts and result tables behind the manuscript section
"The exported unknown families are folded proteins with predicted partners", the
family-sharing paragraph of the preceding section, Fig. 5 and the supplementary
table of supported gene pairs.

All analyses use the common gene calls of the harmonised export comparison
(`../external_export_20261110`): 3,617,011 proteins from 237 Rub' al-Khali and
1,720 external Thermoplasmatota genomes.

| Folder | Analysis | Software |
|---|---|---|
| `families/` | Protein families (30% identity, 80% coverage) and their sharing between desert and external genomes | MMseqs2 v18.8cc5c |
| `localisation/` | Export route, membrane topology, sorting signals and domain content of exported families | SignalP 6.0h, TMbed v1.0.2, HMMER v3.4 (Pfam-A, dbCAN HMMdb v14, TIGRFAM v15.0), DIAMOND v2.2.8 (MEROPS) |
| `structures/` | Mature-protein structure prediction, Foldseek searches, conserved adjacent partners, predicted heterodimers and control pairs | ColabFold v1.5.5 (AlphaFold2-ptm, AlphaFold2-Multimer v3), Foldseek |

Each folder contains `code/` (the scripts and Slurm files as run) and `results/`
(tables and a `receipt.json` with inputs, parameters, seeds and job identifiers).
The Slurm files contain the cluster paths used for the runs; set the paths at the
top of each file to rerun them elsewhere.

## Main tables

- `families/results/dark_exported_widespread_families_id30_c80.tsv`: the 981 exported
  sequence-dark families present in at least ten desert genomes from at least three
  genera or 95% ANI groups, with their presence in each external genome set.
- `families/results/sharing_fractions.tsv.gz`: sharing fractions by family class,
  prevalence filter, detection method and external genome set.
- `localisation/results/family_localisation_dark_exported_widespread_981.tsv`: route,
  anchor and role class of each family with the evidence columns; `rules.txt` states
  the classification rules.
- `localisation/results/summary_class_distributions.tsv`: class counts for the
  sequence-dark families and for the exported families with a UniRef50 match.
- `structures/results/per_family_monomer_metrics.tsv`: representative, alignment depth,
  model confidence and best PDB and AlphaFold Database matches for each family.
- `structures/results/complex_real_pair_metrics.tsv` and
  `complex_control_pair_metrics.tsv`: interface scores of adjacent and control pairs;
  `strict_support` marks the pairs that meet the support criterion of the manuscript.

## Figure and table

In the manuscript repository, `plot_exported_repertoire.py` renders Fig. 5 and
`make_complex_pair_table.py` writes the supplementary table from these tables. Set
`EQ_PROTEOME_REPRO_ROOT` to the `reproducibility` folder of this repository.

Predicted structures are large and are deposited with the bulk archive of the study.
