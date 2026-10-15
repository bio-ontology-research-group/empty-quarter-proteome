# Structure follow-up of the historical 148-family dN/dS screen tail

ESMFold models are present for all 148 historical screen-tail families. Under the deposited
Foldseek/model-confidence fields, 100 have no detected PDB, Swiss-Prot, or AFDB structural
analogue, 40 have an analogue below the combined confidence gate, and eight pass Foldseek
probability at least 0.9, TM-score at least 0.5, and model confidence at least 0.5. These labels are
database-, version-, and threshold-relative.

## Confident structural analogues

| family | members | dN/dS | pLDDT | pocket | nominal sites | TM | analogue-based hypothesis |
|---|---:|---:|---:|---:|---:|---:|---|
| V5PRr1_235427_14 | 4 | 1.026 | 0.51 | 0.392 | 23 | 0.92 | uncharacterised conserved protein |
| F2Sr2_109584_65 | 4 | 1.472 | 0.59 | 0.048 | 1 | 0.86 | toxin-like fold |
| F38PRr1_858279_42 | 3 | 1.392 | 0.71 | 0.156 | 0 | 0.83 | uncharacterised conserved protein |
| 36PRr1_1482479_11 | 3 | 1.162 | 0.64 | 0.478 | 0 | 0.78 | uncharacterised protein |
| 58PRr3_401505_1 | 4 | 1.260 | 0.79 | 0.972 | 34 | 0.74 | lectin-like binding hypothesis |
| V12Dr3_1574070_5 | 3 | 1.007 | 0.72 | 0.655 | 0 | 0.65 | DNA-binding regulator-like hypothesis |
| S11PRr1_690820_66 | 3 | 1.113 | 0.75 | 0.217 | 0 | 0.58 | ABC-transporter-like hypothesis |
| F4PRr2_783354_7 | 16 | 1.046 | 0.69 | 0.978 | 176 | 0.55 | lipid/cell-surface-carrier-like hypothesis |

Five analogue descriptions fall in a broad recognition/binding/surface/transport hypothesis
group and three are uncharacterised/other. This manual categorisation of eight preselected
families is descriptive and is not a category-enrichment test.

## Predicted-pocket overlap

Recomputation from `struct_master.tsv` finds 332 nominal FUBAR-positive sites among families with
P2Rank pocket probability at least 0.5; 23/332 fall inside a predicted pocket. This tally was not a
predeclared interface-enrichment test, has no valid inferential p-value attached, and cannot
localise selected substitutions to binding interfaces. The analysis therefore does not establish
positive selection, interface diversification, or an arms race.

Model confidence is modest, Foldseek matches can be driven by short regions, and most confident
analogues occur in shallow families. Treat every function and ligand label as an experimental
hypothesis. Regenerate the summary tables and figure with `python3 RQ/E2/run.py`.
