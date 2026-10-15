# RQ-C4 — Result

**Question.** Can structure-based annotation rescue function that sequence/KO annotation misses?

## Answer
Structure-based search recovers a **structural** match for **10.76 %** of sequence-dark proteins
(those with no UniRef50 hit), but it recovers a *function hint* — a match to a characterised
(reviewed-SwissProt) target — for only **~1.0 %** of them. Structure rescues *neighbourhood*, not
*function*, for the great majority of the dark matter.

## Numbers
- Sequence-dark reps queried structurally (ProstT5 3Di → FoldSeek vs AFDB): **N_query = 2,839,030**.
- Gain an AFDB **structural** hit despite no UniRef50 hit: **N_afdb_hit = 305,584** → **structure-rescue
  rate = 305,584 / 2,839,030 = 10.76 %** (matches the expected ~11 %).
- Remain structurally novel (no AFDB hit): **N_nohit = 2,533,446 (89.24 %)** — `305,584 + 2,533,446
  = 2,839,030`, an exact partition.
- Of the 80,231,849 hit rows, only **74,631** are to a reviewed SwissProt target.
- Rescued proteins whose AFDB match is a **reviewed SwissProt** entry (curated function = a
  defensible "characterised" call): **28,254 of 305,584 = 9.25 % of the rescued**, i.e.
  **28,254 / 2,839,030 = 0.995 % ≈ 1.0 % of all sequence-dark proteins**.

So of every ~100 sequence-dark proteins, ~11 gain a structural neighbour in AFDB, but only ~1
land next to a protein anyone has actually characterised; the other ~10 hit AFDB models that are
themselves TrEMBL/uncharacterized — a fold neighbour, not a function.

## Caveat
An AFDB structural hit to **another hypothetical/uncharacterized** protein is **not** a function
assignment; it only says "this fold recurs in the predicted-structure universe." We therefore use
reviewed-SwissProt membership as a conservative, accession-level proxy for "characterised" — it is
a **lower bound** (some informative TrEMBL targets are excluded) and, conversely, even a SwissProt
hit at modest TM/lddt is a *hypothesis* that still needs structural inspection, not a final
annotation. This is consistent with the near-saturation picture of AFDB structure space
(yeo2025afesm): once most folds are catalogued, "has a structural neighbour" stops implying "has a
known function." FoldSeek's TMalign mode (vankempen2024foldseek) was used for the search; the dark
reps themselves were predicted with ESMFold (lin2023esmfold), and the existence of a vast reservoir
of recurrent, function-unassigned folds echoes durairaj2023folds.

## Files
- `structure_rescue_summary.tsv` — key/value summary (this directory).
- `methods.md`, `run.py` — commands and the processing script.
- IBEX provenance: `M/04_structure/foldseek_results/{all_query_ids,afdb_hit_ids,struct_nohit_ids}.txt`,
  `afdb_hits.tsv`, job `c4_characterised.sh` → `c4_characterised_summary.tsv`.
