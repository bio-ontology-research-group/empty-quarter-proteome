# Are dark-matter proteins "functionally novel"? — genomic-context test (negative for isolation)

Generated 2026-06-16 (`scripts/ibex/82g_darkmatter_context.sh`, job 47558327). Tests the hypothesis that
structurally-novel (reallydark) proteins resist function transfer not only by sequence (no UniRef50) and fold
(no AFDB) but also by **genomic context** — i.e. that they sit in uncharacterised genomic islands away from
characterized (KO-annotated) genes. If true, that would be a third independent line of evidence for *functional*
(not just structural) novelty.

## Method
For every reallydark protein (orthogroup members, expanded via `orthogroups_cluster.tsv`, across the 2,213
catalog MAGs that carry dark proteins; 100,970 dark proteins) classify its ±5 Bakta neighbours as
**characterized** (neighbour has a KEGG KO in EggNOG) vs not, and compare to a random non-dark gene null from
the same MAGs (874,787 genes). Bridge dark(prodigal)↔Bakta by C-terminal suffix.

## Result — the isolation hypothesis does NOT hold
| metric | dark proteins | random genes |
|---|---|---|
| mean characterized-neighbour fraction | **0.437** | 0.452 |
| proteins with **zero** characterized neighbour | 4.0% | 9.3% |
| uncharacterised island size around the gene | median 2, mean 2.4, max 66 genes | — |
| in island of ≥4 consecutive uncharacterised genes | 18% | — |

**Dark-matter proteins are embedded in characterized genomic context about as much as random genes** (43.7% vs
45.2%), and are *less* likely than random to be neighbourless-of-function (the higher null orphan rate is partly
a short-contig/edge-gene artifact). They are **not** clustered in uncharacterised islands.

## What this means (honest reframing)
- We **cannot** argue functional novelty from genomic isolation — dark proteins live among normal, characterized genes.
- This is in fact **consistent with** (and explains why) the existing context/FBA function-assignment engine works:
  ~34% of dark families *can* be given an EC-level function precisely because they have characterized anchors
  nearby (see `results/darkmatter_function_assignment/`).
- The defensible basis for "functionally novel" therefore rests on **(i)** no UniRef50 sequence homolog,
  **(ii)** no AFDB structural match, and **(iii)** the **~66% of dark families that still resist EC assignment**
  even with all evidence (sequence, structure, genomic context, FBA) — not on physical genomic isolation.
- Separately, dark proteins are **not** specifically enriched in the trace-gas operons (1.5%; `OPERON_DARKMATTER.md`).

## Bottom line
Two clean negatives that *sharpen* the paper by removing tempting overclaims:
1. The trace-gas machinery is characterizable (not built from novel folds).
2. The structural dark matter is **not** genomically isolated; it is interspersed with known genes.
The functional-novelty claim stands on the triple of sequence-novelty + fold-novelty + assignment-resistance,
quantified by the ~66% unassignable fraction — a stronger, data-backed statement than "genomic islands."

Files: `DARK_CONTEXT_summary.txt`, `dark_context_per_protein.tsv`.
