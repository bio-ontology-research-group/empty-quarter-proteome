# RQ-A5 — Methods

**Question.** Do dark protein families encode genuinely novel enzymatic chemistry
(reactions absent from canonical annotation, or completing broken pathways)?

## Input
`results/darkmatter_ecology/family_ec.tsv` — 23,382 decoded dark families. Columns:
`family_id, best_ec, rxn_name, pathways, C_confidence, core_secondary, priority`.
These families are the FESNov-style novel families (no UniRef50 ortholog; carried by
known taxa) decoded by the C·V·T (Coverage · Validity · consensus) engine; `C_confidence`
is the decode confidence. Cf. Rodriguez del Rio et al. 2024 (`rodriguezdelrio2024fesnov`)
for the novel-family selection framing.

## Procedure (`run.py`, stdlib + matplotlib)
1. **Confidence cutoff.** Inspected the `C_confidence` distribution (min 0.045,
   median 0.276, max 0.806). Quartiles: Q1=0.193, Q2=0.276, Q3=0.369, P95=0.486.
   We adopt the **top-quartile cutoff C_confidence ≥ 0.369** as the "confident decode"
   threshold — a defensible, distribution-driven choice that isolates the upper tail
   without assuming a hard ground truth.
2. **EC parsing.** `best_ec` is comma-separated; each token may carry `/`-separated
   alternatives (e.g. `2.7.2.15/2.7.2.1`). A decode is counted as a *specific* enzyme
   call only when at least one fully-resolved 4-level EC (no `-` placeholder levels) is
   present. Confident + specific families are the analysis set.
3. **Landscape.** EC top-level class counts (each family counted once per class),
   full-EC counts, MetaCyc pathway counts (`pathways` parsed by splitting on commas and
   stripping `|`-delimiters), and the `core_secondary` split.
4. **Candidate novel-chemistry shortlist.** Families that are (i) confident, (ii) decoded
   to a fully-specified EC, (iii) carry a human-readable `rxn_name`, and (iv) map to at
   least one pathway (pathway-completion candidate). Ranked by `C_confidence`.
5. **Figure.** `fig_decode_landscape.pdf/png`: EC-class bar chart + top-pathway bar chart.

## What this evidence can and cannot establish
- **CAN:** the decoded-EC/pathway landscape of the dark fraction; enrichment of secondary
  vs core metabolism; pathway-mapped families that are candidate gap-fillers; a ranked
  shortlist of high-confidence, named-reaction candidates.
- **CANNOT (locally):** prove a reaction is "absent from every genome." We have no local
  per-MAG KO/EC coverage matrix here, so true pathway-hole completion and genome-wide
  absence cannot be confirmed; these remain *candidates*. Decode confidence is a model
  score, not experimental validation — biochemical confirmation is out of scope.

## Outputs
`ec_class_landscape.tsv`, `pathway_counts.tsv`, `ec_full_counts.tsv`,
`candidate_novel_chemistry.tsv`, `fig_decode_landscape.pdf/png`, `_summary.txt`.
