# Round-4 corrections

Gemini 4 Argon High identified useful issues alongside multiple incorrect quotations
and numerical claims. Root verified the issues before changing the paper.

IBEX job 53695688 completed in 40 seconds. `scripts/rerun.py` preserves original
upstream scripts and outputs, recomputes the soil-position contrast only on the
547 samples present in the environmental table, propagates its family memberships
into phylum and expedition summaries, and reruns the 224 osmolyte/repair partial
correlations using residual degrees of freedom n - rank(C) - 1 (n - 4 with two
independent covariates plus an intercept). Inputs and source scripts consumed by
that job are identified in `results/input_sha256.json`, including all 5,438 rRNA
GFF hashes. Job log: `reviews/darkmatter-2026-10-13/round4/ibex_corrections.log`.

The corrected root contrast tests 101,897 families: 20,136 enriched and 41,046
depleted. Excluded samples are 62PRr2, 64Dr3, 64PRr2 and 64Sr2. The reviewer supplied
incorrect sample names and a minimum-presence threshold; the implementation uses
25 nonzero observations, consistent with the original method. Updated phylum and
expedition summaries are in `results/compartment_families_summary.tsv`.

Corrected degrees of freedom leave every reported correlation coefficient
unchanged; 57 of 224 tests remain below 5% FDR (previously 58). The sulfur versus
photoprotective-pigment occurrence association changes from q=0.04981 to 0.05117.
The original gradient table is retained under `inputs/`; the generated corrected
script and table are under `results/`. Its inherited final print line names the
old output path, but the actual open/write target is the versioned correction
folder, and the original table was not overwritten.

The rRNA audit compares longest-hit selection against existence of a qualifying
complete hit in all 5,438 GFFs: zero call changes. The corrected parser in
`scripts/mimag_summarise_corrected.py` gives complete hits precedence; its regression
fixture checks that a longer partial hit cannot erase a complete 16S hit. Output
paths point to a separate correction directory if this full summariser is run.
The 126 MIMAG high-quality count is unchanged.

`permanova_blocks_corrected.tsv` changes only the two incorrect logged sample
counts (geochemistry core model 251; climate core model 456). All estimates and
p-values are preserved. `scripts/eco_3_community_corrected.py` contains the source
fix and writes to a separate correction folder; the expensive PERMANOVA was not
rerun because this was a sample-count logging defect. It retains the original
IBEX dependencies. The original source/table remain frozen.

The metabolic summary is frozen under `inputs/`. The matched Halobacteriota
comparison uses 181 pairs, separately from 237-genome lineage prevalences. Module
admission is two thirds of steps, not the review's incorrect 0.75 claim. CRISPR
markers use any of seven KOs, not two of four. Full definitions remain in the
checksum-pinned upstream source.

`../round4_checks.py` is a read-only compact verification step (40 checks).
`../run_current.sh` includes it. It verifies retained corrected outputs; it does
not claim to rerun the remote family abundance or GFF scans locally.
