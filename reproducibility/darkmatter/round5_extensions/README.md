# Round-5 exploratory extensions

Purpose: supply real additional analyses before the fifth Nature Communications
review, then support a bounded, falsifiable extension plan. These post hoc designs
were fixed in DESIGN.md before execution. RESULTS.md is the exact summary sent to
Argon; later root adjudications belong in the review response rather than silently
changing that submitted document.

## Inputs and execution

IBEX workbench: `/ibex/scratch/projects/c2014/rob/metagenomics`.
Scripts: `06_darkmatter/p4_review_20261013/scripts/round5/`.
Outputs: `06_darkmatter/p4_review_20261013/results/round5_extensions/`.
Job logs: `06_darkmatter/p4_review_20261013/logs/round5-*.log`.
Jobs 53697241 and 53697281 completed; 53697242 was cancelled before execution to
use exact original amino-acid lengths in the complete/long subset.

`profiles.sbatch` runs `profile_extensions.py` using metabolism-env (NumPy, pandas,
SciPy; environment manifests retained in ../environments). `culture.sbatch` runs
`culture_caller.py` using the same Python and Pyrodigal 3.7.1 from pyrodigal/bin.
Both are CPU Slurm jobs; no GPU, external API or database download is needed.
The scripts name all source paths explicitly. The CLI configuration is supported
by https://pyrodigal.readthedocs.io/en/stable/guide/cli.html .

Remote replay from the shared workbench (writes the same output directory):

    sbatch 06_darkmatter/p4_review_20261013/scripts/round5/profiles.sbatch
    sbatch 06_darkmatter/p4_review_20261013/scripts/round5/culture.sbatch

Input hashes: results/profile_input_sha256.json and
results/culture_caller/input_sha256.json. The latter also pins the newly produced
22 culture GFFs as evidence read by the locus-comparison stage. Those GFFs and
Pyrodigal stderr logs remain on IBEX. Raw assemblies, the full family presence
matrix, full functional profile matrices and PAFs stay at their documented IBEX
paths. Source snapshots and old results remain preserved. The original full-locus
Bakta classifications are exactly reproduced for all 68,537 eligible DNA-supported
calls before the alternative-caller summaries are accepted.

## Compact verification

From this repository:

    python3 reproducibility/darkmatter/round5_extensions/checks.py

This checks retained-output hashes, population keys, subgroup and pair-level
caller transitions, monotone threshold counts, per-family/null closure, actual
Jaccard/recall denominators, and paired-profile diagonals and permutation tails.
It does not rerun alignments, gene prediction or matrix permutations locally.
`outputs.sha256.tsv` freezes compact result files independently of manuscript
editing. `root_audit.json` adds genome-pair summaries and the number of defined
recall denominators for every threshold (only two at the highest cutoff).
The compact driver retains its narrower guarantees; whole raw-to-final
reproduction and historical claim-registry migration remain open.

## Interpretation boundaries

Caller substitution conditions on DNA alignments, compares overlap and first-base
phase, and does not validate translation. The recurrent-family culture subset is
small and has unchanged aggregate no-CDS counts. The temporal null conditions on
family prevalence/co-occurrence and site-by-position expedition counts, under label
exchangeability. It cannot separate temporal biology from trip-associated assay
effects. The profile null holds common KO identities fixed, includes zero ties and
reports each exchangeability cohort separately; thresholds do not define biological
truth. Every analysis is exploratory, including informative negative results.

After review, job53698522 ran the whole-visit sensitivity in VISIT_DESIGN_ADDENDUM.md
(2m1s, exit0). Code: scripts/temporal_visits.py and visits.sbatch. Results:
temporal_visits_summary.json, temporal_visits_draws.tsv and temporal_visits_design.tsv.
The null mean is80.54%, central95%range78.95–82.11%, lower-tail p0.11; the observed
79.38%is within this range. The original individual-label deficit is sensitive to
the permutation unit. Both results are retained. Read culture partial flags with
pandas dtype={"partial":str} to preserve00/01/10/11.
