# Dark-matter pipeline and evidence boundary

This document distinguishes heavy historical/upstream computation from the lightweight deposited
reproduction package. Cluster steps run through Slurm on IBEX or ORIX. Do not run read mapping,
large database searches, folding, or catalogue clustering on a login node or workstation.

## Analysis object

The ecology/dN/dS/specificity object is a historical set of 23,382 recurrent orthogroups built
from pyrodigal calls over 5,438 dereplicated MAGs. A representative was retained under the
historical darkness rule when it had no UniRef50 match or less than 30% identity, and each retained
orthogroup had at least three members. This object is distinct from both the 2,839,030
sequence-dark representatives in the Bakta ladder and the current 14,426-bin/R232 catalogue.

## Upstream abundance pipeline (IBEX)

1. `code/ibex/90_darkfam_catalog.sh` constructs the dark-family nucleotide catalogue and sample
   manifest.
2. `code/ibex/91_darkfam_coverm_array.sh` maps paired reads with CoverM/minimap2-sr; the documented
   recovery path handles stalled low-complexity samples.
3. `code/ibex/92_darkfam_aggregate.sh` produces the 23,382 × 551 abundance and presence matrices.

The large abundance input deposited for local reproduction is checksum-pinned in
`reproducibility/large_files.sha256` and fetched by `reproducibility/fetch_large_files.sh`.

## Local ecological summaries

- `code/build_env_matrix.py` builds the environmental matrix using the all-trip XRF table.
- `code/eco_stats.py` contains the deposited PERMANOVA/Mantel/variation-partitioning code and
  self-tests.
- `code/phase_b_ecology.py`, `code/phase_d_endemism.py`, and
  `code/phase_d_global.py` regenerate bounded historical ecological/specificity summaries.

Specificity means failure to detect a homolog in a named catalogue under a stated search rule.
For the joint comparison, 5,011/23,382 (21.4%) have no detected homolog in the named MGnify90,
GEM, and five-desert comparison. This does not establish universal absence or geographic
restriction.

## Historical dN/dS pipeline

The original dN/dS workflow predates a member-extraction guard correction. Its exact darkness
definition, withdrawn pre-guard scripts, retained result table, and family metadata are pinned in
`reproducibility/dnds_legacy_provenance.tsv`; the verifier checks every SHA-256, including immutable
Git objects for withdrawn code.

The manuscript therefore uses only bounded dark-only sensitivity:

- 22,767 recurrent orthogroups yield at least one filtered pairwise NG86 value;
- the descriptive median is 0.391 and increases with the minimum family-size floor; and
- the result supports coding-like constraint conditional on recurrent alignable ORFs, not a
  catalogue-wide estimate of selection strength.

The 148-family dN/dS-greater-than-one tail is an exploratory screen. The 23/332 predicted-pocket
overlap is descriptive and does not establish positive selection or interface localisation. Its
reproducible follow-up lives in `RQ/E2/`.

## Current structural calibration

The active structural claim does not come from the historical predicted-3Di no-hit count. Full
ESMFold coordinates were searched against PDB through 30 July 2026 and AFDB v6/UniProt 2025_03.
The deposited analysis reports both query-normalised and alignment-normalised TM, decomposes
missing rows from retrieved-below-threshold rows, and applies identical length/confidence gates to
known-fold controls. Compact results are SHA-256-pinned in
`reproducibility/current_manuscript_artifacts.tsv` and recomputed into manuscript claims by
`reproducibility/verify.py`.

## Reproduce the deposited layer

```bash
bash reproducibility/fetch_large_files.sh
bash reproducibility/run_all.sh
```

The harness checks provenance and table arithmetic, regenerates all figures, and builds the public
manuscript snapshot. It does not claim to rerun the 51.9-million- or 753-million-protein upstream
pipelines locally.
