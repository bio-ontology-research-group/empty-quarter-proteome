# SUPERSEDED: Initial Hypothesis 2 (v1) Pilot Artifacts

These initial v1 artifacts (executed under Slurm Job ID `53701050`) are strictly preserved for auditability and provenance. They are **SUPERSEDED** by v2 due to methodological and technical issues identified in `reviews/darkmatter-2026-10-13/round6/ROOT_AUDIT_NOTES.md`:

1. **Critical Contig-ID Dropping:** `02_extract_neighbourhoods.py` keyed Bakta annotations by `(start, stop, strand)`, dropping contig identity within bins, causing coordinate collisions and misattributions across contigs.
2. **Global Pool Shuffling in Null:** `03_analyse_conservation.py` pooled all KOs globally across the entire dataset rather than preserving contig-specific or host-lineage annotation structures.
3. **Inadequate Species Callability Gate:** Blank species (unbinned contigs) were not strictly excluded from multi-species aggregation, and partial windows were counted toward multi-species callability without distinguishing full observable windows.
4. **Conflated Metrics:** The `exact_coordinate_matches` counter measured discovery of the focal protein in original assembly FASTA headers, not actual cross-caller Pyrodigal-to-Bakta matches.

All corrected logic, rigorous contig-keyed joins, host/contig-preserving null permutations, and audited tables are implemented in the v2 pipeline.
