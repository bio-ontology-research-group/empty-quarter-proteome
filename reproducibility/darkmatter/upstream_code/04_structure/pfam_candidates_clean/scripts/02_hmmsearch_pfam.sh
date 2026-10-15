#!/bin/bash
# Pfam-A 38.2 gathering-threshold scan of the three protein sets.
#
# One hmmsearch over the combined FASTA: the three sets are disjoint, and the
# cost of this search is dominated by walking the 30,134 Pfam profiles, so
# scanning them once for all 9,700 targets is ~3x cheaper than three runs.
#
# Plain bash worker. Runs on a login node or under scripts/02_hmmsearch_pfam.sbatch.
set -euo pipefail

W=/ibex/scratch/projects/c2014/rob/metagenomics
BASE=$W/04_structure/pfam_candidates_clean
PFAM=$W/rebuild_20260729/references/pfam_38_2/hmm/Pfam-A.hmm
# HMMER 3.3.2, the version the workbench Pfam 38.2 reference was built with
# (rebuild_20260729/references/pfam_38_2/hmmer_version.txt).
HMMSEARCH=/ibex/scratch/projects/c2014/rob/miniconda3/envs/bigslice/bin/hmmsearch

CPUS=${SLURM_CPUS_PER_TASK:-8}
QUERY=$BASE/sequences/all_sets.faa
OUTDIR=$BASE/results
mkdir -p "$OUTDIR"

test -s "$PFAM"
test -s "$QUERY"
test -s "$W/rebuild_20260729/references/pfam_38_2/PFAM_38_2_REFERENCE_COMPLETE"

"$HMMSEARCH" -h | sed -n '2p' > "$OUTDIR/hmmsearch_version.txt"
cp "$W/rebuild_20260729/references/pfam_38_2/Pfam.version" "$OUTDIR/pfam_version.txt"

echo "start $(date -Is)  cpus=$CPUS  targets=$(grep -c '^>' "$QUERY")"

"$HMMSEARCH" \
  --cut_ga \
  --cpu "$CPUS" \
  --noali \
  --notextw \
  --domtblout "$OUTDIR/pfam_domtbl.txt" \
  --tblout "$OUTDIR/pfam_tbl.txt" \
  -o "$OUTDIR/pfam_hmmsearch.out" \
  "$PFAM" "$QUERY"

echo "done  $(date -Is)"
grep -vc '^#' "$OUTDIR/pfam_domtbl.txt" | xargs echo "domtbl data rows:"
touch "$OUTDIR/PFAM_SCAN_COMPLETE"
