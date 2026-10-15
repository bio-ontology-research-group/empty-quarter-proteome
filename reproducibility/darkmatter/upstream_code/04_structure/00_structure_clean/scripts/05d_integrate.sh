#!/bin/bash
# 05d_integrate — assemble the validation screen from the banked legs and produce §1.3's numbers.
# Inputs: sc5c_parts/pdb_hits.raw.tsv (banked from job 50603810) + afdb_hits.part1..4.tsv (05c_afdb).
# Output: per_rep_structures.tsv, per_control_structures.tsv, summary_05c.json, SUMMARY_05c.txt
# = K-of-Q confidence-qualified candidate folds, OUR OWN measured control floor (Wilson CI) from the
# 3,684 blinded known-fold controls co-searched under the identical invocation, and the
# alignment-normalised TM sensitivity companion. Logic lives once in 05_integrate.py.
# Submit after the 05c_afdb array: sbatch scripts/05d_integrate.sh
#SBATCH --job-name=sc5d_integrate
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/logs/sc5d_integrate_%j.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/logs/sc5d_integrate_%j.err
set -euo pipefail
M=/ibex/scratch/projects/c2014/rob/metagenomics
V=$M/04_structure/00_structure_clean/05_validate
CALDIR=$M/rebuild_20260729/structural_calibration
MANIFEST=$CALDIR/current_reference_control_manifest_20260802.tsv
PIREPS=$CALDIR/current_structural_novelty_20260730/per_candidate_current_structures.tsv
P=/ibex/user/velazqal/sc5c_parts
NPART=4

# --- gate: every leg banked ---
[ -s "$P/pdb_hits.raw.tsv" ] || { echo "ABORT: missing banked PDB leg"; exit 2; }
for i in $(seq 1 $NPART); do
  [ -s "$P/afdb_hits.part${i}.tsv" ] || { echo "ABORT: AFDB part $i not banked"; exit 2; }
done

# --- assemble; every query must appear in the AFDB union exactly once per part slice ---
mkdir -p "$V/screen"
cp "$P/pdb_hits.raw.tsv" "$V/screen/pdb_hits.raw.tsv"
cat "$P"/afdb_hits.part*.tsv > "$V/screen/afdb_hits.raw.tsv"
echo "pdb rows: $(wc -l < "$V/screen/pdb_hits.raw.tsv")  afdb rows: $(wc -l < "$V/screen/afdb_hits.raw.tsv")"
NQP=$(cut -f1 "$V/screen/afdb_hits.raw.tsv" | sort -u | wc -l)
echo "distinct AFDB-hit queries: $NQP (of 20,905 searched; queries with no hit are absent by design)"

python3 "$M/04_structure/00_structure_clean/scripts/05_integrate.py" "$V" "$PIREPS" "$MANIFEST"
echo "DONE_05d"
