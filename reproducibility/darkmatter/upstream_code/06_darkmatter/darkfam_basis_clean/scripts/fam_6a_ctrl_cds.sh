#!/bin/bash
# Control-arm CDS slicing via the shared parameterized slicer (review fix F5).
#SBATCH --job-name=ale_ctrl_cds
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/darkfam_basis_clean/logs/fam6a_%j.out
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --partition=batch
set -euo pipefail
B=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/darkfam_basis_clean
R=$B/results
mkdir -p "$R/ctrl_cds_parts"
cut -f1 "$R/ctrl_member_fam.tsv" | awk -F'::' '{print $1}' | sort -u > "$R/_ctrl_samples.txt"
echo "[$(date -Iseconds)] control samples: $(wc -l < "$R/_ctrl_samples.txt")"
xargs -P "$SLURM_CPUS_PER_TASK" -I{} python3 "$B/scripts/slice_cds.py" "$R/ctrl_member_fam.tsv" "$R/ctrl_cds_parts" {} < "$R/_ctrl_samples.txt"
cat "$R"/ctrl_cds_parts/*.fna > "$R/ctrl_members.fna.part"
n=$(grep -c '^>' "$R/ctrl_members.fna.part"); exp=$(wc -l < "$R/ctrl_member_fam.tsv")
echo "[$(date -Iseconds)] control CDS: $n / $exp"
[ "$n" -ge $((exp * 99 / 100)) ] || { echo "ERROR: >1% missing"; exit 1; }
mv "$R/ctrl_members.fna.part" "$R/ctrl_members.fna"
echo "[$(date -Iseconds)] fam_6a done"
