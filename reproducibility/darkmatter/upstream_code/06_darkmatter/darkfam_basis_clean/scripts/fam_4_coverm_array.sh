#!/bin/bash
# REBUILD step 4 — member-aware read recruitment (§1.4 leg 3: abundance/prevalence; also §2.5's matrix).
# Maps each sample's reads against the catalogue of ALL ~540k member CDS (not one rep per family),
# so a family is detected whichever of its diverged variants is present. CoverM contig / minimap2-sr,
# TPM + covered_fraction, same design as the June 91_darkfam_coverm_array. Aggregation to family
# level happens afterwards (sum TPM over members; present = any member covered_fraction >= 0.5).
# Reads manifest reused from the June run: darkfam_ecology/samples_darkfam.tsv (sample \t R1 \t R2).
#
# Submit (after fam_3): sbatch --array=1-551%30 scripts/fam_4_coverm_array.sh
#
#SBATCH --job-name=ale_fam_coverm
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/darkfam_basis_clean/logs/fam4_%A_%a.out
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --partition=batch

set -euo pipefail
B=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/darkfam_basis_clean
MAN=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/darkfam_ecology/samples_darkfam.tsv
CANON=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/darkfam_ecology/canon_samples.txt
CAT=$B/results/darkfam_members.fna
OUT=$B/results/coverm
mkdir -p "$OUT"
[ -s "$CAT" ] || { echo "ERROR: run fam_3_cds first"; exit 1; }

s=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$CANON")
[ -n "$s" ] || { echo "no sample for task $SLURM_ARRAY_TASK_ID"; exit 0; }
line=$(grep -P "^${s}\t" "$MAN" || true)
[ -n "$line" ] || { echo "[$s] not in manifest — skip"; exit 0; }
R1=$(echo "$line" | cut -f2); R2=$(echo "$line" | cut -f3)
F=$OUT/${s}.tsv
[ -s "$F" ] && { echo "[$s] done — skip"; exit 0; }
[ -s "$R1" ] || { echo "[$s] missing reads $R1"; exit 1; }

source /ibex/scratch/projects/c2014/rob/miniconda3/etc/profile.d/conda.sh
conda activate eco_probes
echo "[$(date -Iseconds)] [$s] coverm vs $(grep -c '^>' "$CAT") member CDS"
coverm contig -1 "$R1" -2 "$R2" -r "$CAT" -p minimap2-sr \
    -m tpm covered_fraction -t "$SLURM_CPUS_PER_TASK" -o "$F.part" \
    && mv "$F.part" "$F"
echo "[$(date -Iseconds)] [$s] done: $(wc -l < "$F") rows"
