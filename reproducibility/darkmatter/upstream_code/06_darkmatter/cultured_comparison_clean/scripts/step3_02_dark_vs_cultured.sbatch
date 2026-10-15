#!/bin/bash
# step3_02 - do the MAG-dark reps have a homolog in a noise-free cultured genome?
# Query = the 4,292,446 genome-resolved MAG-dark protein reps (20 chunks), target = the
# 3,024,784 cultured Bakta proteins of all 722 genomes. e-5, identity and both coverages reported.
# Query side is the dark rep so that "dark rep with >=1 cultured homolog" is read straight
# off column 1, and -k 500 keeps every cultured genome that hits (the shared-species split
# needs all of them, not just the best). Coverage is NOT filtered here: qcovhsp/scovhsp are
# reported and the reciprocal >= 50 % rule is applied in step3_04, so the relaxed count
# (a full-length cultured protein matching a MAG fragment) stays available as a sensitivity.
#   sbatch scripts/step3_02_dark_vs_cultured.sbatch
#SBATCH --job-name=ale_s3_dia
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/cultured_comparison_clean/logs/step3_dia_%A_%a.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/cultured_comparison_clean/logs/step3_dia_%A_%a.err
#SBATCH --array=0-20%21
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --partition=batch
set -euo pipefail
source /ibex/scratch/projects/c2014/rob/miniconda3/etc/profile.d/conda.sh
conda activate eggnog

P=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/cultured_comparison_clean
OUT=$P/tmp/dia_hits; mkdir -p "$OUT"
# task 20 = the matched control: the 106,042 dark family reps PLUS the 106,042 length- and
# size-matched ANNOTATED family reps of the section 3.1 apparatus, in one file, so the dark
# rate has an annotated-protein ceiling measured on the same instrument.
if [ "$SLURM_ARRAY_TASK_ID" -eq 20 ]; then
    Q=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/specificity_clean/results/query_all.faa
    M8=$OUT/hits_ctrl.m8
else
    printf -v IDX "%02d" "$SLURM_ARRAY_TASK_ID"
    Q=$P/tmp/dark_chunks/chunk_${IDX}.faa
    M8=$OUT/hits_${IDX}.m8
fi

[ -s "$Q" ] || { echo "ERROR: $Q missing"; exit 1; }
if [ -s "$M8" ]; then echo "[$(date -Iseconds)] $(basename "$M8") already done - skip"; exit 0; fi

echo "[$(date -Iseconds)] $(basename "$M8"): $(grep -c '^>' "$Q") queries"
diamond blastp \
    -q "$Q" -d "$P/tmp/cultured_db.dmnd" -o "${M8}.tmp" \
    --very-sensitive -e 1e-5 -k 500 \
    --outfmt 6 qseqid sseqid pident length qlen slen evalue bitscore qcovhsp scovhsp \
    --threads "${SLURM_CPUS_PER_TASK:-32}" -b 6 -c 1 --tmpdir "$OUT"
mv "${M8}.tmp" "$M8"
echo "[$(date -Iseconds)] $(basename "$M8") done - $(wc -l < "$M8") hit rows"
echo "DONE_STEP3_02_${SLURM_ARRAY_TASK_ID}"
