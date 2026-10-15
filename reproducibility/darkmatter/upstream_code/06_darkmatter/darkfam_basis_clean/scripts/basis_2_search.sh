#!/bin/bash
# darkfam_basis_clean step 2 — strict UniRef50 search of the 9 query chunks under the
# paper's locked dark rule: mmseqs -s 7.5 -e 1e-3 --max-seqs 10 (same params + same DB
# as darkclean_3b / 27_uniref50_search_delta). Task 1 = the 23,382 family reps (fast,
# decides the survivor universe); tasks 2-9 = the never-screened gap reps.
# Idempotent per chunk.
#
# Submit (after basis_1): sbatch --array=1-9 scripts/basis_2_search.sh
#
#SBATCH --job-name=ale_basis_search
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/darkfam_basis_clean/logs/search_%A_%a.out
#SBATCH --cpus-per-task=32
#SBATCH --mem=240G
#SBATCH --time=24:00:00
#SBATCH --partition=batch

set -euo pipefail
B=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/darkfam_basis_clean
UNIREF=/ibex/scratch/projects/c2014/rob/metagenomics/04_structure/uniref50/uniref50
R=$B/results
i=$SLURM_ARRAY_TASK_ID
CH=$R/chunks/chunk_${i}.faa
HITS=$R/hits_chunk_${i}.tsv

[ -f "$R/CHUNKS_DONE" ] || { echo "ERROR: run basis_1_prep first"; exit 1; }
[ -s "$CH" ] || { echo "ERROR: $CH missing"; exit 1; }
[ -s "$HITS" ] && { echo "chunk $i already done — skip"; exit 0; }

source /ibex/scratch/projects/c2014/rob/miniconda3/etc/profile.d/conda.sh
conda activate mmseqs2
mmseqs --version

W=$R/work_${i}; mkdir -p "$W"
echo "[$(date -Iseconds)] chunk $i: $(grep -c '^>' "$CH") queries"
mmseqs createdb "$CH" "$W/qdb"
mmseqs search "$W/qdb" "$UNIREF" "$W/res" "$W/tmp" \
    -s 7.5 -e 1e-3 --max-seqs 10 \
    --threads "$SLURM_CPUS_PER_TASK" --split-memory-limit 200G
mmseqs convertalis "$W/qdb" "$UNIREF" "$W/res" "$HITS.part" \
    --format-output 'query,target,pident,evalue,bits' --threads "$SLURM_CPUS_PER_TASK"
mv "$HITS.part" "$HITS"
rm -rf "$W"
echo "[$(date -Iseconds)] chunk $i done: $(cut -f1 "$HITS" | sort -u | wc -l) queries with >=1 hit"
