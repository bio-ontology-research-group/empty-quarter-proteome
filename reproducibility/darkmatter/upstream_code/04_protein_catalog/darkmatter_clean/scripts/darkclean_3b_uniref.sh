#!/bin/bash
# darkclean_3b — UniRef50 sensitive search of the uncovered-rep chunks, as a 120-task array.
# Mirrors 27_uniref50_search_delta.sh (same params: -s 7.5 -e 1e-3 --max-seqs 10). Idempotent per chunk.
#
# Submit (after darkclean_3a): sbatch .../darkmatter_clean/scripts/darkclean_3b_uniref.sh
#
#SBATCH --job-name=velazqal_darkclean3b
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/logs/darkclean3b_%A_%a.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/logs/darkclean3b_%A_%a.err
#SBATCH --array=0-119%120
#SBATCH --cpus-per-task=16
#SBATCH --mem=256G
#SBATCH --time=24:00:00
#SBATCH --partition=batch

set -euo pipefail
source /ibex/scratch/projects/c2014/rob/miniconda3/etc/profile.d/conda.sh
conda activate mmseqs2

M=/ibex/scratch/projects/c2014/rob/metagenomics
CAT=$M/04_protein_catalog
OUT=$CAT/darkmatter_clean/03_uncovered_vs_uniref50
CHUNK_DIR=$OUT/chunks
UNIREF_DB=$M/04_structure/uniref50/uniref50
PER_CHUNK_DIR=$OUT/per_chunk
mkdir -p "$PER_CHUNK_DIR"

printf -v IDX "%03d" "$SLURM_ARRAY_TASK_ID"
CHUNK_FA=$CHUNK_DIR/chunk_${IDX}.faa.gz
HITS_TSV=$PER_CHUNK_DIR/hits_chunk_${IDX}.tsv

[ -s "$CHUNK_FA" ] || { echo "ERROR: chunk $CHUNK_FA missing"; exit 1; }
[ -f "${UNIREF_DB}.dbtype" ] || { echo "ERROR: UniRef50 DB not found at $UNIREF_DB"; exit 2; }
if [ -s "$HITS_TSV" ]; then echo "[$(date -Iseconds)] chunk $IDX already done — skip"; exit 0; fi

WORK=$PER_CHUNK_DIR/chunk_${IDX}
mkdir -p "$WORK/tmp"
QDB=$WORK/query.mmseqsdb
RESULT_DB=$WORK/result.mmseqsdb

if [ ! -f "${QDB}.dbtype" ]; then
    echo "[$(date -Iseconds)] chunk $IDX createdb"
    FIFO=$(mktemp -u --tmpdir=/dev/shm dc_${IDX}_fifo.XXXX); mkfifo "$FIFO"
    pigz -dc "$CHUNK_FA" > "$FIFO" & PID=$!
    mmseqs createdb "$FIFO" "$QDB"; wait $PID; rm -f "$FIFO"
fi

if [ ! -f "${RESULT_DB}.dbtype" ]; then
    echo "[$(date -Iseconds)] chunk $IDX search vs UniRef50"
    mmseqs search "$QDB" "$UNIREF_DB" "$RESULT_DB" "$WORK/tmp" \
        -s 7.5 -e 1e-3 --max-seqs 10 --threads 16 --split-memory-limit 200G
fi

echo "[$(date -Iseconds)] chunk $IDX convertalis"
mmseqs convertalis "$QDB" "$UNIREF_DB" "$RESULT_DB" "${HITS_TSV}.tmp" \
    --format-output 'query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits' \
    --threads 16
mv "${HITS_TSV}.tmp" "$HITS_TSV"
rm -rf "$WORK/tmp"
echo "[$(date -Iseconds)] chunk $IDX done — $(wc -l < "$HITS_TSV") hit rows"
