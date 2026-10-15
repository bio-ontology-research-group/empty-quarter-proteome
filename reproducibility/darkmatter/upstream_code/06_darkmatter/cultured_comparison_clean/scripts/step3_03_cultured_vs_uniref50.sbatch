#!/bin/bash
# step3_03 - clean dark rate: the full cultured proteomes against UniRef50 under the SAME rule that
# defines MAG darkness (darkclean_3b_uniref.sh, verbatim params: -s 7.5 -e 1e-3 --max-seqs 10,
# dark = no hit, no coverage filter). 14 chunks of the 3,024,784 tagged cultured proteins: chunks
# 00-11 of ~226k hold the 651 enrichment and isolate genomes, chunks 12-13 of ~157k the 71
# soil-batch genomes. The hypotheticals are a subset of these ids, so they need no separate search.
# Finished chunks are skipped, so after new genomes are added only the new chunks need searching:
#   sbatch --array=12-13 scripts/step3_03_cultured_vs_uniref50.sbatch
#SBATCH --job-name=ale_s3_uniref
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/cultured_comparison_clean/logs/step3_uniref_%A_%a.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/cultured_comparison_clean/logs/step3_uniref_%A_%a.err
#SBATCH --array=0-13%14
#SBATCH --cpus-per-task=16
#SBATCH --mem=256G
#SBATCH --time=24:00:00
#SBATCH --partition=batch
set -euo pipefail
source /ibex/scratch/projects/c2014/rob/miniconda3/etc/profile.d/conda.sh
conda activate mmseqs2

M=/ibex/scratch/projects/c2014/rob/metagenomics
UNIREF_DB=$M/04_structure/uniref50/uniref50
P=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/cultured_comparison_clean
OUT=$P/tmp/uniref_hits; mkdir -p "$OUT"
printf -v IDX "%02d" "$SLURM_ARRAY_TASK_ID"
CHUNK_FA=$P/tmp/cultured_chunks/chunk_${IDX}.faa
HITS_TSV=$OUT/hits_chunk_${IDX}.tsv

[ -s "$CHUNK_FA" ] || { echo "ERROR: chunk $CHUNK_FA missing"; exit 1; }
[ -f "${UNIREF_DB}.dbtype" ] || { echo "ERROR: UniRef50 DB not found at $UNIREF_DB"; exit 2; }
if [ -s "$HITS_TSV" ]; then echo "[$(date -Iseconds)] chunk $IDX already done - skip"; exit 0; fi

WORK=$OUT/chunk_${IDX}; mkdir -p "$WORK/tmp"
QDB=$WORK/query.mmseqsdb
RESULT_DB=$WORK/result.mmseqsdb

if [ ! -f "${QDB}.dbtype" ]; then
    echo "[$(date -Iseconds)] chunk $IDX createdb"
    mmseqs createdb "$CHUNK_FA" "$QDB"
fi

if [ ! -f "${RESULT_DB}.dbtype" ]; then
    echo "[$(date -Iseconds)] chunk $IDX search vs UniRef50"
    mmseqs search "$QDB" "$UNIREF_DB" "$RESULT_DB" "$WORK/tmp" \
        -s 7.5 -e 1e-3 --max-seqs 10 --threads 16 --split-memory-limit 200G
fi

echo "[$(date -Iseconds)] chunk $IDX convertalis"
mmseqs convertalis "$QDB" "$UNIREF_DB" "$RESULT_DB" "${HITS_TSV}.tmp" \
    --format-output 'query,target,fident,alnlen,evalue,bits' --threads 16
mv "${HITS_TSV}.tmp" "$HITS_TSV"
rm -rf "$WORK/tmp"
echo "[$(date -Iseconds)] chunk $IDX done - $(wc -l < "$HITS_TSV") hit rows"
echo "DONE_STEP3_03_${IDX}"
