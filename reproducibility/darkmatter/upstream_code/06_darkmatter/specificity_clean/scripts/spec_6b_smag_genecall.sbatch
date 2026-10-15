#!/bin/bash
# §3.1 SMAG leg, step b: gene calling. Pyrodigal 3.7.1 (rob/miniconda3 env `pyrodigal`), the same env
# and flags as the project's own assembly gene calls (scripts/ibex/22_pyrodigal_assembly_array.sh):
#   pyrodigal -p meta -j 4 --no-stop-codon
# Headers are prefixed <MAG id>:: so every protein is globally unique and traceable to its MAG.
# Chunk k = rows (k-1)*500+1 .. k*500 of tmp/smag/genome_list.txt -> tmp/smag_proteins/chunk_<k>.faa.
# A chunk file only appears (atomic mv) once all its genomes are done, so a rerun skips it.
#   N=$(( ($(wc -l < tmp/smag/genome_list.txt) + 499) / 500 ))
#   sbatch --array=1-${N} spec_6b_smag_genecall.sbatch
#SBATCH --job-name=ale_spec_smag_genecall
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/specificity_clean/logs/spec_smag_genecall_%A_%a.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/specificity_clean/logs/spec_smag_genecall_%A_%a.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --partition=batch
#SBATCH --account=pi-hohndor
set -euo pipefail
S=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/specificity_clean
LIST=$S/tmp/smag/genome_list.txt
OUTD=$S/tmp/smag_proteins
CHUNK=500
K=${SLURM_ARRAY_TASK_ID:?array task id required}
source /ibex/scratch/projects/c2014/rob/miniconda3/etc/profile.d/conda.sh
conda activate pyrodigal
mkdir -p "$OUTD"
OUT=$(printf '%s/chunk_%03d.faa' "$OUTD" "$K")
[ -s "$OUT" ] && { echo "chunk $K already done"; exit 0; }

W=$(mktemp -d -p "$OUTD" .work_${K}_XXXX)
trap 'rm -rf "$W"' EXIT
echo "[$(date -Iseconds)] chunk $K on $(hostname), pyrodigal $(pyrodigal --version 2>&1 | tail -1)"
n=0
while read -r fa; do
  mag=$(basename "$fa"); mag=${mag%.gz}; mag=${mag%.*}
  in=$fa
  if [[ $fa == *.gz ]]; then zcat "$fa" > "$W/in.fa"; in=$W/in.fa; fi
  pyrodigal -p meta -i "$in" -a "$W/g.faa" -o /dev/null -j 4 --no-stop-codon
  awk -v s="$mag" '/^>/ { sub(/^>/, ">"s"::") } { print }' "$W/g.faa" >> "$W/chunk.faa"
  n=$((n+1))
done < <(sed -n "$(( (K-1)*CHUNK + 1 )),$(( K*CHUNK ))p" "$LIST")
[ "$n" -eq 0 ] && { echo "no genomes for chunk $K"; exit 0; }
mv "$W/chunk.faa" "$OUT"
echo "[$(date -Iseconds)] chunk $K done: $n genomes, $(grep -c '^>' "$OUT") proteins"
