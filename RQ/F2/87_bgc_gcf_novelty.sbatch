#!/bin/bash
# RQ-F2 : Novelty of Empty Quarter BGCs vs known biosynthesis (MIBiG / BiG-FAM).
#   Turn the antiSMASH BGC *census* into a *discovery* claim: co-cluster all
#   119,090 antiSMASH-8.0.4 BGC regions from 15,469 EQ MAGs together with the
#   2,502 MIBiG 3.1 reference BGCs into Gene Cluster Families (GCFs) with
#   BiG-SLICE v2.0.2 (BIRCH centroids, normalized features, threshold T=0.4).
#   A GCF that contains >=1 EQ BGC and 0 MIBiG BGCs = NOVEL family (no
#   experimentally-characterized analog).  MIBiG underpins BiG-FAM, so this is
#   the self-contained version of the BiG-FAM novelty comparison.
#   Refs: kautsar2020bigslice; kautsar2021bigfam; blin2023antismash; terlouw2023mibig.
#
#   Tool choice: 119k BGCs is far beyond a tractable BiG-SCAPE all-vs-all;
#   BiG-SLICE was designed for >1M BGCs (super-linear feature-vector clustering).
#
#   Submit:
#     sbatch /ibex/scratch/projects/c2014/rob/metagenomics/scripts/ibex/87_bgc_gcf_novelty.sbatch
#
#SBATCH --job-name=hohndor_bgc_gcf
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/logs/bgc_gcf_%j.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/logs/bgc_gcf_%j.err
#SBATCH --cpus-per-task=32
#SBATCH --mem=300G
#SBATCH --time=2-00:00:00
#SBATCH --partition=batch
#SBATCH --account=pi-hohndor
set -uo pipefail
export M=/ibex/scratch/projects/c2014/rob/metagenomics
BASE=$M/06_metabolism/bgc_novelty
INPUT=$BASE/bigslice_input
OUTPUT=$BASE/bigslice_out
MODELS=/ibex/scratch/projects/c2014/rob/miniconda3/envs/bigslice/bin/bigslice-models
mkdir -p "$M/logs" "$BASE/results"

source /ibex/scratch/projects/c2014/rob/miniconda3/etc/profile.d/conda.sh
conda activate bigslice
# conda scipy/highs needs the env libstdc++ ahead of the system one
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:${LD_LIBRARY_PATH:-}
THREADS=${SLURM_CPUS_PER_TASK:-32}
echo "[$(date -Iseconds)] host=$(hostname) cpus=$THREADS bigslice=2.0.2 input=$INPUT"

# bigslice writes into a fresh OUTPUT (use --resume to continue a partial run)
rm -rf "$OUTPUT"
bigslice -i "$INPUT" "$OUTPUT" \
    --threshold 0.4 \
    -t "$THREADS" \
    --program_db_folder "$MODELS"
rc=$?
echo "[$(date -Iseconds)] bigslice clustering rc=$rc"
[ "$rc" -eq 0 ] || exit "$rc"

echo "[$(date -Iseconds)] running novelty analysis"
python "$M/scripts/ibex/87_bgc_novelty_analyze.py" \
    "$OUTPUT/result/data.db" "$BASE/results"

echo "[$(date -Iseconds)] RQ-F2 job complete"
