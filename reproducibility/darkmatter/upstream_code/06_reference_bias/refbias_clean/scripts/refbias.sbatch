#!/bin/bash
# §2.1 reference-bias: no-KO ceiling on the 5,438 catalogue (step 1), then the MAG-encoded
# per-sample KO matrix (step 2) and the PICRUSt2 comparison (step 3).
#SBATCH --job-name=ale_refbias
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/06_reference_bias/refbias_clean/logs/refbias_%j.out
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --partition=batch
set -euo pipefail
S=/ibex/scratch/projects/c2014/rob/metagenomics/06_reference_bias/refbias_clean/scripts
cd "$S"
echo "[refbias] 1/3 no-KO ceiling"; python3 refbias_1_ceiling.py 16
echo "[refbias] 2/3 MAG-encoded KO matrix"; python3 refbias_2_magko.py
echo "[refbias] 3/3 PICRUSt2 comparison"; python3 refbias_3_compare.py
echo "[refbias] all done"
