#!/bin/bash
#SBATCH --job-name=ale_refb23
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/06_reference_bias/refbias_clean/logs/refb23_%j.out
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=06:00:00
#SBATCH --partition=batch
set -euo pipefail
cd /ibex/scratch/projects/c2014/rob/metagenomics/06_reference_bias/refbias_clean/scripts
echo "[refbias] 2/3 MAG-encoded KO matrix"; python3 refbias_2_magko.py
echo "[refbias] 3/3 PICRUSt2 comparison";   python3 refbias_3_compare.py
echo "[refbias] done"
