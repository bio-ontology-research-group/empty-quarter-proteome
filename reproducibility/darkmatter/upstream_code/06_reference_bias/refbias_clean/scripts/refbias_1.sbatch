#!/bin/bash
# §2.1 step 1 only: no-KO ceiling on the 5,438 catalogue.
#SBATCH --job-name=ale_ceiling
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/06_reference_bias/refbias_clean/logs/ceiling_%j.out
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --partition=batch
set -euo pipefail
cd /ibex/scratch/projects/c2014/rob/metagenomics/06_reference_bias/refbias_clean/scripts
python3 refbias_1_ceiling.py 16
