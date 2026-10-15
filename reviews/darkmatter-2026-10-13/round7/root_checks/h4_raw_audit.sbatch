#!/bin/bash
#SBATCH --job-name=h4_raw_audit
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p7_analyses_clean/h4_structure/logs/raw_audit_%j.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p7_analyses_clean/h4_structure/logs/raw_audit_%j.err
set -euo pipefail
/ibex/scratch/projects/c2014/rob/miniconda3/envs/metabolism-env/bin/python /ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p7_analyses_clean/h4_structure/scripts_recovery_v1/audit_h4_raw.py
