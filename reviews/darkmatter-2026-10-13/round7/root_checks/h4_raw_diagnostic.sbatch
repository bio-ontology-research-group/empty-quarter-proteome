#!/bin/bash
#SBATCH --job-name=h4_finalize_audit
#SBATCH --account=pi-hohndor
#SBATCH --partition=debug
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p7_analyses_clean/h4_structure/logs/raw_diagnostic_%j.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p7_analyses_clean/h4_structure/logs/raw_diagnostic_%j.err
set -euo pipefail
/ibex/scratch/projects/c2014/rob/miniconda3/envs/metabolism-env/bin/python -u /ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p7_analyses_clean/h4_structure/scripts_recovery_v3/audit_h4_raw_diagnostic.py
