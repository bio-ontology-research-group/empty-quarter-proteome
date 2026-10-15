#!/bin/bash
#SBATCH --job-name=eq-helio-structure-audit
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch,debug
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:30:00
set -euo pipefail
cd /ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/family_followup_20261015_clean
mkdir -p tools
if [ ! -d tools/USalign ]; then git clone --depth 1 https://github.com/pylelab/USalign.git tools/USalign; fi
git -C tools/USalign rev-parse HEAD > results/USalign_commit.txt
if [ ! -x tools/USalign/TMalign ]; then make -C tools/USalign TMalign; fi
/ibex/scratch/projects/c2014/rob/envs/gspa-neural/bin/python scripts/structure_audit.py
