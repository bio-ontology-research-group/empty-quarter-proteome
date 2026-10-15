#!/bin/bash
#SBATCH --job-name=eq-helio-fold
#SBATCH --partition=freecycle-h200
#SBATCH --qos=freecycle
#SBATCH --gpus=nvidia_h200:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --requeue
set -euo pipefail
B=/mnt/data/u/hohndor/eq-family-followup-20261015
export HF_HOME=/mnt/data/u/hohndor/folding/hf
export OMP_NUM_THREADS=4
nvidia-smi --query-gpu=name,memory.total --format=csv
srun /home/hohndor/pixi-envs/esmfold/.pixi/envs/default/bin/python -u "$B/scripts/fold_helio.py" "$B/inputs/${FASTA:-helio_pilot.faa}" "$B/${OUTPUT:-pilot}"
