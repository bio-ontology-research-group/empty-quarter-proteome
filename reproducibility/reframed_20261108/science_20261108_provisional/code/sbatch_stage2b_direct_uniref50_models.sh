#!/bin/bash
#SBATCH --job-name=eq_sci_s2b_u50
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=192
#SBATCH --exclusive
#SBATCH --time=00:45:00
#SBATCH --output=/scratch/project/k10018/eq-science-20261108/logs/eq_sci_s2b_u50_%j.out
#SBATCH --error=/scratch/project/k10018/eq-science-20261108/logs/eq_sci_s2b_u50_%j.err

set -euo pipefail
export EQ_SCIENCE_OUT="/scratch/project/k10018/eq-science-20261108"
cd "$EQ_SCIENCE_OUT"
source /project/k10018/sw/eq-env.sh
conda activate /project/k10018/sw/envs/mapping

python3 -m pytest -o cache_dir="$EQ_SCIENCE_OUT/logs/.pytest_cache" "$EQ_SCIENCE_OUT/tests/test_science_20261108.py" -v

if [[ ! -s "$EQ_SCIENCE_OUT/results/task_b_uniref50_direct_per_gene_484.tsv.gz" ]]; then
    echo "ERROR: Missing or empty required Stage 2 direct per-gene output: $EQ_SCIENCE_OUT/results/task_b_uniref50_direct_per_gene_484.tsv.gz" >&2
    exit 1
fi

python3 "$EQ_SCIENCE_OUT/code/stage2b_join_external_signalp_and_direct_uniref50.py" --allow-pending-signalp
