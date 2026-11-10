#!/bin/bash
#SBATCH --job-name=eq_sci_stg2b
#SBATCH --account=k10018
#SBATCH --partition=workq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=192
#SBATCH --exclusive
#SBATCH --time=00:45:00
#SBATCH --output=/scratch/project/k10018/eq-science-20261108/logs/eq_sci_stg2b_%j.out
#SBATCH --error=/scratch/project/k10018/eq-science-20261108/logs/eq_sci_stg2b_%j.err

set -euo pipefail
export EQ_SCIENCE_OUT="/scratch/project/k10018/eq-science-20261108"
cd "$EQ_SCIENCE_OUT"
source /project/k10018/sw/eq-env.sh
conda activate /project/k10018/sw/envs/mapping

python3 -m pytest -o cache_dir="$EQ_SCIENCE_OUT/logs/.pytest_cache" "$EQ_SCIENCE_OUT/tests/test_science_20261108.py" -v

for req_f in \
    "$EQ_SCIENCE_OUT/results/task_b_uniref50_direct_per_gene_484.tsv.gz" \
    "$EQ_SCIENCE_OUT/results/external247_signalp6_predictions.tsv.gz" \
    "$EQ_SCIENCE_OUT/results/external_signalp_manifest.tsv"; do
    if [[ ! -s "$req_f" ]]; then
        echo "ERROR: Full Stage 2b join requires non-empty $req_f (never falling back to taxonomy-only)." >&2
        exit 1
    fi
done

python3 "$EQ_SCIENCE_OUT/code/stage2b_join_external_signalp_and_direct_uniref50.py"
