#!/bin/bash
#SBATCH --job-name=h3_ref_pilot
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p6_hypotheses_clean/h3_references/logs/pilot_%j.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p6_hypotheses_clean/h3_references/logs/pilot_%j.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=00:30:00
#SBATCH --partition=debug
#SBATCH --account=pi-hohndor

set -euo pipefail

WORKDIR=/ibex/scratch/projects/c2014/rob/metagenomics
H3DIR=$WORKDIR/06_darkmatter/p6_hypotheses_clean/h3_references
PY=/ibex/scratch/projects/c2014/rob/miniconda3/envs/metabolism-env/bin/python
export PATH=/ibex/scratch/projects/c2014/rob/miniconda3/envs/eggnog/bin:$PATH

echo "[$(date -Iseconds)] Starting Hypothesis 3 Subsample Search Pilot on IBEX..."
echo "Node: $(hostname)"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Python: $($PY --version)"
echo "DIAMOND: $(diamond --version 2>&1 | head -n 1)"

# Step 1: Generate query sample if not already done
if [ ! -s "$H3DIR/results/pilot_query_2000.faa" ]; then
    echo "[$(date -Iseconds)] Generating 2,000 query sample..."
    $PY $H3DIR/scripts/02_sample_queries.py
fi

# Step 2: Run empirical subsample search and top-1 bias evaluation
echo "[$(date -Iseconds)] Running 03_pilot_subsample_search.py..."
$PY $H3DIR/scripts/03_pilot_subsample_search.py $SLURM_CPUS_PER_TASK

echo "[$(date -Iseconds)] Pilot job completed successfully."
