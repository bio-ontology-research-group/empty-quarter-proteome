#!/usr/bin/env bash
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p4_review_20261013/logs/round4-verify-%j.log
set -euo pipefail
cd /ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p4_review_20261013/review_package
PY=/ibex/scratch/projects/c2014/rob/miniconda3/envs/metabolism-env/bin/python
"$PY" reproducibility/darkmatter/verify_review_inputs.py
"$PY" reproducibility/darkmatter/current_revision_checks.py
"$PY" reproducibility/darkmatter/round4_checks.py
"$PY" reviews/darkmatter-2026-10-13/round4/root_cross_table_audit.py
"$PY" - <<'PY'
from pathlib import Path
import csv,hashlib
rows=list(csv.DictReader(open('reviews/darkmatter-2026-10-13/round4/final_artifacts.sha256.tsv'),delimiter='\t'))
for r in rows:assert hashlib.sha256(Path(r['path']).read_bytes()).hexdigest()==r['sha256'],r['path']
print('PASS:',len(rows),'final manuscript, figure and analysis-script hashes match local package')
PY
