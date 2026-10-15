#!/bin/bash
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
set -euo pipefail
P=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p4_review_20261013/review_package
PY=/ibex/scratch/projects/c2014/rob/miniconda3/envs/metabolism-env/bin/python
cd "$P"
"$PY" reproducibility/darkmatter/verify_review_inputs.py
"$PY" reproducibility/darkmatter/current_revision_checks.py
"$PY" - <<'PY'
import csv,hashlib
from pathlib import Path
for r in csv.DictReader(open('reviews/darkmatter-2026-10-13/final_artifacts.sha256.tsv'),delimiter='\t'):
 p=Path(r['path']);assert hashlib.sha256(p.read_bytes()).hexdigest()==r['sha256'],str(p)
print('Final manuscript, figures and analysis-script hashes match local manifest')
PY
