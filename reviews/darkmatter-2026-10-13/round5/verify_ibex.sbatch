#!/usr/bin/env bash
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p4_review_20261013/logs/round5-verify-%j.log
set -euo pipefail
cd /ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p4_review_20261013/review_package
PY=/ibex/scratch/projects/c2014/rob/miniconda3/envs/metabolism-env/bin/python
"$PY" reproducibility/darkmatter/verify_review_inputs.py
"$PY" reproducibility/darkmatter/current_revision_checks.py
"$PY" reproducibility/darkmatter/round4_checks.py
"$PY" reproducibility/darkmatter/round5_extensions/checks.py
"$PY" reviews/darkmatter-2026-10-13/round4/root_cross_table_audit.py
"$PY" - <<'PY'
from pathlib import Path
import csv,hashlib,json
for r in csv.DictReader(open('reviews/darkmatter-2026-10-13/round5/final_artifacts.sha256.tsv'),delimiter='\t'):
 assert hashlib.sha256(Path(r['path']).read_bytes()).hexdigest()==r['sha256'],r['path']
print('PASS: final manuscript, figures, code and review hashes match local package')
pins={}
for p in ['profile_input_sha256.json','culture_caller/input_sha256.json','temporal_visits_input_sha256.json']:
 for name,want in json.loads((Path('reproducibility/darkmatter/round5_extensions/results')/p).read_text()).items():
  if name in pins:assert pins[name]==want
  pins[name]=want
for name,want in pins.items():
 h=hashlib.sha256()
 with open(name,'rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 assert h.hexdigest()==want,name
print('PASS:',len(pins),'remote extension input/generated-GFF hashes')
PY
