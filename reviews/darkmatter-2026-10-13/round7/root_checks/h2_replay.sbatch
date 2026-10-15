#!/bin/bash
#SBATCH --job-name=h2_review_replay
#SBATCH --account=pi-hohndor
#SBATCH --partition=debug
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p7_analyses_clean/h2_context/logs/replay_%j.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p7_analyses_clean/h2_context/logs/replay_%j.err
set -euo pipefail
cd /ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p7_analyses_clean/h2_context
PY=/ibex/scratch/projects/c2014/rob/miniconda3/envs/metabolism-env/bin/python
"$PY" -u - <<'PY'
from pathlib import Path
import shutil,subprocess,hashlib,json,sys
src=Path('results_full');dst=Path('results_review_replay');shutil.copytree(src,dst)
subprocess.run([sys.executable,'scripts/03_analyse_conservation.py',str(dst)],check=True)
rows=[]
for p in sorted(src.iterdir()):
 if p.is_file():
  q=dst/p.name;a=hashlib.sha256(p.read_bytes()).hexdigest();b=hashlib.sha256(q.read_bytes()).hexdigest();rows.append({'path':p.name,'before':a,'after':b,'identical':a==b})
out={'passed':all(x['identical'] for x in rows),'files':rows};Path('review_replay_verification.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out));assert out['passed']
PY
