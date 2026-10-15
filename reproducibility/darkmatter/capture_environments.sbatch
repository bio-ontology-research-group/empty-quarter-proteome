#!/bin/bash
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:15:00
set -euo pipefail
OUT=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p4_review_20261013/results/environments
mkdir -p "$OUT"
/ibex/scratch/projects/c2014/rob/miniconda3/envs/metabolism-env/bin/python - "$OUT" <<'PY'
import sys,json,glob,hashlib,platform
from pathlib import Path
out=Path(sys.argv[1]);base=Path('/ibex/scratch/projects/c2014/rob/miniconda3/envs');roots=[base/x for x in ['eggnog','pyrodigal','mmseqs2','selection','foldseek','metabolism-env']]+[Path('/home/velazqal/miniconda3/envs/python-utils')]
rows=[]
for env in roots:
 packages=[];error=''
 try: paths=sorted((env/'conda-meta').glob('*.json'))
 except PermissionError as e: paths=[];error=str(e)
 for p in paths:
  d=json.loads(p.read_text());packages.append({k:d.get(k) for k in ['name','version','build','subdir','sha256','md5']})
 (out/(env.name+'.json')).write_text(json.dumps({'environment':str(env),'packages':packages,'error':error,'scope':'observed environment at audit time, not proof of historical execution'},indent=2)+'\n')
for p in [Path('/home/hohndor/bin/minimap2'),base/'eggnog/bin/hmmsearch',base/'eggnog/bin/diamond',base/'foldseek/bin/foldseek',base/'mmseqs2/bin/mmseqs']:
 if p.exists():rows.append({'path':str(p),'resolved_path':str(p.resolve()),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
(out/'binary_sha256.json').write_text(json.dumps(rows,indent=2)+'\n');(out/'python.json').write_text(json.dumps({'python':sys.version,'platform':platform.platform()},indent=2)+'\n')
PY
