#!/bin/bash
#SBATCH --job-name=smag_part_verify
#SBATCH --account=pi-hohndor
#SBATCH --partition=debug
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:05:00
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p7_analyses_clean/h3_references/logs/verify_transfer_%j.out
set -euo pipefail
cd /ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p7_analyses_clean/h3_references
python3 - <<'PY'
from pathlib import Path
import hashlib,json,os
p=Path('downloads/smag_parts/mag.tar.gz.dw.unimatrix_verified');b=p.read_bytes();h=hashlib.sha256(b).hexdigest();assert len(b)==104857600;assert h=='d77379eacf826adb7105754c26149d843e94cc3aa4ebeade3257c046fdf56325';os.replace(p,p.with_name('mag.tar.gz.dw'))
x={'passed':True,'part':'mag.tar.gz.dw','sha256':h,'bytes':len(b),'origin':'unimatrix01 Slurm9183 direct Zenodo benchmark; source checked publishedMD5'};Path('results/unimatrix_part_transfer_verified.json').write_text(json.dumps(x,indent=2)+'\n');print(json.dumps(x))
PY
