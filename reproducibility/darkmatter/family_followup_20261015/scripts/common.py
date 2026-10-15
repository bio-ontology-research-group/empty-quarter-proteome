import csv, hashlib, os
from pathlib import Path
W=Path('/ibex/scratch/projects/c2014/rob/metagenomics')
D=W/'06_darkmatter/family_discovery_20261015_clean'
N=W/'06_darkmatter/family_followup_20261015_clean'
PY='/ibex/scratch/projects/c2014/rob/envs/gspa-neural/bin/python'
assert os.environ.get('SLURM_JOB_ID'), 'Run analyses on a Slurm worker'
def tab(p):
    with open(p) as f:return list(csv.DictReader(f,delimiter='\t'))
def write(p,rows,fields=None):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    with open(p,'w') as f:
        w=csv.DictWriter(f,fieldnames=fields or list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        while b:=f.read(4*1024*1024):h.update(b)
    return h.hexdigest()
