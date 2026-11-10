import os,socket,json,hashlib,csv
from pathlib import Path
W=Path('/ibex/scratch/projects/c2014/rob/metagenomics')
N=W/'06_darkmatter/heat_stress_20261016'
T=W/'10_thermal'
assert os.environ.get('SLURM_JOB_ID') and 'login' not in socket.gethostname(), 'Slurm worker required'
def tab(p):
 with open(p) as f:return list(csv.DictReader(f,delimiter='\t'))
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
 return h.hexdigest()
def write(p,rows,fields=None):
 p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
 with open(p,'w') as f:
  wr=csv.DictWriter(f,fieldnames=fields or list(rows[0]),delimiter='\t');wr.writeheader();wr.writerows(rows)
